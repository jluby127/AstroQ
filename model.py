import astropy as ap
import astroplan as apl
import numpy as np
from astropy.time import Time
from astropy.time import TimeDelta
from collections import defaultdict
from itertools import permutations
import gurobipy as gp
from gurobipy import *


class TTPModel(object):

    def __init__(self,start,stop,stars,observatory):

        # Will have observatory object
        self.observatory = observatory
        self.stars = stars
        self.nightstarts = start
        self.nightends = stop
        self.create_nodes()
        self.compute_tau_slew()
        self.solve()

    def create_nodes(self):
        dur = np.round((self.nightends-self.nightstarts).jd*24*60,0)*60
        t = np.arange(self.nightstarts.jd,self.nightends.jd,TimeDelta(60,format='sec'))
        t = Time(t,format='jd')

        # Make a lookup table for node index -> star index
        # Also track which targets have degenerate nodes, and where they are
        node_to_star = defaultdict(int)
        multi_visit_ind = defaultdict(list)

        i = 1
        for s in self.stars:
            for j in range(s.visits):
                if s.visits > 1:
                    multi_visit_ind[s.index].append(i)
                node_to_star[i] = s.index
                i += 1

        # N is the total number of nodes + 2 for anchors
        N = sum(s.visits for s in self.stars) + 2
        self.N = N

        te = []
        tl = []
        tau_exp = []
        tau_sep = []

        for i in range(N):
            if i == 0 or i == N-1:
                te.append(0)
                tl.append((self.nightends-t[0]).jd*24*60)
                tau_exp.append(0)
                tau_sep.append(0)
            else:
                s = self.stars[node_to_star[i]]
                tau_exp.append(s.exptime)
                tau_sep.append(s.intra_night_cadence)
                AZ = self.observatory.observer.altaz(t, s.target)
                alt=AZ.alt.deg
                az=AZ.az.deg
                good = np.where(self.observatory.isup(alt,az) == 1)
                if len(good > 0):
                    te.append(np.round(((t[good[0]]+TimeDelta(s.exptime*60,format='sec'))-t[0]).jd*24*60,1))
                    if t[good[-1]].jd < self.nightends.jd:
                        tl.append(np.round((t[good[-1]]-t[0]).jd*24*60,1))
                    elif t[good[-1]].jd >= self.nightends.jd:
                        tl.append(np.round((self.nightends-t[0]).jd*24*60,1))
                else:
                    print('Target {} does not meet observability requirements'.format(s.name))
                    te.append(0)
                    tl.append(0)
        
        self.te = np.array(te)
        self.tl = np.array(tl)
        self.tau_exp = np.array(tau_exp)
        self.tau_sep = np.array(tau_sep)
        self.node_to_star = node_to_star
        self.multi_visit_ind = multi_visit_ind

    def compute_tau_slew(self,M=3,cad=30):
        # Will want to accomodate slot input value
        self.M = M
        night_dur = (self.nightends.jd - self.nightstarts.jd)*24*60 # minutes

        samples_per_slot = max(night_dur/(M*cad),3) # Evaluate at least every thirty minutes

        slot_bounds = Time(np.linspace(self.nightstarts.jd,self.nightends.jd,M+1,endpoint=True),format='jd')
        
        t = []
        for m in range(M):
            times_to_eval = np.linspace(slot_bounds[m],slot_bounds[m+1],samples_per_slot)
            t.append([item for item in times_to_eval])
        times = Time(t,format='jd')

        nodes = self.node_to_star
        N = self.N
        coordinate_matrix = []

        # Dummy coords for start node
        coordinate_matrix.append([(0,0) for item in times])

        # Loop through targets, evaluate coordinates at all times
        for i in range(N)[1:-1]:
            star = self.stars[nodes[i]]
            altaz = self.observatory.observer.altaz(times,star.target)
            coordinate_matrix.append([(item.alt.deg,item.az.deg) for item in altaz])

        # Dummy for end
        coordinate_matrix.append([(0,0) for item in times])

        def to_wrap_frame(angle):
            if self.wrapangle:
                angle += (360-self.wrapangle)
                if angle >= 360:
                    angle -= 360
            return angle
        
        def max_ang_sep(targ1,targ2,slot):

            slot_ind_start = slot*samples_per_slot
            slot_ind_end = (slot+1)*samples_per_slot-1

            coords1 = coordinate_matrix[targ1,slot_ind_start:slot_ind_end].T
            coords2 = coordinate_matrix[targ2,slot_ind_start:slot_ind_end].T

            alt1 = coords1[0]
            alt2 = coords2[0]                         
            az1 = to_wrap_frame(coords1[1])
            az2 = to_wrap_frame(coords2[1])

            return max(max(np.abs(alt1-alt2)),max(np.abs(az1-az2)))

        tau_slew = defaultdict(float) # Holds minutes
        for m in range(M):
            for targ1,targ2 in permutations(range(N)[1:-1],2):
                tau_slew[(targ1,targ2,m)] = np.round(max_ang_sep(targ1,targ2,m)/(60*self.observatory.slewrate),1)

        # Slot start and end times as minutes from start
        w = (np.array(slot_bounds)-slot_bounds[0])*24*60

        self.w = w
        self.tau_slew = tau_slew


    def to_gurobi_model(self,time_limit=300,opt_gap=0,output_flag=True):
        # Translate TTP Model object into a Gurobi Model object

        Mod = gp.Model('TTP')
        Mod.Params.OutputFlag = output_flag

        N = self.N
        M = self.M
        w = self.w
        tau_slew = self.tau_slew
        tau_exp = self.tau_exp
        te = self.te
        tl = self.tl
        tau_sep = self.tau_sep


        Yi = Mod.addVars(range(N),vtype=GRB.BINARY,name='Yi')
        Xijm = Mod.addVars(range(N),range(N),range(M),vtype=GRB.BINARY,name='Xijm')
        tijm = Mod.addVars(range(N),range(N),range(M),vtype=GRB.CONTINUOUS,name='tijm')
        ti = Mod.addVars(range(N),vtype=GRB.CONTINUOUS,lb=0,name='ti')

        tijmdef = Mod.addConstrs((ti[i] == gp.quicksum(tijm[i,j,m] for j in range(N)[1:] for m in range(M)) 
                            for i in range(N)[:-1]),'tijm_def')
        start_origin = Mod.addConstr(gp.quicksum(Xijm[0,j,m] for j in range(N) for m in range(M)) == 1,
                        'start_origin')
        end_origin = Mod.addConstr(gp.quicksum(Xijm[i,N-1,m] for i in range(N) for m in range(M)) == 1,
                        'end_origin')
        visit_once = Mod.addConstrs((gp.quicksum(Xijm[i,j,m] for i in range(N)[:-1] for m in range(M)) == Yi[j]
                        for j in range(R)[1:]), 'visit_once')
        flow_constr = Mod.addConstrs(((gp.quicksum(Xijm[i,k,m] for i in range(N)[:-1] for m in range(M))
                        - gp.quicksum(Xijm[k,j,m] for j in range(N)[1:] for m in range(M)) == 0)
                        for k in range(N)[:-1][1:]), 'flow_constr')
        exp_constr = Mod.addConstrs((ti[j] >= gp.quicksum(tijm[i,j,m] + (tau_slew[(i,j,m)] + tau_exp[j])*Xijm[i,j,m] 
                        for i in range(N)[:-1] for m in range(M)) for j in range(N)[1:])
                        , 'exp_constr')
        t_min = Mod.addConstrs(((tijm[i,j,m] >= w[m]*Xijm[i,j,m]) for j in range(N) for m in range(M)
                        for i in range(N)),'t_min')
        t_max = Mod.addConstrs((tijm[i,j,m] <= w[m+1]*Xijm[i,j,m] for j in range(N) for m in range(M) 
                        for i in range(N)),'t_max')
        rise_constr = Mod.addConstrs((ti[i] >= te[i]*Yi[i] for i in range(N)),'rise_constr')
        set_constr = Mod.addConstrs((ti[i] <= tl[i]*Yi[i] for i in range(N)),'set_constr')

        #Multivisit constraints
        for targ in self.multi_visit_ind.keys():
            indeces = self.multi_visit_ind[targ]
            intra_sep = Mod.addConstrs(((gp.quicksum(tijm[indeces[i],j,m] for j in range(N)[1:] for m in range(M))
                                    >= (gp.quicksum(tijm[indeces[i-1],j,m] for j in range(N)[1:] for m in range(M))
                                    + Yi[indeces[i]]*tau_sep[targ])) 
                                    for i in range(len(indeces))[1:]),'intra_sep_constr')

        # Normalize the weighting in the objective to ensure observations are maximized
        priority_param = 50
        slew_param = 1/100

        Mod.setObjective(priority_param*gp.quicksum(Yi[j] for j in range(N)[1:-1]) 
                            -slew_param *gp.quicksum(tau_slew[(i,j,m)]*Xijm[i,j,m] for i in range(N)[1:-1] 
                                        for j in range(N)[1:-1] for m in range(M))
                            ,GRB.MAXIMIZE)
        print('Building TTP')
        Mod.params.TimeLimit = time_limit
        Mod.params.MIPGap = opt_gap
        Mod.update()

        self.gurobi_model = Mod
        

    def solve(self):
        print('Solving TTP for {} exposures with Gurobi'.format(self.N-2))
        Mod = self.gurobi_model
        Mod.optimize()

        if Mod.SolCount > 0:
            Yi = Mod.getVarByName('Yi')
            ti = Mod.getVarByName('ti')
            num_scheduled = 0
            scheduled_targets = []
            for i in range(self.N)[1:-1]:
                if np.round(Yi[i].X,0) == 1:
                    num_scheduled += 1
                    v = ti[i]
                    scheduled_targets.append((int(v.VarName[3:-1]),int(np.round(v.X))))

            print('{} of {} total exposures scheduled into script.'.format(num_scheduled,self.N-2))

            unordered_times = []

            for i in range(len(scheduled_targets)):
                unordered_times.append(int(scheduled_targets[i][1]))
            order = np.argsort(unordered_times)
            scheduled_targets = [scheduled_targets[i] for i in order]

            names = []
            orders = []
            times = []
            order = 1
            for pair in scheduled_targets:
                node_ind = pair[0]
                s = self.stars[self.node_to_star[node_ind]]
                names.append(s.name)
                times.append(Time(self.nightstarts.jd + pair[1]/(24*60),format='jd'))
                orders.append(order)
                order += 1

            return {'Order': orders, 'Name':names, 'Time':times}
            

            '''if plot_results:
                logger.info('Plotting {}'.format(current_day))

                obs_time=[]
                az_path=[]
                alt_path=[]
                names=[]
                targ_list=[]
                for pair in scheduled_targets:
                    targ=pair[0]
                    time=t[pair[1]]
                    exp=s_i[targ]
                    targ_list.append(targ)
                    names.append(all_targets_frame.loc[ind_to_id[targ],'Starname'])
                    time1=time-TimeDelta(60*exp,format='sec')
                    obs_time.append(time1.jd)
                    obs_time.append(time.jd)
                    az_path.append(get_alt_az(all_targets_frame,ind_to_id[targ],time1,keck)[1])
                    az_path.append(get_alt_az(all_targets_frame,ind_to_id[targ],time,keck)[1])
                    alt_path.append(get_alt_az(all_targets_frame,ind_to_id[targ],time1,keck)[0])
                    alt_path.append(get_alt_az(all_targets_frame,ind_to_id[targ],time,keck)[0])
                    
                time_array = []
                for i in range(len(obs_time)):
                    if i % 2 == 0:
                        time_array.append((obs_time[i] - obs_time[0])*1440)

                start = Time(obs_time[0],format='jd')
                stop = Time(obs_time[-1],format='jd')
                step = TimeDelta(2,format='sec')
                t = np.arange(start.jd,stop.jd,step.jd)
                t = Time(t,format='jd')
                        
                time_counter = Time(obs_time[0],format='jd')
                time_strings = t.isot
                observed_at_time = []
                tel_targs = []
                while time_counter.jd < obs_time[-1]:
                    ind = np.flatnonzero(time_array <= (time_counter.jd - obs_time[0])*1440)[-1]
                    req = ind_to_id[targ_list[ind]]
                    observed_at_time.append(ind)
                    ra,dec = get_ra_dec(all_targets_frame,req)
                    coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
                    target = apl.FixedTarget(name=targ, coord=coords)
                    tel_targs.append(target)
                    time_counter += TimeDelta(2,format='sec')

                AZ1 = keck.altaz(t, tel_targs, grid_times_targets=False)
                tel_az = np.round(AZ1.az.rad,2)
                tel_zen = 90 - np.round(AZ1.alt.deg,2)
                target_list = []
                for targ in targ_list:
                    ra,dec = get_ra_dec(all_targets_frame,ind_to_id[targ])
                    coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
                    target = apl.FixedTarget(name=targ, coord=coords)
                    target_list.append(target)
                    
                AZ = keck.altaz(t, target_list, grid_times_targets=True)

                total_azimuth_list = []
                total_zenith_list = []

                for i in range(len(AZ[0])):
                    total_azimuth_list.append(np.round(AZ[:,i].az.rad,2))
                    total_zenith_list.append(90-np.round(AZ[:,i].alt.deg,2))
                
                plotting.plot_path_2D(obs_time,az_path,alt_path,names,targ_list,
                                    plotpath,current_day)
                plotting.animate_telescope(time_strings,total_azimuth_list,total_zenith_list,
                                        tel_az,tel_zen,observed_at_time,plotpath)
            
            for pair in scheduled_targets:
                ordered_requests.append(ind_to_id[pair[0]])

            for i in range(R)[1:-1]:
                if np.round(yi[i].X,0) == 0:
                    extras.append(ind_to_id[i])'''

        else:
            print('No incumbent solution in time allotted, aborting solve. Try increasing time_limit parameter.')


    