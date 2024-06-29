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

    def __init__(self,start,stop,stars,observatory,runtime=300,optgap=0.01):

        # Will have observatory object
        self.observatory = observatory
        self.stars = stars
        self.nightstarts = start
        self.nightends = stop
        self.runtime = runtime
        self.optgap = optgap
        self.create_nodes()
        self.compute_tau_slew()
        self.solve()

    def create_nodes(self):
        self.dur = np.round((self.nightends-self.nightstarts).jd*24*60,0)
        t = np.arange(self.nightstarts.jd,self.nightends.jd,TimeDelta(60,format='sec').jd)
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
                tau_exp.append(s.expwithreadout)
                tau_sep.append(s.intra_night_cadence * 60) # Hours -> minutes
                AZ = self.observatory.observer.altaz(t, s.target)
                alt=AZ.alt.deg
                az=AZ.az.deg
                # if s.name == '104067':
                #     print("COORDS: ", s.target)
                #     for a in range(len(alt)):
                #         print(alt[a], az[a])
                good = np.where(self.observatory.is_up(alt,az) == 1)[0]
                if len(good > 0):
                    te.append(max(np.round(((t[good[0]]+TimeDelta(s.expwithreadout*60,format='sec'))-t[0]).jd*24*60,1),s.expwithreadout))
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

    def compute_tau_slew(self,cad=30):
        # Will want to accomodate slot input value
        M = self.observatory.nSlots
        self.M = M
        night_dur = (self.nightends.jd - self.nightstarts.jd)*24*60 # minutes

        samples_per_slot = int(max(night_dur/(M*cad),3)) # Evaluate at least every thirty minutes

        slot_bounds = Time(np.linspace(self.nightstarts.jd,self.nightends.jd,M+1,endpoint=True),format='jd')

        # Sample at least 3 times per slot
        t = []
        for m in range(M):
            times_to_eval = np.linspace(slot_bounds[m].jd,slot_bounds[m+1].jd,samples_per_slot)
            for item in times_to_eval:
                t.append(item)
        times = Time(t,format='jd')

        nodes = self.node_to_star
        N = self.N

        # Loop through targets, then evaluate coordinates at all times for all targets
        targets = []
        for i in range(N)[1:-1]:
            star = self.stars[nodes[i]]
            targets.append(star.target)
        coordinate_matrix = self.observatory.observer.altaz(
                                    times,targets,grid_times_targets=True)


        def to_wrap_frame(angle):
            if self.observatory.wrapLimitAngle:
                angle += (360-self.observatory.wrapLimitAngle)
                angle = np.array([x - 360 if x > 360 else x for x in angle])
            return angle

        def max_ang_sep(targ1,targ2,slot):

            slot_ind_start = slot*samples_per_slot
            slot_ind_end = (slot+1)*samples_per_slot-1

            # Subtract 1 since there's no dummy starting node in the coordinate matrix
            coords1 = coordinate_matrix[targ1-1,slot_ind_start:slot_ind_end+1]
            coords2 = coordinate_matrix[targ2-1,slot_ind_start:slot_ind_end+1]

            alt1 = coords1.alt.deg
            alt2 = coords2.alt.deg
            az1 = to_wrap_frame(coords1.az.deg)
            az2 = to_wrap_frame(coords2.az.deg)

            # Take the maximum of all the samples, in both dimensions
            return max(max(np.abs(alt1-alt2)),max(np.abs(az1-az2)))

        tau_slew = defaultdict(float) # Holds minutes
        for m in range(M):
            for targ1,targ2 in permutations(range(N)[1:-1],2):
                tau_slew[(targ1,targ2,m)] = np.round(max_ang_sep(targ1,targ2,m)/(60*self.observatory.slewrate),3)

        # Slot start and end times as minutes from start
        w = (slot_bounds.jd-slot_bounds[0].jd)*24*60

        self.w = w
        self.tau_slew = tau_slew


    def to_gurobi_model(self,output_flag=True):
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

        # Variables and Constraints, all with the same conventions as Handley 2024
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
                        for j in range(N)[1:]), 'visit_once')
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
                                    + Yi[indeces[i]]*tau_sep[indeces[i]]))
                                    for i in range(len(indeces))[1:]),'intra_sep_constr')

        # Normalize the weighting in the objective to ensure observations are maximized
        priority_param = 50
        slew_param = 1/100

        # Maximize number of observations, with a penalty term for long slews
        Mod.setObjective(priority_param*gp.quicksum(Yi[j] for j in range(N)[1:-1])
                            -slew_param *gp.quicksum(tau_slew[(i,j,m)]*Xijm[i,j,m] for i in range(N)[1:-1]
                                        for j in range(N)[1:-1] for m in range(M))
                            ,GRB.MAXIMIZE)
        print('Building TTP')
        Mod.params.TimeLimit = self.runtime
        Mod.params.MIPGap = self.optgap
        Mod.update()

        self.gurobi_model = Mod


    def solve(self):
        print('Solving TTP for {} exposures with Gurobi'.format(self.N-2))
        self.to_gurobi_model()
        Mod = self.gurobi_model
        Mod.optimize()

        
        if Mod.SolCount > 0:
            self.digest_gurobi()
    
            print('{} of {} total exposures scheduled into script.'.format(self.num_scheduled,self.N-2))

        else:
            print('No incumbent solution in time allotted, aborting solve. Try increasing time_limit parameter.')

    def digest_gurobi(self):
        
        N = self.N
        M = self.M
        tau_slew = self.tau_slew
        Mod = self.gurobi_model

        num_scheduled = 0
        scheduled_targets = []
        extras = []
        for i in range(self.N)[1:-1]:
            # Parsing Gurobi variables is difficult. We retrieve each variable
            # by its internal gurobi name using the node index
            Yvar = Mod.getVarByName(f'Yi[{i}]')
            if np.round(Yvar.X,0) == 1:
                num_scheduled += 1
                v = Mod.getVarByName(f'ti[{i}]')
                # Each element will have tuples as (index,time)
                scheduled_targets.append((i,v.X))
            else:# np.round(Yvar.X,0) == 0:
                s = self.stars[self.node_to_star[i]].name
                extras.append(s)
        self.extras = {'Starname':extras}
        self.num_scheduled = num_scheduled

        # Retrieve slew times for statistics
        def to_wrap_frame(angle):
            if self.observatory.wrapLimitAngle:
                angle += (360-self.observatory.wrapLimitAngle)
                if angle > 360:
                    angle -= 360
            return angle

        # Get slew estimate from the tau tensor
        est_slews = []
        for i in range(N)[1:-1]:
            for j in range(N)[1:-1]:
                for m in range(M):
                    var = Mod.getVarByName(f'Xijm[{i},{j},{m}]').X
                    if np.round(var,0) ==1:
                        est_slews.append(tau_slew[i,j,m])

        # Compute the real slew at the time the scheduler chose
        real_slews = []
        for i in range(N)[1:-1]:
            for j in range(N)[1:-1]:
                for m in range(M):
                    t = Mod.getVarByName(f'tijm[{i},{j},{m}]').X
                    if np.round(t,1) != 0:
                        minutes = np.round(t,1)
                        time_of_slew = self.nightstarts + TimeDelta(minutes*60,format='sec')
                        altaz1 = self.observatory.observer.altaz(time_of_slew,self.stars[self.node_to_star[i]].target)
                        altaz2 = self.observatory.observer.altaz(time_of_slew,self.stars[self.node_to_star[j]].target)

                        alt1 = altaz1.alt.deg
                        alt2 = altaz2.alt.deg
                        az1 = to_wrap_frame(altaz1.az.deg)
                        az2 = to_wrap_frame(altaz2.az.deg)

                        separation = max(np.abs(alt1-alt2),np.abs(az1-az2))
                        slew = separation/(60*self.observatory.slewrate)

                        real_slews.append(np.round(slew,3))

        self.estimated_slews = est_slews
        self.real_slews = real_slews
        self.solve_time = Mod.Runtime

        # Reorganize the schedule by time
        unordered_times = []

        for i in range(len(scheduled_targets)):
            unordered_times.append(int(scheduled_targets[i][1]))
        order = np.argsort(unordered_times)
        scheduled_targets = [scheduled_targets[i] for i in order]

        # Get properly ordered metadata
        starnames = []
        orders = []
        t_starts = []
        t_ends = []
        n_shots = []
        exptimes = []
        ordered_target_nodes = []
        all_times = []
        az_path = []
        alt_path = []
        order = 0
        for pair in scheduled_targets:
            node_ind = pair[0]
            ordered_target_nodes.append(node_ind)
            s = self.stars[self.node_to_star[node_ind]]
            starnames.append(s.name)
            n_shots.append(s.shots)
            exptimes.append(s.exptime)

            # Observation start + end times, as both time objects and as
            # minute markers from the start of the night
            t1 = Time(self.nightstarts.jd + pair[1]/(24*60),format='jd') - TimeDelta(60*s.expwithreadout,format='sec')
            t2 = Time(self.nightstarts.jd + pair[1]/(24*60),format='jd')
            t_starts.append(pair[1] - s.expwithreadout)
            t_ends.append(pair[1])

            # Add times, az_path, and alt_path as attributes for easy plotting
            all_times.append(t1)
            all_times.append(t2)
            coords_start_obs = self.observatory.observer.altaz(t1,s.target)
            coords_end_obs = self.observatory.observer.altaz(t2,s.target)
            az_path.append(coords_start_obs.az.deg)
            alt_path.append(coords_start_obs.alt.deg)
            az_path.append(coords_end_obs.az.deg)
            alt_path.append(coords_end_obs.alt.deg)
            orders.append(order)
            order += 1

        t_starts = np.array(t_starts)
        t_ends = np.array(t_ends)
        exptimes = np.array(exptimes)
        n_shots = np.array(n_shots)
        rise_times = (self.te - self.tau_exp)[ordered_target_nodes]
        set_times = self.tl[ordered_target_nodes]


        # Dictionary with simple output for user to view as dataframe
        self.schedule = {'Order': orders, 'Starname':starnames, 'Time':[self.nightstarts.jd + t/(24*60) for t in t_starts]}

        # Verbose dictionary with more metadata for plotting
        # Only include real nodes for this step
        self.plotly = {'Starname':starnames,
                        'First Available':rise_times,
                        'Last Available':set_times,
                        'Start Exposure':t_starts,
                        'Minutes the from Start of the Night':(t_starts + t_ends)/2,
                        'Stop Exposure':t_ends,
                        'N_shots':n_shots,
                        'Exposure Time (min)':exptimes,
                        'Total Exp Time (min)':n_shots*exptimes + (45/60)*(n_shots-1)
                        }


        self.times = all_times
        self.az_path = az_path
        self.alt_path = alt_path

