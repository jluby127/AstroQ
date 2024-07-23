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

    """Object which contains and solves a TTP.

    Args:
        start (astropy time object): Time value of the beginning of the observing interval
        stop (astropy time object): Time value of the end of the observing interval
        stars (list): List of star objects to be scheduled in the model
        observatory (object): Observatory from which observations will take place
        runtime (float): Maximum time (in seconds) to be spent solving the model in the Gurobi kernel
        optgap (float): Percentage gap between current objective and bound objective values at which 
            to terminate optimization

    """

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
        """Transform a target list into a node-based formalism.

        NOTE: All units below are in minutes

        Attributes:
            N (int): Total number of nodes in the model
            dur (float): Duration of observing interval
            te (array): Earliest possible completion times (relative to the start
                of the observing period) of all nodes
            tl (array): Latest possible completion times (relative to the start
                of the observing period) of nodes
            tau_exp (array): Exposure times of nodes
            tau_sep (array): Intranight minimum cadence of nodes
            node_to_star (dict): Map to convert any node index into the index of
                the star object which it represents
            multi_visit_ind (dict): Map which converts a star index into the list 
                of node indeces which comprise all the requested observations
        
        """
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
        """Compute the slew tensor which holds the worst case estimate of the 
        slews between any two targets during any time slot.
        
        Args:
            cad (int): Maximum spacing in minutes at which to sample the slew 
                during every slot

        Attributes:
            M (int): Number of discrete time slots for which slews will be 
                assumed constant
            w (list): Bounds (as minutes from start) of each time slot, including
                the start and end of the night
            tau_slew (array): 3D tensor which holds travel time (minutes) between
                any two nodes during a given slot
        
        """

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
            alt_sep = np.abs(alt1-alt2)
            az_sep = np.abs(az1-az2)
            
            # For telescopes without wraps, no slew greater than 180 deg can exist
            if not self.observatory.wrapLimitAngle:
                az_sep = [360 - x if x > 180 else x for x in az_sep]

            return max(max(alt_sep),max(az_sep))

        tau_slew = defaultdict(float) # Holds minutes
        for m in range(M):
            for targ1,targ2 in permutations(range(N)[1:-1],2):
                tau_slew[(targ1,targ2,m)] = np.round(max_ang_sep(targ1,targ2,m)/(60*self.observatory.slewrate),3)

        # Slot start and end times as minutes from start
        w = (slot_bounds.jd-slot_bounds[0].jd)*24*60

        self.w = w
        self.tau_slew = tau_slew


    def to_gurobi_model(self,output_flag=True):

        """Translate TTP Model object into a Gurobi Model object
        
        Args:
            output_flag (bool): Activates Gurobi console output during solving process if True.
                Console will remain silent during solve if false.
            
        Attributes:
            gurobi_model: An object in the Gurobi package which understands the TTP as a 
                mixed-integer linear program (a matrix problem). See the Gurobi docs 
                at https://www.gurobi.com/documentation/current/refman/py_model.html for
                details. Can be queried for information about the problem setup and the 
                best found solution.
        """
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
        """Convenience function to wrap the Gurobi processes in a separate place. Optimizes
            the gurobi interpretation of the TTP and calls the output routines. Catches
            situations where a solution is not found in the time limit.
        
        """
        print('Solving TTP for {} exposures with Gurobi'.format(self.N-2))
        self.to_gurobi_model()
        Mod = self.gurobi_model
        Mod.optimize()

        
        if Mod.SolCount > 0:
            self.digest_gurobi()
            self.optimization_status()

        else:
            print('No incumbent solution in time allotted, aborting solve. Try increasing time_limit parameter.')

    def digest_gurobi(self):
        """Interprets the solution to the model by parsing the values assigned to
            each variable by Gurobi. Prepares attributes necessary for schedule output
            and plots.
        
        Attributes:
            extras (dict): Holds the name of any stars which were not able to be scheduled
            num_scheduled (int): Total number of observations (>= number of targets)
                scheduled by Gurobi
            estimated_slews (list): The worst case estimates of every slew that was scheduled,
                as queried from the slew tensor
            real_slews (list): The true value of every slew scheduled, by evaluating the 
                angular separation at the exact time of the slew
            schedule (dict): Holds the order of the successfully scheduled targets, their 
                name, and what time they should be observed
            plotly (dict): A more verbose description of the solution which holds additional
                information for the interactive solution plot
            times (list): The times at which the state of the telescope changes (i.e., an
                observation begins or ends) as astropy time objects
            az_path (list): The azimuth direction (degrees) of the telescope at every state 
                change
            alt_path (list): The altitude direction (degrees) of the telescope at every state 
                change
        """
        
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

                        alt_sep = np.abs(alt1-alt2)
                        az_sep = np.abs(az1-az2)
                        
                        # For telescopes without wraps, no slew greater than 180 deg can exist
                        if not self.observatory.wrapLimitAngle:
                            if az_sep >= 180:
                                az_sep = 360 - az_sep

                        separation = max(alt_sep,az_sep)
                        slew = separation/(60*self.observatory.slewrate)

                        real_slews.append(np.round(slew,3))

        self.estimated_slews = est_slews
        self.real_slews = real_slews

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

    def optimization_status(self):
        """Parse the Objective function and output some stats to the console.

        NOTE: Modifying the priorities of different targets will make it
            impossible to interpret the bounds of the objective.

        Attributes:
            obs_bound (int): The maximum number of observations that was deemed
                theoretically possible in the time limit
            slew_bound (float): The lowest slew time to observe all the targets
                deemed theoretically possible in the time limit
            time_idle (float): The total time (minutes) within the observing 
                interval which was not spent slewing or exposing
            time_slewing (float): The total time (minutes) spent slewing
            time_exposing (float): The total time (minutes) spent exposing
            solve_time (float): The runtime (seconds) spent within the Gurobi solving
                kernel. This will be the specified time limit unless the chosen 
                optimality gap was achieved in less than that time limit.
        """

        # The optimization weights set the scale
        slew_param = 1/100
        priority_param = 50

        # Unpack the objective value and bounds
        Mod = self.gurobi_model
        lower = 1 + (Mod.ObjVal // 1)
        slewtime = (lower - Mod.ObjVal)*1/slew_param
        
        upper = 1 + (Mod.ObjBound // 1)
        max_obs = int(upper/priority_param)
        bound_slew = (upper - Mod.ObjBound)*1/slew_param

        time_exposing = 0
        for i in range(self.N)[1:-1]:
            Y = Mod.getVarByName(f'Yi[{i}]').x
            time_exposing += Y*self.tau_exp[i]
        time_idle = self.dur - time_exposing - slewtime

        # Save aspects of the best solution and the tightest bound to get theoretical
        # limits from the solve.
        self.obs_bound = max_obs
        self.slew_bound = bound_slew
        self.time_idle = time_idle
        self.time_slewing = slewtime
        self.time_exposing = time_exposing
        self.solve_time = Mod.Runtime

        print('\n')
        print('------------------------------------')
        print(f'    Model ran for {self.solve_time:.2f} seconds')
        print('------------------------------------')
        print(f'     Observations Requested: {self.N-2}')
        print(f'     Observations Scheduled: {self.num_scheduled}')
        print(f' Maximum Observations Bound: {self.obs_bound}')
        print('------------------------------------')
        print(f'   Observing Duration (min): {self.dur:.2f}')
        print(f'  Time Spent Exposing (min): {self.time_exposing:.2f}')
        print(f'      Time Spent Idle (min): {self.time_idle:.2f}')
        print(f'   Time Spent Slewing (min): {self.time_slewing:.2f}')
        print(f'   Minimum Slew Bound (min): {self.slew_bound:.2f}')
        print('------------------------------------')   