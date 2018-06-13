import numpy as np
import random
import os
import subprocess
import script_utils

def modify_top(top='compound.top',
               include="#include \"/raid6/homes/ahy3nz/Programs/McCabeGroup/atomistic/forcefield.itp\""):
    """ Modify topology file include statement

    Parameters
    ---------
    top : str
        filename of gmx topology
    include : str
        new include statement

    Notes'
    ------
    Top file is overwritten
    Assumes the include statement is in the first line

    Useful includes to remember:
    rahman: /raid6/homes/ahy3nz/Programs/McCabeGroup/atomistic/forcefield.itp
    edison: /global/homes/a/ahy3nz/Programs/McCabeGroup/atomistic/forcefield.itp
    accre: /home/yangah/Programs/McCabeGroup/atomistic/forcefield.itp

    """
    toplines = open(top,'r').readlines()
    toplines[0] = include + "\n"
    with open(top,'w') as f:
        for line in toplines:
            f.write(line)

def submit_job(script, jobid, n_nodes, i):
    """ Submit job to cluster 

    Parameters
    ---------
    script : str
        name of submission script
    n_nodes : int
        number of nodes we are using, equivalently the number
        of jobs running at any one time with others on hold
    jobid : list 
        Constains the last submited job for each node
    i : int
        Index/iteration number of the job we are submitting

    Returns
    -------
    jobid : list
        updated list corresponding to jobids last submitted


    """
    if 'pbs' in script:
        if not jobid[i % n_nodes]:
            p = subprocess.Popen('qsub {}'.format(script), shell=True, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE)
            p.wait()
        elif jobid[i % n_nodes]:
            p = subprocess.Popen('qsub -W depend=afterany:{0} {1}'.format(\
                                 jobid[i % n_nodes], script), shell=True, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE)
            p.wait()
        outs, errs = p.communicate()
        jobid[i % n_nodes] = outs.decode().strip()
    elif "sbatch" in script:
        if not jobid[i % n_nodes]:
            p = subprocess.Popen('sbatch {}'.format(script), shell=True, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE)
            p.wait()
        elif jobid[i % n_nodes]:
            p = subprocess.Popen('sbatch -d afterany:{0} {1}'.format(\
                                 jobid[i % n_nodes], script), shell=True, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE)
            p.wait()
        outs, errs = p.communicate()
        jobid[i % n_nodes] = outs.decode().strip().split()[-1]
    return jobid

###############################
### Equilibration functions ###
###############################

def write_eq_lines(gro='compound.gro', top='compound.top'):
    """ Write EM, NVT, NPT lines for equilibration """
    lines = """
gmx grompp -f em.mdp -c {gro} -p {top} -o em -maxwarn 2 &> em_grompp.log
gmx mdrun -deffnm em 

gmx grompp -f nvt.mdp -c em.gro -p {top} -o nvt -maxwarn 2 &> nvt_grompp.log
gmx mdrun -deffnm nvt

gmx grompp -f npt_500ps.mdp -c nvt.gro -p {top} -o npt_500ps -t nvt.cpt -maxwarn 2 &> npt_500ps_grompp.log
gmx mdrun -deffnm npt_500ps -ntomp 1""".format(**locals())
    return lines

###########################
## Production functions ###
###########################

def write_production_lines(filename='npt'):
    """ Write NPT production """
    lines = "gmx mdrun -ntomp 8 -ntmpi 2 -gpu_id 01 -deffnm {filename} -cpi {filename}_prev.cpt -append".format( **locals())
    return lines


######################
### RWMD functions ###
######################

def write_rwmd_files(components, gro='npt_500ps.gro', top='compound.top', 
                    cooling_rate=1000, t_max=385):
    """ Write files relevant for RWMD

    Parameters
    ---------
    components : []
        List of each component (str) in system 
    gro : str
    top : str
    cooling_rate : float, default 1000
        RWMD ceiling drops by 1 K every `cooling_rate` ps

    Notes
    -----
    Will require a recompiled verion of gromacs, modified for longer string buffers
    """
    # First 25 ns is just the maximum range of heating
    tc_groups = ['non-water', 'water']
    #tc_groups = [i for i in components]
    #tc_groups.append('water')
    t_pairs = [[305, t_max] for _ in tc_groups]
    #t_pairs = [[305,385],[305,385]]
    total_cooling_time = (t_pairs[0][1] - t_pairs[0][0])*cooling_rate 
    with open("heating_phase.mdp", 'w')as f:
        _write_body(f, ref_temp=305, n_steps=15000000, tc_groups=tc_groups)
        times, all_temps, t_pairs = _generate_heating_phase(f, t_pairs=t_pairs, 
                interval=5, time_anneal=30000, 
                dtemp=10, water_thermostat_style='plateau')

        _write_annealing_lines(f, times, all_temps)

    
    # Cooling phase is a series of 25ns simulations
    for i, t_start in enumerate(np.arange(30000, 30000 + total_cooling_time, 20000)):
        #print("Writing cooling_phase{}.mdp".format(i))
        #print("Temperature Pairs: {}".format(t_pairs))
        with open("cooling_phase{}.mdp".format(i),'w') as f:
            _write_body(f, ref_temp=305, n_steps=10000000, t_init=t_start, 
                    tc_groups=tc_groups)
            times, all_temps, t_pairs = _generate_cooling_phase(t_pairs=t_pairs,
                duration=20000, interval=5, cooling_rate=cooling_rate, 
                time_start=t_start, dtemp=10,
                water_thermostat_style='plateau')
            _write_annealing_lines(f, times, all_temps)
    return _write_rahman_rwmd(gro=gro, top=top, n_cooling=i)

def _write_rahman_rwmd(gro='npt_500ps.gro', top='compound.top', n_cooling=4):
    """ Write submission lines

    Parameters
    ---------
    gro : str
    top : str
    n_cooling : int
        last # of cooling cycle
        n_cooling 4 means the last index is 4 (so there are actually 5 cooling phases)
    
    Returns
    lines : sttr
    
    """
    lines = """
mygmx grompp -f heating_phase.mdp -c {gro} -p {top} -o heating_phase &> heating_phase.out
gmx mdrun -ntomp 8 -ntmpi 2 -gpu_id 01 -deffnm heating_phase -cpi heating_phase.cpt -append

mygmx grompp -f cooling_phase0.mdp -c heating_phase.gro -p {top} -t heating_phase.cpt -o cooling_phase0
gmx mdrun -ntomp 8 -ntmpi 2 -gpu_id 01 -deffnm cooling_phase0 -cpi cooling_phase0.cpt -append &> cooling_phase0.out


for ((i=1; i<={n_cooling} ; i++))
do
    mygmx grompp -f cooling_phase${{i}}.mdp -c cooling_phase$((${{i}}-1)).gro -p {top} -t cooling_phase$((${{i}}-1)).cpt -o cooling_phase${{i}}
    gmx mdrun -ntomp 8 -ntmpi 2 -gpu_id 01 -deffnm cooling_phase${{i}} -cpi cooling_phase${{i}}.cpt -append &> cooling_phase${{i}}.out
done

gmx grompp -f npt.mdp -c cooling_phase$((${{i}}-1)).gro -p {top} -o npt > npt_grompp.out
""".format(**locals())
    return lines



def _write_body(f, ref_temp=305, n_steps=50000000,t_init=0,tc_groups=None):
    tc_grps_string = " ".join([thing for thing in tc_groups])
    tau_t_string = " ".join(["1.0" for _ in range(len(tc_groups))])
    ref_t_string = " ".join(["{:5.0f}".format(ref_temp) for _ in range(len(tc_groups))])
    annealing_string = " ".join(["single" for _ in range(len(tc_groups))])
    f.write(""" title                       = RWMD
; Run parameters
integrator                  = md
nsteps                      = {n_steps}     ; 2fs/step * 5e7 steps = 5e7 fs = 50 ns
dt                          = 0.002
tinit                       = {t_init}

; Output control
nstxout                     = 0             ; Don't save coordinates 
nstvout                     = 0             ; Don't save velocities
nstenergy                   = 10000
nstlog                      = 10000
nstxtcout                   = 10000

;bond parameters
continuation                = yes
constraint_algorithm        = lincs
constraints                 = h-bonds
;constraints                 = none
lincs_iter                  = 1
lincs_order                 = 4

; Neighbor searching
cutoff-scheme               = Verlet
ns_type                     = grid
nstlist                     = 10
rcoulomb                    = 1.4
rvdw                        = 1.4

;Electrostatics
coulombtype                 = PME
pme_order                   = 4
fourierspacing              = 0.16

; Temperature coupling
tcoupl                      = nose-hoover
tc-grps                     = {tc_grps_string}    
tau_t                       = {tau_t_string}
ref_t                       = {ref_t_string}

;Pressure coupling
pcoupl                      = Parrinello-Rahman
pcoupltype                  = semiisotropic
tau_p                       = 5.0           ; ps
ref_p                       = 1.0 1.0          ; bar   
compressibility             = 4.5e-5 4.5e-5
refcoord_scaling            = com

;PBC
pbc                         = xyz

;Dispersion correction
DispCorr                    = EnerPres

;Velocity generation
gen_vel                     = no

;Simulated annealing
annealing                   = {annealing_string}
""".format(**locals()))

def _generate_heating_phase(f,t_pairs=None,
        interval=5, time_anneal=25000, final_time=50000, dtemp=20,
        water_thermostat_style=None):
    """Generate annealing points
    f : File to write to 
    t_pairs : list of pairs,
        In each pair, (low temp, high temp). Convention
        is that the first ordered pair is the DSPC and the last 
        ordered pair is the water 
        """
    n_points = int(time_anneal/interval)

    #all_temps = [[] for _ in range(len(temp_pairs))]
    all_temps = np.zeros((len(t_pairs),n_points))

    # Initialize all temperatuers to their low temperature
    all_temps[:,0] = [temp[0] for temp in t_pairs]

    times = [0]

    # Generate scaling factors based on temperatre pairs
    if water_thermostat_style == 'proportional':
        scaling_factors = np.zeros(len(t_pairs))
        for i, t_pair in enumerate(t_pairs): 
            scaling_factors[i] = (t_pairs[i][1] - t_pairs[i][0]) / \
                (t_pairs[0][1] - t_pairs[0][0])
    

    possible_temps = np.arange(t_pairs[0][0], t_pairs[0][1]+1, dtemp)
    # Construct the histogram
    histogram = {val: 0 for val in possible_temps}
    histogram[all_temps[0,0]] += 1
    for i in np.arange(1,n_points):
        # new stuff here
        time = i*interval
        times.append(time)

        # Randomly pick a temperature
        new_temps = np.zeros(len(t_pairs))
        #candidate_temp = random.choice(possible_temps)
        adjacent_temps = [all_temps[0, i-1] - dtemp, all_temps[0, i+1] + dtemp]
        candidate_temp = random.choice(adjacent_temps)

        # Compare it to a criteria
        if histogram[all_temps[0,i-1]] >= histogram[candidate_temp]:
            new_temp = candidate_temp
            # Adjust temperature step for each tc group
            for j, tc_group in enumerate(t_pairs):
                if water_thermostat_style == 'proportional':
                    primary_temp_change = new_temp - t_pairs[0][0]
                    other_temp_change = scaling_factors[j] * primary_temp_change 
                    all_temps[j, i] = tc_group[0] + other_temp_change
            
                if water_thermostat_style == 'plateau':
                    if candidate_temp > tc_group[1]:
                        all_temps[j,i] = tc_group[1]
                    else:
                        all_temps[j,i] = candidate_temp
                    
        else:
            all_temps[:,i] = all_temps[:, i-1]

        histogram[all_temps[0,i]] += 1
        

    ### Saving stuff for reference
    with open('rwmd_histogram.dat', 'w') as thing:
        for k, v in histogram.items():
            thing.write("{}\t{}\n".format(k,v))
    #np.savetxt('timeseries.dat', np.column_stack((times, temps)))
    #### 

    #return times, temps, w_temps, (nonw_temp_low, nonw_temp_high), (w_temp_low, w_temp_high)
    return times, all_temps, t_pairs

    

def _generate_cooling_phase(time_start=25000, duration=50000, interval=5, 
        cooling_rate=500,
        t_pairs=None, dtemp=10,
        water_thermostat_style=None):
    """
    Change temperatures every 5 ps
    Reduce ceiling 2 K every ns
    Reduce ceiling 1 K every 500ps
    time_start : int
        Initial time of simulation (ps)
    duration : int
        Length of cooling phase (ps)
    interval :
        duration of particular temperature before attempting jum (ps)
    cooling_rate : float
        Every `cooling_rate` ps, drop the RWMD temperature ceiling by 1 K
    t_pairs : list of 2-tuples
        Elements of t_pairs are 2-tuples, where each element is a particular
        temperatuer group, and the elements of the tuple are (temp_low and temp_max)
    dtemp : int
        Potential temperature jumps
    water_thermostat_style : str
        Specifies how to thermstat water 
        'plateau' for water to plateau if thermostat exceeds a temperature
        'proportional' for water temp changes to be proportional to lipid temp

    """

    # This is a check to make sure temperatures aren't 
    # excessively higher than the primary temps
    n_points = int(np.floor(duration/interval))
    for i, t_pair in enumerate(t_pairs):
        if t_pair[1] > t_pairs[0][1]:
            t_pairs[i][1] = t_pairs[0][1]
        if t_pair[0] > t_pairs[0][0]:
            t_pairs[i][0] = t_pairs[0][0]
    
    if water_thermostat_style == 'proportional':
        scaling_factors = np.zeros(len(t_pairs))
        for i, t_pair in enumerate(t_pairs): 
            scaling_factors[i] = (t_pairs[i][1] - t_pairs[i][0]) / \
                (t_pairs[0][1] - t_pairs[0][0])

    times = [time_start]

    all_temps = np.zeros((len(t_pairs),n_points+1))

    # Initialize all temperatuers to their high temperature
    all_temps[:,0] = [temp[1] for temp in t_pairs]


    for i,time in enumerate(np.arange(time_start+interval, time_start+duration+interval, interval)):
        times.append(time)
        # If we've hit a cooling step, reduce the temperature range
        # Ensure that there is still a range of temperatures to choose from
        if time % cooling_rate == 0 and t_pairs[0][1] > t_pairs[0][0] + 1:
            t_pairs[0][1] -= 1

        #new_temp = np.random.randint(t_pairs[0][0], t_pairs[0][1])
        new_temp = np.random.choice(all_temps[0,i] - dtemp, all_temps[0,i] + dtemp)

        for j, tc_group in enumerate(t_pairs):
            if water_thermostat_style == 'proportional':
                primary_temp_change = new_temp - t_pairs[0][0]
                other_temp_change = scaling_factors[j] * primary_temp_change 
                all_temps[j, i+1] = tc_group[0] + other_temp_change
        
            if water_thermostat_style == 'plateau':
                if new_temp > tc_group[1]:
                    all_temps[j,i+1] = tc_group[1]
                else:
                    all_temps[j,i+1] = new_temp 

        
    return times, all_temps, t_pairs 


def _write_annealing_lines(f, times, all_temps):
    if len(times) == np.shape(all_temps)[1]:
        temp_string = ""
        small_time_string = "".join(["{} ".format(thing) for thing in times])
        for i in range(all_temps.shape[0]):
            temp_string += "".join(["{:3.0f} ".format(thing) for thing in all_temps[i,:]])


        n_points = "".join(["{} ".format(len(times)) for _ in range(all_temps.shape[0])])
        time_string = " ".join([small_time_string for _ in range(all_temps.shape[0])])


        f.write("annealing-npoints\t={0}\n".format(n_points))
        f.write("annealing-time\t={0}\n".format(time_string))
        f.write("annealing-temp\t={0}\n".format(temp_string))
    else:
        sys.exit("Error, annealing lines do not have matching dimensions")


##########################
### Analysis functions ###
##########################

def analysis_routine(trajfile, grofile, pdbfile):
    traj = mdtraj.load(trajfile, top=grofile)
    traj_pdb = mdtraj.load(trajfile, top=pdbfile)
    topol = traj.topology

    # Compute system information
    lipid_tails, headgroup_dict = bilayer_analysis_functions.identify_groups(traj, 
            forcefield='charmm36')
    n_lipid = len([res for res in traj.topology.residues if not res.is_water])
    n_lipid_tails = len(lipid_tails.keys())
    n_tails_per_lipid = n_lipid_tails/n_lipid



    # Vectorized Calculations start here
    apl_avg, apl_std, apl_list = bilayer_analysis_functions.calc_APL(traj,n_lipid, blocked=True)
    np.savetxt('apl.dat', apl_list)

    angle_avg, angle_std, angle_list = bilayer_analysis_functions.calc_tilt_angle(traj, topol, lipid_tails, blocked=True)
    np.savetxt('angle.dat', angle_list)

    apt_avg, apt_std, apt_list = bilayer_analysis_functions.calc_APT(traj, apl_list, angle_list, n_tails_per_lipid, 
            blocked=True)
    np.savetxt('apt.dat', apt_list)

    s2_ave, s2_std, s2_list = bilayer_analysis_functions.calc_nematic_order(traj, blocked=True)
    np.savetxt('s2.dat', s2_list)

    headgroup_distance_dict = bilayer_analysis_functions.compute_headgroup_distances(traj, topol, headgroup_dict, blocked=True)
    Hpp_ave, Hpp_std, Hpp_list = bilayer_analysis_functions.calc_bilayer_height(traj, headgroup_distance_dict, blocked=True, anchor='DSPC')
    np.savetxt('height.dat', Hpp_list)

    offset_dict = bilayer_analysis_functions.calc_offsets(traj, headgroup_distance_dict, blocked=True, anchor='DSPC')

    d_a, d_t, d_b, bins, interdig_list,interdig_avg, interdig_std = \
        bilayer_analysis_functions.calc_density_profile(traj, topol, 
                                                        blocked=True)
    np.savetxt('idig.dat', interdig_list)
    ##print('Calculating hydrogen bonds...')
    ##hbond_matrix_avg, hbond_matrix_std, hbond_matrix_list, labelmap = bilayer_analysis_functions.calc_hbonds(traj, traj_pdb, topol, lipid_dict, headgroup_dict)
    #
    # Printing properties
    outpdf = PdfPages(('bilayeranalysis.pdf'))
    datafile = OrderedDict()
    datafile['trajectory'] = trajfile
    datafile['structure'] = grofile
    datafile['n_frames'] = traj.n_frames
    datafile['lipids'] = n_lipid
    datafile['tails'] = n_lipid_tails
    datafile['APL'] = OrderedDict()
    datafile['APL']['unit'] = str(apl_avg.unit)
    datafile['APL']['mean'] = float(apl_avg._value)
    datafile['APL']['std'] = float(apl_std._value)
    datafile['APT'] = OrderedDict()
    datafile['APT']['unit'] = str(apt_avg.unit)
    datafile['APT']['mean'] = float(apt_avg._value)
    datafile['APT']['std'] = float(apt_std._value)
    datafile['Bilayer Height'] = OrderedDict()
    datafile['Bilayer Height']['unit'] = str(Hpp_ave.unit)
    datafile['Bilayer Height']['mean'] = float(Hpp_ave._value)
    datafile['Bilayer Height']['std'] = float(Hpp_std._value)
    datafile['Tilt Angle'] = OrderedDict()
    datafile['Tilt Angle']['unit'] = str(angle_avg.unit)
    datafile['Tilt Angle']['Bilayer'] = OrderedDict()
    datafile['Tilt Angle']['Bilayer']['mean'] = float(angle_avg._value)
    datafile['Tilt Angle']['Bilayer']['std'] = float(angle_std._value)
    datafile['S2'] = OrderedDict()
    datafile['S2']['mean'] = s2_ave
    datafile['S2']['std'] = s2_std
    datafile['Interdigitation'] = OrderedDict()
    datafile['Interdigitation']['unit'] = str(interdig_avg.unit)
    datafile['Interdigitation']['mean'] = float(interdig_avg._value)
    datafile['Interdigitation']['std'] = float(interdig_std._value)

    datafile['Offset'] = OrderedDict()
    for key in offset_dict.keys():
        datafile['Offset']['unit'] = str(offset_dict[key][0].unit)
        datafile['Offset'][key] = OrderedDict()
        datafile['Offset'][key]['mean'] = float(offset_dict[key][0]._value )
        datafile['Offset'][key]['std'] = float(offset_dict[key][1]._value )
        #datafile['Offset (A)'][key] = [str(offset_dict[key][0]), str(offset_dict[key][1])]

    datafile['Tilt Angle']['Leaflet 1'] = OrderedDict()
    datafile['Tilt Angle']['Leaflet 1']['mean'] = float(np.mean(angle_list[:, 
                                            0 :int(np.floor(n_lipid_tails/2))])._value)
    datafile['Tilt Angle']['Leaflet 1']['std'] = float(np.std(angle_list[:, 
                                                0 :int(np.floor(n_lipid_tails/2))])._value)

    datafile['Tilt Angle']['Leaflet 2'] = OrderedDict()
    datafile['Tilt Angle']['Leaflet 2']['mean'] = float(np.mean(angle_list[:, 
                                                    int(np.floor(n_lipid_tails/2)):])._value)
    datafile['Tilt Angle']['Leaflet 2']['std'] = float(np.std(angle_list[:, 
                                                int(np.floor(n_lipid_tails/2)):])._value)
    #for row_label in labelmap.keys():
    #    for col_label in labelmap.keys():
    #        row_index = labelmap[row_label]
    #        col_index = labelmap[col_label]
    #        hbond_avg = hbond_matrix_avg[row_index, col_index]
    #        hbond_std = hbond_matrix_std[row_index, col_index]
    #        outfile.write('{:<20s}: {} ({})\n'.format(str(row_label+"-"+ col_label), hbond_avg, hbond_std))


    # Plotting

    fig1 = plt.figure(1)
    plt.subplot(3,2,1)
    plt.plot(apl_list)
    plt.title('APL')

    plt.subplot(3,2,2)
    plt.plot(np.mean(angle_list, axis=1))
    plt.title('Tilt Angle ($^o$)')

    plt.subplot(3,2,3)
    plt.plot(np.mean(apt_list,axis=1))
    plt.title('APT')

    plt.subplot(3,2,4)
    plt.plot(Hpp_list)
    plt.title('H$_{PP}$')

    plt.subplot(3,2,5)
    plt.plot(s2_list)
    plt.title('S2')

    plt.subplot(3,2,6)
    plt.plot(interdig_list)
    plt.title('Interdigitation (A)')

    plt.tight_layout()
    outpdf.savefig(fig1)
    plt.close()

    density_profile_top_avg = np.mean(d_t, axis = 0)
    density_profile_bot_avg = np.mean(d_b, axis = 0)
    density_profile_avg  = np.mean(d_a, axis=0)
    #
    #
    fig2 = plt.figure(2)
    plt.subplot(2,1,1)
    plt.plot(bins,density_profile_avg)
    plt.xlabel('Depth (nm)')
    plt.title('Density Profile (kg m$^{-3}$)')


    plt.subplot(2,1,2)

    #plt.plot(bins,density_profile_bot_avg)
    #plt.plot(bins,density_profile_top_avg)

    plt.hist(np.mean(angle_list[:, 0 : int(np.floor(n_lipid_tails/2))], axis=0)._value, 
                    bins=50,  
                    alpha=0.5, facecolor='blue', normed=True)
    plt.hist(np.mean(angle_list[:, int(np.floor(n_lipid_tails/2)) : ], 
                    axis=0)._value, bins=50,  
                    alpha=0.5, facecolor='red', normed = True)
    plt.title('Angle Distribution by Leaflet')
    plt.xlabel('Angle ($^o$)')

    plt.tight_layout()
    outpdf.savefig(fig2)
    plt.close()
    outpdf.close()
    with open('data.txt', 'w') as f:
        json.dump(datafile, f, indent=2)
    

