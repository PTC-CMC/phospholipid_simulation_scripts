import os
import json
import mdtraj
import numpy as np
import bilayer_analysis_functions
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import plot_ay
plot_ay.setDefaults()

curr_dir = os.getcwd()

index = json.load(open('newindex.txt', 'r'))
frame_step = 10 # one frame is 10 ps
for folder in index.keys():
    print(folder)
    os.chdir(os.path.join(curr_dir, folder))
    traj = mdtraj.load('npt_80-100ns.xtc', top='npt.gro')
    (angles_trajectory, angles_err) = bilayer_analysis_functions.compute_rotations(traj)
    blocked_angles, blocked_err =  bilayer_analysis_functions.block_avg(traj, 
            angles_trajectory) 
    times = frame_step * np.arange(len(angles_trajectory))
    fig, ax = plt.subplots(1,1)
    ax.set_xlabel("Time [ps]")
    ax.set_ylabel("Rotation [rad]")
    ax.set_ylim(-0.1, 0.3)
    l, = ax.plot(times, angles_trajectory)
    ax.fill_between(times, angles_trajectory - angles_err, 
            angles_trajectory + angles_err, alpha=0.4, color=l.get_color())
    plot_ay.tidyUp(fig, ax, gridArgs={}, tightLayoutArgs={})
    fig.savefig('rotations_over_time.png')
    plt.close(fig)
    np.savetxt('rotations_over_time.dat', 
            np.column_stack((times, angles_trajectory,angles_err)))
    with open('rotations.json', 'w') as f:
        stats = {}
        stats['mean'] = np.mean(blocked_angles)
        stats['std'] = np.std(blocked_angles)
        json.dump(stats, f, indent=2)

