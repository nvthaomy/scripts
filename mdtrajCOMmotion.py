import os, sys, re
import numpy as np
import mdtraj as md
       
if __name__ == "__main__":
    import argparse as ap
    import matplotlib.pyplot as plt
    import matplotlib
    showPlots = True
    try:
      os.environ["DISPLAY"] #Detects if display is available
    except KeyError:
      showPlots = False
      matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window
    matplotlib.rc('font', size=7)
    matplotlib.rc('axes', titlesize=7)

    parser = ap.ArgumentParser(description="get COM motion of system")
    parser.add_argument('traj', type=str, help = "AA trajectory")
    parser.add_argument('top', type=str, help = "AA topology")
    parser.add_argument('dt', type=float, help = "dt in simulation in picoseconds")
    parser.add_argument('trajStride', type=float, help = "stride of traj in simulation")
    parser.add_argument('-stride', type=int, help = "stide", default = 1)
    args = parser.parse_args()

    traj = args.traj 
    top = args.top
    dt = args.dt
    trajStride = args.trajStride
    stride = args.stride

    print('...loading trajectory...')
    traj = md.load(traj,top = top, stride=stride)
    top = traj.topology
    print('...done loading...')
    COMs = md.compute_center_of_mass(traj)
    ds = (COMs[0:-1] - COMs[1:])**2
    ds = np.sqrt(np.sum(ds,axis=1))

    dt = dt*trajStride*stride
    vs = ds/dt #nm/picosecond
    v = np.mean(vs)
    v1s = vs * 1000
    v1 = v * 1000 #m/s
    print('Average COM motion is {:3.3e} nm/ps or {:3.3f} m/s'.format(v,v1)) 

    #get water motion:
    atomSelect = top.select('name EPW')
#    print(len(atomSelect),' molecules')
#    print(traj.n_frames, ' frames')
    trajHOH = traj.atom_slice(atomSelect)
    rHOHs = trajHOH.xyz
    dHOHs = (rHOHs[0:-1] - rHOHs[1:])**2
    dHOHs = np.sqrt(np.sum(dHOHs,axis=2))
#    print(dHOHs.shape[0],' frames ',dHOHs.shape[1],' molecules')
    vHOHs = dHOHs/dt
    vHOHs = np.mean(vHOHs,axis=1) #avg over molecule
    vHOH = np.mean(vHOHs,axis=0) #avg over frames
    print('Average velocity of water is {:3.3e} nm/ps'.format(vHOH))

    fig,axs = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
    axs.plot(range(len(v1s)), v1s, marker=None,ls='-',lw=1, c = 'g')
    plt.xlabel('frame')
    plt.ylabel('$v_{COM} (m/s)$')
    title = 'COM motion'
    plt.title(title, loc = 'center')
    plt.savefig('_'.join(re.split(' |=|,',title))+'.png',dpi=500,transparent=True,bbox_inches="tight")
    plt.show()
