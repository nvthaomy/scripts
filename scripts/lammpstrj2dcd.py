import mdtraj,argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t", nargs='+',required=True,help=".lammpstrj files")
parser.add_argument("-top", nargs='+',required=True,help="topology (pdb) files")
args = parser.parse_args()

top = args.top
wtraj = args.t
for i in range(len(wtraj)):
	t=mdtraj.load(wtraj[i],top = top[i])
	out = wtraj[i].split('.lammpstrj')[0]+'.dcd'
	t.save(out)
