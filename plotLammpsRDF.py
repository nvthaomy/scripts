import matplotlib ,os
import matplotlib.pyplot as plt  
import argparse
parser = argparse.ArgumentParser() 
parser.add_argument("-i",required=True,nargs="+", help="rdf.txt files")
parser.add_argument("-p", required = True, help = "pair name") 
parser.add_argument("-l",required=True,nargs="+", help="legends")
parser.add_argument("-t", help="chart title")
args = parser.parse_args() 
files = args.i
pair = args.p
linestyle=('-','--')
legend=args.l
plt.figure()
for k,i in enumerate(files):
	name=i.split('.txt')[0]
	IsLammps = False
	r=[]
	g_r=[]
	f=open(i,'r')
	line=f.readline()
	while len(line):
		if "Time-averaged" in line:
			IsLammps = True
		if not line.startswith('#'):
			if IsLammps:
				if len(line.split())>2:
					r.append(float(line.split()[1]))
					g_r.append(float(line.split()[2]))
			else:
				r.append(10*float(line.split()[0]))
				g_r.append(float(line.split()[1]))
		line=f.readline()
	plt.plot(r,g_r,label=legend[k],linestyle=linestyle[k],linewidth=2)

plt.xlabel('r') 
plt.ylabel('$g_r$(r)')
plt.legend(loc='best')
plt.xlim(0,10) 
plt.ylim(0)
if args.t:
	plt.title(args.t,loc='center')
	plt.savefig('rdf-{}_{}.png'.format(pair,'_'.join(args.t.split())))
else:
	plt.savefig('rdf-{}.png'.format(pair))
plt.show() 


