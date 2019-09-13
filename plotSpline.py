import spline
import numpy as np
""" Plot pair potential with spline, 1 gaussian and 2 gaussians
Output: txt files of pair potential values"""

G1 = False
G2 = False
Spline = True
rcut = float(11)
knots = [2.7835e+02 , 3.3541e+00 , -5.8015e-01, 1.6469e-01 , -1.1965e-01, 5.2720e-02 , -2.3451e-02, 2.6243e-03 ]

#1 Gaussian parameters:
B = 9.2563e+01
K = 7.3395e-01

#2 Gaussian parameters:
B1 = 6.2285e+01
K1 = 6.1858e-01
B2 = -8.5466e-01
K2 = 1.4607e-01 


knots = np.array(knots)
myspline = spline.Spline(rcut,knots)
nknots = len(knots)
dr = 11/nknots
rs = []
for i in range(nknots):
	rs.append(dr*i)
rs = np.array(rs)
#myspline.fitCoeff(rs,knots)
r = np.linspace(0,rcut,1000)
u = []
u_1G = []
u_2G = []
if Spline:
	fs=open("spline.txt",'w')
	fs.write("# r Pair")  
if G1:
	f1 = open("1gauss.txt",'w')
	f1.write("# r Pair")
if G2:
	f2 = open("2gauss.txt",'w')
	f2.write("# r Pair")
for i in r:
	if G1:
		u_1Gi = B*np.exp(-K*i**2)
		u_1G.append(u_1Gi)
		f1.write('\n{}   {}'.format(i,u_1Gi))
	if G2:
		u_2Gi = B1*np.exp(-K1*i**2) + B2*np.exp(-K2*i**2)
		u_2G.append(u_2Gi)
		f2.write('\n{}   {}'.format(i,u_2Gi))
	if Spline:
		u_i = myspline.Val(i)
		u.append(u_i)
		fs.write('\n{}   {}'.format(i,u_i))

import matplotlib.pyplot as plt
if Spline:
	plt.plot(r,u,label="spline")
if G1:
	plt.plot(r,u_1G,label="1 Gaussian")
if G2:
	plt.plot(r,u_2G,label="2 Gaussian")
plt.ylabel('Pair Potential')
plt.xlabel('r')
plt.xlim(0)
plt.legend(loc='best')
plt.savefig('spline.png')

plt.show()


