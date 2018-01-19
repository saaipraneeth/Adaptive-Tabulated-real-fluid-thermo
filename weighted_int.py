'''
########################################################################################################
#  A Program to read the rho, E table and find the interpolated thermodynamic properties for a given point
#  Written: August 2016
#  Author: Sai
########################################################################################################
'''
import numpy as np
import pickle
#import realFluidTableGeneration as rftg
import math
import pdb
import matplotlib.pyplot as plt

#rho_target = input("Enter the Density target in kg/m3 ")
#e_target = input("Enter the energy target in KJ/Kg  ")

# T_target = input("Enter the Temperature target in K ")
# P_target = input("Enter the Pressure target in Pa  ")

f=open("interp.pickle", "rb" )
e_interp,T_interp, P_interp, H_interp, rho_act_int, cv_interp, cp_interp, a_interp = pickle.load(f)

rinterp = rho_act_int
einterp = e_interp

interp = {'T': T_interp, 'P': P_interp, 'H': H_interp, 'Rho': rho_act_int, 'E': e_interp, 'CV': cv_interp, 'CP': cp_interp, 'a': a_interp}
# interp[0] = rftg.T_interp
# interp[1] = rftg.P_interp
# interp[2] = rftg.H_interp

#pdb.set_trace()
# xrinterp = np.amax(P_interp, axis=0)
# yeinterp = np.max(T_interp, axis=1)
xrinterp = np.amax(rinterp, axis=0)
yeinterp = np.max(einterp, axis=1)

def weigh(rho_target, e_target, c):
	def find_position(rho_target, e_target):
		#find the indicies of the immediate upper,lower & left,right values
		idy = np.argmax(yeinterp>e_target)-1
		idx = np.argmax(rinterp[idy-1]>rho_target)
		x3 = rho_act_int[idy][idx]
		x4 = rho_act_int[idy][idx+1]
		x1 = rho_act_int[idy+1][idx]
		x2 = rho_act_int[idy+1][idx+1]
		y3   = e_interp[idy][idx]
		y4   = y3
		y2   = e_interp[idy+1][idx+1]
		y1   = y2
		x31 = x3 - x1
		y31 = y3 - y1
		x42 = x4 - x2
		y42 = y4 - y2
		x21 = x2 - x1
		y21 = y2 - y1
		x43 = x4 - x3
		y43 = y4 - y3

		Px  = rho_target
		Py  = e_target

		if (np.shape(Px) != ()):
			Px = Px[0]
			Py = Py[0]

		xy_m = [[x31, x21], [y31, y21]]
		p_m  = [[Px - x1], [Py - y1]]

		xy_m = np.asarray(xy_m)
		p_m = np.asarray(p_m)
		param = np.linalg.solve(xy_m, p_m)
		t     = param[0]
		s     = param[1]
		return idx, idy, t, s

	def weighted_interp(idx, idy, t, s, c):

		P1 = interp[c][idy+1][idx]
		P2 = interp[c][idy+1][idx+1]
		P3 = interp[c][idy][idx]
		P4 = interp[c][idy][idx+1]

		#print "The P1, P2, P3, P4 are: ", P1, P2, P3, P4
		#print "The values of t, s, are: ", t, s

		P = P1*(1-s)*(1-t) + P2*(s)*(1-t) + P3*(1-s)*t + P4*s*t

		# plt.scatter(rinterp, einterp)
		# plt.scatter(rho_target, e_target, color='green', marker='*')
		# plt.scatter(rinterp[idy][idx], einterp[idy][idx], color='red')
		# plt.scatter(rinterp[idy][idx+1], einterp[idy][idx+1], color='red')
		# plt.scatter(rinterp[idy+1][idx], einterp[idy+1][idx], color='red')
		# plt.scatter(rinterp[idy+1][idx+1], einterp[idy+1][idx+1], color='red')
		# plt.xlabel('Density (kg/m3)')
		# plt.ylabel('Internal Energy (KJ/Kg)')
		# plt.title('Internal Energy vs Density')
		# plt.show()
		# pdb.set_trace()
		return P

	idx, idy, t, s  = find_position(rho_target, e_target)
	if (c=='T'):
		return weighted_interp(idx, idy, t, s, 'T')
	if (c=='P'):
		return weighted_interp(idx, idy, t, s, 'P')
	if (c=='H'):
		return weighted_interp(idx, idy, t, s, 'H')
	if (c=='Rho'):
		return weighted_interp(idx, idy, t, s, 'Rho')
	if (c=='E'):
		return weighted_interp(idx, idy, t, s, 'E')
	if (c=='CV'):
		return weighted_interp(idx, idy, t, s, 'CV')
	if (c=='CP'):
		return weighted_interp(idx, idy, t, s, 'CP')
	if (c=='a'):
		return weighted_interp(idx, idy, t, s, 'a')

'''
Univariate Newton Solver
'''

def derivative(p, e, h):
	return (p(e+h) - p(e-h)) / (2.0*h)

def interpo(e):
	return weigh(rho_target, e, 'P') - P_target

def solver(p, e0, h):
	lastE = e0
	nextE = lastE + 10*h
	while (abs(lastE - nextE) > h):
		newP = p(nextE)
		print "p(", nextE, ") = ", newP
		lastE = nextE
		nextE = lastE - newP / derivative(p, lastE, h)
	return nextE

# Efound = solver(interpo, 230, 0.001)
# print "solution: e = ", Efound

# T_cal = weigh(rho_target, Efound, "T")
# print "The T(",P_target, " , " , rho_target, ") = ", T_cal

'''
Multivariate Newton solver
'''

T_target = input("Enter the target Temperature in K : ")
P_target = input("Enter the target Pressure in Pa : ")

pt = [P_target, T_target]

def jacobian(x1, x2, h):
	dpdrho = (p(x1+h, x2) - p(x1, x2))/h
	dpde   = (p(x1, x2+h) - p(x1, x2))/h
	dtdrho = (t(x1+h, x2) - t(x1, x2))/h
	dtde   = (t(x1, x2+h) - t(x1, x2))/h
	jacobian = np.matrix([[dpdrho[0], dpde[0]], [dtdrho[0], dtde[0]]])
	return jacobian

def twodinterp(x1, x2):
	pt = [weigh(x1, x2, 'P'), weigh(x1, x2, 'T')] 
	pt = np.transpose(np.concatenate((pt[0], pt[1])))- [P_target, T_target]
	return pt

def p(x1, x2):
	return weigh(x1, x2, 'P')
def t(x1, x2):
	return weigh(x1, x2, 'T')

def solver2D(pt, x1, x2, h):
	lastX = [x1, x2]
	nextX = [1,2]
	nextX[0] = lastX[0] + 10*h
	nextX[1] = lastX[1] + 10*h
	while ((abs(lastX[0] - nextX[0]) > h)&(abs(lastX[1] - nextX[1])> h)):
		if (np.shape(nextX[0])!=()):
			nextX[0] = nextX[0][0]
			nextX[1] = nextX[1][0]
		rhoe = pt(nextX[0], nextX[1])
		newX = np.matrix([[rhoe[0]], [rhoe[1]]])
		print "PT(", nextX[0], ",", nextX[1], ") = ", newX
		delX = np.linalg.solve(jacobian(nextX[0], nextX[1], h), newX)
		delX = np.ndarray.tolist(delX)
		lastX = np.copy(nextX)
		nextX[0] = lastX[0] - delX[0]
		nextX[1] = lastX[1] - delX[1]
		nextX = np.asarray(nextX)
		nextX = np.ndarray.tolist(nextX)
	return nextX

XFound = solver2D(twodinterp, 15, 235, 0.01)
print "solution: [rho, e] = ", XFound
Rep = raw_input("Do you wanna proceed to calculate any more parameters ? (Y/N)")
if (Rep == 'Y'):
	c = raw_input("Enter the parameter symbol. For Ex: a- speed of sound, H- Enthalpy ..")
	prop_cal = weigh(XFound[0][0], XFound[1][0], c)
	print "The value is: ", prop_cal





# import matplotlib.pyplot as plt
# fig = plt.figure()
# fig.suptitle('Internal Energy (E) vs Density (Rho)', fontsize = 14, fontweight = 'bold')
# ax = fig.add_subplot(111)
# plt.scatter(rho_target, e_target, color = 'blue')
# #ax.set_title("Internal Energy (E) vs Density (Rho)")
# ax.set_xlabel("Density (Kg / m3)")
# ax.set_ylabel("Internal Energy (KJ / Kg)")
# plt.scatter(val_neigh[0], val_neigh[2], color = 'red')
# plt.scatter(val_neigh[0], val_neigh[3], color = 'red')
# plt.scatter(val_neigh[1], val_neigh[2], color = 'red')
# plt.scatter(val_neigh[1], val_neigh[3], color = 'red')
# plt.scatter(rho_target, e_target, color = 'green', marker = "*")
# #ax.annotate('(rho_target, e_target)', xy = (rho_target, e_target), xytext = (24, 190), arrowprops=dict(facecolor='black', shrink = 0.05))
# plt.show()

#Important
#np.ndarray.flatten(T_interp).shape
#interpolator = NearestNDInterpolator(np.vstack(np.dstack((rho_act_int, e_interp))), np.ndarray.flatten(T_interp))
