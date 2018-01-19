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

rho_target = input("Enter the Density target in kg/m3 ")
e_target = input("Enter the energy target in KJ/Kg  ")

f=open("interp.pickle", "rb" )
e_interp,T_interp, P_interp, H_interp, rho_act_int = pickle.load(f)

#pdb.set_trace()
Nrho= rftg.NRho

rinterp = rho_act_int
#einterp = rftg.e_interp

interp = {'T': rftg.T_interp, 'P': rftg.P_interp, 'H': rftg.H_interp}
# interp[0] = rftg.T_interp
# interp[1] = rftg.P_interp
# interp[2] = rftg.H_interp

#pdb.set_trace()
xrinterp = np.amax(rinterp, axis=0)
yeinterp = np.max(einterp, axis=1)


def find_position(rho_target, e_target):
	#find the indicies of the immediate upper,lower & left,right values
	idy = np.argmax(yeinterp>e_target)-1
	idx = np.argmax(rinterp[idy-1]>rho_target)
	x3 = rftg.rho_act_int[idy][idx]
	x4 = rftg.rho_act_int[idy][idx+1]
	x1 = rftg.rho_act_int[idy+1][idx]
	x2 = rftg.rho_act_int[idy+1][idx+1]
	y3   = rftg.e_interp[idy][idx]
	y4   = y3
	y2   = rftg.e_interp[idy+1][idx+1]
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

	xy_m = [[x31, x21], [y31, y21]]
	p_m  = [Px - x1], [Py - y1]

	xy_m = np.asarray(xy_m)
	p_m = np.asarray(p_m)
	param = np.linalg.solve(xy_m, p_m)
	t     = param[0]
	s     = param[1]
	#pdb.set_trace()
	return idx, idy, t, s

def weighted_interp(idx, idy, t, s, c):

	P1 = interp[c][idy+1][idx]
	P2 = interp[c][idy+1][idx+1]
	P3 = interp[c][idy][idx]
	P4 = interp[c][idy][idx+1]

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
	#pdb.set_trace()
	#print P
	return P


# def get_values_of_neighbours(idx, idy):
# 	#Get the four neighbours for each interpolation point
# 	rho0 = rftg.rho_interp[idy][idx]
# 	rho1 = rftg.rho_interp[idy+1][idx+1]
# 	e0   = rftg.e_interp[idy][idx]
# 	e1   = rftg.e_interp[idy+1][idx+1]
# 	print rho0, rho1, e0, e1
# 	#pdb.set_trace()

# 	T00  = rftg.T_interp[idy][idx]
# 	T01  = rftg.T_interp[idy+1][idx+1]
# 	T10  = rftg.T_interp[idy+1][ idx]
# 	T11  = rftg.T_interp[idy+1][ idx+1]

# 	# idy -= 2
# 	#idx += 1
# 	print idx, idy
# 	P00  = rftg.P_interp[idy][idx]
# 	print P00
# 	P01  = rftg.P_interp[idy+1][idx+1]
# 	P10  = rftg.P_interp[idy+1][ idx]
# 	P11  = rftg.P_interp[idy+1][ idx+1]

# 	print "Temperatures: ", T00, T01, T10, T11
# 	print "Pressures: ", P00, P01, P10, P11

# 	# rho_a0 = rftg.rho_act_int[idy][idx]
# 	# rho_a1  = rftg.rho_act_int[idy+1][idx+1]

# 	return rho0, rho1, e0, e1, T00, T01, T10, T11, P00, P01, P10, P11

# def find_interpolation_weights(rho_target, e_target, val_neigh):
# 	#Coefficients for weighting between lower and upper bounds
# 	x31 = x3 - x1
# 	y31 = y3 - y1
# 	x42 = x4 - x2
# 	y42 = y4 - y2
# 	x21 = x2 - x1
# 	y21 = y2 - y1
# 	x43 = x4 - x3
# 	y43 = y4 - y3
# 	Ax = x1 + x31*t
# 	Ay = y1 + y31*t
# 	bx = x2 + x42*t
# 	by = y2 + y42*t
# 	cx = x1 + x21*s
# 	cy = y1 + y21*s
# 	dx = x3 + x43*s
# 	dy = y3 + y43*s

# 	# Px = Ax + (bx - Ax)*s
# 	# Py = Ay + (by - Ay)*s

# 	Px = rho_target
# 	Py = e_target

# 	xy_m = np.array([x31, x21], [y31, y21])
# 	p_m  = np.array([Px - x1], [Py - y1])

# 	param = np.linalg.solve(xy_m, p_m)
# 	t     = param[0]
# 	s     = param[1]
# 	pdb.set_trace()

# 	P1 = rftg.T_interp[idy+1][idx]
# 	P2 = rftg.T_interp[idy+1][idx+1]
# 	P3 = rftg.T_interp[idy][idx]
# 	P4 = rftg.T_interp[idy][idx+1]


# 	# A = x31*y42 - y31*x42
# 	# B = Py*(x42-x31) - Px*(y42 - y31) + x31*y2 - y31*x2 + x1*y42 - y1*y42
# 	# C = Py*x21 - Px*y21 + x1*y2 - x2*y1

# 	# T = np.roots(np.asarray(A, B, C))
# 	# if 0<T[0]<1:
# 	# 	t = T[0]
# 	# if 0<T[1]<1:
# 	# 	t = T[1]
# 	# s = (Py - Ay)/(By - Ay)
# 	P = P1*(1-s)*(1-t) + P2*(s)*(1-t) + P3*(1-s)*t + P4*s*t




# 	rho0  = val_neigh[0]
# 	rho1  = val_neigh[1]
# 	e0    = val_neigh[2]
# 	e1    = val_neigh[3]
# 	alpha = (rho_target - rho0) / (rho1 - rho0)
# 	beta  = (e_target - e0) / (e1 - e0)
# 	return alpha, beta

idx, idy, t, s  = find_position(rho_target, e_target)
# val_neigh = get_values_of_neighbours(idx, idy)
# alpha, beta = find_interpolation_weights(rho_target, e_target, val_neigh)
T_cal = weighted_interp(idx, idy, t, s, 'T')
P_cal = weighted_interp(idx, idy, t, s, 'P')
H_cal = weighted_interp(idx, idy, t, s, 'H')


print "The input density and energy values in Kg/m3, KJ/Kg K are: ", rho_target, e_target

#new T is T_cal and use e_target


print "The interpolated Temperature and density from the tabulated data are: ", T_cal, "K", P_cal/1E6, "MPa", H_cal, "KJ/Kg"


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
