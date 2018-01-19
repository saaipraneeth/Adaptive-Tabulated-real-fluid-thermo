'''
########################################################################################################
#  A Program to read the rho, E table and find the interpolated thermodynamic properties for a given point
#  Written: August 2016
#  Author: Sai
########################################################################################################
'''
import numpy as np
#import realFluidTableGeneration as rftg
import math
import pdb
import matplotlib.pyplot as plt
import pickle

f=open("interp.pickle", "rb" )
e_interp,T_interp, P_interp, H_interp, rho_interp, cv_interp, cp_interp, a_interp, k_interp, mu_interp, rhohomo, ehomo = pickle.load(f)
interp = {'T': T_interp, 'P': P_interp, 'H': H_interp, 'Rho': rho_interp, 'E': e_interp, 'CV': cv_interp, 'CP': cp_interp, 'a': a_interp, 'k': k_interp, 'mu': mu_interp}

xrinterp = np.amax(rho_interp, axis=0)
yeinterp = np.max(e_interp, axis=1)

rho_target = input("Enter the Density target in kg/m3 ")
e_target = input("Enter the energy target in KJ/Kg  ")

def find_position(rho_target, e_target):
	#find the indicies of the immediate upper,lower & left,right values
	x_u_i = np.argmax(xrinterp>rho_target)
	y_u_i = np.argmax(yeinterp>e_target)
	print x_u_i, y_u_i
	idx = x_u_i
	idy = y_u_i
	#pdb.set_trace()
	return idx, idy


def get_values_of_neighbours(idx, idy):
	#Get the four neighbours for each interpolation point
	rho0 = rho_interp[idy-1][idx-1]
	rho1 = rho_interp[idy][idx]
	e0   = e_interp[idy-1][idx-1]
	e1   = e_interp[idy][idx]
	print rho0, rho1, e0, e1
	#pdb.set_trace()

	T00  = T_interp[idy-1][idx-1]
	T01  = T_interp[idy][idx]
	T10  = T_interp[idy-1][idx-1]
	T11  = T_interp[idy][idx]

	P00  = P_interp[idy-1][idx-1]
	P01  = P_interp[idy][idx]
	P10  = P_interp[idy-1][idx-1]
	P11  = P_interp[idy][idx]

	plt.scatter(rho_interp, e_interp)
	plt.scatter(rho0, e0, color='r')
	plt.scatter(rho1, e1, color='r' )
	plt.scatter(rho0, e1, color='r')
	plt.scatter(rho1, e0, color='r')
	plt.scatter(rho_target, e_target, color='g')
	plt.show()


	return rho0, rho1, e0, e1, T00, T01, T10, T11, P00, P01, P10, P11

def find_interpolation_weights(rho_target, e_target, val_neigh):
	#Coefficients for weighting between lower and upper bounds
	rho0  = val_neigh[0]
	rho1  = val_neigh[1]
	e0    = val_neigh[2]
	e1    = val_neigh[3]
	alpha = (rho_target - rho0) / (rho1 - rho0)
	beta  = (e_target - e0) / (e1 - e0)
	return alpha, beta

def find_target_from_weight_coeff(alpha, beta, val_neigh):
	#Bilinear interpolation of rho and e
	#get the target values for Temp, Pressure, Enthalpy.
	T00 = val_neigh[4]
	T01 = val_neigh[5]
	T10 = val_neigh[6]
	T11 = val_neigh[7]
	P00 = val_neigh[8]
	P01 = val_neigh[9]
	P10 = val_neigh[10]
	P11 = val_neigh[11]
	dtx = T10 - T00
	dty = T01 - T00
	dpx = P10 - P00
	dpy = P01 - P00
	T_cal  = T00 + alpha*dtx + beta*dty + alpha*beta*(T11 - dtx - dty - T00)
	P_cal  = P00 + alpha*dpx + beta*dpy + alpha*beta*(P11 - dpx - dpy - P00)
	return T_cal, P_cal

idx, idy  = find_position(rho_target, e_target)
val_neigh = get_values_of_neighbours(idx, idy)
alpha, beta = find_interpolation_weights(rho_target, e_target, val_neigh)
T_cal, P_cal = find_target_from_weight_coeff(alpha, beta, val_neigh)





def weigh(rho_target, e_target, c):
	def find_position(rho_target, e_target):
		#find the indicies of the immediate upper,lower & left,right values
		x_u_i = np.argmax(xrinterp>rho_target)
		y_u_i = np.argmax(yeinterp>e_target)
		print x_u_i, y_u_i
		idx = x_u_i-1
		idy = y_u_i-1
		#pdb.set_trace()
		return idx, idy

	def weighted_interp(idx, idy, c):

		P1 = interp[c][idy-1][idx-1]
		P2 = interp[c][idy-1][idx]
		P3 = interp[c][idy][idx]
		P4 = interp[c][idy][idx-1]

		dtx = P3 - P1
		dty = P2 - P1
		dpx = P3 - P1
		dpy = P2 - P1

		rho0 = rho_interp[idy-1][idx-1]
		rho1 = rho_interp[idy][idx]
		e0   = e_interp[idy-1][idx-1]
		e1   = e_interp[idy][idx]

		alpha = (rho_target - rho0) / (rho1 - rho0)
		beta  = (e_target - e0) / (e1 - e0)

		T00  = T_interp[idy-1][idx-1]
		T01  = T_interp[idy][idx]
		T10  = T_interp[idy-1][idx-1]
		T11  = T_interp[idy][idx]
		#print "The P1, P2, P3, P4 are: ", P1, P2, P3, P4
		#print "The values of t, s, are: ", t, s

		#P = P1*(1-s)*(1-t) + P2*(s)*(1-t) + P3*(1-s)*t + P4*s*t
		P  = P1 + alpha*dtx + beta*dty + alpha*beta*(P4 - dtx - dty - P1)
		if (raw_input("Do you want to see the plot ? ") == 'Y'):
			rinterp = rho_interp
			einterp = e_interp
			plt.scatter(rinterp, einterp)
			plt.scatter(rho_target, e_target, color='green', marker='*')
			plt.scatter(rinterp[idy][idx], einterp[idy][idx], color='red')
			plt.scatter(rinterp[idy][idx+1], einterp[idy][idx+1], color='red')
			plt.scatter(rinterp[idy+1][idx], einterp[idy+1][idx], color='red')
			plt.scatter(rinterp[idy+1][idx+1], einterp[idy+1][idx+1], color='red')
			plt.xlabel('Density (kg/m3)')
			plt.ylabel('Internal Energy (KJ/Kg)')
			plt.title('Internal Energy vs Density')
			plt.show()
			pdb.set_trace()
		return P

	idx, idy  = find_position(rho_target, e_target)
	if (c=='T'):
		return weighted_interp(idx, idy,'T')
	if (c=='P'):
		return weighted_interp(idx, idy, 'P')
	if (c=='H'):
		return weighted_interp(idx, idy, 'H')
	if (c=='Rho'):
		return weighted_interp(idx, idy, 'Rho')
	if (c=='E'):
		return weighted_interp(idx, idy, 'E')
	if (c=='CV'):
		return weighted_interp(idx, idy, 'CV')
	if (c=='CP'):
		return weighted_interp(idx, idy, 'CP')
	if (c=='a'):
		return weighted_interp(idx, idy, 'a')



# T_cal = weigh(rho_target, e_target, 'T')
# P_cal = weigh(rho_target, e_target, 'P')
print "The input density and energy values in Kg/m3, KJ/Kg K are: ", rho_target, e_target
print "The interpolated Temperature and density from the tabulated data are: ", T_cal, "K", P_cal/1E6, "MPa"