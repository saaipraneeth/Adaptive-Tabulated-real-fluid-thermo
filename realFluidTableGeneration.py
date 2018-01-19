'''
########################################################################################################
#  A Program to read the NIST website, extract isothermal or isobaric data and reformat into a rho-E table
#  Written: April 2013
#  Author: jph
#  Version: 0.23
########################################################################################################
'''

# Generic modules
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import pdb
from pdb import set_trace as keyboard
import requests
from scipy import interpolate
import pickle
# User defined modules
import NIST_reader as NIST
#import netCDF_writer as ncdf
import tecplotOutput as tec
import criticalProperties as cP

'''
##################################################
# User defined quantities
##################################################
'''
fluid  ='O2'         # fluid (H2O, CO, H2, O2 and N2)
isoType='isotherm'   # extract along isobars or isothermal lines (currently hard coded for isotherm in main)
Pmin   =0.8E5        # Minimal pressure to be extracted
Pmax   =4.0E6        # Maximal pressure to be extracted
Tmin   =250.0        # Minimal temperature to be extracted
Tmax   =350.0        # Maximal temperature to be extracted
NbT    =50          # Number of isothermal lines (neglected if isobar used)
NbP    =50          # Number of isobar lines (neglected if isothermal used)
T      =444          #
P      =1.0
NE     =100           # Number of iso energy lines in final table
NRho   =100           # Number of isochoric lines in final table
NVar   =14           # Number of variables (remains unchanged)
ending ='.png'       # Output format of figures (matplotlib)
outpath="../res/"         # Path for figure output
readData=True
'''
##################################################
# User FUNCTIONS (some fcts are not used -> need to clean this up)
##################################################
'''

'''
Interpolates a 1D field (could use scipy.interp1d instead) but this function works with repeated values
'''
def interpolationWeights(Vect, target):
    nmax=np.size(Vect)
    idx=np.argmin(np.abs(Vect - target))
    pt=[-1,-1]
    # for more than 1 identical value, loop over all values
    if np.size(Vect[Vect==Vect[idx]])>1:
        pt=slowSearch(Vect,target)
        weight=(target-Vect[pt[0]])/(Vect[pt[1]]-Vect[pt[0]])
        return weight,pt
    # for exact values, no interpolation needed
    if np.abs(Vect[idx]-target)<1E-6:
        return 0.0,[idx,0]

    # find neighbors
    if idx==0:
        pt[0] = 0
        pt[1] = 1
    elif idx==nmax-1:
        pt[0] = nmax-2
        pt[1] = nmax-1
    elif target>Vect[idx] and target<=Vect[idx+1]:
        pt[0]=idx
        pt[1]=idx+1
    elif target<Vect[idx] and target>=Vect[idx+1]:
        pt[0]=idx
        pt[1]=idx+1
    elif target<=Vect[idx] and target>Vect[idx-1]:
        pt[0]=idx-1
        pt[1]=idx
    elif target>=Vect[idx] and target<Vect[idx-1]:
        pt[0]=idx-1
        pt[1]=idx
    weight=(target-Vect[pt[0]])/(Vect[pt[1]]-Vect[pt[0]])

    #check if target is between neighbors
    if ((target>Vect[pt[0]] and target<=Vect[pt[1]]) or (target<Vect[pt[0]] and target>=Vect[pt[1]])):
        pt=slowSearch(Vect,target)
        weight=(target-Vect[pt[0]])/(Vect[pt[1]]-Vect[pt[0]])
        return weight,pt
    return weight,pt

'''
loops over all possible values in the vector matrix
'''
def slowSearch(Vect,target):
    nmax=np.size(Vect)-1
    found=0
    pt=[-1,-1]
    for idx in range(0,nmax):
        if (target<Vect[idx] and target>=Vect[idx+1]) or (target>Vect[idx] and target<=Vect[idx+1]):
            pt[0]=idx
            pt[1]=idx+1
            found+=1
    if found>1:
       print "Found more than 1 acceptable iteration point...",Vect, target,pt
    elif found==0:
        print "Found no acceptable values..."
    else:
        return pt

def find_closest(A, target):
    #A must be sorted
    return np.argmin(np.abs(AA - target))


'''
##################################################
# MAIN CODE
##################################################
'''
#Defines the start values of the arrays
startValues = [-1] * (NbT+1)

dataIsotherm  =np.zeros([NVar,NbT*NbP*10])
thermo        =cP.solutionVector()
thermo.setThermodynamics(fluid)
RHOarray_maxRho=np.zeros([NbT])
RHOarray_minRho=np.zeros([NbT])
Earray_minRho=np.zeros([NbT])
Earray_maxRho=np.zeros([NbT])
outputName = fluid+"_Pmin"+str(int(Pmin/1E5))+"_Pmax"+str(int(Pmax/1E5))+".pickle"
Tarray=np.zeros([NbT,NbP+1])
Parray=np.zeros([NbT,NbP+1])
RHOarray=np.zeros([NbT,NbP+1])
Harray=np.zeros([NbT,NbP+1])
Earray=np.zeros([NbT,NbP+1])
CVarray=np.zeros([NbT,NbP+1])
CParray=np.zeros([NbT,NbP+1])
Aarray=np.zeros([NbT,NbP+1])
Karray=np.zeros([NbT,NbP+1])
Muarray=np.zeros([NbT,NbP+1])


if readData:
  for ii,T in enumerate(np.linspace(Tmin,Tmax,NbT)):
      dataNIST=NIST.readNIST(isoType, fluid, T, P, Tmin,Tmax,Pmin/1.0E6,Pmax/1.0E6,NbP)
      RHOarray_minRho[ii] = dataNIST[NIST.colNIST('rho'),0]
      RHOarray_maxRho[ii] = dataNIST[NIST.colNIST('rho'),-1]
      Earray_minRho[ii]   = dataNIST[NIST.colNIST('E'),0]
      Earray_maxRho[ii]   = dataNIST[NIST.colNIST('E'),-1]
      RHOarray_minT=dataNIST[NIST.colNIST('rho'),:]
      #keyboard()
      Tarray[ii,:]  =dataNIST[NIST.colNIST('T'),:]
      Parray[ii,:]  =dataNIST[NIST.colNIST('P'),:]
      Harray[ii,:]  =dataNIST[NIST.colNIST('H'),:]
      RHOarray[ii,:]=dataNIST[NIST.colNIST('rho'),:]
      Earray[ii,:]  =dataNIST[NIST.colNIST('E'),:]
      CVarray[ii,:] =dataNIST[NIST.colNIST('Cv'),:]
      Aarray[ii,:]  =dataNIST[NIST.colNIST('C'),:]
      CParray[ii,:] =dataNIST[NIST.colNIST('Cp'),:]
      Karray[ii,:]  =dataNIST[NIST.colNIST('kappa'),:]
      Muarray[ii,:] =dataNIST[NIST.colNIST('mu'),:]
    # Outputs a bunch of isothermal lines and writes figures
      # plt.figure(33)
      # plt.plot(dataNIST[NIST.colNIST('rho'),:],dataNIST[NIST.colNIST('P'),:]/1E6,color='k')
      # plt.figure(34)
      # plt.plot(dataNIST[NIST.colNIST('E'),:],dataNIST[NIST.colNIST('P'),:]/1E6,color='k')
      # plt.figure(35)
      # plt.plot(dataNIST[NIST.colNIST('rho'),:],dataNIST[NIST.cloNIST('E'),:],color='k')
      # plt.figure(36)
      # plt.plot(dataNIST[NIST.colNIST('H'),:],dataNIST[NIST.colNIST('S'),:],color='k')
      # plt.figure(37)
      # plt.plot(dataNIST[NIST.colNIST('P'),:]/1E6,dataNIST[NIST.colNIST('V'),:],color='k')
      # plt.figure(38)
      # plt.plot(dataNIST[NIST.colNIST('rho'),:],dataNIST[NIST.colNIST('T'),:],color='k')
      # f=open(outputName, "wb" )
      # pickle.dump((Tarray,Parray,RHOarray,Earray,Harray, CVarray, CParray, Aarray, Karray, Muarray), f )
      # f.close()
else:
      f=open(outputName, "rb" )
      Tarray,Parray,RHOarray,Earray,Harray, CVarray, CParray, Aarray, Karray, Muarray = pickle.load(f)


# Creates a rho-E table
###########################################
# Get homogeneous energy divisions
###########################################
E_min= min(min(Earray[:,0]),min(Earray[:,-1]))*0.99
E_max=max(max(Earray[:,0]),max(Earray[:,-1]))*1.01
RHO_max = max(max(RHOarray[:,0]),max(RHOarray[:,-1]))
RHO_min = min(min(RHOarray[:,0]),min(RHOarray[:,-1]))
dataHomo=np.zeros([NE,NRho,NVar])
dataTemp=np.zeros([NRho,NbT,NVar])

Ehomo=np.linspace(E_min,E_max,NE)
Rhohomo=np.linspace(RHO_min,RHO_max,NRho)
# Set start and end of the density for each energy level.
from scipy import interpolate

for ii,T in enumerate(np.linspace(Tmin,Tmax,NbT)):
  #print ii
  lrho = RHOarray[:,ii]
  lrho_min=min(lrho)
  lrho_max=max(lrho)
  window= (Rhohomo>lrho_min) &(Rhohomo<lrho_max)
  ind_window=np.where(window)
  #pdb.set_trace()
  sort = np.argsort(lrho)
  lrho = lrho[sort]
  #print ii
  #print window
  #print ind_window
  f_p =   interpolate.interp1d(lrho, Parray[sort,ii])
  f_e =   interpolate.interp1d(lrho, Earray[sort,ii])
  f_T =   interpolate.interp1d(lrho, Tarray[sort,ii])
  f_H =   interpolate.interp1d(lrho, Harray[sort,ii])
  f_R =   interpolate.interp1d(lrho, RHOarray[sort, ii])
  f_CV=   interpolate.interp1d(lrho, CVarray[sort, ii])
  f_CP=   interpolate.interp1d(lrho, CParray[sort, ii])
  f_A =   interpolate.interp1d(lrho, Aarray[sort, ii])
  f_K =   interpolate.interp1d(lrho, Karray[sort, ii])
  f_MU=   interpolate.interp1d(lrho, Muarray[sort, ii])
 # dataTemp[ind_window,ii,0] = Rhohomo
  dataTemp[ind_window,ii,1] = f_p(Rhohomo[window])
  dataTemp[ind_window,ii,2] = f_e(Rhohomo[window])
  dataTemp[ind_window,ii,3] = f_T(Rhohomo[window])
  dataTemp[ind_window,ii,4] = f_H(Rhohomo[window])
  dataTemp[ind_window,ii,5] = Rhohomo[window]
  dataTemp[ind_window,ii,6] = f_R(Rhohomo[window])
  dataTemp[ind_window,ii,7] = f_CV(Rhohomo[window])
  dataTemp[ind_window,ii,8] = f_CP(Rhohomo[window])
  dataTemp[ind_window,ii,9] = f_A(Rhohomo[window])
  dataTemp[ind_window,ii,10] = f_K(Rhohomo[window])
  dataTemp[ind_window,ii,11] = f_MU(Rhohomo[window])
  #print Rhohomo[window]
#pdb.set_trace()
for ii in range(1,NRho-1):
    lenergy = dataTemp[ii,:,2]
    nonZero=lenergy.nonzero()
    count_nonzero= np.count_nonzero(lenergy)
    if count_nonzero>0:
      lenergy_min=min(lenergy[nonZero])
      lenergy_max=max(lenergy[nonZero])
    window= (Ehomo>lenergy_min) &(Ehomo<lenergy_max)
    ind_window=np.where(window)
    if count_nonzero>1:
      #print ii, nonZero,dataTemp[ii,nonZero,3]
      #print count_nonzero
      f_p =   interpolate.interp1d(lenergy[nonZero],dataTemp[ii,nonZero,1])
      f_e =   interpolate.interp1d(lenergy[nonZero],dataTemp[ii,nonZero,2])
      f_T =   interpolate.interp1d(lenergy[nonZero],dataTemp[ii,nonZero,3])
      f_H =   interpolate.interp1d(lenergy[nonZero],dataTemp[ii,nonZero,4])
      f_R =   interpolate.interp1d(lenergy[nonZero],dataTemp[ii,nonZero,5])
      f_r =   interpolate.interp1d(lenergy[nonZero],dataTemp[ii,nonZero,6])
      f_CV=   interpolate.interp1d(lenergy[nonZero],dataTemp[ii,nonZero,7])
      f_CP=   interpolate.interp1d(lenergy[nonZero],dataTemp[ii,nonZero,8])
      f_A =   interpolate.interp1d(lenergy[nonZero],dataTemp[ii,nonZero,9])
      f_K =   interpolate.interp1d(lenergy[nonZero],dataTemp[ii,nonZero,10])
      f_MU=   interpolate.interp1d(lenergy[nonZero],dataTemp[ii,nonZero,11])
      dataHomo[ind_window,ii,1]=f_p(Ehomo[window])
      dataHomo[ind_window,ii,2]=f_e(Ehomo[window])
      dataHomo[ind_window,ii,3]=f_T(Ehomo[window])
      dataHomo[ind_window,ii,4]=f_H(Ehomo[window])
      dataHomo[ind_window,ii,5]=f_R(Ehomo[window])
      dataHomo[ind_window,ii,6]=f_r(Ehomo[window])
      dataHomo[ind_window,ii,7]=f_CV(Ehomo[window])
      dataHomo[ind_window,ii,8]=f_CP(Ehomo[window])
      dataHomo[ind_window,ii,9]=f_A(Ehomo[window])
      dataHomo[ind_window,ii,10]=f_K(Ehomo[window])
      dataHomo[ind_window,ii,11]=f_MU(Ehomo[window])
#pdb.set_trace()
e_interp = dataHomo[:,:,2]
T_interp = dataHomo[:,:,3]
P_interp = dataHomo[:,:,1]
H_interp = dataHomo[:,:,4]
cv_interp = dataHomo[:,:,7]
cp_interp = dataHomo[:,:,8]
a_interp = dataHomo[:,:,9]
rho_interp = dataHomo[:,:,5]
#rho_act_int = dataHomo[:,:,6]
k_interp = dataHomo[:,:,10]
mu_interp = dataHomo[:,:,11]

# f=open("interp.pickle", "wb" )
# pickle.dump((e_interp,T_interp, P_interp, H_interp, rho_interp, cv_interp, cp_interp, a_interp, k_interp, mu_interp, Rhohomo, Ehomo), f )
# f.close()
pdb.set_trace()
# X, Y = np.meshgrid( Rhohomo,Ehomo)
# plt.scatter(X, Y, dataHomo[:,:,3])
# #plt.plot(RHOarray, Earray)
# plt.show()
