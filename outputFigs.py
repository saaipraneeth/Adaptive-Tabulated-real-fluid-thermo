########################################################################################################
#  A Program to read the NIST website and output figures
########################################################################################################
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.mlab as ml
import pdb
import requests
import NIST_reader as NIST
#import netCDF_writer as ncdf
import criticalProperties as cP
from bs4 import BeautifulSoup
from scipy import interpolate
import tecplotOutput as tec
import scipy as sp


#User defined quantities
fluid='O2'
isoType='isotherm'
Pmin=1.0E6
Pmax=5.0E6
Tmin=80
Tmax=840.0
NbT=25
NbP=100
T=444
P=1.0
NVar=14
ending='.svg'

#Defines the start values of the arrays
startValues = [-1] * (NbT+1)  
if isoType=='isotherm':
    rangeNIST=np.linspace(Tmin,Tmax,NbT)
elif isoType=='isobar':
    rangeNIST=np.linspace(Pmin,Pmax,NbP)

dataIsotherm=np.zeros([NVar,NbT*NbP*10])
for ii,valThermo in enumerate(rangeNIST):
    if isoType=='isotherm':
        T=valThermo
    elif isoType=='isobar':
        P=valThermo
    dataNIST=NIST.readNIST(isoType, fluid, T, P/1.0E6, Tmin,Tmax,Pmin/1.0E6,Pmax/1.0E6,NbP)
    if ii==0: 
        Tarray=dataNIST[NIST.colNIST('T'),:]
        Parray=dataNIST[NIST.colNIST('P'),:]
        Harray=dataNIST[NIST.colNIST('H'),:]
        RHOarray=dataNIST[NIST.colNIST('rho'),:]
        Earray=dataNIST[NIST.colNIST('E'),:]
        nPts=np.size(dataNIST[0,:])
    else:     
        Tarray=np.append(Tarray,dataNIST[NIST.colNIST('T'),:])
        Parray=np.append(Parray,dataNIST[NIST.colNIST('P'),:])
        Harray=np.append(Harray,dataNIST[NIST.colNIST('H'),:])
        Earray=np.append(Harray,dataNIST[NIST.colNIST('E'),:])
        RHOarray=np.append(RHOarray,dataNIST[NIST.colNIST('rho'),:])  

    plt.figure(33)
    plt.plot(dataNIST[NIST.colNIST('rho'),:],dataNIST[NIST.colNIST('P'),:]/1E6,color='k')
    plt.figure(34)
    plt.plot(dataNIST[NIST.colNIST('E'),:],dataNIST[NIST.colNIST('P'),:]/1E6,color='k')
    plt.figure(35)
    plt.plot(dataNIST[NIST.colNIST('rho'),:],dataNIST[NIST.colNIST('E'),:],color='k')
    plt.figure(36)
    plt.plot(dataNIST[NIST.colNIST('H'),:],dataNIST[NIST.colNIST('S'),:],color='k')
    plt.figure(37)
    plt.plot(dataNIST[NIST.colNIST('P'),:]/1E6,dataNIST[NIST.colNIST('V'),:],color='k')
    plt.figure(38)
    plt.plot(dataNIST[NIST.colNIST('rho'),:],dataNIST[NIST.colNIST('T'),:],color='k')

prefix=isoType+'_'+fluid
fig=plt.figure(33)
plt.xlabel('rho (kg/m3)')
plt.ylabel('P (MPa)')
plt.savefig(prefix+'_rhoP'+ending)
plt.figure(34)
plt.xlabel('E (kJ/kg)')
plt.ylabel('P (MPa)')
plt.savefig(prefix+'_eP'+ending)
plt.figure(35)
plt.xlabel('rho (kg/m3)')
plt.ylabel('E (kJ/kg)')
plt.savefig(prefix+'_rhoE'+ending)
plt.figure(36)
plt.xlabel('H (kJ/kg)')
plt.ylabel('S (J/g*K)')
plt.savefig(prefix+'_HS'+ending)
plt.figure(37)
plt.xlabel('P (MPa)')
plt.ylabel('V (m3/kg)')
plt.savefig(prefix+'_PV'+ending)
plt.figure(38)
plt.xlabel('rho (kg/m3)')
plt.ylabel('T (K)')
plt.savefig(prefix+'_rhoT'+ending)
pdb.set_trace()
