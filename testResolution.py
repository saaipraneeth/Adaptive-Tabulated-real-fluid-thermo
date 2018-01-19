########################################################################################################
#  A Program to read the NIST website and extract isothermal or isobaric data
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


#User defined quantities
fluid='O2'
isoType='isotherm'
Pmin=0.5E6
Pmax=1.5E6
Tmin=110.0
Tmax=130.0
NbT=40
NbP=40
T=444
P=0.8
NE=20
NRho=20
NVar=14
ending='.svg'

#Defines the start values of the arrays
startValues = [-1] * (NbT+1)  

dataIsotherm=np.zeros([NVar,NbT*NbP*10])
RHOarray_maxRho=np.zeros([NbT])
RHOarray_minRho=np.zeros([NbT])
Earray_minRho=np.zeros([NbT])
Earray_maxRho=np.zeros([NbT])
for ii,T in enumerate(np.linspace(Tmin,Tmax,NbT)):
    dataNIST=NIST.readNIST(isoType, fluid, T, P, Tmin,Tmax,Pmin/1.0E6,Pmax/1.0E6,NbP)
    RHOarray_minRho[ii]=dataNIST[NIST.colNIST('rho'),0]
    RHOarray_maxRho[ii]=dataNIST[NIST.colNIST('rho'),-1]
    Earray_minRho[ii]=dataNIST[NIST.colNIST('E'),0]
    Earray_maxRho[ii]=dataNIST[NIST.colNIST('E'),-1]    
    if ii==0: 
#        RHOarray_minT=dataNIST[NIST.colNIST('rho'),:]
        Tarray=dataNIST[NIST.colNIST('T'),:]
        Parray=dataNIST[NIST.colNIST('P'),:]
        Harray=dataNIST[NIST.colNIST('H'),:]
        RHOarray=dataNIST[NIST.colNIST('rho'),:]
        Nadded=np.size(Tarray)
        Earray=dataNIST[NIST.colNIST('E'),:]
        Earray_minT=dataNIST[NIST.colNIST('E'),:]
        nPts=np.size(dataNIST[0,:])
        dataIsotherm[:,0:nPts]=dataNIST[:,:]
        Rhoarray_minT=dataNIST[NIST.colNIST('rho'),:]
        startValues[0]=0
        startValues[1]=Nadded
    else:     
        Nadded=np.size(dataNIST[0,:])
        Rhoarray_maxT=dataNIST[NIST.colNIST('rho'),:]
        Tarray=np.append(Tarray,dataNIST[NIST.colNIST('T'),:])
        Parray=np.append(Parray,dataNIST[NIST.colNIST('P'),:])
        Harray=np.append(Harray,dataNIST[NIST.colNIST('H'),:])
        Earray=np.append(Harray,dataNIST[NIST.colNIST('E'),:])
        RHOarray=np.append(RHOarray,dataNIST[NIST.colNIST('rho'),:])  
        startValues[ii+1]=startValues[ii]+Nadded
        dataIsotherm[:,startValues[ii]:startValues[ii+1]]=dataNIST[:,:]
    if ii==NbT-1:   Earray_maxT=dataNIST[NIST.colNIST('E'),:]

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

pdb.set_trace()
#Algorithm
###########################################
# Get homogeneous energy divisions
E_min=min(Earray_maxRho)*0.9
E_max=max(Earray_minRho)*0.9
RHO_minmax=np.zeros([NE,2])
dataHomo=np.zeros([NE,NRho,NVar])
dataTemp=np.zeros([NE,NbT,NVar])

Ehomo=np.linspace(E_min,E_max,NE)
# Set start and end of the density for each energy level.

for ii,Elocal in enumerate(Ehomo):    
    #  Start: Minimum density
    if Elocal >= min(Earray_minRho) and Elocal <= max(Earray_minRho):        
        weights,pt = interpolationWeights(Earray_minRho, Elocal)
        k0=startValues[pt[0]]
        k1=startValues[pt[1]]
        for kk in range(NVar):
            if kk==2: print ii,dataIsotherm[kk,k0] + (dataIsotherm[kk,k1]-dataIsotherm[kk,k0])*weights
            dataHomo[ii,0,kk]=dataIsotherm[kk,k0] + (dataIsotherm[kk,k1]-dataIsotherm[kk,k0])*weights
    elif Elocal < min(Earray_minRho):
        weights,pt = interpolationWeights(Earray_minT, Elocal)
        k0=startValues[0]+pt[0]
        k1=startValues[0]+pt[1]      
        for kk in range(NVar):
            if kk==2: print '2',ii,dataIsotherm[kk,k0] + (dataIsotherm[kk,k1]-dataIsotherm[kk,k0])*weights
            dataHomo[ii,0,kk]=dataIsotherm[kk,k0] + (dataIsotherm[kk,k1]-dataIsotherm[kk,k0])*weights

    # End: Maximum density
    if Elocal >= min(Earray_maxT) and Elocal <= max(Earray_maxT):
        weights,pt = interpolationWeights(Earray_maxT, Elocal)
        k0=startValues[NbT-1]+pt[0]
        k1=startValues[NbT-1]+pt[1]  
        for kk in range(NVar):
            dataHomo[ii,-1,kk]=dataIsotherm[kk,k0] + (dataIsotherm[kk,k1]-dataIsotherm[kk,k0])*weights
    elif Elocal > min(Earray_maxRho):
        weights,pt = interpolationWeights(Earray_maxRho, Elocal)
        k0=startValues[pt[0]+1]-1
        k1=startValues[pt[1]+1]-1
        for kk in range(NVar):
            dataHomo[ii,-1,kk]=dataIsotherm[kk,k0] + (dataIsotherm[kk,k1]-dataIsotherm[kk,k0])*weights

    #Interpolate isotherm on the isoenergy grid
    for jj in range(NbT):    
        E1 = dataIsotherm[NIST.colNIST('E'),startValues[jj]]
        E2 = dataIsotherm[NIST.colNIST('E'),startValues[jj+1]-1]
        if (Elocal >= E1 and Elocal <= E2) or (Elocal <= E1 and Elocal >= E2):            
            Etemp = dataIsotherm[NIST.colNIST('E'),startValues[jj]:startValues[jj+1]-1 ]
            weights,pt = interpolationWeights(Etemp, Elocal)            
            k0=startValues[jj]+pt[0]
            k1=startValues[jj]+pt[1]
            for kk in range(NVar):
                dataTemp[ii,jj,kk]=dataIsotherm[kk,k0] + (dataIsotherm[kk,k1]-dataIsotherm[kk,k0])*weights

#Interpolate isoenergy on the homogeneously spaced density field
for ii,Elocal in enumerate(Ehomo):
    rho_division= np.linspace(dataHomo[ii,0,NIST.colNIST('rho')],dataHomo[ii,-1,NIST.colNIST('rho')],NRho)
    allRho=dataTemp[ii,:,NIST.colNIST('rho')]
    allRho_remove0 = np.ma.masked_equal(allRho,0)
    RHOtemp = np.append(np.append(dataHomo[ii,0,NIST.colNIST('rho')],allRho_remove0.compressed()),dataHomo[ii,-1,NIST.colNIST('rho')])
    for jj in range(1,NRho-1):
        weights,pt = interpolationWeights(RHOtemp, rho_division[jj])
        print 'Energy: ',Elocal,'pt: ',pt,weights,jj,RHOtemp,rho_division[jj]
        for kk in range(NVar):
            k0=pt[0]
            k1=pt[1]
            allVar=dataTemp[ii,:,kk]
            allVar_remove0 = np.ma.masked_equal(allVar,0)
            temp = np.append(np.append(dataHomo[ii,0,kk],allVar_remove0.compressed()),dataHomo[ii,-1,kk])
            dataHomo[ii,jj,kk]=temp[k0]+(temp[k1]-temp[k0])*weights

    plt.figure(35)
    plt.scatter(dataHomo[ii,:,2],dataHomo[ii,:,4],color='k')
    plt.scatter(dataHomo[ii,0,2],dataHomo[ii,0,4])
    plt.scatter(dataHomo[ii,-1,2],dataHomo[ii,-1,4],color='r')
    plt.savefig('isotherm_rhoE'+ending)
    plt.figure(36)
    plt.plot(dataHomo[ii,:,NIST.colNIST('H')],dataHomo[ii,:,NIST.colNIST('S')],color='k')
    # Divide the density field in N-1 segments
    
fig=plt.figure(33)
plt.xlabel('rho (kg/m3)')
plt.ylabel('P (MPa)')
plt.savefig('isotherm_rhoP'+ending)
plt.figure(34)
plt.xlabel('E (kJ/kg)')
plt.ylabel('P (MPa)')
plt.savefig('isotherm_eP'+ending)
plt.figure(35)
plt.xlabel('rho (kg/m3)')
plt.ylabel('E (kJ/kg)')
plt.savefig('isotherm_rhoE'+ending)
plt.figure(36)
plt.xlabel('H (kJ/kg)')
plt.ylabel('S (J/g*K)')
plt.savefig('isotherm_HS'+ending)
plt.figure(37)
plt.xlabel('P (MPa)')
plt.ylabel('V (m3/kg)')
plt.savefig('isotherm_PV'+ending)
plt.figure(38)
plt.xlabel('rho (kg/m3)')
plt.ylabel('T (K)')
plt.savefig('isotherm_rhoT'+ending)
pdb.set_trace()




# Write a tecplot file
T=np.reshape(dataHomo[:,:,0],np.size(dataHomo[:,:,0]))
P=np.reshape(dataHomo[:,:,1],np.size(dataHomo[:,:,1]))
RHO=np.reshape(dataHomo[:,:,2],np.size(dataHomo[:,:,2]))
S=np.reshape(dataHomo[:,:,6],np.size(dataHomo[:,:,6]))
phase=np.reshape(dataHomo[:,:,13],np.size(dataHomo[:,:,13]))
# 'rho':2, 'V':3,'E':4,'H':5,'S':6,'Cv':7, 'Cp':8, 'C':9, 'JT':10,'mu':11, 'kappa':12, 'phase':13}
var=(('P',P),('temp',T),('RHO',RHO),('S',S),('phase',phase))
tec.tecplot_WriteRectilinearMesh('tablein_rhoE.dat',rho_division,Ehomo,[],var)


# Write a netcdf file for TAU
#_di = sp.dtype('int32')
#_dd = sp.dtype('float64')
#f = netcdf.netcdf_file('cavitationTable.nc', 'w')
#f.type='Cavitation full thermodynamic description'
#f.NE = NE
#f.NRho= NRho
#f.Ntot= NE*NRho
#f.NVar=NVar
