#####################################################################
#!/usr/bin/env python
# Description: Reads the NIST website (http://webbook.nist.gov/) and
#              extracts the data
# Date       : Mai 2014
__author__ = "jphickey"
__version__ = "0.1"
####################################################################
import numpy as np
import pdb
import requests
#from bs4 import BeautifulSoup

#----------------------------------------------------
# Wrapper function to read the NIST website
#----------------------------------------------------
def readNIST(isoType, fluid, T, P, tmin, tmax, pmin, pmax,N):
    if isoType=='isotherm':
        website=getNISTWebsite_isoTherm(fluid, pmin, pmax, T, N)
    elif isoType=='isobar':
        website=getNISTWebsite_isoBar(fluid, tmin, tmax, P, N)
    else:
        quit('Only isotherm and isobar are implemented')
    return readNISTWebsite(website)


# Sets the column for variables from the NIST website
def colNIST(var):
    colNumber={'T':0, 'P':1, 'rho':2, 'V':3,'E':4,'H':5,'S':6,'Cv':7, 'Cp':8, 'C':9, 'JT':10,'mu':11, 'kappa':12, 'phase':13}
    return colNumber[var]

# Defines the website name for isothermal
def getNISTWebsite_isoTherm(fluid, Plow,Phigh, T, Nb=100,units='kg'):
    if units=='kg':
        unitSuffix='&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm'
    elif units=='mol':
        unitSuffix='&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=mol%2Fl&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm'
    website='http://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID='\
            + fluidNIST(fluid)\
            + '&Type=' + typeNIST() \
            + '&Digits=' + digitNIST() \
            + '&PLow=' + str(Plow)\
            + '&PHigh=' + str(Phigh) \
            + '&PInc='+ str((Phigh-Plow)/float(Nb))\
            + '&T='+str(T) \
            +  unitSuffix
    #print website
    return website

# Defines the website name for isobaric
def getNISTWebsite_isoBar(fluid, Tlow, Thigh, P, Nb=100,units='kg'):
    if units=='kg':
        unitSuffix='RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm'
    elif units=='mol':
        unitSuffix='&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=mol%2Fl&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm'
    website='http://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID='\
            + fluidNIST(fluid)\
            + '&Type=' + typeNIST('IsoBar') \
            + '&Digits=' + digitNIST() \
            + '&P=' + str(P)\
            + '&THigh=' + str(Thigh) \
            + '&TLow=' + str(Tlow) \
            + '&TInc='+ str((Thigh-Tlow)/float(Nb))\
            +  unitSuffix
    return website
#http://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C7782447&Type=IsoTherm&Digits=5&PLow=1&PHigh=2&PInc=&T=300&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm

# NIST uses codes to define the fluids
def fluidNIST(fluid):
    if   fluid=='H2O': return 'C7732185'
    elif fluid=='CO' : return 'C630080'
    elif fluid=='H2' : return 'C1333740'
    elif fluid=='O2' : return 'C7782447'
    elif fluid=='N2' : return 'C7727379'
    else: quit('no fluid found')

def typeNIST(isotype='IsoTherm'):
    if isotype.lower()=='isotherm': return isotype
    elif isotype.lower()=='isobar': return isotype

def digitNIST():
    return '5'

# Read the NIST website
def readNISTWebsite(website):
    session = requests.session()
    doc = session.get(website)
    docstring=doc.content
    return getDatafromString(docstring.split('\n'))

# Converts the state of the fluid to a numeral
def getDatafromString(doc):
    nbVar = 14
    nbData=np.size(doc)-2  #First and last line are skipped
    dataArray=np.zeros([nbVar,nbData])
    removeRow=[]
    for ii in range(nbData):
        doc[ii+1]=doc[ii+1].replace('liquid','1')
        doc[ii+1]=doc[ii+1].replace('vapor','2')
        doc[ii+1]=doc[ii+1].replace('solid','3')
        doc[ii+1]=doc[ii+1].replace('supercritical','4')

        temp=doc[ii+1].split()
        for jj in range(nbVar):
            if temp[jj]=='undefined':
               if ii  not in removeRow: removeRow.append(ii)
            else:
                if jj==colNIST('P'):
                   coef=1.0E6
                else:
                   coef=1.0
                dataArray[jj,ii]=float(temp[jj])*coef
    if removeRow!=[]:
        dataArray = np.delete(dataArray,(removeRow), axis=1)
    return dataArray
