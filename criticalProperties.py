import numpy as np
import pdb


#---------------------------------------------------------------------
#   Generic solution vector
#---------------------------------------------------------------------
class solutionVector():
    def setThermodynamics(self,fluid):
         self.Pcrit,self.Tcrit,self.rhocrit,self.omega,self.dipole = readCritPropDatabase(fluid)
         self.MW=getMW(fluid)
         self.Rcst = 8.3144621	#J/(mol*K)
         self.vcrit   = (self.MW/self.rhocrit)    #m^3/mol
         self.Zcrit   = (self.Pcrit*self.vcrit/(self.Rcst*self.Tcrit))   # Zc=PcVc/RTc


def getMW(fluid):
    MW={'H2O':18.0153,'H2':2.01588,'O2':31.9988,'N2':28.0134}
    return MW[fluid]

def readCritPropDatabase(fluid,dbPath='./critProp.dat'):
    lines = [line.strip() for line in open(dbPath)]
    for line in lines:
        if line.split()[0]==fluid:
            str_critProp=line.split()[2:]
            return [float(val) for val in str_critProp]
            #pCrit[Pa] TCrit[K] rhoCrit[kg/m3] omega(mu) dipole at NBP [debye]





