import NIST_reader as NIST
import numpy as np
from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import PhaseSI
from pdb import set_trace as keyboard
from skimage.transform import ProjectiveTransform
import unicodedata


# def get_dataNIST(x_mid, z_mid):
#     dataNIST=NIST.readNIST(isoType = "isotherm", fluid = 'O2', T=x_mid, P=z_mid/1.0E6, tmin=x_mid, tmax=x_mid, pmin = z_mid/1.0E6, pmax = z_mid/1.0E6, N=1)
#     dataNIST = np.ndarray.tolist(dataNIST)
#     dataNIST = [item for sublist in dataNIST for item in sublist]
#     del dataNIST[:2]
#     return dataNIST
def get_coolprop_TPS(x_mid, z_mid, response, trans, fluid='Oxygen'):
    dataCoolProp =[]

    ####------------------In Temperature-Pressure space ------------####
    if response == "T-P":

        ######------------Re-transforming back to T-P--------##
        [x_mid, z_mid] = np.ndarray.tolist(trans.inverse([x_mid, z_mid])[0])

        data_density = PropsSI('DMASS', 'P', z_mid,'T', x_mid, fluid)
        #data_volume = 0
        data_internal_energy = PropsSI('UMASS', 'P', z_mid,'T', x_mid, fluid)/1000.0
        data_enthalpy= PropsSI('HMASS', 'P', z_mid,'T', x_mid, fluid)/1000.0
        data_entropy= PropsSI('SMASS', 'P', z_mid,'T', x_mid, fluid)/1000.0
        data_cv = PropsSI('CVMASS', 'P', z_mid,'T', x_mid, fluid)/1000.0
        data_cp = PropsSI('CPMASS', 'P', z_mid,'T', x_mid, fluid)/1000.0
        data_a = PropsSI('A', 'P', z_mid,'T', x_mid, fluid)
        #data_joule = 0
        data_mu= PropsSI('VISCOSITY', 'P', z_mid,'T', x_mid, fluid)
        data_k = PropsSI('CONDUCTIVITY', 'P', z_mid,'T', x_mid, fluid)
        data_phase = PropsSI('PHASE', 'P', z_mid,'T', x_mid, fluid)
        #data = [data_density, data_volume, data_internal_energy, data_enthalpy, data_entropy, data_cv, data_cp, data_a, data_joule, data_mu, data_k, data_phase]
        data = [data_density, data_internal_energy, data_enthalpy, data_entropy, data_cv, data_cp, data_a, data_mu, data_k]
        dataCoolProp.append(data)

    ####------------------In Density-Int. energy space---------------####
    elif response == "rho-e":

        [x_mid, z_mid] = np.ndarray.tolist(trans.inverse([x_mid, z_mid])[0])
        data_phase = PhaseSI('UMASS', z_mid,'DMASS', x_mid, fluid)
        data_phase = unicodedata.normalize('NFKD', data_phase).encode('ascii','ignore')
        if data_phase == 'twophase':
            x_mid = 460.5914052903932
            z_mid = 17667.19156915298
        if data_phase == '':
            print "Out of bounds buddy !!"
            keyboard()
 
        data_temp = PropsSI('T', 'UMASS', z_mid,'DMASS', x_mid, fluid)
        data_volume = 0
        data_pressure = PropsSI('P', 'UMASS', z_mid,'DMASS', x_mid, fluid)
        data_enthalpy= PropsSI('HMASS', 'UMASS', z_mid,'DMASS', x_mid, fluid)/1000.0
        data_entropy= PropsSI('SMASS','UMASS', z_mid,'DMASS', x_mid, fluid)/1000.0
        data_cv = PropsSI('CVMASS','UMASS', z_mid,'DMASS', x_mid, fluid)/1000.0
        data_cp = PropsSI('CPMASS', 'UMASS', z_mid,'DMASS', x_mid, fluid)/1000.0
        data_a = PropsSI('A', 'UMASS', z_mid,'DMASS', x_mid, fluid)
        #data_joule = 0
        data_mu= PropsSI('VISCOSITY', 'UMASS', z_mid,'DMASS', x_mid, fluid)
        data_k = PropsSI('CONDUCTIVITY', 'UMASS', z_mid,'DMASS', x_mid, fluid)
        
        #data = [data_density, data_volume, data_internal_energy, data_enthalpy, data_entropy, data_cv, data_cp, data_a, data_joule, data_mu, data_k, data_phase]
        data = [data_temp, data_pressure, data_enthalpy, data_entropy, data_cv, data_cp, data_a, data_mu, data_k]
        dataCoolProp.append(data)

    return dataCoolProp[0]

# def get_dataNIST_I(pointp):
#     point_00 = [pointp[0],pointp[1]]
#     point_10 = [pointp[2],pointp[1]]
#     point_01 = [pointp[0],pointp[3]]
#     point_11 = [pointp[2],pointp[3]]

#     points = [point_00, point_10, point_01, point_11]
#     dataNIST = []
#     for point in points:
#         x_mid = point[0]
#         z_mid = point[1]
#         data = NIST.readNIST(isoType = "isotherm", fluid = 'O2', T=x_mid, P=z_mid/1.0E6, tmin=x_mid, tmax=x_mid, pmin = z_mid/1.0E6, pmax = z_mid/1.0E6, N=1)
#         data= np.ndarray.tolist(data)
#         data= [item for sublist in data for item in sublist]
#         del data[:2]
#         dataNIST.append(data)
#     return dataNIST

def get_coolprop_TP(pointp, response, trans, fluid='Oxygen'):
    point_00 = [pointp[0],pointp[1]]
    point_10 = [pointp[2],pointp[1]]
    point_01 = [pointp[0],pointp[3]]
    point_11 = [pointp[2],pointp[3]]

    points = [point_00, point_10, point_01, point_11]
    dataNIST = []
    dataCoolProp =[]

    ####--------------------In Temperature-Pressure space -----------####
    if response == "T-P":
        for point in points:
            x_mid = point[0]
            z_mid = point[1]
            
            ######------------Re-transforming back to T-P--------##
            x_mid, z_mid = np.ndarray.tolist(trans.inverse([x_mid, z_mid])[0])

            data_density = PropsSI('DMASS', 'P', z_mid,'T', x_mid, fluid)
            #data_volume = 0
            data_internal_energy = PropsSI('UMASS', 'P', z_mid,'T', x_mid, fluid)/1000.0
            data_enthalpy= PropsSI('HMASS', 'P', z_mid,'T', x_mid, fluid)/1000.0
            data_entropy= PropsSI('SMASS', 'P', z_mid,'T', x_mid, fluid)/1000.0
            data_cv = PropsSI('CVMASS', 'P', z_mid,'T', x_mid, fluid)/1000.0
            data_cp = PropsSI('CPMASS', 'P', z_mid,'T', x_mid, fluid)/1000.0
            data_a = PropsSI('A', 'P', z_mid,'T', x_mid, fluid)
            #data_joule = 0
            data_mu= PropsSI('VISCOSITY', 'P', z_mid,'T', x_mid, fluid)
            data_k = PropsSI('CONDUCTIVITY', 'P', z_mid,'T', x_mid, fluid)
            data_phase = PropsSI('PHASE', 'P', z_mid,'T', x_mid, fluid)
            #data = [data_density, data_volume, data_internal_energy, data_enthalpy, data_entropy, data_cv, data_cp, data_a, data_joule, data_mu, data_k, data_phase]
            data = [data_density, data_internal_energy, data_enthalpy, data_entropy, data_cv, data_cp, data_a, data_mu, data_k]
            dataCoolProp.append(data)
            #print dataCoolProp

            #---------------------NIST Database---------------------######
            # data = NIST.readNIST(isoType = "isotherm", fluid = 'O2', T=x_mid, P=z_mid/1.0E6, tmin=x_mid, tmax=x_mid, pmin = z_mid/1.0E6, pmax = z_mid/1.0E6, N=1)
            # data= np.ndarray.tolist(data)
            # data= [item for sublist in data for item in sublist]
            # del data[:2]
            # dataNIST.append(data)
            # print dataNIST

    elif response == "rho-e":
        ####-----------------In Density-Internal Energy space -----------####
        for point in points:
            x_mid = point[0]
            z_mid = point[1]

            [x_mid, z_mid] = np.ndarray.tolist(trans.inverse([x_mid, z_mid])[0])
            data_phase = PhaseSI('UMASS', z_mid,'DMASS', x_mid, fluid)
            data_phase = unicodedata.normalize('NFKD', data_phase).encode('ascii','ignore')
            if data_phase == 'twophase': #bound to critical properties
                x_mid = 460.5914052903932 
                z_mid = 17667.19156915298
            if data_phase == '':
                print "Out of bounds buddy !!"
                keyboard()

            data_temp = PropsSI('T', 'UMASS', z_mid,'DMASS', x_mid, fluid)
            data_volume = 0
            data_pressure = PropsSI('P', 'UMASS', z_mid,'DMASS', x_mid, fluid) #Pa
            data_enthalpy= PropsSI('HMASS', 'UMASS', z_mid,'DMASS', x_mid, fluid)/1000.0 #KJ/kg
            data_entropy= PropsSI('SMASS','UMASS', z_mid,'DMASS', x_mid, fluid)/1000.0 #KJ/kg.k
            data_cv = PropsSI('CVMASS','UMASS', z_mid,'DMASS', x_mid, fluid)/1000.0 #KJ/kg
            data_cp = PropsSI('CPMASS', 'UMASS', z_mid,'DMASS', x_mid, fluid)/1000.0
            data_a = PropsSI('A', 'UMASS', z_mid,'DMASS', x_mid, fluid)
            #data_joule = 0
            data_mu= PropsSI('VISCOSITY', 'UMASS', z_mid,'DMASS', x_mid, fluid)
            data_k = PropsSI('CONDUCTIVITY', 'UMASS', z_mid,'DMASS', x_mid, fluid)
            #data = [data_density, data_volume, data_internal_energy, data_enthalpy, data_entropy, data_cv, data_cp, data_a, data_joule, data_mu, data_k, data_phase]


            data = [data_temp, data_pressure, data_enthalpy, data_entropy, data_cv, data_cp, data_a, data_mu, data_k]
            
            #data = [data_temp]
            dataCoolProp.append(data)

    return dataCoolProp

def check_out_of_bound(x_c, z_c, response, trans, fluid="Oxygen"):

    if response == "T-P":
        [x_mid, z_mid] = np.ndarray.tolist(trans.inverse([x_c, z_c])[0])
        data_phase = PhaseSI('P', z_mid,'T', x_mid, fluid)
        data_phase = unicodedata.normalize('NFKD', data_phase).encode('ascii','ignore')
        if data_phase == '':
            return True
        return False
    if response == "rho-e":
        [x_mid, z_mid] = np.ndarray.tolist(trans.inverse([x_c, z_c])[0])
        data_phase = PhaseSI('UMASS', z_mid,'DMASS', x_mid, fluid)
        data_phase = unicodedata.normalize('NFKD', data_phase).encode('ascii','ignore')
        if data_phase == '':
            return True
        return False
def check_out_of_bound_points(pointp, response, trans, fluid="Oxygen"):
    point_00 = [pointp[0],pointp[1]]
    point_10 = [pointp[2],pointp[1]]
    point_01 = [pointp[0],pointp[3]]
    point_11 = [pointp[2],pointp[3]]

    points = [point_00, point_10, point_01, point_11]
    if response == "T-P":
        for point in points:
            x_mid = point[0]
            z_mid = point[1]
            [x_mid, z_mid] = np.ndarray.tolist(trans.inverse([x_mid, z_mid])[0]) #unit indicies to rho-e or T-P
            data_phase = PhaseSI('UMASS', z_mid,'DMASS', x_mid, fluid)
            data_phase = unicodedata.normalize('NFKD', data_phase).encode('ascii','ignore')
            if data_phase == '':
                return True
        return False
    if response == "rho-e":
        for point in points:
            x_mid = point[0]
            z_mid = point[1]
            [x_mid, z_mid] = np.ndarray.tolist(trans.inverse([x_mid, z_mid])[0]) #unit indicies to rho-e or T-P
            data_phase = PhaseSI('UMASS', z_mid,'DMASS', x_mid, fluid)
            data_phase = unicodedata.normalize('NFKD', data_phase).encode('ascii','ignore')
            if data_phase == '':
                return True
        return False    