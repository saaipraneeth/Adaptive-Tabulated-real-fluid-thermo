import numpy as np
from CoolProp.CoolProp import PropsSI
import pdb
import matplotlib.pyplot as plt
import quad_utilities as utility

Pressure = np.linspace(0.01, 22, 100)
x_mid = 147
error_rho = [None]*100
error_enthalpy =[None]*100
error_entropy=[None]*100
error_cv=[None]*100
error_int_e = [None]*100
error_cp = [None]*100
error_a = [None]*100
error_viscosity =[None]*100
error_cond = [None]*100
for i, item in enumerate(Pressure):
	z_mid = item*1E6
	dataNIST = utility.get_dataNIST(x_mid, z_mid)
	coolprop = utility.get_coolprop_TPS(x_mid, z_mid)
	error_rho[i] = (abs(dataNIST[0] - coolprop[0])/dataNIST[0])*100
	error_int_e[i] = (abs(dataNIST[2] - coolprop[1])/dataNIST[2])*100
	error_enthalpy[i] = (abs(dataNIST[3] - coolprop[2])/dataNIST[3])*100
	error_entropy[i] = (abs(dataNIST[4] - coolprop[3])/dataNIST[4])*100
	error_cv[i] = (abs(dataNIST[5] - coolprop[4])/dataNIST[5])*100
	error_cp[i] = (abs(dataNIST[6] - coolprop[5])/dataNIST[6])*100
	error_a[i] = (abs(dataNIST[7] - coolprop[6])/dataNIST[7])*100
	error_viscosity[i] = (abs(dataNIST[9] - coolprop[7])/dataNIST[9])*100
	error_cond[i] = (abs(dataNIST[10] - coolprop[8])/dataNIST[10])*100

#pdb.set_trace()
plt.figure()
currentAxis = plt.gca()
currentAxis.set_xlim([0.01,22])
currentAxis.set_ylim([0,15])
plt.plot(Pressure,error_rho,'g--', label="Density")
plt.plot(Pressure,error_int_e,'b--', label="Internal Energy")
plt.plot(Pressure,error_enthalpy,'g--', label="Enthalpy")
plt.plot(Pressure,error_a,'g*', label="Speed of Sound")
plt.plot(Pressure,error_viscosity,'r--', label="Viscosity")
plt.plot(Pressure,error_cond,'bs', label="Therm. Conductivity")
plt.xlabel('Pressure[MPa]')
plt.ylabel('Error [%]')
plt.title('CoolProp Error(w.r.t NIST) vs Pressure at 147 K')

plt.legend()
plt.show()