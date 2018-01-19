########################################################################################################
#  A Program to read the NIST website and extract isothermal or isobaric data
########################################################################################################
import numpy as np
from scipy.io import netcdf

def check_io_status(io):
    if io=='w': return True
    elif io=='r': return True
    else: return False

def python_stream_open(fname,io='w'):
    if check_io_status(io): 
        return netcdf.netcdf_file(fname, io)
    else:
        quit('The selected io status is not recognized')

def python_stream_att(fhandle,att,value):
    fhandle.att=value

def python_stream_variable(fhandle,name,dtype)
    time = f.createVariable('time', 'd', ('time',))
    np.
