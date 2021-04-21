"""This model stub allows one to generate simple analytic models and then 
read them into the Analytic3Species model found in model.py which is part of SNEWPY.
"""

from astropy.table import *
import numpy as np

total_energy = (5e52,5e52,2e53)
mean_energy = (15., 15., 15.)
rms_or_pinch = "pinch"
pinch_values = (2.3, 2.3, 2.3)

file_name = "analytic.dat"

table = Table()
table['TIME'] = np.linspace(0,1,2)
table['L_NU_E'] =  np.linspace(1,1,2)*total_energy[0]
table['L_NU_E_BAR'] = np.linspace(1,1,2)*total_energy[1]
table['L_NU_X'] = np.linspace(1,1,2)*total_energy[2]
        
table['E_NU_E'] = np.linspace(1,1,2)*mean_energy[0]
table['E_NU_E_BAR'] = np.linspace(1,1,2)*mean_energy[1]
table['E_NU_X'] = np.linspace(1,1,2)*mean_energy[2]

if rms_or_pinch == "rms":
    table['RMS_NU_E'] = np.linspace(1,1,2)*rms_energy[0]
    table['RMS_NU_E_BAR'] = np.linspace(1,1,2)*rms_energy[1]
    table['RMS_NU_X'] = np.linspace(1,1,2)*rms_energy[2]
    table['ALPHA_NU_E'] = (2.0 * table['E_NU_E'] ** 2 - table['RMS_NU_E'] ** 2) / (
        table['RMS_NU_E'] ** 2 - table['E_NU_E'] ** 2)
    table['ALPHA_NU_E_BAR'] = (2.0 * table['E_NU_E_BAR'] ** 2 - table['RMS_NU_E_BAR'] ** 2) / (
        table['RMS_NU_E_BAR'] ** 2 - table['E_NU_E_BAR'] ** 2)
    table['ALPHA_NU_X'] = (2.0 * table['E_NU_X'] ** 2 - table['RMS_NU_X'] ** 2) / (
        table['RMS_NU_X'] ** 2 - table['E_NU_X'] ** 2)
elif rms_or_pinch == "pinch":
    table['ALPHA_NU_E'] = np.linspace(1,1,2)*pinch_values[0]
    table['ALPHA_NU_E_BAR'] = np.linspace(1,1,2)*pinch_values[1]
    table['ALPHA_NU_X'] = np.linspace(1,1,2)*pinch_values[2]
    table['RMS_NU_E'] = (2.0 + table['ALPHA_NU_E'])/(1.0 + table['ALPHA_NU_E'])*table['E_NU_E']**2
    table['RMS_NU_E_BAR'] =  (2.0 + table['ALPHA_NU_E_BAR'])/(1.0 + table['ALPHA_NU_E_BAR'])*table['E_NU_E_BAR']**2
    table['RMS_NU_X'] = (2.0 + table['ALPHA_NU_X'])/(1.0 + table['ALPHA_NU_X'])*table['E_NU_X']**2 
else:
    print("incorrect second moment method: rms or pinch")
table.write(file_name,format='ascii')
