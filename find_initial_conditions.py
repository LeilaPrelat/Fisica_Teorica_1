#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 22:27:14 2020

@author: leila
"""

from scipy import special #funciones de Bessel
import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
from numpy.linalg import det
from scipy.optimize import minimize

#%%
lista_colores = ['coral','yellowgreen','midnightblue','purple','darkred','orange','aquamarine','hotpink','steelblue','darkgreen'] 

tamfig = (10,8)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%
im_kz = 0
modo = 1
R1 = 70*1e-9
R2 = 85*1e-9

ε_path = r'/home/leila/Desktop/Teo1_definitivo2/dielectric_functions'
os.chdir(ε_path)

det_path = r'/home/leila/Desktop/Teo1_definitivo2' #funcion determinante

#%%

try:
    sys.path.insert(1, ε_path)
    from dielectric_functions import ε_functions
except ModuleNotFoundError as error:
    print(error)
    ε_path = input('path de la carpeta dielectric_functions')
    sys.path.insert(1, ε_path)

ε_SiO2 = ε_functions.ε_SiO2
ε_Ag = ε_functions.ε_Ag
ε_CdS = ε_functions.ε_CdS

try: 
    ε_functions.load_data(ε_path)
    plt.close('all')

except OSError or IOError as error:
    print(error)
    ε_path2 = input('path de los archivos .txt de la carpeta dielectric_functions')
    ε_functions.load_data(ε_path2)
    plt.close('all')
    
try:
    sys.path.insert(1, det_path)
    from determinante import det_M3medios as det_M
except ModuleNotFoundError as error:
    print(error)
    det_path = input('path de det_function.py')
    sys.path.insert(1, det_path)  


ε1,ε2,ε3 = ε_CdS,ε_SiO2,ε_Ag

#%%

if ε1==ε_Ag:
    title1 = 'Ag'
elif ε1==ε_CdS:
    title1 = 'CdS'
elif ε1==ε_SiO2:
    title1 = 'SiO2'
    
if ε2==ε_Ag:
    title2 = 'Ag'
elif ε2==ε_CdS:
    title2 = 'CdS'
elif ε2==ε_SiO2:
    title2 = 'SiO2' 

if ε3==ε_Ag:
    title3 = 'Ag'
elif ε3==ε_CdS:
    title3 = 'CdS'
elif ε3==ε_SiO2:
    title3 = 'SiO2'

title = title1 + '/' + title2 + '/' + title3
    
#%%
print('Minimizacion con Nelder-Mead del determinante hallando Re(kz) y Im(kz)')

Ne = int(1e3)        
Elist_NM = np.linspace(0.8,3.8,Ne)

def det_M_2variable(x):
    Energy = 2.5
    Re_kz = x[0]
    Im_kz = x[1]
    dett = det_M([Re_kz,Im_kz],Energy,modo,R1,R2,ε1,ε2,ε3)
    return np.abs(dett)

x0_0 = [0.8*1e8,0.1*1e8,0.2*1e8]
x0_1 = [1e2,1e3,1e4,1e5,1e6,1e7,1e8,1]
#
#x0_0 = [1e7]
#x0_1 = [1e2,1e3,1e4,1e5]

for j in range(len(x0_0)):  
    for k in range(len(x0_1)):
        
        init_conditions = np.array([x0_0[j],x0_1[k]])
        res = minimize(det_M_2variable, init_conditions, method='Nelder-Mead')

        if res.x[1]>1 and det_M_2variable([res.x[0],res.x[1]])<0.1: 
            print(init_conditions) 
            print(res.message)
            print(res.x[0],',',res.x[1])
            print(det_M_2variable([res.x[0],res.x[1]]))   
            print('')
  