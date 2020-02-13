#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 13:33:49 2020

@author: leila
"""

from scipy import special #funciones de Bessel
import sys
import numpy as np
import os 
import matplotlib.pyplot as plt
from numpy import linalg as LA

#%%
lista_colores = ['coral','yellowgreen','midnightblue','purple','darkred','orange','aquamarine','hotpink','steelblue','darkgreen'] 

tamfig = (10,8)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

R = 85*1e-9

ε_path = r'/home/leila/Desktop/Teo1_definitivo2/dielectric_functions'
os.chdir(ε_path)

det_path = r'/home/leila/Desktop/Teo1_definitivo2'

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
    ε_path2 = input('path de la carpeta dielectric_functions')
    ε_functions.load_data(ε_path2)
    plt.close('all')

ε2,ε3 = ε_SiO2,ε_Ag

#%%
    
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

title = title2 + '/' + title3

#%%

print('Importar las relaciones de dispersion y graficar')

R2 = 85*1e-9
list_R = [85]
list_R = np.array(list_R)*1e-9

list_modos = [0,1]

print('Graficar E vs Re(kz)')
color = 0
plt.figure(figsize=tamfig)
for modo in list_modos:
    for R in list_R:
        color = color +1
    
        folder_mediums = title2 + '_' + title3
        folder_R = 'R' + '_' + str(int(R*1e9)) 
        os.chdir(det_path + '/' + folder_mediums + '/' + folder_R)
        
        try:   
            Elist_NM = np.loadtxt('NM_mod' + str(modo) + '_E_list.txt',delimiter='\t')
            kz_real_min_NM = np.loadtxt('NM_mod' + str(modo) + '_Rekz_list.txt',delimiter='\t')   
        except OSError or IOError as error:
            print(error)
            direc = input('path de los .txt de las relaciones de dispersion')
            os.chdir(direc)
            Elist_NM = np.loadtxt('NM_mod' + str(modo) + '_E_list.txt',delimiter='\t')
            kz_real_min_NM = np.loadtxt('NM_mod' + str(modo) + '_Rekz_list.txt',delimiter='\t')   
            
        kz_real_min_NM2 = np.array(kz_real_min_NM)*1e-8   
        
        plt.title(title +', R = %i nm' %(R*1e9),fontsize=tamtitle)
        plt.plot(np.abs(kz_real_min_NM2),Elist_NM,'.',ms=10,color=lista_colores[color],alpha=0.7,label='modo = %i, R = %i nm' %(modo,int(R*1e9))) 
        plt.xlabel('Re($k_z$) [x 10$^8$ m$^{-1}$]',fontsize=tamletra)
        plt.ylabel('Energy [eV]',fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
        plt.grid(1)

print('Graficar E vs Im(kz)')            
color = 0
plt.figure(figsize=tamfig)
for modo in list_modos:
    for R in list_R:
        color = color +1
    
        folder_mediums = title2 + '_' + title3
        folder_R = 'R' + '_' + str(int(R*1e9)) 
        os.chdir(det_path + '/' + folder_mediums + '/' + folder_R)
        
        try:   
            Elist_NM = np.loadtxt('NM_mod' + str(modo) + '_E_list.txt',delimiter='\t')
            kz_imag_min_NM = np.loadtxt('NM_mod' + str(modo) + '_Imkz_list.txt',delimiter='\t')     
        except OSError or IOError as error:
            print(error)
            direc = input('path de los .txt de las relaciones de dispersion')
            os.chdir(direc)
            Elist_NM = np.loadtxt('NM_mod' + str(modo) + '_E_list.txt',delimiter='\t')
            kz_imag_min_NM = np.loadtxt('NM_mod' + str(modo) + '_Imkz_list.txt',delimiter='\t')  
            
        kz_imag_min_NM2 = np.array(kz_imag_min_NM)*1e-8   
        
        plt.title(title +', R = %i nm' %(R*1e9),fontsize=tamtitle)
        plt.plot(np.abs(kz_imag_min_NM2),Elist_NM,'.',ms=10,color=lista_colores[color],alpha=0.7,label='modo = %i, R = %i nm' %(modo,int(R*1e9))) 
        plt.xlabel('Im($k_z$) [x 10$^8$ m$^{-1}$]',fontsize=tamletra)
        plt.ylabel('Energy [eV]',fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
        plt.grid(1)            


#%%
    
    
