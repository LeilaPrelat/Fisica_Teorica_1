# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 22:51:59 2019

@author: Usuario
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
    
print('Observar el comportamiento de la funcion det_M fijando im_kz = 0 y haciendo un graf por cada E')

n = 7
x = np.logspace(6+np.log10(5),7+np.log10(4),600)
Elist_NM = np.linspace(0.8,3.8,n)
plt.figure(figsize=tamfig)
for Energy in Elist_NM:
    
    def det_M_1variable(Re_kz):
        dett = det_M([Re_kz,im_kz],Energy,modo,R1,R2,ε1,ε2,ε3)
        return np.abs(dett)
    
    detlist=[]
    for real_kz in x:
        detlist.append(det_M_1variable(real_kz))
        
    index_min = int(np.argmin(detlist))    
    print(Energy,x[index_min])
    
    plt.title(title +', modo = %i, Im($k_z$) = %.2e, R1 = %i nm, R2 = %i nm' %(modo,im_kz,R1*1e9,R2*1e9),fontsize=tamtitle)
    plt.plot(x,detlist,'.',ms=10,alpha=0.7,label='E=%.2e eV' %(Energy))
    plt.ylabel('|det(L)|',fontsize=tamletra)
#    plt.ylabel('min{sing. values}',fontsize=tamletra)
    plt.xlabel('Re($k_z$) [$m^{-1}$]',fontsize=tamletra)
    plt.xscale('log')
    plt.yscale('log')
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='upper right',markerscale=2,fontsize=tamlegend)  
    plt.grid(1)
#    plt.tight_layout(1)
    
#%%
print('Minimizacion con Nelder-Mead del determinante hallando Re(kz) y Im(kz)')

Ne = int(1e3)        
Elist_NM = np.linspace(2,3.8,Ne)

det_M_NM = []
kz_real_min_NM = []
kz_imag_min_NM = []

if (ε1,ε2,ε3) == (ε_Ag,ε_SiO2,ε_Ag):
    if modo==0:
        if R1==5*1e-9:
            cond_inicial = np.array([24638187.70156102,350640.09316422])
        
        elif R1==40*1e-9:
            cond_inicial = np.array([1e7,1e5])
        
        elif R1==70*1e-9:
            cond_inicial=np.array([11971949.74033338,68813.08756318])
    #cond_inicial = np.array([1.17*1e7,66240])
            
    elif modo==1:
        if R1==5*1e-9:
            cond_inicial = np.array([203557.97975328687 , -199927851.60810703])
        
        elif R1==40*1e-9:
            cond_inicial = np.array([1.80898938e+06,2.12055049e+08])
    
        elif R1==70*1e-9:
            cond_inicial = np.array([1.80898938e+06,2.12055049e+08])
            

if (ε1,ε2,ε3) == (ε_CdS,ε_SiO2,ε_Ag):
    if modo==0:       
        if R1==30*1e-9:
#            cond_inicial = np.array([8401429.141698383,75501.33358144428])
            cond_inicial = np.array([30691207.23202795,4701587.9489140995])
        
        elif R1==70*1e-9:
            cond_inicial=np.array([8401429.14167578,75501.33357499144])
    #cond_inicial = np.array([1.17*1e7,66240])
            
    elif modo==1:        
        if R1==30*1e-9:
            cond_inicial = np.array([1.80898938e+06,2.12055049e+08])
    
        elif R1==70*1e-9:
            cond_inicial = np.array([8401429.141728252,75501.33361986515])

for Energy in Elist_NM:  
    
    def det_M_2variable(x):
        Re_kz = x[0]
        Im_kz = x[1]
        dett = det_M([Re_kz,Im_kz],Energy,modo,R1,R2,ε1,ε2,ε3)
        return np.abs(dett)
    
    res = minimize(det_M_2variable, cond_inicial, method='Nelder-Mead')
    print(res.message)
    print(res.x[0],',',res.x[1])
    kz_real_min_NM.append(res.x[0])
    kz_imag_min_NM.append(res.x[1])
    det_M_NM.append(det_M_2variable([res.x[0],res.x[1]]))
    cond_inicial = np.array([res.x[0],res.x[1]])

#%%
   
print('Guardar la relacion de dispersión obtenida')
    
os.chdir(det_path)
folder_mediums = title1 + '_' + title2 + '_' + title3

try:
    os.mkdir(folder_mediums)
except FileExistsError as error: 
    print(error)

folder_R1 = 'R1' + '_' + str(int(R1*1e9))   
try:
    os.chdir(det_path + '/' + folder_mediums)
    os.mkdir(folder_R1)
except FileExistsError as error: 
    print(error)
 
os.chdir(det_path + '/' + folder_mediums + '/' + folder_R1)
np.savetxt('NM_mod' + str(modo) + '_E_list.txt',[Elist_NM],delimiter='\t')
np.savetxt('NM_mod' + str(modo) + '_Rekz_list.txt',[kz_real_min_NM],delimiter='\t')   
np.savetxt('NM_mod' + str(modo) + '_Imkz_list.txt',[kz_imag_min_NM],delimiter='\t')   

#%%

kz_real_min_NM2 = np.array(kz_real_min_NM)*1e-8    

#%%
max_kz_real = int(np.argmax(kz_real_min_NM2))
max_energy = Elist_NM[max_kz_real]

marca_x = np.linspace(np.min(kz_real_min_NM2),np.max(kz_real_min_NM2),100)
marca_y = np.ones(100)*max_energy

print('Graficar la minimizacion con Nelder-Mead')

plt.figure(figsize=tamfig)     
plt.title(title +', modo = %i, R1 = %i nm, R2 = %i nm' %(modo,R1*1e9,R2*1e9),fontsize=tamtitle)
plt.plot(np.abs(kz_real_min_NM2),Elist_NM,'.',ms=10,color=lista_colores[7],alpha=0.7,label='Nelder-Mead') 
plt.plot(marca_x,marca_y,'.',ms=10,color=lista_colores[5],alpha=0.7,label='maximo %.2f eV' %(max_energy)) 
plt.xlabel('Re($k_z$) [x 10$^8$ m$^{-1}$]',fontsize=tamletra)
plt.ylabel('Energy [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
plt.grid(1)
#plt.tight_layout(1)

#%%