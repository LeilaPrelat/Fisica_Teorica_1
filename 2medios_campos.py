#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 23:18:27 2020

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

π = np.pi
μ0 = 1
c = 2.99792458*1e8 #m/s (SI units)
h = 4.135667731*1e-15 #eV/s
ħ = h/(2*π)
    
#%%

R = 85*1e-9
modo = 0

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
    ε_path2 = input('path de los archivos .txt de la carpeta dielectric_functions')
    ε_functions.load_data(ε_path2)
    plt.close('all')

ε2,ε3 = ε_Ag,ε_SiO2,ε_Ag

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

print('Importar la relacion de dispersion')

folder_mediums = title2 + '_' + title3
folder_R = 'R' + '_' + str(int(R*1e9))   
os.chdir(det_path + '/' + folder_mediums + '/' + folder_R)

try:   
    Elist_NM = np.loadtxt('NM_mod' + str(modo) + '_E_list.txt',delimiter='\t')
    kz_real_min_NM = np.loadtxt('NM_mod' + str(modo) + '_Rekz_list.txt',delimiter='\t')   
    kz_imag_min_NM = np.loadtxt('NM_mod' + str(modo) + '_Imkz_list.txt',delimiter='\t')     
except OSError or IOError as error:
    print(error)
    direc = input('path de los .txt de las relaciones de dispersion')
    os.chdir(direc)
    Elist_NM = np.loadtxt('NM_mod' + str(modo) + '_E_list.txt',delimiter='\t')
    kz_real_min_NM = np.loadtxt('NM_mod' + str(modo) + '_Rekz_list.txt',delimiter='\t')   
    kz_imag_min_NM = np.loadtxt('NM_mod' + str(modo) + '_Imkz_list.txt',delimiter='\t')  

#%%

print('Se definen los coeficientes del determinante de 4x4')

def coef(indice,R,ε2,ε3):
    
    E = Elist_NM[indice]
    kz_real = kz_real_min_NM[indice]
    kz_imag = kz_imag_min_NM[indice]
    
    kz = kz_real+1j*kz_imag       
        
    ω = (E/ħ)
    k0 = ω/c
    xz_real = kz.real/k0
    xz_imag = kz.imag/k0
    xz = xz_real+1j*xz_imag    

    def epsi(medio):
        if medio==2:
            return ε2
        if medio==3:
            return ε3
        
    def xt2(medio):
        inside = xz**2-μ0*epsi(medio)
        return inside
            
    def xt(medio):
        inside = xt2(medio)+0j
        xtmas=(inside)**(1/2)
        if xtmas.real>=0:
            xtfin=xtmas
        else:
            xtfin=-xtmas       
        return xtfin
    
    def a(j):
        if modo!=0:
            rta1 = xz*modo*xt2(j)
            rta2 = Rbarra
        else: 
            rta1 = 0
            rta2 = 1
        return rta1/rta2
    
    def b(j,k):
        return 1j*μ0*xt(j)*xt2(k)

    def d(j,k):
        return 1j*epsi(j)*xt(j)*xt2(k)
    
    Rbarra = R*k0
    
    def I(k):
        return special.iv(modo,xt(k)*Rbarra)
    
    def derI(n):
        return special.ivp(modo,xt(n)*Rbarra)
    
    def K(n):
        return special.kv(modo,xt(n)*Rbarra)
    
    def derK(n):
        return special.kvp(modo,xt(n)*Rbarra)
    
    A = [I(2),0,-K(3),0]
    B = [0,I(2),0,-K(3)]
    C = [a(3)*I(2),b(2,3)*derI(2),-a(2)*K(3),-b(3,2)*derK(3)]
    D = [d(2,3)*derI(2),-a(3)*I(2),-d(3,2)*derK(3),a(2)*K(3)]
    
    M = np.matrix([A,B,C,D], dtype='complex')
    try:
        λ, v = LA.eig(M)
        mod_λ = []
        for autovalor in λ:
            mod_λ.append(np.abs(autovalor))
        minλ = min(mod_λ)
        index_min = int(np.argmin(mod_λ))        
    except np.linalg.LinAlgError as Error:
        raise Error 
    return [minλ,v[index_min]]

#%%

print('Calcular los coeficientes de los campos')

lambdas = []
autovectores = []

kz_real_min_NM2 = []
kz_imag_min_NM2 = []
Elist_NM2 = []
for j in range(len(Elist_NM)):
    try:
        lambdda, autov = coef(j,R,ε2,ε3)
        lambdas.append(lambdda)
        autovectores.append(autov)
        Elist_NM2.append(Elist_NM[j])
        kz_real_min_NM2.append(kz_real_min_NM[j])
        kz_imag_min_NM2.append(kz_imag_min_NM[j])       
    except np.linalg.LinAlgError as Error:
        print(Error)   

plt.figure(figsize=tamfig)        
plt.title(title +', modo = %i, R = %i nm' %(modo,R*1e9),fontsize=tamtitle)
plt.plot(Elist_NM2,lambdas,'.',ms=10,alpha=0.7)
plt.ylabel('|λ|',fontsize=tamletra)
plt.xlabel('Energy [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.grid(1)

#%%

print('Se define el campo electrico E')

def Etot(indice,medio,ρ,φ,z,ε2,ε3):
    indice = int(indice)
    
    E = Elist_NM2[indice]
    kz_real = kz_real_min_NM2[indice]
    kz_imag = kz_imag_min_NM2[indice]
    
    kz = kz_real+1j*kz_imag     
    
    ω = (E/ħ)
    k0 = ω/c
    xz_real = kz.real/k0
    xz_imag = kz.imag/k0
    xz = xz_real+1j*xz_imag  

    def epsi(medio):
        if medio==2:
            return ε2
        if medio==3:
            return ε3
        
    def xt2(medio):
        inside = xz**2-μ0*epsi(medio)
        return inside
            
    def xt(medio):
        inside = xt2(medio)+0j
        xtmas=(inside)**(1/2)
        if xtmas.real>=0:
            xtfin=xtmas
        else:
            xtfin=-xtmas       
        return xtfin
    
    ρ_barra = ρ*k0
    
    def I(k):
        return special.iv(modo,xt(k)*ρ_barra)
    
    def derI(n):
        return special.ivp(modo,xt(n)*ρ_barra)
    
    def K(n):
        return special.kv(modo,xt(n)*ρ_barra)
    
    def derK(n):
        return special.kvp(modo,xt(n)*ρ_barra)
        
    def a(medio):
        if modo!=0:
            rta1 = xz*modo
            rta2 = ρ_barra*xt2(medio)
        else: 
            rta1 = 0
            rta2 = 1
        return rta1/rta2
    
    def b(medio):
        rta5 = xt(medio)
        return 1j*μ0/rta5
    
    def d(j):
        rta7 = modo*μ0
        rta9 = xt2(j)*ρ_barra
        return rta7/rta9
    
    def e(j):
        return 1j*xz/xt(j)
    
    coeff = np.array(autovectores[indice], dtype='complex')[0]

    coef_A2 = coeff[0]
    
    coef_C2 = coeff[1]
    
    coef_B3 = coeff[2]
    
    coef_D3 = coeff[3]
    
    if medio==3:
    
        Eρ = (-d(3)*coef_D3*K(3)+e(3)*coef_B3*derK(3))*np.e**(1j*modo*φ)*np.e**(-1j*kz*z)
        
        Eφ = (-a(3)*coef_B3*K(3)-b(3)*coef_D3*derK(3))*np.e**(1j*modo*φ)*np.e**(-1j*kz*z)
        
        Ez = coef_B3*K(3)*np.e**(1j*modo*φ)*np.e**(-1j*kz*z)
        
    elif medio==2:
        
        Eρ = (-d(2)*coef_C2*I(2)+e(2)*coef_A2*derI(2))*np.e**(1j*modo*φ)*np.e**(-1j*kz*z)
        
        Eφ = (-a(2)*coef_A2*I(2)-b(2)*coef_C2*derI(2))*np.e**(1j*modo*φ)*np.e**(-1j*kz*z)
        
        Ez = coef_A2*I(2)*np.e**(1j*modo*φ)*np.e**(-1j*kz*z)        
    
    campoE = (np.abs(Eρ))**2+(np.abs(Eφ))**2+(np.abs(Ez))**2
    
    return [campoE,np.abs(Ez)]

#%%
    
lambda_paper = 489 #nm
Energy_paper = h*c/(lambda_paper*10**-9)

print('Graficar el campo E en 3D para cada medio')

indicecerca=np.argmin(np.abs(np.array(Elist_NM2)-Energy_paper*np.ones(len(Elist_NM2)))) #el indice más cerca de Energy_paper
indice = int(indicecerca)
Energy = Elist_NM2[indice]
    
def E_2variable(x,y):
    z = 0     
    φ = np.arctan(y/x)
    ρ = (x**2+y**2)**(1/2)
    if np.abs(ρ)<=R:
        medio = 2
    else: 
        medio = 3
    return Etot(indice,medio,ρ,φ,z,ε2,ε3)[0]

n2 = 50
cota = 100*1e-9
x = np.linspace(-cota,cota,n2)
y = np.linspace(-cota,cota,n2)
X, Y = np.meshgrid(x, y)
f2 = np.vectorize(E_2variable)
Z = f2(X, Y)
Z = Z/np.max(Z)

limits = [min(x) , max(x), min(y) , max(y)]
limits = np.array(limits)*1e9
plt.figure(figsize=tamfig)
plt.xlabel('x [nm]',fontsize=tamletra)
plt.ylabel('y [nm]',fontsize=tamletra)
plt.title(title +', modo = %i, E = %.2f eV, R = %i nm' %(modo,Energy,R*1e9),fontsize=int(tamtitle*0.8)) 
im = plt.imshow(Z, extent = limits, cmap=plt.cm.jet, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.set_label('|E|$^2$',size=tamlegend)

#%%
    
print('Graficar el campo E vs x para cada medio')

n = len(Elist_NM2) - 1
list_indice = [0,int(n*1/4),int(n*1/2),int(n*3/4),n]

indicecerca=np.argmin(np.abs(np.array(Elist_NM2)-Energy_paper*np.ones(len(Elist_NM2)))) #el indice más cerca de Energy_paper
list_indice=[int(indicecerca)]

color = 0
plt.figure(figsize=tamfig)
for indice in list_indice: 
    color = color + 1
    def E_1variable(medio,ρ):
        φ,z = [0,0]
        return np.abs(Etot(indice,medio,ρ,φ,z,ε2,ε3)[0])

    N = int(1e3)
    lista_ρ = np.linspace(-100*1e-9,100*1e-9,N)
    campo_E = []
    for ρ in lista_ρ:
        if np.abs(ρ)<=R:
            campo_E.append(E_1variable(2,np.abs(ρ)))
        else:
            campo_E.append(E_1variable(3,np.abs(ρ)))
            
    campo_E = np.array(campo_E)/np.max(campo_E)        
    lista_ρ = np.array(lista_ρ)*1e9
    Energy = Elist_NM2[indice]
    plt.title(title +', modo = %i, R = %i nm' %(modo,R*1e9),fontsize=tamtitle)
    plt.plot(lista_ρ, campo_E,'.-',ms=10,alpha=0.7,label = 'Energia %.2f eV' %(Energy))
    plt.ylabel('Campo |E|$^2$',fontsize=tamletra)
    plt.xlabel('x[nm]',fontsize=tamletra)
    plt.legend(loc='upper left',markerscale=3,fontsize=tamlegend)
    plt.tick_params(labelsize = tamnum)
    plt.grid(1)
    plt.tight_layout(pad=2.5, w_pad=1, h_pad=1)

#%%
