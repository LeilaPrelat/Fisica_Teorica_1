#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 11:00:11 2020

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

modo = 0
R1 =70*1e-9
R2 = 85*1e-9

ε_path = r'/home/leila/Desktop/Teo1_definitivo2/dielectric_functions'
os.chdir(ε_path)

det_path = r'/home/leila/Desktop/Teo1_definitivo2'

#%%

print('Importar las funciones dielectricas de los medios')

try:
    sys.path.insert(1, ε_path)
    from dielectric_functions import ε_functions
    from diff_dielectric_functions import diff_ε
except ModuleNotFoundError as error:
    print(error)
    ε_path = input('path de la carpeta dielectric_functions')
    sys.path.insert(1, ε_path)

ε_SiO2 = ε_functions.ε_SiO2
ε_Ag = ε_functions.ε_Ag
ε_CdS = ε_functions.ε_CdS

diff_ε_SiO2 = diff_ε.ε_SiO2
diff_ε_Ag = diff_ε.ε_Ag
diff_ε_CdS = diff_ε.ε_CdS

try: 
    ε_functions.load_data(ε_path)
    plt.close('all')

except OSError or IOError as error:
    print(error)
    ε_path2 = input('path de los archivos .txt de la carpeta dielectric_functions')
    ε_functions.load_data(ε_path2)
    plt.close('all')

ε1,ε2,ε3 = ε_Ag,ε_SiO2,ε_Ag
der_ε1,der_ε2,der_ε3 = diff_ε_Ag,diff_ε_SiO2,diff_ε_Ag

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

print('Importar la relacion de dispersion')

folder_mediums = title1 + '_' + title2 + '_' + title3
folder_R1 = 'R1' + '_' + str(int(R1*1e9)) 
os.chdir(det_path + '/' + folder_mediums + '/' + folder_R1)

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
print('Se definen los coeficientes del determinante de 8x8')

def coef(indice,r1,r2,ε1,ε2,ε3):

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
        if medio==1:            
            return ε1(E)
        if medio==2:
            return ε2(E)
        if medio==3:
            return ε3(E)
                
    def xt2(medio):
        inside = xz**2-μ0*epsi(medio)
        return inside
            
    def xt(medio):
        inside = xt2(medio)+0j
        xtmas=(inside)**(1/2)
        if (xtmas).real>=0:
            xtfin=xtmas
        else:
            xtfin=-xtmas       
        return xtfin
                
    def rbarra(medio):
        if medio==1:
            return r1*k0
        if medio==2:
            return r2*k0
    
    def I(k,m):
        return special.iv(modo,xt(k)*rbarra(m))
    
    def derI(n,m):
        return special.ivp(modo,xt(n)*rbarra(m))
    
    def K(n,m):
        return special.kv(modo,xt(n)*rbarra(m))
    
    def derK(n,m):
        return special.kvp(modo,xt(n)*rbarra(m))
    
    def a(j,l):
        if modo!=0:
            rta1 = xz*modo*xt2(j)
            rta2 = rbarra(l)
        else: 
            rta1 = 0
            rta2 = 1
            
        return rta1/rta2

    def b(j,k):
        return 1j*μ0*xt(j)*xt2(k)

    def d(j,k):
        return 1j*epsi(j)*xt(j)*xt2(k)
    
    A = [I(1,1),0,-I(2,1),-K(2,1),0,0,0,0]
    B = [0,I(1,1),0,0,-I(2,1),-K(2,1),0,0]
    C = [0,0,I(2,2),K(2,2),0,0,-K(3,2),0]
    D = [0,0,0,0,I(2,2),K(2,2),0,-K(3,2)]      
       
    F = [a(2,1)*I(1,1),b(1,2)*derI(1,1),-a(1,1)*I(2,1),-a(1,1)*K(2,1),-b(2,1)*derI(2,1),-b(2,1)*derK(2,1),0,0]
    G = [d(1,2)*derI(1,1),-a(2,1)*I(1,1),-d(2,1)*derI(2,1),-d(2,1)*derK(2,1),a(1,1)*I(2,1),a(1,1)*K(2,1),0,0]      
    H = [0,0,a(3,2)*I(2,2),a(3,2)*K(2,2),b(2,3)*derI(2,2),b(2,3)*derK(2,2),-a(2,2)*K(3,2),-b(3,2)*derK(3,2)]
    I = [0,0,d(2,3)*derI(2,2),d(2,3)*derK(2,2),-a(3,2)*I(2,2),-a(3,2)*K(2,2),-d(3,2)*derK(3,2),a(2,2)*K(3,2)]

    M = np.matrix([A,B,C,D,F,G,H,I], dtype='complex')
    try:
        λ, v = LA.eig(M)
        mod_λ = []
        for autovalor in λ:
            mod_λ.append(np.abs(autovalor))
        minλ = min(mod_λ)
        index_min = int(np.argmin(mod_λ))        
    except np.linalg.LinAlgError as Error:
        raise Error 
    return [minλ,v[:,index_min]]

#%%

print('Calcular los coeficientes de los campos')

lambdas = []
autovectores = []

kz_real_min_NM2 = []
kz_imag_min_NM2 = []
Elist_NM2 = []
for j in range(len(Elist_NM)):
    try:
        lambdda, autov = coef(j,R1,R2,ε1,ε2,ε3)
        lambdas.append(lambdda)
        autovectores.append(autov)
        Elist_NM2.append(Elist_NM[j])
        kz_real_min_NM2.append(kz_real_min_NM[j])
        kz_imag_min_NM2.append(kz_imag_min_NM[j])       
    except np.linalg.LinAlgError as Error:
        print(Error)   

plt.figure(figsize=tamfig)        
plt.title(title +', modo = %i, R1 = %i nm, R2 = %i nm' %(modo,R1*1e9,R2*1e9),fontsize=tamtitle)
plt.plot(Elist_NM2,lambdas,'.',ms=10,alpha=0.7)
plt.ylabel('|λ|',fontsize=tamletra)
plt.xlabel('Energy [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.grid(1)

#%%

print('Se define el campo electrico E')

def Etot(indice,medio,ρ,φ,z,ε1,ε2,ε3):
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
       
    def coef2(medio):
        coeff = np.array(autovectores[indice], dtype='complex')
        if medio==1:
            return [coeff[0][0],0,coeff[1][0],0]
        if medio==2:
            return [coeff[2][0],coeff[3][0],coeff[4][0],coeff[5][0]]
        if medio==3:
            return [0,coeff[6][0],0,coeff[7][0]]
    
    def epsi(medio):
        if medio==1:            
            return ε1(E)
        if medio==2:
            return ε2(E)
        if medio==3:
            return ε3(E)
                
    def xt2(medio):
        inside = xz**2-μ0*epsi(medio)
        return inside
            
    def xt(medio):
        inside = xt2(medio)+0j
        xtmas=(inside)**(1/2)
        if (xtmas).real>=0:
            xtfin=xtmas
        else:
            xtfin=-xtmas       
        return xtfin
    
    ρ_barra = ρ*k0
    
    def a(j):
        if modo!=0:
            rta7 = modo*μ0
            rta8 = xt2(j)*ρ_barra
        else: 
            rta7 = 0
            rta8 = 1
        return rta7/rta8
    
    def b(j):
        return 1j*xz/xt(j)
    
    def e(j):
        if modo!=0:
            rta7 = modo*xz
            rta8 = xt2(j)*ρ_barra
        else: 
            rta7 = 0
            rta8 = 1
        return rta7/rta8      
    
    def f(j):
        return 1j*μ0/xt(j)    
      
    def I(n):
        return special.iv(modo,xt(n)*ρ_barra)
    
    def derI(n):
        return special.ivp(modo,xt(n)*ρ_barra)
    
    def K(n):
        return special.kv(modo,xt(n)*ρ_barra)
    
    def derK(n):
        return special.kvp(modo,xt(n)*ρ_barra)
        
    campo1 = coef2(medio)
    campo2 = [b(medio)*derI(medio),b(medio)*derK(medio),-a(medio)*I(medio),-a(medio)*K(medio)]

    campo3 = np.dot(campo1,campo2)
    campo4 = np.e**(1j*modo*φ)
    campo5 = np.e**(-1j*kz*z)
    
    Eρ = campo3*campo4*campo5
        
    campo6 = coef2(medio)
    campo7 = [-e(medio)*I(medio),-e(medio)*K(medio),-f(medio)*derI(medio),-f(medio)*derK(medio)]

    campo8 = np.dot(campo6,campo7)
    campo9 = np.e**(1j*modo*φ)
    campo10 = np.e**(-1j*kz*z)
    
    Eφ = campo8*campo9*campo10
    
    campo11 = coef2(medio)
    campo12 = [I(medio),K(medio),0,0]

    campo13 = np.dot(campo11,campo12)
    campo14 = np.e**(1j*modo*φ)
    campo15 = np.e**(-1j*kz*z)
    
    Ez = campo13*campo14*campo15

    campoE = (np.abs(Eρ))**2+(np.abs(Eφ))**2+(np.abs(Ez))**2
    
    return [campoE,np.abs(Ez)]


#%%

print('Se define el campo magnetico H')

def Htot(indice,medio,ρ,φ,z,ε1,ε2,ε3):
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
       
    def coef3(medio):
        coeff = np.array(autovectores[indice], dtype='complex')
        if medio==1:
            return [coeff[0][0],0,coeff[1][0],0]
        if medio==2:
            return [coeff[2][0],coeff[3][0],coeff[4][0],coeff[5][0]]
        if medio==3:
            return [0,coeff[6][0],0,coeff[7][0]]
    
    def epsi(medio):
        if medio==1:            
            return ε1(E)
        if medio==2:
            return ε2(E)
        if medio==3:
            return ε3(E)
    
    def xt2(medio):
        inside = xz**2-μ0*epsi(medio)
        return inside
            
    def xt(medio):
        inside = xt2(medio)+0j
        xtmas=(inside)**(1/2)
        if (xtmas).real>=0:
            xtfin=xtmas
        else:
            xtfin=-xtmas       
        return xtfin
    
    ρ_barra = ρ*k0
    
    def a(j):
        if modo!=0:
            rta7 = modo*epsi(j)
            rta8 = xt2(j)*ρ_barra
        else: 
            rta7 = 0
            rta8 = 1
        return rta7/rta8
    
    def b(j):
        return 1j*xz/xt(j)

    def e(j):
        if modo!=0:
            rta7 = modo*xz
            rta8 = xt2(j)*ρ_barra
        else: 
            rta7 = 0
            rta8 = 1
        return rta7/rta8      
    
    def f(j):
        return 1j*epsi(j)/xt(j)    
      
    def I(n):
        return special.iv(modo,xt(n)*ρ_barra)
    
    def derI(n):
        return special.ivp(modo,xt(n)*ρ_barra)
    
    def K(n):
        return special.kv(modo,xt(n)*ρ_barra)
    
    def derK(n):
        return special.kvp(modo,xt(n)*ρ_barra)
        
    campo1 = coef3(medio)
    campo2 = [a(medio)*I(medio),a(medio)*K(medio),b(medio)*derI(medio),b(medio)*derK(medio)]

    campo3 = np.dot(campo1,campo2)
    campo4 = np.e**(1j*modo*φ)
    campo5 = np.e**(-1j*kz*z)
    
    Hρ = campo3*campo4*campo5

    campo6 = coef3(medio)
    campo7 = [f(medio)*derI(medio),f(medio)*derK(medio),-e(medio)*I(medio),-e(medio)*K(medio)]

    campo8 = np.dot(campo6,campo7)
    campo9 = np.e**(1j*modo*φ)
    campo10 = np.e**(-1j*kz*z)
    
    Hφ = campo8*campo9*campo10

    campo11 = coef3(medio)
    campo12 = [0,0,I(medio),K(medio)]

    campo13 = np.dot(campo11,campo12)
    campo14 = np.e**(1j*modo*φ)
    campo15 = np.e**(-1j*kz*z)
    
    Hz = campo13*campo14*campo15

    campoH = (np.abs(Hρ))**2+(np.abs(Hφ))**2+(np.abs(Hz))**2
    
    return [campoH,np.abs(Hz)]


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
    if np.abs(ρ)<=R1:
        medio = 1
    elif np.abs(ρ)<=R2:
        medio = 2
    else: 
        medio = 3
    return Etot(indice,medio,ρ,φ,z,ε1,ε2,ε3)[0]

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
plt.title(title +', modo = %i, E = %.2f eV, R1 = %i nm, R2 = %i nm' %(modo,Energy,R1*1e9,R2*1e9),fontsize=int(tamtitle*0.8)) 
im = plt.imshow(Z, extent = limits, cmap=plt.cm.jet, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.set_label('|E|$^2$',size=tamlegend)

#%%
    
lambda_paper = 489 #nm
Energy_paper = h*c/(lambda_paper*10**-9)

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
        return np.abs(Etot(indice,medio,ρ,φ,z,ε1,ε2,ε3)[0])

    N = int(1e3)
    lista_ρ = np.linspace(-100*1e-9,100*1e-9,N)
    campo_E = []
    for ρ in lista_ρ:
        if np.abs(ρ)<=R1:
            campo_E.append(E_1variable(1,np.abs(ρ)))
        elif np.abs(ρ)<=R2:
            campo_E.append(E_1variable(2,np.abs(ρ)))
        else:
            campo_E.append(E_1variable(3,np.abs(ρ)))
            
    campo_E = np.array(campo_E)/np.max(campo_E)        
    lista_ρ = np.array(lista_ρ)*1e9
    Energy = Elist_NM2[indice]
    plt.title(title +', modo = %i, R1 = %i nm, R2 = %i nm' %(modo,R1*1e9,R2*1e9),fontsize=tamtitle)
    plt.plot(lista_ρ, campo_E,'.-',color = lista_colores[color],ms=10,alpha=0.7,label = 'Energia %.2f eV' %(Energy))
    plt.ylabel('Campo |E|$^2$',fontsize=tamletra)
    plt.xlabel('x[nm]',fontsize=tamletra)
    plt.legend(loc='upper left',markerscale=3,fontsize=tamlegend)
    plt.tick_params(labelsize = tamnum)
    plt.grid(1)
    plt.tight_layout(pad=2.5, w_pad=1, h_pad=1)

#%%

print('Definir la funcion densidad de energia')

def energy_density(indice,medio,ρ,φ,z,ε1,ε2,ε3):
    indice = int(indice)
    campo_H = np.abs(Htot(indice,medio,ρ,φ,z,ε1,ε2,ε3)[0])
    campo_E = np.abs(Etot(indice,medio,ρ,φ,z,ε1,ε2,ε3)[0])
    
    E = Elist_NM2[indice]
        
    def epsi(medio):
        if medio==1:            
            return ε1(E)
        if medio==2:
            return ε2(E)
        if medio==3:
            return ε3(E)
    
    def der(E):
        if medio==1:            
            return der_ε1(E)
        if medio==2:
            return der_ε2(E)
        if medio==3:
            return der_ε3(E)  
    
    rta = der(E)*campo_E+μ0*campo_H
    return 0.5*(rta.real)
    
#%%
    
print('Graficar el campo energy density vs x para cada medio')

list_indice=[int(indicecerca)]

color = 0
plt.figure(figsize=tamfig)
for indice in list_indice:
    color = color + 1    
    def energy_density1variable(medio,ρ):
        φ,z = [0,0]
        return energy_density(indice,medio,ρ,φ,z,ε1,ε2,ε3)
    
    N = int(1e3)
    lista_ρ = np.linspace(-100*1e-9,100*1e-9,N)
    e_density = []
    for ρ in lista_ρ:
        if np.abs(ρ)<=R1:
            e_density.append(energy_density1variable(1,np.abs(ρ)))
        elif np.abs(ρ)<=R2:
            e_density.append(energy_density1variable(2,np.abs(ρ)))
        else:
            e_density.append(energy_density1variable(3,np.abs(ρ)))
            
    e_density = np.array(e_density)/np.max(e_density)        
    lista_ρ = np.array(lista_ρ)*1e9
    Energy = Elist_NM2[indice]
    plt.title(title +', modo = %i, R1 = %i nm, R2 = %i nm' %(modo,R1*1e9,R2*1e9),fontsize=tamtitle)
    plt.plot(lista_ρ, e_density,'.-',color = lista_colores[color],ms=10,alpha=0.7,label = 'Energia %.2f eV' %(Energy))
    plt.ylabel('$\omega_e$',fontsize=tamletra)
    plt.xlabel('x[nm]',fontsize=tamletra)
    plt.legend(loc='upper left',markerscale=3,fontsize=tamlegend)
    plt.tick_params(labelsize = tamnum)
    plt.grid(1)
    plt.tight_layout(pad=2.5, w_pad=1, h_pad=1)

#%%