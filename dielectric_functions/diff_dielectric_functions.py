#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 16:07:56 2020

@author: leila
"""

from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import os

#%%

lista_colores = ['coral','yellowgreen','midnightblue','purple','darkred','orange','aquamarine','hotpink','steelblue','darkgreen'] 

tamfig = (12,10)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

graficar = 0 #si se esta importando la funcion ε_functions, este valor tiene que ser 0 (graficar==0)
        
direc = r'/home/leila/Desktop/Teo1_definitivo2/dielectric_functions'

#%%

print('Se importan las funciones dielectricas de los 3 medios')

class diff_ε:    
    
    def load_data2(direc):
        os.chdir(direc)
        try:
            energy_SiO2 = np.loadtxt('energy_Palik.txt',delimiter='\t')
            epsi_real_SiO2 = np.loadtxt('epsi_real_Palik.txt',delimiter='\t')
            epsi_imag_SiO2 = np.loadtxt('epsi_imag_Palik.txt',delimiter='\t')
            epsi_imag_SiO2 = -np.array(epsi_imag_SiO2) #convencion adoptada
            
            energy_Ag = np.loadtxt('energy_JC.txt',delimiter='\t')   
            epsi_real_Ag = np.loadtxt('epsi_real_JC.txt',delimiter='\t')
            epsi_imag_Ag = np.loadtxt('epsi_imag_JC.txt',delimiter='\t')
                 
            return [energy_SiO2,epsi_real_SiO2,epsi_imag_SiO2,energy_Ag,epsi_real_Ag,epsi_imag_Ag]
        except OSError or IOError:
            raise OSError('Los archivos .txt NO estan en la direccion brindada como variable')  
    
    [energy_SiO2,epsi_real_SiO2,epsi_imag_SiO2,energy_Ag,epsi_real_Ag,epsi_imag_Ag] = load_data2(direc)
    min_Elist = np.max([np.min(energy_Ag),np.min(energy_SiO2)])
    max_Elist = np.min([np.max(energy_Ag),np.max(energy_SiO2)])
    
    def ε_SiO2(E):
        delta = 1e-8    
        try:
            f_real2 = interp1d(energy_SiO2, epsi_real_SiO2)
            f_imag2 = interp1d(energy_SiO2, epsi_imag_SiO2)
            rta1 = f_real2(E+delta)+1j*f_imag2(E+delta)
            rta2 = f_real2(E-delta)+1j*f_imag2(E-delta)
            return (rta1-rta2)/(2*delta)
        except ValueError:       
            raise ValueError('El rango valido de energia es:', (min_Elist,max_Elist),'eV')  

    def ε_Ag(E):
        delta = 1e-8
        try:
            f_real1 = interp1d(energy_Ag, epsi_real_Ag)
            f_imag1 = interp1d(energy_Ag, epsi_imag_Ag)
            rta3 = f_real1(E+delta)+1j*f_imag1(E+delta)
            rta4 = f_real1(E-delta)+1j*f_imag1(E-delta)
            return (rta3-rta4)/(2*delta)
        except ValueError:       
            raise ValueError('El rango valido de energia es:', (min_Elist,max_Elist),'eV')  

    def ε_CdS(E):
        delta = 1e-8
        def ε_CdS2(ħω):
    #        ħω=-ħω #convencion nuestra
            εinf_2 = 1.58  
            def Epsilon_0(ħω):
                E0A = 2.482  #eV
                E0B = 2.496  #eV
                E0C = 2.555  #eV
                A0A = 13.0   #eV**1.5
                A0B = 6.6    #eV**1.5
                A0C = 6.4    #eV**1.5
                Γ0A = 0.055  #eV
                Γ0B = 0.055  #eV
                Γ0C = 0.055  #eV
                χ0 = (ħω+1j*Γ0A) / E0A
                fχ0 = χ0**-2 * ( 2-(1+χ0)**0.5-(1-χ0)**0.5 )
                ε = A0A*E0A**-1.5 * fχ0
                χ0 = (ħω+1j*Γ0B) / E0B
                fχ0 = χ0**-2 * ( 2-(1+χ0)**0.5-(1-χ0)**0.5 )
                ε += A0B*E0B**-1.5 * fχ0
                χ0 = (ħω+1j*Γ0C) / E0C
                fχ0 = χ0**-2 * ( 2-(1+χ0)**0.5-(1-χ0)**0.5 )
                ε += A0C*E0C**-1.5 * fχ0
                return ε
            
            def Epsilon_0x(ħω):
                G0  = 0.03   #eV
                A0xA= 0.045  #eV
                A0xB= 0.023  #eV
                A0xC= 0.022  #eV
                Γ0A = 0.055  #eV
                Γ0B = 0.055  #eV
                Γ0C = 0.055  #eV
                E0A = 2.482  #eV
                E0B = 2.496  #eV
                E0C = 2.555  #eV
                ε  = A0xA*(E0A-G0-ħω-1j*Γ0A)**-1
                ε += A0xB*(E0B-G0-ħω-1j*Γ0B)**-1
                ε += A0xC*(E0C-G0-ħω-1j*Γ0C)**-1
                return ε
            
            def Epsilon_1(ħω):
                E1xA= 4.84   #eV
                E1xB= 5.45   #eV
                B1xA= 1.60   #eV
                B1xB= 1.25   #eV
                Γ1A = 0.42   #eV
                Γ1B = 0.38   #eV
                ε  = B1xA*(E1xA-ħω-1j*Γ1A)**-1
                ε += B1xB*(E1xB-ħω-1j*Γ1B)**-1
                return ε
            
            def Epsilon_0pr(ħω):
                E0pr= 6.2    #eV
                C_CdS   = 0.30
                γ   = 0.12
                χ = ħω/E0pr
                return C_CdS / (1 - χ**2 - 1j*χ*γ)
            
            ε0   = Epsilon_0(ħω)
            ε0x  = Epsilon_0x(ħω)
            ε1   = Epsilon_1(ħω)
            ε0pr = Epsilon_0pr(ħω)
            ε = ε0 + ε0x + ε1 + ε0pr + εinf_2
            return ε   
        
        return (ε_CdS2(E+delta)-ε_CdS2(E-delta))/(2*delta)
    
#%%
        
[energy_SiO2,epsi_real_SiO2,epsi_imag_SiO2,energy_Ag,epsi_real_Ag,epsi_imag_Ag] = diff_ε.load_data2(direc)
min_Elist = np.max([np.min(energy_Ag),np.min(energy_SiO2)])
max_Elist = np.min([np.max(energy_Ag),np.max(energy_SiO2)])

n = 100
E = np.linspace(0.8,3.8,n)

ε_SiO2 = diff_ε.ε_SiO2
ε_Ag = diff_ε.ε_Ag
ε_CdS = diff_ε.ε_CdS

if graficar ==1:
    
    print('Palik: Graficar ε de SiO2 vs E')

    plt.figure(figsize=tamfig)
    plt.plot(energy_SiO2,epsi_real_SiO2,'o',ms=10,color=lista_colores[0],label='Re(ε)')
    plt.plot(E,ε_SiO2(E).real,'-',ms=10,color=lista_colores[2],label='interpolacion Re(ε)')
    plt.title('ε de SiO2 de Palik',fontsize=tamtitle)
    plt.xlabel('Energy [eV]',fontsize=tamletra)
    plt.ylabel('dε/dE',fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)
    
    plt.figure(figsize=tamfig)
    plt.plot(energy_SiO2,epsi_imag_SiO2,'o',ms=10,color=lista_colores[1],label='Im(ε)')
    plt.plot(E,ε_SiO2(E).imag,'-',ms=10,color=lista_colores[3],label=' interpolacion Im(ε)')
    plt.title('ε de SiO2 de Palik',fontsize=tamtitle)
    plt.xlabel('Energy [eV]',fontsize=tamletra)
    plt.ylabel('dε/dE',fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)

    print('JC: Graficar ε de Ag vs E')

    plt.figure(figsize=tamfig)
    plt.plot(energy_Ag,epsi_real_Ag,'o',ms=10,color=lista_colores[0],label='Re(ε)')
    plt.plot(E,ε_Ag(E).real,'-',ms=10,color=lista_colores[2],label='interpolacion Re(ε)')
    plt.title('ε de Ag de JC',fontsize=tamtitle)
    plt.xlabel('Energy [eV]',fontsize=tamletra)
    plt.ylabel('dε/dE',fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)
    
    plt.figure(figsize=tamfig)
    plt.plot(energy_Ag,epsi_imag_Ag,'o',ms=10,color=lista_colores[1],label='Im(ε)')
    plt.plot(E,ε_Ag(E).imag,'-',ms=10,color=lista_colores[3],label=' interpolacion Im(ε)')
    plt.title('ε de Ag de JC',fontsize=tamtitle)
    plt.xlabel('Energy [eV]',fontsize=tamletra)
    plt.ylabel('dε/dE',fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)

    print('Graficar ε de CdS vs E')
    
    plt.figure(figsize=tamfig)
    plt.plot(E, ε_CdS(E).real,'.',ms=10,color=lista_colores[0],alpha=0.7,label='Re(ε)')
    plt.plot(E, ε_CdS(E).imag,'.',ms=10,color=lista_colores[1],alpha=0.7,label='Im(ε)')
    plt.title('ε de CdS, modelo DL',fontsize=tamtitle)
    plt.xlabel('Energy (eV)',fontsize=tamletra)
    plt.ylabel('ε',fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
    plt.grid(1)
    
#%%