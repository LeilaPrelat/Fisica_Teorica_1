#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 23:22:25 2020

@author: leila
"""
from scipy import special #funciones de Bessel
from numpy.linalg import det
import numpy as np

#%% 

print('3 medios: Se define el determinante de 8x8')

def det_M3medios(kz,E,modo,r1,r2,ε1,ε2,ε3):
    π = np.pi
    c = 2.99792458*1e8 #m/s (SI units)
    h = 4.135667731*1e-15 #eV/s
    ħ = h/(2*π)
    μ0 = 1
    
    kz_real = kz[0]
    kz_imag = kz[1]   
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
#        return inside**(1/2)
                
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
#        A = np.linalg.svd(M,full_matrices=True,compute_uv=False)
#        rta = np.min(A)
        rta = np.linalg.det(M)
    except np.linalg.LinAlgError as Error:
        raise Error 
    return (np.abs(rta))

#%%

print('2 medios: Se define el determinante de 4x4')


def det_M2medios(kz,E,modo,R,ε2,ε3):
    π = np.pi
    c = 2.99792458*1e8 #m/s (SI units)
    h = 4.135667731*1e-15 #eV/s
    ħ = h/(2*π)
    μ0 = 1
    
    kz_real = kz[0]
    kz_imag = kz[1]   
    kz = kz_real+1j*kz_imag       

    ω = (E/ħ)
    k0 = ω/c
    xz_real = kz.real/k0
    xz_imag = kz.imag/k0
    xz = xz_real+1j*xz_imag    

    def epsi(medio):
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
#        A = np.linalg.svd(M,full_matrices=True,compute_uv=False)
#        rta = np.min(A) #minimo de los valores singulares 
        rta = np.abs(np.linalg.det(M)) #modulo del determinante 
    except np.linalg.LinAlgError as Error:
        raise Error 
    return rta

#%%