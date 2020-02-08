
"""
Created on Fri Feb  7 13:28:30 2020

@author: leila
"""
from scipy import special #funciones de Bessel
import numpy as np
import os 
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import scipy.optimize as opt

#%%
lista_colores = ['coral','yellowgreen','midnightblue','purple','darkred','orange','aquamarine','hotpink','steelblue','darkgreen'] 

tamfig = (10,8)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

direc_Pffeifer = r'/home/leila/Desktop/Teo1_definitivo2/Fisica_Teorica_1/Pffeifer'

modo = 1
modo_real = modo
modo_virtual = modo
alfa = 1

ver_comportamiento = 0

#%%

print('Se define el determinante de 4x4 para los modos reales y para los modos virtuales')

def det_M_real_modes(Kz,Ω,modo,alfa):  #xz es en fracción de k0

    Kz_real = Kz[0]
    Kz_imag = Kz[1]   
    if Kz_imag==0:
        Kz=Kz_real
    else:
        Kz = Kz_real+1j*Kz_imag       

    V1_star = alfa*(-Ω**2+Kz**2+1)**(1/2)
    V2_star = alfa*(-Ω**2+Kz**2)**(1/2)

    X = V1_star*special.ivp(modo,V1_star)/special.iv(modo,V1_star)
    Y = V2_star*special.kvp(modo,V2_star)/special.kv(modo,V2_star)

    rta1 = (Ω**2-1)*(Ω**2-Kz**2)**2*X**2+Ω**2*(Ω**2-Kz**2-1)**2*Y**2
    rta2 = (2*Ω**2-1)*(Ω**2-Kz**2)*(Ω**2-Kz**2-1)*X*Y+modo**2*Kz**2

    return rta1-rta2

def det_M_virtual_modes(Kz,Ω,modo,alfa):    
    Kz_real = Kz[0]
    Kz_imag = Kz[1]  
    Kz = Kz_real+1j*Kz_imag  
    denominador = det_M_real_modes([Kz_real,0],Ω,modo,alfa)
    return np.abs((Ω**2-Kz**2)*(Ω**2-Kz**2-1)/denominador)

#%%

tipo_modo = 'real'

if ver_comportamiento==1:
    from scipy.signal import find_peaks
    
    print('Graficar la funcion det_M fijando im_kz = 0 y haciendo una curva por cada Ω')

    n = 5
    
    if tipo_modo == 'real':
        x = np.linspace(0.70,24,400)
        y = np.linspace(0.64,0.7,n)
    else:
        x = np.linspace(0.05,0.8,400)
        y = np.linspace(0.69,0.71,n)

    plt.figure(figsize=tamfig)
    detlist_tot = []
    color = 0
    for Ω in y:
        color = color + 1
        def det_M_1variable(Re_K):
            if tipo_modo == 'real':
                dett = det_M_real_modes([Re_K,0],Ω,modo,alfa)
                dett = np.abs(dett)
                dett = -np.log10(dett)
            elif tipo_modo == 'virtual':
                dett = det_M_virtual_modes([Re_K,0],Ω,modo,alfa)
            return dett
        
        detlist = []
        xlist = []
        for re_K in x:
            value = det_M_1variable(re_K)
            if np.isnan(value)==False: 
                detlist.append(value)
                xlist.append(re_K)

        detlist_tot.append(detlist)
        
        if tipo_modo == 'virtual':
            index_max = np.argmax(detlist)
            print(Ω, xlist[index_max])
        else:
            peaks, _ = find_peaks(detlist, threshold=0.01/max(detlist))
            print('omega', Ω)
            for j in range(len(peaks)):
                print('K', x[peaks[j]])
            print('')
        
        plt.title('Pffeifer. Modo = %i, Im($k_z$) = 0, alfa = %.1f' %(modo,alfa),fontsize=tamtitle)
        plt.plot(xlist,detlist,'.',color = lista_colores[color],ms=10,alpha=0.7,label='Ω = %.2f' %(Ω))
        for j in range(len(peaks)):
            plt.plot(x[peaks[j]],detlist[peaks[j]],'x',color = lista_colores[color],ms=15,alpha=0.7)
        plt.ylabel('-log$_{10}$(|det(L)|)',fontsize=tamletra)
        plt.xlabel('K',fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.legend(loc='best',markerscale=2,fontsize=tamlegend)  
        plt.grid(1)

#%%
    
print('Hallar la relacion de dispersion del modo REAL %i' %(modo_real))

Omega_real_obt1= []
det_min = []
Kz_obt1 = []

tol = 1e-9
N1 = int(1e3)
if modo_real ==1:
    if alfa==0.1:
        K_list1 = np.linspace(0.71,24,N1)
        x0 = 0.7
    elif alfa==1:
        K_list1 = np.linspace(0.53,18,N1)
        x0 = 0.51
elif modo_real ==5 and alfa==1:
    K_list1 = np.linspace(0.71,18,N1)
    x0 = 0.705   
elif alfa==10:
    if modo_real ==6:
        K_list1 = np.linspace(0.56,18,N1)
        x0 = 0.51   
    elif modo_real == 16:
        K_list1 = np.linspace(0.6755,18,N1)
        x0 = 0.67    

for Re_K in K_list1:
    def det_M_1variable(Ω):
        dett = det_M_real_modes([Re_K,0],Ω,modo_real,alfa)
        return np.abs(dett)
    minimos=opt.fsolve(det_M_1variable,x0,xtol=tol)    

    for j in range(len(minimos)):
        if minimos[j]<Re_K: #real modes
            Omega_real_obt1.append(minimos[j])
            Kz_obt1.append(Re_K)                 
            x0 = minimos[j]

#%%

print('Guardar la relacion de dispersion del modo REAL %i' %(modo_real))

os.chdir(direc_Pffeifer)
folder = 'alfa' + str(alfa)   
try:
    os.mkdir(folder)
except FileExistsError as error: 
    print(error)
os.chdir(direc_Pffeifer + '/' + folder)
np.savetxt('mod_real' + str(modo_real) + '_Omega.txt',[Omega_real_obt1],delimiter='\t')
np.savetxt('mod_real' + str(modo_real) + '_Kz.txt',[Kz_obt1],delimiter='\t')    

#%%

print('Hallar la relacion de dispersion del modo VIRTUAL %i' %(modo_virtual))
     
N2 = int(5*1e3)
#K_list2 = np.linspace(0.001,0.71,N2)

Omega_real_obt2 = []
Omega_imag_obt2 = []
Kz_obt2 = []

#if modo_virtual ==1:
#    if alfa==0.1: 
#        x0 = np.array([0.25,1e-15])
#    elif alfa==1:
#        K_list2 = np.linspace(0.001,0.49,N2)
#        x0 = np.array([0.6105,1e-1])
#    elif alfa == 10:
#        K_list2 = np.linspace(0.01,0.1,N2)
#        x0 = np.array([0.05,1e-1])            
# 
#elif modo_virtual in [2,3,4,5]:
#    if alfa==0.1: 
#        x0 = np.array([0.7075,1e-1])
#    elif alfa==1:
#        K_list2 = np.linspace(0.001,0.7,N2)
#        if modo_virtual ==2:
#            x0 = np.array([0.63,1e-1])
#        else: 
#            x0 = np.array([0.7,1e-1])
         
if alfa==10:  
    if modo_virtual == 2:
        K_list2 = np.linspace(0.001,0.32,N2)
        x0 = np.array([0.153,1e-1])
    elif modo_virtual == 3:
        K_list2 = np.linspace(0.001,0.42,N2)
        x0 = np.array([0.25,1e-1])
    elif modo_virtual == 4:
        K_list2 = np.linspace(0.001,0.48,N2)
        x0 = np.array([0.33,1e-1])
    elif modo_virtual == 5:
        K_list2 = np.linspace(0.001,0.53,N2)
        x0 = np.array([0.4,1e-1])
    elif modo_virtual ==6:
        K_list2 = np.linspace(0.001,0.56,N2)
        x0 = np.array([0.3,1e-1])   
    elif modo_virtual ==7:
        K_list2 = np.linspace(0.001,0.6,N2)
        x0 = np.array([0.51,1e-1])   
    elif modo_virtual ==8:
        K_list2 = np.linspace(0.001,0.61,N2)
        x0 = np.array([0.55,1e-1])   
    elif modo_virtual ==9:
        K_list2 = np.linspace(0.001,0.625,N2)
        x0 = np.array([0.59,1e-1])  
    elif modo_virtual ==10:
        K_list2 = np.linspace(0.001,0.64,N2)
        x0 = np.array([0.61,1e-1])
    elif modo_virtual ==11:
        K_list2 = np.linspace(0.001,0.65,N2)
        x0 = np.array([0.63,1e-1]) 
    elif modo_virtual ==12:
        K_list2 = np.linspace(0.001,0.655,N2)
        x0 = np.array([0.64,1e-1]) 
    elif modo_virtual ==13:
        K_list2 = np.linspace(0.001,0.665,N2)
        x0 = np.array([0.65,1e-1]) 
    elif modo_virtual ==14:
        K_list2 = np.linspace(0.001,0.668,N2)
        x0 = np.array([0.66,1e-1])
    elif modo_virtual ==15:
        K_list2 = np.linspace(0.001,0.673,N2)
        x0 = np.array([0.673,1e-1]) 
    elif modo_virtual == 16:
        K_list2 = np.linspace(0.001,0.6755,N2)
        x0 = np.array([0.67,1e-1])   

            
for Re_K in K_list2:

    def det_M4_2variable(x):
        Ω=x[0]+1j*x[1]
        dett = -det_M_virtual_modes([Re_K,0],Ω,modo_virtual,alfa)
        return (dett)
      
    res = minimize(det_M4_2variable, x0, method='Nelder-Mead')
    print(res.message)
    print(res.x)
    
    if res.x[0]>Re_K: #virtual modes
        Omega_real_obt2.append(res.x[0])
        Omega_imag_obt2.append(res.x[1])
        Kz_obt2.append(Re_K)   
        
        x0 = np.array([res.x[0],res.x[1]])
        
#%%

print('Guardar la relacion de dispersion del modo VIRTUAL %i' %(modo_virtual))

os.chdir(direc_Pffeifer)
folder = 'alfa' + str(alfa)   
try:
    os.mkdir(folder)
except FileExistsError as error: 
    print(error)
os.chdir(direc_Pffeifer + '/' + folder)
np.savetxt('mod_virtual' + str(modo_virtual) + '_Re_Omega.txt',[Omega_real_obt2],delimiter='\t')
np.savetxt('mod_virtual' + str(modo_virtual) + '_Im_Omega.txt',[Omega_imag_obt2],delimiter='\t')
np.savetxt('mod_virtual' + str(modo_virtual) + '_Kz.txt',[Kz_obt2],delimiter='\t')    

#%%
     
print('Graficar la relacion de dispersion obtenida')

x = np.linspace(0.1,0.71,300)
    
plt.figure(figsize=tamfig)
#plt.plot(Kz_obt1,Omega_real_obt1,'.',ms=8,color = lista_colores[0],alpha = 0.8,label='modo real')
plt.plot(Kz_obt2,Omega_real_obt2,'.',ms=8,color = lista_colores[1],alpha = 0.8,label='modo virtual')
plt.plot(x,x,'.-',ms=8,color = lista_colores[2],alpha = 0.8,label='Ω = K')
plt.title('Relación de dispersión. Modo = %i, alfa = %.1f' %(modo,alfa),fontsize=tamtitle)
plt.ylabel('Ω',fontsize=tamletra)
plt.xlabel('K',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
plt.grid(1)
#plt.ylim([0.69,0.71])
#plt.tight_layout(1)

#%%
