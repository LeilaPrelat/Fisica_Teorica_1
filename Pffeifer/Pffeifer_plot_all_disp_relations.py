
"""
Created on Fri Feb  7 13:28:30 2020

@author: leila
"""
import numpy as np
import os 
import matplotlib.pyplot as plt

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

#%%
print('Importar las relaciones de dispersion y graficar')

folder = 'alfa' + str(alfa)  
os.chdir(direc_Pffeifer + '/' + folder)

try:   
    Kz_obt1 = np.loadtxt('mod_real' + str(modo_real) + '_Kz.txt',delimiter='\t')
    Omega_real_obt1 = np.loadtxt('mod_real' + str(modo_real) + '_Omega.txt',delimiter='\t')      
    Kz_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Kz.txt',delimiter='\t')
    Omega_real_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Re_Omega.txt',delimiter='\t')   
    Omega_imag_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Im_Omega.txt',delimiter='\t')         
except OSError or IOError as error:
    print(error)
    direc = input('path de los .txt de las relaciones de dispersion')
    os.chdir(direc)
    Kz_obt1 = np.loadtxt('mod_real' + str(modo_real) + '_Kz.txt',delimiter='\t')
    Omega_real_obt1 = np.loadtxt('mod_real' + str(modo_real) + '_Omega.txt',delimiter='\t')      
    Kz_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Kz.txt',delimiter='\t')
    Omega_real_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Re_Omega.txt',delimiter='\t')   
    Omega_imag_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Im_Omega.txt',delimiter='\t')      

x = np.linspace(0.1,0.71,300)
    
plt.figure(figsize=tamfig)
plt.plot(Kz_obt1,Omega_real_obt1,'.',ms=8,color = lista_colores[0],alpha = 0.8,label='modo real')
plt.plot(Kz_obt2,Omega_real_obt2,'.',ms=8,color = lista_colores[1],alpha = 0.8,label='modo virtual')
plt.plot(x,x,'.-',ms=8,color = lista_colores[2],alpha = 0.8,label='Ω = K')
plt.title('Relación de dispersión. Modo = %i, alfa = %.1f' %(modo,alfa),fontsize=tamtitle)
plt.ylabel('Ω',fontsize=tamletra)
plt.xlabel('K',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
plt.grid(1)

#%%

print('Barrido de los modos virtuales')
print('Importar las relaciones de dispersion y graficar')

list_modos = [1,2,3,5]
x = np.linspace(0.69,0.71,300)

color = 0
plt.figure(figsize=tamfig)
for modo_virtual in list_modos:
    color = color +1

    folder = 'alfa' + str(alfa)  
    os.chdir(direc_Pffeifer + '/' + folder)
    
    try:   
        Kz_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Kz.txt',delimiter='\t')
        Omega_real_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Re_Omega.txt',delimiter='\t')   
        Omega_imag_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Im_Omega.txt',delimiter='\t')     
    except OSError or IOError as error:
        print(error)
        direc = input('path de los .txt de las relaciones de dispersion')
        os.chdir(direc)
        Kz_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Kz.txt',delimiter='\t')
        Omega_real_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Re_Omega.txt',delimiter='\t')   
        Omega_imag_obt2 = np.loadtxt('mod_virtual' + str(modo_virtual) + '_Im_Omega.txt',delimiter='\t')  
    
    
    plt.plot(Kz_obt2,Omega_real_obt2,'.',ms=8,color=lista_colores[color],alpha=0.7,label='modo = %i' %(modo_virtual)) 

plt.title('Relación de dispersión. Modos virtuales. alfa = %.1f' %(alfa),fontsize=tamtitle)
plt.plot(x,x,'.-',ms=8,color=lista_colores[color+1],alpha = 0.7,label='Ω = K')
plt.ylabel('Ω',fontsize=tamletra)
plt.xlabel('K',fontsize=tamletra)
plt.ylim([0.69,0.71])
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
plt.grid(1)

#%%