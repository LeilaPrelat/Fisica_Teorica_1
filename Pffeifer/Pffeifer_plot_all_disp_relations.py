
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

modo = 5
modo_real = modo
modo_virtual = modo
alfa = 10

#%%
print('Importar las relaciones de dispersion y graficar modo real y virtual')

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

if alfa ==0.1:
    list_modos = [1,2,3,5]
    x = np.linspace(0.69,0.71,300)
elif alfa ==1:
    list_modos = [1,2,3,4,5,6]
    x = np.linspace(0.45,0.75,300)
elif alfa ==10:
    list_modos = np.linspace(1,15,15)
    x = np.linspace(0.001,0.75,300)

color = 0
plt.figure(figsize=tamfig)
for modo_virtual in list_modos:
    modo_virtual = int(modo_virtual)
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
    
    if alfa!=10:
        color_graph = lista_colores[color]
        plt.plot(Kz_obt2,Omega_real_obt2,'.-',color = color_graph, ms=8,alpha=0.7,label='modo = %i' %(modo_virtual)) 
    else: 
        plt.plot(Kz_obt2,Omega_real_obt2,'.-', ms=8,alpha=0.7,label='modo = %i' %(modo_virtual)) 
        
plt.title('Relación de dispersión. Modos virtuales. alfa = %.1f' %(alfa),fontsize=tamtitle)
if alfa==10:
    tamlegend = 18
    plt.plot(x,x,'.-',ms=8,color=lista_colores[0],alpha = 0.7,label='Ω = K')
    tamlegend = int(tamlegend*0.8)
else: 
    plt.plot(x,x,'.-',ms=8,color=lista_colores[color+1],alpha = 0.7,label='Ω = K')
plt.ylabel('Ω',fontsize=tamletra)
plt.xlabel('K',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
plt.grid(1)
plt.ylim([x[0],x[-1]])

#%%