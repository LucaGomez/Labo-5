# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 16:49:47 2022

@author: luo
"""

import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython
from scipy.io import loadmat
from scipy.optimize import curve_fit
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from scipy import ndimage

from IPython import get_ipython

get_ipython().run_line_magic('matplotlib', 'qt5')
plt.rcParams.update({'font.size': 12})

#%%

#Varias agitaciones

r0,vt0,er0,evt0=np.loadtxt(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\promedios 1.txt')
r1,vt1,er1,evt1=np.loadtxt(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\promedios 2.txt')
r2,vt2,er2,evt2=np.loadtxt(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\promedios 3.txt')
r3,vt3,er3,evt3=np.loadtxt(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\promedios 4.txt')
r4,vt4,er4,evt4=np.loadtxt(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\promedios 5.txt')
r5,vt5,er5,evt5=np.loadtxt(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\promedios 6.txt')


plt.figure()
plt.scatter(r0,vt0,label='Agitacion 1',color='cornflowerblue',zorder=1,s=10)
plt.scatter(r1,vt1,label='Agitacion 2',color='purple',zorder=1,s=10)
#plt.scatter(r2,vt2,label='Agitacion 3',color='brown',zorder=1,s=10)
plt.scatter(r3,vt3,label='Agitacion 3',color='green',zorder=1,s=10)
plt.scatter(r4,vt4,label='Agitacion 4',color='red',zorder=1,s=10)
#plt.scatter(r5,vt5,label='Agitacion 6',color='yellow',zorder=1,s=10)
plt.legend(loc='upper right')
plt.xlabel("Radios [cm]", fontsize=13)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 13)
plt.title('$V_θ$(r) - 27.2% - varias agitaciones')
plt.show() 

#%%

#Datos de ajuste de dif agitaciones

ajustes=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\Dato ajustes nuevo.txt')
rajuste= np.loadtxt(ajustes,delimiter= ',', skiprows= 0)

chi2s=rajuste[:,8]
chi2reds=rajuste[:,9]
gamma=2*np.pi*rajuste[:,11]
egamma=2*np.pi*rajuste[:,12]
rc=rajuste[:,13]
erc=rajuste[:,14]
r0=rajuste[:,17]
er0=rajuste[:,18]

rcs=np.array((rc[0],rc[1],rc[3],rc[4]))
gammas=np.array((gamma[0],gamma[1],gamma[3],gamma[4]))/(2*np.pi)

Omegas=gammas/(rcs**2)
n=np.array((0,1,2,3))

print(Omegas)

plt.figure()
plt.scatter(n,Omegas)
plt.show()

#%%

popta=np.array((rajuste[4,11]/(2*np.pi),rajuste[4,13],rajuste[4,15]))
popt27=np.array((rajuste[4,11]/(2*np.pi),rajuste[4,13],rajuste[4,15]))
popt50=np.array((rajuste[4,11]/(2*np.pi),rajuste[4,13],rajuste[4,15]))

#%%

#Diferentes viscosidades

r27,vt27,er27,evt27=np.loadtxt(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\promedios 1.txt')
ra,vta,era,evta=np.loadtxt(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\promedios agua.txt')
r50,vt50,er50,evt50=np.loadtxt(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\promedios 50.txt')

xdata=np.linspace(min(r27),6.5,len(r27)*20)

def func(x,a,b,c):
    return ((a/(x-c))*(1 - np.exp(-((x-c)**2/b**2))))

#esta medio cualqeuira esto, medir de nuevo ya fue 
plt.figure()
plt.errorbar(ra,vta,xerr=era,yerr=evta,label='Datos agua',color='r',zorder=1,fmt='.')
plt.errorbar(r27,vt27,xerr=er27,yerr=evt27,label='Glicerina 27.2 %',color='cornflowerblue',zorder=1,fmt='.')
plt.errorbar(r50,vt50,xerr=er50,yerr=evt50,label='Glicerina 50%',color='green',zorder=1,fmt='.')
plt.plot(xdata,func(xdata,*popta),label='Ajuste agua',color='b',zorder=2)
plt.legend(loc='upper right')
plt.xlabel("Radios [cm]", fontsize=13)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 13)
plt.title('$V_θ$(r) - varias soluciones')
plt.xlim(-0.3,7)
plt.show() 

#%%

#analizo vorticidad

vortice0=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\PIVlab.vorticidad.27.2.1.txt')
matrizV0 = np.loadtxt(vortice0,delimiter= ',', skiprows= 3, usecols = (0,1,2,3))

Vx0,Vy0=matrizV0[:,0], matrizV0[:,1] #posicion de vx y vy
Vvx0,Vvy0=matrizV0[:,2], matrizV0[:,3] #valores de vx y vy

plt.figure()
plt.quiver(Vx0,Vy0,Vvx0,Vvy0)

#%%

centro = plt.ginput(1)
x0,y0 = np.array(centro[0])

#%%

#Calibracion espacial

#Calibracion dia 2
#escp=29.93713202275223/5*10 #px por cm
#escp_error=0.13440475595279747*10

#Calibracion dia 3: dos vortices
#escp=21.55473505638159
#escp_error=0.23189202751466623

#Calibracion dia 3: sumidero
escp=31.34482526881721
escp_error=0.30160302831451574

#Conversion a polares

from Analisis_fluidos import converTo_polares
from Analisis_fluidos import homogVtan_R

Mxy,eMxy=np.array([x,y,vx,vy]),np.array([ex,ey,evx,evy])
M_xy,eM_xy=Mxy.T, eMxy.T

M_rt,eM_rt = converTo_polares(M_xy,eM_xy)

#Promedios

matrizProm_Rp_Vtp,matrizProm_Rp_Vtp_sd = homogVtan_R(M_rt, 50)
matrizProm_Rp_Vtp_ev,matrizProm_Rp_Vtp_ev_sd = homogVtan_R(eM_rt, 50)
matrizProm_Rp_Vrp,matrizProm_Rp_Vrp_sd = homogVr_R(M_rt, 50)

#%%

plt.figure()
plt.plot(M_rt[:,0], M_rt[:,3], '.',label='Vorticidad',color='b')
plt.plot(matrizProm_Rp_Vtp[:,0],matrizProm_Rp_Vtp[:,1], '.', label='Promedio vorticidad',color='r')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Vorticidad [1/s]", fontsize = 15)
plt.title('w(r) - {}'.format(archivo))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc='upper right')
plt.show()

#%%

def vortice(x,a,b,c):
    return (a/(np.pi*b**2))*np.exp(-((x-c)**2/b**2))

poptv, pcovv = curve_fit(vortice, x, y,sigma= ,absolute_sigma=True)
perrv=np.diag(pcovv)
rdatav = np.linspace(min( ),max( ),len( ))

plt.figure()
plt.title('w(r) promedio - {}'.format(archivo))
#plt.plot(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], '.',label='datos')
plt.errorbar(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], xerr = matrizProm_Rp_Vtp_sd[:,0], yerr =matrizProm_Rp_Vtp_sd[:,1],fmt='.',label='Datos',zorder=1,color='b')
plt.plot(rdatav,vortice(rdatav,*poptv), label='Ajuste',zorder=2,color='r')
plt.xlabel("Radios [cm]", fontsize=14)
plt.ylabel("Vorticidad [1/s]", fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()


def chi2red(func,x,popt,y,std):
    fitted_function=func(x,*popt)
    dof=len(y)-len(popt)
    chisquare=sum((y-fitted_function)**2/(std**2))
    chi2red=chisquare/dof
    return chi2r


print('Chi-square reducido es', chi2r)
print('Gamma es', poptv[0],'+-', perrv[0])
print('rc o sqrt(4nu/alpha) es', poptv[1],'+-', perrv[1])
print('r0 es', poptv[2],'+-', perrv[2])
