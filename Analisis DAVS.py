# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 21:32:46 2022

@author: luo
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sg
from scipy.optimize import curve_fit
from scipy import stats

from IPython import get_ipython

get_ipython().run_line_magic('matplotlib', 'qt5')
plt.rcParams.update({'font.size': 15})

#%%

ventana=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 2\Con campo\Ventana-DAVSPARAB-conB-22C-1.6mA-0.1kHz-250us-200mV.txt')
vent=np.loadtxt(ventana)

t=vent[0]
V=vent[1]

plt.figure()
plt.plot(t,V)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

#%%
escala=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Resultados escala-frecuencia laser vs tiempo osci.txt')
esc=np.loadtxt(escala)

A=esc[:1][0][1]
B=esc[2:3][0][1]

w=A*t+B

R87=np.array([377.104391,377.105205,377.111226,377.112040])
R85=np.array([377.105909,377.106274,377.108947,377.109307])

w_abs=np.array([R87[3],R85[3],R85[0],R87[0]])

plt.figure()
plt.plot(w,V)
plt.ylabel('Voltaje (V)')
plt.xlabel('Frecuencia (THz)')
plt.axvline(w_abs[0],linestyle='dashed',color='r',label='R87 excitado')
plt.axvline(w_abs[1],linestyle='dashed',color='r',label='R85 excitado')
plt.axvline(w_abs[2],linestyle='dashed',color='r',label='R85 fund')
plt.axvline(w_abs[3],linestyle='dashed',color='r',label='R87 fund')
plt.legend(loc='upper left')
plt.show()

puntos = plt.ginput(2)
puntos = np.array(puntos)

dif=puntos[1][0]-puntos[0][0]

#%%

w_esc=A*t+B-dif

plt.figure()
plt.plot(w_esc,V)
plt.ylabel('Voltaje (V)')
plt.xlabel('Frecuencia (THz)')
plt.axvline(w_abs[0],linestyle='dashed',color='r',label='R87 excitado')
plt.axvline(w_abs[1],linestyle='dashed',color='r',label='R85 excitado')
plt.axvline(w_abs[2],linestyle='dashed',color='r',label='R85 fund')
plt.axvline(w_abs[3],linestyle='dashed',color='r',label='R87 fund')
plt.legend(loc='upper left')
plt.show()

#%%

DAVS=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 2\Con campo\Ventana-DAVS-conB-22C-1.6mA-0.1kHz-250us-200mV.txt')
davs=np.loadtxt(DAVS)

LOC=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 2\Con campo\Ventana-DAVS-conB-22C-1.6mA-0.1kHz-250us-200mV.txt')
loc=np.loadtxt(LOC)


DAVSsinB=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 2\Sin campo\Ventana-DAVS-sinB-22C-1.6mA-0.1kHz-250us-50mV.txt')
davss=np.loadtxt(DAVSsinB)

tdavs=davs[0]
Vdavs=davs[1]

tdavss=davss[0]
Vdavss=davss[1]/4


plt.figure()
plt.plot(tdavs,Vdavs)
plt.plot(tdavss,Vdavss)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

puntos = plt.ginput(2)
puntos = np.array(puntos)


#%%

wdavs=A*tdavs+B-dif
wdavss=A*tdavss+B-dif
pt0,pt1=A*puntos[0][0]+B-dif,puntos[0][1]+0.2
pf0,pf1=A*puntos[1][0]+B-dif,puntos[1][1]+0.2

plt.figure()
plt.plot(wdavs,Vdavs+0.2,label='Con campo',zorder=2)
plt.plot(wdavss,Vdavss,label='Sin campo',zorder=1)
plt.scatter(pt0,pt1,zorder=3,color='r')
plt.scatter(pf0,pf1,zorder=3,color='r')
plt.ylabel('Voltaje (V)')
plt.xlabel('Frecuencia (THz)')
#plt.axvline(w_abs[0],color='r',linestyle='dashed',label='Rb-87: F=2 → F=1')
#plt.axvline(w_abs[1],color='g',linestyle='dashed',label='Rb-85: F=3 → F=2')
#plt.axvline(w_abs[2],color='g',linestyle='dashed',label='Rb-85: F=2 → F=3')
#plt.axvline(w_abs[3],color='r',linestyle='dashed',label='Rb-87: F=1 → F=2')
plt.title('Señal DAVS')
plt.legend(loc='upper left')
plt.show()

#%%


#BUSCO DOS PICOS PARA DEFINIR LA DEIF DE ENERGIA

plt.figure()
plt.plot(wdavs[1100:2000],Vdavs[1100:2000])
plt.ylabel('Voltaje (V)')
plt.xlabel('Frecuencia (THz)')
plt.show()

puntos_davs = plt.ginput(2)
puntos_davs = np.array(puntos_davs)

#%%

plt.figure()
plt.plot(wdavs,Vdavs)
plt.scatter(puntos_davs[0][0],puntos_davs[0][1],color='r')
plt.scatter(puntos_davs[1][0],puntos_davs[1][1],color='r')
plt.ylabel('Voltaje (V)')
plt.xlabel('Frecuencia (THz)')
plt.show()

deltav_zeeman=puntos_davs[1][0]-puntos_davs[0][0]

deltaw_zeeman=2*np.pi*deltav_zeeman*1e12 #=larmor frequency = gamma.B where gamma=g.uB/hslash

g=2
uB=(9.274)*1e-24 #J/Tesla
hslash=(1.05457)*1e-34 #J.s = J/Hz

#w=μB*B/hslash

B=deltaw_zeeman*hslash/uB

print(B)



#%%

puntos1_davs = plt.ginput(2)
puntos1_davs = np.array(puntos1_davs)

#%%

plt.figure()
plt.plot(wdavs,Vdavs)
plt.scatter(puntos1_davs[0][0],puntos1_davs[0][1],color='r')
plt.scatter(puntos1_davs[1][0],puntos1_davs[1][1],color='r')
plt.ylabel('Voltaje (V)')
plt.xlabel('Frecuencia (THz)')
plt.show()

deltav_zeeman1=puntos1_davs[1][0]-puntos1_davs[0][0]

B1=2*np.pi*1e12*deltav_zeeman1*hslash/uB

print(B1)


#%%

puntos2_davs = plt.ginput(2)
puntos2_davs = np.array(puntos2_davs)

plt.figure()
plt.plot(wdavs,Vdavs)
plt.scatter(puntos2_davs[0][0],puntos2_davs[0][1],color='r')
plt.scatter(puntos2_davs[1][0],puntos2_davs[1][1],color='r')
plt.ylabel('Voltaje (V)')
plt.xlabel('Frecuencia (THz)')
plt.show()

deltav_zeeman2=puntos2_davs[1][0]-puntos2_davs[0][0]

B2=2*np.pi*1e12*deltav_zeeman2*hslash/uB

print(B2)


#%%

puntos3_davs = plt.ginput(2)
puntos3_davs = np.array(puntos3_davs)

plt.figure()
plt.plot(wdavs,Vdavs)
plt.scatter(puntos3_davs[0][0],puntos3_davs[0][1],color='r')
plt.scatter(puntos3_davs[1][0],puntos3_davs[1][1],color='r')
plt.ylabel('Voltaje (V)')
plt.xlabel('Frecuencia (THz)')
plt.show()

deltav_zeeman3=puntos3_davs[1][0]-puntos3_davs[0][0]

B3=2*np.pi*1e12*deltav_zeeman3*hslash/uB

print(B3)

Campos=np.array((B,B1,B2,B3))

deltavs=np.array((deltav_zeeman,deltav_zeeman1,deltav_zeeman2,deltav_zeeman3))*1e12
