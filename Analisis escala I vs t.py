# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 12:28:14 2022

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

#CONVERSION ESCALA TEMP:CORRIENTE

ventanaG=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 1\ventas barrido I-22C\Ventana-1.6mA-0.1kHz-1ms-50mV.txt')
ventG=np.loadtxt(ventanaG)

tG=ventG[0]
VG=ventG[1]

plt.figure()
plt.plot(tG,VG)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

puntos = plt.ginput(4)
puntos = np.array(puntos)

inter_t=puntos[1][0]-puntos[0][0] #=0.01 s ok
inter_amp=puntos[3][1]-puntos[2][1] #amplitud_V=0.11 V o  110 mV, o sea Vpp=0.22 V

Esc_tiempos=np.array([puntos[2][0],puntos[3][0]])
Esc_volt=np.array([puntos[2][1],puntos[3][1]])
Esc_amp=np.array([0,0.016])

plt.figure()
plt.plot(tG,VG,zorder=1)
plt.scatter(Esc_tiempos,Esc_volt,color='r',zorder=2)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.title('Señal 1.6 mA')
plt.show()


def f(a,b,x):
    return a*x+b

#Intensidad=A*tiempo+B 
#A=intensidad/tiempo

poptI, pcovI = curve_fit(f, Esc_tiempos,Esc_amp)
perrI = np.sqrt(np.diag(pcovI))

escala_I=tG*poptI[1]+poptI[0]

xd=np.linspace(-0.0025,0.0027,200)

plt.figure()
plt.scatter(Esc_tiempos,Esc_amp)
plt.plot(xd,xd*poptI[1]+poptI[0],color='g')
plt.xlabel('Tiempo (s)')
plt.ylabel('Intensidad (A)')
plt.show()


plt.figure()
plt.plot(escala_I,VG)
#plt.scatter(I_abs,V_abs)
plt.show()

picos = plt.ginput(5)
picos = np.array(picos)

I_abs=np.array([picos[1][0],picos[2][0],picos[3][0],picos[4][0]])
V_absG=np.array([picos[1][1],picos[2][1],picos[3][1],picos[4][1]])
'''
La intensidad de abs son (sacado de la escala de 1.6 mA)

array([0.00604866, 0.00664407, 0.00732067, 0.00761837])

'''

plt.figure()
plt.plot(escala_I,VG,label='Datos medidos',zorder=1)
plt.scatter(I_abs,V_absG,color='r',label='Picos de Intensidades encontrados',zorder=2)
plt.ylabel('Voltaje (V)')
plt.xlabel('Intensidad (A)')
plt.legend(loc='upper left')
plt.show()

#%%

R87=np.array([377.104391,377.105205,377.111226,377.112040])
#R87 = np.array([376.70018737600003, 377.514687376, 376.6993126240001, 377.513812624])
R85=np.array([377.105909,377.106274,377.108947,377.109307])

w_abs=np.array([R87[3],R85[3],R85[0],R87[0]])+0.014

escala_wG=tG*popt[1]+popt[0]

#Paso a escala de frecuencia usando la escala que ya conozco de tiempo vs frecuencia

plt.figure()
plt.plot(escala_wG,VG)
plt.axvline(w_abs[0],linestyle='dashed',color='g')
plt.axvline(w_abs[1],linestyle='dashed',color='g')
plt.axvline(w_abs[2],linestyle='dashed',color='g')
plt.axvline(w_abs[3],linestyle='dashed',color='g')
plt.show()

#Intensidad=frecuencia*A+B

poptIf, pcovIf = curve_fit(f, I_abs,w_abs)
perrIf = np.sqrt(np.diag(pcovIf))

escala_If=escala_I*poptIf[1]+poptIf[0]

xdata=np.linspace(0.0055,0.008,200)

plt.figure()
plt.scatter(I_abs,w_abs, label='Frec. de absorción teóricas\n     vs     \nIntensidades de Picos encontrados')
plt.plot(xdata, xdata*poptIf[1]+poptIf[0],label='Ajuste',color='g')
plt.legend(loc='upper right')
plt.ylabel('Frecuecia del laser (THz)')
plt.xlabel('Intensidad (A)')
plt.plot()


plt.figure()
plt.plot(escala_I,VG,label='Datos medidos',zorder=1)
plt.scatter((w_abs-poptIf[0])/poptIf[1],V_absG,label='Int. de absorción teóricas',zorder=2,color='r')
plt.title('Señal 1.6 mA')
plt.xlabel('Intensidad (A)')
plt.ylabel('Voltaje (V)')
plt.legend(loc='upper left')
plt.show()
