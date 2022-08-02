# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import pyvisa
from scipy.signal import savgol_filter,argrelextrema
import time

from osciloscope import Osciloscope

rm = pyvisa.ResourceManager('@py')
instrumentos = rm.list_resources()
print(instrumentos)

osc = Osciloscope(0)

#%%
#LLevarnos algunas ventanas crudas con los parametros fijos

t, V = osc.getWindow(1)
data_to_save = np.vstack((t, V))
folderpath = '/home/publico/Desktop/temp-G6 L5/dia 3/Bose Eisten/'
#np.savetxt(folderpath  + 'Ventana-1k-sin-250ns-20mV-1050V.txt', data_to_save)

plt.figure()
plt.plot(t,V)
plt.axhline(0.0055)
plt.xlabel('Tiempo')
plt.ylabel('Tension')

#%%
#Guardamos los mìnimos y los picos de todas las ventanas juntas
N=600
#ampCmin=[]
ampCpeak=[]
for i in range(N):
    try:
        t, V = osc.getWindow(1)
        #filtered_col=savgol_filter(V, 15, 3)
        #y_min=filtered_col[argrelextrema(filtered_col, np.less)[0]]
        #ampCmin.append(y_min)
        peaks = find_peaks(V,distance=2)
        amp=V[peaks[0]]
        for j in range(len(amp)):
            ampCpeak.append(amp[j])
            #np.savetxt('/home/publico/Desktop/Grupo 6/Conteo/Ventanas/'  + f'50-sin-50ns-5mv-900V-{i}.txt', amp, fmt = '%.3f')
        print(i)
        time.sleep(50e-3)
    except:
        print('cague')
        pass
    
#np.savetxt('/home/publico/Desktop/Grupo 6/Conteo/Ventanas/'  + f'MinCf-50-sin-100ns-2mvBw-900V.txt', ampCmin, fmt = '%.7f')
np.savetxt('/home/publico/Desktop/temp-G6 L5/dia 3/Bose Eisten/'  + f'umbral-dist2-1k-con-25us-10mV-1050V.txt', ampCpeak, fmt = '%.7f')
#%%
#Probar los histogramas
binsize=np.arange(min(V), max(V), 0.0001)
plt.figure(1)
plt.hist(np.array(ampCpeak), bins=binsize)
plt.yscale('log')

#plt.figure(2)
#plt.hist(ampCmin)
#%%

t, V = osc.getWindow(1)
peaks = find_peaks(V,distance=2)
Amp=V[peaks[0]]
TAmp=t[peaks[0]]


plt.figure()
plt.plot(t,V, color='g')
plt.scatter(TAmp,Amp,color='r')


'''
filtered_col=savgol_filter(V, 15, 3)
y_min=filtered_col[argrelextrema(filtered_col, np.less)[0]]

plt.plot(t,V)
plt.plot(t,filtered_col)
'''

#%%
data=np.loadtxt('Frank.txt')
Amp=data
binsize=np.arange(min(Amp), max(Amp), 0.0002)
plt.figure(1)
plt.hist(Amp, bins=binsize)
plt.yscale('log')

#%%
#Guardamos los mìnimos y los picos de todas las ventanas juntas
N=500

Fotxvent=[]
fotxvent=0
for i in range(N):
    try:
        t, V = osc.getWindow(1)
        peaks = find_peaks(V,distance=2)
        amp=V[peaks[0]]
        for j in range(len(amp)):
            if amp[j] > 0.0055:
                fotxvent+=1
        Fotxvent.append(fotxvent)
        fotxvent=0
        print(i)
        time.sleep(50e-3)
    except:
        print(':c')
        pass
    
np.savetxt('/home/publico/Desktop/temp-G6 L5/dia 3/Bose Eisten/'  + f'ConteoBErapido-dist2-50-con-100ns-10mV-1050V.txt', Fotxvent, fmt = '%.7f')

#%%
mu=2.87
from scipy.stats import poisson
binsize=np.arange(min(Fotxvent), max(Fotxvent)+3,1)
plt.figure()
plt.hist(np.array(Fotxvent), bins=binsize, density=True)
#plt.yscale('log')
#plt.plot(Fotxvent, poisson.pmf(Fotxvent, mu), 'bo', ms=8, label='poisson pmf')

#%%
data=np.loadtxt('/home/publico/Desktop/temp-G6 L5/dia 2/Conteo-dist20-1k-sin-2.5us-10mV-1050V/frankesteindist20.txt')
Amp=data
binsize=np.arange(min(Fotxvent), max(Fotxvent)+1, 1)
plt.figure()
plt.hist(Amp, bins=binsize, density=True)
#%%
mu=24
from scipy.stats import poisson
n, bins, patches = plt.hist(Amp, bins=binsize,density=True)


# Scatter plot
# Now we find the center of each bin from the bin edges
bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(n))]
plt.figure()
plt.scatter(bins_mean, n, color='r')
plt.plot(data, poisson.pmf(data, mu), 'bo', ms=8, label='poisson pmf')

#%%
esc=1e-6 #1 microseg
t, V = osc.getWindow(1)
peaks = find_peaks(V,distance=2)
Numpic=len(peaks)


#%%


N=30

for i in range(N):
    t, V = osc.getWindow(1)
    data_to_save = np.vstack((t, V))
    folderpath = '/home/publico/Desktop/temp-G6 L5/dia 3/Bose Eisten/ModoIntensivo'
    np.savetxt(folderpath  + f'Ventana{i}RAPIDO-215k-con-1ms-500mV-1050V.txt', data_to_save)

#%%
data=np.loadtxt('/home/publico/Desktop/temp-G6 L5/dia 3/Bose Eisten/ModoIntensivo/ModoIntensivoVentana11-215k-con-5ms-2V-1050V.txt')
Amp=data[1]
t=data[0]

plt.figure()
plt.plot(t,Amp)

#%%


