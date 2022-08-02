# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 19:50:55 2022

@author: luo
"""


import matplotlib.pyplot as plt
import numpy as np
from IPython import get_ipython
import pandas as pd
from scipy.optimize import curve_fit
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import summary_table
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from decimal import Decimal



get_ipython().run_line_magic('matplotlib', 'qt5')
plt.rcParams.update({'font.size': 15})

#%%

#ANALISIS DE RESISTENCIA 1k: ventana

ventana=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Conteo\dia 3\1k\Ventana-1k-sin-250ns-20mV-1050V.txt')
data_ventana=np.loadtxt(ventana, delimiter=",") 

t_vent=data_ventana[0]
V_vent=data_ventana[1]

plt.figure()
plt.plot(t_vent,V_vent)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.title('Ventana del osciloscopio')
plt.show()

#%%

#ANALISIS DE RESISTENCIA 1k: ruido

ruido_dia2=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Conteo\dia 2\1k-apagado\Conteoapagado-dist20-1k-sin-2.5us-10mV-1050V.txt')
ruido_dia3=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Conteo\dia 3\1k\Conteoruido-dist20-1k-sin-2.5us-20mV-1050V.txt')

ruxvent_dia2=np.loadtxt(ruido_dia2, delimiter=",") 
ruxvent_dia3=np.loadtxt(ruido_dia3, delimiter=",") 

ruxvent=np.hstack((ruxvent_dia2,ruxvent_dia3))

binsize=np.arange(min(ruxvent)-0.5, max(ruxvent)+3-0.5, 1)

plt.figure()
n, bins, patches =plt.hist(ruxvent, bins=binsize,density=True,zorder=1)
Bins=bins[0:len(bins)-1]+0.5
plt.scatter(Bins,n,s=20,label='Datos',zorder=2,color='red')
plt.xlabel('Número de fotones por pantalla')
plt.ylabel('Frecuencia')
plt.title('Histograma: fuente apagada-resistencia 1k')
plt.show()

print(n)
efotxvent=(n[1]+2*n[2])/500


from scipy.special import factorial

def poisson(t, lamb,A):
    return A*(lamb**t)*np.exp(-lamb)/factorial(t)

popt, pcov = curve_fit(poisson, Bins,n)
perr = np.sqrt(np.diag(pcov))

px=np.linspace(0,3,200)

plt.figure()
plt.scatter(Bins,n,s=10,label='Datos',zorder=1,color='red')
plt.plot(px, poisson(px, popt[0],popt[1]), label='poisson pmf')
#plt.scatter(fxvent, poiss(fxvent, 20), '-', lw=1,label='Ajuste teórico',zorder=2,color='orange') #grafico de la recta de ajuste
plt.legend(loc='upper right')
plt.title('Distribución de Poisson: fuente apagada-resistencia 1k')
plt.xlabel('Número de fotones por pantalla')
plt.ylabel('Frecuencia')
plt.show()


print('el numero de fotones por microsegundo es :')
print(popt[0]/25,'+-', perr[0]/25 )
print('la amplitud de la distribucion de Poisson es :')
print(popt[1],'+-', perr[1] )

from scipy import stats

chi2stat=stats.chisquare(n,poisson(Bins,popt[0],popt[1]),len(n))

print(chi2stat)
dof=len(n)
p=1-stats.chi2.cdf(chi2stat,dof)

print('Chi-square es', chi2stat)
print('p-value es', p)

def R2(yd,yf):
    residuals = n-poisson(Bins,popt[0],popt[1])
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((n-np.mean(n))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

coef_det=R2(n,poisson(Bins,popt[0],popt[1]))
print("coef de determinacion I es",coef_det)

#%%

#ANALISIS DE RESISTENCIA 1k: conteo

conteo_dia3=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Conteo\dia 3\1k\Conteo-dist20-1k-sin-2.5us-20mV-1050V.txt')
fotxvent_dia3=np.loadtxt(conteo_dia3, delimiter=",") 

conteo_dia2=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Conteo\dia 2\1k-prendido\frankenstein conteo.txt')
fotxvent_dia2=np.loadtxt(conteo_dia2, delimiter=",") 

fotxvent=np.hstack((fotxvent_dia2,fotxvent_dia3))

binsize=np.arange(min(fotxvent)-0.5, max(fotxvent)+0.5, 1)

plt.figure()
n, bins, patches = plt.hist(fotxvent, bins=binsize,density=True,zorder=1)
fxvent=bins[0:len(bins)-1]+0.5
plt.scatter(fxvent, n, s=20, color='r',zorder=2)
plt.xlabel('Número de fotones por pantalla')
plt.ylabel('Frecuencia')
plt.title('Histograma: fuente prendida-resistencia 1k')
plt.show()


'''
from scipy.stats import poisson

def poiss(k, lamb):
    return poisson.pmf(k, lamb)

def poiss1(k, lamb, A):
    return A*poisson.pmf(k, lamb)

efxvent=efotxvent*n/len(fotxvent)

popt, pcov = curve_fit(poiss1, fxvent,n)
perr = np.sqrt(np.diag(pcov))

px=np.linspace(8,max(fxvent)+1,20)

plt.figure()
plt.errorbar(fxvent,n,yerr=efxvent,ls='none',label='Datos',zorder=1,color='black')
plt.scatter(fxvent,n,s=10,color='black')
plt.plot(fxvent, popt[1]*poisson.pmf(fxvent, popt[0]), label='poisson pmf')
#plt.scatter(fxvent, poiss(fxvent, 20), '-', lw=1,label='Ajuste teórico',zorder=2,color='orange') #grafico de la recta de ajuste
plt.legend(loc='upper right')
plt.title('Distribución de Poisson')
plt.xlabel('Número de fotones por pantalla')
plt.ylabel('Ocurrencia')
plt.show()
'''

#popt, pcov = curve_fit(poisson, fxvent,n)
popt, pcov = curve_fit(poisson, fxvent,n)
perr = np.sqrt(np.diag(pcov))

px=np.linspace(8,max(fxvent)+1,500)

plt.figure()
#plt.scatter(fxvent,n,s=10,label='Datos',zorder=1,color='red')
plt.scatter(fxvent,n,s=10,label='Datos',zorder=1,color='red')
plt.plot(px, poisson(px, popt[0],popt[1]), label='poisson pmf')
#plt.scatter(fxvent, poiss(fxvent, 20), '-', lw=1,label='Ajuste teórico',zorder=2,color='orange') #grafico de la recta de ajuste
plt.legend(loc='upper right')
plt.title('Distribución de Poisson: fuente prendida-resistencia 1k')
plt.xlabel('Número de fotones por pantalla')
plt.ylabel('Frecuencia')
plt.show()

print('el numero de fotones por microsegundo es :')
print(popt[0]/25,'+-', perr[0]/25 )
print('la amplitud de la distribucion de Poisson es :')
print(popt[1],'+-', perr[1] )



chi2stat=stats.chisquare(n,poisson(fxvent,popt[0],popt[1]),len(n))

print(chi2stat)
dof=len(n)
p=1-stats.chi2.cdf(chi2stat,dof)

print('Chi-square es', chi2stat)
print('p-value es', p)

def R2(yd,yf):
    residuals = n-poisson(fxvent,popt[0],popt[1])
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((n-np.mean(n))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

coef_det=R2(n,poisson(fxvent,popt[0],popt[1]))
print("coef de determinacion I es",coef_det)


#%%


residual=poisson(fxvent,popt[0],popt[1])-n


#%%

#ANALISIS DE 50OHM


conteo_50=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Conteo\dia 3\50\Conteo-dist20-50-sin-500ns-10mV-1050V.txt')
fotxvent_50=np.loadtxt(conteo_50, delimiter=",") 

conteo_51=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Conteo\dia 3\50\Conteobis-dist20-50-sin-500ns-10mV-1050V.txt')
fotxvent_51=np.loadtxt(conteo_51, delimiter=",") 

fotxvent50=np.hstack((fotxvent_50,fotxvent_51))

binsize=np.arange(min(fotxvent50)-0.5, max(fotxvent50)+0.5, 1)

plt.figure()
n, bins, patches = plt.hist(fotxvent50, bins=binsize,density=True,zorder=1)
fxvent50=bins[0:len(bins)-1]+0.5
plt.scatter(fxvent50, n, s=20, color='r',zorder=2)
plt.xlabel('Número de fotones por pantalla')
plt.ylabel('Frecuencia')
plt.title('Histograma: fuente prendida-resistencia 50')
plt.show()


'''
from scipy.stats import poisson

def poiss(k, lamb):
    return poisson.pmf(k, lamb)

def poiss1(k, lamb, A):
    return A*poisson.pmf(k, lamb)

efxvent=efotxvent*n/len(fotxvent)

popt, pcov = curve_fit(poiss1, fxvent,n)
perr = np.sqrt(np.diag(pcov))

px=np.linspace(8,max(fxvent)+1,20)

plt.figure()
plt.errorbar(fxvent,n,yerr=efxvent,ls='none',label='Datos',zorder=1,color='black')
plt.scatter(fxvent,n,s=10,color='black')
plt.plot(fxvent, popt[1]*poisson.pmf(fxvent, popt[0]), label='poisson pmf')
#plt.scatter(fxvent, poiss(fxvent, 20), '-', lw=1,label='Ajuste teórico',zorder=2,color='orange') #grafico de la recta de ajuste
plt.legend(loc='upper right')
plt.title('Distribución de Poisson')
plt.xlabel('Número de fotones por pantalla')
plt.ylabel('Ocurrencia')
plt.show()
'''

#popt, pcov = curve_fit(poisson, fxvent,n)
popt50, pcov50 = curve_fit(poisson, fxvent50,n)
perr50 = np.sqrt(np.diag(pcov50))

px=np.linspace(min(fxvent50),max(fxvent50)+1,500)

plt.figure()
#plt.scatter(fxvent,n,s=10,label='Datos',zorder=1,color='red')
plt.scatter(fxvent50,n,s=10,label='Datos',zorder=1,color='red')
plt.plot(px, poisson(px, popt50[0],popt50[1]), label='poisson pmf')
#plt.scatter(fxvent, poiss(fxvent, 20), '-', lw=1,label='Ajuste teórico',zorder=2,color='orange') #grafico de la recta de ajuste
plt.legend(loc='upper right')
plt.title('Distribución de Poisson: fuente prendida-resistencia 50')
plt.xlabel('Número de fotones por pantalla')
plt.ylabel('Frecuencia')
plt.show()

print('el numero de fotones por microsegundo es :')
print(popt50[0]/0.5,'+-', perr50[0]/0.5 )
print('la amplitud de la distribucion de Poisson es :')
print(popt50[1],'+-', perr50[1] )

chi2stat50=stats.chisquare(n,poisson(fxvent50,popt50[0],popt50[1]),len(n))

print(chi2stat)
dof=len(n)
p=1-stats.chi2.cdf(chi2stat50,dof)

print('Chi-square es', chi2stat50)

def R2(yd,yf):
    residuals = n-poisson(fxvent50,popt50[0],popt50[1])
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((n-np.mean(n))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

coef_det50=R2(n,poisson(fxvent50,popt50[0],popt50[1]))
print("coef de determinacion I es",coef_det50)




#%%


fxt=fxvent/25
efxvent=efotxvent*n/25
etiempo=0.0025 #VER ESTO MEJOR

#popt, pcov = curve_fit(poisson, fxvent,n)
popt, pcov = curve_fit(poisson, fxt,n)
perr = np.sqrt(np.diag(pcov))

px=np.linspace(0.3,1.4,500)

plt.figure()
plt.errorbar(fxt,n,yerr=efxvent,xerr=etiempo,ls='none',label='Datos',zorder=1,color='black')
plt.scatter(fxt,n,s=10,color='black')
plt.plot(px, poisson(px, popt[0],popt[1]), label='poisson pmf')
#plt.scatter(fxvent, poiss(fxvent, 20), '-', lw=1,label='Ajuste teórico',zorder=2,color='orange') #grafico de la recta de ajuste
plt.legend(loc='upper right')
plt.title('Distribución de Poisson')
plt.xlabel('Número de fotones por segundo (1/s)')
plt.ylabel('Frecuencia')
plt.show()


from scipy import optimize


fitfunc = lambda p, x: p[0]*(p[1]**x)*np.exp(-p[1])/factorial(x) # Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
p0 = [1., 20.]
p1, success = optimize.leastsq(errfunc, p0[:], args=(fxvent, n))

#%%


def poisson1(t, parameters):
    A=parameters[0]
    lamb=parameters[1]
    return A*(lamb**t)*np.exp(-lamb)/factorial(t)


def get_residuals(parameters, ydata, xdata):
  theoretical_function = poisson1(xdata,parameters)
  residuals = np.abs(theoretical_function - ydata)
  return residuals

from scipy.optimize import curve_fit, least_squares
guess_parameters=[1,20]

#Performing the fit
res_lsq = least_squares(get_residuals, guess_parameters, args=(n, fxvent))
#Fit results
best_parameters = res_lsq['x']

# Calculamos la matriz de covarianza "pcov"
def calcular_cov(res,y_datos):
    U, S, V = np.linalg.svd(res.jac, full_matrices=False)
    threshold = np.finfo(float).eps * max(res.jac.shape) * S[0]
    S = S[S > threshold]
    V = V[:S.size]
    pcov = np.dot(V.T / S**2, V)

    s_sq = 2 * res.cost / (y_datos.size - res.x.size)
    pcov = pcov * s_sq
    return pcov

pcov = calcular_cov(res_lsq,n)
# De la matriz de covarinza podemos obtener los valores de desviación estándar
# de los parametros hallados
pstd =np.sqrt(np.diag(pcov))

print('Parámetros hallados (con incertezas):')
   
print('Best A: ',round(best_parameters[0],3),'  ± ',round( pstd[0],3))
print('Best lamb: ',round(best_parameters[1],3),'  ± ',round( pstd[1],3))

#%%


