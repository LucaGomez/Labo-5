# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 21:48:13 2022

@author: luo
"""

#importo las librerias a utilizar
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sg
from scipy.optimize import curve_fit
from scipy import stats

from IPython import get_ipython

get_ipython().run_line_magic('matplotlib', 'qt5')
plt.rcParams.update({'font.size': 12})

#%%'

#ESCALA CON LA CORRIENTE DE 1.6mA

ventana=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 2\Ventana-2can-escala-22C-1.6mA-0.1kHz-250us-50mV.txt')
vent=np.loadtxt(ventana)

t=vent[0]
V=vent[1]

plt.figure()
plt.plot(t,V)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

#%%


#filtro la señal
V_filt = sg.savgol_filter(V, 21, 3)


#busco picos
pic_index = sg.find_peaks(-V, prominence=(0.005, None), distance=82)[0]
print(pic_index)

t[pic_index[3]]=0.000764
V_filt[pic_index[3]]=0.27482

plt.rcParams.update({'font.size': 9})


#grafico
plt.figure()
plt.plot(t,V,label='Señal registrada',zorder=1,color='green')
plt.plot(t,V_filt,'-',label='Señal filtrada',zorder=2,color='blue')
plt.scatter(t[pic_index], V_filt[pic_index],s=10,label='Picos encontrados',zorder=3,color='red')
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.title('Señal 1.6 mA')
#plt.axvline(t[pic_index[0]],color='r',linestyle='dashed')
#plt.axvline(t[pic_index[3]],color='r',linestyle='dashed')
#plt.axvline(t[pic_index[1]],color='g',linestyle='dashed')
#plt.axvline(t[pic_index[2]],color='g',linestyle='dashed')
#plt.legend(loc='upper left')
plt.show()


#%%

R87=np.array([377.104391,377.105205,377.111226,377.112040])
#R87 = np.array([376.70018737600003, 377.514687376, 376.6993126240001, 377.513812624])
R85=np.array([377.105909,377.106274,377.108947,377.109307])

#w_abs=np.array([R87[3],R85[3],R85[1],R87[1]])

w_abs=np.array([R87[3],R85[3],R85[0],R87[0]])

t_abs=np.array([t[pic_index[0]],t[pic_index[1]],t[pic_index[2]],0.000764])

V_abs=np.array([V_filt[pic_index[0]],V_filt[pic_index[1]],V_filt[pic_index[2]],0.27482])

def f(a,b,x):
    return a*x+b

#frecuencia=tiempo*A+B

popt, pcov = curve_fit(f, t_abs,w_abs)
perr = np.sqrt(np.diag(pcov))

escala_w=t*popt[1]+popt[0]

xdata=np.linspace(0.00025,0.00078,200)

plt.figure()
plt.scatter(t_abs,w_abs, label='$v_{abs}$ vs $t_{picos}$')
plt.plot(xdata, xdata*popt[1]+popt[0],label='Ajuste',color='g')
plt.legend(loc='upper right')
plt.ylabel('Frecuencia del láser (THz)')
plt.xlabel('Tiempo (s)')
plt.title('Ajuste lineal')
plt.plot()

plt.figure()
plt.plot(escala_w,V_filt)
#plt.scatter(w_abs,V_abs,color='r',label='Frec. de absorción teóricas')
plt.title('Señal 1.6 mA')
plt.xlabel('Frecuencia del láser (THz)')
plt.ylabel('Voltaje (V)')
plt.axvline(w_abs[0],color='r',linestyle='dashed',label='Rb-87: F=2 → F=1')
plt.axvline(w_abs[1],color='g',linestyle='dashed',label='Rb-85: F=3 → F=2')
plt.axvline(w_abs[2],color='g',linestyle='dashed',label='Rb-85: F=2 → F=3')
plt.axvline(w_abs[3],color='r',linestyle='dashed',label='Rb-87: F=1 → F=2')
plt.legend(loc='upper right')
plt.show()

dof=len(w_abs)-2

fitted_function=f(popt[1],popt[0],t_abs)
chi2stat=stats.chisquare(w_abs,fitted_function,dof)
p=1-stats.chi2.cdf(chi2stat,dof)

print('Chi-square es', chi2stat)


def R2(yd,yf):
    residuals = yd-yf
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((yd-np.mean(yd))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

coef_det=R2(w_abs,fitted_function)

print('La escala THz (frecuencia de laser):s (escala temp del osciloscopio) es de: ')
print(popt[1], '1/s^2')
print('Con la frecuencia corrida en')
print(popt[0], 'THz')
print('chi^2 es:', chi2stat)
print('p-value es', p)
print('R2 es:', coef_det)

#%%

#SOSPECHA DEL ULTIMO PICO

plt.figure()
plt.plot(escala_w,V_filt)
plt.scatter(w_abs[3],V_abs[3],color='r',label='Nivel excitado del fundamental, R87')
plt.scatter(w_abs[3],V_abs[3],color='r',label='Nivel excitado del fundamental, R87')
plt.title('Señal 1.6 mA')
plt.xlabel('Frecuecia del laser (THz)')
plt.ylabel('Voltaje (V)')
plt.show()

#%%

#PRUEBO CON CORRIENTE DE 1.7mA


ventana1=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 2\Ventana-2can-escala-22C-1.7mA-0.1kHz-250us-50mV.txt')
vent1=np.loadtxt(ventana1)

t1=vent1[0]
V1=vent1[1]

plt.figure()
plt.plot(t1,V1)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

#filtro la señal
V_filt1 = sg.savgol_filter(V1, 13, 2)

#busco picos
pic_index1 = sg.find_peaks(-V1, prominence=(0.0045, None), distance=105)[0]
print(pic_index1)

#grafico
plt.figure()
plt.plot(t1,V1,color='g',label='Señal registrada')
plt.plot(t1,V_filt1,'-',label='Señal filtrada')
plt.plot(t1[pic_index1], V_filt1[pic_index1], 'o',label='picos')
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.axvline(t1[pic_index1[0]],color='r',linestyle='dashed')
#plt.axvline(t1[pic_index1[3]],color='r',linestyle='dashed')
plt.axvline(t1[pic_index1[1]],color='g',linestyle='dashed')
plt.axvline(t1[pic_index1[2]],color='g',linestyle='dashed')
plt.legend(loc='upper left')
plt.show()

Cuartopico_t = 0.000779
Cuartopico_v = 0.270797

#%%


t1_abs=np.array([t1[pic_index1[0]],t1[pic_index1[1]],t1[pic_index1[2]],Cuartopico_t])

V1_abs=np.array([V_filt1[pic_index1[0]],V_filt1[pic_index1[1]],V_filt1[pic_index1[2]],Cuartopico_v])

popt1, pcov1 = curve_fit(f, t1_abs,w_abs)
perr1 = np.sqrt(np.diag(pcov1))

escala1_w=t1*popt1[1]+popt1[0]

plt.figure()
plt.plot(escala1_w,V_filt1)
plt.scatter(w_abs,V1_abs)
plt.show()

print('La escala THz (frecuencia de laser):s (escala temp del osciloscopio) es de: ')
print(popt1[1], '1/s^2')
print('Con la frecuencia corrida en')
print(popt1[0], 'THz')


fitted_function1=f(popt1[1],popt1[0],t1_abs)
chi2stat1=stats.chisquare(w_abs,fitted_function1,dof)
coef_det1=R2(w_abs,fitted_function1)
p1=1-stats.chi2.cdf(chi2stat1,dof)

print('chi^2 es:', chi2stat1)
print('p-value es', p1)
print('R2 es:', coef_det1)

print('La escala THz (frecuencia de laser):s (escala temp del osciloscopio) es de: ')
print(popt[1], '1/s^2')
print('Con la frecuencia corrida en')
print(popt[0], 'THz')



#%%

#CORRIENTE 1.8mA: FUERA DE ESCALA YA

ventana2=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 2\Ventana-2can-escala-22C-1.8mA-0.1kHz-250us-50mV.txt')
vent2=np.loadtxt(ventana2)

t2=vent2[2][0:1500]
V2=vent2[3][0:1500]

plt.figure()
plt.plot(t2,V2)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

#filtro la señal
V_filt2 = sg.savgol_filter(V2, 17, 2)

#busco picos
pic_index2 = sg.find_peaks(-V2, prominence=(0.009, None), distance=100)[0]
print(pic_index2)


#grafico
plt.figure()
plt.plot(t2,V2,color='g',label='Señal registrada')
plt.plot(t2,V_filt2,'-',label='Señal filtrada')
plt.plot(t2[pic_index2], V_filt2[pic_index2], 'o',label='picos')
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.axvline(t2[pic_index2[0]],color='r',linestyle='dashed')
#plt.axvline(t1[pic_index1[3]],color='r',linestyle='dashed')
plt.axvline(t2[pic_index2[1]],color='g',linestyle='dashed')
plt.axvline(t2[pic_index2[2]],color='g',linestyle='dashed')
plt.legend(loc='upper left')
plt.show()

primerpico_t2 = 0.00051787
primerpico_v2 = 0.256

t2_abs=np.array([primerpico_t2, t2[pic_index2[0]],t2[pic_index2[1]],t2[pic_index2[2]]])

V2_abs=np.array([primerpico_v2, V_filt2[pic_index2[0]],V_filt2[pic_index2[1]],V_filt2[pic_index2[2]]])

popt2, pcov2 = curve_fit(f, t2_abs,w_abs)
perr2 = np.sqrt(np.diag(pcov2))

escala2_w=t2*popt2[1]+popt2[0]

plt.figure()
plt.plot(escala2_w,V_filt2)
plt.scatter(w_abs,V2_abs)
plt.show()

print('La escala THz (frecuencia de laser):s (escala temp del osciloscopio) es de: ')
print(popt2[1], '1/s^2')
print('Con la frecuencia corrida en')
print(popt2[0], 'THz')


fitted_function2=f(popt2[1],popt2[0],t2_abs)
chi2stat2=stats.chisquare(w_abs,fitted_function2,dof)
coef_det2=R2(w_abs,fitted_function2)
p2=1-stats.chi2.cdf(chi2stat2,dof)

print('chi^2 es:', chi2stat2)
print('p-value es', p2)
print('R2 es:', coef_det2)

#%%


#CORRIENTE 1.5mA

ventana3=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 2\Ventana-2can-escala-22C-1.5mA-0.1kHz-250us-50mV.txt')
vent3=np.loadtxt(ventana3)

t3=vent3[0]
V3=vent3[1]

plt.figure()
plt.plot(t3,V3)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

#filtro la señal
V_filt3 = sg.savgol_filter(V3, 13, 2)

#busco picos
pic_index3 = sg.find_peaks(-V3, prominence=(0.005, None), distance=80)[0]
print(pic_index3)

w_abs=np.array([R87[3],R85[3],R85[1],R87[1]])

#w_abs=np.array([R87[2],R85[2],R85[0],R87[0]])

t3_abs=t3[pic_index3]
V3_abs=V_filt3[pic_index3]

#t3_abs=np.array([primerpico_t2, t2[pic_index2[0]],t2[pic_index2[1]],t2[pic_index2[2]]])
#V3_abs=np.array([primerpico_v2, V_filt2[pic_index2[0]],V_filt2[pic_index2[1]],V_filt2[pic_index2[2]]])

popt3, pcov3 = curve_fit(f, t3_abs,w_abs)
perr3 = np.sqrt(np.diag(pcov3))

escala3_w=t3*popt3[1]+popt3[0]

plt.figure()
plt.plot(escala3_w,V_filt3)
plt.scatter(w_abs,V3_abs)
plt.show()

print('La escala THz (frecuencia de laser):s (escala temp del osciloscopio) es de: ')
print(popt3[1], '1/s^2')
print('Con la frecuencia corrida en')
print(popt3[0], 'THz')


fitted_function3=f(popt3[1],popt3[0],t3_abs)
chi2stat3=stats.chisquare(w_abs,fitted_function3,dof)
coef_det3=R2(w_abs,fitted_function3)
p3=1-stats.chi2.cdf(chi2stat3,dof)

print('chi^2 es:', chi2stat3)
print('p-value es', p3)
print('R2 es:', coef_det3)

#%%


#CORRIENTE 1.4mA

ventana4=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 2\Ventana-2can-escala-22C-1.4mA-0.1kHz-250us-50mV.txt')
vent4=np.loadtxt(ventana4)

t4=vent4[0]
V4=vent4[1]

plt.figure()
plt.plot(t4,V4)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

#filtro la señal
V_filt4 = sg.savgol_filter(V4, 13, 2)

#busco picos
pic_index4 = sg.find_peaks(-V4, prominence=(0.005, None), distance=100)[0]
print(pic_index4)

plt.figure()
plt.plot(t4,V4,color='g',label='Señal registrada')
plt.plot(t4,V_filt4,'-',label='Señal filtrada')
plt.plot(t4[pic_index4], V_filt4[pic_index4], 'o',label='picos')
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

Cuartopico_t4 = 0.000541
Cuartopico_v4 = 0.278

t4_abs=np.array([t4[pic_index4[0]],t4[pic_index4[1]],t4[pic_index4[2]],Cuartopico_t4])
V4_abs=np.array([V_filt4[pic_index4[0]],V_filt4[pic_index4[1]],V_filt4[pic_index4[2]],Cuartopico_v4])

popt4, pcov4 = curve_fit(f, t4_abs,w_abs)
perr4 = np.sqrt(np.diag(pcov4))

escala4_w=t4*popt4[1]+popt4[0]

plt.figure()
#plt.plot(escala4_w,V_filt4)
plt.plot(escala4_w,V4)
plt.scatter(w_abs,V4_abs)
plt.show()

print('La escala THz (frecuencia de laser):s (escala temp del osciloscopio) es de: ')
print(popt4[1], '1/s^2')
print('Con la frecuencia corrida en')
print(popt4[0], 'THz')


fitted_function4=f(popt4[1],popt4[0],t4_abs)
chi2stat4=stats.chisquare(w_abs,fitted_function4,dof)
coef_det4=R2(w_abs,fitted_function4)
p4=1-stats.chi2.cdf(chi2stat4,dof)

print('chi^2 es:', chi2stat4)
print('p-value es', p4)
print('R2 es:', coef_det4)


#%%


#CORRIENTE 1.3mA

ventana5=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 2\Ventana-2can-escala-22C-1.3mA-0.1kHz-250us-50mV.txt')
vent5=np.loadtxt(ventana5)

t5=vent5[0]
V5=vent5[1]

plt.figure()
plt.plot(t5,V5)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

#filtro la señal
V_filt5 = sg.savgol_filter(V5, 13, 2)

#busco picos
pic_index5 = sg.find_peaks(-V5, prominence=(0.005, None), distance=100)[0]
print(pic_index5)


plt.figure()
plt.plot(t5,V5,color='g',label='Señal registrada')
plt.plot(t5,V_filt5,'-',label='Señal filtrada')
plt.plot(t5[pic_index5], V_filt5[pic_index5], 'o',label='picos')
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

#%%
Cuartopico_t5 = 0.000648
Cuartopico_v5 = 0.29

t5_abs=np.array([t5[pic_index5[0]],t5[pic_index5[1]],t5[pic_index5[2]],Cuartopico_t5])
V5_abs=np.array([V_filt5[pic_index5[0]],V_filt5[pic_index5[1]],V_filt5[pic_index5[2]],Cuartopico_v5])

popt5, pcov5 = curve_fit(f, t5_abs,w_abs)
perr5 = np.sqrt(np.diag(pcov5))

escala5_w=t5*popt5[1]+popt5[0]

plt.figure()
#plt.plot(escala4_w,V_filt4)
plt.plot(escala5_w,V5,zorder=1)
plt.scatter(w_abs,V5_abs,color='r',zorder=2)
plt.show()

print('La escala THz (frecuencia de laser):s (escala temp del osciloscopio) es de: ')
print(popt5[1], '1/s^2')
print('Con la frecuencia corrida en')
print(popt5[0], 'THz')


fitted_function5=f(popt5[1],popt5[0],t5_abs)
chi2stat5=stats.chisquare(w_abs,fitted_function5,dof)
coef_det5=R2(w_abs,fitted_function5)
p5=1-stats.chi2.cdf(chi2stat5,dof)

print('chi^2 es:', chi2stat5)
print('p-value es', p5)
print('R2 es:', coef_det5)


#%%


#CORRIENTE 1.2mA

ventana6=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Dia 2\Ventana-2can-escala-22C-1.2mA-0.1kHz-500us-50mV.txt')
vent6=np.loadtxt(ventana6)

t6=vent6[0]
V6=vent6[1]

plt.figure()
plt.plot(t6,V6)
plt.xlabel('Tiempo (s)')
plt.ylabel('Voltaje (V)')
plt.show()

#%%

plt.figure()
#plt.plot(t,V,label='1.6 mA')
plt.plot(escala1_w,V1,label='1.7 mA',zorder=1)
plt.scatter(w_abs,V1_abs,color='r',zorder=2)
#plt.plot(t2,V2,label='1.8 mA')
#plt.plot(escala3_w,V3,label='1.5 mA',zorder=1)
#plt.scatter(w_abs,V3_abs,zorder=2)
#plt.plot(t4,V4,label='1.4 mA')
plt.plot(escala5_w,V5,label='1.3 mA',zorder=1)
plt.scatter(w_abs,V5_abs,color='g',zorder=2)
#plt.plot(t6,V6,label='1.2 mA')
plt.legend(loc='upper right')
plt.show()

#Excluyendo las de 1.8 mA y 1.2 mA

Amplitudes=np.array([popt1[1],popt[1],popt3[1],popt4[1],popt5[1]])
Errores_Amp=np.array([perr1[1],perr[1],perr3[1],perr4[1],perr5[1]])

Offset=np.array([popt1[0],popt[0],popt3[0],popt4[0],popt5[0]])
Errores_Off=np.array([perr1[0],perr[0],perr3[0],perr4[0],perr5[0]])
Corrientes=np.array([1.7, 1.6, 1.5,1.4,1.3])

Chi2s=np.array([chi2stat1[0],chi2stat[0],chi2stat3[0],chi2stat4[0],chi2stat5[0]])
ps=np.array([p,p1,p3,p4,p5])
R2s=np.array([coef_det1,coef_det,coef_det3,coef_det4,coef_det5])

print('De las corrientes:', Corrientes , 'mA, obtengo que')
print('Las amplitudes de la escala son: ', Amplitudes, '+-', Errores_Amp)
print('Los Offsets de la escala son: ', Offset, '+-', Errores_Off)

'''
Se definio los picos sabiendo que existen dos picos de R87 y dos de R85 siendo los mas chicos el de R87 por tener menor abundancia
Usando la estructura hiperfina se define que los dos picos del respectivo isotopo son el estado fundamental fundamental y el exitado
En el cado de R87 esta transicion tiene  una frecuencia de 3.035 GHz y para  el R85 esta esta transicion es resonancia a una desintonia de 6.834 GHz
De esta forma mirando la tabla de los 4 niveles de energia para cada isotopo y definiendo una relacion lineal de frecuencia vs absorcion, se nota que la opcion con mas probabilidad es la siguiente

De esta forma puedo hacer la conversion frecuencia del laser: tiempo medido como


Parece convenir usar la de 1.5 mA porque es la que lleva menos error, y tambien fue mas facil de detectar el cuarto pico por algun motivo
Pero su escala difiere un poco del resto, asi que capaz conviene usar la de 1.6 mA para estudiar las cosas

Como conozco la relacion tiempo medido vs Amplitud de la modulacion, puedo sacar la escala frecuencia de laser vs intensidad de corriente

Una vez obtenido eso puedo analizar las frecuencias (ver si lo hace luca)

Y despues a analizar el campo magnetico y DAVS

'''

#%%


data_to_save = np.vstack((Amplitudes, Errores_Amp,Offset,Errores_Off,Chi2s,R2s))

np.savetxt(r'C:\Users\luo\OneDrive\文档\Labo 5\Espectrocopia laser\Resultados escala-frecuencia laser vs tiempo osci.txt', data_to_save)



'''


#%%

#cargo las frecuencias de transicion del rubidio (isotopos 85 y 87)
R87=np.array([377.104391,377.105205,377.111226,377.112040])
#R87 = np.array([376.70018737600003, 377.514687376, 376.6993126240001, 377.513812624])
R85=np.array([377.105909,377.106274,377.108947,377.109307])

#Frecuencia de excitacion estado fundamental al excitado de R85=3.035*0.001
#Frecuencia de excitacion estado fundamental al excitado de R87=6.834*0.001
print(R87[3]-R87[1])
print(R85[3]-R85[1])



#%%

variantes1 = [0,1,2,3]
variantes2 = [0,1,2,3,4,5, 6, 7]

def buscarCercania(fs_posibles, fs_objetivo):
    result = []
    for f_obj in fs_objetivo:
        result.append(min(abs(f_obj*np.ones(len(fs_posibles)) - fs_posibles)))
    return np.array(result)

def Rescala(a,b,x):
        t=a*x+b
        return t
            
def escalarPicos(i,j,k,l):
    FsTotales = np.concatenate([R85, R87])
    
    T1=t[pic_index[i]]
    T2=t[pic_index[j]]
    F1=FsTotales[k]
    F2=FsTotales[l]
    
    
    a=(T2-T1)/(F2-F1)
    b=T1-a*F1
    
    
    FR87=Rescala(a,b,R87)
    FR85=Rescala(a,b,R85)
    #grafico
    return (FR85, FR87)       


for i in variantes1:
    for j in variantes1:
        for k in variantes2:
            for l in variantes2:
                
                #hago la escala
                FR85, FR87 = escalarPicos(i,j,k,l)
                fs_posibles = np.concatenate([FR85, FR87])
                fs_objetivo = t[pic_index]
                
                #si matchea lo grafico
                min_band = buscarCercania(fs_posibles, fs_objetivo)
                
                
                if((min_band < 0.000005).all()):
                    print(min_band)
                    print('pico ', i, 'match con frecuencia', F1 )
                    print('pico ', j, 'match con frecuencia', F2 )
                    #grafico
                    plt.figure()
                    plt.plot(t,V_filt,'-')
                    plt.plot(t[pic_index], V_filt[pic_index], 'o')
    
                    for i,f in enumerate(FR87):
                        plt.plot(np.ones(100)*f,np.linspace(min(V_filt), max(V_filt), 100), '--g')
                        plt.plot(np.ones(100)*FR85[i],np.linspace(min(V_filt), max(V_filt), 100), '--r')
 '''