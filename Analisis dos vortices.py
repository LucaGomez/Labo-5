# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 20:24:48 2022

@author: luo
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage,signal
from IPython import get_ipython
from scipy.io import loadmat
from scipy.optimize import curve_fit
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from scipy import ndimage

from IPython import get_ipython

get_ipython().run_line_magic('matplotlib', 'qt5')
plt.rcParams.update({'font.size': 14})

#%%


#Calibracion dia 3: sumidero
escp=31.34482526881721
escp_error=0.30160302831451574

#Conversion a polares

def converTo_polares(MXY):
    MRT= np.empty((len(MXY), 4), float)
    for i in range(len(MXY)):
        x = MXY[i][0]
        y = MXY[i][1]
        u = MXY[i][2]
        v = MXY[i][3]
        
        #el modulo de (x,y)
        r = np.sqrt(x**2 + y**2)
        
        #el angulo en radianes de (x,y)
        t = np.angle(x + y*1j)
        
        MRT[i][0] = r
        MRT[i][1] = t
        
        #defino las velocidades como 0 inicialmente
        vrad = 0
        vtan = 0
        
        #hay que tener cuidado con la posicion de la matriz que esta en el origen, pues ahi r = 0
        if r != 0 :
            #vel radial
            vrad = (1/r)*(x*u + y*v) 
            #vel tangencial r*(titapunto)
            vtan = (1/r)*(x*v - y*u)
            
        MRT[i][2] = vrad
        MRT[i][3] = vtan
    return MRT


#promedio y desviacion estandar de la dispersion 

def homogVtan_R(MRT, p):
    heightMRT = len(MRT)
    #MRT es la matriz que contiene en la primer y cuarta columna los datos a promediar
    #p es la cantidad de datos que quiero obtener del promediado, si p->inf se reobtiene los datos originales
    
    #cuantos puntos quiero para la curva?, se define al definir la cantidad de bloques en los que se promedia
    blockBorders = np.linspace(min(MRT[:,0]), max(MRT[:,0]), p)
    #se devuelve esta lista transformada a array
    r_prom_Vtan = []
    r_prom_Vtan_ds=[]
    
    for j in range(p):
        rs_toprom = []
        vt_toprom = []
        for i in range(heightMRT):
            if j!= p - 1 and blockBorders[j] < MRT[i][0] < blockBorders[j+1]:     
                vt_toprom.append(MRT[i][3])
                rs_toprom.append(MRT[i][0])     
            #print(blockBorders[j])
            if j == p - 1 and blockBorders[j-2] < MRT[i][0] :
                vt_toprom.append(MRT[i][3])
                rs_toprom.append(MRT[i][0])               
        
        if len(rs_toprom) != 0:
            r_prom_Vtan.append([np.mean(rs_toprom),np.mean(vt_toprom)])
            r_prom_Vtan_ds.append([np.std(rs_toprom),np.std(vt_toprom)])
        
        matrizRp_Vtp = np.array(r_prom_Vtan)
        matrizRp_Vtp_ds = np.array(r_prom_Vtan_ds)
    return matrizRp_Vtp, matrizRp_Vtp_ds

def homogVr_R(MRT, p):
    heightMRT = len(MRT)
    #MRT es la matriz que contiene en la primer y cuarta columna los datos a promediar
    #p es la cantidad de datos que quiero obtener del promediado, si p->inf se reobtiene los datos originales
    
    #cuantos puntos quiero para la curva?, se define al definir la cantidad de bloques en los que se promedia
    blockBorders = np.linspace(min(MRT[:,0]), max(MRT[:,0]), p)
    #se devuelve esta lista transformada a array
    r_prom_Vr = []
    r_prom_Vr_ds=[]
    
    for j in range(p):
        rs_toprom = []
        vr_toprom = []
        for i in range(heightMRT):
            if j!= p - 1 and blockBorders[j] < MRT[i][0] < blockBorders[j+1]:     
                vr_toprom.append(MRT[i][2])
                rs_toprom.append(MRT[i][0])     
            #print(blockBorders[j])
            if j == p - 1 and blockBorders[j-2] < MRT[i][0] :
                vr_toprom.append(MRT[i][2])
                rs_toprom.append(MRT[i][0])               
        
        if len(rs_toprom) != 0:
            r_prom_Vr.append([np.mean(rs_toprom),np.mean(vr_toprom)])
            r_prom_Vr_ds.append([np.std(rs_toprom)/np.sqrt(len(rs_toprom)),np.std(vr_toprom)/np.sqrt(len(vr_toprom))])
        
        matrizRp_Vrp = np.array(r_prom_Vr)
        matrizRp_Vrp_ds = np.array(r_prom_Vr_ds)
    return matrizRp_Vrp, matrizRp_Vrp_ds

#%%

#promedio0=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\PIVlab promedio agua.txt')
promedio2=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\PIVlab.27.2.1.txt')
archivo='27.2 % - dos vórtices'

matrizXY2 = np.loadtxt(promedio2,delimiter= ',', skiprows= 3, usecols = (0,1,2,3))
matrizXY2[np.isnan(matrizXY2)] = 0

px2,py2=matrizXY2[:,0], matrizXY2[:,1] #posicion de vx y vy
pvx2,pvy2=matrizXY2[:,2], matrizXY2[:,3] #valores de vx y vy

#std=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\PIVlab.std.27.2.1.txt')
#matrizSTD = np.loadtxt(std,delimiter= ',', skiprows= 3, usecols = (0,1,2,3))
#epx0,epy0=matrizSTD[:,0], matrizSTD[:,1] #posicion de vx y vy
#epvx0,epvy0=matrizSTD[:,2], matrizSTD[:,3] #valores de vx y vy

plt.figure()
plt.quiver(px2,py2,pvx2,pvy2)
plt.xlabel('X [pixel]',fontsize=14)
plt.ylabel('Y [pixel]',fontsize=14)
plt.title('{}'.format(archivo))

#%%

pxc=500

px_filt,py_filt = px2[px2<pxc],py2[px2<pxc]
pvx_filt,pvy_filt = pvx2[px2<pxc],pvy2[px2<pxc]

px_filt1,py_filt1 = px2[px2>pxc],py2[px2>pxc]
pvx_filt1,pvy_filt1 = px2[px2>pxc],py2[px2>pxc]

#%%


plt.figure()
plt.quiver(px_filt,py_filt,pvx_filt,pvy_filt)
plt.xlabel('X [pixel]',fontsize=14)
plt.ylabel('Y [pixel]',fontsize=14)
plt.title('Vórtice 1 - {}'.format(archivo))

centro = plt.ginput(1)
x0,y0 = np.array(centro[0])

#%%

plt.figure()
plt.quiver(px_filt1,py_filt1,pvx_filt1,pvy_filt1)
plt.xlabel('X [pixel]',fontsize=14)
plt.ylabel('Y [pixel]',fontsize=14)
plt.title('Vórtice 2 - {}'.format(archivo))

centro1 = plt.ginput(1)
x1,y1 = np.array(centro1[0])

#%%

esct=34.4321 #cambiar esto con la calibracion tiempo vs frame despues, por ahora tengo la automatica de 53 fps
esct_error=34.4321-60/(5.16-2.85)

x1,y1,vx1,vy1=(px_filt-x0)/escp, (py_filt-y0)/escp, pvx_filt*esct/escp, pvy_filt*esct/escp
x2,y2,vx2,vy2=(px_filt1-x0)/escp, (py_filt1-y0)/escp, pvx_filt1*esct/escp, pvy_filt1*esct/escp


plt.figure()
plt.quiver(x1,y1,vx1,vy1)
plt.ylabel('Posicion Y (cm)',fontsize=14)
plt.xlabel('Posicion X (cm)',fontsize=14)
plt.title('Vórtice 1 - {}'.format(archivo))

plt.figure()
plt.quiver(x2,y2,vx2,vy2)
plt.ylabel('Posicion Y (cm)',fontsize=14)
plt.xlabel('Posicion X (cm)',fontsize=14)
plt.title('Vórtice 2 - {}'.format(archivo))


#%%

#armo la matriz centrada y calibrada
Mxy1,Mxy2=np.array([x1,y1,vx1,vy1]),np.array([x2,y2,vx2,vy2])
M_xy1,M_xy2=Mxy1.T,Mxy2.T

#obtengo las dos matrices a ajustar, la OG en polares y la promediada de la OG
M_rt1,M_rt2 = converTo_polares(M_xy1),converTo_polares(M_xy2)

mPromVt1,mPromVt_sd1 = homogVtan_R(M_rt1, 110)
mPromVr1,mPromVr_sd1 = homogVr_R(M_rt1, 110)

mPromVt2,mPromVt_sd2 = homogVtan_R(M_rt2, 110)
mPromVr2,mPromVr_sd2 = homogVr_R(M_rt2, 110)

#%%

#Estan negativos por el sentido de los vectores de v, no afecta en los resultados excepto el signo de Gamma (el vortice)
#para dejarlo bonito lo multiplico por -1

M_rt1[:,3],M_rt2[:,3]=-1*M_rt1[:,3],-1*M_rt2[:,3]
mPromVt1[:,1],mPromVt2[:,1]=-1*mPromVt1[:,1], -1*mPromVt2[:,1]

def func(x,a,b,c):
    return ((a/(x-c))*(1 - np.exp(-((x-c)**2/b**2))))

def vortice(x,a,b,c):
    return (a/(np.pi*b**2))*np.exp(-((x-c)**2/b**2))

#%%

#grafico la OG y la promediada

plt.figure()
plt.plot(M_rt1[:,0], M_rt1[:,3], '.',label='Velocidad tangencial',color='b')
plt.plot(mPromVt1[:,0],mPromVt1[:,1], '.', label='Promedio-velocidad tangencial',color='r')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 15)
plt.title('$V_θ$(r) - vórtice 1 - {}'.format(archivo))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc='upper right')


plt.figure()
plt.plot(M_rt2[:,0], M_rt2[:,3], '.',label='Velocidad tangencial',color='b')
plt.plot(mPromVt2[:,0],mPromVt2[:,1], '.', label='Promedio-velocidad tangencial',color='r')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 15)
plt.title('$V_θ$(r) - vórtice 2 - {}'.format(archivo))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc='upper right')


#%%

#VER DE GRAFICAR VORTICIDAD DIRECTAMENTE

rmin= 0.0 #max o menos a ojo 
rmax= 6.5

#filtro y tomo el rango entre rmin y rmax para el ajuste
subMatrix0 = matrizProm_Rp_Vtp[matrizProm_Rp_Vtp[:,0] >rmin]
subMatrix = subMatrix0[subMatrix0[:,0] <rmax]
subMatrix_sd0=matrizProm_Rp_Vtp_sd[matrizProm_Rp_Vtp[:,0] >rmin]
subMatrix_sd = subMatrix_sd0[subMatrix0[:,0] <rmax]

#ajuste burger al promedio de la dispersion con error en Y=desviacion estandar del promedio
popt_prom, pcov_prom = curve_fit(func, subMatrix[:,0], subMatrix[:,1],sigma= subMatrix_sd[:,1],absolute_sigma=True)
rD_prom = np.linspace(min(subMatrix[:,0]),max(subMatrix[:,0]),len(subMatrix[:,0]))
perr_prom=np.diag(pcov_prom)

plt.figure()
plt.title('$V_θ$(r) promedio - {}'.format(archivo))
#plt.plot(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], '.',label='datos')
plt.errorbar(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], xerr = matrizProm_Rp_Vtp_sd[:,0], yerr =matrizProm_Rp_Vtp_sd[:,1],fmt='.',label='Datos',zorder=1,color='b')
plt.plot(rD_prom,func(rD_prom,*popt_prom), label='Ajuste',zorder=2,color='r')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()

from scipy import stats

fitted_function=func(subMatrix[:,0],*popt_prom)
dof=len(subMatrix[:,1])-len(popt_prom)
subMatrix_sd[subMatrix_sd == 0] = 0.1
chisquare=sum((subMatrix[:,1]-fitted_function)**2/(subMatrix_sd[:,1]**2))
chi2red=chisquare/dof
p=1-stats.chi2.cdf(chisquare,dof)


gamma=popt_prom[0]*2*np.pi
egamma=perr_prom[0]*2*np.pi

print('Chi-square es', chisquare)
print('Chi-square reducido es', chi2red)
print('Gamma es', gamma,'+-', egamma)
print('rc o sqrt(4nu/alpha) es', popt_prom[1],'+-', perr_prom[1])
print('r0 es', popt_prom[2],'+-', perr_prom[2])
#print('CTE es', popt_prom[3],'+-', perr_prom[3])
