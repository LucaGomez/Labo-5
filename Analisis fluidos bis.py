# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 19:55:39 2022

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
plt.rcParams.update({'font.size': 11})

#%%

#Calibracion escala espacial px vs cm

imagen_escala = plt.imread(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\escala perladoasd2.jpg')
plt.figure()
plt.imshow(imagen_escala)
plt.xlabel('Posicion X (px)')
plt.ylabel('Posicion Y (px)')

ptos= plt.ginput(2)
ptos=np.array(ptos)

#la roto para que quede mas o menos bien 

p_x=ptos[:,0]
p_y=ptos[:,1]
deltax=abs(p_x[1]-p_x[0])
deltay=abs(p_y[1]-p_y[0])
ang=abs(90-abs(180*np.arctan(deltay/deltax)/np.pi))

esc=ndimage.rotate(imagen_escala, ang)

#grafico la figura para tomar los puntos
mesc=esc[200:600,450:720]
plt.figure()
plt.imshow(mesc)

#%%

puntos_esc = plt.ginput(8)
puntos_esc = np.array(puntos_esc)

coord_esc=puntos_esc[:,1]
diff_esc= np.diff(coord_esc)

# Promedio de los valores
escala=np.mean(diff_esc)
escala_error=np.std(diff_esc)/np.sqrt(len(diff_esc)) 

print('la escala de px/5mm es ',escala ,'+-',escala_error)

#%%

escp=31.118324986915233/5*10 #px por cm
escp_error=0.39669604731659686*10


#%%

#Promedio centro 100 frames

#archivo agua
promedio_agua=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\PIVlab promedio agua.txt')
matrizXYa = np.loadtxt(promedio_agua,delimiter= ',', skiprows= 3, usecols = (0,1,2,3))

#archivo 5050
promedio_centro0=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\PIVlab promedio lejos y centro.txt')
matrizXY0 = np.loadtxt(promedio_centro0,delimiter= ',', skiprows= 3, usecols = (0,1,2,3))

promedio_centro1=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\PIVlab promedio lejos y centro 2.txt')
matrizXY1 = np.loadtxt(promedio_centro1,delimiter= ',', skiprows= 3, usecols = (0,1,2,3))

#matrizXY = np.vstack((matrizXY0, matrizXY1))

pxa,pya=matrizXYa[:,0], matrizXYa[:,1] #posicion de vx y vy
pvxa,pvya=matrizXYa[:,2], matrizXYa[:,3] #valores de vx y vy

px0,py0=matrizXY0[:,0], matrizXY0[:,1] #posicion de vx y vy
pvx0,pvy0=matrizXY0[:,2], matrizXY0[:,3] #valores de vx y vy

px1,py1=matrizXY1[:,0], matrizXY1[:,1] #posicion de vx y vy
pvx1,pvy1=matrizXY1[:,2], matrizXY1[:,3] #valores de vx y vy

#plt.figure()
#plt.quiver(px,py,pvx,pvy)

#%%

#centro = plt.ginput(1)
#x0,y0 = np.array(centro[0])

xa,ya=(443.37130376344095, 329.1645915007058)
x0,y0=(439.07526881720446, 312.81527363964256)
x1,y1=(447.91666666666686, 323.1649678532226)

#%%

esct0=53 #cambiar esto con la calibracion tiempo vs frame despues, por ahora tengo la automatica de 53 fps
esct1=20

x0,y0=(px0-x0)/escp, (py0-y0)/escp
vx0,vy0=pvx0*esct0/escp, pvy0*esct0/escp

x1,y1=(px1-x1)/escp, (py1-y1)/escp
vx1,vy1=pvx1*esct1/escp, pvy1*esct1/escp

#junto los dos archivos de 5050 bien centradas y calibradas separadamente  (tienen distinto fps)

#x,y,vx,vy=np.concatenate((x0, x1)), np.concatenate((y0, y1)), np.concatenate((vx0, vx1)), np.concatenate((vy0, vy1))

x,y,vx,vy=(pxa-xa)/escp, (pya-ya)/escp, pvxa*esct0/escp, pvya*esct0/escp

#%%

#veo en la escala correcta

def vn(vx,vy):
 return (vx**2+vy**2)**(1/2)

vm=vn(vx,vy) #modulo de las velocidades

plt.figure()
plt.quiver(x,y,vx,vy,vm)
plt.ylabel('Posicion Y (cm)')
plt.xlabel('Posicion X (cm)')
plt.title('50:50')
plt.show()

#%%

#Funciones corro esto una vez y listo

#convertir en polares

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
            r_prom_Vtan_ds.append([np.std(rs_toprom)/np.sqrt(len(rs_toprom)),np.std(vt_toprom)/np.sqrt(len(vt_toprom))])
        
        matrizRp_Vtp = np.array(r_prom_Vtan)
        matrizRp_Vtp_ds = np.array(r_prom_Vtan_ds)
    return matrizRp_Vtp, matrizRp_Vtp_ds


#burgerking
def func(x,a,b,c,d):
    return ((a/(x-d))*(1 - np.exp(-((x-d)**2/b**2))))+c

#%%

#armo la matriz centrada y calibrada
Mxy=np.array([x,y,vx,vy])
M_xy=Mxy.T

#obtengo las dos matrices a ajustar, la OG en polares y la promediada de la OG
M_rt = converTo_polares(M_xy)
matrizProm_Rp_Vtp,matrizProm_Rp_Vtp_sd = homogVtan_R(M_rt, 100)

#%%

#Estan negativos por el sentido de los vectores de v, no afecta en los resultados excepto el signo de Gamma (el vortice)
#para dejarlo bonito lo multiplico por -1

M_rt[:,3]=-1*M_rt[:,3]
matrizProm_Rp_Vtp[:,1]=-1*matrizProm_Rp_Vtp[:,1]

#%%

#grafico la OG y la promediada

plt.figure()
plt.plot(M_rt[:,0], M_rt[:,3], '.',label='En polar y centrado')
plt.plot(matrizProm_Rp_Vtp[:,0],matrizProm_Rp_Vtp[:,1], '.', label='centrado promediado angular')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 15)
plt.title('agua')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()

#%%

#esto es porque habia tirado error en un momento con los ceros, pero generalmente no importa
inf, nan = np.inf, np.nan 
a_nan = np.where(np.isnan(M_rt[:,3])) 
M_rt[a_nan,3]=0

#Grafico y ajuste burger de Velocidad Tangencial en funcion del radio
#initialGuess = [0,(0.002)*100,(0.02)*100,(0.01)*100] 
popt, pcov = curve_fit(func, M_rt[:,0], M_rt[:,3])
rD = np.linspace(min(M_rt[:,0]),max(M_rt[:,0]),len(M_rt[:,0]))

plt.figure()
plt.plot(M_rt[:,0], M_rt[:,3], '.',label='centrado')
plt.plot(rD,func(rD,*popt), label='ajuste')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 15)
plt.title('50:50')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()

#%%

#Grafico de Vel. Tangencial en funcion del radio PERO PARA EL PROMEDIADO ESPACIAL
#fijo rmin para deshacerme del hueco alrededor del vortice (si es que hay, sino rmin=0)
#fijo rmax para deshacerme de efecto de borde

rmin= 0.5 #max o menos a ojo 
rmax= 8.9

#filtro y tomo el rango entre rmin y rmax para el ajuste
subMatrix0 = matrizProm_Rp_Vtp[matrizProm_Rp_Vtp[:,0] >rmin]
subMatrix = subMatrix0[subMatrix0[:,0] <rmax]
subMatrix_sd0=matrizProm_Rp_Vtp_sd[:,0][matrizProm_Rp_Vtp[:,0] >rmin]
subMatrix_sd = subMatrix_sd0[subMatrix0[:,0] <rmax]

#lo mismo de recien con los ceros
b_nan = np.where(np.isnan(subMatrix[:,1])) 
subMatrix[b_nan,1]=0

#ajuste burger al promedio de la dispersion con error en Y=desviacion estandar del promedio
popt_prom, pcov_prom = curve_fit(func, subMatrix[:,0], subMatrix[:,1],sigma= subMatrix_sd,absolute_sigma=True)
rD_prom = np.linspace(min(subMatrix[:,0]),max(subMatrix[:,0]),len(subMatrix[:,0]))
perr_prom=np.diag(pcov_prom)

plt.figure()
plt.title('V_tangencial(r) agua')
#plt.plot(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], '.',label='datos')
plt.errorbar(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], xerr = matrizProm_Rp_Vtp_sd[:,0], yerr =matrizProm_Rp_Vtp_sd[:,1],fmt='.',label='datos')
plt.plot(rD_prom,func(rD_prom,*popt_prom), label='ajuste')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()

modelPredictions = func(subMatrix[:,0], *popt) 

absError = modelPredictions -  subMatrix[:,1]

SE = np.square(absError) # squared errors
MSE = np.mean(SE) # mean squared errors
RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (np.var(absError) / np.var(subMatrix[:,1]))
print('RMSE:', RMSE)
print('R-squared:', Rsquared)

from scipy import stats
fitted_function=func(matrizProm_Rp_Vtp[:,0],*popt_prom)
dof=len(matrizProm_Rp_Vtp[:,0])-len(popt_prom)
chi2stat=stats.chisquare(matrizProm_Rp_Vtp[:,1],fitted_function,dof)
p=1-stats.chi2.cdf(chi2stat,dof)

print('Chi-square es', chi2stat)

print('Gamma/2pi es', popt_prom[0],'+-', perr_prom[0])
print('sqrt(4nu/alpha) es', popt_prom[1],'+-', perr_prom[1])
print('CTE es', popt_prom[2],'+-', perr_prom[2])
print('r0 es', popt_prom[3],'+-', perr_prom[3])

#%%

#probe quitando el r0 y no cambia mucho

def func1(x,a,b,c):
    return ((a/x)*(1 - np.exp(-(x**2/b**2))))+c

popt_prom1, pcov_prom1 = curve_fit(func1, subMatrix[:,0], subMatrix[:,1])

plt.figure()
plt.title('V_tangencial(r)')
plt.plot(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], '.',label='datos')
plt.plot(rD_prom,func1(rD_prom,*popt_prom1), label='ajuste')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()

#%%

#Veo el vortice

def vortice(x,a,b,c):
    return 2*(a/b)*np.exp(-((x-c)**2/b**2))

w=vortice(matrizProm_Rp_Vtp[:,0],popt_prom[0],popt_prom[1],popt_prom[3])


plt.figure()
plt.title('Vortice en z')
plt.plot(matrizProm_Rp_Vtp[:,0], w, '.',label='datos')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Vorticidad [1/s]", fontsize = 15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
