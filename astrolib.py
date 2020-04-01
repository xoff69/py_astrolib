import sys 
import numpy as np
import matplotlib.pyplot as plt

h=6.62e-34 #m2 kg / s
c=3e8
k=1.38e-23
G= 6.67e-11 #si
MSUN=2e30 #kg
LSUN=3.846e26 #w
pc=3.0856e16 #m

def ua2pc(ua):
    return 4.85e-6*ua;

def pc2ua(pc):
    return pc/4.85e-6;

def pc2m(var_pc):
    return var_pc*pc;

def rad2arcsec(r):
    return  206264.80625*r


# corps noir en fonction de la temperature et de la longueur d'onde l
def B_l(l,T):
    bel=1/(-1+np.exp((h*c)/(l*k*T)))
    return bel * (2 * h* c**2)/(l**5)
def B_v(v,T):
    bel=2*h*v**3
    bas=(c**2)*(np.exp(h*v/(k*T))-1)
    return bel/bas

def arcsec2rad(arcs):
    return (arcs / ((180.0 / np.pi) * 3600 ))

def arcmin2rad(arcs):
    return (arcs / ((180.0 / np.pi) * 60 ))

#conversion arcversparsec
def arcsec2parsec(arcsec,distance):
    return arcsec2rad(arcsec)*distance

# convertit des m2 en parsec2
def m2topc2(m2):
    return m2*1.050264757570000e-33

#donne un ordre d'idee de la valeur en puissance de 10
def quick_ordre_id(num):
    increment=0
    while num>10:
        num=num/10
        increment=increment+1
    return str(num)+" 10^"+str(increment)

# determination du z
def zcalc(lemis,lrepos):
    return (lemis-lrepos)/lrepos

# vitesse de la galaxie en m/s en fonction de z
def vitessegalaxie(z):
    return c*((z+1)**2-1)/((z+1)**2+1)

# donne la distance en AL 
# le H0 doit est etre mis en m
def distanceAL(H0,z):
    return vitessegalaxie(z)*3.26e6/H0

def lambdadebroglie(m,v):
    return h/(m*v)

zc=zcalc(6970,6562)
print("z=",zc)
v=vitessegalaxie(zc)
print(quick_ordre_id(v))
print((distanceAL(70*1000,zc)))
print(quick_ordre_id(distanceAL(70*1000,zc)))
print(B_v(5,5))