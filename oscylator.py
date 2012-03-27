#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''
Created on 13-03-2012 13:19:22
@author: mc
'''

import csv
import os
from math import sqrt

#t = []       # czas w s
#x = []       # położenie w m
#v = []       # prędkość w m*s**-1

def txv_poczatkowe():
    t = [0.0]
    x = [1.0]
    v = [0.0]
    return t,x,v

k = 1.0         # wsp. sprężystości
m = 1.0         # masa w kg
t_max = 1e1     # koniec symulacji

def rownomierny_wykladnik(i,w1,w2):
    # dla i z zakresu 0..1 zwraca liczbę z zakresu 10**w1..10**w2
    return 10**(w1+(w2-w1)*i)

def sila(x):
    return -k*x

def energia(x,v):
    return m*v**2+k*x**2

def calkuj_euler(t,x,v,dt):
    a = sila(x[-1])/m
    v.append( v[-1]+a*dt )
    x.append( x[-1]+v[-2]*dt )
    t.append( t[-1]+dt )
    return 0,0,0

def calkuj_verlet_klasyczny(t,x,v,dt):
    if len(x)<2: x.append(x[-1])
    a = sila(x[-1])/m
    x.append( 2*x[-1]-x[-2]+a*dt**2 )
    v.append( (x[-1]-x[-3])/(2*dt) )
    t.append( t[-1]+dt )
    return 0,-1,0   # przesunięcie jest liczone z wyprzedzeniem

def calkuj_verlet_predkoscowy(t,x,v,dt):
    if len(x)<2: x.append(x[-1])
    a = sila(x[-1])/m
    x.append( x[-1]+dt*v[-1]+0.5*(dt**2)*a )
    v.append( v[-1]+0.5*dt*(a+sila(x[-1])/m) )
    t.append( t[-1]+dt )
    return 0,0,0

def calkuj_rk_2(t,x,v,dt):
    v1 = v[-1]
    xp = x[-1]+v1*dt
    
    a1 = sila(x[-1])/m
    a2 = sila(xp)/m
    
    v.append( v[-1]+dt*(a1+a2)/2.0 )
    x.append( x[-1]+dt*(v[-2]+v[-1])/2.0 )
    t.append( t[-1]+dt )
    
    return 0,0,0



def testuj_algorytm(algorytm, nazwa_pliku_danych):
    t,x,v = txv_poczatkowe()
    dt = 0.05
    H0 = energia(x[-1],v[-1])
    
    skryba = csv.writer(open('dane/'+nazwa_pliku_danych+".csv", 'w',
                             newline=''), delimiter=' ')
    
    SdH2 = 0.0
    l_dH2 = 0
    
    while t[-1] <= t_max:
        t_slip, x_slip, v_slip = algorytm(t,x,v,dt)
        SdH2 += (H0-energia(x[-1+x_slip],v[-1+v_slip]))**2
        l_dH2 += 1
        skryba.writerow([t[-1+t_slip], x[-1+x_slip], v[-1+v_slip],
                         sqrt(SdH2/l_dH2) ])
    
    skryba = csv.writer(open('dane/'+nazwa_pliku_danych+"_krok.csv", 'w',
                             newline=''), delimiter=' ')
    
    i_max = 20
    for i in range(i_max+1):
        dt = rownomierny_wykladnik(float(i)/i_max, -4, 0)
        t,x,v = txv_poczatkowe()
        E0 = energia(x[-1],v[-1])
        
        num_of_points = 0
        sum_of_en_err_sq = 0
        
        while t[-1] <= t_max:
            algorytm(t,x,v,dt)
            num_of_points += 1
            sum_of_en_err_sq += (energia(x[-1+x_slip],v[-1+v_slip])-E0)**2
            
        skryba.writerow([ dt, sqrt(sum_of_en_err_sq)/E0 ])
    
    skryba = None
    
    os.system("cd wykresy && bash rysuj_razem "+nazwa_pliku_danych)


#testuj_algorytm(calkuj_euler, "Euler")
#testuj_algorytm(calkuj_verlet_klasyczny, "Verlet_klas")
#testuj_algorytm(calkuj_verlet_predkoscowy, "Verlet_pred")
#testuj_algorytm(calkuj_rk_2, "RK_2")



