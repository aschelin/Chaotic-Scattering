# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 17:25:09 2021

@author: asche
"""
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, 100, 1000)

E=.01
w=0
dx = 5
dy=0
r = dx

x0 = 0
y0 = 1

U1=-np.exp(-((x0+dx)**2+(y0+dy)**2));
U2=-np.exp(-((x0-dx)**2+(y0+dy)**2));
U = U1+U2

v = np.sqrt(np.abs(2*(E-U)))   

theta = np.arctan2(y0,x0);

vx0 = v*np.cos(theta);
vy0 = v*np.sin(theta);

vx0 = v
vy0=0

y0 = [vx0, x0, vy0, y0] 

# Equations of Motion 
def gauss2pot(y, t, w,r):
    y1, y2, y3, y4 = y

    dpx = r*(np.cos(2*np.pi*w*t))
    dpy = r*(np.sin(2*np.pi*w*t)); 
    
    Fy1 = -2*(y4+dpy)*np.exp(-((y2+dpx)**2+(y4+dpy)**2))   
    Fx1 = -2*(y2+dpx)*np.exp(-((y2+dpx)**2+(y4+dpy)**2))
    Fy2 = -2*(y4-dpy)*np.exp(-((y2-dpx)**2+(y4-dpy)**2)) 
    Fx2 = -2*(y2-dpx)*np.exp(-((y2-dpx)**2+(y4-dpy)**2))

    dp1 = Fx1+Fx2
    dp2 = Fy1+Fy2
    
    dydt = [dp1,y1,dp2,y3]
    return dydt



sol = odeint(gauss2pot, y0, t,args=(w,r))

plt.plot(sol[:,1],sol[:,3], 'b', label='theta(t)')
