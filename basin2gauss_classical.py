# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 17:25:09 2021
This program calculates the escape basin for the elastic pendulum
@author: aschelin
"""
#importing libraries 
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import imagesc as imagesc

# define integration time
tf = 1000
t = np.linspace(0, tf, 1000)

# Parameters used: E (energy), r (distance between each potential D=2r), w (type angular velocity)
dx=5
dy=0
r= dx
E = 0.01
w=0.00
# Number of pixels in the basin
nx = 400
ny = 400

# declaring size of basin variables 
yvxfim = np.zeros((nx,ny))
yvyfim= np.zeros((nx,ny))
yyfim= np.zeros((nx,ny))
yxfim= np.zeros((nx,ny))


# phase space limits 
p100 = np.linspace(-3,3,nx)
q100 = np.linspace(-3,3,nx)

# star loop to calculate each pixel
for ix in range(nx): 
    print(ix)
    for iy in range(ny):
        x0 = q100[ix]
        y0 = p100[iy]
        
        
        U1=-np.exp(-((x0+dx)**2+(y0+dy)**2));
        U2=-np.exp(-((x0-dx)**2+(y0+dy)**2));
        U = U1+U2

        v = np.sqrt(np.abs(2*(E-U)))   

        theta = np.arctan2(y0,x0);

        vx0 = v*np.cos(theta);
        vy0 = v*np.sin(theta);

#        vx0 = v
#        vy0=0 

        z0 = [vx0, x0, vy0, y0] 

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
# Integrating the ODE:
        sol = odeint(gauss2pot, z0, t,args=(w,r))
        yvxfim[ix,iy] = sol[-1,0]
        yvyfim[ix,iy] = sol[-1,2]
        yxfim[ix,iy] = sol[-1,1]
        yyfim[ix,iy] = sol[-1,3]
                    

#imagesc.plot(yvyfim)
imagesc.plot(yvxfim)