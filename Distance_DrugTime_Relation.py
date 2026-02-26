#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 12:23:50 2026

@author: huima
"""
import numpy as np
import random
import math
import matplotlib.pyplot as plt
from scipy import interpolate
plt.close('all')


x = np.loadtxt('x_N2.txt')
cdf = np.loadtxt('cdf_N2.txt')
def cat2Angle(x, cdf):
    u = np.random.rand()
    samples = np.interp(u, cdf, x)
    return samples  


def one_worm():
    width = 60
    #################
    #sedimentation velocity = 0.52 mm/s
    vg = 0.52
    vswim = 0.3

    #Stokes radius r = 0.37 mm
    #Swimming velocity = 0.3 mm/s

    #time intervals for direction change
    mu = 4.65#s
    sigma = 1.95#s
    N = 500
    X = np.zeros(N)
    Y = np.zeros(N)
    T = np.zeros(N)
    X[0] = width/2

    WallPause = 29.7 # walltime = [51, 15, 24, 65, 33, 20, 17, 27, 15]

    for i in range(1,N):
        dt = random.gauss(mu, sigma)
        if X[i-1] < 1:
            X[i], Y[i] = 1, Y[i-1]
            T[i] = T[i-1] + WallPause
        elif X[i-1] > width - 1:
            X[i], Y[i] = width - 1, Y[i-1]
            T[i] = T[i-1] + WallPause
        else:
            #direction = random.random()*math.pi*2
            #vx = vswim*math.sin(direction); vy = vswim*math.cos(direction) + vg
            direction = cat2Angle(x, cdf)/180*math.pi
            LR = np.random.rand()-0.5
            vx = abs(LR)/LR*vswim*abs(math.sin(direction))
            vy = vswim*abs(math.cos(direction)) + vg
            X[i] = X[i-1] + vx*dt;
            Y[i] = Y[i-1] + vy*dt;
            T[i] = T[i-1] + dt
            if vy < 0:
                print(vy)
        if Y[i] > 200:
            Y[i:] = 200
            T[i:] = T[i]
            break

    #plt.plot(Y, T)
    f = interpolate.interp1d(Y, T, kind='linear')
    xnew = np.arange(0, 200, 1)
    ynew = f(xnew)
    #plt.plot(xnew,ynew)
    return xnew, ynew

WormNum = 5000
MeanTime = np.zeros(200)
ErrorBar = np.zeros(200)
lower, upper = np.zeros(200), np.zeros(200)
Y = np.zeros((WormNum,200))
for i in range(WormNum):
    Cy,Y[i,:] = one_worm()
    
for i in range(200):
    MeanTime[i] = np.mean(Y[:,i])
    ErrorBar[i] = np.std(Y[:,i])
    lower[i] = np.percentile(Y[:,i], 2.5, axis=0)
    upper[i] = np.percentile(Y[:,i], 97.5, axis=0)

plt.errorbar(Cy, MeanTime, yerr=ErrorBar, fmt='o')
plt.xlabel('Distance beneath water (mm)', fontsize = 16)
plt.ylabel('Time exposed to drug (s)', fontsize = 16)

plt.figure()
plt.plot(Cy, MeanTime,label = 'Simulated mean time')
plt.fill_between(Cy, lower, upper, alpha=0.3, label = '95% confidence')
Assm = Cy/0.52
plt.plot(Cy, Assm, label = 'Experiment')
plt.xlabel('Distance beneath water (mm)', fontsize = 16)
plt.ylabel('Time exposed to drug (s)', fontsize = 16)
plt.legend(fontsize=12)
plt.show()


plt.figure()
TimeWin = upper - lower
RLB = TimeWin/MeanTime
plt.plot(Cy, RLB*100)
plt.xlabel('Distance beneath water (mm)', fontsize = 16)
plt.ylabel('Error level (%)', fontsize = 16)











