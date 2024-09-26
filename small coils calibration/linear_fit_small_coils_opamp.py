#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 08:52:35 2021

@author: u
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from scipy.stats import linregress
 
#B_mt_teslameter_1009 .txt
#y=np.asarray(np.loadtxt('/home/u/Downloads/data_2409-small_coils_opamp_cld.txt'))[:,1]*1000
#x=np.asarray(np.loadtxt('/home/u/Downloads/data_2409-small_coils_opamp_cld.txt'))[:,0]

y=np.asarray(np.loadtxt('B_mt_teslameter_1009 .txt'))*1000
x=np.asarray(np.linspace(0,1,21))
#x=np.asarray(np.loadtxt('/home/u/Downloads/data_2409-small_coils_opamp_cld.txt'))[:,0]
#print(linregress(x[0:13],y[0:13]))
#slope, intercept, r_value, p_value, std_err = linregress(x[0:13], y[0:13])
print(linregress(x,y))
slope, intercept, r_value, p_value, std_err = linregress(x, y)
plt.plot(x,y,'ro',label='data')
print(x)
print(y)
plt.plot(x, intercept + slope*x, 'black', label=f'$B = {slope:.3f}*I {intercept:+.3f}$')
#plt.plot(x, intercept + slope*x, 'black', label='fitted line: slope={:.3f}'.format(slope))
plt.legend()
plt.xlabel('I(A)')
plt.ylabel('B(mT)') #max voltage 20

plt.savefig('B_vs_I_smallcoils.pdf')
plt.show()