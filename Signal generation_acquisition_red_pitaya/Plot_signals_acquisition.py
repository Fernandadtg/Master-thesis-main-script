# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 14:31:25 2021

@author: aqs
"""

import matplotlib.pyplot as plt
import numpy as np

y1 = np.array(np.loadtxt('freq100_Voltage0.0721831.txt'))[:,1]
x1 = np.array(np.loadtxt('freq100_Voltage0.0721831.txt'))[:,0]

y2 = np.array(np.loadtxt('freq100_Voltage0.17816901.txt'))[:,1]
x2 = np.array(np.loadtxt('freq100_Voltage0.17816901.txt'))[:,0]

y3 = np.array(np.loadtxt('freq100_Voltage0.355.txt'))[:,1]
x3 = np.array(np.loadtxt('freq100_Voltage0.355.txt'))[:,0]

y4 = np.array(np.loadtxt('freq100_Voltage1.txt'))[:,1]
x4 = np.array(np.loadtxt('freq100_Voltage1.txt'))[:,0]

list_max=[]

list_max.append([np.max(y1),np.max(y2),np.max(y3),np.max(y4)])
print(list_max)

list_min=[]

list_min.append([np.min(y1),np.min(y2),np.min(y3),np.min(y4)])
print(list_min)

c=[list_max,list_min]
with open("freq_2_AMax_AMin.txt", 'w') as file:
  for m in zip(*c):
    file.write("{0}\t{1}\n".format(*m))
plt.style.use('seaborn-notebook')
fig, axs = plt.subplots(2, 2)
fig.suptitle('Frequency =100 Hz \n B=(2,5,10,28.292)')
axs[0, 0].plot(x1, y1)
axs[0,0].set_xlim(0, 0.03)
#axs[0, 0].set_title('Axis [0, 0]')
axs[0, 1].plot(x2, y2, 'tab:orange')
axs[0,1].set_xlim(0, 0.03)
#axs[0, 1].set_title('Axis [0, 1]')
axs[1, 0].plot(x3, y3, 'tab:green')
axs[1,0].set_xlim(0, 0.03)
#axs[1, 0].set_title('Axis [1, 0]')
axs[1, 1].plot(x4, y4, 'tab:red')
axs[1,1].set_xlim(0, 0.03)
axs[1,1].set_ylim(0, 30)
#axs[1, 1].set_title('Axis [1, 1]')

for ax in axs.flat:
    ax.set(xlabel='time(s)', ylabel='B(mT)')

# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
#    ax.label_outer()

fig.savefig('Frequency =100 Hz  B=(2,5,10,28.292).pdf')