# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 17:18:02 2021

@author: aqs
"""

import matplotlib.pyplot as plt
import numpy as np
"""
y1 = np.array(np.loadtxt('freq10_Voltage0.0721831.txt'))[:,1]
x1 = np.array(np.loadtxt('freq10_Voltage0.0721831.txt'))[:,0]

y2 = np.array(np.loadtxt('freq10_Voltage0.17816901.txt'))[:,1]
x2 = np.array(np.loadtxt('freq10_Voltage0.17816901.txt'))[:,0]

y3 = np.array(np.loadtxt('freq10_Voltage0.355.txt'))[:,1]
x3 = np.array(np.loadtxt('freq10_Voltage0.355.txt'))[:,0]

y4 = np.array(np.loadtxt('freq10_Voltage1.txt'))[:,1]
x4 = np.array(np.loadtxt('freq10_Voltage1.txt'))[:,0]

list_max=[]

list_max.append([np.max(y1),np.max(y2),np.max(y3),np.max(y4)])
print(list_max)

list_min=[]

list_min.append([np.min(y1),np.min(y2),np.min(y3),np.min(y4)])
print(list_min)

c=[list_max,list_min]
with open("freq_10_AMax_AMin.txt", 'w') as file:
  for m in zip(*c):
    file.write("{0}\t{1}\n".format(*m))
 """

maximums= np.array(np.loadtxt('frequencies_AMax_AMin_last.txt'))[:,1]
minimums= np.array(np.loadtxt('frequencies_AMax_AMin_last.txt'))[:,2]

y1=maximums[0:4]
x1=minimums[0:4]
y2=maximums[4:8]
x2=minimums[4:8]
y3=maximums[8:12]
x3=minimums[8:12]
y4=maximums[12:16]
x4=minimums[12:16]


plt.plot(x1,'*',y1,'*',label='100 Hz',color='blue')
plt.plot(x2,'o',y2,'o',label='50 Hz',color='green')
plt.plot(x3,'D',y3,'D',label='10 Hz',color='red')
plt.plot(x4,'H',y4,'H',label='2 Hz',color='black')
plt.legend()
plt.title('Max(B) and Min(B)')
plt.xlabel('points')
plt.ylabel('B(mT)')
#plt.plot(x_200,minimums,'.',)
plt.savefig('ww.pdf')
plt.show()

c2=[maximums,minimums]
with open("frequencies_AMax_AMin.txt", 'w') as file:
  for m in zip(*c2):
    file.write("{0}\t{1}\n".format(*m))
"""
plt.style.use('seaborn-notebook')
fig, axs = plt.subplots(2, 2)
#fig.suptitle('Frequency =100,50,10,2')
axs[0, 0].plot(x1,'o')
axs[0, 0].plot(y1,'*')
axs[0,0].set_ylim(-1, 30)
#axs[0,0].set_xlim(0, 0.03)
axs[0, 0].set_title('Frequency=100Hz')
axs[0, 1].plot(x2,'o')
axs[0, 1].plot(y2,'*')
axs[0,1].set_ylim(-1, 30)
axs[0, 1].set_title('Frequency=50Hz')
#axs[0,1].set_xlim(0, 0.03)
#axs[0, 1].set_title('Axis [0, 1]')
axs[1, 0].set_title('Frequency=10Hz')
axs[1, 0].plot(x3,'o')
axs[1, 0].plot(y3,'*')
axs[1,0].set_ylim(-1, 30)

#axs[1,0].set_xlim(0, 0.03)
#axs[1, 0].set_title('Axis [1, 0]')
axs[1, 1].plot(x4,'o')
axs[1, 1].plot(y4,'*')
#axs[1,1].set_xlim(0, 0.03)
axs[1, 1].set_title('Frequency=2Hz')
axs[1,1].set_ylim(-1, 30)

#axs[1, 1].set_title('Axis [1, 1]')

for ax in axs.flat:
    ax.set(ylabel='B(mT)',xlabel='points')

# Hide x labels and tick labels for top plots and y ticks for right plots.
#for ax in axs.flat:
#    ax.label_outer()

fig.savefig('max_min.pdf')

"""
