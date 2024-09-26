#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 12:04:53 2022

@author: u
"""


from __future__ import print_function
import sys
sys.path.insert(0,'/home/u/pyNMR')
import math
import scipy as sp
import re
from operator import itemgetter
import numpy as np
from numpy import mean
from scipy import fft
from scipy.fftpack import fftshift
import nmrDataMod as ndm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import funct_ig as funct_ig
import glob
from glob import iglob
from itertools import groupby
from statistics import mean 
from scipy.stats import linregress
import os
from scipy import optimize
import pylab as plt
import natsort
import pickle5 as pickle
import collections, functools, operator
from numpy import trapz
from scipy.integrate import simps
from collections import Counter
#import scipy.integrate.trapezoid as trapezoid

def flatten(t):
    return [item for sublist in t for item in sublist]


#path=  '/home/u/pyNMR/id277/id277_13c_01'        #'/home/u/pyNMR/id280_13c_01' #'#/home/u/pyNMR/07_12/id273_13c_8

    
    
path= '/home/u/pyNMR/nmr_bruker_400/marzo25/id273_13c_19'        #'/home/u/pyNMR/id280_13c_01' #'#/home/u/pyNMR/07_12/id273_13c_8
info_path = glob.glob(os.path.join(path + str('/202203*')))[0] #info_path
print(info_path)

file_title = path.split('/home/u/pyNMR/nmr_bruker_400/marzo25/', maxsplit=1)[1]
print("file_title")
print(file_title)
script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, f'{file_title} Output /')
print("results_dir")
print(results_dir)
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)

off_field=  info_path + str('/param_dyn.p') #off_field
print(off_field)
my_dict = pickle.load(open(off_field,'rb'))

new_arr = [my_dict[key] for key in my_dict]

new_arr = np.array(new_arr)
mag_field_sweep_values =list(flatten(new_arr))


b=[]
b = mag_field_sweep_values


chunk_size = [b.count(i) for i in b]

print("chunk_size[0]", chunk_size[0])

num_group = chunk_size[0]

def chunkify(lst,n):
  return [lst[i::n] for i in range(n)]



files_i_care_about=[]
filenameslist=[]


filenames_list=[]

all_filenames = [fn for fn in glob.iglob(path+ str('*/'))]
ordered2 = natsort.natsorted(all_filenames)
for filenamed in ordered2:
  if os.path.isfile(filenamed+ "acqus"):
     filenames_list.append(filenamed)
   

filtered = [f for f in glob.glob(path + '/*/') if not re.match(info_path, f) ]


ordered = natsort.natsorted(filtered)[1:]

phase_val_list =[]

# take second element for sort
def takeSecond(elem):
    return elem[0]


class Preprocessing_subfolders:
    def __init__(self, path_to_magnetic_field_list,name_dic_file=[],mag_sweep_list_ent=[]):
        self.name_dic_file = pickle.load(open(path_to_magnetic_field_list,'rb'))
        self.mag_sweep_list_ent= [self.name_dic_file[key] for key in self.name_dic_file]
        self.mag_sweep_list_flat = list(flatten(self.mag_sweep_list_ent))


sweepWidthTD2=100000
sizeTD2 = 512  #int(len(data) / self.sizeTD1 / 2)

#0,5.11 x10""-3}} step 0.01e-03 , widthtd2 100000
dwellTime = 1./sweepWidthTD2
#dwellTime = 1./sweepWidthTD2
# self.fidTime = np.linspace(0, (self.sizeTD2-1)*dwellTime, num = self.sizeTD2)
fidTime = np.linspace(0, (sizeTD2-1)*dwellTime, num = sizeTD2)
print("FID time")
print(fidTime)



class pynmr_set_fourier_transform:   
    def __init__(self,filename, fid_sum,data_or=0 ,data_freq=[],ft_sum=0, ft_list=[], fid_inicial= 0, phase_values=0, filenameslist=[], integration_element=0): #ft_sum=0, ft_list=[], phase_values=0, filenameslist=[], integration_element=[]
      self.fid_sum=0
      self.ft_sum=0
      self.ft_list=[]
      self.fid_inicial = 0
      self.integration_element=[]
      self.data = ndm.nmrData(filename, "TopSpin4")
      self.data.leftShift(0,1,0)      #20
      self.data.lineBroadening(1,2,300) #100
      #self.fid_sig=self.data.allFid[2][0]
      self.data_or=self.data.allFid[2][0]
      
      self.data.fourierTransform(2,3)
      self.data.baselineCorrection(3, 4, [[-6000, 0]], scale = "Hz", order = 1, applyLocally = True)
      #self.data.baselineCorrection(3,4,[[-9000,0]], scale = "Hz", order = 1, applyLocally = True)
      self.data_freq = self.data.frequency
      
      #elf.phase_values=self.data.autoPhase0(3,-1, -9000, -3000, scale ='Hz')
     # self.data.phase(3,4,319)#self.phase_values
      self.ft_sum =self.data.allFid[1][0]
      
      #self.integration_element=self.data.integrateAllRealPart(4, -10000, -5000, scale="Hz", part = "real")
      
     

list_ft_sum_inicial = [pynmr_set_fourier_transform(filename,data_freq=[],data_or=0, fid_sum=0,ft_sum=0,ft_list=[],fid_inicial= 0,filenameslist=[], integration_element=0).ft_sum for filename in ordered]
data_freq_p = [pynmr_set_fourier_transform(filename,data_freq=[],data_or=0, fid_sum=0,ft_sum=0,ft_list=[],fid_inicial= 0,filenameslist=[], integration_element=0).data_freq for filename in ordered]
data_o=[pynmr_set_fourier_transform(filename,data_freq=[],data_or=0, fid_sum=0,ft_sum=0,ft_list=[],fid_inicial= 0,filenameslist=[], integration_element=0).data_or for filename in ordered]
#stuff=[filename for filename in ordered]
#print("sttttttttttttttt")
#print(stuff)



list3= list(map(lambda x,y: [x,y], b,data_o))
list4 = list(sorted(list3, key=lambda x:x[0]))
print("c2")
#print(list4[1][1])

list5=[row[1] for row in list4]
list6=[row[0] for row in list4]

#list_ft = list(zip(*chunkify(list_ft_sum_inicial, 32)))
list_data_o = list(zip(*chunkify(list5, chunk_size[0])))
list_data_r = list(zip(*chunkify(list6, chunk_size[0])))


total_data_o= [sum(x) for x in list_data_o]

print("total")
print(len(list_data_o))

LB=300
leftshift=20

data_o_lb=[sp.multiply(total_data_o[k][leftshift:],sp.exp(-fidTime[:len(total_data_o[k][leftshift:])]*LB*np.pi)) for k in range(len(total_data_o))] 
data_freq_p[0]=data_freq_p[0][leftshift:]



size_freq=data_freq_p[0][1]-data_freq_p[0][0]

plt.plot(total_data_o)
plt.show()
#list_ft_sum2=list_ft[::1]
i1= 200#200#50#150#200#209#int(8000/size_freq) #190
i2= 218#220#250#250#250#222#int(5000/size_freq) #260
no_repeat_B =sorted([*set(list6)])
print(no_repeat_B)
for el_2 in data_o_lb:
   plt.plot(data_freq_p[0],el_2)
#plt.title(str(file_title)+str("NMR signal or"), loc='center')
plt.xlabel('frequency')
plt.ylabel('NMR signal (arb. units)')
#plt.xlim([190,230])
#plt.xlim([-10000,0])
plt.show()

for el44 in data_o_lb:
   plt.plot(np.real(el44))
#plt.title(str(file_title)+str("NMR signal prom fft"), loc='center')
plt.xlabel('frequency(Hz)')
plt.ylabel('NMR signal (arb. units)')
#plt.xlim([200,250])
plt.show()


for elem in data_o_lb:
   plt.plot(elem)
#plt.title(str(file_title)+str("NMR signal"), loc='center')
plt.xlabel('frequency')
plt.ylabel('NMR signal (arb. units)')
plt.savefig(results_dir + f"FID average {file_title}.pdf")
#plt.xlim([-10000,0])
plt.show()


def fourierTransform(self, fromPos, toPos, only = []):
     self.checkToPos(toPos)
     if len(only) > 0:
        self.allFid[toPos] = np.array([fftshift(fft(self.allFid[fromPos][fidIndex])) for fidIndex in only])
     else:
        self.allFid[toPos] = np.array([fftshift(fft(fid)) for fid in self.allFid[fromPos]])
     self.frequency = np.linspace(-self.sweepWidthTD2/2,self.sweepWidthTD2/2,len(self.allFid[fromPos][0]))
     


ft_from_fid=[np.array(fftshift(fft(elem))) for elem in data_o_lb]

for elem in ft_from_fid:
   plt.plot(elem)
plt.title(str(file_title)+str("NMR signal from prom fid"), loc='center')
plt.xlabel('frequency')
plt.ylabel('NMR signal (arb. units)')
plt.xlim([180,230])
plt.show()



step=size_freq
def integral_which_max(totalft, i1, i2):
  integrals = np.zeros(len(totalft))

  for m in range(len(totalft)):
      integrals[m] = np.sum(np.real(totalft[m][i1:i2]))*step#np.sum(np.real(totalft)[m][i1:i2])
  return np.argmax(integrals),np.max(integrals),integrals #integrals[np.argmax(integrals)]




def autoPhase000(elemento, i1, i2):
  phiTest = np.linspace(0, 359, num = 360)
  integrals = np.zeros(np.size(phiTest))
  for k in range(len(integrals)):
     # print(k)
      integrals[k] =np.sum(np.real(elemento[i1:i2]*np.exp(-1j*float(phiTest[k])/180.*np.pi)))
  return phiTest[np.argmax(integrals)]

def autoPhase001(elemento22):
  phiTest = np.linspace(0, 359, num = 360)
  integrals = np.zeros(np.size(phiTest))
  for k in range(len(integrals)):
      integrals[k] =np.sum(np.real(elemento22*np.exp(-1j*float(phiTest[k])/180*np.pi)))
  return phiTest[np.argmax(integrals)]


auto_list_values= [autoPhase000(elemento,i1, i2) for elemento in ft_from_fid]

phimax=auto_list_values
print(auto_list_values)

fft100=[ft_from_fid[k]*np.exp(-1j*float(auto_list_values[k])/180.*np.pi) for k in range(len(auto_list_values))] 


auto_list_values2= [autoPhase000(elemento,i1,i2) for elemento in ft_from_fid]
phimax2 = auto_list_values2
phimax_de_fft200 = auto_list_values[integral_which_max(ft_from_fid, i1, i2)[0]]
fft109 = [ft_from_fid[k]*np.exp(-1j*float(phimax_de_fft200)/180.*np.pi) for k in range(len(ft_from_fid))]

fft200=[fft109[k]*np.exp(-1j*float(auto_list_values[k])/180.*np.pi) for k in range(len(auto_list_values))] #*np.exp(-1j*float(phimax[m]))
#fft400=[ft_from_fid[k] for k in range(len(auto_list_values)) ] 
fft200_integral=[np.sum(np.real(fft200[m][i1:i2])) for m in range(len(fft200))] #*np.exp(-1j*float(phimax[m]))


freq2 = np.linspace(i1*size_freq,i2*size_freq,i2-i1)

#
for i in range(len(fft200)): #(len(fft101)):

  plt.plot(data_freq_p[0],np.real(fft200[i]))#[i1:i2]
 # plt.title(str(file_title)+str("NMR signal fft200"), loc='center')
 # plt.xlim([-10000,0])
 # plt.xlim([200,240]) #205,225 #210,225
 # plt.ylim([-0.5*10**8,1.5*10**8])
plt.show()


for i in range(len(fft109)): #(len(fft101)):

  plt.plot(data_freq_p[0],np.real(fft200[i]))#[i1:i2]
 # plt.title(str(file_title)+str("") + str("FFT average"), loc='center')
# plt.xlim([i1,i2])
plt.xlim([-9000,-1000]) #205,225 #210,225
plt.ylim([-1*10**7,7*10**7])
#plt.title(str(file_title)+str("_average"), loc='center')
plt.xlabel('frequency (HZ)') #'sweep field values (mT/ms)'
plt.ylabel('NMR signal (arb. units)')
  #plt.grid()
plt.savefig(results_dir + f"FFT average {file_title}.pdf")
plt.show()



plt.scatter(no_repeat_B,fft200_integral)
plt.plot(no_repeat_B,fft200_integral,'red')
plt.xticks(no_repeat_B)
#\plt.title(str(file_title)+str("_average"), loc='center')
plt.xlabel(r'Offset field  $B_{0}$(mT)') #'sweep field values (mT/ms)'
plt.ylabel('NMR signal (arb. units)')
#plt.grid()
plt.savefig(results_dir + f" average {file_title}.pdf")
plt.show()

with open("integ_id273.txt", "w") as output:
    output.write(str(fft200_integral))
    output.write("field")
    output.write(str(no_repeat_B))