# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 13:14:27 2021

@author: aqs
"""

import sys
import redpitaya_scpi as scpi
import matplotlib.pyplot as plot
import numpy as np
#import visa
import time
import json
import pandas as pd
rp_s = scpi.scpi('172.26.29.224')


decimation=8192
wave_form  = 'SAWU'
freq = 10
ampl =0.5
offset = 0.0

rp_s.tx_txt('GEN:RST')
rp_s.tx_txt('SOUR1:FUNC ' + str(wave_form).upper())
rp_s.tx_txt('SOUR1:VOLT:OFFS' + str(offset))
rp_s.tx_txt('SOUR1:VOLT ' + str(ampl))
rp_s.tx_txt('SOUR1:FREQ:FIX ' + str(freq))



rp_s.tx_txt('SOUR1:BURS:NCYC 2')
rp_s.tx_txt('SOUR1:BURS:STAT ON')
rp_s.tx_txt('SOUR1:TRIG:SOUR INT')
rp_s.tx_txt('SOUR1:TRIG:IMM')
rp_s.tx_txt('OUTPUT1:STATE ON')

#Enable output
rp_s.tx_txt('ACQ:DEC 8192' )
rp_s.tx_txt('ACQ:TRIG:DLY 8192')
rp_s.tx_txt('ACQ:SOUR1:GAIN HV')


rp_s.tx_txt('ACQ:START')
time.sleep(10) 
rp_s.tx_txt('ACQ:TRIG NOW') 
rp_s.tx_txt('ACQ:TRIG CH1_PE')


while 1:
  rp_s.tx_txt('ACQ:TRIG:STAT?')
  if rp_s.rx_txt() == 'TD':
     break
rp_s.tx_txt('ACQ:SOUR1:DATA?')
rp_s.tx_txt('ACQ:STOP')
rp_s.tx_txt('OUTPUT1:STATE OFF')
buff_string = rp_s.rx_txt()
buff_string = buff_string.strip('{}\n\r').replace(" ", "").split(',')
buff = list(map(float,buff_string))
buff2=np.array(buff)

B=(buff2-2.128)/31.25
x1 = np.linspace(0, 1.074, len(buff)) #1.074


c = [x1,buff]


with open("freq"+str(freq)+"_"+ "Voltage"+str(ampl)+".txt", 'w') as file:
  for m in zip(*c):
    file.write("{0}\t{1}\n".format(*m))
  


print(np.max(buff))
print(np.min(buff))

amp_obtained = np.max(buff)-np.min(buff)
d=[freq, ampl, amp_obtained]

with open("Diff_freq"+str(freq)+"_"+ "Voltage"+str(ampl)+ ".txt", 'a') as file2:
     #file2.write("Frequency   Voltage   Amplitude \n ")
    #file2.write('{}\n'.format(amplitude))
     file2.write("{0}\t{1}\t{2}\n".format(*d))
plot.plot(x1,B)
plot.xlabel('time')
plot.ylabel('B')
plot.savefig("freq"+str(freq)+"_"+ "Voltage"+str(ampl)+".pdf")
plot.show()





