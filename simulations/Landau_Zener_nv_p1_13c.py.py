#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 14:07:13 2022

@author: toledo
"""


import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import time
from functools import reduce
from class_aux_oper_3 import aux_oper_1 as a_o
import os
import pandas as pd
import scipy
from scipy.constants import Avogadro
from scipy.constants import pi
import scipy.special


theta = 0/(180/np.pi)
delta_ang =0/(180/np.pi)
delta_ang_2 =0/(180/np.pi)


h_bar =6.626070*10**(-34)/(2*np.pi) 
sweep_f=1000#mT/s #51,52,---10**(-10) 0.4*10**(-3) , 51.1,51.3 ---10**(-7) 0.4*10**(-3)
B_g = 51
B_fin =51.5
B_y = B_x = 0
a_13c_nv = 0.92*10**6 #2*10**6#0.06*10**6
C_dip_c_ = 1*10**6#'0.35*10**6
B_range = B_fin - B_g# Henshaw apaper
period_sweep = B_range/sweep_f
psi0 = basis(2*3*2,3)
steps_LZ=200000
# Johannes fct
def average_defect_distance(concentration):
    return (-np.log(0.5)*12.011 * 1000000 * 3 * 1E21/(concentration * scipy.constants.Avogadro * 3.51 * 4 * scipy.constants.pi))**(1/3)

concentration=30
distance_defect = 4.8*10**(-9)# average_defect_distance(concentration)


h_bar =6.626070*10**(-34)/(2*np.pi) #j.s
#0.5*10**6#0.5*10**6   #MHz      

gyrom_13c =10.705*10**3

gyrom=-28.03*10**6 #MHZ/mT
g_n =2
u_0 = 1.2566*10**(-6)
D_xx = 0
D_yy =0
D_zz =2870*10**6
A_xx = 85*10**6
A_yy =85*10**6
A_zz =114*10**6
u_Bohr = 9.274*10**(-24 ) #J/T

  #h_bar = 1#6.626070*10**(-34)/(2*np.pi) 
A_o_no = np.diag([A_xx,A_yy,A_zz])
  #a_13c_nv_m = np.diag([a_13c_nv,a_13c_nv,a_13c_nv])
gyrom_arr = [gyrom,gyrom,gyrom]
a_13c_nv_m = np.diag([a_13c_nv,a_13c_nv,a_13c_nv])
D_matrix =  np.diag([D_xx,D_yy,D_zz])
D_matrix_2 = a_o.matrixrot(D_matrix,theta,0)

Sx_nv= qutip.tensor(qeye(2),qutip.jmat(1,"x"),qeye(2))
Sy_nv= qutip.tensor(qeye(2),qutip.jmat(1,"y"),qeye(2))
Sz_nv= qutip.tensor(qeye(2),qutip.jmat(1,"z"),qeye(2))
Sx_p1=qutip.tensor(qutip.jmat(0.5,"x"),qeye(3),qeye(2))
Sy_p1= qutip.tensor(qutip.jmat(0.5,"y"),qeye(3),qeye(2)) 
Sz_p1=qutip.tensor(qutip.jmat(0.5,"z"),qeye(3),qeye(2))
Ix_13c=qutip.tensor(qeye(2),qeye(3),qutip.jmat(0.5,"x"))
Iy_13c=qutip.tensor(qeye(2),qeye(3),qutip.jmat(0.5,"y")) 
Iz_13c=qutip.tensor(qeye(2),qeye(3),qutip.jmat(0.5,"z"))
S_nv_vector = [Sx_nv,Sy_nv,Sz_nv]
S_p1_vector = [Sx_p1,Sy_p1,Sz_p1]
I_13c_vector = [Ix_13c,Iy_13c,Iz_13c]

S_p1_minus = Sx_p1 - 1.0j*Sy_p1
S_p1_plus = Sx_p1 + 1.0j*Sy_p1
S_nv_minus = Sx_nv - 1.0j*Sy_nv
S_nv_plus = Sx_nv + 1.0j*Sy_nv
I_13c_minus = Ix_13c - 1.0j*Iy_13c
I_13c_plus = Ix_13c + 1.0j*Iy_13c

rx=ry=rz=1#np.sqrt(1/2) #m 3.39411255*10**(-9)
r = [rx,0,rz]
d_val = np.sqrt(sum([np.dot(r[i],r[i]) for i in range(len(r))])) #2*10**(-9) # meters
print(d_val)
r_unit = r/d_val
print(r_unit)
rx2=ry2=rz2=-1#np.sqrt(1/2) #m 3.39411255*10**(-9)
r2 = [rx2,0,rz2]
d_val2 = np.sqrt(sum([np.dot(r2[i],r2[i]) for i in range(len(r2))])) #2*10**(-9) # meters
print(d_val2)
r_unit2 = r2/d_val2
print(r_unit2)
  
#H_dipolar_NV_p1_1 = reduce(lambda x,y : x + y , [S_p1_vector[i]*S_nv_vector[i] for i in range(3)])
#H_dipolar_NV_p1_2= -(3*(np.cos(delta_ang))**2)*sum([S_p1_vector[i]*r_unit[i] for i in range(3)])*sum([S_nv_vector[i]*r_unit[i] for i in range(3)])
#H_total_dd = C_dip_c_*H_dipolar_NV_p1_1 + C_dip_c_*H_dipolar_NV_p1_2
hyp_nv_p1= 3*a_o.vector_vector_product(S_nv_vector,r_unit)*a_o.vector_vector_product(S_p1_vector,r_unit)-a_o.vector_vector_product(S_nv_vector,S_p1_vector)
hyp_nv_13c = 3*a_o.vector_vector_product(I_13c_vector,r_unit2)*a_o.vector_vector_product(S_nv_vector,r_unit2) -a_o.vector_vector_product(I_13c_vector,S_nv_vector)
#H_hyp_nv_13c = a_13c_nv*(I_13c_vector[0]*S_nv_vector[0]+I_13c_vector[1]*S_nv_vector[1] + I_13c_vector[2]*S_nv_vector[2])
#H_dip_nv_p1 = (1-3*(np.cos(theta))**2)*(Sz_nv*Sz_p1-0.25*(S_nv_plus*S_p1_minus + S_nv_minus*S_p1_minus))
#H_dip_nv_13c = (1-3*(np.cos(theta))**2)*(Iz_13c*Sz_nv-0.25*(I_13c_minus*S_nv_plus + I_13c_plus*S_nv_minus ))
def qubit_integrate(B_x,B_y,B_g, gamma1, gamma2, psi0,sweep_f,tlist):
  ps=[]
  B_z=B_g
  C_dip = C_dip_c_
  B_vector = [B_x,B_y,B_z] 
  H_zeeman_z_p1 = gyrom*a_o.vector_vector_product(B_vector,S_p1_vector)  
  H_zeeman_13c = gyrom_13c*a_o.vector_vector_product(B_vector,I_13c_vector)
  H_zeeman_z =gyrom*a_o.vector_vector_product(B_vector,S_nv_vector)  #utip.Qobj(gyrom*np.dot(B_z,Sz_nv)    a_o.vector_matrix_vector_product(S_nv_vector,D_matrix_2,S_nv_vector)
  S_D_S_product_nv =  D_zz*S_nv_vector[2]*S_nv_vector[2] 
  
#  hyp_nv_p1= 3*a_o.vector_vector_product(S_nv_vector,r_unit2)*a_o.vector_vector_product(S_p1_vector,r_unit2)-a_o.vector_vector_product(S_nv_vector,S_p1_vector)
 # hyp_nv_13c = 3*a_o.vector_vector_product(I_13c_vector,r_unit)*a_o.vector_vector_product(S_nv_vector,r_unit) -a_o.vector_vector_product(I_13c_vector,S_nv_vector)    
   #
  H0 = -H_zeeman_z_p1 -H_zeeman_13c-H_zeeman_z + S_D_S_product_nv +C_dip*hyp_nv_p1+a_13c_nv*hyp_nv_13c
  #-H_zeeman_z_p1 -H_zeeman_13c-H_zeeman_z + S_D_S_product_nv -a_13c_nv*hyp_nv_13c-C_dip*hyp_nv_p1
  #+ C_dip*H_dip_nv_p1-a_13c_nv*H_dip_nv_13c 
  
  #- H_total_dd+ H_hyp_nv_13c # -C_dip*hyp_nv_p1-a_13c_nv*hyp_nv_13c#+  H_total_dd +H_hyp_nv_13c   
  H1 =(sweep_f)*(-gyrom*Sz_nv -gyrom*Sz_p1-gyrom_13c*Iz_13c ) #---, --+, no funciona -gyrom_13c*Iz_13c
  
  H = [H0, [H1, 't']]

  output = mesolve(H, psi0.unit(), tlist, [], [Sz_p1,Sz_nv, Iz_13c],{})  
  states_result = mesolve(H, psi0.unit(), tlist, {})
  #outputstate = mesolve(H,psi0,tlist,{})
  out_t = list(output.expect[0]),list(output.expect[1]),list(output.expect[2]),list(states_result.states)
  return out_t
 
# T/s


Sz_p1=qutip.tensor(qutip.jmat(0.5,"z"),qeye(3),qeye(2))
Iz_13c=qutip.tensor(qeye(2),qeye(3),qutip.jmat(0.5,"z"))
Sz_nv= qutip.tensor(qeye(2),qutip.jmat(1,"z"),qeye(2))


start_time = time.time()

print("period sweep")
print(period_sweep)


time_c = []
proj_z = []

tlist = np.linspace(0,period_sweep,steps_LZ)#np.linspace(0.0,100*period_sweep,1000000)
gamma1=0
gamma2=0 # which values...

p_total=qubit_integrate(0,0,B_g,gamma1, gamma2,psi0.unit(),sweep_f,tlist)

p_ex1 =p_total[0] #p1
p_ex2 =p_total[1] #nv
p_ex3 =p_total[2] #13c
 #14N
states_f = p_total[3]

script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, f'{B_g}+{B_fin}+{sweep_f}+{period_sweep}_c_basis6_dist_calc+{distance_defect}a_13c+{a_13c_nv}C_dip_c_ +{C_dip_c_}minus/')
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)
list_dict = {'p1':p_ex1, 'nv':p_ex2,'13c':p_ex3} 
df = pd.DataFrame(list_dict) 
df.to_csv(results_dir + f'data.csv', index=False)

script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, f'{B_g}+{B_fin}+{sweep_f}+{period_sweep}_c_basis6_dist_calc+{distance_defect}a_13c+{a_13c_nv}C_dip_c_ +{C_dip_c_}minus/')
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)
list_dict = {'states_f':states_f} 
df = pd.DataFrame(list_dict) 
df.to_csv(results_dir + f'eigenstates.txt', index=False)



print('time elapsed = ' + str(time.time() - start_time))

plt.figure(figsize=(20,10))
plt.rcParams.update({'font.size': 12})

blist=B_g+(sweep_f)*tlist

co_min_1 = np.real(p_ex2)#-np.real(p_ex2[0])
c13_val =np.real(p_ex3)#-np.real(p_ex3[0])
p1_val = np.real(p_ex1)#-np.real(p_ex1[0])

plt.plot(blist, co_min_1, 'b') #/no_r

plt.xlabel('B(mT)')
plt.ylabel(r'$\langle m_{NV} \rangle$')#plt.ylabel('Occupation probability ')
plt.title('NV')

plt.savefig(results_dir + f"nv_{B_g}_{B_fin}.pdf")
plt.show()

plt.figure(figsize=(20,10))
plt.rcParams.update({'font.size': 12})

plt.plot(blist,p1_val, 'g')

plt.xlabel('B(mT)')
plt.ylabel(r'$\langle m_{P1} \rangle$')#plt.ylabel('Occupation probability ')
plt.title(r'$P1$')
plt.savefig(results_dir + f"p1_{B_g}_{B_fin}.pdf")
plt.show()


plt.figure(figsize=(20,10))
plt.rcParams.update({'font.size': 12})
plt.plot(blist,c13_val, 'c')
#plt.xlim(51.1,51.3)
plt.xlabel('B(mT)')
plt.ylabel(r'$\langle m_{13C} \rangle$')
plt.title(r'$^{13}C$')
plt.savefig(results_dir + f"c13_{B_g}_{B_fin}.pdf")
plt.show()