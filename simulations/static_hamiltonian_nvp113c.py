#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 15:24:54 2022

@author: toledo
"""

import os
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import time
from functools import reduce
from class_aux_oper_1 import aux_oper_1 as a_o
from scipy.constants import Avogadro
from scipy.constants import pi
import scipy.special

script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, ' Output_qutip_dd_n_d_henshaw_a13c_0.2_C_dip_c_0.08_t/')
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)
 
theta = 45/(180/np.pi)
theta2 = 180-45/(180/np.pi)
delta_ang =0/(180/np.pi)
delta_ang_2 =0/(180/np.pi)
v =45/(180/np.pi)
v2 =180-45/(180/np.pi)


#Johannes function
def average_defect_distance(concentration):
    return (-np.log(0.5)*12.011 * 1000000 * 3 * 1E21/(concentration * scipy.constants.Avogadro * 3.51 * 4 * scipy.constants.pi))**(1/3)


def compute_hamiltonian_dipolar(B_list):
  distance_defect = 4.8*10**(-9)#average_defect_distance(10)*(10**(-9))
  ps=[]
  hp = 6.626070*10**(-34)
  h_bar =6.626070*10**(-34)/(2*np.pi) #J*s
  a_13c_nv =2#12.5*10**6/(52*np.pi) #12.5*10**6/(2*np.pi) 
  #a_13c_p1 =-0.5*10**6
  sm = destroy(3*2*3*2)
  gyrom_13c =10.705*10**(-3)
  gyrom= -28.03#Hz/T
  D_xx = 0
  D_yy =0
  D_zz =2870
  A_xx = 85*10**6 
  A_yy =85*10**6  
  A_zz =114*10**6  
  g_n = 2
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
  
  #Ix_13c= qutip.tensor(qeye(3),qutip.jmat(0.5,"x"),qeye(2))  
  #Iy_13c= qutip.tensor(qeye(3),qutip.jmat(0.5,"y"),qeye(2)) 
  #Iz_13c= qutip.tensor(qeye(3),qutip.jmat(0.5,"z"),qeye(2)) 


  S_nv_vector = [Sx_nv,Sy_nv,Sz_nv]
  S_p1_vector = [Sx_p1,Sy_p1,Sz_p1]
  I_13c_vector = [Ix_13c,Iy_13c,Iz_13c]
  S_p1_minus = Sx_p1 - 1.0j*Sy_p1
  S_p1_plus = Sx_p1 + 1.0j*Sy_p1
  S_nv_minus = Sx_nv - 1.0j*Sy_nv
  S_nv_plus = Sx_nv + 1.0j*Sy_nv
  I_13c_minus = Ix_13c - 1.0j*Iy_13c
  I_13c_plus = Ix_13c + 1.0j*Iy_13c
  


  A_o_no_2 = a_o.matrixrot(A_o_no,theta, 0)
  #D_matrix_2 = matrixrot(D_matrix,9,0)
  print(A_o_no_2)
  rx=ry=rz=1#np.sqrt(1/2) #m 3.39411255*10**(-9)
  r = [0,0,rz]
  d_val = np.sqrt(sum([np.dot(r[i],r[i]) for i in range(len(r))])) #2*10**(-9) # meters
  print(d_val)
  r_unit = r/d_val
  print(r_unit)
  rx2=ry2=rz2=-1#np.sqrt(1/2) #m 3.39411255*10**(-9)
  r2 = [0,0,rz2]
  d_val2 = np.sqrt(sum([np.dot(r2[i],r2[i]) for i in range(len(r2))])) #2*10**(-9) # meters
  print(d_val2)
  r_unit2 = r2/d_val2
  print(r_unit2)
  
  #np.sqrt(sum([np.dot(r_unit[i],r_unit[i]) for i in range(len(r))]))

  u_Bohr = 9.274*10**(-24 ) #9.274E-24 J/T 
  u_0 = 1.2566*10**(-6) #4 pi*10-7 N/A**2
  #J.s
  #gyrom_rp1 = 2
  idx=0 
  C_dip_c =0.08#0.5*10**6 #2*10**6/(2*np.pi) 0.145 # 1 u_0*gyrom*gyrom*h_bar/(2*distance_defect**3)
  #C_dip_c = u_0*4*u_Bohr**2/(4*h_bar*np.pi*(d_val**3))          # u_0*h_bar*np.pi*gyrom*gyrom/(d_val**3) u_0*gyrom*gyrom*h_bar*np.pi/(d_val**3)
  C_dip = C_dip_c #*10**(-6)

  print("CDIP")
  print(C_dip)
  #evals_mat = n 
  evals_mat = np.zeros((len(B_list),3*2*2))
  expectSz_p1=dict()
  expectIz_14N=dict()
  expectSz_nv=dict()
  expectIz_13c=dict()
  for i in range(12):
   expectSz_p1[i]=[]
   expectIz_14N[i]=[] 
   expectSz_nv[i]=[]
   expectIz_13c[i]=[] 
  for B_z in B_list :
      B_x = 0
      B_y = 0
      B_vector = [B_x,B_y,B_z] 
      #3*a_o.vector_vector_product(I_13c_vector,r)*a_o.vector_vector_product(S_nv_vector,r
      H_zeeman_13c = gyrom_13c*a_o.vector_vector_product(B_vector,I_13c_vector)
      H_zeeman_z =gyrom*a_o.vector_vector_product(B_vector,S_nv_vector)  #utip.Qobj(gyrom*np.dot(B_z,Sz_nv)    a_o.vector_matrix_vector_product(S_nv_vector,D_matrix_2,S_nv_vector)
      S_D_S_product_nv =  D_zz*S_nv_vector[2]*S_nv_vector[2] 
      hyp_nv_p1= 3*a_o.vector_vector_product(S_nv_vector,r_unit)*a_o.vector_vector_product(S_p1_vector,r_unit)-a_o.vector_vector_product(S_nv_vector,S_p1_vector)
      hyp_nv_13c = 3*a_o.vector_vector_product(I_13c_vector,r_unit2)*a_o.vector_vector_product(S_nv_vector,r_unit2) -a_o.vector_vector_product(I_13c_vector,S_nv_vector)
      hyp_nv_13c_2=a_o.vector_vector_product(S_nv_vector,I_13c_vector) 
      hyp_p1_13c = 3*a_o.vector_vector_product(I_13c_vector,r)*a_o.vector_vector_product(S_p1_vector,r)-a_o.vector_vector_product(I_13c_vector,S_p1_vector)
      
      H_dip_nv_p1 = (1-3*(np.cos(theta))**2)*(Sz_nv*Sz_p1-0.25*(S_nv_plus*S_p1_minus + S_nv_minus*S_p1_minus))-1.5*np.sin(theta)*np.cos(theta)*np.exp(-1j*v)*(S_nv_plus*Sz_p1+Sz_nv*S_p1_plus)-1.5*np.sin(theta)*np.cos(theta)*np.exp(1j*v)*(S_nv_minus*Sz_p1+Sz_nv*S_p1_minus)
      -0.75*np.exp(-2j*v)*np.sin(theta)**2*(S_nv_plus*S_p1_plus ) -0.75*np.exp(2j*v)*np.sin(theta)**2*(S_nv_minus*S_p1_minus )
      H_dip_nv_13c = (1-3*(np.cos(theta2))**2)*(Sz_nv*Iz_13c-0.25*(S_nv_plus*I_13c_minus + S_nv_minus*I_13c_plus ))-1.5*np.sin(theta2)*np.cos(theta2)*np.exp(-1j*v2)*(S_nv_plus*Iz_13c+Sz_nv*I_13c_plus)-1.5*np.sin(theta2)*np.cos(theta2)*np.exp(1j*v2)*(S_nv_minus*Iz_13c+Sz_nv*I_13c_minus)
      -0.75*np.exp(-2j*v2)*np.sin(theta2)**2*(S_nv_plus*I_13c_plus) -0.75*np.exp(2j*v2)*np.sin(theta2)**2*(S_nv_minus*I_13c_minus )
     
      
  
     # H_dip_nv_p1_m = (1-3*(np.cos(theta))**2)*(np.matmul(Sz_nv,Sz_p1)-0.25*(np.matmul(S_nv_minus,S_p1_plus) + np.matmul(S_nv_plus,S_p1_minus)))
      H_dipolar_NV_p1_1 = reduce(lambda x,y : x + y , [S_p1_vector[i]*S_nv_vector[i] for i in range(len(B_vector))])
      H_dipolar_NV_p1_2= -(3*(np.cos(delta_ang))**2)*sum([S_p1_vector[i]*r_unit[i] for i in range(3)])*sum([S_nv_vector[i]*r_unit[i] for i in range(3)])
      H_total_dd = -C_dip*(hyp_nv_p1)#*H_dipolar_NV_p1_1 + C_dip*H_dipolar_NV_p1_2
      #a_o.vector_matrix_vector_product(S_nv_vector,D_matrix_2,S_nv_vector) #sum([np.dot(S_nv_vector[k],D_S_product_nv_13c[k]) for k in range(3)])
      H_hyp_nv_13c = a_13c_nv*(I_13c_vector[0]*S_nv_vector[0]+I_13c_vector[1]*S_nv_vector[1] + I_13c_vector[2]*S_nv_vector[2])#reduce(lambda x,y : x + y , [a_13c_nv_m[i][i]*S_nv_vector[i]*I_13c_vector[i] for i in range(len(B_vector))])
    
      H_zeeman_z_p1 = gyrom*a_o.vector_vector_product(B_vector,S_p1_vector)    #
      H_dipolar_NV_p1_1 = S_p1_vector[0]*S_nv_vector[0]+S_p1_vector[1]*S_nv_vector[1]+S_p1_vector[2]*S_nv_vector[2]   #reduce(lambda x,y : x + y , [S_p1_vector[i]*S_nv_vector[i] for i in range(len(B_vector))])
      H_dipolar_p1_r= S_p1_vector[0]*r_unit[0]+S_p1_vector[1]*r_unit[1]+S_p1_vector[2]*r_unit[2]#-3*sum([S_p1_vector[i]*r_unit[i] for i in range(3)])*sum([S_nv_vector[i]*r_unit[i] for i in range(3)])
      H_dipolar_NV_r= S_nv_vector[0]*r_unit[0]+S_nv_vector[1]*r_unit[1]+S_nv_vector[2]*r_unit[2]#
    #  H_hyp_p1_n14 = a_o.vector_matrix_vector_product(I_14N_vector,A_o_no_2,S_p1_vector)     
       #+ H_NV_zeeman +H_hyp_nv_13c #qutip.Qobj(H_NV_zfs)
     
      H_total =-H_zeeman_z_p1 -H_zeeman_13c-H_zeeman_z + S_D_S_product_nv -C_dip*hyp_nv_p1-a_13c_nv*hyp_nv_13c
      evals, ekets = H_total.eigenstates() #H_total_nv_13c.eigenstates()
     
      evals_mat[idx,:] = np.real(evals)
      n_of_ew=len(evals)
      for i in range(n_of_ew):
        #  print(ekets[i])
          #expectSz_p1[i].append(ekets[i].overlap(Sz_p1))
          expectSz_p1[i].append(ekets[i].overlap(Sz_p1))
          expectSz_nv[i].append(ekets[i].overlap(Sz_nv))
          expectIz_13c[i].append(ekets[i].overlap(Iz_13c))
         # expectIz_14N[i].append(ekets[i].overlap(Iz_14N))
      idx += 1
     
     
      
  return evals_mat,H_total,ekets,Sz_p1,Iz_13c,Sz_nv,expectSz_p1,expectSz_nv,expectIz_13c

xticks= np.linspace(48,54,60)
#xticks2=np.linspace(0.0,0.060,1800)
B_list = np.linspace(48, 54,600)  # atom 1 frequency range

#evals_mat_p = compute_hamiltonian(B_list)
evals_mat_p2 = compute_hamiltonian_dipolar(B_list)

#7,8,9,10

fig, ax = plt.subplots(figsize=(20,10))
#ax.plot(B_list , evals_mat_p , 'b')

ax.plot(B_list, evals_mat_p2[0],'blue')

#ax.plot(B_list, evals_mat_p2[:,1],'red')
#ax.plot(B_list, evals_mat_p2 + evals_mat_p,'g')
#ax.set_xticks(xticks2)
ax.set_xlabel('Magnetic field (mT)')
#ax.set_ylim([500,1000])
#ax.set_xlim([0.0495,0.052]) #ax.set_xlim(np.linspace(0.0495,0.052,300))
ax.set_ylabel('Eigenenergies (Hz)')
ax.set_title('Energy diagram');
plt.savefig(results_dir + f"out_energy_levels_tilted_{theta}_zoom.pdf")
plt.show()
#expectSz_p1,expectIz_14N,expectSz_nv,expectIz_13c
fig0, ax0 = plt.subplots(figsize=(20,10))
plt.rcParams.update({'font.size': 14})
n_of_ew=12
evals_mat_p3=evals_mat_p2[0]
evals_mat_p6=evals_mat_p2[6]
evals_mat_p7=evals_mat_p2[7]
evals_mat_p8=evals_mat_p2[8]

Ham_ = evals_mat_p2[1]

state2= basis(2*3*2,2) #0,1/2,1/2
state3= basis(2*3*2,3)#0,1/2,-1/2
state10= basis(2*3*2,10)#-1,-1/2,+1/2
state11= basis(2*3*2,11)#-1,-1/2,-1/2


for n in range(n_of_ew):
    ax0.plot(B_list, evals_mat_p3[:,n], lw=0.1)
    colorencodedplot=ax0.scatter(B_list, evals_mat_p3[:,n],marker="o",s=100, c=evals_mat_p6[n],cmap='winter',vmin=-0.5,vmax=0.5)
#ax.set_xlim([46,56]) [0.62e+09,0.86e+09])
#ax.set_ylim([1000,-1000])
ax0.set_xticks(xticks)
#ax.xticks(np.linspace(48,54,60))
#ax0.set_ylim([0.6e+09,0.8e+09])
ax0.set_ylim([716,718.5]) #[0.7175e+9,0.7177e+9]
ax0.set_xlim([51.1,51.3])
ax0.set_xlabel('magnetic field (mT)')
ax0.set_ylabel('Eigenenergies (Hz)')
ax0.set_title('P1')
cbar=fig0.colorbar(colorencodedplot)
cbar.set_label("m_z", labelpad=+1)
plt.savefig(results_dir + f"m_p1{theta}_zoom.pdf")
plt.show()

fig1, ax1 = plt.subplots(figsize=(20,10))
plt.rcParams.update({'font.size': 14})
n_of_ew=12
for n in range(n_of_ew):
    ax1.plot(B_list, evals_mat_p3[:,n], lw=0.1)
    colorencodedplot=ax1.scatter(B_list, evals_mat_p3[:,n],marker="o",s=100, c=evals_mat_p7[n],cmap='spring',vmin=-1.0,vmax=1.0)
#ax.set_xlim([46,56])
#ax.set_ylim([1000,0])
ax1.set_xticks(xticks)
ax1.set_ylim([716,78.5]) 
ax1.set_xlim([51.1,51.3])
ax1.set_xlabel('magnetic field (mT)')
ax1.set_ylabel('Eigenenergies (Hz)')
ax1.set_title('NV')
cbar=fig1.colorbar(colorencodedplot)
cbar.set_label("m_z", labelpad=+1)
plt.savefig(results_dir + f"m_NV{theta}_zoom.pdf")
plt.show()

fig2, ax2 = plt.subplots(figsize=(20,10))
plt.rcParams.update({'font.size': 14})
n_of_ew=12
for n in range(n_of_ew):
    ax2.plot(B_list, evals_mat_p3[:,n], lw=0.1)
    colorencodedplot=ax2.scatter(B_list, evals_mat_p3[:,n],marker="o",s=100, c=evals_mat_p8[n],cmap='spring',vmin=-0.5,vmax=0.5)
#ax.set_xlim([46,56])
#ax.set_ylim([1000,-1000])
#ax.set_ylim([0.62e+09,0.86e+09])
#ax.set_xlim([0.048,0.054])
ax2.set_xticks(xticks) #0.717325e+9,0.71734e+9
ax2.set_ylim([716,718.5])#0.71726e+9,0.71733e+9]
ax2.set_xlim([51.0,51.3])
ax2.set_xlabel('magnetic field (mT)')
ax2.set_ylabel('Eigenenergies (Hz)')
ax2.set_title('13C')
cbar=fig2.colorbar(colorencodedplot)
cbar.set_label("m_z", labelpad=+1)
plt.savefig(results_dir + f"m_13c{theta}_zoom.pdf")
plt.show()


fig3, ax3 = plt.subplots(figsize=(20,10))
plt.rcParams.update({'font.size': 14})
n_of_ew=12
for n in range(n_of_ew):
    ax3.plot(B_list, evals_mat_p3[:,n], lw=2.0)
    colorencodedplot=ax3.scatter(B_list, evals_mat_p3[:,n],marker="o",s=60, c=evals_mat_p6[n],cmap='cool',vmin=-0.5,vmax=0.5)
#ax5.set_xticks(xticks)
#ax4.set_ylim([620,860])
#ax5.set_xlim([48,54])
ax3.set_xlabel('magnetic field (mT)')
ax3.set_ylabel('S(P1)')
plt.savefig(results_dir + f"m_p1{theta}_zoom_expect.pdf")
plt.show()

fig4, ax4 = plt.subplots(figsize=(20,10))
plt.rcParams.update({'font.size': 14})
n_of_ew=12
for n in range(n_of_ew):
    ax4.plot(B_list, evals_mat_p7[n], lw=2.0)
#ax4.set_xticks(xticks)
#ax5.set_ylim([620,860])
ax4.set_xlim([48,54])
ax4.set_xlabel('magnetic field (mT)')
ax4.set_ylabel('S(NV)')
plt.savefig(results_dir + f"m_NV{theta}_zoom_expect.pdf")
plt.show()


fig5, ax5 = plt.subplots(figsize=(20,10))
plt.rcParams.update({'font.size': 14})
n_of_ew=12
for n in range(n_of_ew):
    ax5.plot(B_list, evals_mat_p8[n], lw=2.0)
#ax4.setticks(xticks)
#ax4.set_ylim([620,860])
ax5.set_xlim([48,54])
ax5.set_xlabel('magnetic field (mT)')
ax5.set_ylabel('S(13C)')
plt.savefig(results_dir + f"m_13C{theta}_zoom_expect.pdf")
plt.show()








        
