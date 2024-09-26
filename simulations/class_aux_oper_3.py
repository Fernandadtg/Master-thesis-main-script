#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 12:45:53 2022

@author: u
"""

import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import time
from functools import reduce

class aux_oper_1:

 
   def vector_vector_product(vect1,vect2):
     producto = 0
     for m in range(len(vect1)):         
         producto += vect1[m]*vect2[m]
     return producto
    
    
   def vector_matrix_product(mat1, vect1):  
     result = [] 
     for i in range(len(vect1)): 
         entry = 0
         for j in range(len(mat1[0])): 
             entry +=  mat1[i][j]*vect1[j]
         result.append(entry) 
     return result


   def vector_matrix_vector_product(vector1,matrix1,vector2):
     return(aux_oper_1.vector_vector_product(vector1,aux_oper_1.vector_matrix_product(matrix1,vector2)))



   def matrixrot(matrix, theta, phi):
      theta = theta / (180/np.pi)
      phi = phi / (180/np.pi)
      D=matrix
      RotX = np.array(((1,0,0),(0,np.cos(theta),-np.sin(theta)),(0,np.sin(theta),np.cos(theta))))
      RotY = np.array(((np.cos(theta),0,np.sin(theta)),(0,1,0),(-np.sin(theta),0,np.cos(theta)))) 
      RotZ = np.array(((np.cos(phi),-np.sin(phi),0),(np.sin(phi),np.cos(phi),0),(0,0,1)))    
      DrotY = np.dot(np.linalg.inv(RotY),np.dot(D,RotY)) 
      Drot = np.dot(np.linalg.inv(RotZ),np.dot(DrotY,RotZ))
      return Drot


    
   def diff_with_error_3(list1,list2,B_some_list,error,orientation):
      diff_list_with_error = []
      B_list_vals = []
      spin_13c = []
      for i in range(len(list1)):
        if abs(list1[i]-list2[i])<=error:
            diff_list_with_error.append(abs(list1[i]-list2[i]))
            B_list_vals.append(B_some_list[i])
            spin_13c.append(orientation)
        else:
            diff_list_with_error.append(abs(list1[i])-list2[i])
            B_list_vals.append(B_some_list[i])
            spin_13c.append(0)
      return (diff_list_with_error,B_list_vals,spin_13c)



   def diff_with_error(list1,list2,B_some_list,error,orientation):
      diff_list_with_error_2 = []
      B_list_vals_2= []
      spin_13c_2 = []
      for i in range(len(list1)):
        if abs(list1[i]-list2[i])<=error:
            diff_list_with_error_2.append(abs(list1[i]-list2[i]))
            B_list_vals_2.append(B_some_list[i])  
            spin_13c_2.append(orientation)
      return (diff_list_with_error_2,B_list_vals_2,spin_13c_2)

        
