#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 14:25:43 2019

@author: peter
"""

import numpy as np
from C2P_EOS import EOS_pressFromWzTCons
from C2P_Constraints import Constraints_vector, Constraints_Jacobain_SemiAnalytic
from C2P_MatrixInv import C2P_MatrixInverse
from C2P_InitialGuess import GetWzTGuess
from C2P_Setup import log_T_array, log_rho_array, Ye_array

T_min = 10**log_T_array[1]
T_max = 10**log_T_array[-2]

def RelDiff(a,b):
    return np.abs(2*(a-b)/(a+b))

def Con2Prim(Cons_target,guesses=None,press_rel_diff=1e-9,its_max=100,verbose=False,var_rel_diff=None):
   
    WzT = GetWzTGuess(Cons_target,guesses)
    W_guess = WzT["W"]
    z_guess = WzT["z"]
    T_guess = WzT["T"]

    WzT_delta = {"W": W_guess*1e-3,
                 "z": z_guess*1e-3,
                 "T": T_guess*1e-3}
    
    p_old = EOS_pressFromWzTCons(WzT,Cons_target)
    Ws = np.zeros(its_max+1)
    zs = np.zeros(its_max+1)
    Ts = np.zeros(its_max+1)
    Ws[0] = W_guess
    zs[0] = z_guess
    Ts[0] = T_guess
    
    if verbose: print("Starting Con2Prim")
    it=0
    done=False
    while done==False:
        WzT_old = WzT.copy()
        f_vector = Constraints_vector(WzT,Cons_target)

        Jacobian = Constraints_Jacobain_SemiAnalytic(WzT,Cons_target,WzT_delta)
        J_inv = C2P_MatrixInverse(Jacobian)

        WzT_array = np.array([WzT["W"],WzT["z"],WzT["T"]])
        delta_WzT_array = np.dot(J_inv,f_vector)
        new_WzT_array = WzT_array - delta_WzT_array
    
        if new_WzT_array[0]<1.0:
            new_WzT_array[0]=1.0
        if new_WzT_array[2]<T_min:
            new_WzT_array[2] = T_min
        elif new_WzT_array[2]>T_max:
            new_WzT_array[2] = T_max

        WzT = {"W": new_WzT_array[0],
               "z": new_WzT_array[1],
               "T": new_WzT_array[2]}
        
        Ws[it+1] = new_WzT_array[0]
        zs[it+1] = new_WzT_array[1]
        Ts[it+1] = new_WzT_array[2]
    
        p_current = EOS_pressFromWzTCons(WzT,Cons_target)
        
        if var_rel_diff:
            converged = np.zeros(3,dtype=bool)
            if RelDiff(WzT["W"],WzT_old["W"])<=var_rel_diff:
                converged[0] = True
            if RelDiff(WzT["z"],WzT_old["z"])<=var_rel_diff:
                converged[1] = True
            if RelDiff(WzT["T"],WzT_old["T"])<=var_rel_diff:
                converged[2] = True
            if np.all(converged):
                done = True
        elif abs(2*(p_current-p_old)/(p_current+p_old))<=press_rel_diff:
            done=True
            if verbose: print("Successful convergence after " + str(it+1) + " iterations.")
        else:
            p_old=p_current
    
        it+=1
        if it>=its_max:
            done=True
            if verbose: print("Max it reached, terminating early")
            
        if it>=20 and it%20==0:
            WzT["W"] = np.average(Ws[10:it])
            WzT["z"] = np.average(zs[10:it])
            WzT["T"] = np.average(Ts[10:it])
            
#            min_frac = 0.9
#            max_frac = 1.1
#            del_frac = max_frac-min_frac
#            WzT = GetWzTGuess(Cons_target)
#            #WzT["W"]*= min_frac + del_frac*random()
#            #WzT["z"]*= min_frac + del_frac*random()
#            #WzT["T"]*= min_frac + del_frac*random()
#            WzT["W"] = new_WzT_array[0]
#            WzT["z"] = new_WzT_array[1]
#            WzT["T"] = new_WzT_array[2]*(min_frac + del_frac*random())
    
    return WzT, it, Ws[:it], zs[:it], Ts[:it]
