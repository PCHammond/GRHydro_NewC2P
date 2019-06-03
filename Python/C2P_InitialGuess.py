#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 15:17:31 2019

@author: peter
"""

import copy
import numpy as np
from random import betavariate, random
from C2P_EOS import EOS_pressFromPrim, EOS_epsFromPrim
from C2P_Setup import log_T_array, log_rho_array, Ye_array

class Targets(object):
    def copy(self):
        obj_copy = copy.deepcopy(self)
        return obj_copy 

def GetWzTGuess(Cons_target,guesses):
    D = Cons_target["D"]
    DYe = Cons_target["DYe"]
    Ye_found = DYe/D
    
    if guesses == None:
        W_guess = 1.0
        T_guess = 1.0
    else:
        vsq_guess = guesses.vsq
        W_guess = 1/np.sqrt(1-vsq_guess)
        T_guess = guesses.T

    rho_guess = D/W_guess

    press_guess = EOS_pressFromPrim(T_guess,rho_guess,Ye_found)
    eps_guess = EOS_epsFromPrim(T_guess,rho_guess,Ye_found)
    h_guess = 1 + eps_guess + press_guess/rho_guess
    z_guess = rho_guess*h_guess*(W_guess**2)

    WzT_guess = {"W": W_guess,
                 "z": z_guess,
                 "T": T_guess}
    
    return WzT_guess

def GetInitialGuess(rho=None,velocity=[],T=None,Ye=None):
    targets = Targets()
    if rho:
        targets.rho = rho
    else:
        ### Density ###
        log_rho_min = -11
        log_rho_max = -3
        alpha = 2
        beta = 7
        log_rho = log_rho_min + ((1-betavariate(alpha,beta))*(log_rho_max-log_rho_min))
        targets.rho = 10**log_rho
    
    if velocity!=[]:
        targets.v1 = velocity[0]
        targets.v2 = velocity[1]
        targets.v3 = velocity[2]
    else:
        ### Velocity ###
        alpha = 2
        beta = 3
        speed = 0.5*betavariate(alpha,beta)
        theta = 2*np.pi*random()
        phi   = np.arccos(-1.0 + 2*random())
        targets.v1 = speed*np.cos(theta)*np.sin(phi)
        targets.v2 = speed*np.sin(theta)*np.sin(phi)
        targets.v3 = speed*np.cos(phi)
    
    if T:
        targets.T = T
    else:
        ### Temperature ###
        log_T_min = log_T_array[0] + 0.5
        log_T_max = log_T_array[-1] - 0.5
        alpha = 2
        beta = 7
        log_T = log_T_min + (betavariate(alpha,beta)*(log_T_max-log_T_min))
        targets.T = 10**log_T
        
    if Ye:
        targets.Ye = Ye
    else:
        ### Ye ###
        Ye_min = Ye_array[2]
        Ye_max = Ye_array[-3]
        alpha = 2
        beta = 11
        targets.Ye = Ye_min + (betavariate(alpha,beta)*(Ye_max-Ye_min))

    ### Primitive ###
    targets.vsq = np.linalg.norm([targets.v1,targets.v2,targets.v3])**2
    targets.W = 1/(np.sqrt(1-targets.vsq))
    
    targets.press = EOS_pressFromPrim(targets.T,targets.rho,targets.Ye)
    targets.eps = EOS_epsFromPrim(targets.T,targets.rho,targets.Ye)
    targets.h = 1 + targets.eps + targets.press/targets.rho
    targets.z = targets.h * targets.rho * (targets.W**2)

    ##### Conserved #####
    D   = targets.rho*targets.W
    S_1 = (targets.rho*targets.h)*(targets.W**2)*targets.v1
    S_2 = (targets.rho*targets.h)*(targets.W**2)*targets.v2
    S_3 = (targets.rho*targets.h)*(targets.W**2)*targets.v3
    tau = (targets.rho*targets.h)*(targets.W**2) - targets.press - D
    DYe = D*targets.Ye
    S2  = ((targets.rho*targets.h)**2)*(targets.W**4)*(targets.vsq)

    Cons_target = {"D"   : D,
                   "S_1"  : S_1,
                   "S_2"  : S_2,
                   "S_3"  : S_3,
                   "tau" : tau,
                   "DYe" : DYe,
                   "S2"  : S2}    
    
    return targets, Cons_target

def GuessFromTargets(targets,min_frac=np.ones(3)*0.98,max_frac=np.ones(3)*1.02):
    guesses = Targets()
    
    rho_target = targets.rho
    rand = min_frac[0]+(max_frac[0]-min_frac[0])*random()
    guesses.rho = rho_target*rand
    
    vsq_target = targets.vsq
    rand = min_frac[1]+(max_frac[1]-min_frac[1])*random()
    guesses.vsq = vsq_target*rand
    guesses.W = 1/(np.sqrt(1-guesses.vsq))
    
    T_target = targets.T
    rand = min_frac[2]+(max_frac[2]-min_frac[2])*random()
    guesses.T = T_target*rand
    
    guesses.press = EOS_pressFromPrim(guesses.T,guesses.rho,targets.Ye)
    guesses.eps = EOS_epsFromPrim(guesses.T,guesses.rho,targets.Ye)
    guesses.h = 1 + guesses.eps + guesses.press/guesses.rho
    guesses.z = guesses.h * guesses.rho * (guesses.W**2)
    
    return guesses

##### Primitives #####
rho_target = 2.8358299540575783e-09
v1_target  = 0.1
v2_target  = -0.2
v3_target  = 0.15
Ye_target  = 0.315
T_target   = 4.166582513489608

vsq_target = np.linalg.norm([v1_target,v2_target,v3_target])**2
W_target = 1/(np.sqrt(1-vsq_target))

press_target = EOS_pressFromPrim(T_target,rho_target,Ye_target)
eps_target = EOS_epsFromPrim(T_target,rho_target,Ye_target)
h_target = 1 + eps_target + press_target/rho_target

##### Conserved #####
D   = rho_target*W_target
S_1  = (rho_target*h_target)*(W_target**2)*v1_target
S_2  = (rho_target*h_target)*(W_target**2)*v2_target
S_3  = (rho_target*h_target)*(W_target**2)*v3_target
tau = (rho_target*h_target)*(W_target**2) - press_target - D
DYe = D*Ye_target
S2 = ((rho_target*h_target)**2)*(W_target**4)*(v2_target)

Cons_target = {"D"   : D,
               "S_1"  : S_1,
               "S_2"  : S_2,
               "S_3"  : S_3,
               "tau" : tau,
               "DYe" : DYe,
               "S2"  : S2}