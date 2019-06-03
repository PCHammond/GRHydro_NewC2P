# -*- coding: utf-8 -*-
"""
Created on Fri May 24 17:46:27 2019

@author: User
"""

import numpy as np
from C2P_GlobalVars import eos_table_name, eos_table_location
from C2P_GlobalVars import T_samples, rho_samples, Ye_samples
from C2P_GlobalVars import EOS_Eval_Order
from C2P_GlobalVars import global_verbose

def Setup_ReadTable(table_name,
                    table_location,
                    T_samples,
                    rho_samples,
                    Ye_samples,
                    verbose=False or global_verbose):
    if verbose:
        print("Reading EoS Table " + table_name + " from " + table_location + table_name)
    EoS_table_full_raw = np.genfromtxt(table_location + table_name).reshape((T_samples,rho_samples,Ye_samples,8))
    
    EoS_temp_raw = EoS_table_full_raw[:,0,0,0].copy()
    EoS_rho_raw  = EoS_table_full_raw[0,:,0,1].copy() * 2.7081965813424477e-03
    EoS_ye_raw   = EoS_table_full_raw[0,0,:,2].copy()

    EoS_press_raw   = EoS_table_full_raw[:,:,:,3].copy() * 2.8864091366088927e-06
    EoS_eps_raw     = EoS_table_full_raw[:,:,:,4].copy()
    #EoS_dP_dRho_raw = EoS_table_full_raw[:,:,:,5].copy() * 1.0658048815563104e-03
    #EoS_dP_deps_raw = EoS_table_full_raw[:,:,:,6].copy() * 2.7082599646973916e-03
    #EoS_cs2_raw     = EoS_table_full_raw[:,:,:,7].copy()
    
    T_array = EoS_temp_raw
    rho_array = EoS_rho_raw
    Ye_array = EoS_ye_raw
    
    pressTable = EoS_press_raw
    epsTable = EoS_eps_raw
    
    return T_array, rho_array, Ye_array, pressTable, epsTable

if EOS_Eval_Order == 3:
    def Setup_GetCoeffs(T,rho,Ye,table):
        coeffs = np.array([[[[np.inf]]]])
        return coeffs

T_array, rho_array, Ye_array, pressTable, epsTable = Setup_ReadTable(eos_table_name,eos_table_location,
                                                                     T_samples,rho_samples,Ye_samples)

log_T_array   = np.log10(T_array)
log_rho_array = np.log10(rho_array)

if EOS_Eval_Order == 1:
    press_eval_kwargs = {"table": pressTable}
    eps_eval_kwargs = {"table": epsTable}

if EOS_Eval_Order == 3:
    pressCoeffs = Setup_GetCoeffs(T_array,rho_array,Ye_array,pressTable)
    epsCoeffs = Setup_GetCoeffs(T_array,rho_array,Ye_array,epsTable)