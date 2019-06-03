# -*- coding: utf-8 -*-
"""
Created on Fri May 24 16:10:14 2019

@author: User
"""

import numpy as np
from C2P_EOS import EOS_pressFromWzTCons, EOS_epsFromWzTCons
from C2P_EOS import EOS_dPdWFromWzTCons, EOS_dEdWFromWzTCons, EOS_dPdTFromWzTCons, EOS_dEdTFromWzTCons

def Constraints_1(WzT,Cons):
    W = WzT["W"]
    z = WzT["z"]
    D = Cons["D"]
    tau = Cons["tau"]
    press = EOS_pressFromWzTCons(WzT,Cons)
    return (tau + D - z + press)*(W**2)

def Constraints_2(WzT,Cons):
    W = WzT["W"]
    z = WzT["z"]
    D = Cons["D"]
    rho = D/W
    h = z/(rho*W*W)
    v2 = 1-(1/(W**2))
    #S2 = ((rho*h)**2)*(W**4)*(v2)
    #S2 = EOS_S2FromWzTCons(WzT,Cons)
    S2 = Cons["S2"]
    press_eos = EOS_pressFromWzTCons(WzT,Cons)
    eps_eos = EOS_epsFromWzTCons(WzT,Cons)
    h_eos = 1 + eps_eos + press_eos/rho
    z_eos = rho * h_eos * (W**2)
    #print("W:",W,", z:",z,", S2:",S2)
    #return (((z**2) - S2)*(W**2)) - (z**2)
    #print((((z**2) - S2)*(W**2)),z**2,z_eos**2)
    return (((z**2) - S2)*(W**2)) - (z**2)

def Constraints_3(WzT,Cons):
    W = WzT["W"]
    z = WzT["z"]
    D = Cons["D"]
    press = EOS_pressFromWzTCons(WzT,Cons)
    eps = EOS_epsFromWzTCons(WzT,Cons)
    #print(WzT)
    #print(((z - D*W - press*(W**2))/(D*W)),eps,((z - D*W - press*(W**2))/(D*W)) - eps)
    return ((z - D*W - press*(W**2))/(D*W)) - eps

def Constraints_Jacobian(WzT,Cons,delta_WzT):
    jacobian = np.zeros((3,3))
    WzT_array = np.array([WzT["W"],WzT["z"],WzT["T"]])
    d_WzT_array = np.array([delta_WzT["W"],delta_WzT["z"],delta_WzT["T"]])
    for i in range(3):
        offset = np.zeros(3)
        offset[i] = 1.0
        WzT_plus_array = WzT_array.copy()
        WzT_minus_array = WzT_array.copy()
        WzT_plus_array[i] += d_WzT_array[i]
        WzT_minus_array[i] -= d_WzT_array[i]
        WzT_plus = {"W":WzT_plus_array[0],"z":WzT_plus_array[1],"T":WzT_plus_array[2]}
        WzT_minus = {"W":WzT_minus_array[0],"z":WzT_minus_array[1],"T":WzT_minus_array[2]}
        #print(WzT_plus)
        #print(WzT_minus)
        #First Constraint
        Constraint_plus = Constraints_1(WzT_plus,Cons)
        Constraint_minus = Constraints_1(WzT_minus,Cons)
        jacobian[0,i] = (Constraint_plus-Constraint_minus)/(2*d_WzT_array[i])
        #print(Constraint_plus)
        #print(Constraint_minus)
        #print(d_WzT_array[i])
        #print((Constraint_plus-Constraint_minus)/(2*d_WzT_array[i]))
        
        #Second Constraint
        Constraint_plus = Constraints_2(WzT_plus,Cons)
        Constraint_minus = Constraints_2(WzT_minus,Cons)
        jacobian[1,i] = (Constraint_plus-Constraint_minus)/(2*d_WzT_array[i])
        #print(Constraint_plus)
        #print(Constraint_minus)
        #print(d_WzT_array[i])
        #print((Constraint_plus-Constraint_minus)/(2*d_WzT_array[i]))
        
        #Third Constraint
        Constraint_plus = Constraints_3(WzT_plus,Cons)
        Constraint_minus = Constraints_3(WzT_minus,Cons)
        jacobian[2,i] = (Constraint_plus-Constraint_minus)/(2*d_WzT_array[i])
        #print(Constraint_plus)
        #print(Constraint_minus)
        #print(d_WzT_array[i])
        #print((Constraint_plus-Constraint_minus)/(2*d_WzT_array[i]))        

    return jacobian

def Constraints_Jacobain_SemiAnalytic(WzT,Cons,delta_WzT):
    jacobian = np.zeros((3,3))
    WzT_array = np.array([WzT["W"],WzT["z"],WzT["T"]])
    d_WzT_array = np.array([delta_WzT["W"],delta_WzT["z"],delta_WzT["T"]])
    tau = Cons["tau"]
    D = Cons["D"]
    S2 = Cons["S2"]
    W = WzT["W"]
    z = WzT["z"]
    T = WzT["T"]
    press = EOS_pressFromWzTCons(WzT,Cons)
    
    ##### Constraints_1 #####
    # df1/dw
    #WzT_p = WzT.copy(); WzT_p["W"] += delta_WzT["W"]
    #WzT_m = WzT.copy(); WzT_m["W"] -= delta_WzT["W"]
    #press_p = EOS_pressFromWzTCons(WzT_p,Cons)
    #press_m = EOS_pressFromWzTCons(WzT_m,Cons)
    #dpdW = (press_p-press_m)/(2*delta_WzT["W"])
    dpdW = EOS_dPdWFromWzTCons(WzT,Cons)
    jacobian[0,0] = 2*W*(tau+D-z+press) + (W**2)*dpdW
    # df1/dz
    jacobian[0,1] = -(WzT["W"])**2
    # df1/dT
    #WzT_p = WzT.copy(); WzT_p["T"] += delta_WzT["T"]
    #WzT_m = WzT.copy(); WzT_m["T"] -= delta_WzT["T"]
    #press_p = EOS_pressFromWzTCons(WzT_p,Cons)
    #press_m = EOS_pressFromWzTCons(WzT_m,Cons)
    #dpdT = (press_p-press_m)/(2*delta_WzT["T"])
    dpdT = EOS_dPdTFromWzTCons(WzT,Cons)
    jacobian[0,2] = (W**2)*dpdT
    
    ##### Constraints_2 #####
    # df2/dw
    jacobian[1,0] = 2*W*(z**2 - S2)
    # df2/dz
    jacobian[1,1] = 2*z*((W**2)-1)
    # df2/dT
    jacobian[1,2] = 0.0
    
    ##### Constraints_3 #####
    # df3/dw
    #WzT_p = WzT.copy(); WzT_p["W"] += delta_WzT["W"]
    #WzT_m = WzT.copy(); WzT_m["W"] -= delta_WzT["W"]
    #press_p = EOS_pressFromWzTCons(WzT_p,Cons)
    #press_m = EOS_pressFromWzTCons(WzT_m,Cons)
    #dpdW = (press_p-press_m)/(2*delta_WzT["W"])
    #eps_p = EOS_epsFromWzTCons(WzT_p,Cons)
    #eps_m = EOS_epsFromWzTCons(WzT_m,Cons)
    #dedW = (eps_p-eps_m)/(2*delta_WzT["W"])
    dpdW = EOS_dPdWFromWzTCons(WzT,Cons)
    dedW = EOS_dEdWFromWzTCons(WzT,Cons)
    jacobian[2,0] = ((-z)/(D*(W**2))) - ((press + W*dpdW)/D) - dedW
    
    # df3/dz
    jacobian[2,1] = 1/(D*W)
    
    # df3/dt
    #WzT_p = WzT.copy(); WzT_p["T"] += delta_WzT["T"]
    #WzT_m = WzT.copy(); WzT_m["T"] -= delta_WzT["T"]
    #press_p = EOS_pressFromWzTCons(WzT_p,Cons)
    #press_m = EOS_pressFromWzTCons(WzT_m,Cons)
    #dpdT = (press_p-press_m)/(2*delta_WzT["T"])
    #eps_p = EOS_epsFromWzTCons(WzT_p,Cons)
    #eps_m = EOS_epsFromWzTCons(WzT_m,Cons)
    #dedT = (eps_p-eps_m)/(2*delta_WzT["T"])
    dpdT = EOS_dPdTFromWzTCons(WzT,Cons)
    dedT = EOS_dEdTFromWzTCons(WzT,Cons)
    jacobian[2,2] = -dedT - (W/D)*dpdT

    return jacobian

def Constraints_vector(WzT,Cons):
    f_vector = np.array([Constraints_1(WzT,Cons),
                         Constraints_2(WzT,Cons),
                         Constraints_3(WzT,Cons)])
    return f_vector