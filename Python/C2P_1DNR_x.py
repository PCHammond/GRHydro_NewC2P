#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 15:57:04 2019

@author: peter
"""

import numpy as np
from C2P_EOS import EOS_pressFromPrim, EOS_epsFromPrim, EOS_dEdTFromPrim, EOS_dPdrhoFromPrim, EOS_dPdTFromPrim
from C2P_Setup import log_T_array, log_rho_array, Ye_array
from C2P_InvertT import TInversionDriver as TInversionDriver
from C2P_InitialGuess import GetInitialGuess

def RelDiff(a,b):
    if a==b:
        return 0.0
    return np.abs(2*(a-b)/(a+b))

def Root(x,Cons):
    D = Cons["D"]
    tau = Cons["tau"]
    S2 = Cons["S2"]
    DYe = Cons["DYe"]
    
    Ye = DYe/D
    
    q = tau/D
    r = S2/(D**2)

    W_guess_m2 = 1 - r/(x**2)
    W_guess = np.sqrt(1/W_guess_m2)
    rho_guess = D/W_guess
    if rho_guess<(10**log_rho_array[0]):
        W_guess = D/(10**log_rho_array[0])
        W_guess_m2 = 1/(W_guess**2)
    eps_guess = -1 + (x*(1-(W_guess**2))/W_guess) + W_guess*(1+q)
    T_guess, its = TInversionDriver(eps_guess,rho_guess,Ye)
    press_guess = EOS_pressFromPrim(T_guess,rho_guess,Ye)
    h_guess = 1 + eps_guess + press_guess/rho_guess
    return x - h_guess*W_guess

def drdx_1DNR(x,r,D,W):
    dWdx = dWdx_1DNR(x,r)
    return -D*dWdx/(W**2)

def dWdx_1DNR(x,r):
    return -r*(x**(-3))*((1-r*(x**(-2)))**(-3/2))

def dedx_1DNR(x,r,W,q):
    dWdx = dWdx_1DNR(x,r)
    p1 = ((1/W) - x*dWdx/(W**2))*(1-(W**2))
    p2 = -2*x + (1+q)*dWdx
    return p1 - p2

def dTdx_1DNR(x,r,D,W,q,T,rho,Ye,eps):
    rho_plu = rho*(1.0+1e-6)
    rho_min = rho*(1.0-1e-6)
    T_plu, its = TInversionDriver(eps,rho_plu,Ye)
    T_min, its = TInversionDriver(eps,rho_min,Ye)
    dTdr = (T_plu-T_min)/(rho_plu-rho_min)
    drdx = drdx_1DNR(x,r,D,W)
    #delta_eps = np.max([np.abs(1e-6*eps),1e-9])
    #eps_plu = eps + delta_eps
    #eps_min = eps - delta_eps
    #T_plu = TInversionDriver(eps_plu,rho,Ye)
    #T_min = TInversionDriver(eps_min,rho,Ye)
    #dTde = (T_plu-T_min)/(eps_plu-eps_min)
    dedT = EOS_dEdTFromPrim(T,rho,Ye)
    dTde = 1/dedT
    dedx = dedx_1DNR(x,r,W,q)
    return dTdr*drdx + dTde*dedx

def dpdx_1DNR(x,r,q,D,W,T,rho,Ye,eps):
    dpdr = EOS_dPdrhoFromPrim(T,rho,Ye)
    drdx = drdx_1DNR(x,r,D,W)
    dpdT = EOS_dPdTFromPrim(T,rho,Ye)
    dTdx = dTdx_1DNR(x,r,D,W,q,T,rho,Ye,eps)
    return dpdr*drdx + dpdT*dTdx

def dhdx_1DNR(x,r,q,D,W,T,rho,Ye,press,eps):
    dedx = dedx_1DNR(x,r,W,q)
    drdx = drdx_1DNR(x,r,D,W)
    dpdx = dpdx_1DNR(x,r,q,D,W,T,rho,Ye,eps)
    return dedx + (rho*dpdx - press*drdx)/(rho**2)

def dRootdx_1DNR(x,Cons):
    D = Cons["D"]
    tau = Cons["tau"]
    S2 = Cons["S2"]
    DYe = Cons["DYe"]
    
    Ye = DYe/D
    
    q = tau/D
    r = S2/(D**2)

    W_guess_m2 = 1 - r/(x**2)
    W_guess = np.sqrt(1/W_guess_m2)
    rho_guess = D/W_guess
    if rho_guess<(10**log_rho_array[0]):
        W_guess = D/(10**log_rho_array[0])
        W_guess_m2 = 1/(W_guess**2)
    eps_guess = -1 + (x*(1-(W_guess**2))/W_guess) + W_guess*(1+q)
    T_guess, its = TInversionDriver(eps_guess,rho_guess,Ye)
    press_guess = EOS_pressFromPrim(T_guess,rho_guess,Ye)
    h_guess = 1 + eps_guess + press_guess/rho_guess

    dhdx = dhdx_1DNR(x,r,q,D,W_guess,T_guess,rho_guess,Ye,press_guess,eps_guess)
    dWdx = dWdx_1DNR(x,r)
    return 1 - W_guess*dhdx - h_guess*dWdx

def ConsX2Prim(Cons,x):
    D = Cons["D"]
    tau = Cons["tau"]
    S2 = Cons["S2"]
    DYe = Cons["DYe"]
    S_1 = Cons["S_1"]
    S_2 = Cons["S_2"]
    S_3 = Cons["S_3"]
    
    Ye = DYe/D
    q = tau/D
    r = S2/(D**2)
    
    W_m2 = 1 - r/(x**2)
    W = np.sqrt(1/W_m2)
    
    rho = D/W
    eps = -1 + (x*(1-(W**2))/W) + W*(1+q)
    T, its = TInversionDriver(eps,rho,Ye)
    
    z = x*rho*W
    
    v1 = S_1/z
    v2 = S_2/z
    v3 = S_3/z
    
    return np.array([T, rho, Ye, v1, v2, v3])

def Con2Prim(Cons,x_guess=None,rel_diff=1e-12):
    D = Cons["D"]
    tau = Cons["tau"]
    
    q = tau/D
    
    if x_guess==None:
        x_current = 1.5*(1.0+q)
    else:
        x_current = x_guess
    
    done=False
    max_its = 500
    it = 0
    
    while done==False:
        root = Root(x_current,Cons)
        drootdx = dRootdx_1DNR(x_current,Cons)
        delta_x = root/drootdx
        
        x_current = x_current + delta_x

        if x_current<(1+q):
            x_current=1+q
            delta_x=1+q
        elif x_current>2*(1+q):
            x_current=2*(1+q)
            delta_x=2*(1+q)
        
        it += 1
        if np.abs(delta_x/x_current)<rel_diff:
            done=True
        
        if it>=max_its:
            done=True
            it=-1
    
    return x_current, it

def Con2Prim_Dekker(Cons,rel_diff=1e-15,max_it=100):    
    D = Cons["D"]
    tau = Cons["tau"]
    S2 = Cons["S2"]
    DYe = Cons["DYe"]
    
    Ye = DYe/D
    
    q = tau/D
    r = S2/(D**2)

    a0 =    1.0 + q
    b0 = 2*(1.0 + q)
    
    epsilon_0 = np.abs(a0-b0)
    epsilon_target = 0.5*(a0+b0)*rel_diff
    n_its = int((np.log(epsilon_0) - np.log(epsilon_target))/np.log(2))
    
    fa0 = Root(a0,Cons)
    fb0 = Root(b0,Cons)
    
    if np.abs(fa0)<np.abs(fb0):
        a0,  b0  = b0,  a0
        fa0, fb0 = fb0, fa0
    
    ak = a0    
    
    bk   = b0
    bkm1 = a0
    fbkm1 = Root(bkm1,Cons)

    fak = Root(ak,Cons)
    fbk = Root(bk,Cons)
    
    it = 0
    max_its = int(np.max([max_it,n_its]))
    done = False
    
    if fak==0:
        bk = ak
        done = True
    elif fbk==0:
        done = True
    
    while done == False:        
        m = (ak + bk)/2.0
        
        if fbk!=fbkm1 and 2*np.abs((bk-bkm1)/(bk+bkm1))>1e-9:
            s = bk - fbk*(bk-bkm1)/(fbk-fbkm1)
        else:
            s = m
        
        if (s-bk)*(bk-m)<0:
            bkp1 = s
        else:
            bkp1 = m
            
        fbkp1 = Root(bkp1,Cons)
        
        if fbkp1*fbk<0:
            akp1 = bk
            fakp1 = fbk
        else:
            akp1 = ak
            fakp1 = fak
            
        if np.abs(fakp1)<np.abs(fbkp1):
            akp1,  bkp1  = bkp1,  akp1
            fakp1, fbkp1 = fbkp1, fakp1
        
        ### update ###
        bkm1  = bk
        fbkm1 = fbk
        
        ak  = akp1
        fak = fakp1
        
        bk  = bkp1
        fbk = fbkp1
        
        #if np.abs(2*(bk-bkm1)/(bk+bkm1))<rel_diff:
        #    done = True
        
        #
        if np.abs(fbk/bk)<rel_diff:
            done = True
            
        it += 1
        if it>=max_its:
            done = True
            it = -1
        
    return bk, ak, bk, it

def Con2Prim_Bisection(Cons,a0=None,b0=None,rel_diff=1e-15,max_it=100):    
    D = Cons["D"]
    tau = Cons["tau"]
    S2 = Cons["S2"]
    DYe = Cons["DYe"]
    
    Ye = DYe/D
    
    q = tau/D
    r = S2/(D**2)
    
    if a0==None:
        a0 =    1.0 + q
    if b0==None:
        b0 = 2*(1.0 + q)
    
    epsilon_0 = np.abs(a0-b0)
    epsilon_target = 0.5*(a0+b0)*rel_diff
    n_its = int((np.log(epsilon_0) - np.log(epsilon_target))/np.log(2))
    
    fa0 = Root(a0,Cons)
    fb0 = Root(b0,Cons)
    
    if fa0*fb0>0.0:
        raise ValueError("fa0, fb0 have same sign!")
    
    ak = a0
    bk = b0
    
    fak = fa0
    fbk = fb0
    
    it = 0
    max_its = int(np.max([max_it,n_its]))
    done = False
    
    if fak==0:
        result = ak
        done = True
    elif fbk==0:
        result = bk
        done = True
    
    while done == False:        
        m = (ak + bk)/2.0
        
        fm = Root(m,Cons)
        
        ### update ###
        if fm==0:
            result = m
            done = True
        elif fm*fak<0:
            bk = m
            fbk = fm
        else:
            ak = m
            fak = fm
        
        if np.abs(fak/ak)<rel_diff:
            result = ak
            done = True
        
        if np.abs(fbk/bk)<rel_diff:
            result = bk
            done = True
            
        it += 1
        if it>=max_its:
            if np.abs(fak)<np.abs(fbk):
                result = ak
            else:
                result = bk
            done = True
            it = -1
        
    return result, ak, bk, it


def table_string(value):
    raw_string = "{:.5E}".format(value)
    if not raw_string[0] == "-":
        raw_string = "+" + raw_string
    return raw_string + "  "

def GetxGuess(targets):
    rho_target = targets.rho
    vsq_target = targets.vsq
    T_target = targets.T
    W = 1/np.sqrt(1-vsq)
    return None

def ConvergenceTest(Cons_target,Cons_found,printTable=True):
    cons_errors = np.array([RelDiff(Cons_target["D"],Cons_found["D"]),
                            RelDiff(Cons_target["DYe"],Cons_found["DYe"]),
                            RelDiff(Cons_target["tau"],Cons_found["tau"]),
                            RelDiff(Cons_target["S2"],Cons_found["S2"])])
    
    con_err = np.max(np.abs(cons_errors))
    
    if printTable:
        print("Conservative Vars")
        print("var      target        found         relative_error")
        print("D:      " + table_string(Cons_target["D"]) + table_string(Cons_found["D"]) + table_string(cons_errors[0]))
        print("DYe:    " + table_string(Cons_target["DYe"]) + table_string(Cons_found["DYe"]) + table_string(cons_errors[1]))
        print("tau:    " + table_string(Cons_target["tau"]) + table_string(Cons_found["tau"]) + table_string(cons_errors[2]))
        print("S2:     " + table_string(Cons_target["S2"]) + table_string(Cons_found["S2"]) + table_string(cons_errors[3]))
    
    return con_err

def MakeGridFig(iterationsForConvergence,failed,savename,vel_idx,fig_title,it_min=None,it_max=None):
    fig = plt.figure(figsize=(16,16))
    subplots = fig.subplots(nrows=4,ncols=4,sharex=True,sharey=True)
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    ax_list = subplots.copy()
    ax_list = ax_list.flatten()

    rho_min = np.log10(test_rho_array[0])
    rho_max = np.log10(test_rho_array[-1])
    delta_rho = np.log10(test_rho_array[1]) - np.log10(test_rho_array[0])

    Ye_min = test_Ye_array[0]
    Ye_max = test_Ye_array[-1]
    delta_Ye = test_Ye_array[1] - test_Ye_array[0]

    imshow_extent = np.array([rho_min - 0.5*delta_rho,
                              rho_max + 0.5*delta_rho,
                              Ye_min  - 0.5*delta_Ye,
                              Ye_max  + 0.5*delta_Ye])
    
    if it_min==None:
        it_min = 1
    if it_max==None:
        it_max = np.max(iterationsForConvergence)
    
    label_x = imshow_extent[0] + 0.1*(imshow_extent[1]-imshow_extent[0])
    label_y = imshow_extent[2] + 0.9*(imshow_extent[3]-imshow_extent[2])

    for T_idx in range(len(test_T_array)):
        ax = ax_list[T_idx]
        ax.text(label_x,
                label_y,
                "T = " + "{:.3E}".format(test_T_array[T_idx]),
                fontdict={"size":12,
                          "color":"w"},
                bbox=dict(boxstyle="square",
                   ec=(1., 1., 1.),
                   fc=(0., 0., 0.),
                   ))
        image = ax.imshow(iterationsForConvergence[T_idx,:,:].T,
                          cmap="inferno",
                          origin="lower",
                          extent=imshow_extent,
                          aspect="auto",
                          vmin=it_min,
                          vmax=it_max)

    for fail in failed:
        if fail[-1]==vel_idx:
            T_fail_idx = fail[0]
            ax = ax_list[T_fail_idx]
            ax.scatter(np.log10(test_rho_array[fail[1]]),test_Ye_array[fail[2]],s=40,marker="x",c="w")
    
    for ax in subplots[:,0]:
        ax.set_ylabel("Ye")

    for ax in subplots[-1,:]:
        ax.set_xlabel("log10(rho)")
    
    fig.tight_layout()

    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.925, 0.1, 0.025, 0.8]) #left bottom width height
    cbar = fig.colorbar(image, cax=cbar_ax)
    cbar.set_label("log10(Iterations)")
    


    fig.subplots_adjust(top=0.95)
    fig.suptitle(fig_title)
    fig.savefig(savename,dpi=480)
    
    return fig, subplots

if __name__=="__main__":
    import datetime
    print("Starting at " + str(datetime.datetime.now()))
#    rho_target = 2.8358299540575783e-05
    v1_target  = 0.1
    v2_target  = -0.2
    v3_target  = 0.15
#    Ye_target  = 0.315
#    T_target   = 4.166582513489608
#    targets, Cons_target = GetInitialGuess(rho=rho_target,Ye=Ye_target,T=T_target,velocity=np.array([v1_target,v2_target,v3_target])) 
#    #x_conv, its_taken = Con2Prim(Cons_target,x_guess=None)
#    x_conv, its_taken = Con2Prim_Dekker(Cons_target)
#    result = ConsX2Prim(Cons_target,x_conv)
#    gubs, Cons_result = GetInitialGuess(rho=result[1],velocity=result[3:],T=result[0],Ye=result[2])
#    #print(targets.T, result[0])
#    #print(targets.rho, result[1])
#    #print(targets.Ye, result[2])
#    #print(targets.v1, result[3])
#    #print(targets.v2, result[4])
#    #print(targets.v3, result[5])
#    
#    Cons_err = ConvergenceTest(Cons_target,Cons_result,printTable=True)

    T_samples = 16
    rho_samples = 32
    Ye_samples = 32
    velocity_samples = 9
    
    test_T_array = 10**(np.linspace(log_T_array[0],log_T_array[-1],num=T_samples+2)[1:-1])
    test_rho_array = 10**(np.linspace((log_rho_array[log_rho_array>-12])[0],(log_rho_array[log_rho_array<-2.5])[-1],num=rho_samples+2)[1:-1])
    test_Ye_array = np.linspace(Ye_array[0],Ye_array[-1],num=Ye_samples+2)[1:-1]
    
    velocity_test_master = np.array([v1_target,v2_target,v3_target])
    velocity_fracs = np.linspace(0.0,2.0,num=velocity_samples)
    
    iterationsForConvergence_full = np.zeros((len(test_T_array),len(test_rho_array),len(test_Ye_array),velocity_samples))
    consErrors_full = np.zeros((len(test_T_array),len(test_rho_array),len(test_Ye_array),velocity_samples))
    failed = []
    
    used_bisection = np.zeros((len(test_T_array),len(test_rho_array),len(test_Ye_array),velocity_samples),dtype=bool)
    
    start = datetime.datetime.now()
    for T_test_idx in range(len(test_T_array)):
        print(T_test_idx)
        T_test = test_T_array[T_test_idx]
        for rho_test_idx in range(len(test_rho_array)):
            rho_test = test_rho_array[rho_test_idx]
            for Ye_test_idx in range(len(test_Ye_array)):
                Ye_test = test_Ye_array[Ye_test_idx]
                for vel_test_idx in range(velocity_samples):
                    velocity_test = velocity_fracs[vel_test_idx]*velocity_test_master.copy()
                    targets, Cons_target = GetInitialGuess(rho=rho_test,Ye=Ye_test,T=T_test,velocity=velocity_test) 
                    #x_guess = GetxGuess()
                    #x_conv, its_taken = Con2Prim(Cons_target,x_guess=None)
                    x_conv, ak, bk, its_taken = Con2Prim_Dekker(Cons_target,rel_diff=1e-15,max_it=100)
                    if its_taken==-1:
                        used_bisection[T_test_idx,rho_test_idx,Ye_test_idx,vel_test_idx] = True
                        x_conv, ak, bk, its_taken = Con2Prim_Bisection(Cons_target,a0=ak,b0=bk,rel_diff=1e-15,max_it=100)
                        if its_taken!=-1:
                            its_taken+=100
                    result = ConsX2Prim(Cons_target,x_conv)
                    gubs, Cons_result = GetInitialGuess(rho=result[1],velocity=result[3:],T=result[0],Ye=result[2])
                    Cons_err = ConvergenceTest(Cons_target,Cons_result,printTable=False)
                    consErrors_full[T_test_idx,rho_test_idx,Ye_test_idx,vel_test_idx] = Cons_err
                    if its_taken == -1:
                        failed.append(np.array([T_test_idx,rho_test_idx,Ye_test_idx,vel_test_idx]))
                        its_taken = 1
                    iterationsForConvergence_full[T_test_idx,rho_test_idx,Ye_test_idx,vel_test_idx] = its_taken    
    
    end = datetime.datetime.now()
    
    print("Finishing at " + str(datetime.datetime.now()))
    
    print(str((end-start)/(T_samples*rho_samples*Ye_samples*velocity_samples)) + " per sample")
    print(str(len(failed))+" failed")
    print(str(np.sum(used_bisection)) + " used bisection")
    print(np.max(iterationsForConvergence_full),np.average(iterationsForConvergence_full))
    print(np.max(consErrors_full),np.average(consErrors_full))
    
    iterationsForConvergence = np.max(iterationsForConvergence_full,axis=-1)
    consErrors = np.max(consErrors_full,axis=-1)
    
    it_min = 0
    it_max = np.log10(np.max(iterationsForConvergence_full))
    iterationsForConvergence_full = np.log10(iterationsForConvergence_full)
    
    plot = True
    if plot == True:
        import matplotlib as mpl
        mpl.use("Agg")
        mpl.rcParams['mathtext.fontset'] = 'stix'
        mpl.rcParams['font.family'] = 'STIXGeneral'

        
        import matplotlib.pyplot as plt
        plt.ioff()
        
        for vel_idx in range(velocity_samples):
            title = r"1D Newton Raphson Con2Prim, $|v| = " + str(np.linalg.norm(velocity_fracs[vel_idx]*velocity_test_master.copy())) + r"$"
            MakeGridFig(iterationsForConvergence_full[:,:,:,vel_idx],failed,"DekkerConvergencev" + str(vel_idx) + ".png",vel_idx,title,it_min=it_min,it_max=it_max)