#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 11:48:26 2019

@author: peter
"""

import numpy as np
from C2P_EOS import EOS_pressFromPrim, EOS_epsFromPrim, EOS_dEdTFromPrim
from C2P_Setup import log_T_array, log_rho_array, Ye_array

def RelDiff(a,b):
    return np.abs(2*(a-b)/(a+b))

def NR_TFromEpsRhoYe(eps_target,rho,Ye,max_its=100,rel_diff=1e-14,abs_diff=1e-15,T_guess=None,NR_fac=1.0):
    if T_guess == None:
        T_guess = 1.0
    it=0
    done=False
    T_current = T_guess
    eps_current = EOS_pressFromPrim(T_current,rho,Ye)
    
    T_record = []
    T_record.append(T_current)
    
    success = False
    
    while done == False:
        depsdT = -EOS_dEdTFromPrim(T_current,rho,Ye)
        root = eps_target-eps_current
        delta_T = root/depsdT
        T_current = T_current - NR_fac*delta_T
        if T_current<10**log_T_array[0]:
            T_current = 10**log_T_array[0]
        elif T_current>10**log_T_array[-1]:
            T_current = 10**log_T_array[-1]
            
        eps_current = EOS_epsFromPrim(T_current,rho,Ye)
        
        T_record.append(T_current)
        
        it+=1
        #if it>20 and it%20==0:
            #T_current = np.average(np.array(T_record)[it-10:it])
        if it>=max_its:
            done = True
            success = False
            T_current = np.average(np.array(T_record)[-20:])
        if np.abs(eps_target-eps_current)<abs_diff:
            done = True
            success = True
            #print("Absolute convergence in " + str(it) + " iterations")
        elif RelDiff(eps_target,eps_current)<rel_diff:
            done = True
            success = True
            #print("Relative convergence in " + str(it) + " iterations")
    return T_current, success, it

def Root(eps,T,rho,Ye):
    return eps - EOS_epsFromPrim(T,rho,Ye)

def Bisection_TFromEpsRho(eps_target,rho,Ye,max_it=100,rel_diff=1e-14,abs_diff=1e-15,T_guess_low=None,T_guess_high=None):
    if T_guess_low==None:
        T_guess_low = 10**log_T_array[0]
    if T_guess_high==None:
        T_guess_high = 10**log_T_array[-1]
        
    a0 = T_guess_low
    b0 = T_guess_high
    
    fa0 = Root(eps_target,a0,rho,Ye)
    fb0 = Root(eps_target,b0,rho,Ye)
        
    epsilon_0 = np.abs(a0-b0)
    epsilon_target = 0.5*(a0+b0)*rel_diff
    n_its = int((np.log(epsilon_0) - np.log(epsilon_target))/np.log(2))
    
    it = 0
    max_its = int(np.max([max_it,n_its]))
    done = False
    
    if fa0==0:
        T_result = a0
        success=True
        done = True
    elif fb0==0:
        T_result = b0
        success=True
        done = True
    
    eps_a0 = EOS_epsFromPrim(a0,rho,Ye)
    eps_b0 = EOS_epsFromPrim(b0,rho,Ye)
    
    if eps_target<eps_a0:
        T_result = a0
        success = False
        done = True
    
    
    if eps_target>eps_b0:
        T_result = b0
        success = False
        done = True
    
    if fa0*fb0>0.0 and done==False:
        print("eps_target",eps_target)
        print("rho",rho)
        print("Ye",Ye)
        print("a0",a0)
        print("fa0",fa0)
        print("eps(a0)",EOS_epsFromPrim(a0,rho,Ye))
        print("b0",b0)
        print("fb0",fb0)
        print("eps(b0)",EOS_epsFromPrim(b0,rho,Ye))
        raise ValueError("fa0, fb0 have same sign!")

    
    ak = a0
    bk = b0
    
    fak = fa0
    fbk = fb0
    
    while done==False:
        m = (ak + bk)/2.0
        
        fm = Root(eps_target,m,rho,Ye)
        
        ### update ###
        if fm==0:
            T_result = m
            success=True
            done = True
        elif fm*fak<0:
            bk = m
            fbk = fm
        else:
            ak = m
            fak = fm
        
        if np.abs(fak/ak)<rel_diff:
            T_result = ak
            success=True
            done = True
        
        if np.abs(fbk/bk)<rel_diff:
            T_result = bk
            success=True
            done = True
        
        it += 1
        if it>=max_its:
            if np.abs(fak)<np.abs(fbk):
                T_result = ak
            else:
                T_result = bk
            done = True
            success=False
    
    return T_result, success, it

def TInversionDriver(eps_target,rho,Ye,max_its=100,rel_diff=1e-14,abs_diff=1e-15,T_guess=None):
    its_taken = 0
    T_conv, success, its = NR_TFromEpsRhoYe(eps_target,rho,Ye,max_its=max_its,rel_diff=rel_diff,abs_diff=abs_diff,T_guess=T_guess,NR_fac=1.0)
    its_taken += its
    if success==True:
        return T_conv, its_taken
    else:
        #print("First attempt failed, using smaller NR_fac")
        T_conv, success, its = NR_TFromEpsRhoYe(eps_target,rho,Ye,max_its=2*max_its,rel_diff=1e-14,abs_diff=1e-15,T_guess=T_conv,NR_fac=0.5)
        its_taken += its    
    if success==True:
        return T_conv, its_taken
    else:
        #print("Second attempt failed, using smaller NR_fac")
        T_conv, success, its = NR_TFromEpsRhoYe(eps_target,rho,Ye,max_its=10*max_its,rel_diff=1e-14,abs_diff=1e-15,T_guess=T_conv,NR_fac=0.1)
        its_taken += its
    if success==True:
        return T_conv, its_taken
    else:
        #raise ValueError("T inversion failed")
        return T_conv, -1
    return None

def TInversionDriverBis(eps_target,rho,Ye,max_its=100,rel_diff=1e-15,abs_diff=1e-15,T_guess_low=None,T_guess_high=None):
    T_conv, success, its = Bisection_TFromEpsRho(eps_target,rho,Ye,max_it=max_its,rel_diff=rel_diff,abs_diff=abs_diff,T_guess_low=T_guess_low,T_guess_high=T_guess_high)
    return T_conv, its

def MakeGridFig(iterationsForConvergence,failed,savename,fig_title,it_min=None,it_max=None):
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
    
    failed = []
    
    T_samples = 16
    rho_samples = 64
    Ye_samples = 64

    test_T_array = 10**(np.linspace(log_T_array[0],log_T_array[-1],num=T_samples+2)[1:-1])
    test_rho_array = 10**(np.linspace((log_rho_array[log_rho_array>-12])[0],(log_rho_array[log_rho_array<-2.5])[-1],num=rho_samples+2)[1:-1])
    test_Ye_array = np.linspace(Ye_array[0],Ye_array[-1],num=Ye_samples+2)[1:-1]
    
    iterationsForConvergence = np.zeros((len(test_T_array),len(test_rho_array),len(test_Ye_array)))
    errors = np.zeros((len(test_T_array),len(test_rho_array),len(test_Ye_array)))

    for T_test_idx in range(len(test_T_array)):
        print(T_test_idx)
        T_test = test_T_array[T_test_idx]
        for rho_test_idx in range(len(test_rho_array)):
            rho_test = test_rho_array[rho_test_idx]
            for Ye_test_idx in range(len(test_Ye_array)):
                Ye_test = test_Ye_array[Ye_test_idx]
                eps_test = EOS_epsFromPrim(T_test,rho_test,Ye_test)
                T_conv, its_taken = TInversionDriverBis(eps_test,rho_test,Ye_test)            
                errors[T_test_idx,rho_test_idx,Ye_test_idx] = 0.5*(T_test-T_conv)/(T_test+T_conv)
                if its_taken == -1:
                    failed.append(np.array([T_test_idx,rho_test_idx,Ye_test_idx]))
                    its_taken = 1
                iterationsForConvergence[T_test_idx,rho_test_idx,Ye_test_idx] = its_taken
    
    print(np.max(errors))
    
    iterationsForConvergence = np.log10(iterationsForConvergence)
    it_min = 0
    it_max = np.max(iterationsForConvergence)
    
    
    plot = True
    if plot == True:
        import matplotlib as mpl
        mpl.use("Agg")
        mpl.rcParams['mathtext.fontset'] = 'stix'
        mpl.rcParams['font.family'] = 'STIXGeneral'

        
        import matplotlib.pyplot as plt
        plt.ioff()
        
        title = r"Bisection T Inversion"
        MakeGridFig(iterationsForConvergence,failed,"BisTConvergence.png",title,it_min=it_min,it_max=it_max)
    
#    import matplotlib.pyplot as plt
#
#    fig = plt.figure(figsize=(16,16))
#    subplots = fig.subplots(nrows=4,ncols=4,sharex=True,sharey=True)
#    fig.subplots_adjust(hspace=0)
#    fig.subplots_adjust(wspace=0)
#    ax_list = subplots.copy()
#    ax_list = ax_list.flatten()
#
#    rho_min = np.log10(test_rho_array[0])
#    rho_max = np.log10(test_rho_array[-1])
#    delta_rho = np.log10(test_rho_array[1]) - np.log10(test_rho_array[0])
#
#    Ye_min = test_Ye_array[0]
#    Ye_max = test_Ye_array[-1]
#    delta_Ye = test_Ye_array[1] - test_Ye_array[0]
#
#    imshow_extent = np.array([rho_min - 0.5*delta_rho,
#                              rho_max + 0.5*delta_rho,
#                              Ye_min  - 0.5*delta_Ye,
#                              Ye_max  + 0.5*delta_Ye])
#
#    it_max = np.max(iterationsForConvergence)
#    it_min = 0
#
#    label_x = imshow_extent[0] + 0.1*(imshow_extent[1]-imshow_extent[0])
#    label_y = imshow_extent[2] + 0.9*(imshow_extent[3]-imshow_extent[2])
#    
#    for T_idx in range(len(test_T_array)):
#        ax = ax_list[T_idx]
#        ax.text(label_x,label_y,"T = " + "{:.3E}".format(test_T_array[T_idx]),fontdict={"size":12,"color":"w"})
#        image = ax.imshow(iterationsForConvergence[T_idx,:,:].T,
#                          cmap="inferno",
#                          origin="lower",
#                          extent=imshow_extent,
#                          aspect="auto",
#                          vmin=it_min,
#                          vmax=it_max)
#
#    for fail in failed:
#        T_fail_idx = fail[0]
#        ax = ax_list[T_fail_idx]
#        ax.scatter(np.log10(test_rho_array[fail[1]]),test_Ye_array[fail[2]],s=40,marker="x",c="w")
#        
#    for ax in subplots[:,0]:
#        ax.set_ylabel("Ye")
#
#    for ax in subplots[-1,:]:
#        ax.set_xlabel("log10(rho)")
#
#    fig.tight_layout()
#    fig.savefig("test_T_NR.png",dpi=480)