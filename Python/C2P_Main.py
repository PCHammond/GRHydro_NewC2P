# -*- coding: utf-8 -*-
"""
Created on Fri May 24 16:09:15 2019

@author: User
"""

import numpy as np
from C2P_InitialGuess import GetInitialGuess, GuessFromTargets
from C2P_EOS import EOS_pressFromPrim, EOS_epsFromPrim
from C2P_Setup import log_T_array, log_rho_array, Ye_array
from C2P_3DNR import Con2Prim

T_min = 10**log_T_array[1]
T_max = 10**log_T_array[-2]

def table_string(value):
    raw_string = "{:.5E}".format(value)
    if not raw_string[0] == "-":
        raw_string = "+" + raw_string
    return raw_string + "  "

def RelDiff(a,b):
    return np.abs(2*(a-b)/(a+b))

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

def ConvergenceTest(targets,WzT,Cons_target,printTable=True):
    W_target = targets.W
    T_target = targets.T
    h_target = targets.h
    rho_target = targets.rho
    press_target = targets.press
    eps_target = targets.eps
    v1_target = targets.v1
    v2_target = targets.v2
    v3_target = targets.v3
    
    Ye_found = Cons_target["DYe"]/Cons_target["D"]
    rho_found = Cons_target["D"]/WzT["W"]
    T_found = WzT["T"]
    v1_found = Cons_target["S_1"]/WzT["z"]
    v2_found = Cons_target["S_2"]/WzT["z"]
    v3_found = Cons_target["S_3"]/WzT["z"]
    press_found = EOS_pressFromPrim(T_found,rho_found,Ye_found)
    eps_found = EOS_epsFromPrim(T_found,rho_found,Ye_found)
    
    gubs, Cons_found = GetInitialGuess(rho=rho_found,velocity=np.array([v1_found,v2_found,v3_found]),T=T_found,Ye=Ye_found)
    cons_errors = np.array([RelDiff(Cons_target["D"],Cons_found["D"]),
                            RelDiff(Cons_target["DYe"],Cons_found["DYe"]),
                            RelDiff(Cons_target["tau"],Cons_found["tau"]),
                            RelDiff(Cons_target["S2"],Cons_found["S2"])])
    
    errors = np.array([2*(rho_found-rho_target)/(rho_found+rho_target),
                       2*(T_found-T_target)/(T_found+T_target),
                       2*(v1_found-v1_target)/(v1_found+v1_target),
                       2*(v2_found-v2_target)/(v2_found+v2_target),
                       2*(v3_found-v3_target)/(v3_found+v3_target),
                       2*(press_found-press_target)/(press_found+press_target),
                       2*(eps_found-eps_target)/(eps_found+eps_target)])
    
    prim_err = np.sqrt(np.sum(errors**2)/len(errors))
    con_err = np.max(np.abs(cons_errors))
    
    if printTable:
        print("Primitive Vars")
        print("var      target        found         relative_error")
        print("rho:    " + table_string(rho_target)   + table_string(rho_found)   + table_string(np.abs(2*(rho_found-rho_target)/(rho_found+rho_target))))
        print("T:      " + table_string(T_target)     + table_string(T_found)     + table_string(np.abs(2*(T_found-T_target)/(T_found+T_target))))
        print("v1:     " + table_string(v1_target)    + table_string(v1_found)    + table_string(np.abs(2*(v1_found-v1_target)/(v1_found+v1_target))))
        print("v2:     " + table_string(v2_target)    + table_string(v2_found)    + table_string(np.abs(2*(v2_found-v2_target)/(v2_found+v2_target))))
        print("v3:     " + table_string(v3_target)    + table_string(v3_found)    + table_string(np.abs(2*(v3_found-v3_target)/(v3_found+v3_target))))
        print("press:  " + table_string(press_target) + table_string(press_found) + table_string(np.abs(2*(press_found-press_target)/(press_found+press_target))))
        print("eps:    " + table_string(eps_target)   + table_string(eps_found)   + table_string(np.abs(2*(eps_found-eps_target)/(eps_found+eps_target))))
        print("")
        print("Conservative Vars")
        print("var      target        found         relative_error")
        print("D:      " + table_string(Cons_target["D"]) + table_string(Cons_found["D"]) + table_string(cons_errors[0]))
        print("DYe:    " + table_string(Cons_target["DYe"]) + table_string(Cons_found["DYe"]) + table_string(cons_errors[1]))
        print("tau:    " + table_string(Cons_target["tau"]) + table_string(Cons_found["tau"]) + table_string(cons_errors[2]))
        print("S2:     " + table_string(Cons_target["S2"]) + table_string(Cons_found["S2"]) + table_string(cons_errors[3]))
    
    return prim_err, con_err

samples = 10000
errors = np.zeros(samples)
con_errors = np.zeros(samples)
rhos = np.zeros(samples)
Ts = np.zeros(samples)
Yes = np.zeros(samples)
presss = np.zeros(samples)
epss = np.zeros(samples)
speeds = np.zeros(samples)
iterations = np.zeros(samples)
guesses_array = np.zeros((samples,3))
targets_array = np.zeros((samples,3))
results_array = np.zeros((samples,3))
guess_press = np.zeros(samples) 
err_max = 0.0
max_it = 0

print("Sample  W             z             T             average error")
for i in range(samples):
    targets, Cons_target = GetInitialGuess()
    #min_fracs = np.array([0.95,0.95,0.95]) #rho, vsq, T
    min_fracs = 0.90*np.ones(3)
    #min_fracs[0] = 0.98
    max_fracs = 1/min_fracs
    guesses = GuessFromTargets(targets,min_frac=min_fracs,max_frac=max_fracs)
    guesses_array[i] = np.array([guesses.W,guesses.z,guesses.T])
    guess_press[i] = guesses.press
    targets_array[i] = np.array([targets.W,targets.z,targets.T])
    WzT_found, its_taken, W_record, z_record, T_record = Con2Prim(Cons_target,guesses=guesses,its_max=100,press_rel_diff=1e-14,var_rel_diff=1e-12)
    results_array[i] = np.array([WzT_found["W"],WzT_found["z"],WzT_found["T"]])
    iterations[i] = its_taken
    error, con_err = ConvergenceTest(targets,WzT_found,Cons_target,printTable=False)
    con_errors[i] = con_err
    if error>=np.max(errors):
        WzT_max_err = WzT_found.copy()
        Cons_max_err = Cons_target.copy()
        targets_max_err = targets.copy()
    if error > err_max:
        err_max = error
        W_record_maxerr = W_record
        z_record_maxerr = z_record
        T_record_maxerr = T_record
        
    if its_taken > max_it:
        max_it = its_taken
        W_record_maxit = W_record
        z_record_maxit = z_record
        T_record_maxit = T_record
        
    errors[i] = error
    rhos[i] = targets.rho
    Ts[i] = targets.T
    Yes[i] = targets.Ye
    presss[i] = targets.press
    epss[i] = targets.eps
    speeds[i] = np.sqrt(targets.vsq)
    WzT_array = np.array(list(WzT_found.values()))
    if (i+1)%100==0: 
        print("{:06d}  ".format(i+1) + table_string(WzT_array[0]) + table_string(WzT_array[1]) + table_string(WzT_array[2]) + table_string(error))

errors = np.clip(errors,a_min=np.finfo(float).eps,a_max=1000.0)
con_errors = np.clip(con_errors,a_min=np.finfo(float).eps,a_max=1000.0)

print("Worst primitive error: " + str(np.max(errors)), "Iteration: " + str(np.argmax(errors)))
print("Worst conservative error: " + str(np.max(con_errors)), "Iteration: " + str(np.argmax(con_errors)))

import matplotlib.pyplot as plt

def plot_scatter_log_error(var,err,var_col=None,log_ax="",labx=None,laby=None,title=None,fig=None,ax=None):
    try:
        var_col[0]
    except:
        var_col=np.ones(len(var))
    if "x" in log_ax:
        x = np.log10(var)
    else:
        x = var
    if "y" in log_ax:
        y = np.log10(err)
    else:
        y = err
    y_clone = y.copy()
    y=y[np.isfinite(y_clone)]
    x=x[np.isfinite(y_clone)]
    if not fig:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
    ax.scatter(x,y,marker="x",c=var_col[np.isfinite(y_clone)],cmap="inferno")
    ax.set_xlim(np.min(x),np.max(x))
    ax.set_ylim(np.min(y),np.max(y))
    if labx: ax.set_xlabel(labx)
    if laby: ax.set_ylabel(laby)
    if title: ax.set_title(title)
    return fig, ax


err_sort = np.argsort(con_errors)
errors = errors[err_sort]
con_errors = con_errors[err_sort]
rhos = rhos[err_sort]
presss = presss[err_sort]
Yes = Yes[err_sort]
iterations = iterations[err_sort]
Ts = Ts[err_sort]
epss = epss[err_sort]
speeds = speeds[err_sort]
guesses_array = guesses_array[err_sort]
targets_array = targets_array[err_sort]
results_array = results_array[err_sort]
guess_press = guess_press[err_sort]

print("Failure count: " + str(np.sum(np.log10(con_errors)>-6)))

plot_any=False
if plot_any==True:
    fig_TRY = plt.figure(figsize=(8,8))
    ax1 = fig_TRY.add_subplot(2,2,1)
    ax2 = fig_TRY.add_subplot(2,2,2)
    ax3 = fig_TRY.add_subplot(2,2,3)
    ax4 = fig_TRY.add_subplot(2,2,4)
    fig_TRY, ax1 = plot_scatter_log_error(iterations,con_errors,log_ax="y",labx="Iterations",laby="log10(err)",fig=fig_TRY,ax=ax1)
    fig_TRY, ax2 = plot_scatter_log_error(rhos,Ts,np.log10(con_errors),log_ax="xy",labx="log10(rho)",laby="log10(T)",fig=fig_TRY,ax=ax2)
    fig_TRY, ax3 = plot_scatter_log_error(rhos,Yes,np.log10(con_errors),log_ax="x",labx="log10(rho)",laby="Ye",fig=fig_TRY,ax=ax3)
    fig_TRY, ax4 = plot_scatter_log_error(Yes,Ts,np.log10(con_errors),log_ax="y",labx="Ye",laby="log10(T)",fig=fig_TRY,ax=ax4)
    ax2.scatter(np.log10(rhos)[np.log10(con_errors)>-6],np.log10(Ts)[np.log10(con_errors)>-3],s=120,edgecolors="k",marker="o",facecolors="k")
    ax3.scatter(np.log10(rhos)[np.log10(con_errors)>-6],Yes[np.log10(con_errors)>-6],s=120,edgecolors="k",marker="o",facecolors="k")
    ax4.scatter(Yes[np.log10(con_errors)>-6],np.log10(Ts)[np.log10(con_errors)>-6],s=120,edgecolors="k",marker="o",facecolors="k")
    ax2.scatter(np.log10(rhos)[np.log10(con_errors)>-6],np.log10(Ts)[np.log10(con_errors)>-6],edgecolors="w",marker="+",facecolors="w")
    ax3.scatter(np.log10(rhos)[np.log10(con_errors)>-6],Yes[np.log10(con_errors)>-6],edgecolors="w",marker="+",facecolors="w")
    ax4.scatter(Yes[np.log10(con_errors)>-6],np.log10(Ts)[np.log10(con_errors)>-6],edgecolors="w",marker="+",facecolors="w")
    fig_TRY.tight_layout()

    fig = plt.figure(figsize=(8,4.5))
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)
    ax1.hist(np.log10(Ts),range=(log_T_array.min(),log_T_array.max()),bins=50,density=True)
    ax1.set_xlabel("log10(T)")
    ax1.set_ylabel("Probability Density")
    ax2.hist(np.log10(rhos),range=(log_rho_array.min(),log_rho_array.max()),bins=50,density=True)
    ax2.set_xlabel("log10(rho)")
    ax2.set_ylabel("Probability Density")
    ax3.hist(Yes,bins=50,range=(Ye_array.min(),Ye_array.max()),density=True)
    ax3.set_xlabel("Ye")
    ax3.set_ylabel("Probability Density")
    fig.tight_layout()

all_plots = False
if all_plots == True:
    plot_scatter_log_error(rhos,errors,log_ax="xy",labx="log10(rho)",laby="log10(err)")
    plot_scatter_log_error(presss,RelDiff(presss,guess_press),np.log10(errors),log_ax="xy",labx="log10(press)",laby="log10(abs delta press)")
    plot_scatter_log_error(rhos,Yes,np.log10(errors),log_ax="x",labx="log10(rho)",laby="Ye")
    plot_scatter_log_error(Yes,np.log10(errors),log_ax="",labx="Ye",laby="log10(err)")
    plot_scatter_log_error(rhos,errors,Ts,log_ax="xy",labx="log10(rho)",laby="log10(err)")    
    plot_scatter_log_error(presss,errors,log_ax="xy",labx="log10(press)",laby="log10(err)")
    plot_scatter_log_error(epss,errors,log_ax="y",labx="eps",laby="log10(err)")
    plot_scatter_log_error(speeds,errors,log_ax="y",labx="speed",laby="log10(err)")    
    plot_scatter_log_error(rhos,Ts,iterations,log_ax="xy",labx="log10(rho)",laby="log10(T)")
    plot_scatter_log_error(rhos,speeds,iterations,log_ax="x",labx="log10(rho)",laby="speed")
    plot_scatter_log_error(speeds,Ts,iterations,log_ax="y",labx="speed",laby="log10(T)")
    plot_scatter_log_error(RelDiff(targets_array[:,0],guesses_array[:,0]),errors,log_ax="y",labx="W Guess Err",laby="log10(err)")
    plot_scatter_log_error(RelDiff(targets_array[:,1],guesses_array[:,1]),errors,log_ax="y",labx="z Guess Err",laby="log10(err)")
    plot_scatter_log_error(RelDiff(targets_array[:,2],guesses_array[:,2]),errors,log_ax="y",labx="T Guess Err",laby="log10(err)")
    plot_scatter_log_error(RelDiff(targets_array[:,0],guesses_array[:,0]),RelDiff(targets_array[:,1],guesses_array[:,1]),np.log10(errors),log_ax="",labx="W Guess Err",laby="z Guess Err")
    plot_scatter_log_error(RelDiff(targets_array[:,1],guesses_array[:,1]),RelDiff(targets_array[:,2],guesses_array[:,2]),np.log10(errors),log_ax="",labx="z Guess Err",laby="T Guess Err")
    plot_scatter_log_error(RelDiff(targets_array[:,2],guesses_array[:,2]),RelDiff(targets_array[:,0],guesses_array[:,0]),np.log10(errors),log_ax="",labx="T Guess Err",laby="W Guess Err")
#from C2P_InitialGuess import W_target, T_target, h_target, rho_target, press_target, eps_target
#from C2P_InitialGuess import v1_target, v2_target, v3_target


#z_target = rho_target*h_target*(W_target**2)

#fig_wzt = plt.figure(1)
#ax_wzt = fig_wzt.add_subplot(1,1,1)
#ax_wzt.plot((W_array[:it+1]-W_target)/W_target)
#ax_wzt.plot((z_array[:it+1]-z_target)/z_target)
#ax_wzt.plot((T_array[:it+1]-T_target)/T_target)

#fig_f = plt.figure(2)
#ax_f = fig_f.add_subplot(1,1,1)
#ax_f.plot(f1_array[:it+1]/np.max(np.abs(f1_array)))
#ax_f.plot(f2_array[:it+1]/np.max(np.abs(f2_array)))
#ax_f.plot(f3_array[:it+1]/np.max(np.abs(f3_array)))

#fig_prim = plt.figure(3)
#ax_press = fig_prim.add_subplot(1,1,1)
#ax_press.plot(press_array[:it+1])