# -*- coding: utf-8 -*-
"""
Created on Fri May 24 16:18:38 2019

@author: User
"""

import numpy as np
from scipy.optimize import minimize
from C2P_GlobalVars import EOS_Eval_Order
from C2P_Setup import press_eval_kwargs, eps_eval_kwargs
from C2P_Setup import log_T_array, log_rho_array, Ye_array

### EoS primitive ###

def GetEOSHandle(name):
    if name=="Table3D":
        return 0
    else:
        raise ValueError("EoS name unknown")
        return None

EOS_Type = 0

def EOS_pressFromPrim(T,rho,Ye):
    press = np.inf
    if EOS_Type == 0:
        press = EOS_Eval(np.log10(T),np.log10(rho),Ye,log_T_array,log_rho_array,Ye_array,**press_eval_kwargs)
    return press

def EOS_epsFromPrim(T,rho,Ye):
    eps = np.inf
    if EOS_Type == 0:
        eps = EOS_Eval(np.log10(T),np.log10(rho),Ye,log_T_array,log_rho_array,Ye_array,**eps_eval_kwargs)
    return eps

def EOS_dPdrhoFromPrim(T,rho,Ye):
    dPdrho = np.inf
    if EOS_Type==0:
        dPdlrho = EOS_dVdlrho(np.log10(T),np.log10(rho),Ye,log_T_array,log_rho_array,Ye_array,**press_eval_kwargs)
        dPdrho  = (1/(rho*np.log(10))) * dPdlrho
    return dPdrho

def EOS_dEdrhoFromPrim(T,rho,Ye):
    dEdrho = np.inf
    if EOS_Type==0:
        dEdlrho = EOS_dVdlrho(np.log10(T),np.log10(rho),Ye,log_T_array,log_rho_array,Ye_array,**eps_eval_kwargs)
        dEdrho  = (1/(rho*np.log(10))) * dEdlrho
    return dEdrho

def EOS_dPdTFromPrim(T,rho,Ye):
    dPdT = np.inf
    if EOS_Type==0:
        dPdlT = EOS_dVdlT(np.log10(T),np.log10(rho),Ye,log_T_array,log_rho_array,Ye_array,**press_eval_kwargs)
        dPdT  = (1/(T*np.log(10))) * dPdlT
    return dPdT

def EOS_dEdTFromPrim(T,rho,Ye):
    dEdT = np.inf
    if EOS_Type==0:
        dEdlT = EOS_dVdlT(np.log10(T),np.log10(rho),Ye,log_T_array,log_rho_array,Ye_array,**eps_eval_kwargs)
        dEdT  = (1/(T*np.log(10))) * dEdlT
    return dEdT

### EoS conservative ###

def EOS_pressFromWzTCons(WzT,Cons):
    W = WzT["W"]
    T = WzT["T"]
    D = Cons["D"]
    DYe = Cons["DYe"]
    rho = D/W
    Ye = DYe/D
    return EOS_pressFromPrim(T,rho,Ye)

def EOS_epsFromWzTCons(WzT,Cons):
    W = WzT["W"]
    T = WzT["T"]
    D = Cons["D"]
    DYe = Cons["DYe"]
    rho = D/W
    Ye = DYe/D
    return EOS_epsFromPrim(T,rho,Ye)

def EOS_dPdWFromWzTCons(WzT,Cons):
    W = WzT["W"]
    T = WzT["T"]
    D = Cons["D"]
    DYe = Cons["DYe"]
    rho = D/W
    Ye = DYe/D
    return -D*(W**(-2))*EOS_dPdrhoFromPrim(T,rho,Ye)

def EOS_dEdWFromWzTCons(WzT,Cons):
    W = WzT["W"]
    T = WzT["T"]
    D = Cons["D"]
    DYe = Cons["DYe"]
    rho = D/W
    Ye = DYe/D
    return -D*(W**(-2))*EOS_dEdrhoFromPrim(T,rho,Ye)

def EOS_dPdTFromWzTCons(WzT,Cons):
    W = WzT["W"]
    T = WzT["T"]
    D = Cons["D"]
    DYe = Cons["DYe"]
    rho = D/W
    Ye = DYe/D
    return EOS_dPdTFromPrim(T,rho,Ye)

def EOS_dEdTFromWzTCons(WzT,Cons):
    W = WzT["W"]
    T = WzT["T"]
    D = Cons["D"]
    DYe = Cons["DYe"]
    rho = D/W
    Ye = DYe/D
    return EOS_dEdTFromPrim(T,rho,Ye)

def EOS_S2FromWzTCons(WzT,Cons):
    W = WzT["W"]
    z = WzT["z"]
    D = Cons["D"]
    rho = D/W
    #press = EOS_pressFromWzTCons(WzT,Cons)
    #eps = EOS_epsFromWzTCons(WzT,Cons)
    #h = 1 + eps + press/rho
    h = z/(rho*W*W)
    v2 = 1-(W**(-2))
    return ((rho*h)**2)*(W**4)*(v2)

### 3-D Tabulated EoS ###

def EOS_Eval_Lin(T,rho,Ye,T_array,rho_array,Ye_array,table=None,coeffs=None):
    d_T   = T_array[1]   - T_array[0]
    d_rho = rho_array[1] - rho_array[0]
    d_Ye  = Ye_array[1]  - Ye_array[0]
    
    T_bin   = int((T   - T_array[0])//d_T)
    rho_bin = int((rho - rho_array[0])//d_rho)
    Ye_bin  = int((Ye  - Ye_array[0])//d_Ye)
    
    #print(T_bin,rho_bin,Ye_bin)
    
    xd = (T   - T_array[T_bin])/d_T
    yd = (rho - rho_array[rho_bin])/d_rho
    zd = (Ye  - Ye_array[Ye_bin])/d_Ye
    
    c_000 = table[T_bin  , rho_bin  , Ye_bin  ]
    c_001 = table[T_bin  , rho_bin  , Ye_bin+1]
    c_010 = table[T_bin  , rho_bin+1, Ye_bin  ]
    c_011 = table[T_bin  , rho_bin+1, Ye_bin+1]
    c_100 = table[T_bin+1, rho_bin  , Ye_bin  ]
    c_101 = table[T_bin+1, rho_bin  , Ye_bin+1]
    c_110 = table[T_bin+1, rho_bin+1, Ye_bin  ]
    c_111 = table[T_bin+1, rho_bin+1, Ye_bin+1]
    
    c_x00 = c_000*(1-xd) + c_100*xd
    c_x01 = c_001*(1-xd) + c_101*xd
    c_x10 = c_010*(1-xd) + c_110*xd
    c_x11 = c_011*(1-xd) + c_111*xd
    
    c_xy0 = c_x00*(1-yd) + c_x10*yd
    c_xy1 = c_x01*(1-yd) + c_x11*yd
    
    c_xyz = c_xy0*(1-zd) + c_xy1*zd
    
    return c_xyz

def EOS_dVdlrho_Lin(lT,lrho,Ye,lT_array,lrho_array,Ye_array,table=None,coeffs=None):
    d_lT   = lT_array[1]   - lT_array[0]
    d_lrho = lrho_array[1] - lrho_array[0]
    d_Ye  = Ye_array[1]  - Ye_array[0]
    
    lT_bin   = int((lT   - lT_array[0])//d_lT)
    lrho_bin = int((lrho - lrho_array[0])//d_lrho)
    Ye_bin  = int((Ye  - Ye_array[0])//d_Ye)
    
    #print(T_bin,rho_bin,Ye_bin)
    
    xd = (lT   - lT_array[lT_bin])/d_lT
    yd = (lrho - lrho_array[lrho_bin])/d_lrho
    zd = (Ye  - Ye_array[Ye_bin])/d_Ye
    
    c_000 = table[lT_bin  , lrho_bin  , Ye_bin  ]
    c_001 = table[lT_bin  , lrho_bin  , Ye_bin+1]
    c_010 = table[lT_bin  , lrho_bin+1, Ye_bin  ]
    c_011 = table[lT_bin  , lrho_bin+1, Ye_bin+1]
    c_100 = table[lT_bin+1, lrho_bin  , Ye_bin  ]
    c_101 = table[lT_bin+1, lrho_bin  , Ye_bin+1]
    c_110 = table[lT_bin+1, lrho_bin+1, Ye_bin  ]
    c_111 = table[lT_bin+1, lrho_bin+1, Ye_bin+1]
    
    c_x00 = c_000*(1-xd) + c_100*xd
    c_x01 = c_001*(1-xd) + c_101*xd
    c_x10 = c_010*(1-xd) + c_110*xd
    c_x11 = c_011*(1-xd) + c_111*xd
    
    c_x0z = c_x00*(1-zd) + c_x01*zd
    c_x1z = c_x10*(1-zd) + c_x11*zd
    
    delta_c_xYz = c_x1z - c_x0z
    
    return delta_c_xYz/d_lrho
    
def EOS_dVdlT_Lin(lT,lrho,Ye,lT_array,lrho_array,Ye_array,table=None,coeffs=None):
    d_lT   = lT_array[1]   - lT_array[0]
    d_lrho = lrho_array[1] - lrho_array[0]
    d_Ye  = Ye_array[1]  - Ye_array[0]
    
    lT_bin   = int((lT   - lT_array[0])//d_lT)
    lrho_bin = int((lrho - lrho_array[0])//d_lrho)
    Ye_bin  = int((Ye  - Ye_array[0])//d_Ye)
    
    #print(T_bin,rho_bin,Ye_bin)
    
    xd = (lT   - lT_array[lT_bin])/d_lT
    yd = (lrho - lrho_array[lrho_bin])/d_lrho
    zd = (Ye  - Ye_array[Ye_bin])/d_Ye
    
    c_000 = table[lT_bin  , lrho_bin  , Ye_bin  ]
    c_001 = table[lT_bin  , lrho_bin  , Ye_bin+1]
    c_010 = table[lT_bin  , lrho_bin+1, Ye_bin  ]
    c_011 = table[lT_bin  , lrho_bin+1, Ye_bin+1]
    c_100 = table[lT_bin+1, lrho_bin  , Ye_bin  ]
    c_101 = table[lT_bin+1, lrho_bin  , Ye_bin+1]
    c_110 = table[lT_bin+1, lrho_bin+1, Ye_bin  ]
    c_111 = table[lT_bin+1, lrho_bin+1, Ye_bin+1]
    
    c_0y0 = c_000*(1-yd) + c_010*yd
    c_0y1 = c_001*(1-yd) + c_011*yd
    c_1y0 = c_100*(1-yd) + c_110*yd
    c_1y1 = c_101*(1-yd) + c_111*yd
    
    c_0yz = c_0y0*(1-zd) + c_0y1*zd
    c_1yz = c_1y0*(1-zd) + c_1y1*zd
    
    delta_c_Xyz = c_1yz - c_0yz
    
    return delta_c_Xyz/d_lT

if EOS_Eval_Order == 1:
    EOS_Eval = EOS_Eval_Lin
    EOS_dVdlT = EOS_dVdlT_Lin
    EOS_dVdlrho = EOS_dVdlrho_Lin

def Residual(eps_cold,epss,Ts):
    eps_thermal = epss-eps_cold
    if np.any(eps_thermal<=0.0):
        return 1e100
    xs = np.log10(Ts)
    ys = np.log10(eps_thermal)
    coeffs = np.polyfit(xs,ys,deg=1)
    ys_polyfit = coeffs[0]*xs + coeffs[1]
    residuals = (ys-ys_polyfit)**2
    return np.sum(residuals), coeffs[1], coeffs[0]

#def Residual(eps_cold,epss,Ts):
#    eps_thermal = epss-eps_cold
#    if np.any(eps_thermal<=0.0):
#        return 1e100
#    xs = np.log10(Ts)
#    ys = np.log10(eps_thermal)
#    slope = 2.0
#    ys_temp = ys-slope*xs
#    intercept = np.average(ys_temp)
#    ys_fit = slope*xs + intercept
#    residuals = (ys-ys_fit)**2
#    return np.sum(residuals), intercept, 2.0

def GoldenSectionSearch(a,b,epss,Ts,tol=1e-15):        
    invphi = (np.sqrt(5) - 1) / 2 # 1/phi
    invphi2 = (3 - np.sqrt(5)) / 2 # 1/phi^2

    (a,b)=(min(a,b),max(a,b))
    h = b - a
    if h <= tol: return (a,b)

    # required steps to achieve tolerance                                                                                                                   
    n = int(np.ceil(np.log(tol/h)/np.log(invphi)))

    c = a + invphi2 * h
    d = a + invphi * h
    yc, intercept, power = Residual(c,epss,Ts)
    yd, intercept, power = Residual(d,epss,Ts)

    for k in np.arange(n-1):
        if yc < yd:
            b = d
            d = c
            yd = yc
            h = invphi*h
            c = a + invphi2 * h
            yc, intercept, power = Residual(c,epss,Ts)
        else:
            a = c
            c = d
            yc = yd
            h = invphi*h
            d = a + invphi * h
            yd, intercept, power = Residual(d,epss,Ts)

    if yc < yd:
        return (a,d)
    else:
        return (c,b)

def GetEpsCold(Ts_full,rho,Ye):
    #Ts = Ts_full.copy()
    Ts = Ts_full[Ts_full<=0.5]
    epss = np.zeros(len(Ts))
    for T_idx in range(len(Ts)):
        epss[T_idx] = EOS_epsFromPrim(Ts[T_idx],rho,Ye)
    eps_cold_max = epss[0]# - 1e-12*np.abs(epss[0])
    #eps_cold_min = epss[0] - 2e-16*np.abs(epss[0])
    eps_1 = EOS_epsFromPrim(Ts[0],rho,Ye)
    eps_2 = EOS_epsFromPrim(2*Ts[0],rho,Ye)
    eps_cold_min = 2*eps_1 - eps_2
    eps_cold_interval = GoldenSectionSearch(eps_cold_min,eps_cold_max,epss,Ts,tol=1e-14)
    eps_cold_interval_mid = 0.5*(eps_cold_interval[0] + eps_cold_interval[1])
    #eps_cold_soln = minimize(Residual,[eps_cold_interval_mid],args=(epss,Ts),bounds=[(eps_cold_interval[0],eps_cold_interval[1])],tol=1e-20)
    #eps_cold = eps_cold_soln.x[0]
    #alpha = Residual(eps_cold_interval[0],epss,Ts)
    #beta  = Residual(eps_cold_interval_mid,epss,Ts)
    #gamma = Residual(eps_cold_interval[1],epss,Ts)
    
    #if (alpha - 2*beta + gamma)==0.0:
    #    print(eps_cold_max,eps_cold_min)
    #    print(Residual(eps_cold_min,epss,Ts))
    #    print(Residual(eps_cold_max,epss,Ts))
    #    print(alpha,beta,gamma)
    #    print(eps_cold_interval)
    #    print(eps_cold_interval_mid)
    #    for i in range(100):
    #        print(Ts[i],epss[i])
    #    print(epss)
    #    import matplotlib.pyplot as plt
    #    plt.loglog(Ts,epss-eps_cold_min)
    #    plt.loglog(Ts,epss-eps_cold_interval_mid)
    #    plt.loglog(Ts,epss-eps_cold_max)
    #
    #p = 0.5*(alpha-gamma)/(alpha - 2*beta + gamma)
    #
    #eps_cold = eps_cold_interval_mid + p*(eps_cold_interval_mid-eps_cold_interval[0])
    #if (alpha - 2*beta + gamma)==0.0:
    #    plt.loglog(Ts,epss-eps_cold)
    
    eps_cold = eps_cold_interval_mid
    residual, intercept, power = Residual(eps_cold,epss,Ts)
    return eps_cold, residual, intercept, power

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '#'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

if __name__=="__main__":
    T_samples = 1001
    rho_samples = 512
    Ye_samples = 8
    
    test_T_array = 10**(np.linspace(log_T_array[0],log_T_array[-1],num=T_samples+2)[1:-1])
    test_rho_array = 10**(np.linspace((log_rho_array[log_rho_array>-12])[0],(log_rho_array[log_rho_array<-2.5])[-1],num=rho_samples+2)[1:-1])
    test_Ye_array = np.linspace(Ye_array[0],Ye_array[-1],num=Ye_samples+2)[1:-1]
    
    eps_colds = np.zeros((rho_samples,Ye_samples))
    residuals = np.zeros((rho_samples,Ye_samples))
    intercepts = np.zeros((rho_samples,Ye_samples))
    powers = np.zeros((rho_samples,Ye_samples))
    
    import matplotlib.pyplot as plt
    
    print("Calculating eps_cold array")
    for rho_idx in range(rho_samples):
        #print(rho_idx)
        rho_current = test_rho_array[rho_idx]
        for Ye_idx in range(Ye_samples):
            Ye_current = test_Ye_array[Ye_idx]
            eps_cold_current, residual_current, intercept_current, power_current = GetEpsCold(test_T_array,rho_current,Ye_current)
            eps_colds[rho_idx,Ye_idx] = eps_cold_current
            residuals[rho_idx,Ye_idx] = residual_current
            intercepts[rho_idx,Ye_idx] = intercept_current
            powers[rho_idx,Ye_idx] = power_current
            
            #epss = np.zeros(T_samples)
            #for T_idx in range(T_samples):
            #    epss[T_idx] = EOS_epsFromPrim(test_T_array[T_idx],rho_current,Ye_current)
            
            #fig = plt.figure(1,clear=True)
            #ax = fig.add_subplot(1,1,1)
            #ax.set_title("eps_thermal: rho=" + str(rho_current) + " Ye=" + str(Ye_current))
            #ax.set_xlabel("log10(T)")
            #ax.set_ylabel("eps_thermal")
            #
            #ax.loglog(test_T_array,epss-eps_cold_current)
            #ax.loglog(test_T_array,10**(intercept_current)*(test_T_array**power_current))
            #fig.savefig("./Figs/epsilon_" + "{:04d}".format(rho_idx) + "_" + "{:04d}".format(Ye_idx) + ".png",dpi=480)
            printProgressBar(rho_idx*Ye_samples + Ye_idx + 1,rho_samples*Ye_samples)
    
    T_extrap_samples = 1000
    log_T_extrap_array = np.linspace(-5,-1,num=T_extrap_samples,endpoint=False)
    d_lT = log_T_extrap_array[1]-log_T_extrap_array[0]
    log_T_EOS_array = np.arange(-1.0,log_T_array[-1],d_lT)
    
    T_extrap_array = 10**log_T_extrap_array
    T_EOS_array = 10**log_T_EOS_array
    T_extrap_rho_array = test_rho_array
    
    Ye_idx = np.argmin(np.abs(test_Ye_array-0.46198265933111149))
    Ye = test_Ye_array[Ye_idx]
    test_eps_th_array = np.zeros((rho_samples,T_extrap_samples+len(T_EOS_array)))
    for rho_idx in range(rho_samples):
        for T_idx in range(T_extrap_samples):
            test_eps_th_array[rho_idx,T_idx] = 10**(intercepts[rho_idx,Ye_idx])*(T_extrap_array[T_idx]**powers[rho_idx,Ye_idx])# + eps_colds[rho_idx,Ye_idx]
        for T_idx in range(T_extrap_samples,T_extrap_samples+len(T_EOS_array)):
            test_eps_th_array[rho_idx,T_idx] = EOS_epsFromPrim(T_EOS_array[T_idx-T_extrap_samples],T_extrap_rho_array[rho_idx],Ye) - eps_colds[rho_idx,Ye_idx]
    
    plot=True
    save_plot = True
    if plot==True:        
        import matplotlib.pyplot as plt
        fig1 = plt.figure()
        ax_1 = fig1.add_subplot(1,1,1)
        ax_1.set_title("eps_cold")
        ax_1.set_xlabel("log10(rho)")
        ax_1.set_ylabel("Ye")
        dlr = np.log10(test_rho_array[1]) - np.log10(test_rho_array[0])
        dYe = test_Ye_array[1] - test_Ye_array[0]
        imshow_extent = [np.log10(test_rho_array[0])-0.5*dlr,
                         np.log10(test_rho_array[-1])+0.5*dlr,
                         test_Ye_array[0] - 0.5*dYe,
                         test_Ye_array[-1] + 0.5*dYe]
        im1 = ax_1.imshow(eps_colds.T,origin="lower",cmap="inferno",extent=imshow_extent,aspect="auto")
        fig1.colorbar(im1,ax=ax_1)
        
        fig2 = plt.figure()
        ax_2 = fig2.add_subplot(1,1,1)
        ax_2.set_title("k")
        ax_2.set_xlabel("log10(rho)")
        ax_2.set_ylabel("Ye")
        im2 = ax_2.imshow(intercepts.T,origin="lower",cmap="inferno",extent=imshow_extent,aspect="auto")
        fig2.colorbar(im2,ax=ax_2)
        
        fig3 = plt.figure()
        ax_3 = fig3.add_subplot(1,1,1)
        ax_3.set_title("err")
        ax_3.set_xlabel("log10(rho)")
        ax_3.set_ylabel("Ye")
        im3 = ax_3.imshow(residuals.T,origin="lower",cmap="inferno",extent=imshow_extent,aspect="auto")
        fig3.colorbar(im3,ax=ax_3)
        
        fig4 = plt.figure()
        ax_4 = fig4.add_subplot(1,1,1)
        ax_4.set_title("power")
        ax_4.set_xlabel("log10(rho)")
        ax_4.set_ylabel("Ye")
        im4 = ax_4.imshow(powers.T,origin="lower",cmap="inferno",extent=imshow_extent,aspect="auto")
        fig4.colorbar(im4,ax=ax_4)
        
        d_lr = np.log10(T_extrap_rho_array[1]) - np.log10(T_extrap_rho_array[0])
        d_lT = log_T_extrap_array[1] - log_T_extrap_array[0]
        imshow_extent_eps_th = [np.log10(T_extrap_rho_array[0])  - d_lr,
                                np.log10(T_extrap_rho_array[-1]) + d_lr,
                                log_T_extrap_array[0] - d_lT,
                                log_T_EOS_array[-1] + d_lT]
        fig5 = plt.figure()
        ax_5 = fig5.add_subplot(1,1,1)
        ax_5.set_title("eps_thermal")
        ax_5.set_xlabel("log10(rho)")
        ax_5.set_ylabel("log10(T)")
        im5 = ax_5.imshow(np.log10(test_eps_th_array.T),origin="lower",cmap="inferno",extent=imshow_extent_eps_th,aspect="auto")
        fig5.colorbar(im5,ax=ax_5)
        
        if save_plot==True: 
            fig1.savefig("eps_cold.png",dpi=480)
            fig2.savefig("k.png",dpi=480)
            fig3.savefig("err.png",dpi=480)
            fig4.savefig("power.png",dpi=480)
            fig5.savefig("eps_th.png",dpi=480)
            
#    epss = np.zeros(len(Ts))
#    for T_idx in range(len(Ts)):
#        epss[T_idx] = EOS_epsFromPrim(Ts[T_idx],rho,Ye)
#    
#    eps_cold_max = epss[0] - 1e-16*np.abs(epss[0])
#    eps_cold_min = -1
#    
#    eps_cold_guess = 0.5*(eps_cold_min+eps_cold_max)
#    
#    
#    print(Residual(eps_cold_max,epss,Ts))
#    print(Residual(eps_cold_guess,epss,Ts))
#    print(Residual(eps_cold_min,epss,Ts))
#    
#    eps_cold_interval = GoldenSectionSearch(eps_cold_min,eps_cold_max,epss,Ts,tol=1e-11)
#    int_minus = eps_cold_interval[0]
#    int_plus  = eps_cold_interval[1]
#    res_int_minus = Residual(int_minus,epss,Ts)
#    res_int_plus = Residual(int_plus,epss,Ts)
#    eps_cold_interval_mid = 0.5*(eps_cold_interval[0] + eps_cold_interval[1])
#    eps_cold_soln = minimize(Residual,[eps_cold_interval_mid],args=(epss,Ts),bounds=[(eps_cold_interval[0],eps_cold_interval[1])],tol=1e-20)
#    eps_cold = eps_cold_soln.x[0]
#    
#    xs = np.log10(Ts)
#    ys = np.log10(epss-eps_cold)
#    coeffs = np.polyfit(xs,ys,deg=1)
#    ys_polyfit = coeffs[0]*xs + coeffs[1] + eps_cold
#    
#    print(Residual(eps_cold,epss,Ts))
#    print(eps_cold)
#    
#    import matplotlib.pyplot as plt
#    plt.ion()
#    plt.figure()
#    plt.loglog(Ts,epss-eps_cold_max)
#    plt.loglog(Ts,epss-eps_cold_min)
#    plt.loglog(Ts,epss-eps_cold)
#    plt.loglog(10**xs,10**ys_polyfit)
#    
#    eps_colds = np.linspace(eps_cold_min,eps_cold_max,num=10000)
#    resids = np.zeros(len(eps_colds))
#    for i in range(len(resids)):
#        resids[i] = Residual(eps_colds[i],epss,Ts)
#    plt.figure()
#    plt.plot(eps_colds,resids)
#    plt.scatter(eps_cold,Residual(eps_cold,epss,Ts))
#    
#    
#    eps_colds = np.linspace(eps_cold_interval[0],eps_cold_interval[1],num=10000)
#    resids = np.zeros(len(eps_colds))
#    for i in range(len(resids)):
#        resids[i] = Residual(eps_colds[i],epss,Ts)
#    plt.figure()
#    plt.plot(eps_colds,resids)
#    plt.scatter(int_minus,res_int_minus)
#    plt.scatter(int_plus,res_int_plus)
#    plt.scatter(eps_cold,Residual(eps_cold,epss,Ts))
#    
#    print(Residual(eps_colds[np.argmax(resids)],epss,Ts))
#    