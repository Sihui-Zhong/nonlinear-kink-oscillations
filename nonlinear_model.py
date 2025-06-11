def nonlinear_model(t, R, roi, period_k, period_A, Vi, delta_m, zeta, C1=0.3):
    """nonlinear damping function by KHI-induced turbulence 
    t: time
    R: loop radius in km
    roi: internal density compared to external density roe=1 ;dimensionless
    zeta: unit mass of mixing layer used to calculate the M_bar, ;if roi=3,zeta=0.24
    period_k: in seconds, oscillation period of kink mode; wk=2*pi/period_k angular frequency of kink mode
    period_A: -->wA: Alfven frequency in the mixing layer; in second
    Vi: initial velocity perturbation, in km/s
    delta_m: dimensionless, mass loss/gain rate in then area defined by rhoT
    C1: a constant describing the physics of mixing in hydrodynamics, dimensionless, default=0.3
    Calling sequence:
    t=np.arange(0,130,0.5)
    y = nonlinear_model(t,R=1000,roi=3,period_k=300,period_A=280,Vi=30,delta_m=-0.85,zeta=0.24,C1=0.3)
    """
    import math
    import numpy as np
    from scipy.special import sici
    
    pi = math.pi
    roe = 1.0
    C2 = np.sqrt(2)*C1
    wk = 2*pi/period_k
    wA = 2*pi/period_A
    wD = (wA-wk)/2.0
    wS = (wA+wk)/2.0
    dV = Vi*1.1
    H = pi*R/4 #depth
    E = 4*R*H*roi #mass in the core
    A = E*Vi
    M_bar = (Vi - dV*np.sqrt(roe)/(np.sqrt(roi)+np.sqrt(roe)))*zeta
    #B--momentum loss rate of the tube core
    B = 1/(2*H)*C2/np.sqrt(2)*np.sqrt(roe)/(np.sqrt(roi)+np.sqrt(roe))*(roi*roe)**(0.25)/(np.sqrt(roi)+np.sqrt(roe))*dV
    #print('B=',B,'1/B (time at zero M_core)=',1/B)
    D = 2*M_bar*C2/np.sqrt(2)*(roi*roe)**(0.25)/(np.sqrt(roi)+np.sqrt(roe))*dV*R/wD
    #G --mass gain/loss rate in the area defined by rho_T-->delta_m
    G = delta_m/(2*H*roi)*C2/np.sqrt(2)*(roi*roe)**(0.25)/(np.sqrt(roi)+np.sqrt(roe))*dV
    A1 = A/(E*(G**2))*(-B*G*np.sin(wk*t)/wk+(B+G)*(np.cos(wk/G)*sici(wk*(t+1/G))[1]+np.sin(wk/G)*sici(wk*(t+1/G))[0]))
    A2 = D/(2*E*G)*(-np.cos(wk/G)*sici(wk*(t+1/G))[0]+np.sin(wk/G)*sici(wk*(t+1/G))[1]
               +np.cos(wA/G)*sici(wA*(t+1/G))[0]-np.sin(wA/G)*sici(wA*(t+1/G))[1])
    A1 -= A1[0]
    A2 -= A2[0]
    return A1+A2
