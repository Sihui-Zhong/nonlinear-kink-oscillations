;SoBAT MCMC fitting with nonlinear damping function by turbulence
;Author: Sihui Zhong (s.zhong3@exeter.ac.uk, sihui.zhong@outlook.com)

function nonlinear_model,t,pars,_extra=_extra
  compile_opt idl2
  ;9-parameter function describing the transverse oscillations nonlinearily damped by turbulence
  
  R = pars[0] ;loop radius 
  H = !dpi*R/4
  roi = pars[1] ;internal density, in unit of the external background density
  roe = 1. ;external
  C1 = pars[2]
  C2 = sqrt(2)*C1
  Vi = pars[3] ;initial velocity
  dV = 1.1*Vi ;delta_V
  period_k = pars[4] ;period of kink oscillations
  wk = 2*!dpi/period_k ;wk, angular frequency, =2*!pi/P
  period_A= pars[5]
  wA = 2*!dpi/period_A 
  wD = (wA-wk)/2.0
  wS = (wA+wk)/2.0
  delta_m = pars[6]
  zeta = pars[7] ;0.9
  t0 = pars[8]
  
  t1 = t-t0
  M_bar = zeta*(Vi - dV*sqrt(roe)/(sqrt(roi)+sqrt(roe))) ;0.24
  
  E = 4*R*H*roi ;mass
  A = E*Vi ;momentum
  B = 1/(2*H)*C2/sqrt(2)*sqrt(roe)/(sqrt(roi)+sqrt(roe))*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV
  D = 2*M_bar*C2/sqrt(2)*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV*R/wD
  G = delta_m/(2*H*roi)*C2/sqrt(2)*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV
  A1 = A/(E*(G^2))*(-B*G*sin(wk*t1)/wk+(B+G)*(cos(wk/G)*(sici(wk*(t1+1/G)))[1,*]+sin(wk/G)*(sici(wk*(t1+1/G)))[0,*]))
  A2 = D/(2*E*G)*(-cos(wk/G)*(sici(wk*(t1+1/G)))[0,*]+sin(wk/G)*(sici(wk*(t1+1/G)))[1,*]+$
    cos(wA/G)*(sici(wA*(t1+1/G)))[0,*]-sin(wA/G)*(sici(wA*(t1+1/G)))[1,*])  
  func = A1+A2-mean(A1+A2) ;A1+A2-A1[0]-A2[0]
  return,func
  
end

function nonlinear_model_v2,t,pars,_extra=_extra
  compile_opt idl2
  ;7-parameter nonlinear function
  ;.compile -v 'sici.pro'
  R = pars[0] ;loop radius in km
  H = !dpi*R/4
  roi = pars[1]
  roe = 1. ;external
  C1 = pars[2]
  C2 = sqrt(2)*C1
  Vi = pars[3] ;initial velocity ; in km/s
  dV = 1.1*Vi ;delta_V
  period_k = pars[4]*1d ;period of kink oscillations
  wk = 2*!dpi/period_k ;wk, angular frequency, =2*!pi/P

  wA = wk * sqrt(roi/roe + 1d)/(sqrt(2)*(roi/roe)^0.25) ;approximate the mean density in the mixing layer as sqrt(roi*roe)
  wD = (wA-wk)/2.0
  wS = (wA+wk)/2.0
  percent = pars[5]

  zeta_max = sqrt(roi*1d) ;mass convervation
  zeta = zeta_max * percent
  max_y1 = 1./(sqrt(roi)+1d) ;length /2.
  delta_m = zeta - max_y1*roi ;don't have to /2. any more

  t0 = pars[6]

  t1 = t-t0
  M_bar = zeta*(Vi - dV*sqrt(roe)/(sqrt(roi)+sqrt(roe))) ;0.24

  E = 4*R*H*roi ;mass
  A = E*Vi ;momentum
  B = 1/(2*H)*C2/sqrt(2)*sqrt(roe)/(sqrt(roi)+sqrt(roe))*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV
  D = 2*M_bar*C2/sqrt(2)*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV*R/wD
  G = delta_m/(2*H*roi)*C2/sqrt(2)*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV
  A1 = A/(E*(G^2))*(-B*G*sin(wk*t1)/wk+(B+G)*(cos(wk/G)*(sici(wk*(t1+1/G)))[1,*]+sin(wk/G)*(sici(wk*(t1+1/G)))[0,*]))
  A2 = D/(2*E*G)*(-cos(wk/G)*(sici(wk*(t1+1/G)))[0,*]+sin(wk/G)*(sici(wk*(t1+1/G)))[1,*]+$
    cos(wA/G)*(sici(wA*(t1+1/G)))[0,*]-sin(wA/G)*(sici(wA*(t1+1/G)))[1,*])
  func = A1+A2-mean(A1+A2) ;A1+A2-A1[0]-A2[0]
  return,func/1d3 ;in Mm

end

function nonlinear_model_v3,t,pars,_extra=_extra
  compile_opt idl2
  ;7-parameter nonlinear function, the density profile is provided in a pre-computed table.
  ;.compile -v 'sici.pro'
  R = pars[0] ;loop radius ;in km
  H = !dpi*R/4
  roi = pars[1]
  roe = 1. ;external
  C1 = pars[2]
  C2 = sqrt(2)*C1
  Vi = pars[3] ;initial velocity ;in km/s
  dV = 1.1*Vi ;delta_V
  period_k = pars[4]*1d ;period of kink oscillations
  wk = 2*!dpi/period_k ;wk, angular frequency, =2*!pi/P
  
  restore,'/home/sihuizhong/Documents/KHI/mcmc/w_A_factor_table.sav' ;rhoi,f ; provide a look-up table to obatin a more accurate ratio of wA/wk
  f_ip = interpol(f,rhoi,roi) 
  wA = wk * f_ip
  wD = (wA-wk)/2.0
  wS = (wA+wk)/2.0
  percent = pars[5]

  zeta_max = sqrt(roi*1d) ;mass convervation
  zeta = zeta_max * percent
  max_y1 = 1./(sqrt(roi)+1d) ;length /2.
  delta_m = zeta - max_y1*roi ;don't have to /2. any more

  t0 = pars[6]

  t1 = t-t0
  M_bar = zeta*(Vi - dV*sqrt(roe)/(sqrt(roi)+sqrt(roe))) ;0.24

  E = 4*R*H*roi ;mass
  A = E*Vi ;momentum
  B = 1/(2*H)*C2/sqrt(2)*sqrt(roe)/(sqrt(roi)+sqrt(roe))*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV
  D = 2*M_bar*C2/sqrt(2)*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV*R/wD
  G = delta_m/(2*H*roi)*C2/sqrt(2)*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV
  A1 = A/(E*(G^2))*(-B*G*sin(wk*t1)/wk+(B+G)*(cos(wk/G)*(sici(wk*(t1+1/G)))[1,*]+sin(wk/G)*(sici(wk*(t1+1/G)))[0,*]))
  A2 = D/(2*E*G)*(-cos(wk/G)*(sici(wk*(t1+1/G)))[0,*]+sin(wk/G)*(sici(wk*(t1+1/G)))[1,*]+$
    cos(wA/G)*(sici(wA*(t1+1/G)))[0,*]-sin(wA/G)*(sici(wA*(t1+1/G)))[1,*])
  func = A1+A2-mean(A1+A2) ;A1+A2-A1[0]-A2[0]
  return,func/1d3 ;in Mm

end

pro mcmc_nonlinear_v3
  ;This procedure fit the data wih 7-parameter nonlinear function
  
  restore,'nl_data_30percentnoise_v7.sav',/v
  y/=1d3
  cgplot,x,y,psym=1

  ;set priors
  priors = [prior_uniform(0.1d3,10d3),$; loop radius in km ;
    prior_uniform(1d,10d),$ ;roi
    prior_normal(0.3d,0.1d),$ ;C1, dimensionless ;prior_uniform(0.01d,0.5d) 
    prior_uniform(30d,100d),$ ;Vi, in km/s
    prior_uniform(280d,320d),$ ; kink period
    prior_uniform(0.00001d,1d),$ ;percent for zeta
    prior_uniform(-20d,20d)] ;t0

  pars = [2d3,3.,0.3,60,300,0.14,0d] ;
  ;define the number of samples
  n_samples = 100000l
  ;define the number of burn in samples
  burn_in = 50000l
  ;run MCMC fitting
  fit = mcmc_fit(x, y, pars, "nonlinear_model_v3", priors = priors, burn_in = burn_in,$ 
    n_samples = n_samples, samples = samples, ppd_samples=ppd_samples, credible_intervals = credible_intervals)
  print,'best fitted params are:',pars;output pars that best fit the data with the given function
  print,'R',credible_intervals[0,*]
  print,'roi',credible_intervals[1,*]
  print,'C1',credible_intervals[2,*]
  print,'Vi',credible_intervals[3,*]
  print,'period_k',credible_intervals[4,*]
  print,'eta',credible_intervals[5,*]
  print,'t0',credible_intervals[6,*]
 
  wdef,1,1400,700
  title=['R [km]','roi','C1','Vi [km/s]','P_k','eta','t0'] ;'Delta_m'
  !p.multi=[0,4,2]
  for k=0,6 do begin
    cgHistoplot,samples[k,*],title=title[k],charsize=3
    cgoplot,[credible_intervals[k,0],credible_intervals[k,0]],[0,100000],color='green'
    cgoplot,[credible_intervals[k,1],credible_intervals[k,1]],[0,100000],color='green'
    cgoplot,[pars[k],pars[k]],[0,100000],color='red'
  endfor
  !p.multi=0

  evidence = mcmc_fit_evidence(samples, x, y, priors, 'nonlinear_model_v3', n_iterations = 1d5) 
  print, 'evidence of nonlinear model:', evidence


end

pro mcmc_nonlinear
  ;This procedure fit the data wih 9-parameter nonlinear function
  
  ;data ;to sort it in ascending order  
  x = findgen(200)*5
  yerr = randomn(1,n_elements(x)) ;default is uniform distribution=[0,1]
  y = nonlinear_model(x,[2000,3.0,0.3,60,300,280,-0.85,0.24,0d])+yerr
; y = nonlinear_model(x,[2500,3.0,0.3,60,180,150,-0.85,0.24,-20d]) + yerr ;  longer time series, 6 cycles

  cgplot,x,y,psym=1

  ;define priors

  priors = [prior_uniform(1d2,10d3),$; loop radius in km
    prior_uniform(1d,10d),$ ;roi, dimensionless
    prior_uniform(0.1d,0.5d),$ ;C1, dimensionless
    prior_uniform(30d,150d),$ ;Vi, in km/s
    prior_uniform(100d,max(x)/2d),$ ; kink period
    prior_uniform(100d,300d),$ ; period of alfven wave, make it smaller than P_k
    prior_uniform(-3d,0d),$ ;delta_m
    prior_uniform(0d,3d),$ ;;thre for M_bar;max=total(model_denisty*dl)
    prior_uniform(-10d,10d)] ;t0  
  
  ;define the initial guess
  pars = [2000,3.0,0.3,60,300,280,-0.85,0.24,0d] ; for synthetic data

  ;define the number of samples
  n_samples = 500000l 
  ;define the number of burn in samples
  burn_in = 100000l
  ;run MCMC fitting
  fit = mcmc_fit(x, y, pars, "nonlinear_model", priors = priors, burn_in = burn_in, $ ;errors=yerr,
    n_samples = n_samples, samples = samples, ppd_samples=ppd_samples, credible_intervals = credible_intervals)
  print,'best fitted params are:',pars;output pars that best fit the data with the given function

  evidence = mcmc_fit_evidence(samples, x, y, priors, 'nonlinear_model', n_iterations = 1d5,errors=yerr) ;errors=yerr
  print, 'evidence of nonlinear model:', evidence
  
  wdef,0,700,400
  cgoplot,x,y,psym=16,symsize=0.6,color='blue'
  cgoplot,x,fit,color='red'
    
  wdef,1,2200,1000
  title=['R [km]','ro_i','C1','Vi [km/s]','P_k','P_A','Delta_m','zeta','t0']
  !p.multi=[0,5,2]
  for k=0,8 do begin
     cgHistoplot,samples[k,*],title=title[k],charsize=3
     cgoplot,[credible_intervals[k,0],credible_intervals[k,0]],[0,100000],color='green'
     cgoplot,[credible_intervals[k,1],credible_intervals[k,1]],[0,100000],color='green'
     cgoplot,[pars[k],pars[k]],[0,100000],color='red'
  endfor
  !p.multi=0

end
