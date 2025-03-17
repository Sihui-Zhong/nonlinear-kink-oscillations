;SoBAT MCMC fitting with nonlinear damping function by turbulence

function nonlinear_model,t,pars,_extra=_extra
  compile_opt idl2

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
  threshold = pars[7] ;0.9
  t0 = pars[8]
  
  t1 = t-t0
  M_bar = threshold*(Vi - dV*sqrt(roe)/(sqrt(roi)+sqrt(roe))) ;0.24
  
  E = 2*R*H*roi ;mass
  A = E*Vi ;momentum
  B = 1/(2*H)*C2/sqrt(2)*sqrt(roe)/(sqrt(roi)+sqrt(roe))*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV
  D = M_bar*C2/sqrt(2)*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV*R/wD
  G = delta_m/(2*H*roi)*C2/sqrt(2)*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV
  A1 = A/(E*(G^2))*(-B*G*sin(wk*t1)/wk+(B+G)*(cos(wk/G)*(sici(wk*(t1+1/G)))[1,*]+sin(wk/G)*(sici(wk*(t1+1/G)))[0,*]))
  A2 = D/(2*E*G)*(-cos(wk/G)*(sici(wk*(t1+1/G)))[0,*]+sin(wk/G)*(sici(wk*(t1+1/G)))[1,*]+$
    cos(wA/G)*(sici(wA*(t1+1/G)))[0,*]-sin(wA/G)*(sici(wA*(t1+1/G)))[1,*])  
  func = A1+A2-mean(A1+A2) ;A1+A2-A1[0]-A2[0]
  return,func
  
end

pro mcmc_nonlinear
  ;data ;remember to sort it in ascending order
  
  x = findgen(200)*5
  yerr = randomn(1,n_elements(x)) ;default is uniform distribution=[0,1]
  y = nonlinear_model(x,[2000,3.0,0.3,60,300,280,-0.85,0.24,0d])+yerr
; y = nonlinear_model(x,[2500,3.0,0.3,60,180,150,-0.85,0.24,-20d]) + yerr ;  longer time series, 6 cycles
  yerr = abs(yerr)

  cgplot,x,y,psym=1
  ;in observation order  
  ;x2 = findgen(1140)
  ;y2 = nonlinear_model(x2,[1*1d3,3.0,0.3,0.12*215,323,288,-0.85,0.24,1d])

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
  title=['R [km]','ro_i','C1','Vi [km/s]','P_k','P_A','Delta_m','thre','t0']
  !p.multi=[0,5,2]
  for k=0,8 do begin
     cgHistoplot,samples[k,*],title=title[k],charsize=3
     cgoplot,[credible_intervals[k,0],credible_intervals[k,0]],[0,100000],color='green'
     cgoplot,[credible_intervals[k,1],credible_intervals[k,1]],[0,100000],color='green'
     cgoplot,[pars[k],pars[k]],[0,100000],color='red'
  endfor
  !p.multi=0

end
