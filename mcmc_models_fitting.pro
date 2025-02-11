;SoBAT MCMC fitting with exponential damping function and Hillier+2024 formula (nonlinear)
function model_exp_decay,x,a,n_trend=n_trend
  t_start = a[0] ;oscillation start time
  period = a[1] ;period
  ;q_factor = a[2]
  tau = a[2] ;q_factor*period ;decay time
  amp = a[3] ;initial amplitude
  displ = a[4]
  ref_y = a[5:5+n_trend-1]
  ref_x = findgen(n_trend)/n_trend*abs(x[0]-x[-1])+x[0]
  tosc = x-t_start
  omega = 1.d*!dpi/period
  phi = asin((displ)) ;initial phase
  ;decay profile
  damp = amp*exp(-(tosc/tau)^1)*(x ge t_start) ;
  oscillation = damp * sin(omega*(tosc>0d)+phi)
  trend = spline(ref_x,ref_y,x)
  return, trend + oscillation
end

function exp_model, x, pars, _extra=_extra
  compile_opt idl2
  a = pars[0]
  tau = pars[1]
  p = pars[2]
  phi = pars[3]
  t0 = pars[4]
  t = x-t0
  return, a*exp(-t/tau)*sin(2*!dpi*t/p+phi)
end

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
  s = pars[8] ;scale factor
  t0 = pars[9]
  
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
  return,func*s;+bg
  
end

pro mcmc_exp
  ;data ;remember to sort it in ascending order
;  x = findgen(200)*0.5
;  noise = randomn(1,n_elements(x))
;  yerr = noise/stdev(noise)
;  y = 10*exp(-x/80)*cos(!dpi*2*x/20.0+!dpi/4)+yerr
;  save,filename='exp_decay_data.sav',x,y,yerr,/verbose
;  restore,'exp_decay_data.sav',/v
;  yerr = abs(yerr)
   ;restore,'nl_data_long_0percentnoise.sav',/v
;   restore,'nl_data_30percentnoise.sav',/v
   
;  restore,'/Users/sihuizhong/Documents/Exeter/KHI/20170907_td_slit19_GaussFit.sav',/v
;  x = (t_loop-t_loop[0])*12.0
;  centre=smooth(centre,3,/edge_tr)
;;  smoothy = smooth(centre,130,/edge_tr)
;;  result = poly_fit(x,centre,2,measure_errors=centre_err,sigma=sigma,yfit=trend)
;  trend = smooth(centre,130,/edge_tr)
;  y = centre-trend
;  yerr = centre_err
  restore,'/home/sihuizhong/Documents/KHI/20120530/20120530_td_GaussFit.sav',/v
  x = (t_loop-t_loop[0])*12.0
  ;result = poly_fit(x,centre,3,measure_errors=centre_err,sigma=sigma,yfit=trend)
  trend = smooth(centre,70,/edge_tr)
  cgplot,x,centre
  cgoplot,x,trend,color='red'
  y = centre - trend
  yerr = centre_err
  
;  x = findgen(200)*5
;  yerr = 10*randomu(seed,n_elements(x),/normal) ;10% noise
;  y = nonlinear_model(x,[2000,3.0,0.3,60,300,280,-0.85,0.24,1d,0d])+yerr
;  yerr=abs(yerr)
;   restore,'data_5percentnoise.sav',/v
  ;x = findgen(240)/2.0
  ;y = nonlinear_model(x,[0.5,3.0,0.3,0.12,0.18,0.2,-0.85,0.24,1d])+randomu(seed,n_elements(x))
  ;y -=mean(y)
  ;err = randomu(seed,n_elements(x))
  ;define priors
  priors = objarr(5) ;distribution of params
  priors[0] = prior_uniform(max(y)*0.1d,max(y)*2d);max(y), amplitude
  priors[1] = prior_uniform(100d,max(x)*2d);tau
  priors[2] = prior_uniform(100d,max(x)/2.);period--intial guess = 0d ;1d-5,max(x)
  priors[3] = prior_uniform(!dpi*(-2d),!dpi*2d) ;phase
  priors[4] = prior_uniform(-200d,200d)
  ;priors[4] = prior_uniform(-max(x)/3d,200d) ;t0
  ;define the initial guess
  ;ppars = [1900, 1000, 180, 0d, 20d]
  ;ppars = [max(y)*1.5, max(x)/2.0, max(x)/4.0, 0d, 0d];
  ;ppars = [max(y), 1100, 300, 0d,0d];
  ppars = [max(y), 1100, 300, 0d,0d];20120530
  ;define the number of samples
  n_samples = 500000l;5000000l
  ;define the number of burn in samples
  burn_in = 100000l
  ;run MCMC fitting
  fit = mcmc_fit(x, y, ppars, "exp_model", priors = priors, burn_in = burn_in,$ ;errors=yerr
    n_samples = n_samples, samples = samples, ppd_samples=ppd_samples, credible_intervals = credible_intervals)
  print,'best fitted params are:',ppars;output pars that best fit the data with the given function
  ;print,credible_intervals
  print,'amplitude',credible_intervals[0,*]
  print,'tau',credible_intervals[1,*]
  print,'period',credible_intervals[2,*]
  print,'phase',credible_intervals[3,*]
  print,'t0',credible_intervals[4,*]

  ;calculate evidence
  evidence = mcmc_fit_evidence(samples, x, y, priors, 'exp_model', n_iterations = 1d5) ;errors=yerr
  print, 'evidence of exp model:', evidence 

  wdef,0,700,400
  nbins=500
  ;ppd=mcmc_ppd_histogram_add(x,y,ppd_samples, nbins = nbins, xnbins = 200,hist_x = hist_x, hist_y = hist_y) & help,ppd
  ppd=mcmc_ppd_histogram(x,y,ppd_samples, nbins = nbins, hist_x = hist_x, hist_y = hist_y)
  range = minmax(hist_y) & yrange=range[1]-range[0]
  range = minmax(hist_x) & xrange=range[1]-range[0]
  help,ppd
  ;plot_image,ppd
  ;plot_image,ppd,origin=[min(hist_x),min(hist_y)],scale=[xrange/n_elements(hist_x),yrange/nbins],position=[0.1,0.1,0.9,0.9],/nor
  loadct,0
  contour,sigrange(1.- ppd/max(ppd),frac=0.999), hist_x, hist_y,/fill, nlevels = 255,/xsty,/ysty,charsize=2,$
    xtitle='!5 time [s]', ytitle='!5 Displacement [px]!X',color=255, xrange = [0,max(x)],yrange=[min(y)*1.2,max(y)*1.2];,/xlog,/ylog
  ;wdef,5,500,400
  cgoplot,x,y,psym=16,color='blue'
  cgoplot,x,fit,color='red'
  ;l=findgen(150)/100.
  ;cgoplot,l,pars[1]*l^pars[0]+pars[2],color='red',thick=2.

  wdef,1,1300,400
  title = ['amplitude','tau','period','phase','t0']
  !p.multi=[0,5,1]
  for k=0,4 do begin
     cgHistoplot,samples[k,*],title=title[k],charsize=2 ;binsize=0.5
     cgoplot,[credible_intervals[k,0],credible_intervals[k,0]],[0,100000],color='green'
     cgoplot,[credible_intervals[k,1],credible_intervals[k,1]],[0,100000],color='green'
  endfor


  !p.multi=0
  ;save,filename='20120530_mcmc_exp_fit.sav',x,y,yerr,fit, ppars,evidnece,samples,credible_intevals,ppd,hist_x,hist_y,priors,/verbose
end

pro mcmc_nonlinear
  ;data ;remember to sort it in ascending order
;  restore,'exp_decay_data.sav',/v
;  yerr = abs(yerr)
  
;  x = findgen(200)*5
;  yerr = randomn(1,n_elements(x)) ;default is uniform distribution=[0,1]
;  ;y = nonlinear_model(x,[2500,3.0,0.3,60,200,180,-0.85,0.24,1d])+yerr
  ;y = nonlinear_model(x,[2000,3.0,0.3,60,300,280,-0.85,0.24,1d,0d])+yerr
;  y = nonlinear_model(x,[2500,3.0,0.3,60,180,150,-0.85,0.24,1d,-20d]) + yerr ;  longer time series, 6 cycles
;  yerr = abs(yerr)
;  save,filename='nl_data_long_0percentnoise.sav',x,y,yerr,/verbose
;  restore,'nl_data_30percentnoise.sav',/v

   restore,'/home/sihuizhong/Documents/KHI/20120530/20120530_td_GaussFit.sav',/v
   x = (t_loop-t_loop[0])*12.0
   ;result = poly_fit(x,centre,3,measure_errors=centre_err,sigma=sigma,yfit=trend)
   centre[22]=centre[22]-1
   centre[23:28]=centre[23:28]-2
   trend = smooth(centre,70,/edge_tr)
   cgplot,x,centre
   cgoplot,x,trend,color='red'
   y = centre - trend
   yerr = centre_err
  cgplot,x,y,psym=1
  ;in observation order  
  ;x2 = findgen(1140);330*4
  ;y2 = nonlinear_model(x2,[1*1d3,3.0,0.3,0.12*215,323,288,-0.85,0.24,1d])
;  restore,'20170907_td_slit19_GaussFit.sav',/v
;  x = (t_loop-t_loop[0])*12.0
;  centre=smooth(centre,3,/edge_tr)
;  trend = smooth(centre,130,/edge_tr)
;  y = (centre-trend) ;*440
;  yerr = centre_err
  ;detrending polynomial
  ;define priors
  priors = [prior_uniform(1d3,5d3),$; loop radius in km ;for 201205030
    prior_uniform(2d,9d),$ ;roi, dimensionless Pk/PA=1.12
    prior_uniform(0.2d,0.6d),$ ;C1, dimensionless
    prior_uniform(10d,100d),$ ;Vi, in km/s
    prior_uniform(200d,350d),$ ; kink period
    prior_uniform(180d,280d),$ ; period of alfven wave, make it smaller than P_k
    prior_uniform(-1d,0.1d),$ ;delta_m
    prior_uniform(0d,5d),$ ;;thre for M_bar;max=total(model_denisty*dl)
    prior_uniform(0d,0.1d),$ ;scale factor for overall amplitude
    prior_uniform(-100d,100d)] ;t0
;  priors = [prior_uniform(1d2,10d3),$; loop radius in km
;    prior_uniform(1d,10d),$ ;roi, dimensionless
;    prior_uniform(0.1d,0.5d),$ ;C1, dimensionless
;    prior_uniform(30d,150d),$ ;Vi, in km/s
;    prior_uniform(250d,max(x)/2d),$ ; kink period
;    prior_uniform(220d,300d),$ ; period of alfven wave, make it smaller than P_k
;    prior_uniform(-3d,0d),$ ;delta_m
;    prior_uniform(0d,3d),$ ;;thre for M_bar;max=total(model_denisty*dl)
;    prior_uniform(0d,2d),$ ;scale factor for overall amplitude
;    prior_uniform(-10d,10d)] ;t0
;  priors = [prior_uniform(1d3,5d3),$; loop radius in km
;    prior_uniform(2d,4d),$ ;roi, dimensionless
;    prior_uniform(0.1d,0.5d),$ ;C1, dimensionless
;    prior_uniform(30d,90d),$ ;Vi, in km/s
;    prior_uniform(100d,max(x)/2d),$ ; kink period
;    prior_uniform(100d,180d),$ ; period of alfven wave, make it smaller than P_k
;    prior_uniform(-1d,0d),$ ;delta_m
;    prior_uniform(0d,1d),$ ;;thre for M_bar;max=total(model_denisty*dl)
;    prior_uniform(0.5d,1.5d),$ ;scale factor for overall amplitude
;    prior_uniform(-50d,0d)] ;t0
  
  ;prior for exp_decay data
;  priors = [prior_uniform(10d,3d3),$; loop radius in km
;    prior_uniform(1d,5d),$ ;roi, dimensionless
;    prior_uniform(0.1d,0.5d),$ ;C1, dimensionless
;    prior_uniform(1d,100d),$ ;Vi, in km/s
;    prior_uniform(5d,max(x)/2d),$ ; kink period
;    prior_uniform(5d,30d),$ ; period of alfven wave
;    prior_uniform(-1d,0d),$ ;delta_m
;    prior_uniform(0d,0.3d),$ ;;thre for M_bar;max=total(model_denisty*dl)
;    prior_uniform(0d,2d),$ ;scale factor for overall amplitude
;    prior_uniform(-20d,20d)] ;t0
  ;define the initial guess
  ;pars = [0.5,3.0,0.3,0.15,2*!pi/0.18d,2*!pi/0.2d,-0.85,0.24,1d];
  ;pars = [1d3,3.0,0.3,0.11*215,331,290,-0.85,0.24,1d];
  ;pars = [2*1d3,4.0,0.3,30,480,450,-0.85,0.24,10d-3,0d]; for 20170907 event
  pars = [2d3,4.0,0.3,20,250,220,-0.85,0.24,10d-3,0d]; for 20120530 event
  ;pars = [2000,3.0,0.3,60,300,280,-0.85,0.24,1d,0d] ; for synthetic data
  ;pars = [2500,3.0,0.3,60,180,150,-0.85,0.24,1d,-20d] ; for synthetic data_long
  ;pars = [500,5.0,0.3,60,20,15,-0.85,0.24,3d-2,-7d] ; for exp_decay data
  ;define the number of samples
  n_samples = 500000l 
  ;define the number of burn in samples
  burn_in = 100000l
  ;run MCMC fitting
  fit = mcmc_fit(x, y, pars, "nonlinear_model", priors = priors, burn_in = burn_in, $ ;errors=yerr,
    n_samples = n_samples, samples = samples, ppd_samples=ppd_samples, credible_intervals = credible_intervals)
  print,'best fitted params are:',pars;output pars that best fit the data with the given function
  print,'R',credible_intervals[0,*]
  print,'roi',credible_intervals[1,*]
  print,'C1',credible_intervals[2,*]
  print,'Vi',credible_intervals[3,*]
  print,'period_k',credible_intervals[4,*]
  print,'period_A',credible_intervals[5,*]
  print,'delta_m',credible_intervals[6,*]
  print,'thre',credible_intervals[7,*]
  print,'scale',credible_intervals[8,*]
  print,'t0',credible_intervals[9,*]
  
  ;stop
  evidence = mcmc_fit_evidence(samples, x, y, priors, 'nonlinear_model', n_iterations = 1d5,errors=yerr) ;errors=yerr
  print, 'evidence of nonlinear model:', evidence
  
  wdef,0,700,400
  nbins=500
  ;ppd=mcmc_ppd_histogram_add(x,y,ppd_samples, nbins = nbins, xnbins = 200,hist_x = hist_x, hist_y = hist_y) & help,ppd
  ppd=mcmc_ppd_histogram(x,y,ppd_samples, nbins = nbins, hist_x = hist_x, hist_y = hist_y)
  range = minmax(hist_y) & yrange=range[1]-range[0]
  range = minmax(hist_x) & xrange=range[1]-range[0]
  ;help,ppd
  ;plot_image,ppd
  ;plot_image,ppd,origin=[min(hist_x),min(hist_y)],scale=[xrange/n_elements(hist_x),yrange/nbins],position=[0.1,0.1,0.9,0.9],/nor
  loadct,0
  contour,sigrange(1.- ppd/max(ppd),frac=0.999), hist_x, hist_y,/fill, nlevels = 255,/xsty,/ysty,charsize=2,$
    xtitle='!X Time [s]', ytitle='Displacement ',color=255, xrange = [0,max(x)],yrange=[min(y)*1.2,max(y)*1.2];,/xlog,/ylog
  ;wdef,5,500,400
  cgoplot,x,y,psym=16,symsize=0.6,color='blue'
  cgoplot,x,fit,color='red'
  ;write_png,'fit_long_34.png',tvrd(/true)
  ;write_png,'fit_00322.png',tvrd(/true)
  write_png,'mcmcfit_corrected.png',tvrd(/true)
  
  ;save,filename='synthetic_mcmc_fit_00322.sav',x,y,yerr,fit,pars,evidence,samples,credible_intervals,hist_x,hist_y,ppd,priors,/verbose
  ;save,filename='20170907_mcmc_fit.sav',x,y,yerr,fit, pars,evidence,samples,credible_intervals,ppd,hist_x,hist_y,priors,/verbose
  save,filename='20120530_mcmc_nl_fit_corrected.sav',x,y,yerr,fit, pars,evidence,samples,credible_intervals,ppd,hist_x,hist_y,priors,/verbose
  
  wdef,1,2200,1000
  title=['R [km]','ro_i','C1','Vi [km/s]','P_k','P_A','Delta_m','thre','scale','t0']
  !p.multi=[0,5,2]
  for k=0,9 do begin
     cgHistoplot,samples[k,*],title=title[k],charsize=3
     cgoplot,[credible_intervals[k,0],credible_intervals[k,0]],[0,100000],color='green'
     cgoplot,[credible_intervals[k,1],credible_intervals[k,1]],[0,100000],color='green'
     cgoplot,[pars[k],pars[k]],[0,100000],color='red'
  endfor
  !p.multi=0
  write_png,'params_samples_nl_corrected.png',tvrd(/true)
  stop
end
