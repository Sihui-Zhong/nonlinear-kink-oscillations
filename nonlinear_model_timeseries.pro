PRO nonlinear_model_timeseries,Pk,Vi,roi,R,C1=C1,percent=percent,t,Amp
;+
;oscillation amplitude as a function of time, in the case of nonlinear wave damping by KHI-induced turbulence 
;length unit of Vi and R should be same, Vi [km/s], R[km]-->amplitude in km
;unit of t1 is same as period_k
;Author: Sihui Zhong (s.zhong3@exeter.ac.uk, sihui.zhong@outlook.com)
;-
  compile_opt idl2

  H = !dpi*R/4.0
  roe = 1.0 ;external
  if not keyword_set(C1) then C1 = 0.3
  C2 = sqrt(2)*C1
  dV = 1.1*Vi ;delta_V

  wk = 2*!dpi/Pk ;wk, angular frequency, =2*!pi/P
  den_profile_v5,roe,roi,velocity1,y1,density1

  wA = wk*sqrt(roe)*mean(1.0/sqrt(density1))/sqrt(2*roe/(roi+roe))
  wD = (wA-wk)/2.0
  wS = (wA+wk)/2.0

  m_t = total(reverse(density1),/cumulative)*0.0001-max(y1/2.0*roi)
  if not keyword_set(percent) then percent=0.9
  ind = max(where(reverse(density1) ge roi*percent))
  delta_m = m_t[ind]
  ind1 = min(where(density1 ge roi*percent))
  threshold = total(density1[ind1:*]*0.0001)

  M_bar = threshold*(Vi - dV*sqrt(roe)/(sqrt(roi)+sqrt(roe))) ;0.24
  E = 2*R*H*roi ;mass
  A = E*Vi ;momentum
  B = 1/(2*H)*C2/sqrt(2)*sqrt(roe)/(sqrt(roi)+sqrt(roe))*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV
  D = M_bar*C2/sqrt(2)*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV*R/wD
  G = delta_m/(2*H*roi)*C2/sqrt(2)*(roi*roe)^(0.25)/(sqrt(roi)+sqrt(roe))*dV
  
  ;print,1.0/B,-1.0/G
  tc = min([1.0/B,-1/G])
  t = round(findgen(250)/250.*tc)
  if max(t1) le 0.25*Pk then print, 't < 1/4P.'

  A1 = A/(E*(G^2))*(-B*G*sin(wk*t1)/wk+(B+G)*(cos(wk/G)*(sici(wk*(t1+1/G)))[1,*]+sin(wk/G)*(sici(wk*(t1+1/G)))[0,*]))
  A2 = D/(2*E*G)*(-cos(wk/G)*(sici(wk*(t1+1/G)))[0,*]+sin(wk/G)*(sici(wk*(t1+1/G)))[1,*]+$
    cos(wA/G)*(sici(wA*(t1+1/G)))[0,*]-sin(wA/G)*(sici(wA*(t1+1/G)))[1,*])
  
  Amp = A1+A2-A1[0]-A2[0] ;A1+A2-mean(A1+A2) 

end
