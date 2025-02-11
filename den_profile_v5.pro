PRO den_profile_v5,rho_1,rho_2,velocity,y,density
;Calculate the spatial-averaged density and velocity profile in the mixing layer
;Author: Andrew Hillier (a.s.hillier@exeter.ac.uk)
;rho_1=1.
;rho_2=10.0
print,rho_2

alpha_1=rho_1/(rho_1+rho_2)
alpha_2=rho_2/(rho_1+rho_2)


params_array=fltarr(1,4)
the_array=fltarr(4,4)
y=findgen(10001)/5000.-1.+(sqrt(alpha_1)-sqrt(alpha_2))/(sqrt(alpha_1)+sqrt(alpha_2))

y_1=-2.*sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2))
y_2=2.*sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2))


KE_min=1000.0
m=0
l=0

;density -- 3-order polynomials
;velocity -- 5-order polynomials
for k=5000,9999 do begin ;vx=0 is at the second hafl
y_match1=k


the_array(0,0)=1.
the_array(1,0)=-2.*sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2)) ;y_1
the_array(2,0)=4.*alpha_2/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array(3,0)=-8.*alpha_2^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array(0,1)=1.
the_array(1,1)=2.*sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2)) ;y_2
the_array(2,1)=4.*alpha_1/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array(3,1)=8.*alpha_1^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array(0,2)=2. 
the_array(1,2)=2.*(sqrt(alpha_1)-sqrt(alpha_2))/(sqrt(alpha_1)+sqrt(alpha_2)) ;y_2^2 - y_1^2
the_array(2,2)=8./3.*((alpha_1)^1.5+(alpha_2)^1.5)/(sqrt(alpha_1)+sqrt(alpha_2))^3 ;y_2^3 - y_1^3
the_array(3,2)=4.*(sqrt(alpha_1)^4.-sqrt(alpha_2)^4.)/(sqrt(alpha_1)+sqrt(alpha_2))^4 ;;y_2^4 - y_1^4
the_array(0,3)=1.
the_array(1,3)=y[y_match1]
the_array(2,3)=y[y_match1]^2
the_array(3,3)=y[y_match1]^3


invert_array=invert(the_array,/double)

const_matrix=fltarr(1,4)

const_matrix=[[rho_1],[rho_2],[2.*sqrt(rho_1*rho_2)],[sqrt(rho_1*rho_2)]]

params_array=const_matrix # invert_array

density=params_array[0]+params_array[1]*y+params_array[2]*y^2+params_array[3]*y^3


if (min(3.*params_array[3]*y^2+2.*params_array[2]*y+params_array[1]) ge 0.0 ) then begin

const_matrix1=fltarr(1,6)

the_array1=fltarr(6,6)

the_array1(0,0)=1
the_array1(1,0)=-2.*sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2))
the_array1(2,0)=4.*alpha_2/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array1(3,0)=-8.*alpha_2^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array1(4,0)=16.*alpha_2^2/(sqrt(alpha_1)+sqrt(alpha_2))^4
the_array1(5,0)=(-2.*sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2)))^5

the_array1(0,1)=1
the_array1(1,1)=2.*sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2))
the_array1(2,1)=4.*alpha_1/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array1(3,1)=8.*alpha_1^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array1(4,1)=16.*alpha_1^2/(sqrt(alpha_1)+sqrt(alpha_2))^4
the_array1(5,1)=(2.*sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2)))^5

the_array1(0,2)=params_array[0]*(y_2-y_1)+0.5*params_array[1]*(y_2^2-y_1^2)+params_array[2]/3.*(y_2^3-y_1^3) $
     +params_array[3]/4.*(y_2^4-y_1^4) 
the_array1(1,2)=params_array[0]/2.*(y_2^2.-y_1^2.)+params_array[1]/3.*(y_2^3-y_1^3)+params_array[2]/4.*(y_2^4-y_1^4) $
     +params_array[3]/5.*(y_2^5-y_1^5)
the_array1(2,2)=params_array[0]/3.*(y_2^3.-y_1^3.)+params_array[1]/4.*(y_2^4-y_1^4)+params_array[2]/5.*(y_2^5-y_1^5) $
     +params_array[3]/6.*(y_2^6-y_1^6)
the_array1(3,2)=params_array[0]/4.*(y_2^4.-y_1^4.)+params_array[1]/5.*(y_2^5-y_1^5)+params_array[2]/6.*(y_2^6-y_1^6) $
     +params_array[3]/7.*(y_2^7-y_1^7)
the_array1(4,2)=params_array[0]/5.*(y_2^5.-y_1^5.)+params_array[1]/6.*(y_2^6-y_1^6)+params_array[2]/7.*(y_2^7-y_1^7) $
     +params_array[3]/8.*(y_2^8-y_1^8)
the_array1(5,2)=params_array[0]/6.*(y_2^6.-y_1^6.)+params_array[1]/7.*(y_2^7-y_1^7)+params_array[2]/8.*(y_2^8-y_1^8) $
     +params_array[3]/9.*(y_2^9-y_1^9)     

the_array1(0,3)=0.
the_array1(1,3)=1.
the_array1(2,3)=-4.*sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2))
the_array1(3,3)=12.*alpha_2/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array1(4,3)=-32.*alpha_2^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array1(5,3)=5.*(-2.*sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2)))^4

the_array1(0,4)=0.0
the_array1(1,4)=1.
the_array1(2,4)=4.*sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2))
the_array1(3,4)=12.*alpha_1/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array1(4,4)=32.*alpha_1^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array1(4,4)=5.*(2.*sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2)))^4

the_array1(0,5)=1.
the_array1(1,5)=y[y_match1]
the_array1(2,5)=y[y_match1]^2
the_array1(3,5)=y[y_match1]^3
the_array1(4,5)=y[y_match1]^4
the_array1(5,5)=y[y_match1]^5



const_matrix1=[[sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2))],[-sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2))],[0.0],[0.0],[0.0],[0.0]]

invert_array1=invert(the_array1,/double)

params_array1=const_matrix1 # invert_array1

velocity=params_array1[0]+params_array1[1]*y+params_array1[2]*y^2+params_array1[3]*y^3+params_array1[4]*y^4+params_array1[5]*y^5
;wait,0.2
;plot,velocity

KE_1=total(0.5*velocity^2*density)*0.002

y_grad=max(where(y lt -params_array1[2]/(3.*params_array1[3])))

;if ((KE_1 lt KE_min) and $
;     (3.*params_array1[3]*y[y_grad]^2+2.*params_array1[2]*y[y_grad]+params_array1[1] eq min(3.*params_array1[3]*y^2+2.*params_array1[2]*y;+params_array1[1])) and (max(3.*params_array1[3]*y^2+2.*params_array1[2]*y+params_array1[1]) le 0.0)) then begin

;print,max(params_array1[1]+params_array1[2]*y*2.+params_array1[3]*y^2.*3.+params_array1[4]*y^3.*4.+params_array1[5]*y^4*5.)


if ((KE_1 lt KE_min) ) and (max(params_array1[1]+params_array1[2]*y*2.+params_array1[3]*y^2.*3.+params_array1[4]*y^3.*4.+params_array1[5]*y^4*5.) le 10.0^(-6.)) then begin

KE_min=KE_1

l=k
endif

endif
endfor

if (l gt 0) then begin

print,'enter'
y_match=m
y_match1=l

params_array=fltarr(1,4)
the_array=fltarr(4,4)

the_array(0,0)=1.
the_array(1,0)=-2.*sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2))
the_array(2,0)=4.*alpha_2/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array(3,0)=-8.*alpha_2^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array(0,1)=1.
the_array(1,1)=2.*sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2))
the_array(2,1)=4.*alpha_1/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array(3,1)=8.*alpha_1^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array(0,2)=2.
the_array(1,2)=2.*(sqrt(alpha_1)-sqrt(alpha_2))/(sqrt(alpha_1)+sqrt(alpha_2))
the_array(2,2)=8./3.*((alpha_1)^1.5+(alpha_2)^1.5)/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array(3,2)=4.*((alpha_1)^2.-(alpha_2)^2)/(sqrt(alpha_1)+sqrt(alpha_2))^4
the_array(0,3)=1.
the_array(1,3)=y[y_match1]
the_array(2,3)=y[y_match1]^2
the_array(3,3)=y[y_match1]^3

;print,the_array
invert_array=invert(the_array,/double)
;print,'invert A-1',invert_array
const_matrix=[[rho_1],[rho_2],[2.*sqrt(rho_1*rho_2)],[sqrt(rho_1*rho_2)]]
;print,'const_matrix',const_matrix
params_array=const_matrix # invert_array

density=params_array[0]+params_array[1]*y+params_array[2]*y^2+params_array[3]*y^3

;print,'coefficients:',params_array
print,(params_array[1]+params_array[2]*y[0]+2.*params_array[3]*y[0]^2)/rho_1
print,(params_array[1]+params_array[2]*y[10000]+2.*params_array[3]*y[10000]^2)/rho_2

const_matrix1=fltarr(1,6)

the_array1=fltarr(6,6)


the_array1(0,0)=1
the_array1(1,0)=-2.*sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2))
the_array1(2,0)=4.*alpha_2/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array1(3,0)=-8.*alpha_2^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array1(4,0)=16.*alpha_2^2/(sqrt(alpha_1)+sqrt(alpha_2))^4
the_array1(5,0)=(-2.*sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2)))^5

the_array1(0,1)=1
the_array1(1,1)=2.*sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2))
the_array1(2,1)=4.*alpha_1/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array1(3,1)=8.*alpha_1^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array1(4,1)=16.*alpha_1^2/(sqrt(alpha_1)+sqrt(alpha_2))^4
the_array1(5,1)=(2.*sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2)))^5

the_array1(0,2)=params_array[0]*(y_2-y_1)+0.5*params_array[1]*(y_2^2-y_1^2)+params_array[2]/3.*(y_2^3-y_1^3) $
     +params_array[3]/4.*(y_2^4-y_1^4) 
the_array1(1,2)=params_array[0]/2.*(y_2^2.-y_1^2.)+params_array[1]/3.*(y_2^3-y_1^3)+params_array[2]/4.*(y_2^4-y_1^4) $
     +params_array[3]/5.*(y_2^5-y_1^5)
the_array1(2,2)=params_array[0]/3.*(y_2^3.-y_1^3.)+params_array[1]/4.*(y_2^4-y_1^4)+params_array[2]/5.*(y_2^5-y_1^5) $
     +params_array[3]/6.*(y_2^6-y_1^6)
the_array1(3,2)=params_array[0]/4.*(y_2^4.-y_1^4.)+params_array[1]/5.*(y_2^5-y_1^5)+params_array[2]/6.*(y_2^6-y_1^6) $
     +params_array[3]/7.*(y_2^7-y_1^7)
the_array1(4,2)=params_array[0]/5.*(y_2^5.-y_1^5.)+params_array[1]/6.*(y_2^6-y_1^6)+params_array[2]/7.*(y_2^7-y_1^7) $
     +params_array[3]/8.*(y_2^8-y_1^8)
the_array1(5,2)=params_array[0]/6.*(y_2^6.-y_1^6.)+params_array[1]/7.*(y_2^7-y_1^7)+params_array[2]/8.*(y_2^8-y_1^8) $
     +params_array[3]/9.*(y_2^9-y_1^9)     

the_array1(0,3)=0.
the_array1(1,3)=1.
the_array1(2,3)=-4.*sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2))
the_array1(3,3)=12.*alpha_2/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array1(4,3)=-32.*alpha_2^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array1(5,3)=5.*(-2.*sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2)))^4

the_array1(0,4)=0.0
the_array1(1,4)=1.
the_array1(2,4)=4.*sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2))
the_array1(3,4)=12.*alpha_1/(sqrt(alpha_1)+sqrt(alpha_2))^2
the_array1(4,4)=32.*alpha_1^1.5/(sqrt(alpha_1)+sqrt(alpha_2))^3
the_array1(4,4)=5.*(2.*sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2)))^4

the_array1(0,5)=1.
the_array1(1,5)=y[y_match1]
the_array1(2,5)=y[y_match1]^2
the_array1(3,5)=y[y_match1]^3
the_array1(4,5)=y[y_match1]^4
the_array1(5,5)=y[y_match1]^5


const_matrix1=[[sqrt(alpha_2)/(sqrt(alpha_1)+sqrt(alpha_2))],[-sqrt(alpha_1)/(sqrt(alpha_1)+sqrt(alpha_2))],[0.0],[0.0],[0.0],[0.0]]

invert_array1=invert(the_array1,/double)

params_array1=const_matrix1 # invert_array1

velocity=params_array1[0]+params_array1[1]*y+params_array1[2]*y^2+params_array1[3]*y^3+params_array1[4]*y^4+params_array1[5]*y^5
;print,params_array1
;print,max(params_array1[1]+params_array1[2]*y*2.+params_array1[3]*y^2.*3.+params_array1[4]*y^3.*4.+params_array1[5]*y^4*5.)
endif else begin

density=fltarr(10001)
velocity=fltarr(10001)


endelse


;RETURN velocity,y

end
