pro plot_fressin

  common DATA, period, radius, period_boundary, radius_boundary, rate, rate_unc

; Fressin frequencies - for now limited to one planet per star
  
  period_boundary = [0.8, 2.0, 3.4, 5.9, 10.0, 17.0, 29.0, 50.0, 85.0, 145.0, 245.0, 418.0]
  radius_boundary = [0.8, 1.25, 2.0, 4.0, 6.0, 15.0] ; Fressin gives 22 as upper limit but we lower it here to 15
  planet_type = ['Earths', 'Super-Earths', 'Small Neptunes', 'Large Neptunes', 'Giants']
  
  rate_fressin = dblarr(11,5)   ; period bin, radius bin
  rate_fressin[0,*] = [0.18, 0.17, 0.035, 0.004, 0.015]
  rate_fressin[1,*] = [0.61, 0.74, 0.18,  0.006, 0.067]
  rate_fressin[2,*] = [1.72, 1.49, 0.73,  0.11,  0.17]
  rate_fressin[3,*] = [2.70, 2.90, 1.93,  0.091, 0.18]
  rate_fressin[4,*] = [2.70, 4.30, 3.67,  0.29,  0.27]
  rate_fressin[5,*] = [2.93, 4.49, 5.29,  0.32,  0.23]
  rate_fressin[6,*] = [4.08, 5.29, 6.45,  0.49,  0.35]
  rate_fressin[7,*] = [3.46, 3.66, 5.25,  0.66,  0.71]
  rate_fressin[8,*] = [0.0,  6.54, 4.31,  0.43,  1.25]
  rate_fressin[9,*] = [0.0,  0.0,  3.09,  0.53,  0.94]
  rate_fressin[10,*] =[0.0,  0.0,  0.0,   0.24,  1.05]
  rate_fressin = rate_fressin/100.

  unc_fressin = dblarr(11,5)   ; period bin, radius bin
  unc_fressin[0,*] = [0.04, 0.03, 0.011, 0.003, 0.007]
  unc_fressin[1,*] = [0.15, 0.13, 0.03, 0.006, 0.018]
  unc_fressin[2,*] = [0.43, 0.23, 0.09, 0.03, 0.03]
  unc_fressin[3,*] = [0.60, 0.56, 0.19, 0.03, 0.04]
  unc_fressin[4,*] = [0.83, 0.73, 0.39, 0.07, 0.06]
  unc_fressin[5,*] = [1.05, 1.00, 0.64, 0.08, 0.06]
  unc_fressin[6,*] = [1.88, 1.48, 1.01, 0.12, 0.10]
  unc_fressin[7,*] = [2.81, 1.21, 1.05, 0.16, 0.17]
  unc_fressin[8,*] = [0.0,  2.2, 1.03, 0.17, 0.29]
  unc_fressin[9,*] = [0.0, 0.0, 0.9,  0.21, 0.28]
  unc_fressin[10,*] = [0.0, 0.0, 0.0, 0.15, 0.30]
  unc_fressin = unc_fressin/100.

  q = where(unc_fressin lt 0.1*rate_fressin)
  if (q[0] ne -1) then unc_fressin[q] = 0.1*rate_fressin[q]
  q = where(unc_fressin eq 0)
  if (q[0] ne -1) then unc_fressin[q] = 0.1

; set up the arrays we want to fit

  rate = rate_fressin[*,0:2]
  rate_unc = unc_fressin[*,0:2]
  radius = [1.0, 1.63, 3.0]

  per = dblarr(11)
  for i=0,10 do begin
     per[i] = sqrt(period_boundary[i+1]*period_boundary[i])
  endfor

;;;;;; fit [deprecated once the best fit was established]

;  p_guess = [0.305, 1.029, 0.375, -0.0189, 0.0598]
;  p_scale = [0.01, 0.1, 0.1,  0.1, 0.1]
;
;  p = amoeba(1E-5, function_name='chisq', p0=p_guess, scale=p_scale)
;  p = p_guess
;  print, p

;;;;;; plot

  !p.multi=[0,3,2]
  title = ['Earths','Super-Earths','Small Neptunes']
  !p.charsize=3

  c = fsc_color(['Red','Yellow','Cyan'])
  c_calc = fsc_color(['Pink', 'Red','Yellow','Orange','Green', 'Cyan', 'Blue', 'Violet'])
  rx = 1.05

; for finely sampled smooth functions
  dlogp = 0.1
  log_per_fine = -1. + dlogp*dindgen(50)
  per_fine = 10.^log_per_fine
;;;;;;

  for k=0,1 do begin
     for j=0,2 do begin
        
        plot, per, rate_fressin[*,j], /nodata, $
              title=title[j], $
              xtit='period [days]', xra=[0.5,100], xsty=1, /xlog, $
              ytit='n!Dpla!N in Fressin P/R bin', yra=[0.0001,0.08], ysty=1, ylog=k
        
     plotsym,0,/fill
     oploterror, per, rate_fressin[*,j], unc_fressin[*,j], psym=8, /nohat

     rate_calc = 0.*per

     for i=0,n_elements(per)-1 do begin
;        rate_calc[i] = int_planet_occurrence_fit(p, $
        rate_calc[i] = int_planet_occurrence( $
                       period_boundary[i], $
                       period_boundary[i+1], $
                       radius_boundary[j], $
                       radius_boundary[j+1])
     endfor
     print, rate_calc
     plotsym,8,/fill
     oplot, per, rate_calc, psym=8, color=c[j]

     rate_fine = 0.*per_fine
     
     for i=0,n_elements(per_fine)-1 do begin
        p1 = 10.0^( log_per_fine[i] - 0.5*dlogp )
        p2 = 10.0^( log_per_fine[i] + 0.5*dlogp )
        rate_fine[i] = int_planet_occurrence(p1, p2, radius_boundary[j], radius_boundary[j+1]) * 2.3
     endfor

     oplot, per_fine, rate_fine, color=c[j]

  endfor
  endfor

;;;;;; now plot smooth functions

;  dlogp = 0.1
;  log_per_calc = -1. + dlogp*dindgen(50)
;  per_calc = 10.^log_per_calc
;  rad_vals = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
;
;  plot, per_calc, 0.*per_calc, /nodata, $
;        xra=[0.3,100], xsty=1, /xlog, $
;        yra=[0.0001,0.08], ysty=1, /ylog
;
;  for j=0,n_elements(rad_vals)-1 do begin
;     rad = rad_vals[j]
;     rate_calc = 0.*per_calc
;     
;     for i=0,n_elements(per_calc)-1 do begin
;        p1 = 10.0^( log_per_calc[i] - 0.5*dlogp )
;        p2 = 10.0^( log_per_calc[i] + 0.5*dlogp )
;        rate_calc[i] = int_planet_occurrence(p, p1, p2, rad_vals[j]-0.1, rad_vals[j]+0.1) * 2.3
;     endfor
;     
;     oplot, per_calc, rate_calc, color=c_calc[j]
;     
;  endfor
;  
end

function planet_occurrence_fit, logp, r
; gives you dN/(dr*dlogp) as a function of logp and r
; p is in days, r in earth radii
; output is number per factor-of-ten in p and unit interval in r

  common PARAMETERS, p, r1, r2

  n_calc = 0.*r

  for i=0,n_elements(r)-1 do begin
     c = p[0] + p[3]*r[i]^2
     logp0 = p[1] + p[4]*r[i]^2
     dlogp = p[2]
     if (logp le logp0) then n_calc[i] = c*exp( -0.5*( (logp-logp0)/dlogp)^2. ) else n_calc[i]=c
  endfor

  return, n_calc

end

function int_planet_occurrence_fit, params, p1, p2, r1, r2

  common PARAMETERS, p, radius1, radius2

  p = params
  radius1 = r1
  radius2 = r2

  ans = int_2d('planet_occurrence_fit', [alog10(p1), alog10(p2)], 'radius_limits_fit', 48, order=0)

  return, ans

end

function radius_limits_fit, logp

  common PARAMETERS, p, radius1, radius2

  return, [radius1, radius2]

end

function chisq, p

  common DATA, per, period_boundary, radius_boundary, rate, rate_unc

  rate_calc = 0.*per

  chisq = 0.0

  for j=0,2 do begin
     
     r1 = radius_boundary[j]
     r2 = radius_boundary[j+1]

;     if (j eq 1) then weight = 5.0+ 0.*per else weight = 1.0 + 0.*per
     weight = 1.0 + 0.*per
     weight[8:10] = 0.0

     for i=0,n_elements(per)-1 do begin
        p1 = period_boundary[i]
        p2 = period_boundary[i+1]
        rate_calc[i] = int_planet_occurrence_fit(p, p1, p2, r1, r2)
     endfor
     
     dr = (rate[*,j]-rate_calc)/(rate_unc[*,j])
     chisq = chisq + total( weight*dr^2. )
     
  endfor

  return, chisq

end
