PRO fullStarStruct__define
  template_star = {
                  starstruct, $
                  dart: 0,  $    ; from the Dartmouth over-write?
                  r: 0.0, $     ; Radius (solar)
                  ;oldr: 0.0, $     ; Radius (solar)
                  m: 0.0, $     ; Mass (solar)
                  teff: 0.0, $  ; Effective T (Kelvin)
                  ;oldteff: 0.0, $  ; Effective T (Kelvin)
                  cosi: 0.0, $  ; cos inclination of planets
                  ;var: 0.0, $   ; 3 hr. variability (ppt)
                  age: 0.0, $   ; in Gyr
                  mini: 0.0, $
                  feh: 0.0, $
                  logg: 0.0, $
                  mag: {magstruct}, $
                  companion: {compstruct}, $
                  pri: 0, $        ; Primary of binary?
                  sec: 0, $        ; Secondary of binary? 
                  spl: 0, $        ; Split (for triples and quadruples)
                  ffi: 0, $         ; Full-frame only (default is postage stamp)
                  gc: 0, $	 ; galactic component (1=thin, 2=thick, 3=halo, 4=bulge)
                  starID: 0., $ ; unique star ID
                  coord: {coordstruct} $ ; give the star unique coordinates
                  }
end
