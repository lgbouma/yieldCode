function make_eclipse, sstruct, bkstruct, estruct, frac, ph_p, dartstruct, tefftic, $
  eclass, tband, min_depth=min_depth, max_depth=max_depth, ps_only=ps_only, pla_err=pla_err
  ecliplen = 0L
  if (keyword_set(min_depth)) then min_depth=min_depth else min_depth = 0.0
  if (keyword_set(max_depth)) then max_depth=max_depth else max_depth = 1.0 

  print, 'Adding eclipses to ', n_elements(sstruct), ' stars.'
  
  ; Add HEBs
  if (eclass[3]) then begin
    gd = add_hebs(sstruct, heb_eclip, frac, ph_p, dartstruct, tefftic, ps_only=ps_only)
    if (gd gt 0) then begin
      if (ecliplen gt 0) then estruct = struct_append(estruct, heb_eclip) $
      else estruct = heb_eclip   
    endif
    ecliplen = ecliplen + gd
  endif
  
  ; Add Planets to all target stars
  if (eclass[0]) then begin
    gd = add_planets(sstruct, p_eclip, frac, ph_p, tband, err=pla_err, min_depth=min_depth, dressing=1, ps_only=ps_only)
    if (gd gt 0) then begin
      if (ecliplen gt 0) then estruct = struct_append(estruct, p_eclip) $
      else estruct = p_eclip
    endif
    ecliplen = ecliplen + gd
  endif

  ; Add EBs (identify EBs among target stars)
  if (eclass[1]) then begin
    gd = add_ebs(sstruct, eb_eclip, frac, ph_p, max_depth=max_depth)
    if (gd gt 0) then begin
      if (ecliplen gt 0) then estruct = struct_append(estruct, eb_eclip) $
      else estruct = eb_eclip   
    endif
    ecliplen = ecliplen + gd
  endif
  
  ; Add BEBs (identify EBs among background stars, attach to target stars)
  if (eclass[2]) then begin
    gd = add_bebs(sstruct, bkstruct, beb_eclip, frac, ph_p, 100, ps_only=ps_only)
    if (gd gt 0) then begin
      if (ecliplen gt 0) then estruct = struct_append(estruct, beb_eclip) $
      else estruct = beb_eclip   
    endif
    ecliplen = ecliplen + gd
  endif
  
  ; Add BTPs (identify EBs among background stars, attach to target stars)
  if (eclass[4]) then begin
    gd = add_btps(sstruct, bkstruct, btp_eclip, frac, ph_p, tband, 100, dressing=1, min_depth=min_depth, ps_only=ps_only)
    if (gd gt 0) then begin
      if (ecliplen gt 0) then estruct = struct_append(estruct, btp_eclip) $
      else estruct = btp_eclip   
    endif
    ecliplen = ecliplen + gd
  endif
  
;  if (keyword_set(ps_only)) then begin
;    egd = where(sstruct[estruct.hostid].ffi lt 1)
;    if (egd[0] ne -1) then begin
;      estruct = estruct[egd]
;      ecliplen = n_elements(egd)
;    endif else ecliplen = 0  
;  endif
  return, ecliplen
END
