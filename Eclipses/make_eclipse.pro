function make_eclipse, sstruct, bkstruct, estruct, frac, ph_p, dartstruct, tefftic, $
  eclass, tband, noEclComp, min_depth=min_depth, max_depth=max_depth, ps_only=ps_only, pla_err=pla_err, $
  extMission=extMission, burtCatalog=burtCatalog
;+
;NAME: make_eclipse
;PURPOSE: wrap around add_hebs, add_planets, add_ebs, add_bebs, which are the
; routines that assign planets and generate "eclip" objects (data struct with
; all eclipse information). Give these eclip objects back to tile_wrapper.
;INPUTS:
; 1. sstruct. An outStarStruct object (targets), selected for postage stamps and FFIs.
; 2. bkstruct. Passed "bkgnds" a starStruct object for K>15 stars (dimmer->diluting)
; 3. estruct. Passed "eclip_trial", undefined at the call in tile_wrapper (i.e., at
;		time of call, has no data).
; 4. frac. PRF data from Deb Woods.
; 5. ph_p. Photon fluxes for T=10 vs T_eff.
; 6. dartstruct. dartStruct object from "dartmouth_grid.sav". 
; 7. tefftic. "tic_fits", from "tic_teff.fits". (?)
; 8. eclass. Include planets / EBs / BEBs / HEBs / BTPs.
; 9. tband. TESS transmission bandpass vs wavelength (cf Sullivan+ Fig 2).
; 10. min_depth. Minimum transit depth to retain from eclipses (nominally 1e-6)
; 11. max_depth. Max '' (nominally 1).
; 12. ps_only. Bool value for whether or not only postage stamps.
; 13. pla_err. 0: (standard) nominal occurrence rates / +1: upper bounds / -1: lower bounds.
; 14. noEclComp. `eclip` data structure containing the data from planets who have planets
;		in their system that transit, although they do not (Burt Catalog relevant).
;OUTPUTS:
; 1. estruct, an array of eclip_struct objects. E.g., `make_eclipse` will make 3
;	  transiting planets out of a total of 78 created planets about 700 stars on
;	  a single tile. Each transiting planet gets its own eclip_struct entry in the 
;     array.
;
;-
  ecliplen = 0L
  if (keyword_set(min_depth)) then min_depth=min_depth else min_depth = 0.0
  if (keyword_set(max_depth)) then max_depth=max_depth else max_depth = 1.0 

  print, 'Make_eclipse: adding eclipses to ', n_elements(sstruct), ' stars (targets).'
  
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
    gd = add_planets(sstruct, p_eclip, frac, ph_p, tband, noEclComp, err=pla_err, $ 
					 min_depth=min_depth, dressing=1, ps_only=ps_only, extMission=extMission,$
					 burtCatalog=burtCatalog)
    if (gd gt 0) then begin
      if (ecliplen gt 0) then estruct = struct_append(estruct, p_eclip) $ ; append eclipsing planets
      else estruct = p_eclip
    endif
    ecliplen += gd
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
