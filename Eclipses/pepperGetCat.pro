function pepperGetCat, targets, selpri, selsing
;+
; PURPOSE: produce a catalog (of TRILEGAL star properties)to be saved to 
; 	CSV good for filtergraph comparison.
;-

	; Some day, I may get around to giving stars coordinates.
    v = [targets[selpri].mag.v, targets[selsing].mag.v]
    ic = [targets[selpri].mag.ic, targets[selsing].mag.ic]
    j = [targets[selpri].mag.j, targets[selsing].mag.j]
    kp = [targets[selpri].mag.kp, targets[selsing].mag.kp]
    t = [targets[selpri].mag.t, targets[selsing].mag.t]
    teff = [targets[selpri].teff, targets[selsing].teff]
    logg = [targets[selpri].logg, targets[selsing].logg]
    rad = [targets[selpri].r, targets[selsing].r]
    m = [targets[selpri].m, targets[selsing].m]

    pepperDatThisTile = [[v], [ic], [j], [kp], [t], [teff], [logg], [rad], [m]]

	RETURN, pepperDatThisTile
END
