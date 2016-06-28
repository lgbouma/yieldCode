PRO eclipstruct__define

  eclip = {eclipstruct,  $
        coord: {coordstruct}, $
        pri: {missionstruct}, $
        ext: {missionstruct}, $
        totpointings: 0, $ ; total pointings throughout primary & ext (.pri.npointings + .ext.npointings)
        class: 0, $ ;1 planet around field star, 2=planet around binary, 3=heb, 4=beb/eb
        trial: 0, $
        r1: 0.0, $ ; radius of primary
        r2: 0.0, $ ; radius of secondary
        m1: 0.0, $ ; mass of primary
        m2: 0.0, $ ; mass of secondary
        k: 0.0, $
        teff1: 0.0, $ ; temp of primary
        teff2: 0.0, $ ; temp of secondary
        p: 0.0, $ ; Period (days)
        a: 0.0, $ ; Semimajor axis (AU)
        s: 0.0, $ ; Insolation of secondary
        w: 0.0, $ ; Argument of perigee
        f: 0.0, $ ; True anomaly
        mult: 0.0, $ ; Multi-planet system?
        tmult: 0.0, $ ; Multi-planet system?
        pr: 0.0, $ ; Period ratio for multi
        ecc: 0.0, $ ; Eccentricity
        cosi: 0.0, $ ; -1 to 1
        b: 0.0, $ ; Impact parameter (0-1)
        dep1: 0.0, $ ; Eclipse depth (0-1)
        dep2: 0.0, $ ; Eclipse depth (0-1)
        dur1: 0.0, $ ; Primary eclipse duration (days)
        dur2: 0.0, $ ; Primary eclipse duration (days)
        censhift1: 0.0, $ ; centroid shift (BEBs and HEBs only)
        censhift2: 0.0, $ ; centroid shift of sec eclipse (pixels)
        cenerr1: 0.0, $ ; centroid error in shift direction (pixels)
        cenerr2: 0.0, $ ; error in sec eclipse direction
        gress1: 0.0, $    ; duration of in+engress (days)
        gress2: 0.0, $    ; duration of in+engress (days)
        snrgress1: 0.0, $   ; SNR of ingress/egress
        snrgress2: 0.0, $   ; SNR of ingress/egress
        snrf: 0., $ ; SNR of phase folded light curve, accounting for transits&occultations from pri+ext
        detf: 0, $  ; Detected based on snrf? (and a threshold)
        npix: 0, $       ; Optimal number of pix in aperture
        tsys: 0.0, $      ; system tmag (for binaries)
        kpsys: 0.0, $      ; system tmag (for binaries)
        icsys: 0.0, $      ; system tmag (for binaries)
        jsys: 0.0, $      ; system tmag (for binaries)
        star_ph: 0.0, $      ; photons/s/cm^2 from star
        bin_ph: 0.0, $      ; photons/s/cm^2 from within 0.5 pix
        bk_ph: 0.0, $      ; photons/s/cm^2/pix from other stars
        zodi_ph: 0.0, $      ; photons/s/cm^2/pix
        sep: 0.0, $         ; separation (in pixels) of eclipse from host
        sat: 0, $        ; saturation flag
        dil: 0.0, $        ; dilution ratio
        var: 0.0, $        ; host star variability
        hostid: 0L, $
        tileNum: 0, $ ; tile number of primary
        ffiClass: 0, $; 0: none, 1: ffi, 2: PS only in primary, 3: PS only in extended, 4: PS in both
        uniqEclipID: 0UL $; integer unique to this eclip. New trials get diff, but primary/ext keep same.
	    }
end
