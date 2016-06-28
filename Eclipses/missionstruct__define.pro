PRO missionstruct__define

 mission = {missionstruct, $
            npointings: 0, $
            dep1_eff: 0., $ ; effective depth (note std depth is fixed btwn pri/ext)
            dep2_eff: 0., $ ; effective secondary eclipse depth
            dur1_eff: 0., $ ; effective duration (dur1 fixed btwn pri/ext; _eff is not)
            dur2_eff: 0., $ ; effective duration
            neclip_obs1: 0, $ ; Number of primary eclipses observed
            neclip_obs2: 0, $ ; Number of secondary eclipses observed
            snr: 0., $ ; SNR of primary+secondary eclipses in phase-folded lightcurve over mission
            snrhr: 0., $ ; SNR of target per hour
            snr1: 0., $ ; SNR of primary eclipses in phase-folded lightcurve
            snr2: 0., $ ; SNR of secondary eclipses in phase-folded lightcurve
            snreclp1: 0., $ ; SNR per primary eclipse
            snreclp2: 0., $ ; SNR per secondary eclipse
            det: 0, $  ; Detected (in this mission)?
            det1: 0, $  ; Detected primary?
            det2: 0 $  ; Detected secondary?
           }
END
