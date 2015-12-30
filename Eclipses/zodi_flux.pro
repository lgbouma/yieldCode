PRO zodi_flux, elat, aspix, zodi_ph
  ; e/pix from zodi
  dlat = (abs(elat)-90.)/90.
  vmag_zodi = 23.345 - 1.148*dlat^2.
  zodi_ph = 10.0^(-0.4*(vmag_zodi-22.8)) * 2.56D-3 * aspix^2.
END
