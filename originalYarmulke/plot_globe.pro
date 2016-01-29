pro plot_globe

  restore, filen='survey.sav'

  vdec = 45
  vlong = 0
  zrot = 30

  MAP_SET, vdec, vlong, zrot, /ORTHO, /ISOTROPIC, /HORIZON, xmargin=0,ymargin=0,/noborder
  img_Number_of_Pointings  = MAP_IMAGE(Number_of_Pointings,Startx,Starty, COMPRESS=1)
  atv,img_Number_of_Pointings

end
