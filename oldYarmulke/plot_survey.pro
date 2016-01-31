pro plot_survey, filen=filen

  if (keyword_set(filen)) then filen=filen else filen='survey.sav'
  restore, filen

  !P.multi=0

  image = Number_of_Pointings
  minValue = Floor(Min(image))
  maxValue = Ceil(Max(image))
  nLevels = 10
  xtitle = 'ecliptic longitude'
  ytitle = 'ecliptic latitude'
  position =   [0.125, 0.125, 0.9, 0.800]
  cbposition = [0.125, 0.865, 0.9, 0.895]
  cbTitle = 'Data Value'
  
  cgDisplay, 1200, 600, Title='Image Plot with Contours'
  cgLoadCT, 33, CLIP=[30,255]
  cgImage, image, Stretch=1, MinValue=minValue, MaxValue=maxValue, $
           /Axes, XTitle=xtitle, YTitle=ytitle, Position=position, $
           XRange=[-180,180], YRange=[-90, 90], /Keep_Aspect

end
