function delfosse_mass_given_mj, mj

  logm = 1d-3*(1.6 + 6.01*mj + 14.888*mj^2. - 5.3557*mj^3. + 0.285181*mj^4.)
  m=10.^logm
  return, m

end
