function get_pm,RA=ra,DE=de,PLX=plx,PM_RA=pm_ra,PM_DE=pm_de,RHO=rho,THETA=theta,MJD0=mjd0,N_MJD=n_mjd,N_BEFORE=n_before

  ;;Function to get change in RA and Dec as a function of MJD for a
  ;;given coordinate, proper motion, and parallax. This assumes you
  ;;are measuring delta RA, delta DE, rho and pa relative to the central
  ;;star. This is important to consider, as the proper motion vectors
  ;;will have to be multiplied by -1 because background objects will
  ;;have a proper motion of (-pmra, -pmde), and the central star will
  ;;be stationary if you have aligned the images on the central star.

  ;; Input:
  ;; ra = target RA (decimal degrees)
  ;; de = target DEC (decimal degrees)
  ;; plx = target parallax (arcseconds)
  ;; pm_ra = proper motion in RA (arcseconds /yr)
  ;; pm_de = proper motion in DEC (arcseconds /yr)
  ;; rho - separation of companion in FIRST epoch (in arcseconds)
  ;; theta - position angle of companion in FIRST epoch (degrees)
  ;; mjd0 = starting MJD
  ;; n_mjd = number of days for which proper motion curve is desired
  ;; (e.g. for 2 years would be 365.25 x 2 = 730.5 days)
  ;; n_before = number of days before to plot too.

  ;;Output:
  ;; result - an (n_mjd X 5) array containing:
  ;;      MJD, delta_ra, delta_de, delta_rho, delta_theta
  ;; for each MJD required.
  ;; delta_ra - change in the RA offset between star and background object
  ;; delta_de - change in the DEC offset between star and background object
  ;; delta_rho - change in the rho of a companion measured relative to the star if it is an unnasociated background object
  ;; delta_theta - ditto for the position angle

  ;;Usage example:
  ;; result = GET_PM(RA=308.82723181, DE=14.67421319, PLX=0.01482, PM_RA=0.04552, PM_DE=0.01174, RHO=13.51, THETA=143.35, MJD0=54717.2398, N_MJD=1000.0, N_BEFORE=100.0)

  ;;Should check if variables are declared

  ;;Define path to JPL ephemiris
  path='/Users/Kim/IDL/pro/kdw/'

  ;;Convert to radians
  ra_rad = ra * (!dpi/180.0)
  de_rad = de * (!dpi/180.0)

  IF mjd0 gt 1000000.0 THEN BEGIN
     print,'Use MJD not JD! (MJD = JD-2400000), unless of course this is being used in the year 8000AD...'
     stop
  ENDIF

  mjd_arr = FINDGEN(n_mjd+n_before+1) + (mjd0-n_before)

  ;;Load the position of the Earth at each MJD requested
  
  JPLEPHREAD,path+'JPLEPH.405',pinfo,pdata
  JPLEPHINTERP,pinfo,pdata,mjd_arr+2400000.0,X,Y,Z,/EARTH,posunits='AU'

  ;;Also load the parallax of the object at mjd0
  
  JPLEPHREAD,path+'JPLEPH.405',pinfo,pdata
  JPLEPHINTERP,pinfo,pdata,mjd0+2400000.0,X0,Y0,Z0,/EARTH,posunits='AU'

  ;;Now determine how much of a shift in RA and DEC the parallax is
  ;;causing at each MJD requested

  plx_ra = plx * (1/COS(de_rad)) * ((X*SIN(ra_rad)) - (Y*COS(ra_rad)))
  plx_de = plx * ((X*COS(ra_rad)*SIN(de_rad)) + (Y*SIN(ra_rad)*SIN(de_rad))-(Z*COS(de_rad)))

  ;;And do the same for MJD0

  plx_ra0 = plx * (1/COS(de_rad)) * ((X0*SIN(ra_rad)) - (Y0*COS(ra_rad)))
  plx_de0 = plx * ((X0*COS(ra_rad)*SIN(de_rad)) + (Y0*SIN(ra_rad)*SIN(de_rad))-(Z0*COS(de_rad)))
  
  ;; Determine the cumulative proper motion in RA and DEC

  pm_ra_c = (((mjd_arr-mjd0)/365.25)*pm_ra)
  pm_de_c = (((mjd_arr-mjd0)/365.25)*pm_de)
  ;;Now add the proper motion to the parallax to get the change in RA
  ;;and DEC as a function of time

  delta_ra = (pm_ra_c + (plx_ra - plx_ra0))
  delta_de = (pm_de_c + (plx_de - plx_de0))

  ;;Multiply by -1 since the measurmeents are done relative to the moving star, so background objects are not stationary, they move in the opposite direction to the measured proper motion of the primary.

  delta_ra *= (-1.0)
  delta_de *= (-1.0)
  
  ;;This can be converted into a change in rho and theta. First determine
  ;;the RA and DEC offset of the companion in the first epoch

  offset_ra0 = rho * SIN(theta*(!dpi/180.0))
  offset_de0 = rho * COS(theta*(!dpi/180.0))

  ;;Now work how the RA and DEC offsets would change as a function of time if the object were a background object. Do this by adding the offset within the first epcoh (offset_ra0, offset_de0) and the change in RA and DEC as a function of time:

  offset_ra = offset_ra0 + delta_ra
  offset_de = offset_de0 + delta_de

  ;;Then convert this into a rho and theta for a background object

  delta_rho = SQRT( (offset_ra^2.0) + (offset_de^2.0) )
  delta_theta = ((ATAN(offset_de,-offset_ra)*(180.0/!dpi))+270.0) MOD 360.0

  ;;Errors - a proper treatment would be to run this routine many times drawing values of plx, pmra, pmde, rho, theta from a gaussian distribution centred on the measured value for these parameters, of 1sigma width equal to the stated uncertainties on the measured values. Then at each MJD, we find what the 1sigma range of delta_rho and delta_pa values are, and these are the upper and lower error bounds for the delta_rho and delta_pa curves.

  ;;I will work on this (might want to include it on my paper), but thought I'd send the code without errors for now.

  result = [[mjd_arr],[delta_ra],[delta_de],[delta_rho],[delta_theta]]

  RETURN,result
END
  
