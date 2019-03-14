function orbital_motion,PLX=plx,M1=m1,M2=m2,RHO=rho,E_RHO=e_rho,N_DAYS=n_days,FACTOR=factor

  ;;Usage example:
  ;;  result=ORBITAL_MOTION(PLX=97.84,M1=0.3387,M2=0.2280,RHO=0.16,E_RHO=0.1,N_DAYS=1001)

  ;;Making sure the variables are doubles
  plx=DOUBLE(plx)
  m1=DOUBLE(m1)
  m2=DOUBLE(m2)
  rho=DOUBLE(rho)
  e_rho=DOUBLE(e_rho)
  n_days=DOUBLE(CEIL(n_days))

  ;;Plot every nth day, make sure n_days is divisible by n_plot
  n_plot = 15.0
  WHILE (n_days MOD n_plot) ne 0 DO n_days+=1.0

  ;;Define how finely sampled omega and a are sampled
  n_w = 51
  n_a = 251
  
  ;;Define some constants
  au_to_m = 149597870700.0d
  sec_to_yr = 3.16888d-8
  msol_to_kg = 1.9891d30
  g=6.67384d-11
  
  ;;Convert days to years
  n_years = n_days/365.25

  ;;Create array to store min and max seps
  minmax = DBLARR(n_days/n_plot,2)
  minmax[*,1]=1d10 ;;Making sure that the max by setting the first value really high.

  ;;Looping over two rho values
  old_rho=rho

  ;plot,[0],XRANGE=[0,n_days],YRANGE=[0,2]

  FOR k=0,1 DO BEGIN
     IF k eq 0 THEN rho=old_rho-e_rho
     IF k eq 1 THEN rho=old_rho+e_rho

     ;;Range of a values to search, going to do it equally spaced in log(a)
     ;;Going from rho+0.01 to (rho+0.01)*1000
     a_min = ALOG10(rho+0.01)
     a_max = ALOG10((rho+0.01)*factor)

     ;;Now create omega and a arrays
     a = (10^(((DINDGEN(n_a)/(n_a-1.))*(a_max-a_min))+a_min))
     w = (DINDGEN(n_w)/(n_w-1.))*(2.0*!dpi)
     

     FOR i=0,n_a-1 DO BEGIN
        ;;Work out au separation and period from masses
        au = a[i] * (1.0/plx)
        P = SQRT((((au*au_to_m)^3.0)*4.0*!dpi*!dpi)/(g*(m1+m2)*msol_to_kg)) * sec_to_yr
        
        FOR j=0,n_w-1 DO BEGIN
           ;;Calculate true anomaly (which equals mean anomaly for
           ;;circular orbits) where the separation is equal to the
           ;;measured value. There should be four of these
           nu1=((-1.0)*ACOS((-1.0)*SQRT((rho^2.0)/(a[i]^2.0))))-w[j]
           nu2=((1.0)*ACOS((-1.0)*SQRT((rho^2.0)/(a[i]^2.0))))-w[j]
           nu3=((-1.0)*ACOS((1.0)*SQRT((rho^2.0)/(a[i]^2.0))))-w[j]
           nu4=((1.0)*ACOS((1.0)*SQRT((rho^2.0)/(a[i]^2.0))))-w[j]
           
           ;;Making sure nu values aren't negative here
           IF nu1 lt 0 THEN nu1+=(2.0*!dpi)
           IF nu2 lt 0 THEN nu2+=(2.0*!dpi)
           IF nu3 lt 0 THEN nu3+=(2.0*!dpi)
           IF nu4 lt 0 THEN nu4+=(2.0*!dpi)
           
           IF nu1 lt 0 THEN nu1+=(2.0*!dpi)
           IF nu2 lt 0 THEN nu2+=(2.0*!dpi)
           IF nu3 lt 0 THEN nu3+=(2.0*!dpi)
           IF nu4 lt 0 THEN nu4+=(2.0*!dpi)
           
           IF nu1 lt 0 THEN nu1+=(2.0*!dpi)
           IF nu2 lt 0 THEN nu2+=(2.0*!dpi)
           IF nu3 lt 0 THEN nu3+=(2.0*!dpi)
           IF nu4 lt 0 THEN nu4+=(2.0*!dpi)
           
           nu_arr=[nu1,nu2,nu3,nu4]
           nu_arr=nu_arr[SORT(nu_arr)]
           
           ;;How many periods can we fit in to n_years?
           n_p = (n_years/P)
           
           ;;Calculate the tru anomaly for each day requested. nu1=nu3,
           ;;nu2=nu4 so they are not plotted
           this_nu0 = ((((DINDGEN(n_days/n_plot)/((n_days/n_plot)-1.))*(2.0D*!dpi*n_p)) + nu_arr[0]) MOD (2.0D*!dpi))
           this_nu1 = ((((DINDGEN(n_days/n_plot)/((n_days/n_plot)-1.))*(2.0D*!dpi*n_p)) + nu_arr[1]) MOD (2.0D*!dpi))
           
           ;;Then work out rho on each day
           rho0 = (a[i] * ABS(COS(this_nu0+w[j])))
           rho1 = (a[i] * ABS(COS(this_nu1+w[j])))
           
           ;;Replace maximal values if they are bigger in the current orbit, and minimal
           ;;ones if smaller ones in the current orbit
           minmax[*,0]=minmax[*,0]>rho0
           minmax[*,0]=minmax[*,0]>rho1
           minmax[*,1]=minmax[*,1]<rho0
           minmax[*,1]=minmax[*,1]<rho1
        ENDFOR
     ENDFOR
  ENDFOR

  result=[[DINDGEN(n_days/n_plot)*n_plot],[minmax[*,0]],[minmax[*,1]]]
  RETURN,result
END
