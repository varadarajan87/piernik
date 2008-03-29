 $DOMAIN_SIZES
    nxd = 25
    nyd = 25
    nzd = 25
    nb  = 5
 /  

 $MPI_BLOCKS
    pxsize = 1
    pysize = 1
    pzsize = 1
 /

 $START_CONTROL
 /

 $RESTART_CONTROL
    restart  = 'last'
    new_id   = ''
    nrestart = 0
 /

 $END_CONTROL
    tend   = 1.0e0
    nend   = 100000 
 /

 $OUTPUT_CONTROL
    dt_hdf  = 1.0e0
    dt_res  = 1.0e0
    dt_log  = 0.001
    dt_tsl  = 0.00001
    vars(1) = 'dens'
    vars(2) = 'velx'
    vars(3) = 'vely'
    vars(4) = 'velz'
    vars(5) = 'ener'
    vars(6) = 'magx'
    vars(7) = 'magy'
    vars(8) = 'magz'
    vars(9) = 'encr'
    domain  = 'full_domain'
 / 

 $BOUNDARIES
    bnd_xl = 'cor'  
    bnd_xr = 'out'  
    bnd_yl = 'cor'  
    bnd_yr = 'out'  
    bnd_zl = 'outd'  
    bnd_zr = 'outd'
 /

 $DOMAIN_LIMITS
    xmin   = 0.0
    xmax   = 1.6e4
    ymin   = 0.0
    ymax   = 1.6e4
    zmin   =-0.5e3
    zmax   = 0.5e3
 /

 $EQUATION_OF_STATE
    c_si   = 10.0
    alpha  = 1.0
 /  

 $NUMERICAL_SETUP
    cfl    = 0.9
    smalld = 1.e-6
    smallei= 1.e-4
    integration_order = 2
    magnetic = 'yes'
    dimensions = '3d'
 /
 
 $GRAVITY
    nsub   = 20
    h_grav = 300.0
    n_gravh = 2
 /
 
 $THERMAL
 /
 
 $RESISTIVITY
    cfl_resist   =    0.9   ! liczba Couranta dla opornosci, (wsp. 0.5) w kodzie 
    eta_0   	 =  100.0   ! opornosc jednorodna
    eta_1        =    0.0   ! opornosc zlokalizowana, dziala powyzej j_crit
    j_crit       = 1000.0   ! prad krytyczny
 /
 
  $COSMIC_RAYS
    cr_active    = 1.0      ! gdy 0 wtedy grad p_cr nie jest brany pod uwage w r.ruchu
    cfl_cr       = 0.8      ! liczba couranta dla dyfuzji CR (wsp 0.5 w kodzie)
    cr_eff       = 0.1      ! czesc energii SN (=10^51 erg) zamieniona w CR
    gamma_cr     = 1.555555 ! wykl. adiabatyczny CR
    K_cr_paral   = 1000.   ! wsp. dyfuzji CR wzdluz B
    K_cr_perp    = 100.    ! wsp. dyfuzji CR wpoprzek B
    beta_cr      = 1.0      ! okresla wklad CR w poczatkowej rownowadze (podobnie jak alfa dla B)
    amp_cr       = 0.0      ! wybuch w init_problem w x0,y0,z0
    smallecr     = 1.0e-6    ! lower limit & outflow boundary value for ecr
 /
 
 $GALACTIC_PARAMS
 /

 $SN_DISTR
 snenerg = 1.0
 sn1time = 445.e2
 sn2time = 52.e2
 /

 $PROBLEM_CONTROL
    problem_name ='galdisk'
    run_id  =  'ts0'
    rhoa    = 1.0e-4
    d0      = 1.0
    r_max   = 1.2e4
    mtr      = 10
    mf_orient = 'toroidal'
 /

