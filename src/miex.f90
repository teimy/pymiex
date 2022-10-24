! ****************************************************************************************************
! MIEX: MIE SCATTERING CODE FOR LARGE GRAINS
!                                                _____________________________________________________
!                                                Contact information:   swolf@mpia.de (Sebastian Wolf)
! ====================================================================================================
! GENERAL CODE DESCRIPTION
! ------------------------
!
! Based on Mie scattering theorie, the following quantities for
!
!       (a) single grain sizes / chemical components
!   and (b) mixtures of chemically different grains with a size distributions
!
! can be calculated:
!
!   * Scattering matrix elements S11, S12, S33, and S34
!   * Extinction effiency factor         (Qext) & Extinction cross section     (Cext)
!   * Scattering effiency factor         (Qsca) & Scattering cross section     (Csca)
!   * Absorption effiency factor         (Qabs) & Absorption cross section     (Cabs)
!   * Backscattering effiency factor     (Qbk)  & Backscattering cross section (Cbk)
!   * Radiation pressure effiency factor (Qpr)
!   * Albedo
!   * Scattering assymetry factor (g).
!
! ____________________________________________________________________________________________________
! The optical data of the grains have to be provided in files with the following tabular form
!   * first  row: wavelength [micron]
!   * second row: n (=real[ri])
!   * second row: k (=imag[ri]).
! 
! Rem.: For astrophysical applications, such tables can be found, e.g., at
!       http://www.astro.uni-jena.de/Users/database/entry.html
!       ("Jena-Petersburg Database of Optical Constants"). See, also, 
!       N.V.Voshchinnikov: "Optics of Cosmic Dust",
!                           Astrophysics and Space Physics Review 12,  1 (2002)
!       for further references.
!
! ====================================================================================================
! SUBROUTINES / MODULES
! ---------------------
!
! (1) Module 'datatype'     : data type definitions
! (2) Module 'mie_routines'
!     (2.1) aa2      : calculations of the ratio of derivative to the function for Bessel functions
!                      of half order with complex argument: J'(n)/J(n)
!     (2.2) shexqnn2 : derive quantities (listed above) for a single size parameter
!                      and chemical composition`
! ****************************************************************************************************

! ====================================================================================================
! Definition of the data type to be used for floating point operations: 
! - r1: real*4  (FORTRAN 77: "real")
! - r2: real*8  (FORTRAN 77: "double precision")
! ====================================================================================================
module datatype
  implicit none
  integer, parameter, public ::  r1=selected_real_kind(1)  ! real*4
  integer, parameter, public ::  r2=selected_real_kind(9)  ! real*8 (double precision)
end module datatype


! ====================================================================================================
! Collection of subroutines
! ====================================================================================================
module mie_routines
  private :: aa2
  public  :: shexqnn2
contains
  ! ==================================================================================================
  ! Subroutine for calculations of the ratio of derivative to the function for Bessel functions
  ! of half order with complex argument: J'(n)/J(n). The calculations are given by the recursive
  ! expression ``from top to bottom'' beginning from n=num.
  ! *  a=1/x (a=2*pi*a(particle radius)/lambda - size parameter).
  ! *  ri - complex refractive index.
  ! *  ru-array of results.
  ! - this routine is based on the routine 'aa' published by
  !       N.V.Voshchinnikov: "Optics of Cosmic Dust",
  !                           Astrophysics and Space Physics Review 12,  1 (2002)
  ! ==================================================================================================
  subroutine aa2( a, ri, num, ru )
    use datatype

    implicit none

    ! variables for data exchange.....................................................................
    real(kind=r2), intent(in)                   :: a
    complex(kind=r2), intent(in)                :: ri
    integer, intent(in)                         :: num
    complex(kind=r2), dimension(:), intent(out) :: ru

    ! local variables.................................................................................
    integer :: i, i1, j, num1
    complex(kind=r2) :: s, s1
    !-------------------------------------------------------------------------------------------------
    ! initialisierung: not necessary (+ slowes the code down remarkably)
    ! ru(:) = (0.0, 0.0)

    s       = a / ri
    ru(num) = real(num+1,kind=r2) * s
    num1    = num - 1
    do j=1, num1
       i     = num - j
       i1    = i + 1
       s1    = i1 * s
       ru(i) = s1 - 1.0_r2 / (ru(i1) + s1)
    end do
  end subroutine aa2


  !===================================================================================================
  ! shexqnn2
  ! --------
  ! - for a given size parameter 'x' and (complex) refractive index 'ri' the following quantities
  !   are determined:
  !   * Qext     - extinction effiency
  !   * Qsca     - scattering effiency
  !   * Qabs     - absorption effiency
  !   * Qbk      - backscattering effiency
  !   * Qpr      - radiation pressure effiency
  !   * albedo   - Albedo
  !   * g        - g scattering assymetry factor
  !   * SA1, SA2 - scattering amplitude function
  ! - further input parameters
  !   * doSA = .true.  ->  calculation of the scattering amplitudes
  !   * nang ... half number of scattering angles theta in the intervall 0...PI/2
  !              (equidistantly distributed)
  ! - this routine is based on the routine 'shexqnn' published by
  !       N.V.Voshchinnikov: "Optics of Cosmic Dust",
  !                           Astrophysics and Space Physics Review 12,  1 (2002)
  !===================================================================================================
  subroutine shexqnn2( ri, x, Qext, Qsca, Qabs, Qbk, Qpr, albedo, g, ier, SA1, SA2, doSA, nang )
    use datatype

    implicit none
    
    ! variables for data exchange.....................................................................
    complex(kind=r2), intent(in)                :: ri
    real(kind=r2), intent(in)                   :: x
    real(kind=r2), intent(out)                  :: Qext, Qsca, Qabs, Qbk, Qpr, albedo, g
    integer, intent(out)                        :: ier
    complex(kind=r2), dimension(:), intent(out) :: SA1, SA2
    logical, intent(in)                         :: doSA
    integer, intent(in)                         :: nang

    ! local variables.................................................................................
    integer       :: iterm, nterms, num, iu0, iu1, iu2, iang2, iang
    real(kind=r2) :: r_iterm, factor, eps, pi, ax, besJ0, besJ1, besJ2, besY0, besY1, besY2, b, an, &
         y, ass, w1, qq, fac, an2, P, T, Si, Co, z, xmin
    complex(kind=r2) :: ra0, rb0, ra1, rb1, r, ss, s1, s2, s3, s, rr

    real(kind=r2), dimension(0:1)               :: fact
    real(kind=r2), dimension(:), allocatable    :: mu, fpi, fpi0, fpi1, ftau
    complex(kind=r2), allocatable, dimension(:) :: ru
    !-------------------------------------------------------------------------------------------------
    ! Maximum number of terms to be considered
    nterms = 20000000

    ! Accuracy to be achieved
    eps    = 1.0e-20_r2

    ! Minimum size parameter
    xmin   = 1.0e-6_r2

    !-------------------------------------------------------------------------------------------------
    ! initialization
    allocate( ru(1:nterms), mu(1:nang), fpi(1:nang), fpi0(1:nang), fpi1(1:nang), ftau(1:nang) )
    ier     = 0
    Qext    = 0.0_r2
    Qsca    = 0.0_r2
    Qabs    = 0.0_r2
    Qbk     = 0.0_r2
    Qpr     = 0.0_r2
    albedo  = 0.0_r2
    g       = 0.0_r2
    fact(0) = 1.0_r2
    fact(1) = 1.0e+250_r2
    factor  = 1.0e+250_r2

    ! null argument
    if (x <= xmin) then
       ier = 1
       print *, "<!> Error in subroutine shexqnn2:"
       print *, "    - Mie scattering limit exceeded:"
       print *, "      current size parameter: ", x
    else
       pi = 4.0_r2 * atan(1.0_r2) ! PI = 3.14...
       ax = 1.0_r2 / x
       b  = 2.0_r2 * ax**2
       ss = (0.0_r2, 0.0_r2)
       s3 = (0.0_r2,-1.0_r2)
       an = 3.0_r2
       
       ! define the number for subroutine aa2 [Loskutov (1971)]
       y   = sqrt( RI * conjg(ri) )  *  x
       num = 1.25 * y + 15.5
       
       if      ( y<1.0_r2 ) then
          num = 7.5 * y + 9.0
       else if ( (y>100.0_r2) .and. (y<50000.0_r2) ) then
          num = 1.0625 * y + 28.5
       else if ( y>=50000.0_r2 ) then
          num=1.005*y+50.5
       end if

       if(num > nterms) then
          ier = 2
          print *, "<!> Error in subroutine shexqnn2:"
          print *, "    - Maximum number of terms  : ", nterms
          print *, "    - Number of terms required : ", num
          print *, "    ** Solution: Increase default value of the variable 'nterm' **"
       else
          ! logarithmic derivative to Bessel function (complex argument)
          call aa2(ax,ri,num,ru)

          ! ------------------------------------------------------------------------------------------
          ! FIRST TERM
          ! ------------------------------------------------------------------------------------------
          ! initialize term counter
          iterm = 1
          
          ! Bessel functions
          ass = sqrt( pi / 2.0_r2 * ax )
          w1  = 2.0_r2/pi * ax
          Si  = sin(x)/x
          Co  = cos(x)/x
          
          ! n=0
          besJ0 =  Si / ass
          besY0 = -Co / ass
          iu0   = 0
          
          ! n=1
          besJ1 = ( Si * ax - Co) / ass
          besY1 = (-Co * ax - Si) / ass
          iu1   = 0
          iu2   = 0
          
          ! Mie coefficients
          s   = ru(1) / ri + ax
          s1  = s * besJ1 - besJ0
          s2  = s * besY1 - besY0
          ra0 = s1 / (s1 - s3 * s2)   ! coefficient a_1
          
          s   = ru(1) * ri + ax
          s1  = s * besJ1 - besJ0
          s2  = s * besY1 - besY0
          rb0 = s1 / (s1 - s3 * s2)   ! coefficient b_1
          
          ! efficiency factors
          r    = -1.5_r2 * (ra0-rb0)
          Qext = an * (ra0 + rb0)
          Qsca = an * (ra0 * conjg(ra0)  +  rb0 * conjg(rb0))
          
          ! scattering amplitude functions
          if (doSA) then
             do iang=1, nang
                mu(iang) = cos( (real(iang,kind=r2)-1.0_r2) * (pi/2.0_r2)/real(nang-1,kind=r2) )
             end do

             fpi0(:) = 0.0_r2
             fpi1(:) = 1.0_r2
             SA1(:)  = cmplx( 0.0_r2, 0.0_r2 )
             SA2(:)  = cmplx( 0.0_r2, 0.0_r2 )
             
             r_iterm = real(iterm,kind=r2)  ! double precision             
             fac     = (2.0*r_iterm + 1.0_r2) / (r_iterm * (r_iterm+1.0_r2))
             
             do iang=1, nang
                iang2      = 2 * nang - iang
                
                fpi(iang)  = fpi1(iang)
                ftau(iang) = r_iterm * mu(iang) * fpi(iang)  -  (r_iterm+1.0) * fpi0(iang)
                
                P          = (-1.0)**(iterm-1)
                SA1(iang)  = SA1(iang)   +   fac * (ra0*fpi(iang)  + rb0*ftau(iang))      
                
                T          = (-1.0)**iterm
                SA2(iang)  = SA2(iang)   +   fac * (ra0*ftau(iang) + rb0*fpi(iang) )
                
                if  ( iang /= iang2 )  then
                   SA1(iang2) = SA1(iang2)   +   fac * (ra0*fpi( iang)*P + rb0*ftau(iang)*T)
                   SA2(iang2) = SA2(iang2)   +   fac * (ra0*ftau(iang)*T + rb0*fpi( iang)*P)
                end if
             end do
             
             iterm   = iterm + 1
             r_iterm = real(iterm, kind=r2)
             
             do iang=1, nang
                fpi1(iang) = ((2.0*r_iterm-1.0) / (r_iterm-1.0))   *   mu(iang)  *  fpi(iang)
                fpi1(iang) = fpi1(iang)   -   r_iterm * fpi0(iang)/(r_iterm-1.0)
                fpi0(iang) = fpi(iang)
             end do
          else
             ! start value for the next terms
             iterm = 2
          end if

          ! ------------------------------------------------------------------------------------------
          ! 2., 3., ... num
          ! ------------------------------------------------------------------------------------------
          z = -1.0_r2

          do        
             an  = an + 2.0_r2
             an2 = an - 2.0_r2
             
             ! Bessel functions
             if(iu1 == iu0) then
                besY2 = an2 * ax * besY1 - besY0
             else
                besY2 = an2 * ax * besY1 - besY0 / factor
             end if
             if(dabs(besY2) > 1.0e+300_r2) then
                besY2 = besY2 / factor
                iu2   = iu1 + 1
             end if
             besJ2 = (w1 + besY2 * besJ1) / besY1
             
             ! Mie coefficients
             r_iterm = real(iterm,kind=r2)
             
             s   = ru(iterm) / ri + r_iterm * ax

             s1  = s * besJ2 / fact(iu2) - besJ1 / fact(iu1)
             s2  = s * besY2 * fact(iu2) - besY1 * fact(iu1)
             ra1 = s1 / (s1 - s3 * s2)                        ! coefficient a_n, (n=iterm)
             
             s   = ru(iterm) * ri + r_iterm * ax
             s1  = s * besJ2 / fact(iu2) - besJ1 / fact(iu1)
             s2  = s * besY2 * fact(iu2) - besY1 * fact(iu1)
             rb1 = s1 / (s1 - s3 * s2)                        ! coefficient b_n, (n=iterm)
             
             ! efficiency factors
             z  = -z
             rr = z * (r_iterm + 0.5_r2) * (ra1 - rb1)
             r  = r + rr
             ss = ss + (r_iterm - 1.0_r2) * (r_iterm + 1.0_r2) / r_iterm * (ra0 * conjg(ra1)  &
                  + rb0 * conjg(rb1)) &
                  + an2 / r_iterm / (r_iterm - 1.0_r2) * (ra0 * conjg(rb0))
             qq   = an * (ra1 + rb1)
             Qext = Qext + qq
             Qsca = Qsca + an * (ra1 * conjg(ra1) + rb1 * conjg(rb1))
             
             ! leaving-the-loop criterion
             if ( dabs(qq / qext) < eps ) then
                exit
             end if
             
             ! Bessel functions
             besJ0 = besJ1
             besJ1 = besJ2
             besY0 = besY1
             besY1 = besY2
             iu0   = iu1
             iu1   = iu2
             ra0   = ra1
             rb0   = rb1
             
             ! scattering amplitude functions
             if (doSA) then
                r_iterm = real(iterm,kind=r2)       
                fac      = (2.0 * r_iterm+1.0) / (r_iterm * (r_iterm+1.0))
                
                do iang=1, nang
                   iang2      = 2 * nang - iang
                   
                   fpi(iang)  = fpi1(iang)
                   ftau(iang) = r_iterm * mu(iang) * fpi(iang)  -  (r_iterm+1.0) * fpi0(iang)
                   
                   P          = (-1.0)**(iterm-1)
                   SA1(iang)  = SA1(iang)   +   fac * (ra0*fpi(iang) + rb0*ftau(iang))      
                   
                   T          = (-1.0)**iterm
                   SA2(iang)  = SA2(iang)   +   fac * (ra0*ftau(iang) + rb0*fpi(iang))
                   
                   if  ( iang /= iang2 ) then
                      SA1(iang2) = SA1(iang2)   +   fac * (ra0*fpi(iang)*P  + rb0*ftau(iang)*T)
                      SA2(iang2) = SA2(iang2)   +   fac * (ra0*ftau(iang)*T + rb0*fpi( iang)*P)
                   end if
                end do
                
                iterm   = iterm + 1
                r_iterm = real(iterm,kind=r2)
                
                do iang=1, nang
                   fpi1(iang) = ((2.0*r_iterm-1.0) / (r_iterm-1.0))   *   mu(iang)  *  fpi(iang)
                   fpi1(iang) = fpi1(iang)   -   r_iterm * fpi0(iang)/(r_iterm-1.0)
                   fpi0(iang) = fpi(iang)
                end do
             else
                iterm = iterm + 1
             endif

             if ( iterm==num ) then
                exit
             else
                cycle
             end if
          end do
          
          ! efficiency factors (final calculations)
          Qext   = b * Qext
          Qsca   = b * Qsca
          Qbk    = 2.0_r2 * b * r * conjg(r)
          Qpr    = Qext - 2.0_r2 * b * ss
          Qabs   = Qext - Qsca
          albedo = Qsca / Qext
          g      = (Qext - Qpr) / Qsca
       end if
    end if
    deallocate( ru, mu, fpi, fpi0, fpi1, ftau )
  end subroutine shexqnn2
end module mie_routines


! ====================================================================================================
! Main Program
! ====================================================================================================
program mie
  use datatype
  use mie_routines

  implicit none

  logical :: doSA, svsep
  integer :: icomp, ilam, nlam, ncomp, nrad, ask1, ask2, ask3, nang, nang2, irad, ier, iang

  real(kind=r2) :: radmin, radmax, alpha, radminlog, radmaxlog, steplog, refmed, rad, rad1, delrad, &
       pi, x, qextx, qscax, qabsx, qbkx, qprx, albedox, gscattx, wqsc, wqscx, weisum, wrad, weight, &
       wradx, angx
  real(kind=r2), allocatable, dimension(:) :: lambda, abun, albedo, gscatt, S11x, S12x, S33x, S34x, &
       qext, qsca, qabs, qbk, qpr, &
       cext, csca, cabs, cbk, cpr
  real(kind=r2), allocatable, dimension(:,:) :: n, k, S11, S12, S33, S34

  complex(kind=r2) :: ri
  complex(kind=r2), allocatable, dimension(:) :: S1x, S2x

  character(len=8) :: fresult
  character(len=8), allocatable, dimension(:) :: fname
  ! --------------------------------------------------------------------------------------------------
  print *
  print *, "=========================================================================="
  print *, "                                                                          "
  print *, "     MIE SCATTERING     MIEX V1.04                                        "
  print *, "                                                _________________________ "
  print *, "                                                Contact:    swolf@mpia.de "
  print *, "=========================================================================="
  print *

  !---------------------------------------------------------------------------------------------------
  ! 0. General settings
  !---------------------------------------------------------------------------------------------------
  ! Real refractive index of the surrouding medium
  print *, "Real refractive index of the surrouding medium    : "; read *, refmed


  !---------------------------------------------------------------------------------------------------
  ! 1. Get main parameters
  !---------------------------------------------------------------------------------------------------
  print *, "Number of wavelengths                         : "; read *, nlam
  print *, "Number of chemical components                 : "; read *, ncomp

  allocate( fname(1:ncomp), lambda(1:nlam), n(1:ncomp,1:nlam), k(1:ncomp,1:nlam), abun(1:ncomp), &
       qext(1:nlam), qsca(1:nlam), qabs(1:nlam), qbk(1:nlam), qpr(1:nlam), &
       cext(1:nlam), csca(1:nlam), cabs(1:nlam), cbk(1:nlam), cpr(1:nlam), &
       albedo(1:nlam), gscatt(1:nlam) )

  lambda(:) = 0.0_r2
  n(:,:)    = 0.0_r2
  k(:,:)    = 0.0_r2
  abun(:)   = 0.0_r2

  qext(:)   = 0.0_r2
  qsca(:)   = 0.0_r2
  qabs(:)   = 0.0_r2
  qbk(:)    = 0.0_r2
  qpr(:)    = 0.0_r2

  cext(:)   = 0.0_r2
  csca(:)   = 0.0_r2
  cabs(:)   = 0.0_r2
  cbk(:)    = 0.0_r2
  cpr(:)    = 0.0_r2

  albedo(:) = 0.0_r2
  gscatt(:) = 0.0_r2

  !---------------------------------------------------------------------------------------------------
  print *, "Name of the dust data files (lambda/n/k data); 8 characters"
  print *, " [all data files have to contain the refractive]"
  print *, " [index the same wavelength distribution       ]"
  do icomp=1,ncomp
     print *, "    ", icomp, ". component : "; read  *, fname(icomp)
  end do

  if (ncomp>1) then
     print *, "<?> Relative abundances of the different components [%]"
     do icomp=1,ncomp
        print *, "    ", icomp, ". component : "; read  *, abun(icomp)
     end do
  else
     abun(:) = 100.0_r2
  end if
  abun(:) = abun(:)/100.0_r2

  print *, "-1- Single grain size"
  print *, "-2- Grain size distribution                   : "
  read *, ask1
  if (ask1==1) then
     print *, "    Grain radius [micron]                     : "
     read  *,radmin
     radmax = radmin
     nrad   = 1
     alpha  = 0.0_r2
  else
     print *, "    Minimum grain size [micron]               : "
     read  *,radmin
     print *, "    Maximum grain size [micron]               : "
     read  *,radmax
     print *,              "    Size distribution exponent                  "
     print *, "    alpha [n(a) .propto. a^alpha]             : "
     read  *,alpha
     print *, "    Number of size bins                       : "
     read  *,nrad
  end if

  print *, "Calculate scattering matrix elements (0=n/1=y): "
  read *, ask2
  if (ask2==1) then
     print *, "Number of scattering angles in the interval"
     print *, " [0°,180°];  odd number!"
     print *, " [example: '181' -> step width = 1°]          : "
     read *, nang2
     nang = (nang2-1)/2 + 1
     doSA = .true.
  else
     nang  = 1
     nang2 = 1
     doSA = .false.
  end if

  print *, "Project name (8 characters)                   : "
  read *,fresult
  
  print *, "Save results in separate files (0=n/1=y)      : "
  read *, ask3
  if (ask3==0) then 
     svsep = .false.
  else
     svsep = .true.
  end if


  !---------------------------------------------------------------------------------------------------
  ! 2. Read data files & Prepare the calculations
  !---------------------------------------------------------------------------------------------------
  print *
  print *,"<i> Calculation started ..."

  ! read lambda/n/k database
  do icomp=1, ncomp
     open(unit=1, file="./ri-data/"//fname(icomp), action="read", status="unknown", form="formatted")
     do ilam=1,nlam
        read(unit=1,fmt=*) lambda(ilam), n(icomp,ilam), k(icomp,ilam)
     end do
     close(unit=1)
  end do

  ! define radial step width
  radminlog = log10(radmin)
  radmaxlog = log10(radmax)
  if  (nrad>1) then
     steplog = (radmaxlog - radminlog) / real(nrad-1,kind=r2)
  else
     steplog = 0.0_r2
  endif

  pi = 4.0_r2 * atan(1.0_r2) ! PI = 3.14...
 
  allocate( S1x(1:nang2), S2x(1:nang2), &
       S11(1:nang2,1:nlam), S12(1:nang2,1:nlam), S33(1:nang2,1:nlam), S34(1:nang2,1:nlam), &
       S11x(1:nang2),       S12x(1:nang2),       S33x(1:nang2),       S34x(1:nang2) )
  
  S11(:,:) = 0.0_r2
  S12(:,:) = 0.0_r2
  S33(:,:) = 0.0_r2
  S34(:,:) = 0.0_r2

  !---------------------------------------------------------------------------------------------------
  ! 3. Run the Mie scattering routines
  !---------------------------------------------------------------------------------------------------
  do ilam=1, nlam
     weisum = 0.0_r2
     wrad   = 0.0_r2
     wqsc   = 0.0_r2

     do icomp=1, ncomp
        do irad=1, nrad
           ! initialize some arrays
           S1x(:)  = (0.0_r2,0.0_r2)
           S2x(:)  = (0.0_r2,0.0_r2)
           S11x(:) = 0.0_r2
           S12x(:) = 0.0_r2
           S33x(:) = 0.0_r2
           S34x(:) = 0.0_r2

           ! current radius / radius interval
           rad  = 10.0**(radminlog + (irad-1)*steplog) 
           rad1 = 10.0**(radminlog +  irad   *steplog)
           if  (nrad>1) then
              delrad = rad1 - rad
           else
              delrad = 1.0_r2
           endif
           
           ! size parameter
           x = 2.0_r2*pi * rad * refmed / lambda(ilam)

           ! complex refractive index
           ri = cmplx( n(icomp,ilam), k(icomp,ilam) )  /  refmed

           ! derive the scattering parameters
           call shexqnn2( ri, x, qextx, qscax, qabsx, qbkx, qprx, albedox, gscattx, ier, S1x, S2x, &
                doSA, nang )

           if (ier==1 .or. ier==2) then
              print *,"Program stopped."
              stop
           end if

           ! update average values
           weight = abun(icomp) * rad**alpha * delrad
           weisum = weisum + weight

           wradx  = pi*(rad/1.0e+6_r2)**2         * weight
           wqscx  = pi*(rad/1.0e+6_r2)**2 * qscax * weight

           wrad   = wrad + wradx
           wqsc   = wqsc + wqscx

           cext(  ilam)  =  cext(  ilam)  +  qextx * wradx
           csca(  ilam)  =  csca(  ilam)  +  qscax * wradx
           cbk(   ilam)  =  cbk(   ilam)  +  qbkx  * wradx
           cabs(  ilam)  =  cabs(  ilam)  +  qabsx * wradx

           qext(  ilam)  =  qext(  ilam)  +  qextx * wradx
           qsca(  ilam)  =  qsca(  ilam)  +  qscax * wradx
           qbk(   ilam)  =  qbk(   ilam)  +  qbkx  * wradx
           qabs(  ilam)  =  qabs(  ilam)  +  qabsx * wradx

           gscatt(ilam)  =  gscatt(ilam)  +  gscattx * wqscx

           S11x(:) =           0.5_r2 * abs(S2x(:)) * abs(S2x(:))
           S11x(:) = S11x(:) + 0.5_r2 * abs(S1x(:)) * abs(S1x(:))
           S12x(:) =           0.5_r2 * abs(S2x(:)) * abs(S2x(:))
           S12x(:) = S12x(:) - 0.5_r2 * abs(S1x(:)) * abs(S1x(:))
           S33x(:) = real(  S2x(:) * conjg(S1x(:)), kind=r2 )
           S34x(:) = aimag( S2x(:) * conjg(S1x(:)) )

           S11(:,ilam) = S11(:,ilam) + S11x(:)*weight
           S12(:,ilam) = S12(:,ilam) + S12x(:)*weight
           S33(:,ilam) = S33(:,ilam) + S33x(:)*weight
           S34(:,ilam) = S34(:,ilam) + S34x(:)*weight
        end do
     end do

     cext(  ilam) = cext(  ilam) / weisum
     csca(  ilam) = csca(  ilam) / weisum
     cbk(   ilam) = cbk(   ilam) / weisum
     cabs(  ilam) = cabs(  ilam) / weisum

     qext(  ilam) = qext(  ilam) / wrad
     qsca(  ilam) = qsca(  ilam) / wrad
     qbk(   ilam) = qbk(   ilam) / wrad
     qabs(  ilam) = qabs(  ilam) / wrad

     S11( :,ilam) = S11( :,ilam) / weisum
     S12( :,ilam) = S12( :,ilam) / weisum
     S33( :,ilam) = S33( :,ilam) / weisum
     S34( :,ilam) = S34( :,ilam) / weisum

     albedo(ilam) = csca(ilam)/cext(ilam)
     gscatt(ilam) = gscatt(ilam) / wqsc
  end do


  !---------------------------------------------------------------------------------------------------
  ! 4. Save the results
  !---------------------------------------------------------------------------------------------------
  open(unit=1, file="./results/"//fresult, action="write", status="unknown", form="formatted")

  write(unit=1,fmt=*) "# *** PROJECT PARAMETERS ***"
  write(unit=1,fmt=*)
  write(unit=1,fmt=*) "# Number of wavelengths            : ", nlam
  write(unit=1,fmt=*) "# Number of chemical components    : ", ncomp
  write(unit=1,fmt=*) "# Relative abundances [%]          : "
  do icomp=1,ncomp
     write(unit=1,fmt=*) "# ", icomp, ". component: ", real(abun(icomp)*100.0_r2,kind=r1)
  end do
  write(unit=1,fmt=*) "# Name(s) of the dust data file(s) :"
  do icomp=1,ncomp
     write(unit=1,fmt=*) "# ", icomp, ". component: ", fname(icomp)
  end do
  write(unit=1,fmt=*) "# Minimum grain radius [micron]    : ", real(radmin,kind=r1)
  write(unit=1,fmt=*) "# Maximum grain radius [micron]    : ", real(radmax,kind=r1)
  write(unit=1,fmt=*) "# Size distribution exponent       : ", real( alpha,kind=r1)
  write(unit=1,fmt=*) "# Number of size bins              : ", nrad
  if (doSA) then
     write(unit=1,fmt=*) "# Number of scattering angles      : ", nang2
  end if
  write(unit=1,fmt=*) 
  write(unit=1,fmt=*) 
  !---------------------------------------------------------------------------------------------------
  write(unit=1,fmt=*) "# *** RESULTS ***"
  write(unit=1,fmt=*)
  write(unit=1,fmt=*) "# 1. Wavelength [micron], Extinction efficiency factor / cross section [m^2]"
  do ilam=1,nlam
      write(unit=1,fmt=*) lambda(ilam), qext(ilam), cext(ilam)
  end do
  write(unit=1,fmt=*)

  write(unit=1,fmt=*) "# 2. Wavelength [micron], Scattering efficiency factor / cross section [m^2]"
  do ilam=1,nlam
      write(unit=1,fmt=*) lambda(ilam), qsca(ilam), csca(ilam)
  end do
  write(unit=1,fmt=*)

  write(unit=1,fmt=*) "# 3. Wavelength [micron], Backscattering efficiency factor/cross section [m^2]"
  do ilam=1,nlam
      write(unit=1,fmt=*) lambda(ilam), qbk(ilam), cbk(ilam)
  end do
  write(unit=1,fmt=*)

  write(unit=1,fmt=*) "# 4. Wavelength [micron], Absorption efficiency factor / cross section [m^2]"
  do ilam=1,nlam
      write(unit=1,fmt=*) lambda(ilam), qabs(ilam), cabs(ilam)
  end do
  write(unit=1,fmt=*)

  write(unit=1,fmt=*) "# 5. Wavelength [micron], Albedo"
  do ilam=1,nlam
      write(unit=1,fmt=*) lambda(ilam), albedo(ilam)
  end do
  write(unit=1,fmt=*)

  write(unit=1,fmt=*) "# 6. Wavelength [micron], Scattering asymmetry factor g"
  do ilam=1,nlam
      write(unit=1,fmt=*) lambda(ilam), gscatt(ilam)
  end do
  write(unit=1,fmt=*)

  write(unit=1,fmt=*) "# 7. Radiation pressure efficiency factor Qpr"
  do ilam=1,nlam
     write(unit=1,fmt=*) lambda(ilam), qext(ilam) - gscatt(ilam)*qsca(ilam)
  end do
  write(unit=1,fmt=*)

  write(unit=1,fmt=*) "# 8. Wavelength [micron], theta [degree], F11-F12-F33-F34 "
  if (doSA) then
     do ilam=1,nlam
        do iang=1,nang2
           ! angx: scattering angle [°]
           angx = real(iang-1,kind=r2) * 180.0_r2/real(nang2-1,kind=r2)
           write(unit=1,fmt=*) lambda(ilam), angx, S11(iang,ilam)
           write(unit=1,fmt=*) lambda(ilam), angx, S12(iang,ilam)
           write(unit=1,fmt=*) lambda(ilam), angx, S33(iang,ilam)
           write(unit=1,fmt=*) lambda(ilam), angx, S34(iang,ilam)
        end do
     end do
  else
     write(unit=1,fmt=*) "# not calculated."
  end if
  close(unit=1)

  !---------------------------------------------------------------------------------------------------
  if (svsep) then
     ! 1. Extinktion efficiency factor
     open(unit=1, file="./results/"//fresult//".qext", action="write", status="unknown", &
          form="formatted")
     do ilam=1,nlam
        write(unit=1,fmt=*) lambda(ilam), qext(ilam)
     end do
     close(unit=1)
     
     ! 2. Extinktion cross section [m^2]
     open(unit=1, file="./results/"//fresult//".cext", action="write", status="unknown", &
          form="formatted")
     do ilam=1,nlam
        write(unit=1,fmt=*) lambda(ilam), cext(ilam)
     end do
     close(unit=1)
     
     ! 3. Scattering efficiency factor
     open(unit=1, file="./results/"//fresult//".qsca", action="write", status="unknown", &
          form="formatted")
     do ilam=1,nlam
        write(unit=1,fmt=*) lambda(ilam), qsca(ilam)
     end do
     close(unit=1)

     ! 4. Scattering cross section [m^2]
     open(unit=1, file="./results/"//fresult//".csca", action="write", status="unknown", &
          form="formatted")
     do ilam=1,nlam
        write(unit=1,fmt=*) lambda(ilam), csca(ilam)
     end do
     close(unit=1)
     
     ! 5. Backscattering effiency factor
     open(unit=1, file="./results/"//fresult//".qbk", action="write", status="unknown", &
          form="formatted")
     do ilam=1,nlam
        write(unit=1,fmt=*) lambda(ilam), qbk(ilam)
     end do
     close(unit=1)
     
     ! 6. Backscattering cross section [m^2]
     open(unit=1, file="./results/"//fresult//".cbk", action="write", status="unknown", &
          form="formatted")
     do ilam=1,nlam
        write(unit=1,fmt=*) lambda(ilam), cbk(ilam)
     end do
     close(unit=1)
     
     ! 7. Absorption efficieny factor
     open(unit=1, file="./results/"//fresult//".qabs", action="write", status="unknown", &
          form="formatted")
     do ilam=1,nlam
        write(unit=1,fmt=*) lambda(ilam), qabs(ilam)
     end do
     close(unit=1)
     
     ! 8. Absorption cross section [m^2]
     open(unit=1, file="./results/"//fresult//".cabs", action="write", status="unknown", &
          form="formatted")
     do ilam=1,nlam
        write(unit=1,fmt=*) lambda(ilam), cabs(ilam)
     end do
     close(unit=1)
     
     ! 9. Albedo
     open(unit=1, file="./results/"//fresult//".alb", action="write", status="unknown", &
          form="formatted")
     do ilam=1,nlam
        write(unit=1,fmt=*) lambda(ilam), albedo(ilam)
     end do
     close(unit=1)
     
     ! 10. Scattering asymmetry factor g
     open(unit=1, file="./results/"//fresult//".g", action="write", status="unknown", &
          form="formatted")
     do ilam=1,nlam
        write(unit=1,fmt=*) lambda(ilam), gscatt(ilam)
     end do
     close(unit=1)
     
     ! 11. Qpr
     open(unit=1, file="./results/"//fresult//".qpr", action="write", status="unknown", &
          form="formatted")
     do ilam=1,nlam
        write(unit=1,fmt=*) lambda(ilam), qext(ilam) - gscatt(ilam)*qsca(ilam)
     end do
     close(unit=1)

     ! 12. Scattering Matrix elements S11, S12, S33, and S34
     if (doSA) then
        ! 11.1. S11
        open(unit=1, file="./results/"//fresult//".f11", action="write", status="unknown", &
             form="formatted")
        do ilam=1,nlam
           write(unit=1,fmt=*) lambda(ilam)
           do iang=1,nang2
              ! angx: scattering angle [°]
              angx = real(iang-1,kind=r2) * 180.0_r2/real(nang2-1,kind=r2)
              write(unit=1,fmt=*) angx, S11(iang,ilam)
           end do
        end do
        close(unit=1)

        ! 11.2. S12
        open(unit=1, file="./results/"//fresult//".f12", action="write", status="unknown", &
             form="formatted")
        do ilam=1,nlam
           write(unit=1,fmt=*) lambda(ilam)
           do iang=1,nang2
              ! angx: scattering angle [°]
              angx = real(iang-1,kind=r2) * 180.0_r2/real(nang2-1,kind=r2)
              write(unit=1,fmt=*) angx, S12(iang,ilam)
            end do
        end do
        close(unit=1)

        ! 11.3. S33
        open(unit=1, file="./results/"//fresult//".f33", action="write", status="unknown", &
             form="formatted")
        do ilam=1,nlam
           write(unit=1,fmt=*) lambda(ilam)
           do iang=1,nang2
              ! angx: scattering angle [°]
              angx = real(iang-1,kind=r2) * 180.0_r2/real(nang2-1,kind=r2)
              write(unit=1,fmt=*) angx, S33(iang,ilam)
           end do
        end do
        close(unit=1)

        ! 11.4. S34
        open(unit=1, file="./results/"//fresult//".f34", action="write", status="unknown", &
             form="formatted")
        do ilam=1,nlam
           write(unit=1,fmt=*) lambda(ilam)
           do iang=1,nang2
              ! angx: scattering angle [°]
              angx = real(iang-1,kind=r2) * 180.0_r2/real(nang2-1,kind=r2)
              write(unit=1,fmt=*) angx, S34(iang,ilam)
           end do
        end do
        close(unit=1)
     end if
  end if

  !---------------------------------------------------------------------------------------------------
  ! 5. Clean-up
  !---------------------------------------------------------------------------------------------------
  deallocate( fname, lambda, n, k, abun, qext, qsca, qabs, qbk, qpr, cext, csca, cabs, cbk, cpr, &
       albedo, gscatt, S1x, S2x, S11, S12, S33, S34, S11x, S12x, S33x, S34x)
  
  print *,"    ... done."
  print *

end program mie
