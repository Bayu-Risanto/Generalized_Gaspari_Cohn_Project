! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module cov_cutoff_mod


use     types_mod, only : r8, PI
use utilities_mod, only : error_handler, E_ERR, E_MSG, &
                          do_output, do_nml_file, do_nml_term, nmlfileunit, &
                          find_namelist_in_file, check_namelist_read
use location_mod,  only : location_type
use gengc_function, only : gengc, f1, f2, f3, f4, f5, f6, f7
use linspace_mod, only: linspace

implicit none
private

public :: comp_cov_factor

character(len=*), parameter :: source = 'cov_cutoff_mod.f90'


!============================================================================

!---- namelist with default values
logical :: namelist_initialized = .false.

integer :: select_localization = 4
! Value 1 selects default Gaspari-Cohn cutoff
! Value 2 selects boxcar
! Value 3 selects ramped boxcar
! Value 4 selects Generalized Gaspari Cohn

namelist / cov_cutoff_nml / select_localization

! Fallback point count used only if num_close is not supplied by the
! caller (keeps old behavior/callers working without modification).
integer, parameter :: default_n_prescb = 6371

!============================================================================

contains

!======================================================================



function comp_cov_factor(z_in, c, obs_loc, obs_type, target_loc, target_kind, &
   localization_override, num_close)
!----------------------------------------------------------------------
! function comp_cov_factor(z_in, c)
!
! Computes a covariance cutoff function from Gaspari and Cohn
! QJRMS, 125, 723-757.  (their eqn. 4.10)
!
! z_in is the distance while c is the cutoff distance. 
! For distances greater than 2c, the cov_factor returned goes to 0.

! Other ramping shapes are also available and can be selected by a namelist
! parameter. At present, these include a boxcar with the given halfwidth
! and a ramp in which the weight is set at 1.0 up to the half-width 
! distance and then decreases linearly to 0 at twice the half-width 
! distance.

! Additional information is passed in about the location and specific type of the
! observation and the location and generic kind of the variable being targeted for
! increments. These can be used for more refined algorithms that want to 
! make the cutoff a function of these additional arguments. 
!
! num_close is the scalar returned by get_close_obs() in assim_tools_mod.f90
! -- the number of close locations found for the current observation. When
! present, it is used (localization_type == 4 only) to set the number of
! points in the prescribed z profile, so the increment/spacing of
! z_prescb adapts to however many close locations this observation has,
! rather than using a fixed point count. When absent, default_n_prescb
! is used instead, preserving prior behavior for existing callers.

implicit none

real(r8),                      intent(in) :: z_in, c
type(location_type), optional, intent(in) :: obs_loc, target_loc
integer,             optional, intent(in) :: obs_type, target_kind
integer,             optional, intent(in) :: localization_override
integer,             optional, intent(in) :: num_close
real(r8)                                  :: comp_cov_factor

real(r8) :: z, r, Re, z_top, z_cen, rateredc1 
real(r8) :: BETA1, BETA2 
real(r8), allocatable :: avals(:), cvals(:), z_prescb(:)
real(r8), allocatable :: out(:)
integer  :: n_prescb, idx
integer  :: iunit, io
integer  :: localization_type

!--------------------------------------------------------
! Initialize namelist if not already done
if(.not. namelist_initialized) then


   namelist_initialized = .true.

   ! Read the namelist entry
   call find_namelist_in_file("input.nml", "cov_cutoff_nml", iunit)
   read(iunit, nml = cov_cutoff_nml, iostat = io)
   call check_namelist_read(iunit, io, "cov_cutoff_nml")

   if (do_nml_file()) write(nmlfileunit,nml=cov_cutoff_nml)
   if (do_nml_term()) write(     *     ,nml=cov_cutoff_nml)


   if (do_output()) then
      select case (select_localization)
         case (1)
            call error_handler(E_MSG,'comp_cov_factor:', &
               'Standard Gaspari Cohn localization selected')
         case (2)
            call error_handler(E_MSG,'comp_cov_factor:', &
               'Boxcar localization selected')
         case (3)
            call error_handler(E_MSG,'comp_cov_factor:', &
               'Ramped localization selected')
         case (4)
            call error_handler(E_MSG,'comp_cov_factor:', &
                'Generalized Gaspari Cohn localization selected')
         case default
            call error_handler(E_ERR,'comp_cov_factor', &
               'Illegal value of "select_localization" in cov_cutoff_mod namelist', source)
      end select
   endif

endif
!---------------------------------------------------------

if(present(localization_override)) then
   localization_type = localization_override
else
   localization_type = select_localization
endif

z = abs(z_in)

!----------------------------------------------------------

if(localization_type == 1) then ! Standard Gaspari Cohn localization

   if( z >= c*2.0_r8 ) then

      comp_cov_factor = 0.0_r8

   else if( z <= c ) then
      r = z / c
      comp_cov_factor = &
           ( ( ( -0.25_r8*r +0.5_r8 )*r +0.625_r8 )*r -5.0_r8/3.0_r8 )*r**2 + 1.0_r8
!!$           r**5 * (-0.25_r8 ) + &
!!$           r**4 / 2.0_r8 +              &
!!$           r**3 * 5.0_r8/8.0_r8 -       &
!!$           r**2 * 5.0_r8/3.0_r8 + 1.0_r8
   else

      r = z / c
      comp_cov_factor = &
           ( ( ( ( r/12.0_r8 -0.5_r8 )*r +0.625_r8 )*r +5.0_r8/3.0_r8 )*r -5.0_r8 )*r &
!!$           r**5 / 12.0_r8  -  &
!!$           r**4 / 2.0_r8   +  &
!!$           r**3 * 5.0_r8 / 8.0_r8 + &
!!$           r**2 * 5.0_r8 / 3.0_r8 - 5.0_r8*r &
           + 4.0_r8 - 2.0_r8 / (3.0_r8 * r) 
   endif

else if(localization_type == 2) then ! BOXCAR localization

   if(z < 2.0_r8 * c) then
      comp_cov_factor = 1.0_r8
   else
      comp_cov_factor = 0.0_r8
   endif

else if(localization_type == 3) then ! Ramped localization

   if(z >= 2.0_r8 * c) then
      comp_cov_factor = 0.0_r8
   else if(z >= c .and. z < 2.0_r8 * c) then
      comp_cov_factor = (2.0_r8 * c - z) / c
   else
      comp_cov_factor = 1.0_r8
   endif

else if(localization_type == 4) then ! Generalized Gaspari Cohn localization

   !! define ak,al,ck,cl here
   ! ***********Testing GeGC with a varying and c varying **********************
   ! ***************************************************************************
   ! TEST 6 ------ testing_forecastQ_a_c_gengc_06.ipynb---------
   Re     = 6371.0_r8   ! Earth Radius (in km)
   
   z_top  = 170.0_r8      ! the top level where impact goes to zero
   z_cen  = 100_r8        ! central point where value goes to zero
   rateredc1 = 7.5_r8     ! reduction rate from unity to 0.8 radians
   BETA1  = 0.06_r8       ! curvature of top end (avals)
   BETA2  = 0.06_r8       ! curvature of top end (cvals)

   ! Number of points in the prescribed z profile now tracks num_close,
   ! the count of close locations get_close_obs() found for this
   ! observation (assim_tools_mod.f90), so the spacing/increment of
   ! z_prescb adapts per observation. Falls back to a fixed default
   ! if num_close was not passed in (e.g. older call sites).
   if (present(num_close)) then
      n_prescb = num_close
   else
      n_prescb = default_n_prescb
   endif

   if (n_prescb <= 0) then
      ! no close locations for this observation -- nothing contributes
      comp_cov_factor = 0.0_r8
      return
   endif

   ! z_prescb sized/spaced by n_prescb (derived from num_close above)
   ! note the fartherst distance is half circumference = PI * Re
   call linspace(0.0_r8, PI * Re, n_prescb, z_prescb, endpoint=.false.)

   ! --- avals --- (allocatable arrays auto-(re)allocate on assignment
   ! to match the shape of z_prescb, so no explicit allocate() needed)
   avals = (1.0_r8 / 1000.0_r8) * tanh(BETA1 * (z_prescb - z_cen))**2
   
   ! --- cvals ---
   cvals = (0.4_r8 * z_top - 100.0_r8) * tanh(BETA2 * (z_prescb - z_cen)) * rateredc1 + 1.5_r8

   !! call subroutine gengc.f90 . Convert from km to meters . We need gengc to spit out one by one!
   call gengc(z_prescb, avals(1) * 1000.0_r8, avals * 1000.0_r8, cvals(1) * 1000.0_r8, cvals * 1000.0_r8, n_prescb, out) 
   
   if (z == 0.0_r8) then
      comp_cov_factor = out(1)                 ! alpha_d[0] -> alpha_d(1)
   else if (z > 0.0_r8 .and. z < PI * Re) then
      ! nearest-index lookup into the precomputed alpha_d profile
      idx = nint(z_in / (Re*1000.0_r8) * real(n_prescb - 1, r8)) + 1
      idx = max(1, min(n_prescb, idx))
      comp_cov_factor = out(idx)   !! need to be corrected!
   else
      comp_cov_factor = 0.0_r8
   end if
   

else ! Otherwise namelist parameter is illegal; this is an error

     call error_handler(E_ERR,'comp_cov_factor', &
              'Illegal value of "localization" in cov_cutoff_mod namelist', source)

endif

end function comp_cov_factor

end module cov_cutoff_mod
