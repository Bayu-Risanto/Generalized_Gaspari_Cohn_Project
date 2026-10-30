! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module cov_cutoff_mod


use     types_mod, only : r8, PI, earth_radius
use utilities_mod, only : error_handler, E_ERR, E_MSG, &
                          do_output, do_nml_file, do_nml_term, nmlfileunit, &
                          find_namelist_in_file, check_namelist_read
use location_mod,  only : location_type, get_location
use gengc_function, only : gengc, f1, f2, f3, f4, f5, f6, f7

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

!============================================================================

contains

!======================================================================



function comp_cov_factor(z_in, c, obs_loc, obs_type, target_loc, target_kind, &
   separ_km, separ_idx, localization_override)
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
! CBR: separ_km is the horizontal distance between the obs and the state vars
! as calculated in assim_tools_mod.f90. It is an input array here. separ_idx is
! the index j in the loop STATE_UPDATE and OBS_UPDATE.

implicit none

real(r8),                      intent(in) :: z_in, c
type(location_type), optional, intent(in) :: obs_loc, target_loc
integer,             optional, intent(in) :: obs_type, target_kind
integer,             optional, intent(in) :: localization_override  
real(r8)                                  :: comp_cov_factor

real(r8) :: z, r, z_top, z_cen, rateredc1 
real(r8) :: BETA1, BETA2 
real(r8), allocatable :: avals(:), cvals(:)
real(r8),           optional, intent(in) :: separ_km(:)
integer,            optional, intent(in) :: separ_idx
real(r8), allocatable :: out(:)
integer  :: idx, n_prescb
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
      if(.not.present(separ_km)) then
           call error_handler(E_ERR,'comp_cov_factor', &
                   'Generalized Gaspari Cohn (select_localization = 4) requires separ_km, ' // &
              'but it was not passed by the caller', source)
      endif     

   !! CBR: define ak,al,ck,cl here
   ! ***********Testing GeGC with a varying and c varying **********************
   ! ***************************************************************************
   ! TEST 8 ------ testing_forecastQ_a_c_gengc_07.ipynb--------- what if we call my_state_loc(:) in radians
   z_top  = 170.0_r8      ! the top level where impact goes to zero
   z_cen  = 100_r8        ! central point where value goes to zero
   rateredc1 = 7.5_r8     ! reduction rate from unity to 0.8 radians
   BETA1  = 0.06_r8       ! curvature of top end (avals)
   BETA2  = 0.06_r8       ! curvature of top end (cvals)

   ! calculate the distance in degrees using my_state_loc(:) and my_obs_log(:) in the assim_tools_mod.f90

   ! convert degrees to distance in km since earth_radius is in km -- in the assim_tools_mod.f90

   ! get length of separation in km. And call it here. Maybe horizontal only is ok because the length of lat, lon, vloc
   ! in my_state_loc is the same. I am not sure yet.
   n_prescb = size(separ_km)

   ! --- avals ---
   avals = (1.0_r8 / 1000.0_r8) * tanh(BETA1 * (separ_km - z_cen))**2
   
   ! --- cvals ---
   cvals = (0.4_r8 * z_top - 100.0_r8) * tanh(BETA2 * (separ_km - z_cen)) * rateredc1 + 1.5_r8

   !! call subroutine gengc.f90. It is to create a profile of localization in meters
   call gengc(separ_km * 1000.0_r8, avals(1) * 1000.0_r8, avals * 1000.0_r8, cvals(1) * 1000.0_r8, cvals * 1000.0_r8, n_prescb, out) 
   
   ! determine covariance factors based on the profile. Note: 0 to PI is the farthest distance on a great circle.
   if (z < PI) then
      comp_cov_factor = out(separ_idx)   
   else
      comp_cov_factor = 0.0_r8
   end if
   

else ! Otherwise namelist parameter is illegal; this is an error

     call error_handler(E_ERR,'comp_cov_factor', &
              'Illegal value of "localization" in cov_cutoff_mod namelist', source)

endif

end function comp_cov_factor

end module cov_cutoff_mod
