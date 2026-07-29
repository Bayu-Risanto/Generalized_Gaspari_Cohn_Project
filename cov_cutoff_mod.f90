! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module cov_cutoff_mod


use     types_mod, only : r8
use utilities_mod, only : error_handler, E_ERR, E_MSG, &
                          do_output, do_nml_file, do_nml_term, nmlfileunit, &
                          find_namelist_in_file, check_namelist_read
use location_mod,  only : location_type
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
   localization_override)
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

implicit none

real(r8),                      intent(in) :: z_in, c
type(location_type), optional, intent(in) :: obs_loc, target_loc
integer,             optional, intent(in) :: obs_type, target_kind
integer,             optional, intent(in) :: localization_override
real(r8)                                  :: comp_cov_factor

real(r8) :: z, z_ch, r, Re, z_top, z_cen, z_trans 
real(r8) :: BETA1, BETA2, dist, cbtm
real(r8), allocatable :: avals(:), cvals(:), mz(:)
real(r8), allocatable :: out(:), z_arr(:)
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
   ! test 21 ------Good for TUS and FGZ ---------
   Re     = 6378000.0_r8   ! Earth Radius (meters)
   
   z_top  = 17.0_r8      ! the top level where impact goes to zero
   z_cen  = 7.2_r8
   z_trans = 10.0_r8
   BETA1  = 0.1_r8
   BETA2  = 0.8_r8
   
   ! --- Derived variables (declare these with appropriate dimensions) ---
   
   dist     = ABS(z_in)
   z_arr    = 2.0_r8 * Re * SIN(dist / 2.0_r8)   ! Convert great circle distance to chordal
   
   ! --- Mean of ens_z along first axis (axis=1 in Python = dim=1 in Fortran) ---
   ! ens_z shape assumed: (nz, ens_size) NEED TO CALL THE COMPLETE COLUMN OF ALTITUDE!!!
   ! for now, let's do hardcoded in km. Then change to meter in the call below. 

   mz = (/&
      0.73481584_r8,  0.7611376_r8,   0.8209757_r8,   0.8969528_r8,   0.99284023_r8, &
      1.1129234_r8,   1.2618954_r8,   1.444625_r8,    1.6657923_r8,   1.9294372_r8,  &
      2.2385476_r8,   2.5945466_r8,   2.9967809_r8,   3.4426503_r8,   3.9289868_r8,  &
      4.4563317_r8,   5.0289235_r8,   5.651413_r8,    6.313133_r8,    6.9845695_r8,  &
      7.6496086_r8,   8.308618_r8,    8.961394_r8,    9.607662_r8,   10.247328_r8,   &
     10.880244_r8,   11.505017_r8,   12.1181555_r8,   12.717808_r8,  13.304913_r8,    &
     13.880817_r8,   14.447325_r8,   15.006821_r8,   15.562005_r8,   16.115314_r8,    &
     16.668756_r8,   17.224102_r8,   17.78292_r8,    18.346317_r8,   18.914207_r8,    &
     19.48559_r8,    20.060015_r8,   20.638084_r8 /)

   cbtm = ABS(mz(1) - z_trans)
   ! --- avals ---
   avals(:) = -( (0.36_r8 * z_top - cbtm) / 20.0_r8 * TANH(BETA1 * (mz(:) - z_cen))**2 &
             - z_top / 50.0_r8 )
   ! --- cvals ---
   cvals(:) = ( (0.36_r8 * z_top - cbtm) * TANH(BETA2 * (mz(:) - z_cen)) * 8.0_r8 &
             + (1.4_r8 * z_top) )
  

   !! call subroutine gengc.f90 . Convert from km to meters . We need gengc to spit out one by one!
   call gengc(z_arr, avals(1) * 1000.0_r8, avals * 1000.0_r8, cvals(1) * 1000.0_r8, cvals * 1000.0_r8, 43, out) 
   comp_cov_factor = out(1)
   

else ! Otherwise namelist parameter is illegal; this is an error

     call error_handler(E_ERR,'comp_cov_factor', &
              'Illegal value of "localization" in cov_cutoff_mod namelist', source)

endif

end function comp_cov_factor

end module cov_cutoff_mod

