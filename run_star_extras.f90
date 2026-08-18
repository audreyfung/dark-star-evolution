! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff

      implicit none

      include 'test_suite_extras_def.inc'
      include 'xtra_coeff_os/xtra_coeff_os_def.inc'

      ! these routines are called by the standard run_star check_model

      contains


      include 'test_suite_extras.inc'
      include 'xtra_coeff_os/xtra_coeff_os.inc'


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         include 'xtra_coeff_os/xtra_coeff_os_controls.inc'
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         ! s% other_momentum => default_other_momentum
         s% other_momentum_implicit => default_other_momentum
      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 1
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         !  new column
         names(1) = 'new_G_at_core'
         vals(1) = s%extra_grav(s%nz)%val
         ierr = 0
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 3
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n), mass_ds_in, mass_ds_out
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         if (n /= 3) stop 'data_for_extra_profile_columns'
         names(1) = 'extra_gravity'
         do k = 1, nz
            vals(k,1) = s% extra_grav(k)%val
         end do
         names(2) = 'mass_DS'
         do k = 1, nz
            vals(k,2) = -s%extra_grav(k)%val * s%r(k)**2 / standard_cgrav
         end do
         names(3) = 'rho_DS'

         do k = 1, nz
            if (k .eq. nz) then 
               mass_ds_in = 0.d0 
               mass_ds_out = -s%extra_grav(k-1)%val * s%r(k-1)**2 / standard_cgrav
               vals(k,3) = (mass_ds_out - mass_ds_in)/(4*pi*s%r(k)**2* (s%r(k-1)))
            else if (k .eq. 1) then 
               vals(k,3) = 0.d0
            else 
               mass_ds_in = -s%extra_grav(k+1)%val * s%r(k+1)**2 / standard_cgrav
               mass_ds_out = -s%extra_grav(k-1)%val * s%r(k-1)**2 / standard_cgrav
               vals(k,3) = (mass_ds_out - mass_ds_in)/(4*pi*s%r(k)**2* (s%r(k-1) - s%r(k+1)))
            end if 
         end do
      end subroutine data_for_extra_profile_columns


      subroutine default_other_momentum(id, ierr)
         use dark_star_support
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         call compute_G_prime(id, ierr)

      end subroutine default_other_momentum

      ! returns either keep_going, retry, or terminate.
      integer function extras_finish_step(id)
         use chem_def
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step


      end module run_star_extras

