module dark_star_support

    use star_def
    use const_def
    implicit none

    real(dp), parameter :: target_mass = log10(5.d-10 * Msun)
    real(dp), parameter :: f_c_at_equil = log10(5.724206d-10)
    real(dp), parameter :: r_DS_at_equil = log10(1.594146d0*Rsun)
    real(dp), parameter :: f_floor = log10(1.d-14)
    real(dp), parameter :: kappa = 5.222835d-14
    
    real(dp), save :: previous_successful_f_c = -50



contains

    function get_dmdx(x, f, m_d) result(dmdx)
        real(dp), intent(in) :: x, f, m_d
        real(dp) :: dmdx
        dmdx = 4*pi * 10**(3*x+f-m_d)
    end function get_dmdx

    function get_dfdx(x, m_s, m_d, f) result(dfdx)
        real(dp), intent(in) :: x, m_s, m_d, f
        real(dp) :: dfdx

        dfdx = -standard_cgrav * (10**(m_s) + 10**(m_d))/(2 * kappa * 10**(x+f))
    end function get_dfdx

    function interp_mass(x, xp, fp) result(f)
        real(dp), intent(in) :: x, xp(:), fp(:)
        real(dp), allocatable :: fp_flip(:), xp_flip(:)
        real(dp) :: f, f1, f2, x1, x2
        integer :: len, idx
        

        len = size(xp)

        allocate(fp_flip(len))
        allocate(xp_flip(len))

        xp_flip = xp(len:1:-1)
        fp_flip = fp(len:1:-1)
        
        if (x <= xp_flip(len) .and. x >= xp_flip(1)) then 
            idx = minloc(abs(xp_flip - x), dim=1)
            if (xp_flip(idx) > x) then 
                x1 = xp_flip(idx-1)
                x2 = xp_flip(idx)

                f1 = fp_flip(idx-1)
                f2 = fp_flip(idx)

            elseif (xp_flip(idx) < x) then 
                x1 = xp_flip(idx)
                x2 = xp_flip(idx+1)

                f1 = fp_flip(idx)
                f2 = fp_flip(idx+1)
            elseif (xp_flip(idx) == x) then  
                f = fp_flip(idx)
                return 
            end if

            f = f1 + (f2-f1)/(x2-x1) * (x-x1)
            return
        elseif (x < xp_flip(1)) then 
            f = fp_flip(1) + 3*x - 3*xp_flip(1)
            return
        elseif (x > xp_flip(len)) then 
            f = fp_flip(len)
            return

        end if 
            

    end function interp_mass

    function interp_dark_rho(x, xp, fp) result(f)
        real(dp), intent(in) :: x, xp(:), fp(:)
        real(dp) :: f, f1, f2, x1, x2
        real(dp), allocatable :: xp_flip(:), fp_flip(:)
        integer :: len, idx
        
        len = size(xp)

        allocate(xp_flip(len))
        allocate(fp_flip(len))

        xp_flip = xp(len:1:-1)
        fp_flip = fp(len:1:-1)
        
        if (x <= xp_flip(len) .and. x >= xp_flip(1)) then 
            idx = minloc(abs(xp_flip - x), dim=1)
            if (xp_flip(idx) > x) then 
                x1 = xp_flip(idx-1)
                x2 = xp_flip(idx)

                f1 = fp_flip(idx-1)
                f2 = fp_flip(idx)

            elseif (xp_flip(idx) < x) then 
                x1 = xp_flip(idx)
                x2 = xp_flip(idx+1)

                f1 = fp_flip(idx)
                f2 = fp_flip(idx+1)
            elseif (xp_flip(idx) == x) then  
                f = fp_flip(idx)
                return 
            end if

            f = f1 + (f2-f1)/(x2-x1) * (x-x1)
            return
        elseif (x < xp_flip(1)) then 
            f = fp_flip(1)
            return
        elseif (x > xp_flip(len)) then 
            f = -50.0
            return

        end if 
            

    end function interp_dark_rho

    function interp_array(x, xp, yp, left, right) result(y)
        real(dp), intent(in) :: x(:)        ! A 1D array of points you want to evaluate
        real(dp), intent(in) :: xp(:)       ! Your reference X lookup grid
        real(dp), intent(in) :: yp(:)       ! Your reference Y lookup grid
        real(dp), intent(in), optional :: left
        real(dp), intent(in), optional :: right
        
        real(dp) :: y(size(x))              ! Output array matching the exact size of x
        
        integer :: n_xp, n_x, i, j
        real(dp) :: weight

        n_xp = size(xp)
        n_x  = size(x)

        do j = 1, n_x
            
            ! 1. Handle Out-of-Bounds (Left Side)
            if (x(j) < xp(1)) then
                if (present(left)) then
                    y(j) = left
                else
                    y(j) = yp(1)
                end if
                cycle
            end if

            ! 2. Handle Out-of-Bounds (Right Side)
            if (x(j) > xp(n_xp)) then
                if (present(right)) then
                    y(j) = right
                else
                    y(j) = yp(n_xp)
                end if
                cycle
            end if

            ! 3. Find the interval for x(j)
            i = 1
            do while (i < n_xp .and. xp(i+1) < x(j))
                i = i + 1
            end do

            ! 4. Linear Interpolation for x(j)
            weight = (x(j) - xp(i)) / (xp(i+1) - xp(i))
            y(j) = yp(i) + weight * (yp(i+1) - yp(i))

        end do

    end function interp_array

    subroutine print_data(data, filename, log_number)
        real(dp), intent(in) :: data(:)
        integer, intent(in) :: log_number
        integer :: k
        character(len=*), intent(in) :: filename
        character(len=3) :: log_string
        character(len=100) :: full_filename
        character(len=52) :: path
        
        write(log_string, '(I0)') log_number
        path = '/Users/audreyfung/mesa-HDM/dark_star_log/'
        full_filename = trim(path) // trim(filename) // trim(log_string) // '.txt'
        
        open(unit=10, file=full_filename, status='replace', action='write')
        do k = 1, size(data)
            write(10, *) data(k)
        end do

        close(10)

    end subroutine print_data

    subroutine print_DS_profile(data1, data2, data3, lognumber, len, extra_name)
        real(dp), intent(in) :: data1(:), data2(:), data3(:)
        integer, intent(in) :: lognumber, len
        character(len=*), intent(in) :: extra_name

        integer :: i, io_unit, ios
        character(len=10) :: log_string
        character(len=100) :: full_filename
        character(len=52) :: path
        

        path = '/Users/audreyfung/mesa-HDM/dark_star_log/'
        write(log_string, '(I0)') lognumber
        full_filename = trim(path) // 'DS_profile_mn' // extra_name // trim(log_string) // '.txt'

        open(newunit=io_unit, file=trim(full_filename), status='replace', action='write', iostat=ios)

        if (ios /= 0) then
            return 
        end if 

        write(io_unit, '(3A24)') &
            'r(cm)                   ', &
            'rho (g/cm^3)            ', &
            'm_enclosed (g)          '

        do i = 1, len
            write(io_unit, '(3(1PE24.16))') data1(i), data2(i), data3(i)
        end do
        flush(io_unit)
        close(io_unit)

    end subroutine print_DS_profile

    subroutine get_xrange(star_radius, x_range)
        real(dp), intent(in) :: star_radius(:)
        real(dp), intent(out) :: x_range(:)
        real(dp) :: xstart, xend, x_outer(50)
        integer :: i
        
        xstart = star_radius(1) + 5
        xend = star_radius(1)
        x_outer = [(xstart + (i-1) * (xend - xstart) / 50, i = 1, 50)]

        x_range = [x_outer, star_radius]

    end subroutine get_xrange

    subroutine get_DS_xrange(x_range)
        real(dp), intent(out) :: x_range(:)
        real(dp) :: xstart, xend
        integer :: i, len

        len = size(x_range)
        xstart = r_DS_at_equil + 2
        xend = -5
        x_range = [(xstart + (i-1) * (xend - xstart) / len, i = 1, len)]

    end subroutine get_DS_xrange

    subroutine integrate_profile(f_c, star_radius, star_enclosed_mass, record_arrays, m_profile, f_profile, &
                                final_mass, final_radius, x_range)
        real(dp), intent(in)      :: star_radius(:), star_enclosed_mass(:), f_c
        logical                   :: record_arrays
        real(dp), intent(out)     :: m_profile(:), f_profile(:), x_range(:), final_mass, final_radius
        integer  :: i, len
        real(dp) :: x_previous, x_k, x_mid, f_previous, f_k, f_mid, md_previous, md_k, md_mid, ms_mid
        real(dp), allocatable :: ms(:)

        call get_DS_xrange(x_range)
        len = size(x_range)
        allocate(ms(len))

        ! here this is mesa star info 
        ! writing a function above called interp_star
        do i = 1, len
            ms(i) = interp_mass(x_range(i), star_radius(:), star_enclosed_mass(:))
        end do

        if (record_arrays) then 
                m_profile(:) = -50
                f_profile(:) = -50
        end if
        
        do i = len, 2, -1
            if (i .eq. len) then 
                f_k = f_c
                md_k = log10(4./3. * pi) + 3 * x_range(i) +f_k
                f_previous = f_k
                md_previous = md_k
            end if

            if (i .ne. len) then 
                x_previous = x_range(i+1)
                x_k = x_range(i)
                x_mid = (x_previous + x_k)/2.
                ms_mid = (ms(i+1) + ms(i))/2.

                ! midpoint
                md_mid = md_previous + 0.5 * (x_k - x_previous) * get_dmdx(x_previous, f_previous, md_previous)
                f_mid = f_previous + 0.5 * (x_k - x_previous) * get_dfdx(x_previous, ms(i+1), md_previous, f_previous)

                ! update
                md_k = md_previous + (x_k - x_previous) * get_dmdx(x_mid, f_mid, md_mid)
                f_k = f_previous + (x_k - x_previous) * get_dfdx(x_mid, ms_mid, md_mid, f_mid)

                md_previous = md_k
                f_previous = f_k
                end if
            
            if (record_arrays) then 
                m_profile(i) = md_k
                f_profile(i) = f_k
            end if

            if (f_k <= f_floor) exit

        end do

        if (i .ge. 2) m_profile(:i-1) = md_k

        final_mass = md_k ! This is passed back to the root-finder
        final_radius = x_range(i)

    end subroutine integrate_profile

    subroutine compute_DS_profile(star_radius, star_enclosed_mass, m_profile, f_profile, x_range)
        
        real(dp) :: f_c_guess, mass_error, tolerance, dr_outer
        real(dp) :: old_f_c, old_mass_error, temp_f_c, final_mass, final_radius
        real(dp), intent(in) :: star_radius(:), star_enclosed_mass(:)
        real(dp), intent(out) :: m_profile(:), f_profile(:), x_range(:)
        integer  :: iteration, i
        logical  :: converged        

        tolerance = 1.0d-4
        converged = .false.

        ! Setup initial guess
        if (previous_successful_f_c == -50) then
            f_c_guess = f_c_at_equil
        else
            f_c_guess = previous_successful_f_c
        end if
        
        ! --- Shooting Loop ---
        do iteration = 1, 100
            ! Call the profile integrator (set record_arrays = .false. to save speed while guessing)
            call integrate_profile(f_c_guess, star_radius, star_enclosed_mass, .false., m_profile, & 
                                f_profile, final_mass, final_radius, x_range)
            
            mass_error = (final_mass - target_mass) / target_mass
            
            if (abs(mass_error) < tolerance) then
                converged = .true.
                previous_successful_f_c = f_c_guess
                exit
            end if
            
            ! Secant Root Finder Engine
            if (iteration == 1) then
                old_mass_error = mass_error
                old_f_c = f_c_guess
                f_c_guess = f_c_guess +  log10(1.5)
            else
                temp_f_c = f_c_guess
                f_c_guess = f_c_guess - mass_error * (f_c_guess - old_f_c) / (mass_error - old_mass_error)
                old_f_c = temp_f_c
                old_mass_error = mass_error
            end if
        end do
        if (.not. converged) stop "CRITICAL: Shooting method failed!"

        ! --- Final Storage Run ---
        ! Run it one last time with record_arrays = .true. to save the final profile
        call integrate_profile(previous_successful_f_c, star_radius, star_enclosed_mass, .true., m_profile, &
                                f_profile, final_mass, final_radius, x_range)

    end subroutine compute_DS_profile

    subroutine compute_G_prime(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type(star_info), pointer :: s
        real(dp), allocatable :: m_profile(:), f_profile(:), x_range(:), star_radius(:), star_enclosed_mass(:)
        real(dp), allocatable :: m_profile_std(:), f_profile_std(:), x_range_std(:)
        real(dp) :: g_extra, start_year, end_year
        integer :: i, len, trigger_model
        logical :: threshold

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr/= 0) return

        allocate(star_radius(s%nz))
        allocate(star_enclosed_mass(s%nz))
        star_radius = log10(s%r(:s%nz))
        star_enclosed_mass = log10(s%m_grav(:s%nz))

        len = s%nz

        allocate(m_profile(5000))
        allocate(f_profile(5000))
        allocate(x_range(5000))


        allocate(m_profile_std(len))
        allocate(f_profile_std(len))
        allocate(x_range_std(len))


        call compute_DS_profile(star_radius, star_enclosed_mass, m_profile, f_profile, x_range)

        do i = 1, s%nz
            m_profile_std(i) = interp_mass(star_radius(i), x_range, m_profile)
            f_profile_std(i) = interp_dark_rho(star_radius(i), x_range, f_profile)
            x_range_std(i) = star_radius(i)
        end do

        start_year = 1.d1
        end_year = 1.d9

        if ((s%star_age) < end_year) then 
            threshold = .true.
            do i = 1, s%nz
                ! s%cgrav(i) = standard_cgrav * (1.d0 + (10**m_profile_std(i) / s%m_grav(i)))
                g_extra = -standard_cgrav *(10**m_profile_std(i))*s%rho(i) / s%r(i)**2
                s% extra_grav(i)%val = -abs(g_extra) * abs((s%star_age - start_year)/(end_year-start_year))
                s% extra_grav(i)%d1Array(:) = 0.d0
                if (i .eq. 1) then 
                    s% extra_grav(i)%d1Array(i_lnR_00) = (-2.d0 * abs(g_extra) - standard_cgrav / s%r(i) * &
                    (10**m_profile_std(1) - 10**m_profile_std(2))/(x_range_std(1) - x_range_std(2))) &
                    * abs((s%star_age - start_year)/(end_year-start_year))
                else if (i .eq. s%nz) then 
                    s% extra_grav(i)%d1Array(i_lnR_00) = (-2.d0 * abs(g_extra) - standard_cgrav / s%r(i) * & 
                    (10**m_profile_std(i))/(x_range_std(i))) & 
                    * abs((s%star_age - start_year)/(end_year-start_year))
                else
                    s% extra_grav(i)%d1Array(i_lnR_00) = (-2.d0 * abs(g_extra) - standard_cgrav / s%r(i) * &
                     (10**m_profile_std(i-1) - 10**m_profile_std(i+1))/(x_range_std(i-1) - x_range_std(i+1))) &
                     * abs((s%star_age - start_year)/(end_year-start_year))
                end if 
            end do

        else if ((s%star_age) > end_year) then

            do i = 1, s%nz
                
                ! s%cgrav(i) = standard_cgrav * (1.d0 + (10**m_profile_std(i) / s%m_grav(i)))
                g_extra = -standard_cgrav *(10**m_profile_std(i))*s%rho(i) / s%r(i)**2
                s% extra_grav(i)%val = -abs(g_extra)
                s% extra_grav(i)%d1Array(:) = 0.d0
                if (i .eq. 1) then 
                    s% extra_grav(i)%d1Array(i_lnR_00) = -2.d0 * g_extra - standard_cgrav / s%r(i) * &
                    (10**m_profile_std(1) - 10**m_profile_std(2))/(x_range_std(1) - x_range_std(2))
                else if (i .eq. s%nz) then 
                    s% extra_grav(i)%d1Array(i_lnR_00) = -2.d0 * g_extra - standard_cgrav / s%r(i) * & 
                    (10**m_profile_std(i))/(x_range_std(i))
                else
                    s% extra_grav(i)%d1Array(i_lnR_00) = -2.d0 * g_extra - standard_cgrav / s%r(i) * &
                     (10**m_profile_std(i-1) - 10**m_profile_std(i+1))/(x_range_std(i-1) - x_range_std(i+1))
                end if 
            end do
        end if 
        
        ! call print_DS_profile(10**x_range, 10**f_profile, 10**m_profile, s%cgrav, s%model_number, s%nz)
        ! call print_DS_profile(10**x_range, 10**f_profile, 10**m_profile, s%cgrav, s%model_number, s%nz, ramp_factor)
        if (mod(s% model_number, 50) == 0) then 
            call print_DS_profile(10**x_range, 10**f_profile, 10**m_profile, s%model_number, 5000, '')
            call print_DS_profile(10**x_range_std, 10**f_profile_std, 10**m_profile_std, s%model_number, s%nz, '_mesagrid')
        end if 
        
    end subroutine compute_G_prime





end module dark_star_support