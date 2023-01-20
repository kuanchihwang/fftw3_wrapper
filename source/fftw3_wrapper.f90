module fftw3_wrapper
    use, intrinsic :: iso_fortran_env
    use :: fftw3
    implicit none
    private
    save
    ! Forward fast Fourier transform family
    interface fft
        module procedure &
            dcfft1, drfft1, &
            dcfft2, drfft2, &
            dcfft3, drfft3, &
            dcfft4, drfft4
        module procedure &
            scfft1, srfft1, &
            scfft2, srfft2, &
            scfft3, srfft3, &
            scfft4, srfft4
    end interface fft
    ! Forward complex fast Fourier transform family
    interface cfft
        module procedure &
            dcfft1, &
            dcfft2, &
            dcfft3, &
            dcfft4
        module procedure &
            scfft1, &
            scfft2, &
            scfft3, &
            scfft4
    end interface cfft
    ! Forward real fast Fourier transform family
    interface rfft
        module procedure &
            drfft1, &
            drfft2, &
            drfft3, &
            drfft4
        module procedure &
            srfft1, &
            srfft2, &
            srfft3, &
            srfft4
    end interface rfft
    ! Inverse fast Fourier transform family
    interface ifft
        module procedure &
            dicfft1, dirfft1, &
            dicfft2, dirfft2, &
            dicfft3, dirfft3, &
            dicfft4, dirfft4
        module procedure &
            sicfft1, sirfft1, &
            sicfft2, sirfft2, &
            sicfft3, sirfft3, &
            sicfft4, sirfft4
    end interface ifft
    ! Inverse complex fast Fourier transform family
    interface icfft
        module procedure &
            dicfft1, &
            dicfft2, &
            dicfft3, &
            dicfft4
        module procedure &
            sicfft1, &
            sicfft2, &
            sicfft3, &
            sicfft4
    end interface icfft
    ! Inverse real fast Fourier transform family
    interface irfft
        module procedure &
            dirfft1, &
            dirfft2, &
            dirfft3, &
            dirfft4
        module procedure &
            sirfft1, &
            sirfft2, &
            sirfft3, &
            sirfft4
    end interface irfft
    ! Forward zero-frequency shift family
    interface fftshift
        module procedure &
            dfftshift1, &
            dfftshift2, &
            dfftshift3, &
            dfftshift4
        module procedure &
            sfftshift1, &
            sfftshift2, &
            sfftshift3, &
            sfftshift4
    end interface fftshift
    ! Inverse zero-frequency shift family
    interface ifftshift
        module procedure &
            difftshift1, &
            difftshift2, &
            difftshift3, &
            difftshift4
        module procedure &
            sifftshift1, &
            sifftshift2, &
            sifftshift3, &
            sifftshift4
    end interface ifftshift
    ! Generic public procedures
    public :: &
        fft, cfft, rfft, fftshift, &
        ifft, icfft, irfft, ifftshift
    ! Specific public procedures
    public :: &
        dcfft1, drfft1, dicfft1, dirfft1, &
        dcfft2, drfft2, dicfft2, dirfft2, &
        dcfft3, drfft3, dicfft3, dirfft3, &
        dcfft4, drfft4, dicfft4, dirfft4
    public :: &
        scfft1, srfft1, sicfft1, sirfft1, &
        scfft2, srfft2, sicfft2, sirfft2, &
        scfft3, srfft3, sicfft3, sirfft3, &
        scfft4, srfft4, sicfft4, sirfft4
contains
subroutine dcfft1(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:)
    complex(real64), intent(in) :: inp(:)
    integer, optional, intent(in) :: num(1)
    integer :: &
        xsize, &
        xinpsize
    complex(c_double_complex), pointer :: &
        inpptr(:) => null(), &
        outptr(:) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
    end if
    fftwdata = fftw_alloc_complex(int(xsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize])
    call c_f_pointer(fftwdata, outptr, [xsize])
    fftwplan = fftw_plan_dft_1d( &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_forward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize) = cmplx( &
        real( inp(1:xinpsize)), &
        aimag(inp(1:xinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft(fftwplan, inpptr, outptr)
    out(:) = cmplx(real(outptr), aimag(outptr), real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dcfft1
subroutine drfft1(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:)
    real(real64), intent(in) :: inp(:)
    integer, optional, intent(in) :: num(1)
    integer :: &
        xsize, &
        xhalfsize, &
        xinpsize
    real(c_double), pointer :: &
        inpptr(:) => null()
    complex(c_double_complex), pointer :: &
        outptr(:) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xsize, size(inp, 1))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xsize
    end if
    ! For in-place, real-input and complex-output ("r2c") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftw_alloc_complex(int(xhalfsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xhalfsize * 2])
    call c_f_pointer(fftwdata, outptr, [xhalfsize])
    fftwplan = fftw_plan_dft_r2c_1d( &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! Carefully specify the range of assignments because "inpptr" contains extra padding.
    inpptr = real(0.0, c_double)
    inpptr( 1:xinpsize) = real( &
        inp(1:xinpsize), &
        c_double &
    )
    call fftw_execute_dft_r2c(fftwplan, inpptr, outptr)
    ! The output array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "outptr" only contains the first half by design.
    out(1:xhalfsize) = cmplx(real(outptr), aimag(outptr), real64)
    if (mod(xsize, 2) == 0) then
        ! Even number of elements.
        ! The element at "xhalfsize" is the Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1) = &
            conjg(out(2:(xhalfsize - 1)))
    else
        ! Odd number of elements.
        ! No Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1) = &
            conjg(out(2:(xhalfsize)))
    end if
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine drfft1
subroutine dfftshift1(out, inp)
    complex(real64), allocatable, intent(out) :: out(:)
    complex(real64), intent(in) :: inp(:)
    allocate(out, source=inp)
    out(:) = cshift(out, -size(out, 1) / 2, 1)
end subroutine dfftshift1
subroutine dicfft1(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:)
    complex(real64), intent(in) :: inp(:)
    integer, optional, intent(in) :: num(1)
    integer :: &
        xsize, &
        xinpsize
    complex(c_double_complex), pointer :: &
        inpptr(:) => null(), &
        outptr(:) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
    end if
    fftwdata = fftw_alloc_complex(int(xsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize])
    call c_f_pointer(fftwdata, outptr, [xsize])
    fftwplan = fftw_plan_dft_1d( &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_backward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize) = cmplx( &
        real( inp(1:xinpsize)), &
        aimag(inp(1:xinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft(fftwplan, inpptr, outptr)
    out(:) = cmplx(real(outptr), aimag(outptr), real64)
    out = out / real(xsize, real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dicfft1
subroutine dirfft1(out, inp, num)
    real(real64), allocatable, intent(out) :: out(:)
    complex(real64), intent(in) :: inp(:)
    integer, optional, intent(in) :: num(1)
    integer :: &
        xsize, &
        xhalfsize, &
        xinpsize
    real(c_double), pointer :: &
        outptr(:) => null()
    complex(c_double_complex), pointer :: &
        inpptr(:) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xhalfsize, size(inp, 1) / 2 + 1)
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xhalfsize
    end if
    ! For in-place, complex-input and real-output ("c2r") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftw_alloc_complex(int(xhalfsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize])
    fftwplan = fftw_plan_dft_c2r_1d( &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! The input array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "inpptr" only needs the first half by design.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize) = cmplx( &
        real( inp(1:xinpsize)), &
        aimag(inp(1:xinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft_c2r(fftwplan, inpptr, outptr)
    ! Carefully specify the range of assignments because "outptr" contains extra padding.
    out(:) = real(outptr(1:xsize), real64)
    out = out / real(xsize, real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dirfft1
subroutine difftshift1(out, inp)
    complex(real64), allocatable, intent(out) :: out(:)
    complex(real64), intent(in) :: inp(:)
    allocate(out, source=inp)
    out(:) = cshift(out, size(out, 1) / 2, 1)
end subroutine difftshift1
subroutine dcfft2(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:, :)
    complex(real64), intent(in) :: inp(:, :)
    integer, optional, intent(in) :: num(2)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize
    complex(c_double_complex), pointer :: &
        inpptr(:, :) => null(), &
        outptr(:, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
    end if
    fftwdata = fftw_alloc_complex(int(xsize * ysize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize])
    fftwplan = fftw_plan_dft_2d( &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_forward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize, 1:yinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft(fftwplan, inpptr, outptr)
    out(:, :) = cmplx(real(outptr), aimag(outptr), real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dcfft2
subroutine drfft2(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:, :)
    real(real64), intent(in) :: inp(:, :)
    integer, optional, intent(in) :: num(2)
    integer :: &
        i, &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize
    real(c_double), pointer :: &
        inpptr(:, :) => null()
    complex(c_double_complex), pointer :: &
        outptr(:, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
    end if
    ! For in-place, real-input and complex-output ("r2c") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftw_alloc_complex(int(xhalfsize * ysize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xhalfsize * 2, ysize])
    call c_f_pointer(fftwdata, outptr, [xhalfsize, ysize])
    fftwplan = fftw_plan_dft_r2c_2d( &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! Carefully specify the range of assignments because "inpptr" contains extra padding.
    inpptr = real(0.0, c_double)
    inpptr( 1:xinpsize, 1:yinpsize) = real( &
        inp(1:xinpsize, 1:yinpsize), &
        c_double &
    )
    call fftw_execute_dft_r2c(fftwplan, inpptr, outptr)
    ! The output array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "outptr" only contains the first half by design.
    out(1:xhalfsize, :) = cmplx(real(outptr), aimag(outptr), real64)
    if (mod(xsize, 2) == 0) then
        ! Even number of elements.
        ! The element at "xhalfsize" is the Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)]) = &
            conjg(out(2:(xhalfsize - 1), :))
    else
        ! Odd number of elements.
        ! No Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)]) = &
            conjg(out(2:(xhalfsize), :))
    end if
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine drfft2
subroutine dfftshift2(out, inp)
    complex(real64), allocatable, intent(out) :: out(:, :)
    complex(real64), intent(in) :: inp(:, :)
    allocate(out, source=inp)
    out(:, :) = cshift(out, -size(out, 1) / 2, 1)
    out(:, :) = cshift(out, -size(out, 2) / 2, 2)
end subroutine dfftshift2
subroutine dicfft2(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:, :)
    complex(real64), intent(in) :: inp(:, :)
    integer, optional, intent(in) :: num(2)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize
    complex(c_double_complex), pointer :: &
        inpptr(:, :) => null(), &
        outptr(:, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
    end if
    fftwdata = fftw_alloc_complex(int(xsize * ysize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize])
    fftwplan = fftw_plan_dft_2d( &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_backward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize, 1:yinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft(fftwplan, inpptr, outptr)
    out(:, :) = cmplx(real(outptr), aimag(outptr), real64)
    out = out / real(xsize * ysize, real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dicfft2
subroutine dirfft2(out, inp, num)
    real(real64), allocatable, intent(out) :: out(:, :)
    complex(real64), intent(in) :: inp(:, :)
    integer, optional, intent(in) :: num(2)
    integer :: &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize
    real(c_double), pointer :: &
        outptr(:, :) => null()
    complex(c_double_complex), pointer :: &
        inpptr(:, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xhalfsize, size(inp, 1) / 2 + 1)
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xhalfsize
        ysize = size(inp, 2)
        yinpsize = ysize
    end if
    ! For in-place, complex-input and real-output ("c2r") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftw_alloc_complex(int(xhalfsize * ysize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2, ysize])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize, ysize])
    fftwplan = fftw_plan_dft_c2r_2d( &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! The input array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "inpptr" only needs the first half by design.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize, 1:yinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft_c2r(fftwplan, inpptr, outptr)
    ! Carefully specify the range of assignments because "outptr" contains extra padding.
    out(:, :) = real(outptr(1:xsize, :), real64)
    out = out / real(xsize * ysize, real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dirfft2
subroutine difftshift2(out, inp)
    complex(real64), allocatable, intent(out) :: out(:, :)
    complex(real64), intent(in) :: inp(:, :)
    allocate(out, source=inp)
    out(:, :) = cshift(out, size(out, 1) / 2, 1)
    out(:, :) = cshift(out, size(out, 2) / 2, 2)
end subroutine difftshift2
subroutine dcfft3(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:, :, :)
    complex(real64), intent(in) :: inp(:, :, :)
    integer, optional, intent(in) :: num(3)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize
    complex(c_double_complex), pointer :: &
        inpptr(:, :, :) => null(), &
        outptr(:, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
    end if
    fftwdata = fftw_alloc_complex(int(xsize * ysize * zsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize, zsize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize, zsize])
    fftwplan = fftw_plan_dft_3d( &
        int(zsize, c_int), &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_forward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft(fftwplan, inpptr, outptr)
    out(:, :, :) = cmplx(real(outptr), aimag(outptr), real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dcfft3
subroutine drfft3(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:, :, :)
    real(real64), intent(in) :: inp(:, :, :)
    integer, optional, intent(in) :: num(3)
    integer :: &
        i, &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize
    real(c_double), pointer :: &
        inpptr(:, :, :) => null()
    complex(c_double_complex), pointer :: &
        outptr(:, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
    end if
    ! For in-place, real-input and complex-output ("r2c") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftw_alloc_complex(int(xhalfsize * ysize * zsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xhalfsize * 2, ysize, zsize])
    call c_f_pointer(fftwdata, outptr, [xhalfsize, ysize, zsize])
    fftwplan = fftw_plan_dft_r2c_3d( &
        int(zsize, c_int), &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! Carefully specify the range of assignments because "inpptr" contains extra padding.
    inpptr = real(0.0, c_double)
    inpptr( 1:xinpsize, 1:yinpsize, 1:zinpsize) = real( &
        inp(1:xinpsize, 1:yinpsize, 1:zinpsize), &
        c_double &
    )
    call fftw_execute_dft_r2c(fftwplan, inpptr, outptr)
    ! The output array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "outptr" only contains the first half by design.
    out(1:xhalfsize, :, :) = cmplx(real(outptr), aimag(outptr), real64)
    if (mod(xsize, 2) == 0) then
        ! Even number of elements.
        ! The element at "xhalfsize" is the Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)], &
            [1, (i, i = zsize, 2, -1)]) = &
            conjg(out(2:(xhalfsize - 1), :, :))
    else
        ! Odd number of elements.
        ! No Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)], &
            [1, (i, i = zsize, 2, -1)]) = &
            conjg(out(2:(xhalfsize), :, :))
    end if
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine drfft3
subroutine dfftshift3(out, inp)
    complex(real64), allocatable, intent(out) :: out(:, :, :)
    complex(real64), intent(in) :: inp(:, :, :)
    allocate(out, source=inp)
    out(:, :, :) = cshift(out, -size(out, 1) / 2, 1)
    out(:, :, :) = cshift(out, -size(out, 2) / 2, 2)
    out(:, :, :) = cshift(out, -size(out, 3) / 2, 3)
end subroutine dfftshift3
subroutine dicfft3(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:, :, :)
    complex(real64), intent(in) :: inp(:, :, :)
    integer, optional, intent(in) :: num(3)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize
    complex(c_double_complex), pointer :: &
        inpptr(:, :, :) => null(), &
        outptr(:, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
    end if
    fftwdata = fftw_alloc_complex(int(xsize * ysize * zsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize, zsize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize, zsize])
    fftwplan = fftw_plan_dft_3d( &
        int(zsize, c_int), &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_backward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft(fftwplan, inpptr, outptr)
    out(:, :, :) = cmplx(real(outptr), aimag(outptr), real64)
    out = out / real(xsize * ysize * zsize, real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dicfft3
subroutine dirfft3(out, inp, num)
    real(real64), allocatable, intent(out) :: out(:, :, :)
    complex(real64), intent(in) :: inp(:, :, :)
    integer, optional, intent(in) :: num(3)
    integer :: &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize
    real(c_double), pointer :: &
        outptr(:, :, :) => null()
    complex(c_double_complex), pointer :: &
        inpptr(:, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xhalfsize, size(inp, 1) / 2 + 1)
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xhalfsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
    end if
    ! For in-place, complex-input and real-output ("c2r") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftw_alloc_complex(int(xhalfsize * ysize * zsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2, ysize, zsize])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize, ysize, zsize])
    fftwplan = fftw_plan_dft_c2r_3d( &
        int(zsize, c_int), &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! The input array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "inpptr" only needs the first half by design.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft_c2r(fftwplan, inpptr, outptr)
    ! Carefully specify the range of assignments because "outptr" contains extra padding.
    out(:, :, :) = real(outptr(1:xsize, :, :), real64)
    out = out / real(xsize * ysize * zsize, real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dirfft3
subroutine difftshift3(out, inp)
    complex(real64), allocatable, intent(out) :: out(:, :, :)
    complex(real64), intent(in) :: inp(:, :, :)
    allocate(out, source=inp)
    out(:, :, :) = cshift(out, size(out, 1) / 2, 1)
    out(:, :, :) = cshift(out, size(out, 2) / 2, 2)
    out(:, :, :) = cshift(out, size(out, 3) / 2, 3)
end subroutine difftshift3
subroutine dcfft4(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:, :, :, :)
    complex(real64), intent(in) :: inp(:, :, :, :)
    integer, optional, intent(in) :: num(4)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize, &
        tsize, &
        tinpsize
    complex(c_double_complex), pointer :: &
        inpptr(:, :, :, :) => null(), &
        outptr(:, :, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
        tsize = num(4)
        tinpsize = min(tsize, size(inp, 4))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
        tsize = size(inp, 4)
        tinpsize = tsize
    end if
    fftwdata = fftw_alloc_complex(int(xsize * ysize * zsize * tsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize, zsize, tsize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize, zsize, tsize])
    fftwplan = fftw_plan_dft( &
        4_c_int, [ &
            int(tsize, c_int), &
            int(zsize, c_int), &
            int(ysize, c_int), &
            int(xsize, c_int) &
        ], &
        inpptr, &
        outptr, &
        fftw_forward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize, tsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft(fftwplan, inpptr, outptr)
    out(:, :, :, :) = cmplx(real(outptr), aimag(outptr), real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dcfft4
subroutine drfft4(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:, :, :, :)
    real(real64), intent(in) :: inp(:, :, :, :)
    integer, optional, intent(in) :: num(4)
    integer :: &
        i, &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize, &
        tsize, &
        tinpsize
    real(c_double), pointer :: &
        inpptr(:, :, :, :) => null()
    complex(c_double_complex), pointer :: &
        outptr(:, :, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
        tsize = num(4)
        tinpsize = min(tsize, size(inp, 4))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
        tsize = size(inp, 4)
        tinpsize = tsize
    end if
    ! For in-place, real-input and complex-output ("r2c") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftw_alloc_complex(int(xhalfsize * ysize * zsize * tsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xhalfsize * 2, ysize, zsize, tsize])
    call c_f_pointer(fftwdata, outptr, [xhalfsize, ysize, zsize, tsize])
    fftwplan = fftw_plan_dft_r2c( &
        4_c_int, [ &
            int(tsize, c_int), &
            int(zsize, c_int), &
            int(ysize, c_int), &
            int(xsize, c_int) &
        ], &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize, tsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! Carefully specify the range of assignments because "inpptr" contains extra padding.
    inpptr = real(0.0, c_double)
    inpptr( 1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize) = real( &
        inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize), &
        c_double &
    )
    call fftw_execute_dft_r2c(fftwplan, inpptr, outptr)
    ! The output array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "outptr" only contains the first half by design.
    out(1:xhalfsize, :, :, :) = cmplx(real(outptr), aimag(outptr), real64)
    if (mod(xsize, 2) == 0) then
        ! Even number of elements.
        ! The element at "xhalfsize" is the Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)], &
            [1, (i, i = zsize, 2, -1)], &
            [1, (i, i = tsize, 2, -1)]) = &
            conjg(out(2:(xhalfsize - 1), :, :, :))
    else
        ! Odd number of elements.
        ! No Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)], &
            [1, (i, i = zsize, 2, -1)], &
            [1, (i, i = tsize, 2, -1)]) = &
            conjg(out(2:(xhalfsize), :, :, :))
    end if
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine drfft4
subroutine dfftshift4(out, inp)
    complex(real64), allocatable, intent(out) :: out(:, :, :, :)
    complex(real64), intent(in) :: inp(:, :, :, :)
    allocate(out, source=inp)
    out(:, :, :, :) = cshift(out, -size(out, 1) / 2, 1)
    out(:, :, :, :) = cshift(out, -size(out, 2) / 2, 2)
    out(:, :, :, :) = cshift(out, -size(out, 3) / 2, 3)
    out(:, :, :, :) = cshift(out, -size(out, 4) / 2, 4)
end subroutine dfftshift4
subroutine dicfft4(out, inp, num)
    complex(real64), allocatable, intent(out) :: out(:, :, :, :)
    complex(real64), intent(in) :: inp(:, :, :, :)
    integer, optional, intent(in) :: num(4)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize, &
        tsize, &
        tinpsize
    complex(c_double_complex), pointer :: &
        inpptr(:, :, :, :) => null(), &
        outptr(:, :, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
        tsize = num(4)
        tinpsize = min(tsize, size(inp, 4))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
        tsize = size(inp, 4)
        tinpsize = tsize
    end if
    fftwdata = fftw_alloc_complex(int(xsize * ysize * zsize * tsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize, zsize, tsize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize, zsize, tsize])
    fftwplan = fftw_plan_dft( &
        4_c_int, [ &
            int(tsize, c_int), &
            int(zsize, c_int), &
            int(ysize, c_int), &
            int(xsize, c_int) &
        ], &
        inpptr, &
        outptr, &
        fftw_backward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize, tsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft(fftwplan, inpptr, outptr)
    out(:, :, :, :) = cmplx(real(outptr), aimag(outptr), real64)
    out = out / real(xsize * ysize * zsize * tsize, real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dicfft4
subroutine dirfft4(out, inp, num)
    real(real64), allocatable, intent(out) :: out(:, :, :, :)
    complex(real64), intent(in) :: inp(:, :, :, :)
    integer, optional, intent(in) :: num(4)
    integer :: &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize, &
        tsize, &
        tinpsize
    real(c_double), pointer :: &
        outptr(:, :, :, :) => null()
    complex(c_double_complex), pointer :: &
        inpptr(:, :, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xhalfsize, size(inp, 1) / 2 + 1)
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
        tsize = num(4)
        tinpsize = min(tsize, size(inp, 4))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xhalfsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
        tsize = size(inp, 4)
        tinpsize = tsize
    end if
    ! For in-place, complex-input and real-output ("c2r") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftw_alloc_complex(int(xhalfsize * ysize * zsize * tsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2, ysize, zsize, tsize])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize, ysize, zsize, tsize])
    fftwplan = fftw_plan_dft_c2r( &
        4_c_int, [ &
            int(tsize, c_int), &
            int(zsize, c_int), &
            int(ysize, c_int), &
            int(xsize, c_int) &
        ], &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize, tsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! The input array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "inpptr" only needs the first half by design.
    inpptr = cmplx(0.0, 0.0, c_double_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        c_double_complex &
    )
    call fftw_execute_dft_c2r(fftwplan, inpptr, outptr)
    ! Carefully specify the range of assignments because "outptr" contains extra padding.
    out(:, :, :, :) = real(outptr(1:xsize, :, :, :), real64)
    out = out / real(xsize * ysize * zsize * tsize, real64)
    nullify(inpptr)
    nullify(outptr)
    call fftw_destroy_plan(fftwplan)
    call fftw_free(fftwdata)
end subroutine dirfft4
subroutine difftshift4(out, inp)
    complex(real64), allocatable, intent(out) :: out(:, :, :, :)
    complex(real64), intent(in) :: inp(:, :, :, :)
    allocate(out, source=inp)
    out(:, :, :, :) = cshift(out, size(out, 1) / 2, 1)
    out(:, :, :, :) = cshift(out, size(out, 2) / 2, 2)
    out(:, :, :, :) = cshift(out, size(out, 3) / 2, 3)
    out(:, :, :, :) = cshift(out, size(out, 4) / 2, 4)
end subroutine difftshift4
subroutine scfft1(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:)
    complex(real32), intent(in) :: inp(:)
    integer, optional, intent(in) :: num(1)
    integer :: &
        xsize, &
        xinpsize
    complex(c_float_complex), pointer :: &
        inpptr(:) => null(), &
        outptr(:) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
    end if
    fftwdata = fftwf_alloc_complex(int(xsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize])
    call c_f_pointer(fftwdata, outptr, [xsize])
    fftwplan = fftwf_plan_dft_1d( &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_forward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize) = cmplx( &
        real( inp(1:xinpsize)), &
        aimag(inp(1:xinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft(fftwplan, inpptr, outptr)
    out(:) = cmplx(real(outptr), aimag(outptr), real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine scfft1
subroutine srfft1(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:)
    real(real32), intent(in) :: inp(:)
    integer, optional, intent(in) :: num(1)
    integer :: &
        xsize, &
        xhalfsize, &
        xinpsize
    real(c_float), pointer :: &
        inpptr(:) => null()
    complex(c_float_complex), pointer :: &
        outptr(:) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xsize, size(inp, 1))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xsize
    end if
    ! For in-place, real-input and complex-output ("r2c") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftwf_alloc_complex(int(xhalfsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xhalfsize * 2])
    call c_f_pointer(fftwdata, outptr, [xhalfsize])
    fftwplan = fftwf_plan_dft_r2c_1d( &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! Carefully specify the range of assignments because "inpptr" contains extra padding.
    inpptr = real(0.0, c_float)
    inpptr( 1:xinpsize) = real( &
        inp(1:xinpsize), &
        c_float &
    )
    call fftwf_execute_dft_r2c(fftwplan, inpptr, outptr)
    ! The output array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "outptr" only contains the first half by design.
    out(1:xhalfsize) = cmplx(real(outptr), aimag(outptr), real32)
    if (mod(xsize, 2) == 0) then
        ! Even number of elements.
        ! The element at "xhalfsize" is the Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1) = &
            conjg(out(2:(xhalfsize - 1)))
    else
        ! Odd number of elements.
        ! No Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1) = &
            conjg(out(2:(xhalfsize)))
    end if
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine srfft1
subroutine sfftshift1(out, inp)
    complex(real32), allocatable, intent(out) :: out(:)
    complex(real32), intent(in) :: inp(:)
    allocate(out, source=inp)
    out(:) = cshift(out, -size(out, 1) / 2, 1)
end subroutine sfftshift1
subroutine sicfft1(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:)
    complex(real32), intent(in) :: inp(:)
    integer, optional, intent(in) :: num(1)
    integer :: &
        xsize, &
        xinpsize
    complex(c_float_complex), pointer :: &
        inpptr(:) => null(), &
        outptr(:) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
    end if
    fftwdata = fftwf_alloc_complex(int(xsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize])
    call c_f_pointer(fftwdata, outptr, [xsize])
    fftwplan = fftwf_plan_dft_1d( &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_backward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize) = cmplx( &
        real( inp(1:xinpsize)), &
        aimag(inp(1:xinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft(fftwplan, inpptr, outptr)
    out(:) = cmplx(real(outptr), aimag(outptr), real32)
    out = out / real(xsize, real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine sicfft1
subroutine sirfft1(out, inp, num)
    real(real32), allocatable, intent(out) :: out(:)
    complex(real32), intent(in) :: inp(:)
    integer, optional, intent(in) :: num(1)
    integer :: &
        xsize, &
        xhalfsize, &
        xinpsize
    real(c_float), pointer :: &
        outptr(:) => null()
    complex(c_float_complex), pointer :: &
        inpptr(:) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xhalfsize, size(inp, 1) / 2 + 1)
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xhalfsize
    end if
    ! For in-place, complex-input and real-output ("c2r") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftwf_alloc_complex(int(xhalfsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize])
    fftwplan = fftwf_plan_dft_c2r_1d( &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! The input array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "inpptr" only needs the first half by design.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize) = cmplx( &
        real( inp(1:xinpsize)), &
        aimag(inp(1:xinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft_c2r(fftwplan, inpptr, outptr)
    ! Carefully specify the range of assignments because "outptr" contains extra padding.
    out(:) = real(outptr(1:xsize), real32)
    out = out / real(xsize, real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine sirfft1
subroutine sifftshift1(out, inp)
    complex(real32), allocatable, intent(out) :: out(:)
    complex(real32), intent(in) :: inp(:)
    allocate(out, source=inp)
    out(:) = cshift(out, size(out, 1) / 2, 1)
end subroutine sifftshift1
subroutine scfft2(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:, :)
    complex(real32), intent(in) :: inp(:, :)
    integer, optional, intent(in) :: num(2)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize
    complex(c_float_complex), pointer :: &
        inpptr(:, :) => null(), &
        outptr(:, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
    end if
    fftwdata = fftwf_alloc_complex(int(xsize * ysize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize])
    fftwplan = fftwf_plan_dft_2d( &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_forward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize, 1:yinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft(fftwplan, inpptr, outptr)
    out(:, :) = cmplx(real(outptr), aimag(outptr), real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine scfft2
subroutine srfft2(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:, :)
    real(real32), intent(in) :: inp(:, :)
    integer, optional, intent(in) :: num(2)
    integer :: &
        i, &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize
    real(c_float), pointer :: &
        inpptr(:, :) => null()
    complex(c_float_complex), pointer :: &
        outptr(:, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
    end if
    ! For in-place, real-input and complex-output ("r2c") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftwf_alloc_complex(int(xhalfsize * ysize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xhalfsize * 2, ysize])
    call c_f_pointer(fftwdata, outptr, [xhalfsize, ysize])
    fftwplan = fftwf_plan_dft_r2c_2d( &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! Carefully specify the range of assignments because "inpptr" contains extra padding.
    inpptr = real(0.0, c_float)
    inpptr( 1:xinpsize, 1:yinpsize) = real( &
        inp(1:xinpsize, 1:yinpsize), &
        c_float &
    )
    call fftwf_execute_dft_r2c(fftwplan, inpptr, outptr)
    ! The output array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "outptr" only contains the first half by design.
    out(1:xhalfsize, :) = cmplx(real(outptr), aimag(outptr), real32)
    if (mod(xsize, 2) == 0) then
        ! Even number of elements.
        ! The element at "xhalfsize" is the Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)]) = &
            conjg(out(2:(xhalfsize - 1), :))
    else
        ! Odd number of elements.
        ! No Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)]) = &
            conjg(out(2:(xhalfsize), :))
    end if
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine srfft2
subroutine sfftshift2(out, inp)
    complex(real32), allocatable, intent(out) :: out(:, :)
    complex(real32), intent(in) :: inp(:, :)
    allocate(out, source=inp)
    out(:, :) = cshift(out, -size(out, 1) / 2, 1)
    out(:, :) = cshift(out, -size(out, 2) / 2, 2)
end subroutine sfftshift2
subroutine sicfft2(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:, :)
    complex(real32), intent(in) :: inp(:, :)
    integer, optional, intent(in) :: num(2)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize
    complex(c_float_complex), pointer :: &
        inpptr(:, :) => null(), &
        outptr(:, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
    end if
    fftwdata = fftwf_alloc_complex(int(xsize * ysize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize])
    fftwplan = fftwf_plan_dft_2d( &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_backward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize, 1:yinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft(fftwplan, inpptr, outptr)
    out(:, :) = cmplx(real(outptr), aimag(outptr), real32)
    out = out / real(xsize * ysize, real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine sicfft2
subroutine sirfft2(out, inp, num)
    real(real32), allocatable, intent(out) :: out(:, :)
    complex(real32), intent(in) :: inp(:, :)
    integer, optional, intent(in) :: num(2)
    integer :: &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize
    real(c_float), pointer :: &
        outptr(:, :) => null()
    complex(c_float_complex), pointer :: &
        inpptr(:, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xhalfsize, size(inp, 1) / 2 + 1)
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xhalfsize
        ysize = size(inp, 2)
        yinpsize = ysize
    end if
    ! For in-place, complex-input and real-output ("c2r") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftwf_alloc_complex(int(xhalfsize * ysize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2, ysize])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize, ysize])
    fftwplan = fftwf_plan_dft_c2r_2d( &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! The input array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "inpptr" only needs the first half by design.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize, 1:yinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft_c2r(fftwplan, inpptr, outptr)
    ! Carefully specify the range of assignments because "outptr" contains extra padding.
    out(:, :) = real(outptr(1:xsize, :), real32)
    out = out / real(xsize * ysize, real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine sirfft2
subroutine sifftshift2(out, inp)
    complex(real32), allocatable, intent(out) :: out(:, :)
    complex(real32), intent(in) :: inp(:, :)
    allocate(out, source=inp)
    out(:, :) = cshift(out, size(out, 1) / 2, 1)
    out(:, :) = cshift(out, size(out, 2) / 2, 2)
end subroutine sifftshift2
subroutine scfft3(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:, :, :)
    complex(real32), intent(in) :: inp(:, :, :)
    integer, optional, intent(in) :: num(3)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize
    complex(c_float_complex), pointer :: &
        inpptr(:, :, :) => null(), &
        outptr(:, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
    end if
    fftwdata = fftwf_alloc_complex(int(xsize * ysize * zsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize, zsize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize, zsize])
    fftwplan = fftwf_plan_dft_3d( &
        int(zsize, c_int), &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_forward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft(fftwplan, inpptr, outptr)
    out(:, :, :) = cmplx(real(outptr), aimag(outptr), real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine scfft3
subroutine srfft3(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:, :, :)
    real(real32), intent(in) :: inp(:, :, :)
    integer, optional, intent(in) :: num(3)
    integer :: &
        i, &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize
    real(c_float), pointer :: &
        inpptr(:, :, :) => null()
    complex(c_float_complex), pointer :: &
        outptr(:, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
    end if
    ! For in-place, real-input and complex-output ("r2c") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftwf_alloc_complex(int(xhalfsize * ysize * zsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xhalfsize * 2, ysize, zsize])
    call c_f_pointer(fftwdata, outptr, [xhalfsize, ysize, zsize])
    fftwplan = fftwf_plan_dft_r2c_3d( &
        int(zsize, c_int), &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! Carefully specify the range of assignments because "inpptr" contains extra padding.
    inpptr = real(0.0, c_float)
    inpptr( 1:xinpsize, 1:yinpsize, 1:zinpsize) = real( &
        inp(1:xinpsize, 1:yinpsize, 1:zinpsize), &
        c_float &
    )
    call fftwf_execute_dft_r2c(fftwplan, inpptr, outptr)
    ! The output array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "outptr" only contains the first half by design.
    out(1:xhalfsize, :, :) = cmplx(real(outptr), aimag(outptr), real32)
    if (mod(xsize, 2) == 0) then
        ! Even number of elements.
        ! The element at "xhalfsize" is the Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)], &
            [1, (i, i = zsize, 2, -1)]) = &
            conjg(out(2:(xhalfsize - 1), :, :))
    else
        ! Odd number of elements.
        ! No Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)], &
            [1, (i, i = zsize, 2, -1)]) = &
            conjg(out(2:(xhalfsize), :, :))
    end if
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine srfft3
subroutine sfftshift3(out, inp)
    complex(real32), allocatable, intent(out) :: out(:, :, :)
    complex(real32), intent(in) :: inp(:, :, :)
    allocate(out, source=inp)
    out(:, :, :) = cshift(out, -size(out, 1) / 2, 1)
    out(:, :, :) = cshift(out, -size(out, 2) / 2, 2)
    out(:, :, :) = cshift(out, -size(out, 3) / 2, 3)
end subroutine sfftshift3
subroutine sicfft3(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:, :, :)
    complex(real32), intent(in) :: inp(:, :, :)
    integer, optional, intent(in) :: num(3)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize
    complex(c_float_complex), pointer :: &
        inpptr(:, :, :) => null(), &
        outptr(:, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
    end if
    fftwdata = fftwf_alloc_complex(int(xsize * ysize * zsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize, zsize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize, zsize])
    fftwplan = fftwf_plan_dft_3d( &
        int(zsize, c_int), &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_backward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft(fftwplan, inpptr, outptr)
    out(:, :, :) = cmplx(real(outptr), aimag(outptr), real32)
    out = out / real(xsize * ysize * zsize, real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine sicfft3
subroutine sirfft3(out, inp, num)
    real(real32), allocatable, intent(out) :: out(:, :, :)
    complex(real32), intent(in) :: inp(:, :, :)
    integer, optional, intent(in) :: num(3)
    integer :: &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize
    real(c_float), pointer :: &
        outptr(:, :, :) => null()
    complex(c_float_complex), pointer :: &
        inpptr(:, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xhalfsize, size(inp, 1) / 2 + 1)
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xhalfsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
    end if
    ! For in-place, complex-input and real-output ("c2r") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftwf_alloc_complex(int(xhalfsize * ysize * zsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2, ysize, zsize])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize, ysize, zsize])
    fftwplan = fftwf_plan_dft_c2r_3d( &
        int(zsize, c_int), &
        int(ysize, c_int), &
        int(xsize, c_int), &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! The input array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "inpptr" only needs the first half by design.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft_c2r(fftwplan, inpptr, outptr)
    ! Carefully specify the range of assignments because "outptr" contains extra padding.
    out(:, :, :) = real(outptr(1:xsize, :, :), real32)
    out = out / real(xsize * ysize * zsize, real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine sirfft3
subroutine sifftshift3(out, inp)
    complex(real32), allocatable, intent(out) :: out(:, :, :)
    complex(real32), intent(in) :: inp(:, :, :)
    allocate(out, source=inp)
    out(:, :, :) = cshift(out, size(out, 1) / 2, 1)
    out(:, :, :) = cshift(out, size(out, 2) / 2, 2)
    out(:, :, :) = cshift(out, size(out, 3) / 2, 3)
end subroutine sifftshift3
subroutine scfft4(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:, :, :, :)
    complex(real32), intent(in) :: inp(:, :, :, :)
    integer, optional, intent(in) :: num(4)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize, &
        tsize, &
        tinpsize
    complex(c_float_complex), pointer :: &
        inpptr(:, :, :, :) => null(), &
        outptr(:, :, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
        tsize = num(4)
        tinpsize = min(tsize, size(inp, 4))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
        tsize = size(inp, 4)
        tinpsize = tsize
    end if
    fftwdata = fftwf_alloc_complex(int(xsize * ysize * zsize * tsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize, zsize, tsize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize, zsize, tsize])
    fftwplan = fftwf_plan_dft( &
        4_c_int, [ &
            int(tsize, c_int), &
            int(zsize, c_int), &
            int(ysize, c_int), &
            int(xsize, c_int) &
        ], &
        inpptr, &
        outptr, &
        fftw_forward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize, tsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft(fftwplan, inpptr, outptr)
    out(:, :, :, :) = cmplx(real(outptr), aimag(outptr), real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine scfft4
subroutine srfft4(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:, :, :, :)
    real(real32), intent(in) :: inp(:, :, :, :)
    integer, optional, intent(in) :: num(4)
    integer :: &
        i, &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize, &
        tsize, &
        tinpsize
    real(c_float), pointer :: &
        inpptr(:, :, :, :) => null()
    complex(c_float_complex), pointer :: &
        outptr(:, :, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
        tsize = num(4)
        tinpsize = min(tsize, size(inp, 4))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
        tsize = size(inp, 4)
        tinpsize = tsize
    end if
    ! For in-place, real-input and complex-output ("r2c") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftwf_alloc_complex(int(xhalfsize * ysize * zsize * tsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xhalfsize * 2, ysize, zsize, tsize])
    call c_f_pointer(fftwdata, outptr, [xhalfsize, ysize, zsize, tsize])
    fftwplan = fftwf_plan_dft_r2c( &
        4_c_int, [ &
            int(tsize, c_int), &
            int(zsize, c_int), &
            int(ysize, c_int), &
            int(xsize, c_int) &
        ], &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize, tsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! Carefully specify the range of assignments because "inpptr" contains extra padding.
    inpptr = real(0.0, c_float)
    inpptr( 1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize) = real( &
        inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize), &
        c_float &
    )
    call fftwf_execute_dft_r2c(fftwplan, inpptr, outptr)
    ! The output array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "outptr" only contains the first half by design.
    out(1:xhalfsize, :, :, :) = cmplx(real(outptr), aimag(outptr), real32)
    if (mod(xsize, 2) == 0) then
        ! Even number of elements.
        ! The element at "xhalfsize" is the Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)], &
            [1, (i, i = zsize, 2, -1)], &
            [1, (i, i = tsize, 2, -1)]) = &
            conjg(out(2:(xhalfsize - 1), :, :, :))
    else
        ! Odd number of elements.
        ! No Nyquist frequency.
        out(xsize:(xhalfsize + 1):-1, &
            [1, (i, i = ysize, 2, -1)], &
            [1, (i, i = zsize, 2, -1)], &
            [1, (i, i = tsize, 2, -1)]) = &
            conjg(out(2:(xhalfsize), :, :, :))
    end if
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine srfft4
subroutine sfftshift4(out, inp)
    complex(real32), allocatable, intent(out) :: out(:, :, :, :)
    complex(real32), intent(in) :: inp(:, :, :, :)
    allocate(out, source=inp)
    out(:, :, :, :) = cshift(out, -size(out, 1) / 2, 1)
    out(:, :, :, :) = cshift(out, -size(out, 2) / 2, 2)
    out(:, :, :, :) = cshift(out, -size(out, 3) / 2, 3)
    out(:, :, :, :) = cshift(out, -size(out, 4) / 2, 4)
end subroutine sfftshift4
subroutine sicfft4(out, inp, num)
    complex(real32), allocatable, intent(out) :: out(:, :, :, :)
    complex(real32), intent(in) :: inp(:, :, :, :)
    integer, optional, intent(in) :: num(4)
    integer :: &
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize, &
        tsize, &
        tinpsize
    complex(c_float_complex), pointer :: &
        inpptr(:, :, :, :) => null(), &
        outptr(:, :, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
        tsize = num(4)
        tinpsize = min(tsize, size(inp, 4))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
        tsize = size(inp, 4)
        tinpsize = tsize
    end if
    fftwdata = fftwf_alloc_complex(int(xsize * ysize * zsize * tsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, inpptr, [xsize, ysize, zsize, tsize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize, zsize, tsize])
    fftwplan = fftwf_plan_dft( &
        4_c_int, [ &
            int(tsize, c_int), &
            int(zsize, c_int), &
            int(ysize, c_int), &
            int(xsize, c_int) &
        ], &
        inpptr, &
        outptr, &
        fftw_backward, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize, tsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft(fftwplan, inpptr, outptr)
    out(:, :, :, :) = cmplx(real(outptr), aimag(outptr), real32)
    out = out / real(xsize * ysize * zsize * tsize, real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine sicfft4
subroutine sirfft4(out, inp, num)
    real(real32), allocatable, intent(out) :: out(:, :, :, :)
    complex(real32), intent(in) :: inp(:, :, :, :)
    integer, optional, intent(in) :: num(4)
    integer :: &
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize, &
        tsize, &
        tinpsize
    real(c_float), pointer :: &
        outptr(:, :, :, :) => null()
    complex(c_float_complex), pointer :: &
        inpptr(:, :, :, :) => null()
    type(c_ptr) :: &
        fftwdata, &
        fftwplan
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xhalfsize, size(inp, 1) / 2 + 1)
        ysize = num(2)
        yinpsize = min(ysize, size(inp, 2))
        zsize = num(3)
        zinpsize = min(zsize, size(inp, 3))
        tsize = num(4)
        tinpsize = min(tsize, size(inp, 4))
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xhalfsize
        ysize = size(inp, 2)
        yinpsize = ysize
        zsize = size(inp, 3)
        zinpsize = zsize
        tsize = size(inp, 4)
        tinpsize = tsize
    end if
    ! For in-place, complex-input and real-output ("c2r") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.
    fftwdata = fftwf_alloc_complex(int(xhalfsize * ysize * zsize * tsize, c_size_t))
    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if
    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2, ysize, zsize, tsize])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize, ysize, zsize, tsize])
    fftwplan = fftwf_plan_dft_c2r( &
        4_c_int, [ &
            int(tsize, c_int), &
            int(zsize, c_int), &
            int(ysize, c_int), &
            int(xsize, c_int) &
        ], &
        inpptr, &
        outptr, &
        fftw_estimate &
    )
    if (.not. c_associated(fftwplan)) then
        print *, 'Error calling fftw_plan_dft_*'
        stop 1
    end if
    allocate(out(xsize, ysize, zsize, tsize))
    ! As stressed by the FFTW documents, initialize the input only after creating a plan.
    ! The input array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "inpptr" only needs the first half by design.
    inpptr = cmplx(0.0, 0.0, c_float_complex)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        c_float_complex &
    )
    call fftwf_execute_dft_c2r(fftwplan, inpptr, outptr)
    ! Carefully specify the range of assignments because "outptr" contains extra padding.
    out(:, :, :, :) = real(outptr(1:xsize, :, :, :), real32)
    out = out / real(xsize * ysize * zsize * tsize, real32)
    nullify(inpptr)
    nullify(outptr)
    call fftwf_destroy_plan(fftwplan)
    call fftwf_free(fftwdata)
end subroutine sirfft4
subroutine sifftshift4(out, inp)
    complex(real32), allocatable, intent(out) :: out(:, :, :, :)
    complex(real32), intent(in) :: inp(:, :, :, :)
    allocate(out, source=inp)
    out(:, :, :, :) = cshift(out, size(out, 1) / 2, 1)
    out(:, :, :, :) = cshift(out, size(out, 2) / 2, 2)
    out(:, :, :, :) = cshift(out, size(out, 3) / 2, 3)
    out(:, :, :, :) = cshift(out, size(out, 4) / 2, 4)
end subroutine sifftshift4
end module fftw3_wrapper
