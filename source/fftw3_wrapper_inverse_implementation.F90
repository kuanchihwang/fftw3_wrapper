#define PP_CAT_I(a, b) a ## b
#define PP_CAT(a, b) PP_CAT_I(a, b)
#define PP_STR_I(a) #a
#define PP_STR(a) PP_STR_I(a)

#if defined(FFT_1D)
#define WRAP_DIMENSION_FAMILY 1
#elif defined(FFT_2D)
#define WRAP_DIMENSION_FAMILY 2
#elif defined(FFT_3D)
#define WRAP_DIMENSION_FAMILY 3
#elif defined(FFT_4D)
#define WRAP_DIMENSION_FAMILY 4
#else
#error "Undefined FFT Dimension"
#endif

#if defined(FFT_DOUBLE)
#define FFTW_PRECISION_FAMILY fftw_
#define WRAP_PRECISION_FAMILY d
#define F_FLOATING_POINT_TYPE real64
#define C_FLOATING_POINT_TYPE c_double
#define F_COMPLEX_TYPE        real64
#define C_COMPLEX_TYPE        c_double_complex
#elif defined(FFT_SINGLE)
#define FFTW_PRECISION_FAMILY fftwf_
#define WRAP_PRECISION_FAMILY s
#define F_FLOATING_POINT_TYPE real32
#define C_FLOATING_POINT_TYPE c_float
#define F_COMPLEX_TYPE        real32
#define C_COMPLEX_TYPE        c_float_complex
#else
#error "Undefined FFT Precision"
#endif

subroutine PP_CAT(PP_CAT(WRAP_PRECISION_FAMILY, icfft), WRAP_DIMENSION_FAMILY)(out, inp, num)
#if defined(FFT_1D)
    complex(F_COMPLEX_TYPE), allocatable, intent(out) :: out(:)
    complex(F_COMPLEX_TYPE),              intent(in)  :: inp(:)
#elif defined(FFT_2D)
    complex(F_COMPLEX_TYPE), allocatable, intent(out) :: out(:, :)
    complex(F_COMPLEX_TYPE),              intent(in)  :: inp(:, :)
#elif defined(FFT_3D)
    complex(F_COMPLEX_TYPE), allocatable, intent(out) :: out(:, :, :)
    complex(F_COMPLEX_TYPE),              intent(in)  :: inp(:, :, :)
#elif defined(FFT_4D)
    complex(F_COMPLEX_TYPE), allocatable, intent(out) :: out(:, :, :, :)
    complex(F_COMPLEX_TYPE),              intent(in)  :: inp(:, :, :, :)
#else
#error "Undefined FFT Dimension"
#endif
    integer, optional,                    intent(in)  :: num(WRAP_DIMENSION_FAMILY)

    integer :: &
#if defined(FFT_1D)
        xsize, &
        xinpsize
#elif defined(FFT_2D)
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize
#elif defined(FFT_3D)
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize
#elif defined(FFT_4D)
        xsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize, &
        tsize, &
        tinpsize
#else
#error "Undefined FFT Dimension"
#endif

    complex(C_COMPLEX_TYPE), pointer :: &
#if defined(FFT_1D)
        inpptr(:)          => null(), &
        outptr(:)          => null()
#elif defined(FFT_2D)
        inpptr(:, :)       => null(), &
        outptr(:, :)       => null()
#elif defined(FFT_3D)
        inpptr(:, :, :)    => null(), &
        outptr(:, :, :)    => null()
#elif defined(FFT_4D)
        inpptr(:, :, :, :) => null(), &
        outptr(:, :, :, :) => null()
#else
#error "Undefined FFT Dimension"
#endif
    type(c_ptr) :: &
        fftwdata, &
        fftwplan

#if defined(FFT_1D)
    if (present(num)) then
        xsize = num(1)
        xinpsize = min(xsize, size(inp, 1))
    else
        xsize = size(inp, 1)
        xinpsize = xsize
    end if
#elif defined(FFT_2D)
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
#elif defined(FFT_3D)
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
#elif defined(FFT_4D)
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
#else
#error "Undefined FFT Dimension"
#endif

#if defined(FFT_1D)
    fftwdata = PP_CAT(FFTW_PRECISION_FAMILY, alloc_complex)(int(xsize, c_size_t))

    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if

    call c_f_pointer(fftwdata, inpptr, [xsize])
    call c_f_pointer(fftwdata, outptr, [xsize])

    fftwplan = PP_CAT(FFTW_PRECISION_FAMILY, plan_dft_1d)( &
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
#elif defined(FFT_2D)
    fftwdata = PP_CAT(FFTW_PRECISION_FAMILY, alloc_complex)(int(xsize * ysize, c_size_t))

    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if

    call c_f_pointer(fftwdata, inpptr, [xsize, ysize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize])

    fftwplan = PP_CAT(FFTW_PRECISION_FAMILY, plan_dft_2d)( &
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
#elif defined(FFT_3D)
    fftwdata = PP_CAT(FFTW_PRECISION_FAMILY, alloc_complex)(int(xsize * ysize * zsize, c_size_t))

    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if

    call c_f_pointer(fftwdata, inpptr, [xsize, ysize, zsize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize, zsize])

    fftwplan = PP_CAT(FFTW_PRECISION_FAMILY, plan_dft_3d)( &
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
#elif defined(FFT_4D)
    fftwdata = PP_CAT(FFTW_PRECISION_FAMILY, alloc_complex)(int(xsize * ysize * zsize * tsize, c_size_t))

    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if

    call c_f_pointer(fftwdata, inpptr, [xsize, ysize, zsize, tsize])
    call c_f_pointer(fftwdata, outptr, [xsize, ysize, zsize, tsize])

    fftwplan = PP_CAT(FFTW_PRECISION_FAMILY, plan_dft)( &
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
#else
#error "Undefined FFT Dimension"
#endif

    ! As stressed by the FFTW documents, initialize the input only after creating a plan.

    inpptr = cmplx(0.0, 0.0, C_COMPLEX_TYPE)
#if defined(FFT_1D)
    inpptr(1:xinpsize) = cmplx( &
        real( inp(1:xinpsize)), &
        aimag(inp(1:xinpsize)), &
        C_COMPLEX_TYPE &
    )
#elif defined(FFT_2D)
    inpptr(1:xinpsize, 1:yinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize)), &
        C_COMPLEX_TYPE &
    )
#elif defined(FFT_3D)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        C_COMPLEX_TYPE &
    )
#elif defined(FFT_4D)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        C_COMPLEX_TYPE &
    )
#else
#error "Undefined FFT Dimension"
#endif

    call PP_CAT(FFTW_PRECISION_FAMILY, execute_dft)(fftwplan, inpptr, outptr)

#if defined(FFT_1D)
    out(:)          = cmplx(real(outptr), aimag(outptr), F_COMPLEX_TYPE)
    out = out / real(xsize,                         F_FLOATING_POINT_TYPE)
#elif defined(FFT_2D)
    out(:, :)       = cmplx(real(outptr), aimag(outptr), F_COMPLEX_TYPE)
    out = out / real(xsize * ysize,                 F_FLOATING_POINT_TYPE)
#elif defined(FFT_3D)
    out(:, :, :)    = cmplx(real(outptr), aimag(outptr), F_COMPLEX_TYPE)
    out = out / real(xsize * ysize * zsize,         F_FLOATING_POINT_TYPE)
#elif defined(FFT_4D)
    out(:, :, :, :) = cmplx(real(outptr), aimag(outptr), F_COMPLEX_TYPE)
    out = out / real(xsize * ysize * zsize * tsize, F_FLOATING_POINT_TYPE)
#else
#error "Undefined FFT Dimension"
#endif

    nullify(inpptr)
    nullify(outptr)

    call PP_CAT(FFTW_PRECISION_FAMILY, destroy_plan)(fftwplan)
    call PP_CAT(FFTW_PRECISION_FAMILY, free)(fftwdata)
end subroutine PP_CAT(PP_CAT(WRAP_PRECISION_FAMILY, icfft), WRAP_DIMENSION_FAMILY)

subroutine PP_CAT(PP_CAT(WRAP_PRECISION_FAMILY, irfft), WRAP_DIMENSION_FAMILY)(out, inp, num)
#if defined(FFT_1D)
    real(F_FLOATING_POINT_TYPE), allocatable, intent(out) :: out(:)
    complex(F_COMPLEX_TYPE),                  intent(in)  :: inp(:)
#elif defined(FFT_2D)
    real(F_FLOATING_POINT_TYPE), allocatable, intent(out) :: out(:, :)
    complex(F_COMPLEX_TYPE),                  intent(in)  :: inp(:, :)
#elif defined(FFT_3D)
    real(F_FLOATING_POINT_TYPE), allocatable, intent(out) :: out(:, :, :)
    complex(F_COMPLEX_TYPE),                  intent(in)  :: inp(:, :, :)
#elif defined(FFT_4D)
    real(F_FLOATING_POINT_TYPE), allocatable, intent(out) :: out(:, :, :, :)
    complex(F_COMPLEX_TYPE),                  intent(in)  :: inp(:, :, :, :)
#else
#error "Undefined FFT Dimension"
#endif
    integer, optional,                    intent(in)  :: num(WRAP_DIMENSION_FAMILY)

    integer :: &
#if defined(FFT_1D)
        xsize, &
        xhalfsize, &
        xinpsize
#elif defined(FFT_2D)
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize
#elif defined(FFT_3D)
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize
#elif defined(FFT_4D)
        xsize, &
        xhalfsize, &
        xinpsize, &
        ysize, &
        yinpsize, &
        zsize, &
        zinpsize, &
        tsize, &
        tinpsize
#else
#error "Undefined FFT Dimension"
#endif

    real(C_FLOATING_POINT_TYPE), pointer :: &
#if defined(FFT_1D)
        outptr(:)          => null()
#elif defined(FFT_2D)
        outptr(:, :)       => null()
#elif defined(FFT_3D)
        outptr(:, :, :)    => null()
#elif defined(FFT_4D)
        outptr(:, :, :, :) => null()
#else
#error "Undefined FFT Dimension"
#endif
    complex(C_COMPLEX_TYPE), pointer :: &
#if defined(FFT_1D)
        inpptr(:)          => null()
#elif defined(FFT_2D)
        inpptr(:, :)       => null()
#elif defined(FFT_3D)
        inpptr(:, :, :)    => null()
#elif defined(FFT_4D)
        inpptr(:, :, :, :) => null()
#else
#error "Undefined FFT Dimension"
#endif
    type(c_ptr) :: &
        fftwdata, &
        fftwplan

#if defined(FFT_1D)
    if (present(num)) then
        xsize = num(1)
        xhalfsize = xsize / 2 + 1
        xinpsize = min(xhalfsize, size(inp, 1) / 2 + 1)
    else
        xsize = size(inp, 1)
        xhalfsize = xsize / 2 + 1
        xinpsize = xhalfsize
    end if
#elif defined(FFT_2D)
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
#elif defined(FFT_3D)
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
#elif defined(FFT_4D)
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
#else
#error "Undefined FFT Dimension"
#endif

    ! For in-place, complex-input and real-output ("c2r") transforms,
    ! the complex array is roughly halved in the last* dimension and
    ! the real array requires extra padding in the last* dimension.
    !
    ! *: "Last" for C but "first" for Fortran.

#if defined(FFT_1D)
    fftwdata = PP_CAT(FFTW_PRECISION_FAMILY, alloc_complex)(int(xhalfsize, c_size_t))

    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if

    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize])

    fftwplan = PP_CAT(FFTW_PRECISION_FAMILY, plan_dft_c2r_1d)( &
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
#elif defined(FFT_2D)
    fftwdata = PP_CAT(FFTW_PRECISION_FAMILY, alloc_complex)(int(xhalfsize * ysize, c_size_t))

    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if

    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2, ysize])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize,     ysize])

    fftwplan = PP_CAT(FFTW_PRECISION_FAMILY, plan_dft_c2r_2d)( &
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
#elif defined(FFT_3D)
    fftwdata = PP_CAT(FFTW_PRECISION_FAMILY, alloc_complex)(int(xhalfsize * ysize * zsize, c_size_t))

    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if

    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2, ysize, zsize])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize,     ysize, zsize])

    fftwplan = PP_CAT(FFTW_PRECISION_FAMILY, plan_dft_c2r_3d)( &
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
#elif defined(FFT_4D)
    fftwdata = PP_CAT(FFTW_PRECISION_FAMILY, alloc_complex)(int(xhalfsize * ysize * zsize * tsize, c_size_t))

    if (.not. c_associated(fftwdata)) then
        print *, 'Error calling fftw_alloc_*'
        stop 1
    end if

    call c_f_pointer(fftwdata, outptr, [xhalfsize * 2, ysize, zsize, tsize])
    call c_f_pointer(fftwdata, inpptr, [xhalfsize,     ysize, zsize, tsize])

    fftwplan = PP_CAT(FFTW_PRECISION_FAMILY, plan_dft_c2r)( &
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
#else
#error "Undefined FFT Dimension"
#endif

    ! As stressed by the FFTW documents, initialize the input only after creating a plan.

    ! The input array has Hermitian symmetry.
    ! The other half of the elements is the complex conjugate of the first half.
    ! Here "inpptr" only needs the first half by design.

    inpptr = cmplx(0.0, 0.0, C_COMPLEX_TYPE)
#if defined(FFT_1D)
    inpptr(1:xinpsize) = cmplx( &
        real( inp(1:xinpsize)), &
        aimag(inp(1:xinpsize)), &
        C_COMPLEX_TYPE &
    )
#elif defined(FFT_2D)
    inpptr(1:xinpsize, 1:yinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize)), &
        C_COMPLEX_TYPE &
    )
#elif defined(FFT_3D)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize)), &
        C_COMPLEX_TYPE &
    )
#elif defined(FFT_4D)
    inpptr(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize) = cmplx( &
        real( inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        aimag(inp(1:xinpsize, 1:yinpsize, 1:zinpsize, 1:tinpsize)), &
        C_COMPLEX_TYPE &
    )
#else
#error "Undefined FFT Dimension"
#endif

    call PP_CAT(FFTW_PRECISION_FAMILY, execute_dft_c2r)(fftwplan, inpptr, outptr)

    ! Carefully specify the range of assignments because "outptr" contains extra padding.

#if defined(FFT_1D)
    out(:)          = real(outptr(1:xsize),          F_FLOATING_POINT_TYPE)
    out = out / real(xsize,                          F_FLOATING_POINT_TYPE)
#elif defined(FFT_2D)
    out(:, :)       = real(outptr(1:xsize, :),       F_FLOATING_POINT_TYPE)
    out = out / real(xsize * ysize,                  F_FLOATING_POINT_TYPE)
#elif defined(FFT_3D)
    out(:, :, :)    = real(outptr(1:xsize, :, :),    F_FLOATING_POINT_TYPE)
    out = out / real(xsize * ysize * zsize,          F_FLOATING_POINT_TYPE)
#elif defined(FFT_4D)
    out(:, :, :, :) = real(outptr(1:xsize, :, :, :), F_FLOATING_POINT_TYPE)
    out = out / real(xsize * ysize * zsize * tsize,  F_FLOATING_POINT_TYPE)
#else
#error "Undefined FFT Dimension"
#endif

    nullify(inpptr)
    nullify(outptr)

    call PP_CAT(FFTW_PRECISION_FAMILY, destroy_plan)(fftwplan)
    call PP_CAT(FFTW_PRECISION_FAMILY, free)(fftwdata)
end subroutine PP_CAT(PP_CAT(WRAP_PRECISION_FAMILY, irfft), WRAP_DIMENSION_FAMILY)

subroutine PP_CAT(PP_CAT(WRAP_PRECISION_FAMILY, ifftshift), WRAP_DIMENSION_FAMILY)(out, inp)
#if defined(FFT_1D)
    complex(F_COMPLEX_TYPE), allocatable, intent(out) :: out(:)
    complex(F_COMPLEX_TYPE),              intent(in)  :: inp(:)
#elif defined(FFT_2D)
    complex(F_COMPLEX_TYPE), allocatable, intent(out) :: out(:, :)
    complex(F_COMPLEX_TYPE),              intent(in)  :: inp(:, :)
#elif defined(FFT_3D)
    complex(F_COMPLEX_TYPE), allocatable, intent(out) :: out(:, :, :)
    complex(F_COMPLEX_TYPE),              intent(in)  :: inp(:, :, :)
#elif defined(FFT_4D)
    complex(F_COMPLEX_TYPE), allocatable, intent(out) :: out(:, :, :, :)
    complex(F_COMPLEX_TYPE),              intent(in)  :: inp(:, :, :, :)
#else
#error "Undefined FFT Dimension"
#endif

    allocate(out, source=inp)

#if defined(FFT_1D)
    out(:)          = cshift(out, size(out, 1) / 2, 1)
#elif defined(FFT_2D)
    out(:, :)       = cshift(out, size(out, 1) / 2, 1)
    out(:, :)       = cshift(out, size(out, 2) / 2, 2)
#elif defined(FFT_3D)
    out(:, :, :)    = cshift(out, size(out, 1) / 2, 1)
    out(:, :, :)    = cshift(out, size(out, 2) / 2, 2)
    out(:, :, :)    = cshift(out, size(out, 3) / 2, 3)
#elif defined(FFT_4D)
    out(:, :, :, :) = cshift(out, size(out, 1) / 2, 1)
    out(:, :, :, :) = cshift(out, size(out, 2) / 2, 2)
    out(:, :, :, :) = cshift(out, size(out, 3) / 2, 3)
    out(:, :, :, :) = cshift(out, size(out, 4) / 2, 4)
#else
#error "Undefined FFT Dimension"
#endif
end subroutine PP_CAT(PP_CAT(WRAP_PRECISION_FAMILY, ifftshift), WRAP_DIMENSION_FAMILY)

#undef WRAP_DIMENSION_FAMILY

#undef FFTW_PRECISION_FAMILY
#undef WRAP_PRECISION_FAMILY
#undef F_FLOATING_POINT_TYPE
#undef C_FLOATING_POINT_TYPE
#undef F_COMPLEX_TYPE
#undef C_COMPLEX_TYPE
