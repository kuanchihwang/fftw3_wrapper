#define PP_CAT_I(a, b) a ## b
#define PP_CAT(a, b) PP_CAT_I(a, b)
#define PP_STR_I(a) #a
#define PP_STR(a) PP_STR_I(a)

#if defined(FFT_COMPLEX)
#define WRAP_DATA_TYPE_FAMILY c
#elif defined(FFT_REAL)
#define WRAP_DATA_TYPE_FAMILY r
#else
#error "Undefined FFT Data Type"
#endif

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

subroutine PP_CAT(PP_CAT(PP_CAT(WRAP_DATA_TYPE_FAMILY, fft), WRAP_DIMENSION_FAMILY), _test)(xshape)
    integer, intent(in) :: xshape(WRAP_DIMENSION_FAMILY)

    character(:), allocatable :: filename
    complex(real32), allocatable :: cbuffer(:)

#if defined(FFT_COMPLEX)
#if defined(FFT_1D)
    complex(real64), allocatable :: &
        xread(:), &
        xtemp(:), &
        xtest(:), &
        yread(:), &
        ytest(:)
    logical, allocatable :: &
        check(:)
#elif defined(FFT_2D)
    complex(real64), allocatable :: &
        xread(:, :), &
        xtemp(:, :), &
        xtest(:, :), &
        yread(:, :), &
        ytest(:, :)
    logical, allocatable :: &
        check(:, :)
#elif defined(FFT_3D)
    complex(real64), allocatable :: &
        xread(:, :, :), &
        xtemp(:, :, :), &
        xtest(:, :, :), &
        yread(:, :, :), &
        ytest(:, :, :)
    logical, allocatable :: &
        check(:, :, :)
#elif defined(FFT_4D)
    complex(real64), allocatable :: &
        xread(:, :, :, :), &
        xtemp(:, :, :, :), &
        xtest(:, :, :, :), &
        yread(:, :, :, :), &
        ytest(:, :, :, :)
    logical, allocatable :: &
        check(:, :, :, :)
#else
#error "Undefined FFT Dimension"
#endif
#elif defined(FFT_REAL)
    real(real32), allocatable :: rbuffer(:)
#if defined(FFT_1D)
    real(real64), allocatable :: &
        xread(:), &
        xtemp(:), &
        xtest(:)
    complex(real64), allocatable :: &
        yread(:), &
        ytest(:)
    logical, allocatable :: &
        check(:)
#elif defined(FFT_2D)
    real(real64), allocatable :: &
        xread(:, :), &
        xtemp(:, :), &
        xtest(:, :)
    complex(real64), allocatable :: &
        yread(:, :), &
        ytest(:, :)
    logical, allocatable :: &
        check(:, :)
#elif defined(FFT_3D)
    real(real64), allocatable :: &
        xread(:, :, :), &
        xtemp(:, :, :), &
        xtest(:, :, :)
    complex(real64), allocatable :: &
        yread(:, :, :), &
        ytest(:, :, :)
    logical, allocatable :: &
        check(:, :, :)
#elif defined(FFT_4D)
    real(real64), allocatable :: &
        xread(:, :, :, :), &
        xtemp(:, :, :, :), &
        xtest(:, :, :, :)
    complex(real64), allocatable :: &
        yread(:, :, :, :), &
        ytest(:, :, :, :)
    logical, allocatable :: &
        check(:, :, :, :)
#else
#error "Undefined FFT Dimension"
#endif
#else
#error "Undefined FFT Data Type"
#endif

    integer :: i
    integer, allocatable :: yshape(:, :)

#if defined(FFT_COMPLEX)
    call get_filename(filename, 'cfft', xshape)
    call read_test_array(cbuffer, filename, product(xshape))

    xread = reshape(cbuffer, xshape)

    call get_filename(filename, 'cfft', xshape, xshape)
    call read_test_array(cbuffer, filename, product(xshape))

    yread = reshape(cbuffer, xshape)
#elif defined(FFT_REAL)
    call get_filename(filename, 'rfft', xshape)
    call read_test_array(rbuffer, filename, product(xshape))

    xread = reshape(rbuffer, xshape)

    call get_filename(filename, 'rfft', xshape, xshape)
    call read_test_array(cbuffer, filename, product(xshape))

    yread = reshape(cbuffer, xshape)
#else
#error "Undefined FFT Data Type"
#endif

    print '(a)', 'Calling fft(y, x)'

    call fft(ytest, xread)

    check = almost_equal(ytest, yread)

    if (.not. all(check)) then
        print '(a)', 'Self-test failed'

#if defined(FFT_COMPLEX)
        print '(a, *("(", es13.6e2, ", ", es13.6e2, ")", :, " "))', &
            'Expected: ', pack(yread, .not. check)
        print '(a, *("(", es13.6e2, ", ", es13.6e2, ")", :, " "))', &
            'Actual:   ', pack(ytest, .not. check)
#elif defined(FFT_REAL)
        print '(a, *(es13.6e2, :, " "))', &
            'Expected: ', pack(yread, .not. check)
        print '(a, *(es13.6e2, :, " "))', &
            'Actual:   ', pack(ytest, .not. check)
#else
#error "Undefined FFT Data Type"
#endif

        stop 1
    end if

    print '(a)', 'Calling ifft(x, y)'

    call ifft(xtest, yread)

    check = almost_equal(xtest, xread)

    if (.not. all(check)) then
        print '(a)', 'Self-test failed'

#if defined(FFT_COMPLEX)
        print '(a, *("(", es13.6e2, ", ", es13.6e2, ")", :, " "))', &
            'Expected: ', pack(xread, .not. check)
        print '(a, *("(", es13.6e2, ", ", es13.6e2, ")", :, " "))', &
            'Actual:   ', pack(xtest, .not. check)
#elif defined(FFT_REAL)
        print '(a, *(es13.6e2, :, " "))', &
            'Expected: ', pack(xread, .not. check)
        print '(a, *(es13.6e2, :, " "))', &
            'Actual:   ', pack(xtest, .not. check)
#else
#error "Undefined FFT Data Type"
#endif

        stop 2
    end if

    call get_transform_size(yshape, xshape)

    do i = 1, size(yshape, 2)
#if defined(FFT_COMPLEX)
        call get_filename(filename, 'cfft', xshape, yshape(:, i))
#elif defined(FFT_REAL)
        call get_filename(filename, 'rfft', xshape, yshape(:, i))
#else
#error "Undefined FFT Data Type"
#endif
        call read_test_array(cbuffer, filename, product(yshape(:, i)))

        yread = reshape(cbuffer, yshape(1:WRAP_DIMENSION_FAMILY, i))

        print '(a)', 'Calling fft(y, x, n)'

        call fft(ytest, xread, yshape(:, i))

        check = almost_equal(ytest, yread)

        if (.not. all(check)) then
            print '(a)', 'Self-test failed'

#if defined(FFT_COMPLEX)
            print '(a, *("(", es13.6e2, ", ", es13.6e2, ")", :, " "))', &
                'Expected: ', pack(yread, .not. check)
            print '(a, *("(", es13.6e2, ", ", es13.6e2, ")", :, " "))', &
                'Actual:   ', pack(ytest, .not. check)
#elif defined(FFT_REAL)
            print '(a, *(es13.6e2, :, " "))', &
                'Expected: ', pack(yread, .not. check)
            print '(a, *(es13.6e2, :, " "))', &
                'Actual:   ', pack(ytest, .not. check)
#else
#error "Undefined FFT Data Type"
#endif

            stop 3
        end if

        print '(a)', 'Calling ifft(x, y)'

        call ifft(xtest, yread)

        if (allocated(xtemp)) deallocate(xtemp)

        allocate(xtemp, mold=xtest)

        xtemp = 0.0_real64
#if defined(FFT_1D)
        xtemp( &
            1:min(size(xread, 1), size(xtemp, 1))  &
        ) = xread( &
            1:min(size(xread, 1), size(xtemp, 1))  &
        )
#elif defined(FFT_2D)
        xtemp( &
            1:min(size(xread, 1), size(xtemp, 1)), &
            1:min(size(xread, 2), size(xtemp, 2))  &
        ) = xread( &
            1:min(size(xread, 1), size(xtemp, 1)), &
            1:min(size(xread, 2), size(xtemp, 2))  &
        )
#elif defined(FFT_3D)
        xtemp( &
            1:min(size(xread, 1), size(xtemp, 1)), &
            1:min(size(xread, 2), size(xtemp, 2)), &
            1:min(size(xread, 3), size(xtemp, 3))  &
        ) = xread( &
            1:min(size(xread, 1), size(xtemp, 1)), &
            1:min(size(xread, 2), size(xtemp, 2)), &
            1:min(size(xread, 3), size(xtemp, 3))  &
        )
#elif defined(FFT_4D)
        xtemp( &
            1:min(size(xread, 1), size(xtemp, 1)), &
            1:min(size(xread, 2), size(xtemp, 2)), &
            1:min(size(xread, 3), size(xtemp, 3)), &
            1:min(size(xread, 4), size(xtemp, 4))  &
        ) = xread( &
            1:min(size(xread, 1), size(xtemp, 1)), &
            1:min(size(xread, 2), size(xtemp, 2)), &
            1:min(size(xread, 3), size(xtemp, 3)), &
            1:min(size(xread, 4), size(xtemp, 4))  &
        )
#else
#error "Undefined FFT Dimension"
#endif

        check = almost_equal(xtest, xtemp)

        if (.not. all(check)) then
            print '(a)', 'Self-test failed'

#if defined(FFT_COMPLEX)
            print '(a, *("(", es13.6e2, ", ", es13.6e2, ")", :, " "))', &
                'Expected: ', pack(xtemp, .not. check)
            print '(a, *("(", es13.6e2, ", ", es13.6e2, ")", :, " "))', &
                'Actual:   ', pack(xtest, .not. check)
#elif defined(FFT_REAL)
            print '(a, *(es13.6e2, :, " "))', &
                'Expected: ', pack(xtemp, .not. check)
            print '(a, *(es13.6e2, :, " "))', &
                'Actual:   ', pack(xtest, .not. check)
#else
#error "Undefined FFT Data Type"
#endif

            stop 4
        end if
    end do
end subroutine PP_CAT(PP_CAT(PP_CAT(WRAP_DATA_TYPE_FAMILY, fft), WRAP_DIMENSION_FAMILY), _test)

#undef WRAP_DATA_TYPE_FAMILY

#undef WRAP_DIMENSION_FAMILY
