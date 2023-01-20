program fftw3_wrapper_test
    use, intrinsic :: iso_fortran_env
    use :: fftw3_wrapper

    implicit none

    interface almost_equal
        procedure almost_equal_complex, almost_equal_real
    end interface almost_equal

    interface read_test_array
        procedure read_test_array_complex, read_test_array_real
    end interface read_test_array

    real(real64), parameter :: errortolerance = 1000.0_real64 * epsilon(1.0_real32)

    print '(a)', 'Testing 1-D FFT'
    call cfft1_test([10])
    call cfft1_test([11])
    call rfft1_test([10])
    call rfft1_test([11])
    print '(a)', 'PASS!'

    print '(a)', 'Testing 2-D FFT'
    call cfft2_test([10, 13])
    call cfft2_test([11, 14])
    call rfft2_test([10, 13])
    call rfft2_test([11, 14])
    print '(a)', 'PASS!'

    print '(a)', 'Testing 3-D FFT'
    call cfft3_test([10, 13, 16])
    call cfft3_test([11, 14, 17])
    call rfft3_test([10, 13, 16])
    call rfft3_test([11, 14, 17])
    print '(a)', 'PASS!'

    print '(a)', 'Testing 4-D FFT'
    call cfft4_test([10, 13, 16, 19])
    call cfft4_test([11, 14, 17, 20])
    call rfft4_test([10, 13, 16, 19])
    call rfft4_test([11, 14, 17, 20])
    print '(a)', 'PASS!'
contains
! AlmostEqualRelativeAndAbs
! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition

pure elemental function almost_equal_complex(a, b) result(r)
    complex(real64), intent(in) :: a, b
    logical :: r

    real(real64) :: d, m

    d = abs(a - b)
    m = max(abs(a), abs(b))

    r = d < errortolerance

    if (r) return

    r = d < errortolerance * m
end function almost_equal_complex

pure elemental function almost_equal_real(a, b) result(r)
    real(real64), intent(in) :: a, b
    logical :: r

    real(real64) :: d, m

    d = abs(a - b)
    m = max(abs(a), abs(b))

    r = d < errortolerance

    if (r) return

    r = d < errortolerance * m
end function almost_equal_real

pure subroutine get_filename(filename, filenameprefix, xshape, yshape)
    character(:), allocatable, intent(out) :: filename
    character(*), intent(in) :: filenameprefix
    integer, intent(in) :: xshape(:)
    integer, optional, intent(in) :: yshape(:)

    character(4) :: filenamefragment
    integer :: i

    if (size(xshape) < 1) return

    if (present(yshape)) then
        if (size(yshape) < 1) return
    end if

    write(filenamefragment, '(i4.1)') size(xshape)
    filename = trim(adjustl(filenameprefix)) // trim(adjustl(filenamefragment))

    filename = filename // '_x_'
    write(filenamefragment, '(i4.2)') xshape(1)
    filename = filename // trim(adjustl(filenamefragment))

    do i = 2, size(xshape)
        filename = filename // '-'
        write(filenamefragment, '(i4.2)') xshape(i)
        filename = filename // trim(adjustl(filenamefragment))
    end do

    if (present(yshape)) then
        filename = filename // '_y_'
        write(filenamefragment, '(i4.2)') yshape(1)
        filename = filename // trim(adjustl(filenamefragment))

        do i = 2, size(yshape)
            filename = filename // '-'
            write(filenamefragment, '(i4.2)') yshape(i)
            filename = filename // trim(adjustl(filenamefragment))
        end do
    end if

    filename = filename // '.txt'
end subroutine get_filename

pure subroutine get_transform_size(yshape, xshape)
    integer, allocatable, intent(out) :: yshape(:, :)
    integer, intent(in) :: xshape(:)

    integer :: i

    allocate(yshape(size(xshape), 9))

    do i = 1, size(yshape, 2)
        yshape(:, i) = xshape
    end do

    yshape(:, 1) = yshape(:, 1) / 2 - 1
    yshape(:, 2) = yshape(:, 2) / 2
    yshape(:, 3) = yshape(:, 3) / 2 + 1
    yshape(:, 4) = yshape(:, 4) - 1
    ! yshape(:, 5) = yshape(:, 5)
    yshape(:, 6) = yshape(:, 6) + 1
    yshape(:, 7) = yshape(:, 7) * 2 - 1
    yshape(:, 8) = yshape(:, 8) * 2
    yshape(:, 9) = yshape(:, 9) * 2 + 1
end subroutine get_transform_size

subroutine read_test_array_complex(buffer, filename, length)
    complex(real32), allocatable, intent(out) :: buffer(:)
    character(*), intent(in) :: filename
    integer, intent(in) :: length

    integer :: i, u

    open(newunit=u, file=filename, action='read', form='formatted', status='old')

    allocate(buffer(length))

    do i = 1, length
        read(u, *) buffer(i)
    end do

    close(u)
end subroutine read_test_array_complex

subroutine read_test_array_real(buffer, filename, length)
    real(real32), allocatable, intent(out) :: buffer(:)
    character(*), intent(in) :: filename
    integer, intent(in) :: length

    integer :: i, u

    open(newunit=u, file=filename, action='read', form='formatted', status='old')

    allocate(buffer(length))

    do i = 1, length
        read(u, *) buffer(i)
    end do

    close(u)
end subroutine read_test_array_real

#define FFT_COMPLEX
#define FFT_1D
#include "fftw3_wrapper_test_implementation.F90"
#undef  FFT_1D

#define FFT_2D
#include "fftw3_wrapper_test_implementation.F90"
#undef  FFT_2D

#define FFT_3D
#include "fftw3_wrapper_test_implementation.F90"
#undef  FFT_3D

#define FFT_4D
#include "fftw3_wrapper_test_implementation.F90"
#undef  FFT_4D
#undef  FFT_COMPLEX

#define FFT_REAL
#define FFT_1D
#include "fftw3_wrapper_test_implementation.F90"
#undef  FFT_1D

#define FFT_2D
#include "fftw3_wrapper_test_implementation.F90"
#undef  FFT_2D

#define FFT_3D
#include "fftw3_wrapper_test_implementation.F90"
#undef  FFT_3D

#define FFT_4D
#include "fftw3_wrapper_test_implementation.F90"
#undef  FFT_4D
#undef  FFT_REAL
end program fftw3_wrapper_test
