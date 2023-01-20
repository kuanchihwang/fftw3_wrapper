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
#if defined(HAVE_LIBFFTW3F)
        module procedure &
            scfft1, srfft1, &
            scfft2, srfft2, &
            scfft3, srfft3, &
            scfft4, srfft4
#endif
    end interface fft

    ! Forward complex fast Fourier transform family

    interface cfft
        module procedure &
            dcfft1, &
            dcfft2, &
            dcfft3, &
            dcfft4
#if defined(HAVE_LIBFFTW3F)
        module procedure &
            scfft1, &
            scfft2, &
            scfft3, &
            scfft4
#endif
    end interface cfft

    ! Forward real fast Fourier transform family

    interface rfft
        module procedure &
            drfft1, &
            drfft2, &
            drfft3, &
            drfft4
#if defined(HAVE_LIBFFTW3F)
        module procedure &
            srfft1, &
            srfft2, &
            srfft3, &
            srfft4
#endif
    end interface rfft

    ! Inverse fast Fourier transform family

    interface ifft
        module procedure &
            dicfft1, dirfft1, &
            dicfft2, dirfft2, &
            dicfft3, dirfft3, &
            dicfft4, dirfft4
#if defined(HAVE_LIBFFTW3F)
        module procedure &
            sicfft1, sirfft1, &
            sicfft2, sirfft2, &
            sicfft3, sirfft3, &
            sicfft4, sirfft4
#endif
    end interface ifft

    ! Inverse complex fast Fourier transform family

    interface icfft
        module procedure &
            dicfft1, &
            dicfft2, &
            dicfft3, &
            dicfft4
#if defined(HAVE_LIBFFTW3F)
        module procedure &
            sicfft1, &
            sicfft2, &
            sicfft3, &
            sicfft4
#endif
    end interface icfft

    ! Inverse real fast Fourier transform family

    interface irfft
        module procedure &
            dirfft1, &
            dirfft2, &
            dirfft3, &
            dirfft4
#if defined(HAVE_LIBFFTW3F)
        module procedure &
            sirfft1, &
            sirfft2, &
            sirfft3, &
            sirfft4
#endif
    end interface irfft

    ! Forward zero-frequency shift family

    interface fftshift
        module procedure &
            dfftshift1, &
            dfftshift2, &
            dfftshift3, &
            dfftshift4
#if defined(HAVE_LIBFFTW3F)
        module procedure &
            sfftshift1, &
            sfftshift2, &
            sfftshift3, &
            sfftshift4
#endif
    end interface fftshift

    ! Inverse zero-frequency shift family

    interface ifftshift
        module procedure &
            difftshift1, &
            difftshift2, &
            difftshift3, &
            difftshift4
#if defined(HAVE_LIBFFTW3F)
        module procedure &
            sifftshift1, &
            sifftshift2, &
            sifftshift3, &
            sifftshift4
#endif
    end interface ifftshift

    ! Generic public procedures

    public :: &
        fft,  cfft,  rfft,  fftshift, &
        ifft, icfft, irfft, ifftshift

    ! Specific public procedures

    public :: &
        dcfft1, drfft1, dicfft1, dirfft1, &
        dcfft2, drfft2, dicfft2, dirfft2, &
        dcfft3, drfft3, dicfft3, dirfft3, &
        dcfft4, drfft4, dicfft4, dirfft4

#if defined(HAVE_LIBFFTW3F)
    public :: &
        scfft1, srfft1, sicfft1, sirfft1, &
        scfft2, srfft2, sicfft2, sirfft2, &
        scfft3, srfft3, sicfft3, sirfft3, &
        scfft4, srfft4, sicfft4, sirfft4
#endif
contains
#define FFT_DOUBLE
#define FFT_1D
#include "fftw3_wrapper_forward_implementation.F90"
#include "fftw3_wrapper_inverse_implementation.F90"
#undef  FFT_1D

#define FFT_2D
#include "fftw3_wrapper_forward_implementation.F90"
#include "fftw3_wrapper_inverse_implementation.F90"
#undef  FFT_2D

#define FFT_3D
#include "fftw3_wrapper_forward_implementation.F90"
#include "fftw3_wrapper_inverse_implementation.F90"
#undef  FFT_3D

#define FFT_4D
#include "fftw3_wrapper_forward_implementation.F90"
#include "fftw3_wrapper_inverse_implementation.F90"
#undef  FFT_4D
#undef  FFT_DOUBLE

#if defined(HAVE_LIBFFTW3F)
#define FFT_SINGLE
#define FFT_1D
#include "fftw3_wrapper_forward_implementation.F90"
#include "fftw3_wrapper_inverse_implementation.F90"
#undef  FFT_1D

#define FFT_2D
#include "fftw3_wrapper_forward_implementation.F90"
#include "fftw3_wrapper_inverse_implementation.F90"
#undef  FFT_2D

#define FFT_3D
#include "fftw3_wrapper_forward_implementation.F90"
#include "fftw3_wrapper_inverse_implementation.F90"
#undef  FFT_3D

#define FFT_4D
#include "fftw3_wrapper_forward_implementation.F90"
#include "fftw3_wrapper_inverse_implementation.F90"
#undef  FFT_4D
#undef  FFT_SINGLE
#endif
end module fftw3_wrapper
