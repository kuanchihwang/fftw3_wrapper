# fftw3_wrapper

- [fftw3_wrapper](#fftw3_wrapper)
  - [Author](#author)
  - [Generic Interfaces](#generic-interfaces)
  - [Specific Interfaces](#specific-interfaces)
  - [Compilation](#compilation)
  - [Testing](#testing)
  - [Usage](#usage)
  - [To Do](#to-do)

`fftw3_wrapper` is a Fortran wrapper library around [FFTW](https://www.fftw.org), whose official introduction states that:

> FFTW is a C subroutine library for computing the discrete Fourier transform (DFT) in one or more dimensions, of arbitrary input size, and of both real and complex data (as well as of even/odd data, i.e. the discrete cosine/sine transforms or DCT/DST).

`fftw3_wrapper` provides easy-to-use Fortran interfaces that take the burdens of cross-language interoperability, memory management, and implementation details of FFTW away from ordinary users. `fftw3_wrapper` tries to mimic the syntax of [`fft()`](https://www.mathworks.com/help/matlab/ref/fft.html), [`ifft()`](https://www.mathworks.com/help/matlab/ref/ifft.html)... from Matlab and [`numpy.fft.fft()`](https://numpy.org/doc/stable/reference/generated/numpy.fft.fft.html), [`numpy.fft.ifft()`](https://numpy.org/doc/stable/reference/generated/numpy.fft.ifft.html)... from Python NumPy. Users familiar with those programming languages should already know how to use this library.

## Author

Kuan-Chih "Stargazer" Wang

## Generic Interfaces

* `subroutine fft(out, inp, num)`

    Compute the one-, two-, three-, or four-dimensional forward discrete Fourier transform (DFT) of input array `inp` into output array `out` by using a fast Fourier transform (FFT) algorithm.

    * `out`: Array; type and kind of **complex(real32)** or **complex(real64)**; rank of **1**, **2**, **3**, or **4**; **allocatable**

        Output array containing the Fourier coefficients. It has the same kind and rank as the input array `inp`. Its shape is determined by the transform length `num`.

    * `inp`: Array; type and kind of **complex(real32)**, **complex(real64)**, **real(real32)**, or **real(real64)**; rank of **1**, **2**, **3**, or **4**

        Input array containing the original data.

    * `num`: Array; type and kind of **integer**; rank of **1**; **optional**

        Transform length along each dimension of the input array `inp`. If `num` is absent, it is equivalent to specifying `shape(inp)`. If a particular dimension length in `num` is smaller than that in `inp`, the `inp` is truncated. If a particular dimension length in `num` is larger than that in `inp`, the `inp` is padded with trailing zeros.

* `subroutine cfft(out, inp, num)`

    Same as `subroutine fft(out, inp, num)`. However, only complex-valued input array `inp` is accepted.

* `subroutine rfft(out, inp, num)`

    Same as `subroutine fft(out, inp, num)`. However, only real-valued input array `inp` is accepted.

* `subroutine fftshift(out, inp)`

    Shift the zero-frequency component to the center of spectrum.

    * `out`: Array; type and kind of **complex(real32)** or **complex(real64)**; rank of **1**, **2**, **3**, or **4**; **allocatable**

        Output array.

    * `inp`: Array; type and kind of **complex(real32)** or **complex(real64)**; rank of **1**, **2**, **3**, or **4**

        Input array.

* `subroutine ifft(out, inp, num)`

    Compute the one-, two-, three-, or four-dimensional inverse discrete Fourier transform (DFT) of input array `inp` into output array `out` by using a fast Fourier transform (FFT) algorithm.

    * `out`: Array; type and kind of **complex(real32)**, **complex(real64)**, **real(real32)**, or **real(real64)**; rank of **1**, **2**, **3**, or **4**; **allocatable**

        Output array containing the original data. It has the same kind and rank as the input array `inp`. Its shape is determined by the transform length `num`.

    * `inp`: Array; type and kind of **complex(real32)** or **complex(real64)**; rank of **1**, **2**, **3**, or **4**

        Input array containing the Fourier coefficients.

    * `num`: Array; type and kind of **integer**; rank of **1**; **optional**

        Transform length along each dimension of the input array `inp`. If `num` is absent, it is equivalent to specifying `shape(inp)`. If a particular dimension length in `num` is smaller than that in `inp`, the `inp` is truncated. If a particular dimension length in `num` is larger than that in `inp`, the `inp` is padded with trailing zeros.

* `subroutine icfft(out, inp, num)`

    Same as `subroutine ifft(out, inp, num)`. However, only complex-valued output array `out` is accepted.

* `subroutine irfft(out, inp, num)`

    Same as `subroutine ifft(out, inp, num)`. However, only real-valued output array `out` is accepted.

* `subroutine ifftshift(out, inp)`

    Inverse the actions carried out by `subroutine fftshift(out, inp)`.

    * `out`: Array; type and kind of **complex(real32)** or **complex(real64)**; rank of **1**, **2**, **3**, or **4**; **allocatable**

        Output array.

    * `inp`: Array; type and kind of **complex(real32)** or **complex(real64)**; rank of **1**, **2**, **3**, or **4**

        Input array.

## Specific Interfaces

Ordinary users should just call generic interfaces instead of these. Anyway, specific interfaces are named according to the following patterns:

`{d,s}{,i}{c,r}fft{1,2,3,4}`

where

* `d` indicates that the arguments are double precision. `s` indicates that the arguments are single precision.
* If `i` is absent, then it is forward FFT. If `i` is present, then it is inverse FFT.
* `c` indicates complex-input forward FFT or complex-output inverse FFT, depending on the presence of `i`. `r` indicates real-input forward FFT or real-output inverse FFT, depending on the presence of `i`.
* `1`, `2`, `3`, and `4` indicate the dimensionality of FFT.

For example, `drfft2` is a subroutine for computing double precision, real-input, and two-dimensional forward FFT. Its inverse would be `dirfft2`.

## Compilation

* Dependencies: FFTW (`libfftw3`: Double precision version; `libfftw3f`: Single precision version)
* Requirements: `gcc`, `gfortran`, `make`

1. Compile and install FFTW.

    ```bash
    # Single precision version of FFTW
    ./configure --prefix=/opt/fftw \
    --enable-shared --enable-static \
    --enable-openmp --enable-threads \
    --enable-avx2 --enable-single
    make
    make install

    # Double precision version of FFTW
    ./configure --prefix=/opt/fftw \
    --enable-shared --enable-static \
    --enable-openmp --enable-threads \
    --enable-avx2
    make
    make install
    ```

2. Type `make`. By default, `Makefile` expects that FFTW is installed at `/opt/fftw`.

    * Supply `FFTWPATH=<PATH>` to override the path to FFTW.
    * Supply `DEBUG=1` to compile with debugging options.

    For example, suppose that you have installed FFTW under your own home directory at `$HOME/fftw`. You would invoke `make` like:

    ```bash
    make FFTWPATH="$HOME/fftw"
    ```

3. Both shared (`libfftw3_wrapper.so`) and static (`libfftw3_wrapper.a`) libraries are generated along with Fortran module files (`fftw3.mod` and `fftw3_wrapper.mod`).

## Testing

* Additional Requirements: Python, NumPy

Type `make check` to run the self-tests, which check the implementation correctness of `fftw3_wrapper` against Python NumPy under different transform types, dimensions, and lengths. All results must show "PASS".

## Usage

To use `fftw3_wrapper` in a Fortran program unit, include the following `use` statement near the beginning:

```fortran
use :: fftw3_wrapper
```

When compiling, remember to specify the include paths (e.g. `-I`...) and link against (e.g. `-L`, `-l`...) the `libfftw3_wrapper`, `libfftw3`, and `libfftw3f` libraries.

## To Do

* Support CMake
* Support FFTW wisdom
* Support multi-threaded FFTW
