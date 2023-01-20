#!/usr/bin/env python3

import numpy  as np
import random as rd

def get_random_array(iscomplex: bool, shape: tuple) -> np.ndarray:
	if (iscomplex):
		v = np.array(
				[complex(rd.uniform(-1.0, 1.0), rd.uniform(-1.0, 1.0)) for i in range(np.prod(shape))],
				dtype=np.cdouble
			).reshape(
				shape,
				order='F'
			)
	else:
		v = np.array(
				[rd.uniform(-1.0, 1.0) for i in range(np.prod(shape))],
				dtype=np.double
			).reshape(
				shape,
				order='F'
			)

	return v

def get_transform_size(shape: tuple) -> list:
	v = [np.array(shape, dtype=np.int_) for i in range(8)]

	v[0] = v[0] // 2 - 1
	v[1] = v[1] // 2
	v[2] = v[2] // 2 + 1
	v[3] = v[3] - 1
	v[4] = v[4] + 1
	v[5] = v[5] * 2 - 1
	v[6] = v[6] * 2
	v[7] = v[7] * 2 + 1

	return v

def write_test_array(filename: str, array: np.ndarray) -> None:
	np.savetxt(
		filename,
		array.flatten(order='F'),
		fmt='(%13.6E,%13.6E)' if np.issubdtype(array.dtype, np.complexfloating) else '%13.6E'
	)

def generate_test_case(iscomplex: bool, shape: tuple) -> None:
	x = get_random_array(iscomplex, shape)

	filename = (
		('cfft' if iscomplex else 'rfft') + '{:d}'.format(x.ndim) +
		'_x_' + ('{:02d}' + '-{:02d}' * (x.ndim - 1)).format(*(x.shape)) +
		'.txt'
	)

	write_test_array(filename, x)

	if (x.ndim == 1):
		y = np.fft.fft(x)
	elif (x.ndim == 2):
		y = np.fft.fft2(x)
	else:
		y = np.fft.fftn(x)

	filename = (
		('cfft' if iscomplex else 'rfft') + '{:d}'.format(x.ndim) +
		'_x_' + ('{:02d}' + '-{:02d}' * (x.ndim - 1)).format(*(x.shape)) +
		'_y_' + ('{:02d}' + '-{:02d}' * (y.ndim - 1)).format(*(y.shape)) +
		'.txt'
	)

	write_test_array(filename, y)

	l = get_transform_size(shape)

	for n in l:
		if (x.ndim == 1):
			y = np.fft.fft(x, n[0])
		elif (x.ndim == 2):
			y = np.fft.fft2(x, n)
		else:
			y = np.fft.fftn(x, n)

		filename = (
			('cfft' if iscomplex else 'rfft') + '{:d}'.format(x.ndim) +
			'_x_' + ('{:02d}' + '-{:02d}' * (x.ndim - 1)).format(*(x.shape)) +
			'_y_' + ('{:02d}' + '-{:02d}' * (y.ndim - 1)).format(*(y.shape)) +
			'.txt'
		)

		write_test_array(filename, y)

print('Generating random FFT test cases')

generate_test_case(True,  (10,))
generate_test_case(True,  (10, 13))
generate_test_case(True,  (10, 13, 16))
generate_test_case(True,  (10, 13, 16, 19))

generate_test_case(True,  (11,))
generate_test_case(True,  (11, 14))
generate_test_case(True,  (11, 14, 17))
generate_test_case(True,  (11, 14, 17, 20))

generate_test_case(False, (10,))
generate_test_case(False, (10, 13))
generate_test_case(False, (10, 13, 16))
generate_test_case(False, (10, 13, 16, 19))

generate_test_case(False, (11,))
generate_test_case(False, (11, 14))
generate_test_case(False, (11, 14, 17))
generate_test_case(False, (11, 14, 17, 20))
