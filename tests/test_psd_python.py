"""Regression tests for an installed psd_python; no source-path injection."""
import contextlib
from collections import defaultdict
from fractions import Fraction
import io
import unittest
from unittest.mock import patch

import numpy as np
from scipy import signal

import psd_python as psd


def exact_small_field_reference(a, b):
    """PSD oracle using rational probabilities and explicit pixel spreading.

    This deliberately uses neither library helpers, FFTs, mean normalization
    nor bisection. Small fixtures have a one-cell search step in production.
    """
    left = {index: Fraction(float(a[index])) for index in np.ndindex(a.shape) if a[index]}
    right = {index: Fraction(float(b[index])) for index in np.ndindex(b.shape) if b[index]}
    total_a, total_b = sum(left.values()), sum(right.values())
    difference = {index: left.get(index, 0) / total_a - right.get(index, 0) / total_b
                  for index in left.keys() | right.keys()}
    difference = {index: value for index, value in difference.items() if value}
    q = sum(abs(value) for value in difference.values()) / 2
    if not q:
        return 0., 0, 0.
    for radius in range(1, 20):
        offsets = [(y, x) for y in range(-radius, radius + 1)
                   for x in range(-radius, radius + 1) if x*x + y*y <= radius*radius]
        spread = defaultdict(Fraction)
        for (y, x), value in difference.items():
            for dy, dx in offsets:
                spread[y + dy, x + dx] += value
        score = 1 - sum(abs(value) for value in spread.values()) / (2 * q * len(offsets))
        if score >= Fraction(1, 2):
            return float(Fraction(808, 1000) * q * radius), radius, float(q)
    raise AssertionError('Small reference fixture did not cross the threshold')


class PSDRegression(unittest.TestCase):
    def call(self, a, b):
        with contextlib.redirect_stdout(io.StringIO()):
            return psd.calculate_PSD(a, b)

    def fields(self):
        a = np.zeros((12, 12)); b = a.copy()
        a[2:4, 2:4] = 1; b[6:8, 7:9] = 1
        return a, b

    def test_documented_example(self):
        a = np.zeros((500, 500)); b = a.copy()
        a[200:220, 200:220] = 1; b[200:220, 300:320] = 1
        self.assertAlmostEqual(self.call(a, b), 100.192, places=10)

    def test_symmetry(self):
        rng = np.random.default_rng(17)
        for _ in range(20):
            a, b = rng.random((9, 8)), rng.random((9, 8))
            self.assertAlmostEqual(self.call(a, b), self.call(b, a), places=12)

    def test_scale_invariance(self):
        a, b = self.fields()
        baseline = self.call(a, b)
        for sa, sb in [(1e308, 1e-300), (5e-324, 5e-324), (3., 17.)]:
            self.assertAlmostEqual(self.call(a * sa, b * sb), baseline, places=12)

    def test_overflow_regression(self):
        a = np.array([[1., 1.], [0., 0.]])
        b = a[::-1].copy()
        self.assertGreater(self.call(a * 1e308, b * 1e308), 0)
        self.assertEqual(self.call(a * 1e308, b * 1e308), self.call(a, b))

    def test_subnormal_regression(self):
        a, b = np.eye(2), np.fliplr(np.eye(2))
        self.assertEqual(self.call(a * 5e-324, b * 5e-324), self.call(a, b))

    def test_near_identical_both_orders(self):
        a = np.ones((2, 2)); b = a.copy()
        b[0, 0] = np.nextafter(1., 2.)
        self.assertEqual(self.call(a, b), 0)
        self.assertEqual(self.call(b, a), 0)

    def test_small_radius_formula(self):
        stripes = np.tile(np.arange(30) % 2, (30, 1))
        for eps in [1., .01, 1e-6]:
            a, b = 1 + eps * stripes, 1 + eps * (1 - stripes)
            self.assertAlmostEqual(self.call(a, b), .808 * eps / (2 + eps), places=14)

    def test_identical_and_scaled(self):
        a, _ = self.fields()
        self.assertEqual(self.call(a, a), 0)
        self.assertEqual(self.call(a, 7 * a), 0)

    def test_invalid_inputs(self):
        good = np.ones((2, 2))
        bad = [np.ones((3, 2)), np.ones(4), np.empty((0, 2)),
               np.zeros((2, 2)), -good, good * np.nan, good * np.inf,
               good.astype(complex), good.astype(str), good.astype(object),
               np.ma.array(good), [[1, 1], [1, 1]], None,
               np.full((2, 2), np.datetime64('2020-01-01'))]
        for value in bad:
            with self.subTest(dtype=getattr(value, 'dtype', type(value))):
                self.assertEqual(self.call(good, value), -1)
                self.assertEqual(self.call(value, good), -1)

    def test_input_unchanged(self):
        a, b = self.fields()
        aa, bb = a.copy(), b.copy()
        a.flags.writeable = False; b.flags.writeable = False
        self.call(a, b)
        np.testing.assert_array_equal(a, aa)
        np.testing.assert_array_equal(b, bb)

    def test_matching_invalid_dimensions_and_empty_fields(self):
        # Matching shapes reach dimension/size validation, rather than failing
        # the earlier shape-mismatch check.
        for shape in [(), (3,), (1, 2, 3), (0, 2), (2, 0), (0, 0)]:
            with self.subTest(shape=shape):
                a, b = np.ones(shape, dtype=float), np.ones(shape, dtype=float)
                self.assertEqual(self.call(a, b), -1)
                with self.assertRaises(ValueError):
                    psd.calculate_PSS(a, b, 1, 1)

    def test_integer_and_boolean_rejected(self):
        a, b = self.fields()
        for dtype in [np.int8, np.uint8, np.int64, np.uint64, bool]:
            with self.subTest(dtype=dtype):
                x, y = a.astype(dtype), b.astype(dtype)
                self.assertEqual(self.call(x, b), -1)
                self.assertEqual(self.call(a, y), -1)
                self.assertEqual(self.call(x, y), -1)
                with self.assertRaises(ValueError):
                    psd.calculate_PSS(x, b, 1, 1)
                with self.assertRaises(ValueError):
                    psd.calculate_PSS(a, y, 1, 1)
                self.assertEqual(self.call(x.astype(float), y.astype(float)), self.call(a, b))

    def test_full_convolution_equivalence(self):
        rng = np.random.default_rng(71)
        a, b = rng.random((6, 7)), rng.random((6, 7))
        for r in [1, 2, 4]:
            kernel = psd.construct_circular_kernel_field(r)
            old = 1 - np.abs(signal.fftconvolve(a, kernel, 'full') -
                             signal.fftconvolve(b, kernel, 'full')).sum() / (2 * a.size * .3)
            self.assertAlmostEqual(psd.calculate_PSS(a, b, r, .3), old, places=14)

    def test_kernel(self):
        for r in [0, 1, 2, 2.5]:
            kernel = psd.construct_circular_kernel_field(r)
            self.assertAlmostEqual(kernel.sum(), 1)
            np.testing.assert_array_equal(kernel, kernel[::-1, ::-1])
        for r in [-1, np.nan, np.inf, 1j, 1+0j, 'bad', None, [1]]:
            with self.assertRaises(ValueError):
                psd.construct_circular_kernel_field(r)

    def test_invalid_q(self):
        a, b = self.fields()
        for q in [0, -1, np.nan, np.inf, 1j, 1+0j, 'bad', None, [1]]:
            with self.assertRaises(ValueError):
                psd.calculate_PSS(a, b, 1, q)

    def test_search_limit(self):
        a, b = self.fields()
        with patch.object(psd, '_pss_from_difference', return_value=0) as mock:
            self.assertEqual(self.call(a, b), -1)
            self.assertLessEqual(mock.call_count, 83)

    def test_nonfinite_pss(self):
        a, b = self.fields()
        for value in [np.nan, np.inf, -np.inf]:
            with patch.object(psd, '_pss_from_difference', return_value=value):
                self.assertEqual(self.call(a, b), -1)

    def test_helper_rejects_nonfinite_convolution(self):
        a, b = self.fields()
        for value in [np.nan, np.inf]:
            with patch.object(psd.signal, 'fftconvolve', return_value=np.array([[value]])):
                with self.assertRaises(ValueError):
                    psd.calculate_PSS(a, b, 1, 1)

    def test_exact_threshold(self):
        a, b = self.fields()
        with patch.object(psd, '_pss_from_difference', return_value=.5):
            self.assertAlmostEqual(self.call(a, b), .808)

    def test_low_precision_inputs_promoted_before_normalization(self):
        for dtype in [np.float16, np.float32]:
            a = np.array([[1, 3], [3, 1]], dtype=dtype)
            self.assertEqual(self.call(a, a.astype(np.float64)), 0)
            self.assertAlmostEqual(self.call(a, np.fliplr(a)), .808, places=14)

    def test_extended_precision_range(self):
        if np.finfo(np.longdouble).max <= np.finfo(float).max:
            self.skipTest('longdouble has no wider range on this platform')
        a = np.array([[1, 3], [3, 1]], dtype=np.longdouble)
        for scale in [np.longdouble('1e400'), np.longdouble('1e-400')]:
            self.assertAlmostEqual(self.call(a * scale, np.fliplr(a) * scale), .808, places=14)

    def test_pss_float_types_and_extreme_scales(self):
        for dtype, value in [(np.float16, 2), (np.float32, 2), (np.float64, 2),
                             (float, 1e308), (float, np.nextafter(0., 1.) * 2)]:
            a = np.array([[value, 0]], dtype=dtype)
            b = np.array([[0, value]], dtype=dtype)
            self.assertAlmostEqual(psd.calculate_PSS(a, b, 0, value / 2), 0, places=14)

    def test_pss_rejects_invalid_arrays(self):
        a = np.ones((2, 2))
        for invalid in [np.ones((1, 2)), np.ones(2), np.empty((0, 2)),
                        np.full((2, 2), np.nan), -a, a.astype(complex),
                        a.astype(object), np.ma.array(a), None, [[1, 1], [1, 1]]]:
            with self.subTest(value=type(invalid)):
                with self.assertRaises(ValueError):
                    psd.calculate_PSS(a, invalid, 1, 1)
        with self.assertRaises(ValueError):
            psd.calculate_PSS(a, np.zeros_like(a), 1, np.nextafter(0., 1.))

    def test_pss_independent_spatial_reference(self):
        rng = np.random.default_rng(72)
        for shape in [(1, 5), (5, 1), (4, 7)]:
            a, b = rng.random(shape), rng.random(shape)
            for r in [0, 1, 2, 3]:
                coordinates = np.arange(-r, r + 1)
                yy, xx = np.meshgrid(coordinates, coordinates, indexing='ij')
                kernel = (xx * xx + yy * yy <= r * r).astype(float)
                kernel /= kernel.sum()
                expected = 1 - np.abs(signal.convolve2d(a, kernel, 'full') -
                                      signal.convolve2d(b, kernel, 'full')).sum() / (2 * a.size * .3)
                self.assertAlmostEqual(psd.calculate_PSS(a, b, r, .3), expected, places=13)

    def test_real_exact_threshold_roundoff(self):
        # Equal-mass binary fields have exactly PSS(1)=0.5 by direct summation.
        a = np.array([[0, 0, 1, 0, 1], [1, 1, 0, 0, 0],
                      [1, 0, 1, 0, 0], [1, 0, 1, 0, 1]], dtype=float)
        b = np.array([[0, 1, 0, 1, 1], [0, 0, 0, 1, 0],
                      [0, 1, 0, 0, 0], [0, 1, 1, 1, 1]], dtype=float)
        for transform in [lambda x: x, np.transpose, np.fliplr, np.flipud]:
            self.assertAlmostEqual(self.call(transform(a), transform(b)), .808 * 2 / 3, places=14)

    def test_threshold_tolerance_is_small(self):
        a, b = self.fields()
        with patch.object(psd, '_pss_from_difference', side_effect=lambda difference, r: .5 - 1e-10 if r == 1 else .6):
            self.assertAlmostEqual(self.call(a, b), .808 * 2)

    def test_large_invalid_scalar_options(self):
        a, b = self.fields()
        with self.assertRaises(ValueError):
            psd.construct_circular_kernel_field(10 ** 400)
        with self.assertRaises(ValueError):
            psd.calculate_PSS(a, b, 1, 10 ** 400)

    def test_exact_end_to_end_reference_fields(self):
        cases = []
        # Isolated events, including an event at the boundary.
        for start, end in [((3, 3), (3, 4)), ((0, 0), (0, 5)), ((1, 1), (5, 6))]:
            a = np.zeros((8, 10)); b = a.copy()
            a[start] = 2.; b[end] = 7.
            cases.append(('isolated-' + str(end), a, b))
        # Partially overlapping rain bands with different intensities.
        a = np.zeros((8, 10)); b = a.copy()
        a[2:5, 1:4] = 2.; b[2:5, 3:6] = 5.
        cases.append(('partial-overlap', a, b))
        # Two regions with different intensities and displacements.
        a = np.zeros((8, 10)); b = a.copy()
        a[1:3, 1:3] = 4.; a[5:7, 6:8] = 1.
        b[1:3, 3:5] = 3.; b[4:6, 5:7] = 2.
        cases.append(('multiple-regions', a, b))
        # A shared background with localized positive and negative residuals.
        a = np.ones((6, 8)); b = a.copy()
        a[1:3, 0:2] = 5.; b[3:5, 5:7] = 3.
        cases.append(('background', a, b))
        cases.append(('proportional', a, 3 * a))
        for name, a, b in cases:
            with self.subTest(field=name):
                expected, radius, q = exact_small_field_reference(a, b)
                self.assertAlmostEqual(self.call(a, b), expected, places=12)
                self.assertAlmostEqual(self.call(b, a), expected, places=12)

    def test_seeded_shapes_and_geometric_invariance(self):
        rng = np.random.default_rng(401)
        for shape in [(1, 1), (1, 7), (7, 1), (4, 6), (7, 9)]:
            for _ in range(2):
                a, b = rng.random(shape), rng.random(shape)
                expected = self.call(a, b)
                self.assertGreaterEqual(expected, 0)
                self.assertTrue(np.isfinite(expected))
                for name, x, y in [('transpose', a.T, b.T),
                                   ('reflection', a[::-1, ::-1], b[::-1, ::-1]),
                                   ('scaled', a * 1e200, b * 1e-200),
                                   ('Fortran', np.asfortranarray(a), np.asfortranarray(b))]:
                    with self.subTest(shape=shape, transform=name):
                        self.assertAlmostEqual(self.call(x, y), expected, places=12)

    def test_strided_readonly_inputs(self):
        a, b = self.fields()
        x, y = np.zeros((24, 24)), np.zeros((24, 24))
        x[::2, ::2], y[::2, ::2] = a, b
        x, y = x[::2, ::2], y[::2, ::2]
        x.flags.writeable = False; y.flags.writeable = False
        self.assertFalse(x.flags.c_contiguous)
        self.assertEqual(self.call(x, y), self.call(a, b))
        np.testing.assert_array_equal(x, a)
        np.testing.assert_array_equal(y, b)

    def test_fast_bisection_multiple_crossing_limitation(self):
        # Preserve the documented behavior of the deliberately fast search.
        a = np.zeros((100, 100)); b = a.copy()
        stripes = np.tile(np.arange(30) % 2, (30, 1)).astype(float)
        a[5:35, 5:35] = stripes / stripes.sum() * .58
        b[5:35, 5:35] = (1 - stripes) / stripes.sum() * .58
        a[80, 60] = .42; b[80, 90] = .42
        scores = [psd.calculate_PSS(a / a.mean(), b / b.mean(), r, 1.) for r in [2, 3, 4]]
        self.assertGreater(scores[0], .5)
        self.assertLess(scores[1], .5)
        self.assertGreater(scores[2], .5)
        self.assertAlmostEqual(self.call(a, b), 3.232, places=12)
        self.assertAlmostEqual(self.call(np.pad(a, 50), np.pad(b, 50)), 1.616, places=12)


if __name__ == '__main__':
    unittest.main(verbosity=2)
