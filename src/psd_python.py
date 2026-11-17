"""Precipitation Smoothing Distance for fields on a planar grid."""

from numbers import Real

import numpy as np
from scipy import signal
from skimage.draw import disk

__all__ = ["calculate_PSD", "calculate_PSS", "construct_circular_kernel_field"]


def _real_scalar(value, name, *, positive=False):
    """Validate scalar options without leaking NumPy type/overflow errors."""
    if not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        finite = np.isfinite(value)
    except TypeError:
        raise ValueError(f"{name} must be a finite real number") from None
    if not finite or (value <= 0 if positive else value < 0):
        qualifier = "positive" if positive else "non-negative"
        raise ValueError(f"{name} must be a finite, {qualifier} real number")


def _validate_fields(fa, fb):
    """Validate shared PSD/PSS array requirements; zero fields are permitted."""
    if isinstance(fa, np.ma.MaskedArray) or isinstance(fb, np.ma.MaskedArray):
        raise ValueError("masked arrays are not allowed")
    if not isinstance(fa, np.ndarray) or not isinstance(fb, np.ndarray):
        raise ValueError("fa and fb must be numpy arrays")
    if fa.shape != fb.shape:
        raise ValueError("fa and fb arrays do not have the same shape")
    if fa.ndim != 2 or fb.ndim != 2:
        raise ValueError("fa and fb arrays are not two-dimensional")
    if fa.size == 0:
        raise ValueError("the input arrays are empty")
    for field in (fa, fb):
        if field.dtype.kind != 'f':
            raise ValueError("fa and fb must be floating-point numpy arrays")
        if not np.isfinite(field).all():
            raise ValueError("fa or fb contain non-finite values")
        if (field < 0).any():
            raise ValueError("negative values are not allowed")


def construct_circular_kernel_field(r):
    """Return a normalized circular kernel, including pixels on its boundary."""
    _real_scalar(r, "r")
    r_int_max = int(np.floor(r))
    shape = (2 * r_int_max + 1, 2 * r_int_max + 1)
    kernel = np.zeros(shape)
    rr, cc = disk((r_int_max, r_int_max), r + 0.0001, shape=shape)
    kernel[rr, cc] = 1
    return kernel / kernel.sum()


def calculate_PSS(fa, fb, r, Q):
    """Calculate PSS for nonnegative floating-point 2-D residuals of mean Q.

    Invalid inputs and non-finite arithmetic raise ValueError. Callers are
    responsible for preparing the residual fields with the stated mean Q.
    """
    _real_scalar(Q, "Q", positive=True)
    _validate_fields(fa, fb)
    # Divide before subtracting to avoid overflow in the residual sum and
    # denominator when Q or the field values have extreme scales.
    dtype = np.result_type(fa.dtype, fb.dtype, np.asarray(Q).dtype, np.float64)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        a = np.asarray(fa, dtype=dtype) / Q
        b = np.asarray(fb, dtype=dtype) / Q
        difference = np.asarray(a - b, dtype=float)
    if not np.isfinite(difference).all():
        raise ValueError("residual fields divided by Q must remain finite")
    return _pss_from_difference(difference, r)


def _pss_from_difference(difference, r):
    """Smooth a validated difference whose residuals each have mean one."""
    kernel = construct_circular_kernel_field(r)
    # Linearity permits one convolution; full mode preserves mass at the edges.
    difference_smooth = signal.fftconvolve(difference, kernel, mode='full')
    score = float(1.0 - np.abs(difference_smooth).sum() / (2.0 * difference.size))
    if not np.isfinite(score):
        raise ValueError("smoothing produced a non-finite PSS")
    return score


def _error(message):
    print("ERROR: " + message + ". Returning -1 as result!")
    return -1


def _normalize(field):
    # Scale before summing to avoid overflow and underflow in the mean.
    # Promote float16/float32 before division, but retain wider floating types
    # until scaling so values outside float64's range remain supported.
    scaled = np.array(field, dtype=np.result_type(field.dtype, np.float64), copy=True)
    scaled /= scaled.max()
    scaled = np.asarray(scaled, dtype=float)
    scaled /= scaled.mean()
    return scaled


def calculate_PSD(fa, fb):
    """Estimate displacement for two nonnegative floating-point 2-D arrays.

    Invalid inputs and numerical search failures print a message and return -1.
    Differences at floating-point roundoff level are treated as zero.
    Input arrays are never modified.
    """
    try:
        _validate_fields(fa, fb)
    except ValueError as exc:
        return _error(str(exc))
    for field in (fa, fb):
        if field.max() == 0:
            print("WARNING: At least one field contains only zero values. "
                  "Returning -1 as result!")
            return -1

    a = _normalize(fa)
    b = _normalize(fb)
    difference = a - b
    positive = np.maximum(difference, 0)
    negative = np.maximum(-difference, 0)
    mass_a = float(positive.sum())
    mass_b = float(negative.sum())
    # Average both residual masses so roundoff cannot break symmetry.
    Q = (mass_a + mass_b) / (2.0 * a.size)
    if Q <= 16 * np.finfo(float).eps:
        return 0.0
    if mass_a == 0 or mass_b == 0 or not np.isfinite(Q):
        return _error("normalization produced invalid residual masses")
    # Balance residual masses and divide out Q once, before the radius search.
    # Each residual now has mean one; only their difference is needed by FFT.
    positive *= a.size / mass_a
    negative *= a.size / mass_b
    difference = positive - negative
    del a, b, positive, negative

    def evaluate(radius):
        value = _pss_from_difference(difference, radius)
        if not np.isfinite(value):
            raise ValueError("smoothing produced a non-finite PSS")
        return value

    # A finite search limit prevents runaway allocations if a crossing is absent.
    # Four domain diagonals is a conservative limit for the circular smoothing.
    diagonal = float(np.hypot(*fa.shape))
    step = max(1, int(np.ceil(diagonal * 0.05)))
    max_radius = max(2, int(np.ceil(4 * diagonal)))
    # FFT and normalization roundoff must not move an exact threshold crossing
    # to the next integer radius, which would materially change the distance.
    threshold = 0.5 - 16 * np.finfo(float).eps
    try:
        lower = 1
        if evaluate(lower) >= threshold:
            return float(0.808 * Q * lower)
        upper = min(lower + step, max_radius)
        while evaluate(upper) < threshold:
            if upper == max_radius:
                return _error("PSS did not reach 0.5 within the radius limit")
            lower = upper
            upper = min(upper + step, max_radius)

        while upper - lower > 1:
            middle = (lower + upper) // 2
            if evaluate(middle) >= threshold:
                upper = middle
            else:
                lower = middle
    except ValueError as exc:
        return _error(str(exc))
    return float(0.808 * Q * upper)
