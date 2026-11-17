# PSD Python: Precipitation Smoothing Distance

Precipitation Smoothing Distance (PSD) estimates spatial displacement between
precipitation fields on a planar grid. This Python library uses circular
smoothing with FFT convolution and a bisection search to return a distance
in grid cells.

The package is available on PyPI and is distributed and imported as `psd_python`.

GitHub repository: [skokg/PSD](https://github.com/skokg/PSD).

## Installation

Python 3.10 or newer is required. Install the package with:

```sh
pip install psd_python
```

Alternatively, download and extract the source from the
[GitHub repository](https://github.com/skokg/PSD). Copy `psd_python.py` from the
`src` folder into the folder containing your Python script, then import it as
shown below. NumPy, SciPy and scikit-image must already be installed in your
Python environment.

## Usage

Import the library with:

```python
import psd_python as psd
```

Inputs are two equally shaped, nonempty 2D floating-point NumPy arrays. They represent the
same rectangular grid, with equal spacing along both axes. Values must be
real, finite and nonnegative, with at least one positive value in each field.
Masked arrays are not supported. Input arrays are never modified.

### Example: displaced precipitation fields

Two 500 × 500 fields contain identical 20 × 20 rainy squares, with the forecast
square shifted 100 grid cells to the right. After installing the library, run
the following complete example:

```python
import numpy as np
import psd_python as psd

observed = np.zeros((500, 500), dtype=float)
forecast = observed.copy()
observed[200:220, 200:220] = 1
forecast[200:220, 300:320] = 1

distance = psd.calculate_PSD(observed, forecast)
print(f"PSD: {distance:.4f} grid cells")
# PSD: 100.1920 grid cells
```

For grid spacing of 2 km, multiply this result by 2 to obtain a distance in km.
The calculation uses grid indices; no geographic coordinates are supplied.

## Function descriptions

### `psd.calculate_PSD(fa, fb)`

**Inputs:** `fa` and `fb` are floating-point NumPy arrays shaped `(y, x)`, with
identical shapes and nonempty dimensions. Each must contain only finite,
nonnegative values and at least one positive value. Integer, Boolean, masked,
complex, string and object arrays are rejected. Convert integer data explicitly
with `field.astype(float)` before calling the function.

**Output:** A scalar PSD distance in grid cells. Invalid inputs or numerical
search failures print a message and return `-1`. An all-zero field returns
`-1`, including when both fields are all zero. Identical normalized fields
return zero.

**Functionality:** Normalizes both fields, removes their overlap and computes
the unmatched fraction `Q`. Circular FFT smoothing and a fast bisection search
estimate the radius `r` where PSS reaches 0.5, returning `0.808 * Q * r`.
Radii start at one grid cell and are resolved to one grid cell; smoothing
preserves contributions beyond the input edges.

Normalization removes differences in total precipitation. Stable scaling and
roundoff tolerances handle extreme values and nearly identical fields.
Search is limited to four domain diagonals; numerical failures return `-1`.

Discrete smoothing can produce multiple threshold crossings, so bisection
may miss the earliest one. Adding dry borders can change the selected crossing;
compare PSD on a fixed domain. Large radii require more memory.

## Author and license

Gregor Skok, Faculty of Mathematics and Physics, University of Ljubljana,
Slovenia.

Email: [Gregor.Skok@fmf.uni-lj.si](mailto:Gregor.Skok@fmf.uni-lj.si).

Licensed under the MIT License. See LICENSE.

## References

Please cite the PSD paper when publishing work using this measure:

- Skok, G. (2022). A New Spatial Distance Metric for Verification of
  Precipitation. Applied Sciences, 12(8), 4048.
  [https://doi.org/10.3390/app12084048](https://doi.org/10.3390/app12084048)
