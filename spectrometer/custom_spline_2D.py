"""Natural cubic splines for evaluating fields on a rectangular grid.

The implementation is the reusable form of the usual ``A @ f'' = b``
derivation: second derivatives are zero at the two outer boundaries, the
tridiagonal system supplies the interior second derivatives, and every
interval is stored as a local cubic polynomial.
"""

import numpy as np


def P_tag_tag_fun(vec, f):
    """Find the natural-spline second derivatives (the user's ``p''``).

    This is the reusable version of ``A @ p_tag_tag = b`` from the original
    one-dimensional code.  Axis 0 is the interpolation axis; any remaining
    axes are independent data series that are solved in the same operation.
    """
    vec = np.asarray(vec, dtype=float)
    f = np.asarray(f, dtype=float)

    if vec.ndim != 1 or len(vec) < 2:
        raise ValueError("vec must be a one-dimensional array of length >= 2")
    if f.shape[0] != len(vec):
        raise ValueError("f.shape[0] must equal len(vec)")

    intervals = len(vec) - 1
    delta = np.diff(vec)
    if np.any(delta <= 0.0):
        raise ValueError("spline nodes must be strictly increasing")

    # Natural boundary conditions: p''(x_0) = p''(x_N) = 0.
    p_tag_tag = np.zeros_like(f, dtype=float)
    matrix_size = intervals - 1
    if matrix_size == 0:
        return p_tag_tag

    # For an equally spaced grid this is exactly the original tridiagonal
    # matrix with 4*Delta_x on its main diagonal and Delta_x beside it.
    main_diag = 2.0 * (delta[:-1] + delta[1:])
    A = np.diag(main_diag)
    if matrix_size > 1:
        off_diag = delta[1:-1]
        A += np.diag(off_diag, 1) + np.diag(off_diag, -1)

    delta_shape = (intervals,) + (1,) * (f.ndim - 1)
    slopes = np.diff(f, axis=0) / delta.reshape(delta_shape)
    b = 6.0 * np.diff(slopes, axis=0)

    # Reshaping lets one matrix solve handle all rows/columns of a 2-D field.
    solved = np.linalg.solve(A, b.reshape(matrix_size, -1))
    p_tag_tag[1:-1] = solved.reshape(p_tag_tag[1:-1].shape)
    return p_tag_tag


def Params_finder(vec, f, p_tag_tag):
    """Return the original ``alpha, beta, gamma, etha`` spline parameters.

    On interval i the spline is written exactly as in the original project::

        alpha[i] * (x - vec[i])**3
      + beta[i]  * (x - vec[i + 1])**3
      + gamma[i] * (x - vec[i])
      + etha[i]  * (x - vec[i + 1])

    The formulas also support unequal node spacing and array-valued ``f``.
    """
    vec = np.asarray(vec, dtype=float)
    f = np.asarray(f, dtype=float)
    p_tag_tag = np.asarray(p_tag_tag, dtype=float)

    if f.shape != p_tag_tag.shape or f.shape[0] != len(vec):
        raise ValueError("f and p_tag_tag must match vec along axis 0")

    delta = np.diff(vec)
    delta_shape = (len(delta),) + (1,) * (f.ndim - 1)
    delta = delta.reshape(delta_shape)

    alpha = p_tag_tag[1:] / (6.0 * delta)
    beta = -p_tag_tag[:-1] / (6.0 * delta)
    gamma = (-p_tag_tag[1:] * delta**2 + 6.0 * f[1:]) / (6.0 * delta)
    etha = (p_tag_tag[:-1] * delta**2 - 6.0 * f[:-1]) / (6.0 * delta)

    # Coefficient index 0..3 means alpha, beta, gamma, etha.
    return np.stack((alpha, beta, gamma, etha), axis=1)


def _natural_cubic_coefficients(nodes, values):
    """Build spline parameters using the two original derivation steps."""
    p_tag_tag = P_tag_tag_fun(nodes, values)
    return Params_finder(nodes, values, p_tag_tag)


def _spline_basis(points, nodes, indices):
    """Return the four factors multiplying alpha, beta, gamma and etha."""
    left_distance = points - nodes[indices]
    right_distance = points - nodes[indices + 1]
    return np.stack(
        (
            left_distance**3,
            right_distance**3,
            left_distance,
            right_distance,
        ),
        axis=1,
    )


class NaturalCubicSpline1D:
    """Natural cubic interpolation of scalar or vector-valued samples."""

    def __init__(self, nodes, values):
        self.nodes = np.asarray(nodes, dtype=float)
        self.coefficients = _natural_cubic_coefficients(self.nodes, values)

    def __call__(self, x, clip=True):
        x_array = np.asarray(x, dtype=float)
        flat_x = x_array.ravel()
        if clip:
            flat_x = np.clip(flat_x, self.nodes[0], self.nodes[-1])

        indices = np.searchsorted(self.nodes, flat_x, side="right") - 1
        indices = np.clip(indices, 0, len(self.nodes) - 2)
        selected = self.coefficients[indices]

        basis = _spline_basis(flat_x, self.nodes, indices)
        basis = basis.reshape(basis.shape + (1,) * (selected.ndim - 2))
        result = np.sum(selected * basis, axis=1)
        output_shape = x_array.shape + self.coefficients.shape[2:]
        result = result.reshape(output_shape)
        return result.item() if result.ndim == 0 else result


class NaturalCubicSpline2D:
    """Tensor-product natural cubic spline for ``field[y_index, x_index]``."""

    def __init__(self, x_axis, y_axis, field):
        self.x = np.asarray(x_axis, dtype=float)
        self.y = np.asarray(y_axis, dtype=float)
        field = np.asarray(field, dtype=float)
        if field.shape != (len(self.y), len(self.x)):
            raise ValueError(
                "field shape must be (len(y_axis), len(x_axis)); "
                f"got {field.shape}"
            )

        # First form the four x-polynomial coefficients at every y node.
        # Shape after moveaxis: (Ny, Nx-1, 4).
        x_coefficients = _natural_cubic_coefficients(self.x, field.T)
        x_coefficients = np.moveaxis(x_coefficients, -1, 0)

        # Spline each x coefficient in y.  Final shape is
        # (Ny-1, 4_y, Nx-1, 4_x), ready for constant-time point queries.
        self.coefficients = _natural_cubic_coefficients(self.y, x_coefficients)

    def __call__(self, x, y, clip=True):
        x_array, y_array = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        )
        flat_x = x_array.ravel()
        flat_y = y_array.ravel()
        if clip:
            flat_x = np.clip(flat_x, self.x[0], self.x[-1])
            flat_y = np.clip(flat_y, self.y[0], self.y[-1])

        ix = np.searchsorted(self.x, flat_x, side="right") - 1
        iy = np.searchsorted(self.y, flat_y, side="right") - 1
        ix = np.clip(ix, 0, len(self.x) - 2)
        iy = np.clip(iy, 0, len(self.y) - 2)

        px = _spline_basis(flat_x, self.x, ix)
        py = _spline_basis(flat_y, self.y, iy)

        local = self.coefficients[iy, :, ix, :]
        result = np.einsum("ni,nij,nj->n", py, local, px)
        result = result.reshape(x_array.shape)
        return result.item() if result.ndim == 0 else result

