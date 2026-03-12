import numpy as np
from scipy.interpolate import interp1d

def make_noise_interpolator(t_span, dt, sigma, size=3):
    """Pre-generate noise on a fixed grid, interpolate during integration."""
    t_grid = np.arange(t_span[0], t_span[1], dt)
    noise  = np.random.normal(0, sigma, size=(len(t_grid), size))
    return interp1d(t_grid, noise, axis=0, bounds_error=False, fill_value=0.0)
