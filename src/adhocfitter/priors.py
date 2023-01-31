import numpy as np
import numpy.linalg as linalg

LOG_TWO_PI = np.log(2 * np.pi)

def calculate_abs_mag(apparent_mag, parallax):
    """Calculate the absolute magnitude.

    Args:
        apparent_mag: Apparent magnitude, in arbitrary photometric system.
        parallax: Measured parallax of the object in milliarcseconds.
    """
    return apparent_mag + 5. * (np.log10(parallax) - 3. + 1.)

def calculate_distance_modulus_err(parallax, parallax_err):
    """Calculate the standard error of the distance modulus from parallax.

    Args:
        parallax: In milliarcseconds.
        parallax_err: In milliarcseconds.
    """
    return 5. * parallax_err / parallax / np.log(10)

def calculate_magnitude_err(flux, flux_err):
    return 2.5 * flux_err / flux / np.log(10)

def calculate_abs_mag_err_squared(flux, flux_err, parallax, parallax_err):
    """Calculate the standard error of the absolute magnitude.

    Args:
        flux: Measured flux of the object.
        flux_err: Standard error of the measured flux.
        parallax: Measured parallax of the object in milliarcseconds.
        parallax_err: Standard error of the measured parallax in
            milliarcseconds.
    """
    return (
        ((5. * parallax_err / parallax) ** 2 + (2.5 * flux_err / flux) ** 2)
        / (np.log(10) ** 2))

def log_prob_gaussian(predict, observe, variance):
    return -0.5 * np.sum(
        (predict - observe) ** 2 / variance
        + np.log(variance) + LOG_TWO_PI)

def multivariate_gaussian_normalize(cov):
    ndim = len(cov)
    return -0.5 * np.log(linalg.det(cov)) - 0.5 * ndim * LOG_TWO_PI

def log_prob_multivariate_gaussian(predict, observe, cov, normalize):
    delta = predict - observe
    ln_like = -0.5 * delta @ linalg.solve(cov, delta)
    return ln_like + normalize

def multivariate_gaussian_prior(mean, cov):
    ndim = len(mean)
    normalize = -0.5 * np.log(linalg.det(cov)) - 0.5 * ndim * LOG_TWO_PI
    def ln_prob(params):
        delta = params - mean
        ln_like = -0.5 * delta @ linalg.solve(cov, delta)
        return ln_like + normalize
    return ln_prob
