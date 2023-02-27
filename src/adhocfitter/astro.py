import astropy.constants as const
import astropy.units as u
import numpy as np
import scipy.special as special

AU_IN_METERS = 149597870700.
SOLAR_RADIUS_IN_METERS = 695700000.
AU_OVER_SOLAR_RADIUS = AU_IN_METERS / SOLAR_RADIUS_IN_METERS
NEWTON_GRAV_CONST = 6.67430e-11  # m^3 kg^-1 s^-2
NOMINAL_SOLAR_MASS_PARAM = 1.3271244e20  # m^3 s^-2
NOMINAL_SOLAR_MASS_CORRECTION = (NOMINAL_SOLAR_MASS_PARAM /
                                 ((2 * np.pi) **2
                                  * AU_IN_METERS ** 3
                                  / (365.25 * 24 * 3600) ** 2))
SOLAR_DENSITY = (NOMINAL_SOLAR_MASS_PARAM / NEWTON_GRAV_CONST
                 / SOLAR_RADIUS_IN_METERS ** 3 / (4 * np.pi / 3))  # kg m^-3
SOLAR_TEFF = 5772.  # Kelvin
EARTH_EQ_RADIUS_IN_METERS = 6378100.
JUPITER_EQ_RADIUS_IN_METERS = 7.1492e7
NOMINAL_EARTH_MASS_PARAM = 3.986004e14  # m^3/s^2
NOMINAL_JUPITER_MASS_PARAM =  1.2668653e17  # m^3/s^2
SOLAR_MASS_IN_EARTH_MASS = NOMINAL_SOLAR_MASS_PARAM / NOMINAL_EARTH_MASS_PARAM
SOLAR_MASS_IN_JUPITER_MASS = NOMINAL_SOLAR_MASS_PARAM / NOMINAL_JUPITER_MASS_PARAM
DAY_IN_SECONDS = 24 * 60 * 60

STEFAN_BOLTZMAN = 5.67037441918442945397099673188923087584012297029130e-8  # J m^-2 s^-1 K^-4
FLUX_CGS_IN_SI = 1000
SUN_ABS_MAG = 4.74

def calculate_stellar_radius(stellar_mass, logg):
    """Calculate stellar radius from scaled semimajor axis, period, and stellar mass.

    Args:
        stellar_mass: Stellar mass in nominal solar mass.
        logg: Base-10 logarithm of surface gravity in cgs units.
    """
    surface_g = 10. ** logg
    stellar_radius = (
        np.sqrt(stellar_mass * NOMINAL_SOLAR_MASS_PARAM * 1e6 / surface_g)
        / (SOLAR_RADIUS_IN_METERS * 1e2))
    return stellar_radius

def calculate_stellar_logg(stellar_mass, stellar_radius):
    """Calculate surface gravity in dex(cgs). Inputs are in solar units."""
    surface_g = (
        NOMINAL_SOLAR_MASS_PARAM * stellar_mass
        / (SOLAR_RADIUS_IN_METERS * stellar_radius)**2
        * 1e2)  # cm/s^2.
    return np.log10(surface_g)

def calculate_semimajor_axis_au(stellar_mass, period):
    """Calculate the semimajor axis in au.

    Args:
        stellar_mass: Stellar mass in units of nominal solar mass.
        period: Planet orbital period in days.
    """
    return (stellar_mass * NOMINAL_SOLAR_MASS_CORRECTION * (period / 365.25) ** 2) ** (1/3.)

def calculate_planet_radius(planet_radius_ratio, stellar_radius):
    return planet_radius_ratio * stellar_radius * SOLAR_RADIUS_IN_METERS / JUPITER_EQ_RADIUS_IN_METERS

def calculate_min_planet_mass(rv_semiamp, ecc, period, stellar_mass):
    """rv_semiamp in m/s, period in days, stellar_mass in solar mass -> m_p sin(i) in Jupiter mass"""
    mp_sini = (rv_semiamp * np.sqrt(1 - ecc**2)
               * NOMINAL_JUPITER_MASS_PARAM ** (-1/3.)
               * (stellar_mass * SOLAR_MASS_IN_JUPITER_MASS) ** (2/3.)
               * (period * DAY_IN_SECONDS / (2 * np.pi)) ** (1/3.))
    return mp_sini

def calculate_min_planet_mass_earth(rv_semiamp, ecc, period, stellar_mass):
    """rv_semiamp in m/s, period in days, stellar_mass in solar mass -> m_p sin(i) in Earth mass"""
    mp_sini = (rv_semiamp * np.sqrt(1 - ecc**2)
               * NOMINAL_EARTH_MASS_PARAM ** (-1/3.)
               * (stellar_mass * SOLAR_MASS_IN_EARTH_MASS) ** (2/3.)
               * (period * DAY_IN_SECONDS / (2 * np.pi)) ** (1/3.))
    return mp_sini

def calculate_rv_semiamp(min_planet_mass, ecc, period, stellar_mass):
    """Calculate the RV semiamplitude in m/s.

    Args:
        min_planet_mass (float): Minimum planet mass (Mp sin i) in units of
            nominal Solar mass.
        ecc (float): Planet's eccentricity.
        period (float): Planet's orbital period.
        stellar_mass (float): Stellar mass in units of nominal Solar mass.
    """
    rv_semiamp = (min_planet_mass * NOMINAL_SOLAR_MASS_PARAM ** (1/3.)
                  * stellar_mass ** (-2/3.)
                  * (period * DAY_IN_SECONDS / (2 * np.pi)) ** (-1/3.)
                  / np.sqrt(1 - ecc ** 2))
    return rv_semiamp

def calculate_stellar_density(stellar_mass, logg):
    """Calculate stellar density in cgs units.

    Args:
        stellar_mass (float): Stellar mass in units of nominal solar mass.
        logg (flot): Base-10 logarithm of surface gravity in cgs units."""
    stellar_radius = calculate_stellar_radius(stellar_mass, logg)
    stellar_density = stellar_mass / (stellar_radius ** 3) * SOLAR_DENSITY * 1e-3
    return stellar_density

def calculate_aor_stellar_density(stellar_density, period):
    aor = (stellar_density * NOMINAL_SOLAR_MASS_CORRECTION
           * (period / 365.25) ** 2) ** (1/3.) * AU_OVER_SOLAR_RADIUS
    return aor

def calculate_aor_stellar_radius(stellar_mass, stellar_radius, period):
    stellar_density = stellar_mass / (stellar_radius ** 3)
    aor = (stellar_density * NOMINAL_SOLAR_MASS_CORRECTION
           * (period / 365.25) ** 2) ** (1/3.) * AU_OVER_SOLAR_RADIUS
    return aor

def calculate_aor(stellar_mass, logg, period):
    """Calculate scaled semimajor axis.

    Args:
        stellar_mass: Stellar mass in units of nominal solar mass.
        logg: Base-10 logarithm of surface gravity in cgs units.
        period: Planet orbital period in days.
    """
    stellar_radius = calculate_stellar_radius(stellar_mass, logg)
    semimajor = calculate_semimajor_axis_au(stellar_mass, period)
    aor = semimajor * AU_OVER_SOLAR_RADIUS / stellar_radius
    return aor

def calculate_aor_err(aor, stellar_mass, stellar_mass_err, stellar_radius, stellar_radius_err, period, period_err):
    aor_var = aor ** 2 * ((stellar_mass_err / stellar_mass) ** 2
                          + 3 * (stellar_radius_err / stellar_radius) ** 2
                          + 2 * (period_err / period) ** 2)
    return np.sqrt(aor_var)

def calculate_transit_duration(period, impact_param, scaled_semimajor_axis, secosw, sesinw, scaled_planet_radius):
    ecc = calculate_eccentricity(secosw, sesinw)
    esinw = sesinw * np.sqrt(ecc)
    return ((period / np.pi) * np.sqrt(1 - ecc ** 2) / (1 + esinw)
            * np.arcsin(np.sqrt(((1 + scaled_planet_radius) ** 2 - impact_param ** 2)
                                / (scaled_semimajor_axis ** 2 - impact_param ** 2))))

def calculate_transit_duration_2(period, impact_param, scaled_semimajor_axis, ecc, sin_omega, scaled_planet_radius):
    esinw = sin_omega * ecc
    return ((period / np.pi) * np.sqrt(1 - ecc ** 2) / (1 + esinw)
            * np.arcsin(np.sqrt(((1 + scaled_planet_radius) ** 2 - impact_param ** 2)
                                / (scaled_semimajor_axis ** 2 - impact_param ** 2))))

def calculate_inclination(impact_param, a_over_r_star):
    return np.degrees(np.arccos(impact_param / a_over_r_star))

def calculate_eccentricity(secosw, sesinw):
    return secosw ** 2 + sesinw ** 2

def calculate_omega(secosw, sesinw):
    return np.degrees(np.arctan2(sesinw, secosw))

def calculate_planet_density(planet_mass, planet_radius):
    mp = planet_mass * NOMINAL_JUPITER_MASS_PARAM / NEWTON_GRAV_CONST
    rp = planet_radius * JUPITER_EQ_RADIUS_IN_METERS
    planet_volume = (4/3 * np.pi) * rp ** 3
    return mp / planet_volume * 1e-3  # Convert to cgs units.

def calculate_temperature_eq(star_teff, aor, bond_albedo, ecc=0.):
    elliptic_correction = (np.sqrt(1 - ecc) * special.ellipe(2 * ecc / (ecc - 1))
                           + np.sqrt(1 + ecc) * special.ellipe(2 * ecc / (ecc + 1))) / np.pi
    return star_teff / np.sqrt(2 * aor) * elliptic_correction * (1 - bond_albedo) ** (1/4.)

def calculate_temperature_eq_flux_avg(star_teff, aor, bond_albedo, ecc=0., heat_dist=1.):
    flux_ecc_correction = np.sqrt(1 - ecc * ecc) ** (-1/8.)
    return star_teff / np.sqrt(2 * aor) * flux_ecc_correction * ((1 - bond_albedo) / heat_dist) ** (1/4.)

def calculate_irradiation(star_teff, aor, ecc=0.):
    irradiation = STEFAN_BOLTZMAN * star_teff ** 4 / (aor ** 2 * np.sqrt(1 - ecc * ecc))
    return irradiation * FLUX_CGS_IN_SI

def calculate_irradiation_err(irradiation, teff, teff_err, aor, aor_err):
    irrad_var = ((4 * (teff_err / teff) ** 2 + 2 * (aor_err / aor) ** 2)
                 * irradiation ** 2)
    return np.sqrt(irrad_var)

def calculate_luminosity(star_teff, stellar_radius):
    return stellar_radius ** 2 * (star_teff / SOLAR_TEFF) ** 4.

def calculate_mag_bolometric(luminosity):
    """Luminosity in solar units"""
    return -2.5 * np.log10(luminosity) + SUN_ABS_MAG
