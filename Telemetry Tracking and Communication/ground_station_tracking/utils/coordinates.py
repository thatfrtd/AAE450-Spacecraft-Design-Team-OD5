"""
utils/coordinates.py

Coordinate conversions for orbital mechanics simulation.
All units are SI: metres, seconds, radians unless otherwise noted.

Frames:
    ECI  - Earth-Centered Inertial (J2000 approximation)
    ECEF - Earth-Centered Earth-Fixed
    PQW  - Perifocal (intermediate, used in Keplerian conversions)

Note on TEME vs J2000:
    SGP4 returns states in TEME (True Equator Mean Equinox). The rotation
    from TEME to J2000 is small (~arcseconds) and is treated as negligible
    here.

Unit note:
    Constants (MU_EARTH, RE_EARTH) are in METRES for consistency with the
    propagator and OrbitalState. The internal Keplerian conversion functions
    preserve whatever units are passed in via mu/a — so passing MU_EARTH in
    m^3/s^2 and a in metres produces positions in metres.
"""

import numpy as np
import math
from datetime import datetime, timezone

# ── Earth constants (SI — metres) ─────────────────────────────────────────────
MU_EARTH    = 3.986004418e14   # m^3/s^2   gravitational parameter
RE_EARTH    = 6_378_137.0      # m         equatorial radius (WGS-84)
J2          = 1.08262668e-3    # -          J2 oblateness coefficient
OMEGA_EARTH = 7.2921150e-5     # rad/s     Earth rotation rate


# ── Time / GMST ───────────────────────────────────────────────────────────────

def julian_date(dt: datetime) -> float:
    """Return the Julian Date for a UTC datetime."""
    dt = dt.replace(tzinfo=timezone.utc) if dt.tzinfo is None else dt
    j2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=timezone.utc)
    return 2_451_545.0 + (dt - j2000).total_seconds() / 86_400.0


def gmst_from_datetime(dt: datetime) -> float:
    """
    Compute Greenwich Mean Sidereal Time (GMST) in radians.
    Uses the IAU 1982 formula accurate to ~0.1 arcsec over several decades.
    """
    jd = julian_date(dt)
    T  = (jd - 2_451_545.0) / 36_525.0
    gmst_deg = (100.4606184
                + 36_000.77004 * T
                + 0.000387933 * T**2
                - T**3 / 38_710_000.0)
    ut1_seconds = (dt.hour * 3600 + dt.minute * 60
                   + dt.second + dt.microsecond * 1e-6)
    gmst_deg += 360.98564724 * ut1_seconds / 86_400.0
    return np.deg2rad(gmst_deg % 360.0)


# ── ECI ↔ ECEF ────────────────────────────────────────────────────────────────

def _rz(theta: float) -> np.ndarray:
    """3×3 rotation matrix about the Z-axis by angle theta (radians)."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[ c,  s, 0],
                     [-s,  c, 0],
                     [ 0,  0, 1]])


def eci_to_ecef(r_eci: np.ndarray, dt: datetime) -> np.ndarray:
    """
    Rotate an ECI position vector to ECEF.

    Args:
        r_eci : (3,) position vector in ECI, metres
        dt    : UTC epoch
    Returns:
        (3,) position vector in ECEF, metres
    """
    return _rz(gmst_from_datetime(dt)) @ r_eci


def ecef_to_eci(r_ecef: np.ndarray, dt: datetime) -> np.ndarray:
    """
    Rotate an ECEF position vector to ECI.

    Args:
        r_ecef : (3,) position vector in ECEF, metres
        dt     : UTC epoch
    Returns:
        (3,) position vector in ECI, metres
    """
    return _rz(gmst_from_datetime(dt)).T @ r_ecef


def eci_vel_to_ecef(r_eci: np.ndarray, v_eci: np.ndarray,
                    dt: datetime) -> np.ndarray:
    """
    Convert an ECI velocity vector to ECEF, accounting for Earth rotation.

    Args:
        r_eci : (3,) ECI position, metres
        v_eci : (3,) ECI velocity, m/s
        dt    : UTC epoch
    Returns:
        (3,) velocity in ECEF, m/s
    """
    theta    = gmst_from_datetime(dt)
    R        = _rz(theta)
    omega_vec = np.array([0.0, 0.0, OMEGA_EARTH])
    return R @ v_eci - np.cross(omega_vec, R @ r_eci)


# ── Geodetic helpers ──────────────────────────────────────────────────────────

def ecef_to_geodetic(r_ecef: np.ndarray) -> tuple:
    """
    Convert ECEF position to geodetic latitude, longitude, altitude.
    Uses Bowring's iterative method (WGS-84).

    Args:
        r_ecef : (3,) ECEF position, metres
    Returns:
        (lat_rad, lon_rad, alt_m)
    """
    x, y, z = r_ecef
    f   = 1 / 298.257223563
    b   = RE_EARTH * (1 - f)
    e2  = 1 - (b / RE_EARTH)**2
    lon = np.arctan2(y, x)
    p   = np.sqrt(x**2 + y**2)
    lat = np.arctan2(z, p * (1 - e2))

    for _ in range(10):
        sin_lat = np.sin(lat)
        N       = RE_EARTH / np.sqrt(1 - e2 * sin_lat**2)
        lat_new = np.arctan2(z + e2 * N * sin_lat, p)
        if abs(lat_new - lat) < 1e-12:
            lat = lat_new
            break
        lat = lat_new

    sin_lat = np.sin(lat)
    N   = RE_EARTH / np.sqrt(1 - e2 * sin_lat**2)
    alt = (p / np.cos(lat) - N
           if abs(np.cos(lat)) > 1e-10
           else abs(z) / sin_lat - N * (1 - e2))
    return lat, lon, alt


def geodetic_to_ecef(lat_rad: float, lon_rad: float,
                     alt_m: float) -> np.ndarray:
    """
    Convert geodetic (lat, lon, alt) to ECEF position.

    Args:
        lat_rad : geodetic latitude, radians
        lon_rad : longitude, radians
        alt_m   : altitude above WGS-84 ellipsoid, metres
    Returns:
        (3,) ECEF position, metres
    """
    f       = 1 / 298.257223563
    e2      = 2 * f - f**2
    sin_lat = np.sin(lat_rad)
    cos_lat = np.cos(lat_rad)
    N       = RE_EARTH / np.sqrt(1 - e2 * sin_lat**2)
    return np.array([
        (N + alt_m) * cos_lat * np.cos(lon_rad),
        (N + alt_m) * cos_lat * np.sin(lon_rad),
        (N * (1 - e2) + alt_m) * sin_lat,
    ])


# ── Anomaly helpers ───────────────────────────────────────────────────────────

def solve_kepler(M: float, e: float, tol: float = 1e-10,
                 max_iter: int = 50) -> float:
    """Solve Kepler's equation M = E - e·sin(E) for eccentric anomaly E."""
    M = np.mod(M, 2 * np.pi)
    E = M
    for _ in range(max_iter):
        dE = -(E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E += dE
        if abs(dE) < tol:
            return E
    raise RuntimeError("Kepler solver did not converge.")


def kepler_eqn(E: float, e: float) -> float:
    """Return mean anomaly M from eccentric anomaly E."""
    return np.mod(E - e * np.sin(E), 2 * np.pi)


def eccentric_to_true(E: float, e: float) -> float:
    """Convert eccentric anomaly E to true anomaly nu."""
    nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2),
                        math.sqrt(1 - e) * math.cos(E / 2))
    return np.mod(nu, 2 * np.pi)


def true_to_eccentric(nu: float, e: float) -> float:
    """Convert true anomaly nu to eccentric anomaly E."""
    E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2),
                       np.sqrt(1 + e) * np.cos(nu / 2))
    return np.mod(E, 2 * np.pi)


# ── Rotation matrices ─────────────────────────────────────────────────────────

def RTN_to_inert_313(RAAN: float, inc: float, aop: float) -> np.ndarray:
    """
    3-1-3 Euler rotation matrix: orbital frame → inertial frame.

    Args:
        RAAN : right ascension of ascending node, rad
        inc  : inclination, rad
        aop  : argument of perigee, rad
    Returns:
        (3,3) rotation matrix
    """
    cO, sO = np.cos(RAAN), np.sin(RAAN)
    ci, si = np.cos(inc),  np.sin(inc)
    cw, sw = np.cos(aop),  np.sin(aop)
    return np.array([
        [ cO*cw - sO*ci*sw,  -cO*sw - sO*ci*cw,  sO*si],
        [ sO*cw + cO*ci*sw,  -sO*sw + cO*ci*cw, -cO*si],
        [ si*sw,              si*cw,              ci   ],
    ])


def inert_to_RTN_313(RAAN: float, inc: float, aop: float) -> np.ndarray:
    """Inertial frame → orbital frame (transpose of RTN_to_inert_313)."""
    return RTN_to_inert_313(RAAN, inc, aop).T


# ── Cartesian ↔ Keplerian (array interface) ───────────────────────────────────
#
# These wrappers use a (6,) array interface so that OrbitalState and the
# propagator can call them as:
#     cartesian_to_keplerian(state_vector, mu)
#     keplerian_to_cartesian(kep_array, mu)
#
# Internally they delegate to the scalar functions below which preserve
# your original logic exactly.

def cartesian_to_keplerian(state: np.ndarray,
                            mu: float = MU_EARTH) -> np.ndarray:
    """
    Convert a (6,) Cartesian ECI state to classical Keplerian elements.

    Args:
        state : (6,) [rx, ry, rz, vx, vy, vz] in metres and m/s
        mu    : gravitational parameter, m^3/s^2

    Returns:
        (6,) [a (m), e, i (rad), RAAN (rad), argp (rad), nu (rad)]
    """
    r_vec = state[:3]
    v_vec = state[3:]
    a, e, inc, RAAN, aop, M = _cartesian_to_keplerian_scalar(r_vec, v_vec, mu)
    # Convert mean anomaly back to true anomaly for the return value
    E  = solve_kepler(M, e)
    nu = eccentric_to_true(E, e)
    return np.array([a, e, inc, RAAN, aop, nu])


def keplerian_to_cartesian(kep: np.ndarray,
                            mu: float = MU_EARTH) -> np.ndarray:
    """
    Convert classical Keplerian elements to a (6,) Cartesian ECI state.

    Args:
        kep : (6,) [a (m), e, i (rad), RAAN (rad), argp (rad), nu (rad)]
              true anomaly expected (not mean anomaly)
        mu  : gravitational parameter, m^3/s^2

    Returns:
        (6,) [rx, ry, rz, vx, vy, vz] in metres and m/s
    """
    a, e, inc, RAAN, aop, nu = kep
    pos, vel = _keplerian_to_cartesian_scalar(a, e, inc, RAAN, aop, nu, mu,
                                              anomaly_type="true")
    return np.concatenate([pos, vel])


# ── Cartesian ↔ Keplerian (scalar interface — your original functions) ─────────

def _keplerian_to_cartesian_scalar(a, e, inc, RAAN, aop, anomaly, mu,
                                   frame="inert", anomaly_type="mean"):
    """
    Convert Keplerian elements to Cartesian coordinates (scalar interface).
    Preserves original logic exactly.
    """
    if anomaly_type == "mean":
        E  = solve_kepler(anomaly, e)
        nu = eccentric_to_true(E, e)
    elif anomaly_type == "true":
        nu = anomaly
        E  = true_to_eccentric(anomaly, e)
    else:
        raise ValueError("anomaly_type must be 'mean' or 'true'")

    r_c   = a * (1 - e * math.cos(E))
    pos_o = r_c * np.array([math.cos(nu), math.sin(nu), 0.0])
    vel_o = (math.sqrt(mu / a) / (1 - e * math.cos(E))) * np.array([
        -math.sin(E),
        math.sqrt(1 - e**2) * math.cos(E),
        0.0,
    ])

    if frame == "inert":
        R     = RTN_to_inert_313(RAAN, inc, aop)
        pos_i = R @ pos_o
        vel_i = R @ vel_o
    else:
        pos_i = pos_o
        vel_i = vel_o

    return pos_i, vel_i


def _cartesian_to_keplerian_scalar(r_vec, v_vec, mu):
    """
    Convert Cartesian ECI state to Keplerian elements (scalar interface).
    Returns (a, e, inc, RAAN, aop, M) — note M is mean anomaly.
    Preserves original logic exactly.
    """
    h_vec = np.cross(r_vec, v_vec)
    r_hat = r_vec / np.linalg.norm(r_vec)
    e_vec = (np.cross(v_vec, h_vec) / mu) - r_hat

    k_vec = np.array([0.0, 0.0, 1.0])
    n_vec = np.cross(k_vec, h_vec)
    n     = np.linalg.norm(n_vec)
    e     = np.linalg.norm(e_vec)

    cos_nu = np.dot(e_vec, r_vec) / (e * np.linalg.norm(r_vec))
    nu     = math.acos(np.clip(cos_nu, -1.0, 1.0))
    if np.dot(r_vec, v_vec) < 0:
        nu = 2 * np.pi - nu

    E    = true_to_eccentric(nu, e)
    M    = kepler_eqn(E, e)
    inc  = np.arccos(h_vec[2] / np.linalg.norm(h_vec))
    RAAN = np.arccos(np.clip(n_vec[0] / n, -1.0, 1.0))
    if n_vec[1] < 0:
        RAAN = 2 * np.pi - RAAN

    cos_aop = np.dot(n_vec, e_vec) / (n * e)
    aop     = np.arccos(np.clip(cos_aop, -1.0, 1.0))
    if e_vec[2] < 0:
        aop = 2 * np.pi - aop

    a = 1.0 / (2.0 / np.linalg.norm(r_vec) - np.linalg.norm(v_vec)**2 / mu)
    return a, e, inc, RAAN, aop, M


# ── Angle utilities ───────────────────────────────────────────────────────────

def wrap_to_2pi(angle: float) -> float:
    """Wrap an angle to [0, 2π)."""
    return angle % (2.0 * np.pi)


def wrap_to_pi(angle: float) -> float:
    """Wrap an angle to (-π, π]."""
    return (angle + np.pi) % (2.0 * np.pi) - np.pi