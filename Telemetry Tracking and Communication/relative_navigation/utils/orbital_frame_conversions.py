import numpy as np
import math
from utils.kepler import *
from utils.constants import *

import numpy as np
from dataclasses import dataclass

@dataclass
class KeplerianElements:
    a     : float   # semi-major axis, km
    e     : float   # eccentricity
    i     : float   # inclination, rad
    raan  : float   # right ascension of ascending node, rad
    argp  : float   # argument of perigee, rad
    nu    : float   # true anomaly, rad
    # ── derived ──────────────────────
    M     : float   # mean anomaly, rad
    n     : float   # mean motion, rad/s

def tle_to_keplerian(line1: str, line2: str) -> KeplerianElements:
    """
    Parse a TLE and return classical Keplerian orbital elements.

    Args:
        line1 : TLE line 1 string
        line2 : TLE line 2 string

    Returns:
        KeplerianElements dataclass (SI units, radians)
    """
    # ── Parse line 2 ──────────────────────────────────────────────────────────
    i_deg    = float(line2[8:16])
    raan_deg = float(line2[17:25])
    e        = float("0." + line2[26:33].strip())   # implicit decimal point
    argp_deg = float(line2[34:42])
    M_deg    = float(line2[43:51])
    n_rev    = float(line2[52:63])                  # mean motion, rev/day

    # ── Unit conversions ──────────────────────────────────────────────────────
    deg2rad = np.pi / 180.0

    i    = i_deg    * deg2rad
    raan = raan_deg * deg2rad
    argp = argp_deg * deg2rad
    M    = M_deg    * deg2rad
    n    = n_rev * 2.0 * np.pi / 86400.0            # rad/s

    # ── Semi-major axis from mean motion (Kepler's 3rd law) ───────────────────
    a = (MU_EARTH / n**2) ** (1.0 / 3.0)

    # ── Mean anomaly → eccentric anomaly (Newton's method) ───────────────────
    E = solve_kepler(M, e)

    # ── Eccentric anomaly → true anomaly ─────────────────────────────────────
    nu = 2.0 * np.arctan2(
        np.sqrt(1.0 + e) * np.sin(E / 2.0),
        np.sqrt(1.0 - e) * np.cos(E / 2.0),
    )
    nu = nu % (2.0 * np.pi)

    return KeplerianElements(a=a, e=e, i=i, raan=raan, argp=argp, nu=nu, M=M, n=n)


def keplerian_to_cartesian(
    elements: KeplerianElements,
    mu: float = MU_EARTH,
    anomaly_type: str = "mean",
    frame: str = "inert",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert Keplerian orbital elements to Cartesian coordinates.

    Args:
        elements     : KeplerianElements dataclass (km, rad)
        mu           : gravitational parameter, km^3/s^2
        anomaly_type : 'mean' — use elements.M  |  'true' — use elements.nu
        frame        : 'inert' — rotate to ECI  |  else — perifocal frame

    Returns:
        pos : (3,) ECI position,  km
        vel : (3,) ECI velocity,  km/s
    """
    a, e, inc = elements.a, elements.e, elements.i
    RAAN, aop = elements.raan, elements.argp

    if anomaly_type == "mean":
        E  = solve_kepler(elements.M, e)
        nu = eccentric_to_true(E, e)
    elif anomaly_type == "true":
        nu = elements.nu
        E  = true_to_eccentric(nu, e)
    else:
        raise ValueError("anomaly_type must be 'mean' or 'true'")

    r_c   = a * (1.0 - e * math.cos(E))
    pos_o = r_c * np.array([math.cos(nu), math.sin(nu), 0.0])
    vel_o = (math.sqrt(mu / a) / (1.0 - e * math.cos(E))) * np.array(
        [-math.sin(E), math.sqrt(1.0 - e**2) * math.cos(E), 0.0]
    )

    if frame == "inert":
        dcm   = RTN_to_inert_313(RAAN, inc, aop)
        pos_i = dcm @ pos_o
        vel_i = dcm @ vel_o
    else:
        pos_i, vel_i = pos_o, vel_o

    return pos_i, vel_i


def cartesian_to_keplerian(
    r_vec: np.ndarray,
    v_vec: np.ndarray,
    mu: float = MU_EARTH,
) -> KeplerianElements:
    """
    Convert Cartesian ECI coordinates to Keplerian orbital elements.

    Args:
        r_vec : (3,) ECI position,  km
        v_vec : (3,) ECI velocity,  km/s
        mu    : gravitational parameter, km^3/s^2

    Returns:
        KeplerianElements dataclass (km, rad)
    """
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)

    h_vec = np.cross(r_vec, v_vec)
    e_vec = (np.cross(v_vec, h_vec) / mu) - (r_vec / r)
    e     = np.linalg.norm(e_vec)

    n_vec = np.cross(np.array([0.0, 0.0, 1.0]), h_vec)
    n     = np.linalg.norm(n_vec)

    cos_nu = np.dot(e_vec, r_vec) / (e * r)
    nu     = math.acos(np.clip(cos_nu, -1.0, 1.0))
    if np.dot(r_vec, v_vec) < 0:
        nu = 2.0 * math.pi - nu

    E   = true_to_eccentric(nu, e)
    M   = kepler_eqn(E, e)
    inc = math.acos(np.clip(h_vec[2] / np.linalg.norm(h_vec), -1.0, 1.0))

    raan = math.acos(np.clip(n_vec[0] / n, -1.0, 1.0))
    if n_vec[1] < 0:
        raan = 2.0 * math.pi - raan

    cos_aop = np.dot(n_vec, e_vec) / (n * e)
    argp    = math.acos(np.clip(cos_aop, -1.0, 1.0))
    if e_vec[2] < 0:
        argp = 2.0 * math.pi - argp

    a        = 1.0 / (2.0 / r - v**2 / mu)
    n_motion = math.sqrt(mu / a**3)

    return KeplerianElements(
        a=a, e=e, i=inc, raan=raan, argp=argp, nu=nu, M=M, n=n_motion
    )

def equinoctial_to_cartesian(p, f, g, h, k, L, mu, frame="inert"):
    """
    Convert equinotical orbital elements to keplerian orbit elments. 
    Parameters : 
    p : float
        semi-latus rectum [km]
    f : float

    g : float

    h : float

    k : float 

    L : float

    mu : float
        Gravitational Parameter
    
        
    """
    # a, e, inc, aop, RAAN, nu = equinoctal_to_keplerian(p, f, g, h, k, L)
    alpha2 = h**2 - k**2
    s2 = 1 + h**2 + k**2
    w = 1 + f * np.cos(L) + g * np.sin(L)
    r = p / w

    # position vector r
    r_vec = np.array([
        (r / s2) * (np.cos(L) + alpha2 * np.cos(L) + 2*h*k*np.sin(L)),
        (r / s2) * (np.sin(L) - alpha2 * np.sin(L) + 2*h*k*np.cos(L)),
        (2*r / s2) * (h*np.sin(L) - k*np.cos(L))
    ])

    # velocity vector v
    v_vec = np.array([
        -(1 / s2) * np.sqrt(mu / p) *
        (np.sin(L) + alpha2*np.sin(L) - 2*h*k*np.cos(L) + g - 2*f*h*k + alpha2*g),

        -(1 / s2) * np.sqrt(mu / p) *
        (-np.cos(L) + alpha2*np.cos(L) + 2*h*k*np.sin(L) - f + 2*g*h*k + alpha2*f),

        (2 / s2) * np.sqrt(mu / p) *
        (h*np.cos(L) + k*np.sin(L) + f*h + g*k)
    ])
    
    return r_vec, v_vec

def cartesian_to_equinoctial(r_vec, v_vec, mu):
    a, e, inc, RAAN, aop, M = cartesian_to_keplerian(r_vec, v_vec, mu)
    # Semi-parameter
    p = a * (1 - e**2)
    f = e * np.cos(aop + RAAN)
    g = e * np.sin(aop + RAAN)
    h = np.tan(inc/2) * np.cos(RAAN)
    k = np.tan(inc/2) * np.sin(RAAN)
    # Convert Mean Anomaly to True Anomaly
    E = solve_kepler(M, e)
    nu = eccentric_to_true(E, e) 
    # True Longitude
    L = RAAN + aop + nu
    L = np.mod(L, 2 * np.pi)
    return p, f, g, h, k, L


def milankovitch_to_cartesian(h_vec, e_vec, L, mu, frame="inert"):
    # Node Vector
    k_vec = np.array([0.0, 0.0, 1.0])
    n_vec = np.cross(k_vec, h_vec)
    n = np.linalg.norm(n_vec)

    # Find orbit eccentricity
    e = np.linalg.norm(e_vec)

    # Angular Momentum
    h = np.linalg.norm(h_vec)

    #  Inclination using inverse dcm mapping
    inc = np.acos(h_vec[2] / h)

    # RAAN
    if n > 1e-12:
        RAAN = np.arctan2(n_vec[1], n_vec[0])
    else:
        RAAN = 0.0
    RAAN = np.mod(RAAN, 2*np.pi)

    # Argument of periapsis
    if e > 1e-10 and n > 1e-10:
        aop = np.atan2(
            np.dot(np.cross(n_vec, e_vec), h_vec) / (n * h),
            np.dot(n_vec, e_vec) / n
        )
    else:
        aop = 0.0
    aop = np.mod(aop, 2*np.pi)


    # True anomaly
    nu = L - RAAN - aop
    nu = np.mod(nu, 2*np.pi)

    # Eccentric Anomaly
    E = true_to_eccentric(nu, e)

    # Mean Anomaly
    M = kepler_eqn(E, e)

    # Semi-latus rectum
    p = (np.linalg.norm(h_vec)**2) / mu

    # Semi-major axis
    a = p / (1 - e**2)

    r_vec, v_vec = keplerian_to_cartesian(a, e, inc, RAAN, aop, M, mu, frame=frame, anomaly_type="mean")

    return r_vec, v_vec


def cartesian_to_milankovitch(r_vec, v_vec, mu): 
    
    a, e, inc, RAAN, aop, M = cartesian_to_keplerian(r_vec, v_vec, mu)
    
    # Calculate Angular Momentum Vector
    h_vec = np.cross(r_vec, v_vec)

    # Calculate eccentricity vector
    r_hat = r_vec / np.linalg.norm(r_vec)
    e_vec = (np.cross(v_vec, h_vec) / mu) - r_hat

    # Convert Mean Anomaly to True Anomaly
    E = solve_kepler(M, e)
    nu = eccentric_to_true(E, e) 
    
    # true longitude
    L = RAAN + aop + nu
    L = np.mod(L, 2 * np.pi)
    
    return h_vec, e_vec, L



def RTN_to_inert_313(RAAN, inc, aop):
    """
    Convert 313 Euler angles (RAAN, inc, aop) to rotation matrix. Rotates orbital frame to intertial frame. 

    Parameters:
    RAN, inc, aop : float
        Euler angles in radians
        RAAN : RAAN, first rotation about z-axis
        inc: Inclination, rotation about x-axis
        aop: AOP, second rotation about z-axis

    Returns:
    R : 3x3 numpy array
        Rotation matrix
    """
    cO = np.cos(RAAN)
    sO = np.sin(RAAN)
    ci = np.cos(inc)
    si = np.sin(inc)
    cw = np.cos(aop)
    sw = np.sin(aop)

    R = np.array([
        [cO*cw - sO*ci*sw, -cO*sw - sO*ci*cw,  sO*si],
        [sO*cw + cO*ci*sw, -sO*sw + cO*ci*cw, -cO*si],
        [si*sw,             si*cw,             ci]
    ])
    return R


def inert_to_RTN_313(RAAN, inc, aop):
    R = RTN_to_inert_313(RAAN, inc, aop).T
    return R

