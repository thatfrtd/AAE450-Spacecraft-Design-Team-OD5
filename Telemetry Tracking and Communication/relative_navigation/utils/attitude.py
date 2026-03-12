import numpy as np
from utils.constants import *

def quat_multiply(p, q):
    # p, q as [x, y, z, w]
    px, py, pz, pw = p
    qx, qy, qz, qw = q
    return np.array([
        pw*qx + px*qw + py*qz - pz*qy,
        pw*qy - px*qz + py*qw + pz*qx,
        pw*qz + px*qy - py*qx + pz*qw,
        pw*qw - px*qx - py*qy - pz*qz
    ])

def kde(q, w):
    omega_quat = np.array([w[0], w[1], w[2], 0])       # pure quaternion
    q_dot = 0.5 * quat_multiply(omega_quat, q)         # your quat multiply
    return q_dot   

def dde(I, w, q, r, noise, u_tau):
    tau_g   = grav_gradient(I, w, q, r)
    w_dot = np.linalg.solve(I, -np.cross(w, I @ w) + tau_g + u_tau) + noise # solve is the same as I^-1 b
    return w_dot.flatten()

def grav_gradient(I, w, q, r):
    """
    Gravity gradient torque in the body frame.

    Args:
        I : (3,3) inertia tensor, body frame
        w : (3,)  angular velocity, body frame  (unused but kept for signature)
        q : (4,)  quaternion [qx, qy, qz, qw]
        r : (3,)  ECI position, km
    """
    # Rotate ECI position into body frame using quaternion
    r_body = eci_to_body(q, r)

    r_norm = np.linalg.norm(r_body)
    r_hat  = r_body / r_norm
    L_g    = (3.0 * MU_EARTH / r_norm**3) * np.cross(r_hat, I @ r_hat)
    return L_g

def eci_to_body(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """
    Rotate a vector from ECI to body frame using quaternion [qx, qy, qz, qw].
    """
    qx, qy, qz, qw = q
    # Rotation matrix from ECI to body  (transpose of body-to-ECI)
    R = np.array([
        [1 - 2*(qy**2 + qz**2),     2*(qx*qy + qw*qz),     2*(qx*qz - qw*qy)],
        [    2*(qx*qy - qw*qz), 1 - 2*(qx**2 + qz**2),     2*(qy*qz + qw*qx)],
        [    2*(qx*qz + qw*qy),     2*(qy*qz - qw*qx), 1 - 2*(qx**2 + qy**2)],
    ])
    return R @ v
import numpy as np

def quat2dcm(q):
    q1, q2, q3, q4 = q[0], q[1], q[2], q[3]
    T = np.array([
        [1 - 2*q2**2 - 2*q3**2,  2*(q1*q2 + q3*q4),      2*(q1*q3 - q2*q4)],
        [2*(q1*q2 - q3*q4),      1 - 2*q1**2 - 2*q3**2,  2*(q2*q3 + q1*q4)],
        [2*(q1*q3 + q2*q4),      2*(q2*q3 - q1*q4),      1 - 2*q1**2 - 2*q2**2]
    ])
    return T

def dcm2quat(R):
    trC = np.trace(R)
    epsilon4 = 0.5 * np.sqrt(1 + trC)
    epsilon123 = (1 / (4 * epsilon4)) * np.array([
        R[1, 2] - R[2, 1],
        R[2, 0] - R[0, 2],
        R[0, 1] - R[1, 0]
    ])
    q = np.concatenate([epsilon123, [epsilon4]])
    return q

