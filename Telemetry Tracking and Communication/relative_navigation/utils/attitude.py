import numpy as np
from utils.dynamics import grav_gradient

def quat_times(q1:np.ndarray, q2:np.ndarray) -> np.ndarray:
    q3 = 1
    return q3

def kde(q, w):
    ep13  = q[0:3]
    ep4   = q[3]
    q_dot_13 = 0.5 * (ep4 * w + np.cross(ep13, w))
    q_dot_4  = -0.5 * np.dot(ep13, w)
    return np.append(q_dot_13, q_dot_4)   # → (4,)

def dde(I, w, q, r, noise):
    tau_g   = grav_gradient(I, w, q, r)
    w_dot = np.linalg.solve(I, -np.cross(w, I @ w) + tau_g) + noise # solve is the same as I^-1 b
    return w_dot.flatten()

