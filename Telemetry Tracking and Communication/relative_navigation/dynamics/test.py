import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from dynamics.attitude import quat2dcm, dcm2quat


def test_quat_dcm_roundtrip():
    # --- Unique starting quaternion (not identity, not axis-aligned) ---
    q_orig = np.array([0.4, 0.5, -0.6, np.sqrt(1 - 0.4**2 - 0.5**2 - 0.6**2)])
    assert abs(np.linalg.norm(q_orig) - 1.0) < 1e-10, "Input quaternion must be unit"

    # Step 1: quat -> DCM
    R = quat2dcm(q_orig)

    # Sanity checks on DCM
    assert R.shape == (3, 3), "DCM must be 3x3"
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-10), "DCM must be orthogonal (R R^T = I)"
    assert abs(np.linalg.det(R) - 1.0) < 1e-10, "DCM determinant must be +1 (proper rotation)"

    # Step 2: DCM -> quat
    q_recovered = dcm2quat(R)

    # Sanity check on recovered quaternion
    assert abs(np.linalg.norm(q_recovered) - 1.0) < 1e-10, "Recovered quaternion must be unit"

    # Round-trip check: q and -q represent the same rotation, so allow sign flip
    same     = np.allclose(q_orig, q_recovered,  atol=1e-10)
    flipped  = np.allclose(q_orig, -q_recovered, atol=1e-10)
    assert same or flipped, (
        f"Round-trip failed!\n  Original:  {q_orig}\n  Recovered: {q_recovered}"
    )
    print("✓ Round-trip passed")
    print(f"  q_orig     = {q_orig}")
    print(f"  q_recovered= {q_recovered}")
    print(f"  Sign-flip  = {flipped}")

test_quat_dcm_roundtrip()
