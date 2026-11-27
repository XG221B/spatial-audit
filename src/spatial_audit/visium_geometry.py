"""Calibrate native Visium capture-array addresses against deposited image coordinates."""

import numpy as np

V1_LAYOUT: dict = {
    "name": "V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz",
    "url": "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz",
    "bytes": 9960960,
    "sha256": "415530a9d31172f63fd358d5d1cefb0798698ef9c0e614f52825894cd45c61bf",
}

PITCH = 100.0
PITCH_TOLERANCE = 0.1
IMAGE_DIAMETER_UM = 65.0
SCALING_SOURCE = (
    "https://www.10xgenomics.com/support/software/space-ranger/2.1/analysis/outputs/spatial-outputs"
)


def checked_visium_lattice(positions, pixels, diameter, vendor_positions=None):
    """Validate the nominal lattice pitch against the independent Space Ranger image calibration."""
    rc = np.asarray(positions, dtype=float)
    pixels = np.asarray(pixels, dtype=float)
    if rc.ndim != 2 or rc.shape[1] != 2 or not len(rc):
        raise ValueError("A nonempty two-column Visium grid is required")
    if not np.isfinite(rc).all() or not np.equal(rc, np.rint(rc)).all() or (rc < 0).any():
        raise ValueError("Invalid deposited array address")
    if (
        pixels.shape != rc.shape
        or not np.isfinite(pixels).all()
        or not np.isfinite(diameter)
        or diameter <= 0
    ):
        raise ValueError("Invalid archived image calibration")
    rc = rc.astype(int)
    if np.any((rc[:, 0] + rc[:, 1]) % 2) or len(np.unique(rc, axis=0)) != len(rc):
        raise ValueError("Unexpected duplicate or off-lattice spot")
    if vendor_positions is not None and not np.array_equal(rc, vendor_positions):
        raise ValueError("Archived array addresses differ from manufacturer barcode layout")
    locations = {tuple(row): i for i, row in enumerate(rc)}
    pairs = [
        (i, locations[key])
        for i, (r, c) in enumerate(rc)
        for key in ((r, c + 2), (r + 1, c + 1), (r + 1, c - 1))
        if key in locations
    ]
    if not pairs:
        raise ValueError("No neighboring spot pair available for independent pitch calibration")
    l, r = np.asarray(pairs).T
    # Space Ranger defines this display scale against a theoretical 65 um
    # diameter. The capture area inside the physical spot has diameter 55 um.
    um = np.linalg.norm(pixels[l] - pixels[r], axis=1) * IMAGE_DIAMETER_UM / diameter
    quantiles = np.quantile(um, [0.05, 0.5, 0.95])
    image_agrees = bool(np.max(abs(quantiles / PITCH - 1)) <= PITCH_TOLERANCE)
    if not image_agrees and vendor_positions is None:
        raise ValueError(f"Image and array pitch calibrations differ: {quantiles.tolist()}")
    xy = np.column_stack((PITCH * rc[:, 1] / 2, PITCH * rc[:, 0] * np.sqrt(3) / 2))
    return xy, {
        "adjacent_pairs": len(pairs),
        "image_pitch_within_tolerance": image_agrees,
        "vendor_addresses_checked": len(rc) if vendor_positions is not None else None,
        "image_pitch_p05_p50_p95_um": quantiles,
        "archived_spot_diameter_fullres_px": float(diameter),
        "spaceranger_theoretical_diameter_um": IMAGE_DIAMETER_UM,
        "spaceranger_scale_source": SCALING_SOURCE,
    }
