# /home/karkodan/agri-stress-index/tests/test_leaf_temp.py

from models.leaf_temp import compute_leaf_temperature

inputs = {
    "pa": 101.3,
    "gh": 0.02,
    "Tair": 30.0,
    "D": 2.0,
    "gv": 0.03,
    "es": 4.2,
    "Qshort": 600,
    "csi": 1.5,
    "cloud_frac": 0.3
}

Tl = compute_leaf_temperature(**inputs)
print(f"🌿 Leaf temperature: {Tl:.2f} °C")
