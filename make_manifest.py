from pathlib import Path
import json

ATLAS_DIR = Path("atlas")  # your local atlas folder
BASE_URL = "https://earth.cs.mcgill.ca/sesame-explorer/atlas"

manifest = {}

for sphere_dir in sorted(p for p in ATLAS_DIR.iterdir() if p.is_dir()):
    sphere = sphere_dir.name
    files = sorted(sphere_dir.rglob("*.nc"))
    manifest[sphere] = [
        f"{BASE_URL}/{sphere}/{f.relative_to(sphere_dir).as_posix()}"
        for f in files
    ]

Path("atlas_manifest.json").write_text(
    json.dumps(manifest, indent=2)
)

print("atlas_manifest.json created")
