from pathlib import Path
import zipfile
import tarfile
import gzip

roots = [Path("data/raw/liver_cell_atlas")] + sorted(Path("data/raw").glob("GSE*"))

for root in roots:
    if not root.exists():
        continue

    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue

        print("\n" + "=" * 80)
        print(path)
        print(f"Size: {path.stat().st_size / 1024 / 1024:.1f} MB")

        if path.suffix == ".zip":
            print("ZIP:", "valid" if zipfile.is_zipfile(path) else "incomplete/corrupt")
        elif path.suffix == ".tar":
            print("TAR:", "valid" if tarfile.is_tarfile(path) else "incomplete/corrupt")
        elif path.suffix == ".gz":
            try:
                with gzip.open(path, "rb") as f:
                    f.read(1024)
                print("GZ: readable")
            except Exception:
                print("GZ: incomplete/corrupt")
        else:
            print("Plain file or unknown format")
