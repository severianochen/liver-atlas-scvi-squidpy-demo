from pathlib import Path
import yaml


def load_yaml(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def ensure_dirs(config):
    for value in config["paths"].values():
        Path(value).mkdir(parents=True, exist_ok=True)


def find_files(root, suffixes):
    root = Path(root)
    out = []
    for suffix in suffixes:
        out.extend(root.rglob(f"*{suffix}"))
    return sorted(set(out))


def print_adata_summary(adata, name="adata"):
    print(f"\n{name}")
    print("=" * len(name))
    print(adata)
    print("obs columns:", list(adata.obs.columns)[:40])
    print("var columns:", list(adata.var.columns)[:40])
    print("obsm keys:", list(adata.obsm.keys()))
    print("layers:", list(adata.layers.keys()))
