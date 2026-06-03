from pathlib import Path
import subprocess

URLS = {
    # Liver Cell Atlas portal
    "lca/rawData_human.zip": "https://www.livercellatlas.org/data_files/toDownload/rawData_human.zip",
    "lca/annot_humanAll.csv": "https://www.livercellatlas.org/data_files/toDownload/annot_humanAll.csv",
    "lca/rawData_humanVisium.zip": "https://www.livercellatlas.org/data_files/toDownload/rawData_humanVisium.zip",
    "lca/annot_humanVisium.csv": "https://www.livercellatlas.org/data_files/toDownload/annot_humanVisium.csv",
    "lca/Human_fullResImg.zip": "https://www.livercellatlas.org/data_files/toDownload/Human_fullResImg.zip",

    # Very good HCC/fetal liver/mouse liver fallback, includes ready h5ad files
    "GSE156625/GSE156625_HCCscanpyobj.h5ad.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156625/suppl/GSE156625_HCCscanpyobj.h5ad.gz",
    "GSE156625/GSE156625_HCCFscanpyobj.h5ad.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156625/suppl/GSE156625_HCCFscanpyobj.h5ad.gz",
    "GSE156625/GSE156625_mousescanpyobj.h5ad.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156625/suppl/GSE156625_mousescanpyobj.h5ad.gz",

    # HCC ecosystem, good cancer fallback
    "GSE149614/GSE149614_HCC.metadata.updated.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz",
    "GSE149614/GSE149614_HCC.scRNAseq.S71915.count.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz",

    # Normal human liver atlas, smaller
    "GSE115469/GSE115469_CellClusterType.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115469/suppl/GSE115469_CellClusterType.txt.gz",
    "GSE115469/GSE115469_Data.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115469/suppl/GSE115469_Data.csv.gz",

    # Aizarani normal liver atlas, small processed files
    "GSE124395/GSE124395_Normalhumanlivercellatlasdata.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395_Normalhumanlivercellatlasdata.txt.gz",
    "GSE124395/GSE124395_Normalhumanliverdata.RData.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395_Normalhumanliverdata.RData.gz",
    "GSE124395/GSE124395_clusterpartition.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395_clusterpartition.txt.gz",

    # Liver cancer T-cell dataset, small and fast
    "GSE98638/GSE98638_HCC.TCell.S5063.count.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98638/suppl/GSE98638_HCC.TCell.S5063.count.txt.gz",
    "GSE98638/GSE98638_HCC.TCell.S5063.TPM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98638/suppl/GSE98638_HCC.TCell.S5063.TPM.txt.gz",

    # Cirrhosis/fibrotic niche, useful if cancer downloads are slow
    "GSE136103/GSE136103_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136103/suppl/GSE136103_RAW.tar",
}

BASE = Path("data/raw/liver_datasets")
LOGS = Path("logs/downloads")
BASE.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)

for rel, url in URLS.items():
    out = BASE / rel
    out.parent.mkdir(parents=True, exist_ok=True)
    log = LOGS / (rel.replace("/", "__") + ".log")

    print(f"\n=== Downloading/resuming {rel} ===")
    cmd = [
        "wget", "-c",
        "--tries=0",
        "--timeout=30",
        "--read-timeout=30",
        "--waitretry=10",
        "-O", str(out),
        url,
    ]

    with open(log, "a") as handle:
        subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=False)

print("\nDone. Check:")
print("find data/raw/liver_datasets -type f -exec ls -lh {} \\;")
