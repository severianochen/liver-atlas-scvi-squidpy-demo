IMAGE=liverdemo:0.1
PROJECT_ROOT=/home/severiano/projects/liver-atlas-scvi-squidpy-demo

CACHE_ENV=-e HOME=/tmp/home -e XDG_CACHE_HOME=/tmp/xdg-cache -e MPLCONFIGDIR=/tmp/matplotlib-cache -e NUMBA_CACHE_DIR=/tmp/numba-cache
BASE_RUN=podman run --rm --userns=keep-id --user "$$(id -u):$$(id -g)" $(CACHE_ENV) -v $(PROJECT_ROOT):$(PROJECT_ROOT):Z -w $(PROJECT_ROOT) $(IMAGE)
BASE_RUN_GPU=podman run --rm --gpus all --userns=keep-id --user "$$(id -u):$$(id -g)" $(CACHE_ENV) -v $(PROJECT_ROOT):$(PROJECT_ROOT):Z -w $(PROJECT_ROOT) $(IMAGE)
BASE_RUN_IT=podman run --rm -it --gpus all --userns=keep-id --user "$$(id -u):$$(id -g)" $(CACHE_ENV) -v $(PROJECT_ROOT):$(PROJECT_ROOT):Z -w $(PROJECT_ROOT) $(IMAGE)

PY=bash -c 'export PATH=/opt/conda/bin:$$PATH; mkdir -p "$$HOME" "$$XDG_CACHE_HOME" "$$MPLCONFIGDIR" "$$NUMBA_CACHE_DIR"; exec python "$$@"' _
SH=bash -c 'export PATH=/opt/conda/bin:$$PATH; mkdir -p "$$HOME" "$$XDG_CACHE_HOME" "$$MPLCONFIGDIR" "$$NUMBA_CACHE_DIR"; exec bash'

build:
	podman build --format=docker -t $(IMAGE) -f Containerfile .

build-clean:
	podman build --no-cache --format=docker -t $(IMAGE) -f Containerfile .

gpu-test:
	$(BASE_RUN_GPU) $(PY) -c "import torch; print(torch.__version__); print('cuda:', torch.cuda.is_available()); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU only')"

import-test:
	$(BASE_RUN_GPU) $(PY) -c "import scanpy as sc, squidpy as sq, scvi, anndata as ad, torch; print('scanpy:', sc.__version__); print('squidpy:', sq.__version__); print('scvi-tools:', scvi.__version__); print('anndata:', ad.__version__); print('torch:', torch.__version__); print('cuda:', torch.cuda.is_available()); print('imports ok')"

shell:
	$(BASE_RUN_IT) $(SH)

download:
	$(BASE_RUN) $(PY) src/00_download_lca.py

inspect:
	$(BASE_RUN) $(PY) src/00_inspect_archives.py

import-sc:
	$(BASE_RUN) $(PY) src/01_import_lca_single_cell.py

import-manual:
	$(BASE_RUN) $(PY) src/01_import_manual_h5ad.py

qc:
	$(BASE_RUN) $(PY) src/02_qc_preprocess.py

markers:
	$(BASE_RUN) $(PY) src/03_score_markers.py

scvi:
	$(BASE_RUN_GPU) $(PY) src/04_train_scvi.py

import-visium:
	$(BASE_RUN) $(PY) src/05_import_lca_visium.py

add-visium-coords:
	$(BASE_RUN) $(PY) src/05b_add_lca_visium_geo_coordinates.py

spatial:
	$(BASE_RUN) $(PY) src/06_spatial_squidpy.py

spatial-smoketest:
	$(BASE_RUN) $(PY) src/06_spatial_squidpy_smoketest.py

holdout:
	$(BASE_RUN) $(PY) src/07_holdout_single_cell.py

singlecell-all: import-sc qc markers scvi

spatial-all: import-visium add-visium-coords spatial

all: download inspect import-sc qc markers scvi import-visium add-visium-coords spatial
