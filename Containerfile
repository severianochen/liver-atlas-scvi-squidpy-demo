FROM docker.io/mambaorg/micromamba:2.0.5

COPY environment.yml /tmp/environment.yml
COPY requirements-pip.txt /tmp/requirements-pip.txt

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

RUN micromamba run -n base python -m pip install --upgrade pip && \
    micromamba run -n base python -m pip install --no-cache-dir torch==2.7.0 --index-url https://download.pytorch.org/whl/cu128 && \
    micromamba run -n base python -m pip install --no-cache-dir -r /tmp/requirements-pip.txt && \
    micromamba run -n base python -m pip uninstall -y torchvision torchaudio || true && \
    micromamba run -n base python -m pip install --no-cache-dir --force-reinstall --no-deps torch==2.7.0 --index-url https://download.pytorch.org/whl/cu128

WORKDIR /home/severiano/projects/liver-atlas-scvi-squidpy-demo

CMD ["bash"]
