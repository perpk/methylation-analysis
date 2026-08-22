# This is to setup a conda environment based on rapidsai with Jupyter

- Get Miniconda for Linux x86(64)
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3-latest-Linux-x86_64.sh
```
- Install Miniconda
```
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```
- Create a conda environment with rapids and python 3.13
```
conda create -n rapids-26.04 -c rapidsai -c conda-forge rapids=26.04 python=3.13 'cuda-version>=12.0,<=12.6'

```
- Activate the environment
```
conda activate rapids-26.04
```
- Install GPU-enabled PyTorch

  RAPIDS environments mix the `rapidsai` and `conda-forge` channels, and conda's
  dependency solver will often resolve `pytorch`/`pytorch-cuda` against a
  **CPU-only** build (e.g. `cpu_mkl_*`) in that context instead of the CUDA
  build, even when `pytorch-cuda` is explicitly requested. To reliably get a
  CUDA-enabled build, install PyTorch via `pip` using the official CUDA wheel
  index instead of conda. Pick the `cuXXX` index matching the `cuda-version`
  used when creating the environment (e.g. `cu124` for CUDA 12.4):
  ```
  pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu124
  ```
  Verify the install actually sees the GPU:
  ```python
  import torch
  print(torch.__version__, torch.cuda.is_available(), torch.cuda.get_device_name(0))
  ```
  If a CPU-only `pytorch`/`libtorch` was previously installed via conda, remove
  it first so it doesn't shadow/conflict with the pip-installed GPU build:
  ```
  conda remove pytorch libtorch -y
  ```
- Install Jupyter 
```
conda install conda-forge::jupyter
```

- Change/Set the password
```
jupyter notebook password
```

- Bootstrap a server listening on port 8080 (Needs port forwarding via an ssh-session)
```
jupyter notebook --ip=0.0.0.0 --port=8080 --allow-root --no-browser
```
- Install a realtime monitoring tool for GPU usage visualization e.g. nvitop
  
```
conda install -c conda-forge nvitop
```

## Docker-based setup (docker-ml-cuda)

The [docker-ml-cuda/Dockerfile](docker-ml-cuda/Dockerfile) automates the steps
above for a fresh container:

1. Starts from `nvidia/cuda:12.4.1-devel-ubuntu22.04`.
2. Installs Miniconda and creates the `rapids-env` conda environment
   (`rapids=26.02`, `python=3.13`, `cuda-version=12.4`, `jupyterlab`).
3. Registers `rapids-env` as a Jupyter kernel.
4. **Installs GPU-enabled PyTorch via `pip`** (`--index-url
   https://download.pytorch.org/whl/cu124`), matching the `cuda-version=12.4`
   used to create the environment. This step exists specifically because
   `conda install pytorch pytorch-cuda -c pytorch -c nvidia` resolves to a
   CPU-only build inside `rapids-env` (see the "Install GPU-enabled PyTorch"
   note above) — pip with the official CUDA wheel index is the reliable way
   to get a CUDA build here.
5. Sets up SSH access and a startup script that runs `sshd` plus JupyterLab
   under `rapids-env` on port `8080`.

To verify PyTorch has GPU access after building/running the container:
```python
import torch
print(torch.__version__, torch.cuda.is_available(), torch.cuda.get_device_name(0))
```

If you change the CUDA version in the Dockerfile (base image, `cuda-version`
in the `conda create` step), update the `cuXXX` suffix in the pip
`--index-url` for PyTorch to match.
