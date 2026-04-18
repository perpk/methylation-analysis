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
conda create -n rapids-26.04 -c rapidsai -c conda-forge rapids=26.04 python=3.13 'cuda-version>=13.0,<=13.1'
```
- Activate the environment
```
conda activate rapids-26.04
```
- Install Jupyter 
```
conda install conda-forge::jupyter
```
- Bootstrap a server listening on port 8080 (Needs port forwarding via an ssh-session)
```
jupyter notebook --ip=0.0.0.0 --port=8080 --allow-root --no-browser
```