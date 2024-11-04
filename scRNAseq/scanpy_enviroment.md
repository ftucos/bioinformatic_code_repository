# Scanpy enviroment

```bash
conda create -n scanpy python=3.7
conda activate scanpy
```

Install packages

```bash
conda install -y numpy scipy pandas matplotlib seaborn pytables scikit-learn statsmodels numba pytables notebook
conda install -y -c conda-forge python-igraph leidenalg louvain loompy h5py jupyterlab multicore-tsne ipywidgets jupyterlab_widgets
conda install -y -c conda-forge -c bioconda cellrank
conda install -y -c anaconda cytoolz ipykernel jupyter

pip install 'scanpy[leiden]' openpyxl fa2 pybind11 hnswlib pyscenic harmonypy bbknn phate wishbone_dev scikit-misc scrublet PhenoGraph
pip install git+https://github.com/theislab/scvelo@develop
pip install git+https://github.com/theislab/diffxpy
pip3 install git+https://github.com/jacoblevine/phenograph.git
pip install git+https://github.com/saezlab/decoupler-py
```

Before installing rpy2 be sure that you don't have other r installations in your environment with `where R`.

Type `R RHOME` and verify it's pointing to the system install of R (`/Library/Frameworks/R.framework/Resources`)

It's not recomended to rely on conda's installations of R since not all packages are available/compatible, and you will often encounter problems installing GitHub packages or packages that require compilation.

```bash
# install rpy2 with pip and not conda otherwise it will install r as a dependency
pip install rpy2
```

# sciVI tool enviroment

```python
conda create -n scvi-env python=3.9
conda activate scvi-env


conda install -y pytorch torchvision torchaudio -c pytorch
conda install -y jax jaxlib -c conda-forge

conda install -y -c conda-forge -c bioconda -c anaconda numpy scipy pandas matplotlib seaborn pytables scikit-learn statsmodels numba pytables notebook python-igraph leidenalg louvain loompy h5py jupyterlab multicore-tsne ipywidgets jupyterlab_widgets cytoolz ipykernel jupyter

conda install -y -c conda-forge scanpy python-igraph leidenalg

pip install 'scanpy[leiden]' openpyxl fa2 pybind11 hnswlib pyscenic harmonypy bbknn phate wishbone_dev scikit-misc scrublet PhenoGraph
pip install git+https://github.com/theislab/scvelo@develop
# pip install git+https://github.com/theislab/diffxpy
# pip3 install git+https://github.com/jacoblevine/phenograph.git
# pip install git+https://github.com/saezlab/decoupler-py

conda install -c conda-forge -c bioconda gseapy

conda install -c conda-forge scikit-misc

pip install rpy2

pip install Cython --install-option="--no-cython-compile"
conda install -c bioconda forceatlas2-python
```

conda install -y -c conda-forge -c bioconda -c anaconda numpy scipy pandas matplotlib seaborn pytables scikit-learn statsmodels numba pytables notebook python-igraph leidenalg louvain loompy h5py jupyterlab multicore-tsne ipywidgets jupyterlab_widgets cytoolz ipykernel jupyter

