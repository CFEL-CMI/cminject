Recommended install commands:

#For python3.9+

python -m pip install -U pip setuptools wheel
# prefer binary wheels for heavy deps
python -m pip install --only-binary=:all: numpy scipy pandas numba h5py matplotlib
# install your package (PEP 517)
python -m pip install -e .
# docs if needed
python -m pip install .[docs]

#For Python 3.6.8 (legacy):

python -m pip install -U "pip<21.4" "setuptools<60" "wheel<0.38"
python -m pip install -U "Cython>=0.29.36,<0.30" "numpy==1.19.5"
python -m pip install --only-binary=:all: "scipy==1.5.4" "pandas==1.0.5" "numba==0.53.1" "matplotlib==3.3.4" "h5py==3.1.0"
python -m pip install -e . --no-build-isolation

