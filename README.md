# GRIP-Tomo: GRaph Identification of Proteins in Tomograms

[![Tests](https://github.com/EMSL-Computing/grip-tomo/actions/workflows/python_build.yml/badge.svg)](https://github.com/EMSL-Computing/grip-tomo/actions/workflows/python_build.yml)
[![Documentation Status](https://readthedocs.org/projects/grip-tomo/badge/?version=latest)](https://grip-tomo.readthedocs.io/en/latest/?badge=latest)


A pipeline to help identify and classify structures in tomography data using global and local topological (network) features. 

The package consists of 3 core modules: 
1. `pdb2graph.py` - converts a PDB structure into a graph network. 
2. `density2graph.py` - converts a 3D density to a graph. 
3. `graph2class.py` - measures graph features and classifies

[Documentation](https://grip-tomo.readthedocs.io/en/latest/)

---

### Quickstart guide:

0. **open terminal / command prompt**
1. **clone the repository:** `git clone https://github.com/EMSL-Computing/grip-tomo`
2. **install the dependencies:** `pip install -r /path/to/requirements.txt`  
3. **run tests:** `cd path/to/grip_tomo/` then `python _test_basic.py`
4. **review example notebook, `grip_tomo/example_notebook.ipynb`**
    - To interact with `.ipynb` files install Jupyter-lab

---

### Citation
Please cite this work if you use it:

> George, A, Kim, DN, Moser, T, Gildea, IT, Evans, JE, Cheung, MS. Graph identification of proteins in tomograms (GRIP-Tomo). Protein Science. 2023; 32( 1):e4538. https://doi.org/10.1002/pro.4538




---

See the included license and disclaimer files

August George, PNNL, 2022
