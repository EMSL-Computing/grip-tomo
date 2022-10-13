# GRIP-Tomo: GRaph Identification of Proteins in Tomograms
A pipeline to help identify and classify structures in tomography data using global and local topological (network) features. 

[under development]

**requires** `python 3.7 or 3.8`

**API documentation** in `/docs`

### Quickstart guide:

0. **open terminal / command prompt**
1. **clone the repository:** `git clone https://github.com/EMSL-Computing/grip-tomo`
2. **install the dependencies:** `pip install -r /path/to/requirements.txt`  
3. **run tests:** `cd path/to/scripts/` then `python test_basic.py`
4. **review example notebook, `scripts/example_notebook.ipynb`**
    - To interact with `.ipynb` files, please install Jupyter notebook
---

### Overview

Please see example notebook and API for usage

The pipeline consists of 3 core modules: 
1. `pdb2graph.py` - converts a PDB structure into a graph network. Can be used to make a 'control' graph for classification.  
2. `density2graph.py` - converts a 3D density to a graph. 
3. `graph2class.py` - measures graph features and classifies

---

### Manuscript and Citation
---
**add citation here**

---

Please see included license and disclaimer files

August George, PNNL, 2022
