.. grip-tomo documentation master file, created by
   sphinx-quickstart on Wed Oct 19 22:27:24 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to grip-tomo's documentation!
=====================================

GRaph Identification of Proteins in TOMOgrams (grip-tomo).

The package consists of 3 core modules: 

* `pdb2graph.py` - converts a PDB structure into a graph network. 

* `density2graph.py` - converts a 3D density to a graph. 

* `graph2class.py` - measures graph features and classifies


For more example usage see the included tutorial nobook.  

If you found this tool helpful or want more information please cite: "Graph identification of proteins in tomograms (GRIP-Tomo)" by August George, et al. Preprint DOI: https://doi.org/10.48550/arXiv.2210.08194 


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   pdb2graph
   density2graph
   graph2class


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
