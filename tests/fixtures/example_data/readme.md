August George, PNNL, 2022

Example datasets to use for testing/validation/tutorial

I used the command line interface to generate the data from `pdb2graph.py` and `density2graph.py`. See the docs/api reference for usage examples.

---

PDB files: 
- `3vjf.pdb`, Alpha helix bundle
- `4rlc.pdb`, Beta sheet barrel
- `4rly.pdb`, Ankyrin Repeats
- `1prw.pdb`, Calmodulin collapsed conformation 
- `1cll.pdb`, Calmodulin extended conformaion 
---

3D density (.mrc) files:
- PDB to ideal density (.mrc) using `pdb2mrc_by_EMAN2.py` (by Doo Nam) 
    - alpha carbon only, 1 A/pixel, 3 A resolution (equivalent to Gaussian lowpass with 1/e width at 1/res), 400 pixel box size
        - `3vjf_CA_apix_1_box_dimension_400_res_3.mrc`
- Reconstructed tomograms using `generate_rawtlt.bash` and `reconstruct.bash` along with density projection .mrcs files made by Doo Nam
    - 180 degree tilt series (no missing wedge), 730A defocus. reconstruction files (.rec) renamed to .mrc file for processing
        - `3vjf_CA_apix_1_box_dimension_400_res_3_reconstructed.mrc`

---

Graph (.gexf) files:
- PDB to graph using `pdb2graph.py` - these are the 'control' graphs 
    - using alpha carbons as nodes, 8A pairwise cutoff distance for edges
        - `3vjf.gexf`
        - `4rlc.gexf`
        - `4rly.gexf`
        - `1prw.gexf`
        - `1cll.pdb`

- PDB to ideal density to graph using `pdb2mrc_by_EMAN2.py` and `density2graph.py`
    - pdb to ideal density: alpha carbon only, 1 A/pixel, 3 A resolution, 400 pixel box size
    - density to graph: variable threshold (see below), cluster radius = 4 pixels, min cluster size = 1, 8A pairwise cutoff distance
        - `3vjf_CA_apix_1_box_dimension_400_res_3.gexf`, unormalized pixel threshold = 0.425
- PDB to ideal density to synthetic tomogram density to graph using `pdb2mrc_by_EMAN2.py`, `generate_rawtlt.bash`,  `reconstruct.bash` and `density2graph.py`
    - pdb to ideal density: alpha carbon only, 1 A/pixel, 3 A resolution, 400 pixel box size
    - ideal density to tomogram: 180 degree tilt series, 730A defocus
    - tomogram density to graph: variable threshold (see below), cluster radius = 4 pixels, min cluster size = 1, 8A pairwise cutoff distance
        - `3vjf_CA_apix_1_box_dimension_400_res_3_reconstructed.gexf`, unormalized pixel threshold = 26985
