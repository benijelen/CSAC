CSAC
====

Coordination Sphere Analysis and Comparison (CSAC) is a Python workflow for analyzing
metal-binding environments in protein structures. This repository snapshot includes the
core CSAC extraction scripts together with the downstream analysis and figure-generation
scripts used for the Fe-binding aerobic-versus-anaerobic manuscript.

Repository contents
-------------------
- `CSAC.py`: main extraction workflow for identifying target metal centers and local
  coordination environments from PDB structures.
- `CS_Freq_Vector.py`: converts CSAC database output into normalized amino-acid
  composition vectors.
- `stats_v2.py`: per-protein and mixed-effects statistical analyses used in the manuscript.
- `compositional_analysis_v2.py`: compositional amino-acid analyses used in the manuscript.
- `plot_figure4_ligand_distribution.py`
- `plot_figure4_ligand_distribution_fe_per_protein.py`
- `plot_figure5_ligand_size_violin.py`
- `plot_figure6_hydropathy_by_oxygentol.py`
- `plot_figure6_hydropathy_by_oxygentol_per_protein.py`
- `plot_figure7ab_pca_scree.py`
- `data/CSAC_sample_input.csv`: manuscript-aligned sample input file.
- `data/Supplementary_Table_1.xlsx`: finalized supplementary dataset used for submission.

Python requirements
-------------------
CSAC requires Python 3.9+ and the packages listed in `requirements.txt`.

Install dependencies with:

```bash
pip install -r requirements.txt
```

Basic CSAC usage
----------------
Run the core extraction workflow from a working directory that contains your input CSV and
any local `.pdb` structure files you want CSAC to use.

```bash
python CSAC.py
```

CSAC will prompt for:
- the input CSV filename
- the target metal element code(s), in uppercase (for example `FE` or `FE,CU`)

Core outputs
------------
- `*_DATABASE.csv`: per-center coordination-environment output
- `*_FREQ_VECTOR.csv`: normalized amino-acid composition vectors derived from the database

The downstream project-specific analysis scripts in this repository were used to generate
the summary statistics and figures reported for the Fe-binding manuscript.

Sample input format
-------------------
The first column should be the PDB identifier. Additional metadata columns may be included
to support downstream grouping and analysis. See `data/CSAC_sample_input.csv` for the
sample manuscript-aligned input file.

```csv
PDB ID,Metabolism,OxygenTOL
1APX,ROS Detox,aerobic
1B4U,Aromatic Dioxygenation,aerobic
1CNO,Nitrification,aerobic
```


Ni/Fe manuscript data artifacts
------------------------------
The Fe/Ni manuscript-supporting CSAC input, CSAC database output, regenerated amino-acid
frequency vectors, PDB-level descriptor table, and shared-family/global comparison summary
tables are provided in `data/ni_fe_manuscript/`.

Key files include:
- `csac_input_fe_ni_combined_local181.csv`: finalized combined Fe/Ni CSAC input set.
- `CSAC_fe_ni_combined_local_2026-03-12_FE_NI_DATABASE.csv`: final per-center CSAC output database.
- `CSAC_fe_ni_combined_local_2026-03-12_FE_NI_FREQ_VECTOR.csv`: regenerated amino-acid frequency vectors.
- `CSAC_fe_ni_combined_local_2026-03-12_pdb_level_descriptors.csv`: PDB-metal-level hydropathy and amino-acid descriptor table.
- `stats_track_ab_2026-03-12_compact.csv`: compact shared-family and global statistical summary.

License
-------
This project is distributed under the MIT License. See `LICENSE` for details.

Citing CSAC
-----------
If you use CSAC, please cite the software archive:

Jelen BI, Moore EK, Christensen B. CSAC (v1.1.4) [Computer software]. Zenodo.
https://doi.org/10.5281/zenodo.19319074

Repository
----------
GitHub: https://github.com/benijelen/CSAC
