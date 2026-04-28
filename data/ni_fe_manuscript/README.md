# Ni/Fe Manuscript CSAC Artifacts

This directory contains the Fe/Ni manuscript-supporting CSAC input and output artifacts referenced by the methods and data-availability text.

## Core Input

- `csac_input_fe_ni_combined_local181.csv`: finalized combined local Fe/Ni CSAC input list. This input produced analyzable Fe and/or Ni output for 181 unique PDB structures.
- `csac_input_fe_ni_combined_local181_excluded_missing_fe_ni.csv`: audit list of PDB IDs removed because the combined CSAC run produced no Fe or Ni output rows.

## Core CSAC Outputs

- `CSAC_fe_ni_combined_local_2026-03-12_FE_NI_DATABASE.csv`: final CSAC per-center output database for the combined Fe/Ni run.
- `CSAC_fe_ni_combined_local_2026-03-12_FE_NI_FREQ_VECTOR.csv`: amino-acid frequency vectors regenerated from the final CSAC database.
- `CSAC_fe_ni_combined_local_2026-03-12_pdb_level_descriptors.csv`: PDB-metal-level descriptor table aggregating average hydropathy, hydropathy SD, amino-acid frequencies, cysteine fraction, nearby-residue count, and CSAC center counts.
- `CSAC_fe_ni_combined_local_2026-03-12_missing_metals.csv`: audit output from the multi-metal CSAC run. For Fe/Ni runs, missing requested metals are expected for Fe-only or Ni-only PDBs and are not automatically failures.

## Analysis Summaries

- `stats_track_ab_2026-03-12_compact.csv`: compact Track A and Track B statistical summary used for shared-family and global comparisons.
- `stats_track_a_family_2026-03-12.csv`: family-stratified shared-family statistics.
- `fig6_trackAB_cys_summary.csv`: Figure 6 cysteine summary table.
- `Supplementary_Table_SX_PDB_family_metal_map.csv`: PDB-family-metal mapping table.
- `Supplementary_Table_SY_family_level_counts.csv`: family-level count table.
- `Supplementary_Table_SZ_shared_families_pdb_list_fig6A.csv`: shared-family PDB list used for Figure 6A.

## Analysis Policy

The frozen analysis policy for these artifacts is:

- CSAC analysis unit: ligand centroid (`B` mode in the local project history)
- NMR entries retained
- first model only for multi-model structures
- PDB-level inference used for statistical comparisons to avoid overweighting proteins with many metal centers

The shared-family Track A comparison is restricted to ACS, CODH, and NiFe hydrogenase families. The global Track B comparison uses the full curated Fe/Ni dataset and is composition-sensitive.
