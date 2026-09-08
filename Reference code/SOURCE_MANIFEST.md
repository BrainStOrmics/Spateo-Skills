# Working-notebook consolidation manifest

The source paths below are relative to `/DATA/User/gaomohan/figures`. Exact duplicates, incomplete scratch notebooks, and historical compatibility versions were consolidated instead of being published repeatedly.

The four historical alignment drafts are relative to `/DATA/User/gaomohan`, and the two supplementary developmental notebooks originate from `/DATA/User/gaomohan/figures2_补充`. Their published paths use the English source alias `figures2_supplement` in notebook provenance.

| Working notebook | Curated destination or disposition |
|---|---|
| `Alignment/Basic Alignment.ipynb` | `00_alignment_methods/01_basic_serial_slice_alignment.ipynb` |
| `Alignment/Nonrigid Alignment.ipynb` | `00_alignment_methods/02_nonrigid_serial_slice_alignment.ipynb` |
| `Alignment/cs13_Alignment.ipynb` | `00_alignment_methods/03_cs13_alignment_parameter_comparison.ipynb` |
| `Alignment/drosophila_Alignment.ipynb` | `00_alignment_methods/04_drosophila_serial_section_alignment.ipynb` |
| `figures2_supplement/7dpa_10dpa_Align_models.ipynb` | `04_alignment_and_morphogenesis/01_7dpa_10dpa_model_alignment.ipynb` |
| `figures2_supplement/morphogenesis_7dpa_10dpa_CNS.ipynb` | `04_alignment_and_morphogenesis/02_7dpa_10dpa_cns_morphogenesis.ipynb` |
| `figure2/a/multi_tech_align.ipynb` | Empty import/data-loading scaffold; omitted |
| `figure2/c/10dpa_14dpa_Align_models.ipynb` | Superseded by `04_alignment_and_morphogenesis/03_10dpa_14dpa_model_alignment.ipynb` |
| `figure2/c/10dpa_14dpa_Align_models.spateo_current.ipynb` | `04_alignment_and_morphogenesis/03_10dpa_14dpa_model_alignment.ipynb` |
| `figure2/c/10dpa_14dpa_lineage_analysis.ipynb` | Optional lineage-label QC merged into the model-alignment notebook |
| `figure2/c/ROD_morph_backbone.ipynb` | `03_axis_and_backbone/01_rod_backbone.ipynb` |
| `figure2/d/ROD_morph_backbone_with_continuous_heatmap.ipynb` | `03_axis_and_backbone/02_rod_backbone_gene_heatmap.ipynb` |
| `figure2/d/ROD_morph_backbone_with_heatmap.ipynb` | Superseded by the continuous-heatmap notebook |
| `figure2/f/morphogenesis_10dpa_14dpa_CNS.ipynb` | Superseded by `04_alignment_and_morphogenesis/04_10dpa_14dpa_cns_morphogenesis.ipynb` |
| `figure2/f/morphogenesis_10dpa_14dpa_CNS.spateo_current.ipynb` | `04_alignment_and_morphogenesis/04_10dpa_14dpa_cns_morphogenesis.ipynb` |
| `figure3/b/cs13_models_downsampling.ipynb` | `01_three_dimensional_reconstruction/04_cs13_downsampling_point_clouds.ipynb` |
| `figure3/b/cs13_models_reconstruction(trn100k).ipynb` | `01_three_dimensional_reconstruction/02_cs13_trn100k_reconstruction.ipynb` |
| `figure3/b/downspamling_mesh_similarity.ipynb` | `02_morphological_features/03_cs13_downsampling_mesh_similarity.ipynb` |
| `figure3/b/mesh_similarity.ipynb` | Incomplete scratch cells merged into the CS13 downsampling mesh-similarity notebook |
| `figure3/b_test/cs13_models_reconstruction(trn100k).ipynb` | Parameter trials consolidated into the final TRN100k reconstruction profile |
| `figure3/b_test/cs13_models_reconstructuon.ipynb` | `01_three_dimensional_reconstruction/03_cs13_point_cloud_mesh_voxel.ipynb` |
| `figure3/b_test/cs13_trn100k_3d_v1/scripts_snapshot/cs13_trn100k_celltype_marker_plotting.ipynb` | Exact duplicate; see `01_three_dimensional_reconstruction/06_cs13_celltype_marker_plotting.ipynb` |
| `figure3/b_test/cs13_trn100k_3d_v1/scripts_snapshot/cs13_trn100k_reference_mesh_reconstruction.ipynb` | Exact duplicate; see `01_three_dimensional_reconstruction/01_cs13_reference_mesh_reconstruction.ipynb` |
| `figure3/b_test/cs13_trn100k_celltype_marker_plotting.ipynb` | `01_three_dimensional_reconstruction/06_cs13_celltype_marker_plotting.ipynb` |
| `figure3/b_test/cs13_trn100k_reference_mesh_reconstruction.ipynb` | `01_three_dimensional_reconstruction/01_cs13_reference_mesh_reconstruction.ipynb` |
| `figure3/b_test/downsampling.ipynb` | Incomplete exploratory subset merged into the CS13 downsampling notebooks |
| `figure3/c/cs13_celltype.ipynb` | `01_three_dimensional_reconstruction/05_cs13_celltype_visualization.ipynb` |
| `figure3/d/3D_Models_Reconstruction(myotome).ipynb` | `01_three_dimensional_reconstruction/07_myotome_bilateral_reconstruction.ipynb` |
| `figure3/d/Myotome_backbone_gene_analysis_20260821/notebooks/backbone(myotome)_left_with_continuous_genes.ipynb` | `03_axis_and_backbone/03_myotome_left_backbone_genes.ipynb` |
| `figure3/d/Myotome_backbone_gene_analysis_20260821/notebooks/backbone(myotome)_right_with_continuous_genes.ipynb` | `03_axis_and_backbone/04_myotome_right_backbone_genes.ipynb` |
| `figure3/e/cell_density.ipynb` | `02_morphological_features/02_tissue_cell_density.ipynb` |
| `figure3/f/backbone(myotome)_left.ipynb` | Superseded by the continuous-gene left-myotome notebook |
| `figure3/f/backbone(myotome)_right.ipynb` | Superseded by the continuous-gene right-myotome notebook |
| `figure3/f/mesh(downsampling).ipynb` | `02_morphological_features/04_myotome_downsampling_mesh_similarity.ipynb` |
| `figure3/f/mesh(myotome).ipynb` | `02_morphological_features/01_myotome_morphology.ipynb` |
| `figure3/g/3D_interpolation(投稿).ipynb` | `05_interpolation/01_myotome_3d_gene_expression_interpolation.ipynb` |
| `figure3/g/3D_interpolation.ipynb` | Neural-interpolation cells merged into the curated interpolation notebook |

## Curation changes applied throughout

- Removed all cached outputs and execution counts.
- Removed historical PyMeshFix and POT monkeypatches now handled by Spateo.
- Removed unused Dynamo imports and replaced required Dynamo calls with Spateo-native APIs.
- Translated markdown, comments, and user-facing messages to English.
- Normalized Python formatting and corrected obvious incomplete or duplicate cells.
- Added tutorial-stage headings and provenance metadata.
- Preserved dataset-specific scientific parameters unless a documented current-Spateo compatibility change was required.
- Replaced Dynamo preprocessing/PCA in historical alignment drafts with `st.pp` normalization, `log1p_layer`, PCA, or `group_pca` as appropriate.
- Removed the unrelated 14 dpa branch from the 7-to-10 dpa alignment notebook and corrected its output path.
- Harmonized shared nonzero genes before each cross-stage CNS cell-mapping workflow.
