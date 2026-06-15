# GOTACH — Ground Truth & Catch-log

Track issues found during testing, fixes applied, and remaining gaps.

## Round 1: Environment + API Discovery (2026-06-15)

| # | Stage | File | Issue | Severity | Status | Fix |
|---|-------|------|-------|----------|--------|-----|
| 1 | env | spateo conda env | spateo not installed — empty env | CRITICAL | FIXED | Use spateo_auto env instead |
| 2 | env | spateo_auto | scipy ABI mismatch | CRITICAL | FIXED | `conda install "scipy>=1.9,<1.11"` |
| 3 | env | spateo_auto | scanpy/anndata import crash | CRITICAL | FIXED | `conda install "scanpy>=1.9" "anndata>=0.8,<0.10"` |
| 4 | env | spateo_auto | POT missing | HIGH | FIXED | `conda install -c conda-forge pot` |
| 5 | env | spateo 0.0.0 | `spateo.io.protocol_io` does NOT exist | CRITICAL | FIXED | Rewrote IO to use `st.io.read_*` |
| 6 | env | spateo 0.0.0 | `spateo.preprocessing.protocol_pipeline` missing | CRITICAL | FIXED | Use `st.pp.normalize_total/scale/log1p` |
| 7 | env | spateo 0.0.0 | `st.tools.morphofield_sparsevfc` missing | HIGH | FIXED | Found in `st.tdr` instead |
| 8 | env | spateo 0.0.0 | `st.tools.morphopath` missing | HIGH | FIXED | Found in `st.tdr` instead |
| 9 | IO | manual_read.py | Uses non-existent `spateo.io.protocol_io` | CRITICAL | FIXED | Rewrite to `st.io.read_*` |
| 10 | IO | auto_read.py | Uses non-existent `spateo.io.protocol_io.spatial.auto` | CRITICAL | FIXED | Implement dir-structure detection |
| 11 | alignment | two_slice.py | Uses non-existent `preprocess_spatial` | CRITICAL | FIXED | Use `st.pp.*` chain |
| 12 | alignment | multi_slice.py | Same as #11 | CRITICAL | FIXED | Same as #11 |
| 13 | morphogenesis | vectorfield.py | Uses non-existent `morphofield_sparsevfc` | HIGH | FIXED | Uses `st.tdr.*` |
| 14 | morphogenesis | feature.py | Depends on vector field output | HIGH | OPEN | Blocked by #13 |

## Round 2: Code Fixes (2026-06-15)

| # | Stage | File | Issue | Severity | Status | Fix |
|---|-------|------|-------|----------|--------|-----|
| 15 | IO | manual_read.py | Rewrote to use `st.io.read_*` directly | CRITICAL | FIXED | Manual CSV readers for MERFISH/seqFISH/STARmap+ |
| 16 | IO | auto_read.py | Rewrote with directory-structure detection | CRITICAL | FIXED | PLATFORM_SIGNATURES dict |
| 17 | alignment | two_slice.py | Rewrote preprocess_slices to use st.pp chain | CRITICAL | FIXED | normalize_total → log1p → select_hvf_seurat |
| 18 | morphogenesis | vectorfield.py | Uses `st.tdr.*` — all exist | INFO | OK | No fix needed |

## Round 3: Import Path Fixes (2026-06-15)

| # | Stage | File | Issue | Severity | Status | Fix |
|---|-------|------|-------|----------|--------|-----|
| 19 | all | *.py | `from ..shared.data_helpers` → `skills.shared` doesn't exist, directory is `00_shared` | CRITICAL | FIXED | Renamed `skills/00_shared/` → `skills/shared/` (numeric prefix invalid in Python import syntax) |
| 20 | all | *.py | `from ..00_shared.data_helpers` is syntax error (numeric module name) | CRITICAL | FIXED | Renamed directory to `shared`, use `from skills.shared.data_helpers` |
| 21 | pipeline | run_all.py | ReconstructionConfig has no `model_type` param | HIGH | FIXED | Use `models=[config.model_type]` list |
| 22 | pipeline | run_all.py | ReconstructionResult uses `output_dir` not `outdir` | HIGH | FIXED | Changed `.outdir` → `.output_dir` |
| 23 | pipeline | run_all.py | VectorFieldConfig has no `groupby` param | HIGH | FIXED | Removed `groupby`, pass `output_path` |
| 24 | pipeline | run_all.py | VectorFieldResult has no `outdir` | HIGH | FIXED | Return `vf_out` path directly from config |

## Round 4: Environment Dependencies (2026-06-15)

| # | Stage | File | Issue | Severity | Status | Fix |
|---|-------|------|-------|----------|--------|-----|
| 25 | env | matplotlib | CXXABI_1.3.15 not found | CRITICAL | FIXED | `conda install -c conda-forge libstdcxx-ng` |
| 26 | env | llvmlite | Same CXXABI issue | CRITICAL | FIXED | Same fix as #25 |
| 27 | env | numpy | pip pulled 2.2.6, scipy needs <1.27 | HIGH | FIXED | `pip install --force-reinstall "numpy==1.26.4"` |
| 28 | env | pynndescent | Numba dispatch error | HIGH | FIXED | `pip install "numba==0.58.1" "llvmlite==0.41.1" "pynndescent==0.5.13"` |
| 29 | reconstruction | spateo | PyMCubes missing | MEDIUM | FIXED | `pip install PyMCubes` |
| 30 | reconstruction | spateo | pymeshfix missing | MEDIUM | FIXED | `pip install pymeshfix==0.17.0` |
| 31 | reconstruction | spateo | pyacvd missing | MEDIUM | FIXED | `pip install pyacvd` |
| 32 | reconstruction | surface | Segfault in surface mesh (pymeshfix native) | HIGH | OPEN | Native library compat issue; point-cloud works |

## Round 5: Spateo Internal Bugs (2026-06-15)

| # | Stage | File | Issue | Severity | Status | Fix |
|---|-------|------|-------|----------|--------|-----|
| 33 | morphogenesis | spateo/alignment/paste.py | `paste_pairwise_align` passes `layer=` to `align_preprocess` which expects `rep_layer=` | CRITICAL | FIXED | Monkey-patch in vectorfield.py `patch_align_preprocess_kwarg()` + direct paste.py edit |
| 34 | morphogenesis | spateo/alignment/paste.py | `paste_pairwise_align` passes `select_high_exp_genes=`, `normalize_c=`, `normalize_g=` to `align_preprocess` which don't accept these | HIGH | FIXED | Patch drops these kwargs |
| 35 | morphogenesis | spateo/alignment/deprecated_utils.py | `align_preprocess` returns 8 values, but paste.py expects 7 | CRITICAL | FIXED | Updated unpacking in paste.py to accept 8 values |
| 36 | morphogenesis | spateo/alignment/paste.py | `exp_matrices[i]` is a list of arrays, not an array | HIGH | FIXED | Changed `exp_matrices[i]` → `exp_matrices[i][0]` |
| 37 | morphogenesis | spateo/alignment/paste.py | Variable renamed to `_nx` breaks later `nx.` usage | HIGH | FIXED | Kept variable as `nx` |

## Confirmed Working APIs (spateo 0.0.0 in spateo_auto env)

**IO:** `st.io.read_10x`, `st.io.read_bgi`, `st.io.read_bgi_agg`, `st.io.read_slideseq`, `st.io.read_nanostring`, `st.io.read_image`, `st.read(.h5ad)`
**PP:** `st.pp.normalize_total`, `st.pp.log1p`, `st.pp.scale`, `st.pp.select_hvf_seurat`, `st.pp.filter`
**Align:** `st.align.morpho_align`, `st.align.morpho_align_ref`, `st.align.group_pca`, `st.align.morpho_align_transformation`, `st.align.split_slice`, `st.align.paste_align`, `st.align.paste_align_ref`
**TDR (3D dev tools):** `st.tdr.morphofield_sparsevfc`, `st.tdr.morphopath`, `st.tdr.morphofield_velocity`, `st.tdr.morphofield_acceleration`, `st.tdr.morphofield_curvature`, `st.tdr.morphofield_curl`, `st.tdr.morphofield_torsion`, `st.tdr.morphofield_jacobian`, `st.tdr.glm_degs`, `st.tdr.pc_KDE`, `st.tdr.model_morphology`, `st.tdr.read_model`, `st.tdr.save_model`, `st.tdr.cell_directions`, `st.tdr.construct_surface`, `st.tdr.construct_point_cloud`

## Test Results Summary

| Stage | Test | Result | Notes |
|-------|------|--------|-------|
| 0: IO | MERFISH read | PASS | 78329 cells, 649 genes |
| 0: IO | STARmap+ read | PASS | 1022 cells, 91246 genes |
| 0: IO | Open-ST read | PASS | 49048 cells, 21609 genes |
| 0: IO | Auto-detect | PASS | Detects xenium, visium, merfish, starmap+ |
| 1: Alignment | Two-slice alignment | PASS | 43205 + 49927 cells, 86s on CPU |
| 2: Reconstruction | Point-cloud model | PASS | VTK output generated |
| 2: Reconstruction | Surface mesh | FAIL | pymeshfix native segfault (env issue, not code) |
| 3: Morphogenesis | Cell mapping + vector field | PASS | 4 models saved, 21s on CPU |

## Remaining Action Items

1. **[HIGH]** Fix surface mesh segfault — pymeshfix/pyvista native library incompatibility with system libstdc++
2. **[MEDIUM]** Add PyMCubes, pymeshfix, pyacvd to environment.yml/requirements.txt
3. **[MEDIUM]** Document spateo internal bugs and patches in SKILL.md for morphogenesis stage
4. **[LOW]** Run multi-slice alignment test (multi_slice.py)
5. **[LOW]** Run two_stage_3d.py test
6. **[LOW]** Run Stage 3 with compute_trajectory=True
