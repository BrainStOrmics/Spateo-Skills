<div align="center">

<img src="_static/logo.png" width="100" alt="Spateo-Skills Logo" />

# 🧬 Spateo-Skills

**面向 Spateo 空间转录组分析流程的标准化、可复用 Skill 模块**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Spateo](https://img.shields.io/badge/spateo-0.0.0-orange.svg)](https://github.com/aristoteleo/spateo-release)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](../LICENSE)

📖 [文档](README.md) · 🌐 [English](../README.md) · 🧪 [测试](../tests/) · 🐛 [问题记录](../tests/GOTACH.md)

</div>

---

## 📋 简介

Spateo-Skills 将 Jupyter Notebook 中的分析流程转化为**可调用、参数化、可通过 CLI 运行**的 Python 模块。每个 Skill 遵循统一模式：

```
Config（类型化 dataclass） → run_pipeline(config) → Result（类型化 dataclass）
```

✨ **懒加载导入** — 导入 Skill 不会加载 Spateo  
✨ **类型化接口** — `@dataclass` Config/Result，参数均有默认值  
✨ **可组合** — 每阶段输出直接作为下一阶段输入  
✨ **CLI 可用** — 每个 Skill 均可通过 `python -m` 独立运行  
✨ **Agent 友好** — 接口可预测、无魔法、无副作用  

---

## 🚀 快速开始

### 1️⃣ 安装环境

```bash
bash skills/00_env_setup/setup.sh   # 自动检测 conda 或 uv
conda activate spateo_auto
```

### 2️⃣ 一键运行全流程

```python
from skills.00_pipeline.run_all import PipelineConfig, run_full_pipeline

result = run_full_pipeline(PipelineConfig(
    data_path="./data/xenium_outs/",
    platform="xenium",
    stage2_path="./data/stage2.h5ad",
    out_dir="./output/",
))
```

### 3️⃣ 分步调用 Skill

```python
import importlib

# 📖 读取数据 — 自动检测平台
reader = importlib.import_module("skills.01_data_io.auto_read")
adata, platform = reader.read_auto_spatial_data("./data/unknown/", return_match=True)
adata.write_h5ad("./data/stage0.h5ad")

# 🔗 对齐 — 两片组织切片
align = importlib.import_module("skills.02_slice_alignment.two_slice")
cfg = align.TwoSliceConfig(
    slice1_path="./data/slice_a.h5ad",
    slice2_path="./data/slice_b.h5ad",
    outdir="./output/alignment/",
)
r1 = align.run_two_slice_pipeline(cfg)

# 🏗️ 重建 — 3D 点云 + 表面网格
recon = importlib.import_module("skills.03_3d_reconstruction.reconstruction")
cfg = recon.ReconstructionConfig(
    input_path="./data/aligned.h5ad",
    out_dir="./output/models/",
    models=["point-cloud", "surface"],
    groupby="tissue_type",
)
r2 = recon.run_reconstruction_pipeline(cfg)
# r2.point_cloud_model → VTK 文件路径

# 🧭 形态发生 — 细胞迁移向量场
vf = importlib.import_module("skills.04_morphogenesis.vectorfield")
cfg = vf.VectorFieldConfig(
    stage1_path="./data/stage1.h5ad",
    stage2_path="./data/stage2.h5ad",
    output_path="./output/morph/",
)
r3 = vf.run_vectorfield_pipeline(cfg)
# r3.saved_model_paths → VTK 文件列表
```

---

## 📦 Skill 索引

| 阶段 | Skill | 模块 | 输入 | 输出 |
|---|---|---|---|---|
| 🔧 | **env_setup** | `skills.00_env_setup` | `setup.sh` | `spateo_auto` 环境 |
| 📖 | **auto_read** | `skills.01_data_io.auto_read` | 数据目录 | `AnnData` |
| 📖 | **manual_read** | `skills.01_data_io.manual_read` | 平台 + 路径 | `AnnData` |
| 🔗 | **two_slice** | `skills.02_slice_alignment.two_slice` | 两个 `.h5ad` | 对齐后 `AnnData` |
| 🔗 | **multi_slice** | `skills.02_slice_alignment.multi_slice` | N 个 `.h5ad` | 对齐后 `AnnData` |
| 🏗️ | **reconstruction** | `skills.03_3d_reconstruction.reconstruction` | 对齐 `.h5ad` | VTK 模型 |
| 🎨 | **interpolation** | `skills.03_3d_reconstruction.interpolation` | `.h5ad` + VTK | 表达量映射网格 |
| 📐 | **morphology** | `skills.03_3d_reconstruction.morphology` | `.h5ad` + VTK | 形态学指标 |
| 🧭 | **vectorfield** | `skills.04_morphogenesis.vectorfield` | 两阶段 `.h5ad` | 向量场 + 模型 |
| 📊 | **feature** | `skills.04_morphogenesis.feature` | vectorfield `.h5ad` | GLM 差异基因表 |
| ▶️ | **run_all** | `skills.00_pipeline.run_all` | 原始数据 | 全阶段输出 |

---

## 🔄 流程架构

```
📖 Stage 0 (IO)          🔗 Stage 1 (Alignment)        🏗️ Stage 2 (Reconstruction)
原始数据 ──→ AnnData      AnnData ──→ 对齐 AnnData       对齐 AnnData ──→ VTK 模型
                │                        │                    ├── 🎨 interpolation
                auto_read                two_slice            ├── 📐 morphology
                manual_read              multi_slice          └── reconstruction
                                                                │
                                                                ↓
🧭 Stage 3 (Morphogenesis)
对齐 AnnData + VTK 模型 ──→ 向量场
         │                            ├── 📊 feature (DEG)
         └── vectorfield              └── vectorfield
```

阶段间数据格式：**AnnData `.h5ad`**，对齐坐标存储于 `.obsm` 中。

---

## 🎯 支持的平台

| 平台 | 识别特征 |
|---|---|
| 🔬 MERFISH | `cell_by_gene.csv` + `cell_metadata.csv` |
| 🔬 Xenium | `cell_boundaries.csv` + `gene_panel.json` |
| 🔬 Visium | `spatial/scalefactors_json.json` |
| 🔬 Visium HD | `spatial/tissue_positions.parquet` |
| 🔬 Slide-seq | `MappedDGEForR.csv` |
| 🔬 Stereo-seq | `.gem` 文件 |
| 🔬 STARmap+ | `processed_expression_pd.csv` |
| 🔬 seqFISH | 文件名含 `CxG` |
| 🔬 Open-ST | `.h5ad` 文件 |

---

## 🧩 Skill 模式

所有 Skill 遵循统一接口：

```python
from dataclasses import dataclass

@dataclass
class SkillConfig:
    """类型化参数 — 全部带默认值。"""
    input_path: str
    out_dir: str = "./output/"
    # ...

@dataclass
class SkillResult:
    """类型化输出 — 路径为字符串，notes 用于诊断。"""
    output_dir: str
    notes: list[str]

def run_skill_pipeline(config: SkillConfig) -> SkillResult:
    """端到端入口。首次调用时懒加载 spateo。"""
```

### 🪶 懒加载导入

Skill 模块导入时不会加载 spateo、torch 或 pyvista。重型依赖延迟到 `run_*_pipeline()` 中加载。安全用于探索和局部流程构建。

```python
import importlib
skill = importlib.import_module("skills.02_slice_alignment.two_slice")
# ✅ spateo 尚未加载
result = skill.run_two_slice_pipeline(cfg)
# ✅ spateo 在此加载
```

---

## ⚠️ 已知问题

| 问题 | 阶段 | 状态 | 解决方式 |
|---|---|---|---|
| 表面网格段错误（pymeshfix 原生库） | 🏗️ 2 | 🟡 OPEN | 仅用 `models=["point-cloud"]` |
| Spateo `paste_pairwise_align` 参数不匹配 | 🧭 3 | 🟢 PATCHED | `patch_align_preprocess_kwarg()` 自动修补 |
| 无 GPU（CUDA 驱动过旧） | 全部 | 🟡 ENV | 使用 CPU 运行（较慢） |

---

## 📁 项目结构

```
Spateo-Skills/
├── skills/
│   ├── 00_shared/          # lazy_imports, data_helpers, configure_logging
│   ├── 00_env_setup/       # setup.sh, environment.yml, uv_requirements.txt
│   ├── 00_pipeline/        # run_all.py — 全流程编排器
│   ├── 01_data_io/         # auto_read.py, manual_read.py
│   ├── 02_slice_alignment/ # two_slice.py, multi_slice.py, two_stage_3d.py
│   ├── 03_3d_reconstruction/ # reconstruction.py, interpolation.py, morphology.py
│   └── 04_morphogenesis/   # vectorfield.py, feature.py
├── docs/                   # PYTHON-PATTERNS.md, CODEMAP.md
├── tests/                  # GOTACH.md（问题记录）
└── pyproject.toml
```

---

## 📚 开发

- 📖 编码规范 → [`README.md`](README.md)
- 🗺️ 文件索引 → [`architecture/CODEMAP.md`](architecture/CODEMAP.md)
- 🧪 测试结果 → [`../tests/GOTACH.md`](../tests/GOTACH.md)

---

<div align="center">

用 ❤️ 为空间转录组而生

</div>
