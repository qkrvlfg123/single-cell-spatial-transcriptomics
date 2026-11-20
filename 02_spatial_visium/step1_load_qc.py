"""
작업 2 - 1단계: Visium 데이터 로드 + QC  (데이터: 10x V1_Adult_Mouse_Brain)

[데이터]
- 성체 마우스 뇌 coronal(관상면) 절편, fresh frozen, H&E 염색. 10x Genomics 공개 데모(Space Ranger 1.1.0).
- squidpy 로더(sq.datasets.visium)가 filtered count matrix + 조직 이미지 + spot 좌표를 내려받는다.

[Visium이 scRNA-seq과 다른 점 — QC 해석이 달라진다]
- 관측 단위가 "세포"가 아니라 "spot"이다. spot 지름 55 µm, 중심 간격 100 µm, 한 spot에 대략 1~10개 세포가 들어간다.
  → 한 spot의 발현값은 여러 세포(종종 여러 세포타입)의 합이다.
- 각 spot은 슬라이드 위 고정된 좌표(array_row, array_col)를 가진다. 조직을 해리하지 않으므로 위치 정보가 보존된다.
- UMI 수가 적다고 곧 "저품질"이 아니다. 세포 밀도가 낮은 부위(백질, 뇌실 주변)는 원래 RNA가 적다.
  scRNA-seq처럼 기계적으로 자르면 특정 해부학 구조가 통째로 사라진다 → 반드시 "공간 위에" QC 지표를 그려 보고 판단한다.

[실행]
    .venv\\Scripts\\python.exe 02_spatial_visium\\step1_load_qc.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import squidpy as sq

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)

sc.settings.verbosity = 1
sc.set_figure_params(dpi=100, facecolor="white")

SAMPLE = "V1_Adult_Mouse_Brain"
# QC 기준 — 아래 분포와 공간 그림을 보고 정한 값 (근거는 README 1단계)
MIN_COUNTS = 5000    # spot당 총 UMI 하한
MAX_PCT_MT = 30.0    # 미토콘드리아 비율 상한(%)
MIN_CELLS = 10       # 유전자가 검출된 최소 spot 수

# ----------------------------------------------------------------------------
# 0. 로드
# ----------------------------------------------------------------------------
adata = sq.datasets.visium(SAMPLE, base_dir=DATA / "visium", include_hires_tiff=False)
adata.var_names_make_unique()
print(f"[로드] spot {adata.n_obs:,}개 x 유전자 {adata.n_vars:,}개")
# obs["in_tissue"]      : 조직 아래에 있는 spot인지 (filtered matrix라 전부 1)
# obs["array_row/col"]  : 슬라이드 격자 좌표 (육각 격자)
# obsm["spatial"]       : 원본 해상도 이미지 위의 픽셀 좌표
# uns["spatial"]        : H&E 이미지(hires/lowres)와 scale factor
sf = adata.uns["spatial"][SAMPLE]["scalefactors"]
print(f"  in_tissue=1인 spot: {int(adata.obs['in_tissue'].sum()):,} / 슬라이드 전체 4,992 spot")
print(f"  spot 지름(원본 픽셀) {sf['spot_diameter_fullres']:.1f}px → 55µm 이므로 1px ≈ {55 / sf['spot_diameter_fullres']:.2f}µm")

# ----------------------------------------------------------------------------
# 1. QC 지표
# ----------------------------------------------------------------------------
# 마우스 유전자는 첫 글자만 대문자 → 미토콘드리아 유전자는 "mt-" (사람은 "MT-")
adata.var["mt"] = adata.var_names.str.startswith("mt-")
# 헤모글로빈: 혈관/출혈 부위 표시. 뇌 절편에서 혈액 오염 정도를 본다.
adata.var["hb"] = adata.var_names.str.match(r"^Hb[ab]-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "hb"], percent_top=None, log1p=False, inplace=True)

print(f"  미토콘드리아 유전자 {int(adata.var['mt'].sum())}개: {', '.join(adata.var_names[adata.var['mt']])}")
print("\n[QC 지표 요약 - 필터 전]")
qc_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_hb"]
print(adata.obs[qc_cols].describe(percentiles=[0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]).round(1))

# 분포
fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))
for ax, col, cut in zip(axes, qc_cols, [MIN_COUNTS, None, MAX_PCT_MT, None]):
    ax.hist(adata.obs[col], bins=60, color="#4c72b0")
    if cut is not None:
        ax.axvline(cut, color="red", ls="--", lw=1)
    ax.set_title(col)
    ax.grid(False)
fig.tight_layout()
fig.savefig(FIG / "01_qc_hist.png", dpi=150, bbox_inches="tight")
plt.close("all")

# 조직 이미지와 QC 지표를 같은 좌표계에 그린다 — 공간 QC의 핵심
sq.pl.spatial_scatter(adata, color=None, title="H&E", figsize=(6, 6))
plt.savefig(FIG / "02_he_image.png", dpi=150, bbox_inches="tight")
plt.close("all")

sq.pl.spatial_scatter(adata, color=qc_cols, ncols=4, cmap="viridis", img_alpha=0.3, figsize=(5, 5))
plt.savefig(FIG / "03_qc_spatial.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 2. 필터링
# ----------------------------------------------------------------------------
low = adata.obs["total_counts"] < MIN_COUNTS
high_mt = adata.obs["pct_counts_mt"] > MAX_PCT_MT
adata.obs["qc_fail"] = np.select([low & high_mt, low, high_mt], ["low_counts+high_mt", "low_counts", "high_mt"], "pass")
adata.obs["qc_fail"] = adata.obs["qc_fail"].astype("category")
print("\n[QC 판정]")
print(adata.obs["qc_fail"].value_counts().to_string())

# 제거되는 spot이 조직 어디에 있는지 확인: 가장자리·접힌 부위에 흩어져 있으면 기술적 문제,
# 한 해부학 구조에 몰려 있으면 생물학일 수 있으므로 기준을 다시 생각해야 한다.
sq.pl.spatial_scatter(adata, color="qc_fail", img_alpha=0.5, figsize=(6, 6),
                      palette=matplotlib.colors.ListedColormap(["#d62728", "#ff7f0e", "#9467bd", "#d9d9d9"][-adata.obs["qc_fail"].cat.categories.size:]))
plt.savefig(FIG / "04_qc_fail_spatial.png", dpi=150, bbox_inches="tight")
plt.close("all")

n0 = adata.shape
adata = adata[adata.obs["qc_fail"] == "pass"].copy()
sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
print(f"\n[필터링] spot {n0[0]:,} → {adata.n_obs:,}  /  유전자 {n0[1]:,} → {adata.n_vars:,}")

adata.layers["counts"] = adata.X.copy()
out = DATA / "visium_step1_qc.h5ad"
adata.write_h5ad(out)
print(f"[저장] {out}")
