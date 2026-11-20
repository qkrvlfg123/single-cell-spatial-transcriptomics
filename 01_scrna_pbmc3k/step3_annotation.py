"""
작업 1 - 3단계: 마커 유전자 기반 세포타입 annotation

입력 : data/pbmc3k_step2_clustered.h5ad  (obs["leiden"] = 9개 클러스터)
출력 : data/pbmc3k_step3_annotated.h5ad, results/*.csv, figures/10~13_*.png

[두 가지 근거를 교차 확인한다]
  (A) data-driven : 클러스터별로 "나머지 세포 대비 높게 발현되는 유전자"를 통계 검정으로 뽑는다 (편견 없음)
  (B) knowledge-driven : 면역학에서 이미 알려진 계통 마커가 어느 클러스터에서 켜지는지 본다
  (A)의 상위 유전자가 (B)의 기대와 일치할 때만 세포타입 이름을 붙인다.

[실행]
    .venv\\Scripts\\python.exe 01_scrna_pbmc3k\\step3_annotation.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
RES = HERE / "results"
RES.mkdir(exist_ok=True)

sc.settings.verbosity = 1
sc.set_figure_params(dpi=100, facecolor="white")

adata = sc.read_h5ad(DATA / "pbmc3k_step2_clustered.h5ad")

# ----------------------------------------------------------------------------
# 1. (A) 클러스터별 마커 유전자: Wilcoxon rank-sum 검정
# ----------------------------------------------------------------------------
# 각 클러스터 vs 나머지 전체. scRNA-seq 발현값은 0이 많고 정규분포가 아니므로 t-test보다
# 분포 가정이 없는 Wilcoxon을 쓴다. 입력은 스케일된 값이 아니라 log-normalized 값이어야 한다
# (z-score에서는 fold change가 의미를 잃는다).
# 주의: 같은 데이터로 클러스터를 만들고 그 클러스터 간 차이를 검정하므로 p값은 과대평가된다(double dipping).
#       → p값은 순위용으로만 쓰고, 판단은 logFC와 "발현 세포 비율"로 한다.
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", layer="lognorm", use_raw=False, pts=True)

frames = []
for cl in adata.obs["leiden"].cat.categories:
    df = sc.get.rank_genes_groups_df(adata, cl)
    df.insert(0, "cluster", cl)
    frames.append(df)
markers = pd.concat(frames, ignore_index=True)
markers = markers.rename(columns={"names": "gene", "pct_nz_group": "pct_in", "pct_nz_reference": "pct_out"})
# 리보솜/미토콘드리아 유전자는 거의 모든 세포에 있어 "마커"로서 정보가 없으므로 상위 목록에서는 제외
informative = ~markers["gene"].str.match(r"^(RPS|RPL|MT-)")
top = markers[informative & (markers["logfoldchanges"] > 1) & (markers["pvals_adj"] < 0.05)]
top = top.groupby("cluster", observed=True).head(25)
top.round(4).to_csv(RES / "cluster_markers_top25.csv", index=False)

print("[클러스터별 상위 마커]  유전자(log2FC, 클러스터 내 발현비율% / 나머지 발현비율%)")
for cl, g in top.groupby("cluster", observed=True):
    items = [f"{r.gene}({r.logfoldchanges:.1f}, {r.pct_in * 100:.0f}/{r.pct_out * 100:.0f})" for r in g.head(8).itertuples()]
    print(f"  {cl}: " + ", ".join(items))

# ----------------------------------------------------------------------------
# 2. (B) 알려진 계통 마커 패널
# ----------------------------------------------------------------------------
# 각 유전자가 그 세포의 마커인 이유(기능)를 함께 적어 둔다 — 면접에서 유전자 선택 근거를 설명하기 위함.
MARKER_PANEL = {
    # T세포 공통: CD3 복합체 = TCR 신호전달 단위. T세포 정체성의 정의에 가까운 마커
    "T (pan)": ["CD3D", "CD3E"],
    # naive/central memory: CCR7·SELL(CD62L) = 림프절 homing 수용체, LEF1/TCF7 = 미분화 상태 유지 전사인자
    "Naive": ["CCR7", "SELL", "LEF1", "TCF7"],
    # memory/활성화: IL7R(생존 신호), S100A4·ANXA1(항원 경험 T세포에서 증가)
    "Memory": ["IL7R", "S100A4", "ANXA1"],
    # CD8 T: CD8 공수용체(MHC-I 인식), GZMK = effector memory 단계의 granzyme
    "CD8 T": ["CD8A", "CD8B", "GZMK", "CCL5"],
    # 세포독성: perforin(PRF1)이 구멍을 내고 granzyme(GZMB)이 들어가 표적세포 사멸 유도
    "Cytotoxic/NK": ["NKG7", "GNLY", "GZMB", "PRF1", "KLRF1"],
    # B세포: CD79A/B = BCR 신호전달 단위, MS4A1 = CD20(rituximab 표적)
    "B": ["MS4A1", "CD79A", "CD79B"],
    # classical monocyte: CD14 = LPS 공수용체, S100A8/9 = calprotectin, LYZ = lysozyme
    "CD14 Mono": ["CD14", "LYZ", "S100A8", "S100A9"],
    # non-classical monocyte: FCGR3A = CD16(Fcγ 수용체 IIIa). NK에도 있으므로 MS4A7과 함께 본다
    "FCGR3A Mono": ["FCGR3A", "MS4A7"],
    # 수지상세포(cDC2): FCER1A, CLEC10A + MHC class II 고발현(전문 항원제시)
    "DC": ["FCER1A", "CLEC10A", "HLA-DQA1"],
    # 혈소판: PPBP·PF4 = α-granule에 저장되는 케모카인
    "Platelet": ["PPBP", "PF4"],
}

sc.pl.dotplot(adata, MARKER_PANEL, groupby="leiden", layer="lognorm", use_raw=False,
              standard_scale="var", show=False)
plt.savefig(FIG / "10_dotplot_markers_by_cluster.png", dpi=150, bbox_inches="tight")
plt.close("all")

sc.pl.umap(adata, color=["CD3D", "CCR7", "S100A4", "CD8A", "NKG7", "GNLY", "MS4A1", "CD14", "FCGR3A",
                         "FCER1A", "PPBP", "LYZ"],
           layer="lognorm", use_raw=False, ncols=4, cmap="viridis", show=False)
plt.savefig(FIG / "11_umap_marker_genes.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 3. annotation (근거는 README 3단계 표 참고)
# ----------------------------------------------------------------------------
CELL_TYPES = {
    "0": "Naive T",            # CD3D+, CCR7/LEF1/SELL+, 활성화 유전자 낮음. CD4·CD8 naive 혼재
    "1": "CD8 T",              # CD3D+, CD8A/CD8B+, GZMK/CCL5/NKG7+ (effector memory 성향)
    "2": "Memory CD4 T",       # CD3D+, IL7R/LTB+, S100A4/ANXA1 ↑, CD8 낮음
    "3": "NK",                 # CD3D−, GNLY/GZMB/PRF1/NKG7 매우 높음, FCGR3A+
    "4": "B",                  # CD79A/CD79B/MS4A1+, MHC-II 높음
    "5": "FCGR3A+ Mono",       # FCGR3A/MS4A7+, CD14 낮음, LST1/AIF1/IFITM3
    "6": "CD14+ Mono",         # CD14/LYZ/S100A8/S100A9/FCN1
    "7": "DC",                 # FCER1A/CLEC10A+, MHC-II 최고
    "8": "Platelet",           # PPBP/PF4/GNG11, 유전자 수 매우 적음(무핵)
}
adata.obs["cell_type"] = adata.obs["leiden"].map(CELL_TYPES).astype("category")

comp = (adata.obs.groupby(["leiden", "cell_type"], observed=True).size().rename("n_cells").reset_index())
comp["pct"] = (comp["n_cells"] / adata.n_obs * 100).round(1)
comp.to_csv(RES / "cell_type_composition.csv", index=False)
print("\n[세포타입 구성]")
print(comp.to_string(index=False))

sc.pl.umap(adata, color="cell_type", legend_loc="on data", legend_fontsize=8, legend_fontoutline=2,
           title="PBMC3k cell types", show=False)
plt.savefig(FIG / "12_umap_cell_types.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 4. 보조 확인: NK 클러스터 안의 소집단 (2단계에서 res=0.8일 때만 분리되던 14개 세포)
# ----------------------------------------------------------------------------
# 사람 NK는 CD56bright(사이토카인 분비형, 림프절 homing)와 CD56dim(세포독성형, CD16+)으로 나뉜다.
nk = adata[adata.obs["leiden"] == "3"].copy()
nk.obs["nk_sub"] = np.where(nk.obs["leiden_r0.8"] == "2", "subset (n=14)", "main")
genes = ["GZMK", "XCL1", "KLRC1", "SELL", "FCGR3A", "FGFBP2", "GZMB", "PRF1", "NCAM1"]
X = pd.DataFrame(nk[:, genes].layers["lognorm"].toarray(), columns=genes, index=nk.obs_names)
nk_tab = (X > 0).groupby(nk.obs["nk_sub"]).mean().mul(100).round(0).T
nk_tab.loc["n_cells"] = nk.obs["nk_sub"].value_counts()
nk_tab.to_csv(RES / "nk_subset_pct_expressing.csv")
print("\n[NK 소집단: 발현 세포 비율(%)]")
print(nk_tab)

sc.pl.dotplot(adata, MARKER_PANEL, groupby="cell_type", layer="lognorm", use_raw=False,
              standard_scale="var", show=False)
plt.savefig(FIG / "13_dotplot_markers_by_cell_type.png", dpi=150, bbox_inches="tight")
plt.close("all")

out = DATA / "pbmc3k_step3_annotated.h5ad"
adata.write_h5ad(out)
print(f"\n[저장] {out}")
