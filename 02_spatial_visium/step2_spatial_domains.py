"""
작업 2 - 2단계: 클러스터링으로 공간 도메인(조직 구조) 식별

입력 : data/visium_step1_qc.h5ad
출력 : data/visium_step2_domains.h5ad, results/*.csv, figures/05~11_*.png

[핵심 아이디어]
  클러스터링에는 좌표를 전혀 쓰지 않는다. "발현이 비슷한 spot끼리" 묶기만 한다(scRNA-seq과 같은 절차).
  그런데 그 결과를 조직 위에 다시 그리면 클러스터가 해부학 구조 모양으로 나타난다.
  → 위치 정보를 넣지 않았는데 공간적으로 뭉친다는 것 자체가 "클러스터가 실제 조직 구조를 반영한다"는 독립적인 검증이다.
  scRNA-seq의 클러스터는 "세포타입"이지만, Visium의 클러스터는 "세포 조성이 비슷한 조직 구역(도메인)"이다.

[실행]
    .venv\\Scripts\\python.exe 02_spatial_visium\\step2_spatial_domains.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
RES = HERE / "results"
RES.mkdir(exist_ok=True)

sc.settings.verbosity = 1
sc.set_figure_params(dpi=100, facecolor="white")

SEED = 0
N_HVG = 2000
N_PCS = 30
N_NEIGHBORS = 15
RESOLUTIONS = [0.4, 0.6, 0.8, 1.0, 1.2]
FINAL_RES = 0.8

adata = sc.read_h5ad(DATA / "visium_step1_qc.h5ad")
print(f"[로드] {adata.shape}")

# ----------------------------------------------------------------------------
# 1. 정규화 → HVG → PCA  (scRNA-seq과 동일한 절차)
# ----------------------------------------------------------------------------
# spot마다 총 UMI가 다른 이유에는 기술적 요인(포획 효율)과 생물학적 요인(세포 밀도)이 섞여 있다.
# library size 정규화는 둘을 구분하지 못하고 함께 지운다 → "세포가 많아서 높다"는 정보는 사라지고
# "spot 안에서의 상대 조성"만 남는다. 도메인을 조성으로 구분하려는 목적에는 맞는 선택.
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()

# seurat_v3 flavor: raw count에서 평균-분산 관계를 직접 모델링해 상위 N개를 고른다(개수를 지정할 수 있어 재현이 쉽다).
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=N_HVG, layer="counts")
print(f"[HVG] {int(adata.var['highly_variable'].sum())}개. 상위: "
      + ", ".join(adata.var.sort_values('highly_variable_rank').index[:15]))

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, svd_solver="arpack", random_state=SEED)
vr = adata.uns["pca"]["variance_ratio"]
print(f"[PCA] 누적 설명 분산: PC1-10 {vr[:10].sum() * 100:.1f}%, PC1-{N_PCS} {vr[:N_PCS].sum() * 100:.1f}%")

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
plt.axvline(N_PCS + 0.5, color="red", ls="--", lw=1)
plt.savefig(FIG / "05_pca_elbow.png", dpi=150, bbox_inches="tight")
plt.close("all")

# PC 점수를 조직 위에 그리면 각 PC가 어떤 해부학적 축인지 바로 보인다
for i in range(4):
    adata.obs[f"PC{i + 1}"] = adata.obsm["X_pca"][:, i]
sq.pl.spatial_scatter(adata, color=[f"PC{i + 1}" for i in range(4)], ncols=4, cmap="RdBu_r", img_alpha=0.2, figsize=(4.5, 4.5))
plt.savefig(FIG / "06_pc_spatial.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 2. 발현 기반 kNN 그래프 → Leiden → UMAP   (좌표 미사용)
# ----------------------------------------------------------------------------
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
print("\n[Leiden] 해상도별 클러스터 수")
for res in RESOLUTIONS:
    key = f"leiden_r{res}"
    sc.tl.leiden(adata, resolution=res, key_added=key, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
    print(f"  res={res:<4} → {adata.obs[key].nunique():>2}개: {adata.obs[key].value_counts().sort_index().tolist()}")
adata.obs["cluster"] = adata.obs[f"leiden_r{FINAL_RES}"].copy()

sq.pl.spatial_scatter(adata, color=[f"leiden_r{r}" for r in RESOLUTIONS], ncols=5, img_alpha=0.2, figsize=(4.5, 4.5),
                      legend_loc=None)
plt.savefig(FIG / "07_spatial_resolutions.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 3. 공간 이웃 그래프와 공간 응집도 (여기서 처음으로 좌표를 쓴다 — 검증 용도)
# ----------------------------------------------------------------------------
# Visium은 육각 격자라 각 spot의 물리적 이웃이 최대 6개다. n_rings=1 → 바로 인접한 한 겹만 연결.
sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6, n_rings=1)
A = adata.obsp["spatial_connectivities"]
print(f"\n[공간 그래프] spot당 평균 물리적 이웃 {A.sum() / adata.n_obs:.2f}개")

# 응집도 = "물리적으로 인접한 spot 쌍 중 같은 클러스터인 비율". 좌표 없이 만든 클러스터가 조직 위에서 뭉치는 정도.
# 비교 기준(null): 클러스터 라벨을 무작위로 섞었을 때의 같은 비율.
labels = adata.obs["cluster"].cat.codes.to_numpy()
rows, cols = A.nonzero()
obs_same = float((labels[rows] == labels[cols]).mean())
rng = np.random.default_rng(SEED)
null = [float((p[rows] == p[cols]).mean()) for p in (rng.permutation(labels) for _ in range(200))]
print(f"[공간 응집도] 인접 spot이 같은 클러스터인 비율: 관측 {obs_same:.3f} vs 무작위 {np.mean(null):.3f} (최대 {np.max(null):.3f})")

per = pd.DataFrame({"cluster": adata.obs["cluster"].cat.categories})
per["n_spots"] = adata.obs["cluster"].value_counts().sort_index().to_numpy()
per["same_cluster_neighbor_frac"] = [
    float((labels[cols[labels[rows] == k]] == k).mean()) for k in range(len(per))]
per["expected_by_chance"] = (per["n_spots"] / adata.n_obs).round(3)
per.round(3).to_csv(RES / "cluster_spatial_coherence.csv", index=False)
print(per.round(3).to_string(index=False))

# 이웃 enrichment: 서로 맞닿아 있는 도메인 쌍 확인 (z-score, 순열 검정)
sq.gr.nhood_enrichment(adata, cluster_key="cluster", seed=SEED, show_progress_bar=False)
sq.pl.nhood_enrichment(adata, cluster_key="cluster", cmap="RdBu_r", vmin=-30, vmax=30, figsize=(5, 5))
plt.savefig(FIG / "10_nhood_enrichment.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 4. 그림 + 클러스터별 마커
# ----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 6))
sc.pl.umap(adata, color="cluster", legend_loc="on data", legend_fontsize=9, legend_fontoutline=2, ax=axes[0], show=False,
           title="UMAP (expression only)")
sq.pl.spatial_scatter(adata, color="cluster", img_alpha=0.3, ax=axes[1], title="Same clusters on tissue")
fig.savefig(FIG / "08_umap_vs_spatial.png", dpi=150, bbox_inches="tight")
plt.close("all")

# 클러스터를 하나씩 따로 그려 모양을 확인 (annotation 근거)
cats = adata.obs["cluster"].cat.categories
ncols = 6
nrows = int(np.ceil(len(cats) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows))
xy = adata.obsm["spatial"]
for ax, c in zip(axes.ravel(), cats):
    m = (adata.obs["cluster"] == c).to_numpy()
    ax.scatter(xy[~m, 0], xy[~m, 1], s=2, c="#e0e0e0", linewidths=0)
    ax.scatter(xy[m, 0], xy[m, 1], s=4, c="#d62728", linewidths=0)
    ax.set_title(f"cluster {c} (n={m.sum()})", fontsize=10)
    ax.invert_yaxis(); ax.set_aspect("equal"); ax.axis("off")
for ax in axes.ravel()[len(cats):]:
    ax.axis("off")
fig.tight_layout()
fig.savefig(FIG / "09_clusters_individual.png", dpi=130, bbox_inches="tight")
plt.close("all")

sc.tl.rank_genes_groups(adata, "cluster", method="wilcoxon", layer="lognorm", use_raw=False, pts=True)
frames = []
for c in cats:
    df = sc.get.rank_genes_groups_df(adata, c)
    df.insert(0, "cluster", c)
    frames.append(df)
markers = pd.concat(frames).rename(columns={"names": "gene", "logfoldchanges": "log2FC", "pct_nz_group": "pct_in", "pct_nz_reference": "pct_out"})
top = markers[(markers["log2FC"] > 1) & (markers["pvals_adj"] < 0.05) & ~markers["gene"].str.match(r"^(Rp[sl]|mt-|Gm\d|.*Rik$)")]
top = top.groupby("cluster", observed=True).head(25)
top.round(4).to_csv(RES / "cluster_markers_top25.csv", index=False)
print("\n[클러스터별 상위 마커]  유전자(log2FC)")
for c, g in top.groupby("cluster", observed=True):
    print(f"  {c}: " + ", ".join(f"{r.gene}({r.log2FC:.1f})" for r in g.head(12).itertuples()))

# ----------------------------------------------------------------------------
# 5. 도메인 annotation: (위치·모양) + (마커 유전자) 두 근거가 모두 맞을 때 이름을 붙인다
# ----------------------------------------------------------------------------
# scRNA-seq annotation과 다른 점: 근거가 하나 더 있다. 클러스터가 조직 위에서 어떤 모양으로 어디에 있는지를
# H&E·뇌 지도(Allen Mouse Brain Atlas coronal)와 대조할 수 있다. 상세 근거는 README 2-4.
DOMAINS = {
    "1": "Cortex L1",                  # 가장 바깥 얇은 띠. 뉴런 세포체가 적고 astrocyte 유전자(Slc1a2, Slc1a3) 비중이 큼
    "5": "Cortex L2/3",                # Lamp5, Cux2, Calb1, Rasgrf2
    "16": "Cortex L4",                 # Rorb 최고(피질 클러스터 중). L2/3와 L5 사이
    "8": "Cortex L5",                  # Fezf2, Bcl11b
    "17": "Cortex L6",                 # Tbr1, Ccn2(=Ctgf, L6b), Nr4a2. 백질 바로 바깥
    "0": "Piriform cortex",            # 복외측 가장자리 호. Slc17a7+ 글루탐산성
    "12": "Hippocampus CA1",           # Fibcd1, Wfs1, Spink8 — 바깥쪽 얇은 호
    "11": "Hippocampus DG/CA3",        # C1ql2, Prox1, Dsp(치아이랑 과립세포) + Bok(CA3)
    "10": "Hippocampus neuropil",      # 세포체층 사이. 수상돌기로 운반되는 mRNA(Camk2a, Ddn, Shank1)
    "15": "Thalamus (lateral)",        # Prkcd, Rora, Tnnt1, Slc17a6 — 감각 중계핵 쪽
    "3": "Thalamus (medial)",          # Tcf7l2, Slc17a6, Calb2, Cbln1/2/4 — 정중선 쪽
    "4": "GABAergic nuclei (TRN/ZI)",  # Gad1/2, Slc32a1, Pvalb + 수초 유전자. 시상 가장자리를 감싸는 띠
    "2": "Hypothalamus",               # Pmch, Hcrt(외측 시상하부 신경펩타이드), Dlk1, Baiap3
    "9": "Amygdala (medial/central)",  # Hap1, Cartpt, Penk, Zcchc12, Gad1
    "14": "Amygdala (basolateral)",    # Lypd1, Lmo3, Slc30a3, Slc17a7+
    "18": "Striatum",                  # Ppp1r1b(DARPP-32), Penk, Gpr88, Adora2a — medium spiny neuron
    "6": "White matter",               # Mbp, Plp1, Mobp, Mag — oligodendrocyte 수초 유전자
    "13": "Ventricle (ependyma/ChP)",  # Ttr(맥락총), Foxj1 및 섬모 유전자(뇌실막세포)
    "7": "Meninges",                   # Col1a2, Dcn, Mgp, Slc6a13, Ptgds — 조직 표면에 흩어짐
}
adata.obs["domain"] = adata.obs["cluster"].astype(str).map(DOMAINS).astype("category")
order = list(DOMAINS.values())
adata.obs["domain"] = adata.obs["domain"].cat.reorder_categories(order)

comp = adata.obs.groupby(["cluster", "domain"], observed=True).size().rename("n_spots").reset_index()
comp = comp.merge(per[["cluster", "same_cluster_neighbor_frac"]], on="cluster").set_index("domain").loc[order].reset_index()
comp.round(3).to_csv(RES / "domain_summary.csv", index=False)
print("\n[도메인 요약]")
print(comp.round(3).to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), gridspec_kw={"width_ratios": [1, 1.25]})
sq.pl.spatial_scatter(adata, color=None, ax=axes[0], title="H&E")
sq.pl.spatial_scatter(adata, color="domain", img_alpha=0.25, ax=axes[1], title="Spatial domains", palette="tab20")
fig.savefig(FIG / "11_spatial_domains.png", dpi=150, bbox_inches="tight")
plt.close("all")

PANEL = {
    "Cortex layers": ["Lamp5", "Cux2", "Rorb", "Fezf2", "Tbr1", "Ccn2"],
    "Hippocampus": ["Fibcd1", "Wfs1", "C1ql2", "Prox1", "Ddn"],
    "Thalamus": ["Prkcd", "Tcf7l2", "Slc17a6"],
    "GABA": ["Gad1", "Pvalb"],
    "Hypothal.": ["Pmch", "Hcrt"],
    "Amygdala": ["Hap1", "Lypd1"],
    "Striatum": ["Ppp1r1b", "Penk"],
    "Glia/other": ["Mbp", "Plp1", "Ttr", "Foxj1", "Col1a2"],
    "Glu": ["Slc17a7"],
}
sc.pl.dotplot(adata, PANEL, groupby="domain", layer="lognorm", use_raw=False, standard_scale="var", show=False)
plt.savefig(FIG / "12_dotplot_domain_markers.png", dpi=150, bbox_inches="tight")
plt.close("all")

sq.pl.spatial_scatter(adata, color=["Lamp5", "Rorb", "Fezf2", "Tbr1", "Fibcd1", "C1ql2", "Prkcd", "Pmch", "Ppp1r1b", "Mbp", "Ttr", "Gad1"],
                      layer="lognorm", use_raw=False, ncols=6, cmap="magma", img_alpha=0.0, figsize=(3.2, 3.2))
plt.savefig(FIG / "13_marker_genes_spatial.png", dpi=120, bbox_inches="tight")
plt.close("all")

adata.write_h5ad(DATA / "visium_step2_domains.h5ad")
print(f"\n[저장] {DATA / 'visium_step2_domains.h5ad'}")
