"""
작업 1 - 2단계: 스케일링 → PCA → 이웃그래프(kNN) → Leiden 클러스터링 → UMAP

입력 : data/pbmc3k_step1_qc_norm_hvg.h5ad  (X = log-normalized, HVG 플래그 포함)
출력 : data/pbmc3k_step2_clustered.h5ad, figures/05~09_*.png

[전체 흐름과 이유]
  13,714차원(유전자) 공간에서는 세포 간 거리가 노이즈에 묻힌다("차원의 저주").
  HVG 1,838개 → PCA로 수십 개 축으로 압축 → 그 공간에서 "비슷한 세포끼리" 연결한 그래프 생성
  → 그래프에서 촘촘히 연결된 덩어리를 찾는 것이 Leiden 클러스터링, 같은 그래프를 2D로 펼친 것이 UMAP.
  ※ 클러스터링과 UMAP은 둘 다 "kNN 그래프"에서 나온 형제 결과다. UMAP 그림을 보고 클러스터를 나누는 게 아니다.

[실행]
    .venv\\Scripts\\python.exe 01_scrna_pbmc3k\\step2_pca_cluster_umap.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = Path(__file__).resolve().parent / "figures"

sc.settings.verbosity = 2
sc.set_figure_params(dpi=100, facecolor="white")

SEED = 0            # PCA/Leiden/UMAP 모두 난수를 쓰므로 고정해야 재현된다
N_PCS = 10          # 이웃그래프에 쓸 PC 개수 (아래 elbow plot으로 근거 확인)
N_NEIGHBORS = 10    # 세포 하나당 연결할 이웃 수 (k). 작을수록 국소 구조, 클수록 큰 구조
RESOLUTIONS = [0.3, 0.5, 0.8, 1.0, 1.5]   # Leiden 해상도 후보: 클수록 클러스터가 잘게 쪼개진다
FINAL_RES = 0.5     # 최종 채택 해상도 (근거는 README 2-3: 0.8부터는 뚜렷한 마커 없이 T세포 덩어리만 더 쪼개짐)

adata = sc.read_h5ad(DATA / "pbmc3k_step1_qc_norm_hvg.h5ad")
print(f"[로드] {adata.shape}, HVG {int(adata.var['highly_variable'].sum())}개")

# ----------------------------------------------------------------------------
# 1. 스케일링 (유전자별 z-score)
# ----------------------------------------------------------------------------
# log-normalized 값은 마커 유전자 확인/차등발현에 계속 써야 하므로 따로 보존한다.
adata.layers["lognorm"] = adata.X.copy()

# PCA는 분산이 큰 방향을 찾는다. 스케일링을 안 하면 원래 발현량이 큰 유전자(예: MALAT1, 리보솜 유전자)가
# 단지 숫자가 크다는 이유로 PC를 지배한다. 유전자마다 평균 0, 분산 1로 맞춰 "모든 유전자에 같은 발언권"을 준다.
# max_value=10: 극소수 세포에서만 튀는 값(z > 10)이 PC를 끌고 가지 않도록 잘라낸다.
sc.pp.scale(adata, max_value=10)

# ----------------------------------------------------------------------------
# 2. PCA
# ----------------------------------------------------------------------------
# highly_variable 플래그가 있으면 scanpy는 HVG만 사용해 PCA를 계산한다.
# svd_solver="arpack": 결정론적(같은 입력 → 같은 결과) 솔버.
sc.tl.pca(adata, n_comps=50, svd_solver="arpack", random_state=SEED)

var_ratio = adata.uns["pca"]["variance_ratio"]
print("\n[PCA] PC별 설명 분산 비율(%) 상위 15개")
print("  " + " ".join(f"{v * 100:.2f}" for v in var_ratio[:15]))
print(f"  누적: PC1-10 = {var_ratio[:10].sum() * 100:.1f}%, PC1-50 = {var_ratio.sum() * 100:.1f}%")

# 각 PC를 끌고 가는 유전자(loading 상위) — PC가 어떤 생물학적 축인지 해석하는 근거
loadings = adata.varm["PCs"]
hvg_names = adata.var_names
print("\n[PC loading 상위 유전자]")
for i in range(5):
    order = np.argsort(loadings[:, i])
    pos = ", ".join(hvg_names[order[-6:][::-1]])
    neg = ", ".join(hvg_names[order[:6]])
    print(f"  PC{i + 1} (+): {pos}")
    print(f"  PC{i + 1} (-): {neg}")

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
plt.axvline(N_PCS + 0.5, color="red", ls="--", lw=1)
plt.savefig(FIG / "05_pca_elbow.png", dpi=150, bbox_inches="tight")
plt.close("all")

sc.pl.pca(adata, color=["CST3", "NKG7", "MS4A1", "PPBP"], layer="lognorm", ncols=4, show=False)
plt.savefig(FIG / "06_pca_markers.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 3. 이웃 그래프 (kNN graph)
# ----------------------------------------------------------------------------
# 상위 N_PCS개 PC 공간에서 각 세포의 가장 가까운 k개 세포를 찾아 연결한다.
# 결과는 adata.obsp["connectivities"](가중치 그래프)에 저장되고, Leiden과 UMAP이 이것을 공유한다.
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=SEED)

# ----------------------------------------------------------------------------
# 4. Leiden 클러스터링
# ----------------------------------------------------------------------------
# 그래프에서 "내부 연결은 촘촘하고 외부 연결은 성긴" 커뮤니티를 찾는 알고리즘(Louvain의 개선판).
# resolution은 정답이 있는 값이 아니라 "얼마나 세밀하게 볼지"를 정하는 손잡이 → 여러 값을 비교한다.
print("\n[Leiden] 해상도별 클러스터 수 / 클러스터 크기")
for res in RESOLUTIONS:
    key = f"leiden_r{res}"
    sc.tl.leiden(adata, resolution=res, key_added=key, flavor="igraph", n_iterations=2,
                 directed=False, random_state=SEED)
    sizes = adata.obs[key].value_counts().sort_index()
    print(f"  res={res:<4} → {sizes.size:>2}개: {sizes.tolist()}")

adata.obs["leiden"] = adata.obs[f"leiden_r{FINAL_RES}"].copy()

# ----------------------------------------------------------------------------
# 5. UMAP (시각화 전용)
# ----------------------------------------------------------------------------
# kNN 그래프의 이웃 관계를 최대한 보존하면서 2D에 배치한다.
# 주의: UMAP에서 "가까운 점은 비슷하다"는 믿을 수 있지만, 덩어리 사이의 거리·덩어리의 크기는 해석하면 안 된다.
sc.tl.umap(adata, random_state=SEED)

sc.pl.umap(adata, color=[f"leiden_r{r}" for r in RESOLUTIONS], ncols=3, legend_loc="on data",
           legend_fontsize=8, show=False)
plt.savefig(FIG / "07_umap_resolutions.png", dpi=150, bbox_inches="tight")
plt.close("all")

sc.pl.umap(adata, color="leiden", legend_loc="on data", title=f"Leiden (res={FINAL_RES})", show=False)
plt.savefig(FIG / "08_umap_leiden.png", dpi=150, bbox_inches="tight")
plt.close("all")

# 클러스터가 생물학이 아니라 품질 차이(UMI 수, MT 비율)로 갈린 것은 아닌지 점검
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "pct_counts_mt"], ncols=3, show=False)
plt.savefig(FIG / "09_umap_qc.png", dpi=150, bbox_inches="tight")
plt.close("all")

print("\n[클러스터별 QC 중앙값] (특정 클러스터만 품질이 나쁘면 '세포타입'이 아니라 '저품질 덩어리'일 수 있음)")
print(adata.obs.groupby("leiden", observed=True)[["total_counts", "n_genes_by_counts", "pct_counts_mt"]]
      .median().round(2).assign(n_cells=adata.obs["leiden"].value_counts().sort_index()))

out = DATA / "pbmc3k_step2_clustered.h5ad"
adata.write_h5ad(out)
print(f"\n[저장] {out}")
