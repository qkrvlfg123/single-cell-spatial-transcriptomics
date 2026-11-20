"""
작업 2 - 4단계: Deconvolution — 단일세포 레퍼런스로 spot별 세포타입 조성 추정 (NNLS)

입력 : data/visium_step3_svg.h5ad (layers["counts"], obs["domain"])
       squidpy 내장 마우스 피질 scRNA-seq 레퍼런스 (sq.datasets.sc_mouse_cortex, 세포타입 라벨 포함)
출력 : data/visium_step4_deconv.h5ad, results/deconv_*.csv, figures/18~22_*.png

[배경]
  Visium spot 하나에는 여러 세포가 들어 있어 발현값이 세포타입들의 합이다.
  3단계의 한계(조성 차이와 세포 내 발현 변화를 구분하지 못함)를 풀려면 spot을 세포타입 비율로 분해해야 한다.

[방법: 비음수 최소제곱(NNLS)]
  spot의 발현 벡터 y(유전자 G개)를 세포타입별 평균 발현 프로파일 S(G x K)의 가중합으로 설명한다.
        y ≈ S · w,   w ≥ 0        →  w를 합이 1이 되도록 나눈 값을 세포타입 비율로 해석
  - S는 단일세포 레퍼런스에서 세포타입별 평균 발현으로 만든다.
  - 유전자는 세포타입 마커만 쓴다(모든 타입에서 비슷한 유전자는 분해에 정보를 주지 않고 노이즈만 더한다).
  - 비율 w가 음수일 수 없으므로 일반 최소제곱이 아니라 NNLS를 쓴다.
  - log 변환한 값이 아니라 선형 스케일(CP10K)에서 푼다. 혼합은 선형 공간에서 더해지기 때문이다.

[검증 순서]
  1) 레퍼런스 세포를 섞어 만든 가상 spot(정답 비율을 아는 데이터)에서 방법 자체의 정확도 확인
  2) 실제 Visium에서 피질 층별 뉴런이 2단계의 층 도메인에 맞게 배치되는지 확인
  3) spot별 적합도(R²)를 조직 위에 그려, 레퍼런스에 없는 세포타입으로 이루어진 영역을 확인

[실행]
    .venv\\Scripts\\python.exe 02_spatial_visium\\step4_deconvolution.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import squidpy as sq
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
RES = HERE / "results"

sc.settings.verbosity = 1
sc.set_figure_params(dpi=100, facecolor="white")

SEED = 0
LABEL = "cell_subclass"   # 레퍼런스의 세포타입 라벨 열
MIN_CELLS_PER_TYPE = 50   # 세포 수가 이보다 적은 타입은 평균 프로파일이 불안정하므로 제외
N_MARKERS = 60            # 세포타입당 마커 유전자 수
TEST_FRAC = 0.3           # 가상 spot 검증용으로 떼어 두는 레퍼런스 세포 비율
N_SYNTH = 1000            # 가상 spot 개수
CELLS_PER_SYNTH = 8       # 가상 spot 하나에 섞는 세포 수 (Visium spot의 1~10개 범위)

rng = np.random.default_rng(SEED)


def cp10k(X):
    """세포(또는 spot)별 총합을 10,000으로 맞춘 선형 스케일 값."""
    X = sp.csr_matrix(X, dtype=np.float64)
    tot = np.asarray(X.sum(axis=1)).ravel()
    return sp.diags(1e4 / np.maximum(tot, 1)) @ X


def deconvolve(Y, S):
    """각 행(spot)에 대해 NNLS를 풀어 (비율, R²)를 돌려준다. Y: spot x gene, S: gene x type."""
    W = np.zeros((Y.shape[0], S.shape[1]))
    r2 = np.zeros(Y.shape[0])
    for i in range(Y.shape[0]):
        y = Y[i]
        w, resid = nnls(S, y)
        W[i] = w
        ss_tot = ((y - y.mean()) ** 2).sum()
        r2[i] = 1 - resid ** 2 / ss_tot if ss_tot > 0 else np.nan
    prop = W / np.maximum(W.sum(axis=1, keepdims=True), 1e-12)
    return prop, r2


# ----------------------------------------------------------------------------
# 0. 로드
# ----------------------------------------------------------------------------
adata = sc.read_h5ad(DATA / "visium_step3_svg.h5ad")
ref = sq.datasets.sc_mouse_cortex(path=DATA / "sc_mouse_cortex.h5ad")
# 레퍼런스의 X는 가공된 값이므로 raw count가 있으면 그것을 쓴다
if ref.raw is not None:
    ref = ref.raw.to_adata()
ref.var_names_make_unique()
print(f"[Visium] {adata.shape}   [레퍼런스] {ref.shape}")

counts = ref.obs[LABEL].value_counts()
keep_types = counts[counts >= MIN_CELLS_PER_TYPE].index.tolist()
ref = ref[ref.obs[LABEL].isin(keep_types)].copy()
ref.obs[LABEL] = ref.obs[LABEL].astype(str).astype("category")
types = ref.obs[LABEL].cat.categories.tolist()
safe = {t: t.replace("/", "-") for t in types}   # h5ad 키에는 "/"를 쓸 수 없다 (L2/3 IT → L2-3 IT)
print(f"[세포타입] {len(types)}개 (세포 {MIN_CELLS_PER_TYPE}개 미만 타입 제외): {dict(counts[keep_types])}")

common = ref.var_names.intersection(adata.var_names)
ref = ref[:, common].copy()
print(f"[공통 유전자] {len(common):,}개")

# 레퍼런스를 signature용(train)과 가상 spot용(test)으로 나눈다. 같은 세포로 만들고 평가하면 정확도가 과대평가된다.
is_test = np.zeros(ref.n_obs, dtype=bool)
for t in types:
    idx = np.where(ref.obs[LABEL].to_numpy() == t)[0]
    is_test[rng.choice(idx, size=int(len(idx) * TEST_FRAC), replace=False)] = True
train, test = ref[~is_test].copy(), ref[is_test].copy()
print(f"[분할] signature용 {train.n_obs:,}개 / 가상 spot용 {test.n_obs:,}개")

# ----------------------------------------------------------------------------
# 1. 마커 유전자 선택 + signature 행렬
# ----------------------------------------------------------------------------
train.layers["cp10k"] = cp10k(train.X)
train.X = train.layers["cp10k"].copy()
sc.pp.log1p(train)
sc.tl.rank_genes_groups(train, LABEL, method="t-test_overestim_var", n_genes=N_MARKERS)
markers = pd.DataFrame(train.uns["rank_genes_groups"]["names"])
genes = pd.unique(markers.to_numpy().ravel())
print(f"[마커] 타입당 {N_MARKERS}개 → 합집합 {len(genes):,}개")

S = np.column_stack([
    np.asarray(train[train.obs[LABEL].to_numpy() == t, genes].layers["cp10k"].mean(axis=0)).ravel() for t in types
])
# 유전자별로 크기를 맞춘다. 맞추지 않으면 발현량이 큰 몇몇 유전자의 오차만 줄이는 해가 나온다.
gene_scale = S.max(axis=1)
gene_scale[gene_scale == 0] = 1
S_scaled = S / gene_scale[:, None]
pd.DataFrame(S, index=genes, columns=types).round(3).to_csv(RES / "deconv_signature_cp10k.csv")

# ----------------------------------------------------------------------------
# 2. 검증 1: 가상 spot (정답을 아는 혼합)
# ----------------------------------------------------------------------------
test_labels = test.obs[LABEL].to_numpy()
test_counts = sp.csr_matrix(test[:, genes].X)
test_total = np.asarray(test.X.sum(axis=1)).ravel()
truth = np.zeros((N_SYNTH, len(types)))
Y_syn = np.zeros((N_SYNTH, len(genes)))
for i in range(N_SYNTH):
    # 실제 조직처럼 spot마다 2~4개 타입만 섞이도록, 타입을 먼저 고르고 그 안에서 세포를 뽑는다
    chosen = rng.choice(types, size=rng.integers(2, 5), replace=False)
    pool = np.where(np.isin(test_labels, chosen))[0]
    cells = rng.choice(pool, size=CELLS_PER_SYNTH, replace=True)
    for c in cells:
        truth[i, types.index(test_labels[c])] += 1
    # raw count를 더한 뒤 CP10K: 실제 spot에서 mRNA가 섞이는 방식과 같다
    Y_syn[i] = np.asarray(test_counts[cells].sum(axis=0)).ravel() / test_total[cells].sum() * 1e4
truth /= CELLS_PER_SYNTH

prop_syn, r2_syn = deconvolve(Y_syn / gene_scale, S_scaled)
per_type = pd.DataFrame({
    "cell_type": types,
    "n_ref_cells": [int(counts[t]) for t in types],
    "pearson_r": [np.corrcoef(truth[:, k], prop_syn[:, k])[0, 1] for k in range(len(types))],
    "rmse": np.sqrt(((truth - prop_syn) ** 2).mean(axis=0)),
    "mean_true": truth.mean(axis=0),
    "mean_est": prop_syn.mean(axis=0),
})
per_type.round(3).to_csv(RES / "deconv_synthetic_validation.csv", index=False)
overall_r = np.corrcoef(truth.ravel(), prop_syn.ravel())[0, 1]
overall_rmse = np.sqrt(((truth - prop_syn) ** 2).mean())
print(f"\n[가상 spot 검증] 전체 Pearson r = {overall_r:.3f}, RMSE = {overall_rmse:.3f}, 적합도 R² 중앙값 = {np.nanmedian(r2_syn):.3f}")
print(per_type.round(3).to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
axes[0].scatter(truth.ravel(), prop_syn.ravel(), s=3, alpha=0.2, c="#4c72b0", linewidths=0)
axes[0].plot([0, 1], [0, 1], color="red", ls="--", lw=1)
axes[0].set_xlabel("true proportion"); axes[0].set_ylabel("estimated proportion"); axes[0].grid(False)
axes[0].set_title(f"Synthetic spots (r = {overall_r:.2f}, RMSE = {overall_rmse:.3f})")
o = per_type.sort_values("pearson_r")
axes[1].barh(o["cell_type"], o["pearson_r"], color="#4c72b0")
axes[1].set_xlabel("Pearson r (true vs estimated)"); axes[1].grid(False); axes[1].set_xlim(0, 1)
axes[1].set_title("Per cell type")
fig.tight_layout()
fig.savefig(FIG / "18_deconv_synthetic_validation.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 3. 실제 Visium에 적용
# ----------------------------------------------------------------------------
Y = np.asarray(cp10k(adata.layers["counts"])[:, adata.var_names.get_indexer(genes)].todense())
prop, r2 = deconvolve(Y / gene_scale, S_scaled)
prop_df = pd.DataFrame(prop, index=adata.obs_names, columns=[safe[t] for t in types])
adata.obsm["deconv_prop"] = prop_df
adata.obs["deconv_r2"] = r2
for t in types:
    adata.obs[f"prop_{safe[t]}"] = prop_df[safe[t]].to_numpy()
print(f"\n[Visium] 적합도 R² 중앙값 {np.nanmedian(r2):.3f} (가상 spot에서는 {np.nanmedian(r2_syn):.3f})")
print(f"  spot당 비율 5% 이상인 세포타입 수: 중앙값 {int(np.median((prop >= 0.05).sum(axis=1)))}개")

# 도메인별 평균 조성과 적합도
by_dom = prop_df.groupby(adata.obs["domain"], observed=True).mean()
by_dom.round(3).to_csv(RES / "deconv_mean_proportion_by_domain.csv")
r2_dom = adata.obs.groupby("domain", observed=True)["deconv_r2"].median().round(3)
r2_dom.to_csv(RES / "deconv_fit_r2_by_domain.csv")
print("\n[도메인별 적합도 R² 중앙값 / 비율 상위 3개 세포타입]")
for d in by_dom.index:
    top = by_dom.loc[d].sort_values(ascending=False).head(3)
    print(f"  {d:<28} R²={r2_dom[d]:.2f}  " + ", ".join(f"{k} {v:.2f}" for k, v in top.items()))

fig, ax = plt.subplots(figsize=(0.45 * len(types) + 3, 0.38 * len(by_dom) + 2))
im = ax.imshow(by_dom.to_numpy(), aspect="auto", cmap="viridis", vmin=0, vmax=min(1.0, float(by_dom.to_numpy().max())))
ax.set_xticks(range(len(types))); ax.set_xticklabels(types, rotation=90)
ax.set_yticks(range(len(by_dom))); ax.set_yticklabels(by_dom.index)
ax.grid(False)
fig.colorbar(im, ax=ax, label="mean estimated proportion", fraction=0.03)
ax.set_title("Estimated cell-type composition by spatial domain")
fig.savefig(FIG / "19_deconv_domain_heatmap.png", dpi=150, bbox_inches="tight")
plt.close("all")

sq.pl.spatial_scatter(adata, color=[f"prop_{safe[t]}" for t in types], ncols=6, cmap="magma", img_alpha=0.0,
                      figsize=(3.2, 3.2), vmin=0)
plt.savefig(FIG / "20_deconv_all_types_spatial.png", dpi=110, bbox_inches="tight")
plt.close("all")

sq.pl.spatial_scatter(adata, color="deconv_r2", cmap="viridis", img_alpha=0.2, figsize=(6, 6), title="NNLS fit (R²) per spot")
plt.savefig(FIG / "21_deconv_fit_r2_spatial.png", dpi=150, bbox_inches="tight")
plt.close("all")

adata.write_h5ad(DATA / "visium_step4_deconv.h5ad")
print(f"\n[저장] {DATA / 'visium_step4_deconv.h5ad'}")
