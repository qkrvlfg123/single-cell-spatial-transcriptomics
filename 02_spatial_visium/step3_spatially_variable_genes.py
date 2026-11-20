"""
작업 2 - 3단계: 공간 변이 유전자(SVG, spatially variable genes) 탐색 — Moran's I

입력 : data/visium_step2_domains.h5ad
출력 : data/visium_step3_svg.h5ad, results/svg_*.csv, figures/14~18_*.png

[SVG란]
  발현이 조직 위에서 무작위가 아니라 "공간적 패턴"을 이루는 유전자.
  HVG는 spot 간 변동의 크기만 보고 위치는 보지 않는다. SVG는 가까운 spot끼리 값이 비슷한지를 본다.
  → 변동이 커도 위치와 무관하게 흩어져 있으면 HVG이지만 SVG는 아니다.

[Moran's I]  공간 자기상관 통계량
       I = (N / W) * Σ_ij w_ij (x_i - x̄)(x_j - x̄) / Σ_i (x_i - x̄)²
  w_ij = spot i와 j가 인접하면 1 (2단계에서 만든 육각 격자 그래프), W = 가중치 총합
  직관: "이웃한 spot 쌍의 발현 상관계수". +1에 가까우면 이웃끼리 비슷(뭉친 패턴), 0이면 무작위, 음수면 체스판 무늬.
  유의성은 순열 검정: spot 위치를 무작위로 섞어 I의 귀무분포를 만들고 관측값과 비교한다.

[실행]
    .venv\\Scripts\\python.exe 02_spatial_visium\\step3_spatially_variable_genes.py
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

sc.settings.verbosity = 1
sc.set_figure_params(dpi=100, facecolor="white")

SEED = 0
MIN_SPOT_FRAC = 0.05   # 전체 spot의 5% 이상에서 검출된 유전자만 검정 (너무 희소한 유전자는 I가 불안정)
N_PERMS = 100
I_CUT = 0.3            # "뚜렷한 공간 패턴"으로 볼 Moran's I 기준 (분포를 보고 정함)

adata = sc.read_h5ad(DATA / "visium_step2_domains.h5ad")
adata.X = adata.layers["lognorm"].copy()   # 스케일된 값이 아니라 log-normalized 값으로 계산

# ----------------------------------------------------------------------------
# 1. Moran's I 계산
# ----------------------------------------------------------------------------
frac = np.asarray((adata.layers["counts"] > 0).mean(axis=0)).ravel()
keep = (frac >= MIN_SPOT_FRAC) & ~adata.var_names.str.match(r"^(mt-|Rp[sl])")
genes = adata.var_names[keep].tolist()
print(f"[검정 대상] {len(genes):,}개 유전자 (전체 {adata.n_vars:,}개 중, spot의 {MIN_SPOT_FRAC:.0%} 이상에서 검출)")

sq.gr.spatial_autocorr(adata, mode="moran", genes=genes, n_perms=N_PERMS, n_jobs=1, seed=SEED, show_progress_bar=False)
svg = adata.uns["moranI"].copy()
svg.index.name = "gene"
svg["frac_spots"] = pd.Series(frac, index=adata.var_names).loc[svg.index]
svg["is_hvg"] = adata.var["highly_variable"].loc[svg.index].to_numpy()
svg["hvg_rank"] = adata.var["highly_variable_rank"].loc[svg.index].to_numpy()

# 각 유전자가 어느 도메인에서 가장 높은지 + 그 도메인에 얼마나 국한되는지
X = pd.DataFrame(adata[:, svg.index].X.toarray(), columns=svg.index)
dom_mean = X.groupby(adata.obs["domain"].to_numpy()).mean()
svg["peak_domain"] = dom_mean.idxmax()
# 도메인이 유전자 발현 분산을 설명하는 정도 (η² = 집단 간 제곱합 / 전체 제곱합).
# 1에 가까우면 "도메인 마커형" 패턴, 낮으면 도메인 경계와 무관한 패턴(구배, 도메인 내부 구조 등).
n_k = adata.obs["domain"].value_counts().reindex(dom_mean.index).to_numpy()[:, None]
ss_between = (n_k * (dom_mean - X.mean()) ** 2).sum()
ss_total = ((X - X.mean()) ** 2).sum()
svg["domain_eta2"] = (ss_between / ss_total).round(3)

svg = svg.sort_values("I", ascending=False)
svg.round(5).to_csv(RES / "svg_moran_all.csv")

sig = svg["pval_sim_fdr_bh"] < 0.05
strong = svg["I"] >= I_CUT
print(f"[결과] 순열 검정 FDR<0.05: {int(sig.sum()):,}개 ({sig.mean():.0%})  |  I>={I_CUT}: {int(strong.sum()):,}개")
print(f"  I 분포: 중앙값 {svg['I'].median():.3f}, 90% {svg['I'].quantile(0.9):.3f}, 99% {svg['I'].quantile(0.99):.3f}, 최대 {svg['I'].max():.3f}")
print(f"  순열 검정의 최소 p값 = 1/({N_PERMS}+1) = {1 / (N_PERMS + 1):.4f} → p값으로는 상위 유전자 간 순위를 가릴 수 없어 I 값으로 순위를 매긴다")

print("\n[상위 30개 SVG]  유전자(I, 최고 발현 도메인)")
for r in svg.head(30).itertuples():
    print(f"  {r.Index:<10} I={r.I:.3f}  eta2={r.domain_eta2:.2f}  {r.peak_domain}")

# ----------------------------------------------------------------------------
# 2. SVG vs HVG
# ----------------------------------------------------------------------------
top_n = 500
top_svg = set(svg.head(top_n).index)
hvg_ranked = adata.var.loc[adata.var["highly_variable"]].sort_values("highly_variable_rank")
top_hvg = set(hvg_ranked.index[:top_n])
print(f"\n[SVG vs HVG] 상위 {top_n}개끼리 겹침: {len(top_svg & top_hvg)}개")

hvg_not_svg = svg[svg["is_hvg"] & (svg["I"] < 0.1)].sort_values("hvg_rank").head(15)
svg_not_hvg = svg[~svg["is_hvg"] & strong].head(15)
print("  HVG인데 공간 패턴이 약한 유전자(I<0.1): " + ", ".join(f"{g}({i:.2f})" for g, i in zip(hvg_not_svg.index, hvg_not_svg["I"])))
print(f"  HVG가 아닌데 I>={I_CUT}인 유전자({int((~svg['is_hvg'] & strong).sum())}개): "
      + ", ".join(f"{g}({i:.2f})" for g, i in zip(svg_not_hvg.index, svg_not_hvg["I"])))

# ----------------------------------------------------------------------------
# 3. 도메인과의 관계
# ----------------------------------------------------------------------------
s = svg[strong]
by_dom = s.groupby("peak_domain").agg(n_svg=("I", "size"), top_genes=("I", lambda x: ", ".join(x.index[:6])))
by_dom = by_dom.reindex(adata.obs["domain"].cat.categories).fillna({"n_svg": 0, "top_genes": ""})
by_dom["n_svg"] = by_dom["n_svg"].astype(int)
by_dom.to_csv(RES / "svg_by_peak_domain.csv")
print(f"\n[I>={I_CUT}인 SVG {len(s)}개의 최고 발현 도메인]")
print(by_dom.to_string())

low_eta = s[s["domain_eta2"] < 0.5].head(15)
print(f"\n[도메인으로 잘 설명되지 않는 SVG] I>={I_CUT} 이면서 eta2<0.5: {int((s['domain_eta2'] < 0.5).sum())}개")
print("  " + ", ".join(f"{g}(I={r.I:.2f}, eta2={r.domain_eta2:.2f})" for g, r in low_eta.iterrows()))

# ----------------------------------------------------------------------------
# 4. 그림
# ----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].hist(svg["I"], bins=80, color="#4c72b0")
axes[0].axvline(I_CUT, color="red", ls="--", lw=1)
axes[0].set_xlabel("Moran's I"); axes[0].set_ylabel("genes"); axes[0].set_yscale("log"); axes[0].grid(False)
axes[0].set_title(f"Moran's I of {len(svg):,} genes")
axes[1].scatter(svg["domain_eta2"], svg["I"], s=3, c="#4c72b0", alpha=0.4, linewidths=0)
axes[1].axhline(I_CUT, color="red", ls="--", lw=1)
axes[1].set_xlabel("variance explained by domain (eta²)"); axes[1].set_ylabel("Moran's I"); axes[1].grid(False)
axes[1].set_title("Spatial autocorrelation vs domain structure")
fig.tight_layout()
fig.savefig(FIG / "14_moran_distribution.png", dpi=150, bbox_inches="tight")
plt.close("all")

sq.pl.spatial_scatter(adata, color=svg.index[:12].tolist(), ncols=6, cmap="magma", img_alpha=0.0, figsize=(3.2, 3.2))
plt.savefig(FIG / "15_top_svg.png", dpi=120, bbox_inches="tight")
plt.close("all")

if len(hvg_not_svg) >= 6:
    sq.pl.spatial_scatter(adata, color=hvg_not_svg.index[:6].tolist(), ncols=6, cmap="magma", img_alpha=0.0, figsize=(3.2, 3.2))
    plt.savefig(FIG / "16_hvg_not_svg.png", dpi=120, bbox_inches="tight")
    plt.close("all")

if len(low_eta) >= 6:
    sq.pl.spatial_scatter(adata, color=low_eta.index[:6].tolist(), ncols=6, cmap="magma", img_alpha=0.0, figsize=(3.2, 3.2))
    plt.savefig(FIG / "17_svg_not_domain.png", dpi=120, bbox_inches="tight")
    plt.close("all")

adata.write_h5ad(DATA / "visium_step3_svg.h5ad")
print(f"\n[저장] {DATA / 'visium_step3_svg.h5ad'}")
