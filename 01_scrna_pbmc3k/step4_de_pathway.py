"""
작업 1 - 4단계(선택): 차등발현(DE) + 경로분석
  비교: CD14+ (classical) 단핵구  vs  FCGR3A+ (non-classical) 단핵구

입력 : data/pbmc3k_step3_annotated.h5ad
출력 : results/de_*.csv (깊이 맞춤 = 본 결과, naive = 비교용), results/enrich_*.csv, results/gsea_hallmark.csv, figures/14~17_*.png

[왜 이 비교인가]
  3단계의 마커 검정은 "한 클러스터 vs 나머지 전부"라서 계통 차이(골수계 vs 림프계)가 결과를 지배한다.
  같은 계통 안의 두 집단을 직접 비교해야 "둘을 구분 짓는 기능 차이"가 드러난다.
  단핵구는 classical(CD14++CD16−) → intermediate → non-classical(CD14+CD16++)로 성숙하며
  기능이 다르다고 잘 알려져 있어, 분석 결과를 기존 지식과 대조해 검증하기 좋다.

[두 가지 경로분석을 함께 쓴다]
  ORA  (over-representation): 유의한 DE 유전자 "목록"이 특정 gene set에 우연보다 많이 겹치는가 (Fisher 검정)
  GSEA (prerank)           : 컷오프 없이 "전체 유전자 순위"에서 gene set이 한쪽 끝에 몰려 있는가

[실행]
    .venv\\Scripts\\python.exe 01_scrna_pbmc3k\\step4_de_pathway.py
"""
from pathlib import Path

import gseapy as gp
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

sc.settings.verbosity = 1
sc.set_figure_params(dpi=100, facecolor="white")

GROUP, REF = "FCGR3A+ Mono", "CD14+ Mono"   # logFC > 0 이면 FCGR3A+ 단핵구에서 높음
LFC_CUT, PADJ_CUT, MIN_PCT = 1.0, 0.05, 0.10
SEED = 0

adata = sc.read_h5ad(DATA / "pbmc3k_step3_annotated.h5ad")
mono = adata[adata.obs["cell_type"].isin([GROUP, REF])].copy()
print(f"[대상] {GROUP} {int((mono.obs['cell_type'] == GROUP).sum())}개 vs {REF} {int((mono.obs['cell_type'] == REF).sum())}개")

# ----------------------------------------------------------------------------
# 1. 차등발현 (Wilcoxon rank-sum) — 시퀀싱 깊이를 맞춘 뒤 비교
# ----------------------------------------------------------------------------
# 한계 1 (pseudoreplication): 공여자가 1명이라 "세포"를 반복으로 취급한다. 세포들은 독립 표본이 아니므로
#   p값은 실제보다 훨씬 작게 나온다. 여러 공여자가 있다면 공여자별 pseudobulk + DESeq2/edgeR이 표준.
#   → p값은 엄밀한 유의성이 아니라 순위로 보고, 효과크기(log2FC)와 발현 비율을 함께 본다.
#
# 한계 2 (깊이 편향): FCGR3A+ 단핵구는 세포당 UMI 중앙값이 3,756, CD14+는 2,208이다.
#   CP10K 정규화는 "총합"만 맞출 뿐 "0의 개수"는 못 맞춘다. UMI가 적은 세포는 저발현 유전자가 0으로 찍힐
#   확률이 높고, log1p 평균과 Wilcoxon 순위는 0이 많은 쪽을 체계적으로 낮게 본다.
#   → 그대로 비교하면 깊이가 깊은 FCGR3A+ 쪽에서 "거의 모든 유전자가 높다"는 가짜 신호가 생긴다.
#   해결: 모든 세포의 UMI를 같은 수(DEPTH)로 무작위 downsampling한 뒤 정규화·검정한다.
DEPTH = 1500


def run_de(ad):
    """ad.X(log-normalized)로 GROUP vs REF Wilcoxon 검정 + 집단별 발현 세포 비율."""
    sc.tl.rank_genes_groups(ad, "cell_type", groups=[GROUP], reference=REF, method="wilcoxon", use_raw=False)
    df = sc.get.rank_genes_groups_df(ad, GROUP).rename(columns={"names": "gene", "logfoldchanges": "log2FC"})
    is_g = (ad.obs["cell_type"] == GROUP).to_numpy()
    expressed = ad.X > 0
    pct = pd.DataFrame({
        "pct_fcgr3a": np.asarray(expressed[is_g].mean(axis=0)).ravel(),
        "pct_cd14": np.asarray(expressed[~is_g].mean(axis=0)).ravel(),
    }, index=ad.var_names)
    df = df.join(pct, on="gene")
    # 두 집단 어디서도 거의 발현되지 않는 유전자는 logFC가 불안정하므로 제외. 리보솜/MT 유전자도 제외.
    df = df[(df["pct_fcgr3a"] >= MIN_PCT) | (df["pct_cd14"] >= MIN_PCT)]
    df = df[~df["gene"].str.match(r"^(RPS|RPL|MT-)")].copy()
    df["direction"] = "ns"
    sig = df["pvals_adj"] < PADJ_CUT
    df.loc[sig & (df["log2FC"] >= LFC_CUT), "direction"] = "up_FCGR3A"
    df.loc[sig & (df["log2FC"] <= -LFC_CUT), "direction"] = "up_CD14"
    return df


def summarize(df):
    return (f"검정 {len(df):,}개, FCGR3A+↑ {(df['direction'] == 'up_FCGR3A').sum()}개, "
            f"CD14+↑ {(df['direction'] == 'up_CD14').sum()}개, "
            f"통계량>0인 유전자 {(df['scores'] > 0).mean() * 100:.0f}%")


print(mono.obs.groupby("cell_type", observed=True)[["total_counts", "n_genes_by_counts"]].median())

# (a) 보정 없이: 편향을 보여주기 위한 비교용
naive = mono.copy()
naive.X = naive.layers["lognorm"].copy()
de_naive = run_de(naive)
de_naive.sort_values("scores", ascending=False).round(5).to_csv(RES / "de_mono_naive_depth_biased.csv", index=False)
print(f"[DE 보정 전]      {summarize(de_naive)}")

# (b) 깊이 맞춤: 본 분석
matched = mono.copy()
matched.X = matched.layers["counts"].copy()
matched = matched[np.asarray(matched.X.sum(axis=1)).ravel() >= DEPTH].copy()   # DEPTH 미만 세포는 제외
sc.pp.downsample_counts(matched, counts_per_cell=DEPTH, random_state=SEED)
sc.pp.normalize_total(matched, target_sum=1e4)
sc.pp.log1p(matched)
print(f"[깊이 맞춤] 세포당 {DEPTH} UMI → {matched.obs['cell_type'].value_counts().to_dict()}")
de = run_de(matched)
de.sort_values("scores", ascending=False).round(5).to_csv(RES / "de_fcgr3a_vs_cd14_mono.csv", index=False)
print(f"[DE 깊이 맞춤 후] {summarize(de)}")
n_keep = len(set(de.loc[de["direction"] == "up_FCGR3A", "gene"]) & set(de_naive.loc[de_naive["direction"] == "up_FCGR3A", "gene"]))
print(f"  보정 전 FCGR3A+↑ {(de_naive['direction'] == 'up_FCGR3A').sum()}개 중 깊이 맞춤 후에도 남는 것: {n_keep}개")
background = de["gene"].tolist()   # 경로분석의 배경 = 실제로 검정한 유전자

up_f = de[de["direction"] == "up_FCGR3A"].sort_values("scores", ascending=False)
up_c = de[de["direction"] == "up_CD14"].sort_values("scores")
fmt = lambda d: ", ".join(f"{r.gene}({r.log2FC:.1f})" for r in d.head(20).itertuples())
print("  FCGR3A+ ↑ :", fmt(up_f))
print("  CD14+   ↑ :", fmt(up_c))

# volcano plot
fig, ax = plt.subplots(figsize=(7, 5.5))
y = -np.log10(de["pvals_adj"].clip(lower=1e-300))
colors = de["direction"].map({"ns": "#c8c8c8", "up_FCGR3A": "#8c564b", "up_CD14": "#e377c2"})
# 한쪽 집단에서 발현이 0인 유전자는 log2FC가 ±20 이상으로 발산하므로 그림에서는 ±8로 잘라 표시
ax.scatter(de["log2FC"].clip(-8, 8), y, c=colors, s=6, linewidths=0)
ax.set_xlim(-8.5, 8.5)
ax.axvline(LFC_CUT, color="k", ls="--", lw=0.6)
ax.axvline(-LFC_CUT, color="k", ls="--", lw=0.6)
ax.axhline(-np.log10(PADJ_CUT), color="k", ls="--", lw=0.6)
for r in pd.concat([up_f.head(12), up_c.head(12)]).itertuples():
    ax.annotate(r.gene, (r.log2FC, -np.log10(max(r.pvals_adj, 1e-300))), fontsize=7,
                xytext=(2, 2), textcoords="offset points")
ax.set_xlabel("log2 fold change  (← CD14+ Mono | FCGR3A+ Mono →)")
ax.set_ylabel("-log10 adjusted p")
ax.set_title(f"FCGR3A+ vs CD14+ monocytes  ({len(up_f)} up / {len(up_c)} down)")
ax.grid(False)
fig.savefig(FIG / "14_volcano_mono.png", dpi=150, bbox_inches="tight")
plt.close("all")

top_genes = up_f.head(12)["gene"].tolist() + up_c.head(12)["gene"].tolist()
sc.pl.dotplot(mono, top_genes, groupby="cell_type", layer="lognorm", use_raw=False, show=False)  # 그림은 원래 log-normalized 값으로
plt.savefig(FIG / "15_dotplot_mono_de.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 2. ORA: GO Biological Process
# ----------------------------------------------------------------------------
# gene set 라이브러리는 Enrichr에서 내려받아 로컬에서 Fisher 검정을 한다(유전자 목록을 외부로 보내지 않음).
# background를 "사람 전체 유전자"가 아니라 "이 비교에서 실제로 검정한 유전자"로 줘야 한다.
# 전체 유전자를 배경으로 쓰면 단핵구에서 발현되는 유전자라는 이유만으로 면역 관련 term이 전부 유의하게 나온다.
go = gp.get_library("GO_Biological_Process_2023", organism="Human")
ora_tables = {}
for name, genes in [("up_FCGR3A", up_f["gene"].tolist()), ("up_CD14", up_c["gene"].tolist())]:
    enr = gp.enrich(gene_list=genes, gene_sets=go, background=background, outdir=None, no_plot=True)
    res = enr.results.sort_values("Adjusted P-value")
    res = res[res["Overlap"].str.split("/").str[0].astype(int) >= 3]   # 유전자 1~2개로 생긴 term 제외
    res.to_csv(RES / f"enrich_GO_BP_{name}.csv", index=False)
    ora_tables[name] = res
    print(f"\n[ORA GO BP] {name}  (입력 {len(genes)}개, padj<0.05 term {int((res['Adjusted P-value'] < 0.05).sum())}개)")
    for r in res.head(10).itertuples():
        print(f"  {r.Term[:62]:<62} {r.Overlap:>7}  padj={r._5:.1e}  {r.Genes[:70]}")

fig, axes = plt.subplots(1, 2, figsize=(15, 5))
for ax, (name, color) in zip(axes, [("up_FCGR3A", "#8c564b"), ("up_CD14", "#e377c2")]):
    t = ora_tables[name].head(10).iloc[::-1]
    ax.barh([s.split(" (GO")[0][:48] for s in t["Term"]], -np.log10(t["Adjusted P-value"]), color=color)
    ax.axvline(-np.log10(0.05), color="k", ls="--", lw=0.6)
    ax.set_xlabel("-log10 adjusted p")
    ax.set_title(f"GO BP: {name}")
    ax.grid(False)
    ax.tick_params(axis="y", labelsize=8)
fig.tight_layout()
fig.savefig(FIG / "16_go_bp_ora.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 3. GSEA prerank: MSigDB Hallmark (50개 대표 경로)
# ----------------------------------------------------------------------------
# 순위 지표로 Wilcoxon 통계량(scores)을 쓴다: 방향과 일관성을 함께 반영한다.
# NES > 0 : FCGR3A+ 쪽에 몰린 경로,  NES < 0 : CD14+ 쪽에 몰린 경로
rnk = de.set_index("gene")["scores"].sort_values(ascending=False)
hallmark = gp.get_library("MSigDB_Hallmark_2020", organism="Human")
pre = gp.prerank(rnk=rnk, gene_sets=hallmark, min_size=10, max_size=500, permutation_num=1000,
                 seed=SEED, outdir=None, no_plot=True, threads=4, verbose=False)
gsea = pre.res2d.copy()
for c in ["NES", "FDR q-val", "NOM p-val"]:
    gsea[c] = pd.to_numeric(gsea[c])
gsea = gsea.sort_values("NES", ascending=False)
gsea[["Term", "NES", "NOM p-val", "FDR q-val", "Tag %", "Lead_genes"]].to_csv(RES / "gsea_hallmark.csv", index=False)
print("\n[GSEA Hallmark]  FDR < 0.25")
for _, r in gsea[gsea["FDR q-val"] < 0.25].iterrows():
    print(f"  NES={r['NES']:+.2f}  FDR={r['FDR q-val']:.3f}  {r['Term']:<36} {r['Lead_genes'][:80]}")

