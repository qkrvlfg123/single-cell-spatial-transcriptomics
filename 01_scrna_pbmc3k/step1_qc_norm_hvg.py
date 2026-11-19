"""
작업 1 - 1단계: QC → 정규화 → HVG 선택  (데이터: 10x PBMC3k)

[데이터]
- 건강한 공여자 1명의 말초혈액 단핵세포(PBMC) 약 2,700개, 10x Chromium v1, Illumina NextSeq 500.
- scanpy 내장 로더(sc.datasets.pbmc3k)가 10x 공개 서버에서 "filtered" count matrix를 내려받는다.
  filtered = Cell Ranger가 이미 빈 droplet(세포 없는 바코드)을 걸러낸 행렬.
  따라서 여기서의 QC는 "빈 droplet 제거"가 아니라 "죽어가는 세포 / doublet 의심 세포 제거"가 목적.

[실행]  (프로젝트 루트에서)
    .venv\\Scripts\\python.exe 01_scrna_pbmc3k\\step1_qc_norm_hvg.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 창을 띄우지 않고 파일로만 저장
import matplotlib.pyplot as plt
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = Path(__file__).resolve().parent / "figures"
DATA.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

sc.settings.datasetdir = DATA  # 내려받은 원본을 data/ 에 캐시
sc.settings.verbosity = 2
sc.set_figure_params(dpi=100, facecolor="white")

# QC 기준값 — 아래 분포 그림(violin/scatter)을 보고 정한 값. 한곳에 모아 둬야 근거를 설명하기 쉽다.
MIN_GENES = 200      # 유전자가 200개 미만 검출된 바코드 = 빈 droplet/파편에 가까움
MIN_CELLS = 3        # 3개 미만 세포에서만 잡힌 유전자 = 통계적으로 쓸 수 없는 노이즈
MAX_GENES = 2500     # 유전자 수가 비정상적으로 많음 = 세포 2개가 한 droplet에 들어간 doublet 의심
MAX_PCT_MT = 5.0     # 미토콘드리아 비율 5% 초과 = 손상/사멸 세포 의심

# ----------------------------------------------------------------------------
# 0. 로드
# ----------------------------------------------------------------------------
adata = sc.datasets.pbmc3k()          # obs = 세포(바코드), var = 유전자, X = raw UMI count (희소행렬)
adata.var_names_make_unique()         # 같은 gene symbol이 중복되면 뒤에 -1 등을 붙여 구분
n_cells_raw, n_genes_raw = adata.shape
print(f"[로드] 세포 {n_cells_raw:,}개 x 유전자 {n_genes_raw:,}개")

# ----------------------------------------------------------------------------
# 1. QC 지표 계산
# ----------------------------------------------------------------------------
# 사람 미토콘드리아 유전체에 암호화된 유전자는 symbol이 "MT-"로 시작한다 (MT-CO1, MT-ND1 ... 13개 단백질 유전자).
#   세포막이 손상되면 세포질 mRNA는 빠져나가지만, 미토콘드리아는 막에 싸여 있어 그 안의 전사체는 남는다.
#   → 손상/사멸 세포일수록 전체 UMI 중 MT 비율이 올라간다. 그래서 pct_counts_mt가 세포 품질 지표.
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# 리보솜 단백질 유전자(RPS*/RPL*)는 필터 기준으로는 쓰지 않고 참고용으로만 계산.
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)
# 생성되는 주요 컬럼
#   n_genes_by_counts : 세포당 검출된 유전자 수 (count > 0)
#   total_counts      : 세포당 총 UMI 수 (= library size, 시퀀싱 깊이)
#   pct_counts_mt     : 총 UMI 중 미토콘드리아 유전자 비율(%)

print("\n[QC 지표 요약 - 필터 전]")
print(adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]].describe().round(2))

# 필터 전 분포 그림
sc.pl.violin(
    adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4, multi_panel=True, show=False,
)
plt.savefig(FIG / "01_qc_violin_before.png", dpi=150, bbox_inches="tight")
plt.close("all")

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", ax=axes[0], show=False)
axes[0].axhline(MAX_PCT_MT, color="red", ls="--", lw=1)
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", ax=axes[1], show=False)
axes[1].axhline(MAX_GENES, color="red", ls="--", lw=1)
axes[1].axhline(MIN_GENES, color="red", ls="--", lw=1)
fig.tight_layout()
fig.savefig(FIG / "02_qc_scatter_before.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 2. 필터링 (단계별로 몇 개가 빠지는지 기록 → README에 그대로 인용)
# ----------------------------------------------------------------------------
log = {"시작": adata.shape}

sc.pp.filter_cells(adata, min_genes=MIN_GENES)
sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
log[f"min_genes>={MIN_GENES}, min_cells>={MIN_CELLS}"] = adata.shape

n_doublet_like = int((adata.obs["n_genes_by_counts"] >= MAX_GENES).sum())
n_high_mt = int((adata.obs["pct_counts_mt"] >= MAX_PCT_MT).sum())
adata = adata[adata.obs["n_genes_by_counts"] < MAX_GENES].copy()
log[f"n_genes<{MAX_GENES}"] = adata.shape
adata = adata[adata.obs["pct_counts_mt"] < MAX_PCT_MT].copy()
log[f"pct_mt<{MAX_PCT_MT}"] = adata.shape

print("\n[필터링 경과]  (세포 수, 유전자 수)")
for k, v in log.items():
    print(f"  {k:<32} {v[0]:>6,} cells  {v[1]:>7,} genes")
print(f"  - 유전자 수 과다(doublet 의심)로 제거: {n_doublet_like}개")
print(f"  - 미토콘드리아 비율 과다로 제거(위 조건과 중복 가능): {n_high_mt}개")

sc.pl.violin(
    adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4, multi_panel=True, show=False,
)
plt.savefig(FIG / "03_qc_violin_after.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ----------------------------------------------------------------------------
# 3. 정규화
# ----------------------------------------------------------------------------
# raw count는 따로 보관. (나중에 차등발현의 일부 방법·deconvolution 레퍼런스는 raw count가 필요)
adata.layers["counts"] = adata.X.copy()

# (1) library size 정규화: 세포마다 총 UMI가 다르다(포획 효율·시퀀싱 깊이 차이 = 기술적 변동).
#     각 세포의 총합을 10,000으로 맞춰 "세포 내 상대 발현량"으로 바꾼다 (CP10K).
sc.pp.normalize_total(adata, target_sum=1e4)
# (2) log1p 변환: count 데이터는 평균이 클수록 분산도 커진다(소수의 고발현 유전자가 거리를 지배).
#     log(1+x)로 분산을 안정화하고, 발현 차이를 fold-change 스케일로 해석할 수 있게 한다.
sc.pp.log1p(adata)

# ----------------------------------------------------------------------------
# 4. HVG (highly variable genes) 선택
# ----------------------------------------------------------------------------
# 대부분의 유전자는 모든 세포에서 비슷하게 발현(housekeeping)되거나 거의 0 → 세포타입 구분에 정보가 없다.
# "평균 발현량이 비슷한 유전자들끼리 비교했을 때" 유독 세포 간 변동이 큰 유전자만 골라
# 이후 PCA/클러스터링의 입력으로 쓴다. (노이즈 감소 + 계산량 감소)
# flavor="seurat": log 정규화된 값에서 평균 구간(bin)별로 dispersion을 z-score화하는 방식.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, flavor="seurat")
n_hvg = int(adata.var["highly_variable"].sum())
print(f"\n[HVG] {n_hvg:,}개 선택 / 전체 {adata.n_vars:,}개")

top_hvg = adata.var[adata.var["highly_variable"]].sort_values("dispersions_norm", ascending=False).head(20)
print("[HVG 상위 20개 (정규화 dispersion 기준)]")
print("  " + ", ".join(top_hvg.index))

sc.pl.highly_variable_genes(adata, show=False)
plt.savefig(FIG / "04_hvg.png", dpi=150, bbox_inches="tight")
plt.close("all")

# 유전자를 HVG로 잘라내지는 않는다(subset 안 함). highly_variable 플래그만 남겨 두면
# PCA는 HVG만 쓰고, 마커 유전자 확인·차등발현은 전체 유전자로 할 수 있다.
out = DATA / "pbmc3k_step1_qc_norm_hvg.h5ad"
adata.write_h5ad(out)
print(f"\n[저장] {out}  (shape={adata.shape})")
