# 단일세포 · 공간전사체 분석 포트폴리오

공개 데이터로 scRNA-seq와 공간전사체(Visium)를 직접 분석하고 생물학적으로 해석한 프로젝트.
각 단계의 스크립트에는 수행 내용, 이유, 생물학적 의미를 한국어 주석으로 남겼고, 해석과 한계는 폴더별 README에 정리했다.

| 작업 | 데이터 | 내용 | 문서 |
|---|---|---|---|
| 1. scRNA-seq | 10x PBMC3k (2,700 cells) | QC → 정규화 → 클러스터링 → 세포타입 annotation → 차등발현·경로분석 | [01_scrna_pbmc3k/README.md](01_scrna_pbmc3k/README.md) |
| 2. 공간전사체 | 10x Visium 성체 마우스 뇌 (2,702 spots) | QC → 공간 도메인 → 공간 변이 유전자 → deconvolution | [02_spatial_visium/README.md](02_spatial_visium/README.md) |

진행 상황은 [PROGRESS.txt](PROGRESS.txt) 참고.

## 작업 1 주요 결과

- QC 후 2,638개 세포를 9개 세포타입으로 annotation (T 3종, NK, B, 단핵구 2종, DC, 혈소판)
- Leiden 해상도는 여러 값을 비교해, 뚜렷한 마커 없이 쪼개지는 과분할을 피하는 0.5를 채택
- 단핵구 아형 차등발현에서 시퀀싱 깊이 차이로 인한 편향(DE 유전자 277개 → 깊이 맞춤 후 50개)을 발견하고 downsampling으로 보정
- classical 단핵구는 염증·항균(S100A8/A9, LYZ, CD14), non-classical 단핵구는 Fcγ 수용체 신호·interferon 반응(FCGR3A, IFITM1/2/3)이 특징

## 작업 2 주요 결과

- QC 후 2,683개 spot을 좌표 없이 클러스터링해 19개 공간 도메인을 식별 (피질 L1~L6, 해마 3구역, 시상, 시상하부, 편도체, 줄무늬체, 백질, 뇌실, 수막)
- 인접한 spot이 같은 클러스터에 속하는 비율 0.724 (라벨을 무작위로 섞었을 때 0.064) → 발현만으로 조직 구조가 재현됨을 정량 확인
- Moran's I로 공간 변이 유전자 1,070개(I ≥ 0.3)를 선별. 상위 500개 SVG와 HVG의 겹침은 39%로, 두 기준이 서로 다른 유전자를 고른다는 것을 확인
- NNLS deconvolution으로 spot별 세포타입 조성을 추정. 가상 spot 검증에서 r = 0.989, 실제 데이터에서 피질 네 층 모두 해당 층의 흥분성 뉴런이 1위
- 실제 Visium의 적합도(R² 0.17)가 낮은 원인을 플랫폼 차이와 레퍼런스 범위(피질 한정)로 분석하고, 해석 가능한 영역을 피질·백질로 제한

## 환경

```
uv venv --python 3.12 .venv
uv pip install --python .venv\Scripts\python.exe scanpy anndata leidenalg igraph matplotlib seaborn squidpy gseapy
.venv\Scripts\python.exe 01_scrna_pbmc3k\step1_qc_norm_hvg.py   # 이후 step2 ~ step4 순서대로
```

데이터는 스크립트가 자동으로 내려받아 `data/`에 저장한다(저장소에는 포함하지 않음).
