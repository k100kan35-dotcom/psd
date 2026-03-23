# PSD (Power Spectral Density) 계산 프로그램 사양서

노면 거칠기 프로파일(1D)에서 2D 등방성 PSD C(q)를 계산하기 위한 표준 파이프라인 사양.

---

## 1. 입력 데이터

| 항목 | 설명 |
|------|------|
| 포맷 | 2열 텍스트 (x[m], h[m]) 또는 1열 (h[m]만, dx 별도 지정) |
| 단위 | SI (m) |
| 프로파일 | 등간격 1D 높이 프로파일 h(x) |

---

## 2. 전처리 (Preprocessing)

### 2.1 Detrend
- **기본값: `mean`** (평균 제거)
- 선택지: `none`, `mean`, `linear`, `quadratic`
- `mean`: h(x) → h(x) - mean(h)
- `linear`: 최소제곱 직선 피팅 후 제거
- `quadratic`: 2차 다항식 피팅 후 제거

### 2.2 Top PSD (선택)
- 프로파일에서 상위 접촉면(bearing area) 영역만 추출하여 PSD 계산
- 기본값: **OFF**

---

## 3. 1D PSD 계산

### 3.1 Window / Method
- **기본값: `multitaper`**
- 선택지: `none`, `multitaper`, `welch`, `hanning`, `hamming`, `blackman`

#### multitaper (권장)
- DPSS (Discrete Prolate Spheroidal Sequences) 테이퍼 사용
- NW (half-bandwidth parameter) = 4, 테이퍼 수 K = 2×NW - 1 = 7
- 각 테이퍼로 FFT 후 PSD를 가중 평균
- 스펙트럼 누출 억제에 가장 효과적

#### 단일 윈도우 (hanning, hamming, blackman)
- 해당 윈도우 함수를 프로파일에 적용 후 단일 FFT
- 윈도우 에너지 보상: C1D /= mean(w²)

#### welch
- 프로파일을 겹치는 세그먼트로 분할 (기본 nperseg=N//8, overlap=50%)
- 각 세그먼트에 Hanning 윈도우 적용 후 평균

#### none
- 윈도우 없이 직접 FFT

### 3.2 1D PSD 공식

```
C1D(q) = (dx / (2π·N)) × |FFT(h)|²
```

여기서:
- q = 2π·f (공간 파수, rad/m)
- dx = 샘플 간격 (m)
- N = 프로파일 포인트 수
- 양의 q 영역만 사용 (q > 0)

### 3.3 sinc² 보정
- **기본값: ON (적용)**
- 이산 샘플링에 의한 sinc²(f·dx) 감쇠를 보상
- Nyquist 근처에서 PSD가 인위적으로 떨어지는 것을 방지

```
x = q / (2·q_max)
sinc(x) = sin(πx) / (πx)
C1D_corrected = C1D / sinc²(x)
```

---

## 4. 1D → 2D 변환

1D 라인 스캔 PSD를 2D 등방성(isotropic) 표면 PSD로 변환.

### 4.1 변환 방법

#### sqrt 방법 (기본값, 권장)

```
C2D(q) = C1D(q) / (π·q) × √(1 + 3H)
```

- 보정 계수(Correction Factor): **√(1 + 3H)**
- H=0.80일 때: √(1 + 2.4) = √3.4 ≈ **1.8439**

#### gamma 방법

```
C2D(q) = (C1D(q) / q) × Γ(1+H) / [√π · Γ(H+½)]
```

- Γ: 감마 함수
- H=0.80일 때: 보정 계수 ≈ 0.5765

#### standard 방법

```
C2D(q) = C1D(q) / (π·q) × corr
```

- corr: 사용자 직접 지정

#### none

```
C2D(q) = C1D(q)   (변환 없음, 1D PSD 그대로)
```

### 4.2 Hurst 지수 (H)
- **기본값: 0.80**
- 범위: 0.1 ~ 1.0
- 자기 친화(self-affine) 거칠기의 스케일링 지수
- sqrt 및 gamma 방법의 보정 계수 계산에 사용

### 4.3 보정 계수 자동 계산
- `sqrt` → corr = √(1 + 3H), 자동 계산, 사용자 수정 불가
- `gamma` → corr = Γ(1+H)/[√π·Γ(H+½)], 자동 계산, 사용자 수정 불가
- `standard` → corr = 사용자 입력
- H 값 변경 시 보정 계수 자동 업데이트

---

## 5. 출력 비닝 (Output Binning)

### 5.1 Log-bin (Persson 방식)
- **기본 빈 수: 88**
- log10(q) 공간에서 등간격 빈 생성
- 각 빈 내 기하 평균(geometric mean) 사용 (log 공간에서 산술 평균)
- 저주파에서 FFT 주파수 간격이 빈 폭보다 클 경우 개별 빈 할당

---

## 6. 권장 기본 설정 요약

| 파라미터 | 기본값 | 비고 |
|----------|--------|------|
| Detrend | `mean` | 평균 제거 |
| Window | `multitaper` | DPSS 기반 |
| Top PSD | OFF | |
| sinc² correction | **ON** | Nyquist 근처 보정 |
| Conversion Method | `sqrt` | √(1+3H) 보정 |
| Hurst H | 0.80 | |
| Correction Factor | 1.8439 | √(1+3×0.8), 자동 계산 |
| Log bins | 88 | Persson 방식 |

---

## 7. 출력

### 7.1 CSV (linear)
- 2열: `q [rad/m]`, `C2D [m⁴]`
- 선형 스케일

### 7.2 Persson (log10)
- 2열: `log10(q)`, `log10(C2D)`
- Persson의 2D.RSG 포맷과 호환

---

## 8. 앙상블 (Ensemble) 모드

여러 프로파일에서 PSD 앙상블을 구성하고 PCA 기반 샘플링 수행.

| 파라미터 | 값 |
|----------|-----|
| 개별 PSD 계산 | 위 §2–§5와 동일한 파이프라인 |
| PCA 분산 임계값 | 90% (0.90) |
| 샘플 수 | 1000 (기본) |

앙상블 psd_params 예시:
```python
psd_params = {
    'detrend': 'mean',
    'window': 'multitaper',
    'use_top_psd': False,
    'sinc2_correct': True,
    'conversion_method': 'sqrt',
    'correction_factor': 1.8439,  # sqrt(1+3*0.8) for H=0.8
    'n_bins': 88,
}
```
