# 코드 공개 — GitHub 저장소 + Zenodo DOI

세 리뷰어가 모두 요구한 필수 항목입니다. **심사 중에 공개되어 있어야** 하므로
"게재 시 공개"가 아니라 지금 올리셔야 합니다. 전 과정 15분.

**업로드할 폴더는 이미 준비되어 있습니다:**
`CYB5R3.../cyb5r3-marc1-lihc-upload/`
(원본 작업폴더에서 `cache/`와 임시파일을 뺀 깨끗한 사본, 3.2 MB, 92개 파일 —
원고용 그림 PNG와 각 그림의 수치 원자료 CSV 포함. TIFF는 용량 때문에 제외했고
MDPI 에는 별도 제출합니다.)

---

## 1단계 · GitHub 저장소 만들기 (5분)

1. https://github.com/new 접속
2. **Repository name:** `cyb5r3-marc1-lihc`
3. **Public** 선택 — Private이면 리뷰어가 볼 수 없습니다
4. "Add a README file", ".gitignore", "license" 는 **모두 체크 해제** (이미 들어 있습니다)
5. **Create repository**
6. 다음 화면 중간의 **"uploading an existing file"** 링크 클릭
7. Finder에서 `cyb5r3-marc1-lihc-upload` 폴더를 열고, **폴더 자체가 아니라 그 안의
   항목 전부**(R, docs, results, figures, README.md, RESULTS_SUMMARY.md, LICENSE,
   CITATION.cff, RUN_*.R)를 선택해 브라우저 창으로 드래그
8. 아래 Commit message 에 `Analysis pipeline for IJMS ijms-4481300 revision 1` 입력
9. **Commit changes**

업로드 후 저장소 첫 화면에 README가 표로 보이면 성공입니다.

---

## 2단계 · Zenodo 연동 (3분)

1. https://zenodo.org 접속 → 우상단 **Log in** → **Log in with GitHub**
2. 로그인 후 우상단 계정 메뉴 → **GitHub**
3. 저장소 목록에서 `cyb5r3-marc1-lihc` 옆 스위치를 **ON**

**순서가 중요합니다.** 스위치를 켠 *뒤에* release를 만들어야 Zenodo가 잡아냅니다.

---

## 3단계 · Release 발행 → DOI 자동 발급 (3분)

1. GitHub 저장소 → 우측 **Releases** → **Create a new release**
2. **Choose a tag** → `v1.0.0` 입력 → "Create new tag on publish"
3. **Release title:** `IJMS ijms-4481300 revision 1 — analysis pipeline`
4. **Description:** 아래를 붙여넣기
5. **Publish release**
6. 1~2분 뒤 Zenodo Uploads 목록에 나타나고 **DOI**가 부여됩니다 (`10.5281/zenodo.XXXXXXX`)

```
Reproducible R pipeline for "Divergence of the CYB5R3–mARC1 redox axis across the
human MASLD-to-hepatocellular-carcinoma continuum: an exploratory analysis of
liver-biopsy and tumour cohorts" (International Journal of Molecular Sciences,
manuscript ijms-4481300, revision 1).

Contents: analysis scripts (R/), the numerical source data underlying every figure
and table (results/, figures/*_source.csv), figures at 300 dpi (figures/), run logs,
package versions and session information.

Data sources: UCSC Xena (TCGA, GTEx, PanCanAtlas), GEO (GSE14520, GSE130970,
GSE135251, GSE162694) and DepMap. No raw patient data are redistributed; all inputs
are downloaded by the scripts from their public repositories.

Licence: MIT (code).
```

---

## 4단계 · 원고에 반영 (2분)

Zenodo 기록 페이지의 **"Cite all versions"** 아래에 있는 **concept DOI**를 쓰십시오
(버전이 올라가도 안 바뀝니다).

원고에서 바꿀 곳은 **두 군데뿐**입니다.

**Data Availability Statement:**
> …is available at **https://github.com/drswnam69-debug/cyb5r3-marc1-lihc** and archived
> with the citable DOI **⟨Zenodo DOI⟩**.

→ GitHub 아이디(`drswnam69-debug`)는 이미 채워 넣었습니다. 남은 것은 `⟨Zenodo DOI⟩` 하나로,
   Zenodo 가 발급한 `10.5281/zenodo.XXXXXXX` 로 바꾸면 됩니다.

**응답서**에도 같은 URL·DOI를 한 번 적어야 리뷰어가 확인할 수 있습니다.
`CITATION.cff`의 `repository-code` 줄도 함께 고쳐 두면 깔끔합니다.

---

## 흔한 문제

| 증상 | 해결 |
|---|---|
| Zenodo에 저장소가 안 보임 | 스위치를 켠 **뒤에** release를 만들어야 합니다. 순서 확인 |
| DOI가 안 나옴 | Zenodo → GitHub 페이지에서 **Sync now** |
| 저장소가 Private | Settings → General → 맨 아래 Change visibility → Public |
| 드래그가 안 됨 | 폴더째가 아니라 **안의 항목들**을 선택해서 드래그 |

## 이후 수정이 생기면

새 release(`v1.0.1`)를 만들면 Zenodo가 자동으로 새 버전을 만듭니다.
원고에는 concept DOI를 적어 두었으므로 고칠 필요가 없습니다.
