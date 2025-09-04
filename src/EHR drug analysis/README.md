# EHR Drug–Cancer Re-Purposing 

> Survival analysis of incident cancer following drug exposure using EHR data.  
> Cohorts are built from first and recurrent exposures, covariates are balanced with propensity scores + IPTW, and hazards are estimated via Cox PH with weighted KM curves.


---

## TL;DR (quick start)

```bash
# 1) Create a clean environment
conda create -n ehr-cancer python=3.10 -y
conda activate ehr-cancer

# 2) Install dependencies
pip install -r requirements.txt
# If you plan to use Parquet (recommended):
pip install pyarrow fastparquet

# 3) Run the main analysis script
python iptw_drug_analysis.py
```

Outputs will land alongside the script as timestamped CSVs, e.g.:

- `output_YYYY-MM-DD_HH-MM-SS.csv`
- `KM_curve_treated_YYYY-MM-DD_HH-MM-SS.csv`
- `KM_curve_controlled_YYYY-MM-DD_HH-MM-SS.csv`

---

## Repository layout

```
.
├── iptw_drug_analysis.py                        # main pipeline
├── requirements.txt                             # pinned deps (see below)
├── data/                                        # input data (see schemas)
│   ├── WHO ATC-DDD 2021-12-03.csv
│   ├── valid_treat_control_drug_ingredient.csv
│   ├── drug_table.csv
│   ├── pat_cancer.parquet | .csv
│   ├── pat_record.parquet | .csv
│   ├── pat_demo.parquet   | .csv
│   ├── drug_exposure_by_drug/
│   │   ├── 1001.parquet | .csv
│   │   ├── 1002.parquet | .csv
│   │   └── ...
│   └── confounder_by_phecode_detail/
│       ├── 250.2.parquet | .csv
│       ├── 401.1.parquet | .csv
│       └── ...
└── scripts/
    └── csv_to_parquet.py                        # helper to convert CSV → Parquet
```

> The main script expects **Parquet**. CSVs are fine for explaining the schema; convert them before running (helper below).

---

## Data contracts (schemas)

### 1) ATC reference

**`data/WHO ATC-DDD 2021-12-03.csv`**

| column     | type   | notes                       |
|------------|--------|-----------------------------|
| atc_code   | str    | 7-char ATC code (e.g., C10AA05) |
| atc_name   | str    | human name (e.g., Atorvastatin) |

Used to find **ATC-similar** controls by prefix (level 1–4).

---

### 2) Treated/Control catalog

**`data/valid_treat_control_drug_ingredient.csv`**

| column                | type | notes                                               |
|----------------------|------|-----------------------------------------------------|
| MolecularBrandName   | str  | drug label used across tables                       |
| concept_name         | str  | synonym; not required by code logic                 |
| ingredient_concept_id| int  | numeric key used to locate exposure files           |
| type                 | str  | `"treat"` or `"control"`                            |

---

### 3) Analysis pairs (what to run)

**`data/drug_table.csv`**

| column              | type | notes                  |
|---------------------|------|------------------------|
| MolecularBrandName  | str  | treated drug to test   |
| Cancer              | str  | target cancer code     |

> The script loops rows here: for each treated drug, it finds ATC-similar controls and runs the pipeline.

---

### 4) Drug exposure logs (per ingredient)

**`data/drug_exposure_by_drug/{ingredient_concept_id}.parquet`**

| column                         | type      | notes                                       |
|--------------------------------|-----------|---------------------------------------------|
| person_id                      | int       | patient key                                 |
| drug_exposure_start_datetime   | timestamp | exposure start; daily resolution is fine     |

Rules applied by the script:

- Keep exposures in **[2000-01-01, 2024-02-29]**.  
- For each person, take **first exposure** and the **first recurrence ≥ 30 days** later.  
- Patients without a qualifying recurrence are excluded.

---

### 5) Cancer registry

**`data/pat_cancer.parquet`**

| column            | type      | notes                              |
|-------------------|-----------|------------------------------------|
| person_id         | int       |                                    |
| cancer            | str       | must match the `Cancer` values     |
| cancer_start_date | timestamp | first observed cancer date         |

---

### 6) Record span

**`data/pat_record.parquet`**

| column          | type      | notes                                  |
|-----------------|-----------|----------------------------------------|
| person_id       | int       |                                        |
| min_record_date | timestamp | first seen in system                   |
| max_record_date | timestamp | last seen; capped at N_YEAR in script  |

---

### 7) Demographics

**`data/pat_demo.parquet`**

| column          | type      | notes                                                     |
|-----------------|-----------|-----------------------------------------------------------|
| person_id       | int       |                                                           |
| birth_datetime  | timestamp | used for age; age ≥ 40 at event/censor required           |
| gender          | str       | `"MALE"`/`"FEMALE"`; encoded to 1/0                       |
| race            | str       | harmonized to `WHITE`, `BLACK`, `ASIAN`, `UNKNOWN`, `OTHER` |
| ethnicity       | str       | `HISPANIC/LATINO`, `NOT HISPANIC/LATINO`, or `UNKNOWN`   |

---

### 8) Confounder catalog

**`data/confounding_phecode.csv`**

| column   | type | notes                                        |
|----------|------|----------------------------------------------|
| Phecode  | str  | like `"401.1"`                               |
| Category | str  | `"G"` (general) or cancer-specific (e.g., `"LUNG"`) |

The script selects **general + cancer-specific** phecodes per run.

---

### 9) Confounder detail (per phecode)

**`data/confounder_by_phecode_detail/{phecode}.parquet`**

| column       | type      | notes                                        |
|--------------|-----------|----------------------------------------------|
| person_id    | int       |                                              |
| phecode      | str       |                                              |
| concept_date | timestamp | only occurrences **on/before drug_start_date** count |

The script builds a person×phecode 0/1 one-hot matrix from these files.

---

## Study design (implemented in code)

- **Window:** fixed **N_YEAR = 10** years after **drug_start_date**.  
- **Eligibility:**  
  - Must observe a **recurrence ≥ 30 days** after first exposure.  
  - If cancer occurs, it must be **≥ 3 months after recurrence** (reduces protopathic bias).  
  - If no cancer, follow-up must extend **through recurrence**.  
  - **Baseline ≥ 12 months** prior to drug start.  
  - **Age ≥ 40** at event/censor.  
  - Sex-restricted cancers (e.g., BRCA/PRAD) are filtered appropriately.

- **Treatment assignment:** treated drug vs ATC-similar controls (ATC prefix match at level **n_level=2** by default).

- **Confounding control:** Logistic regression **propensity score** with **GridSearchCV**; custom scoring rewards **AUC** and **weighted balance** (share of covariates with SMD ≤ 0.1). Stabilized **IPTW** with truncation.

- **Outcomes:**  
  - Weighted **Cox PH** for HR and 95% CI (robust).  
  - Weighted **KM** curves with CIs for treated and control.  
  - **PH assumption** via `proportional_hazard_test`; **linearity** check for age via a Box–Tidwell style term.

---

## Outputs (what to expect)

### Primary table (`output_*.csv`)

| column                              | description                                             |
|-------------------------------------|---------------------------------------------------------|
| Cancer                               | target cancer (from input row)                         |
| treated drug / controlled drug       | the pair evaluated                                     |
| HR, CI lower, CI upper, P value      | hazard ratio (treatment effect)                        |
| unbalanced_covariate_percentage      | share of covariates with SMD > 0.1 after weighting     |
| treated / controlled group size      | N in each arm after all filters                        |
| treated / controlled cancer group size| events per arm                                         |
| mean / std IPTW weight               | weight diagnostics                                     |
| missing proportion by column/person  | demographic missingness                                 |
| median follow up, IQR follow up      | months                                                 |
| cencor proportion                    | crude censoring proxy                                  |
| treatment patient year               | person-months/12                                       |
| proportional hazards assumption      | “violate” if any covariate fails                       |
| linearity assumption                 | “violate” if age violates logistic linearity           |

### KM curve data

- `KM_curve_treated_*.csv` and `KM_curve_controlled_*.csv`  
  Survival estimates and confidence intervals over time (months).

---

## Reproducibility notes

- Seed: `random_seed = 341` (set in the script; used for CV splits).  
- Weight truncation: `np.clip(weights, 5e-06, 5e1)`.  
- Balance metric: SMD with IPTW means/variances; threshold 0.1.

---

## Dependencies

Minimal set (pin in `requirements.txt` as you prefer):

```
pandas>=2.0
numpy>=1.24
scikit-learn>=1.2
lifelines>=0.27
statsmodels>=0.14
patsy>=0.5
python-dateutil>=2.8
tqdm>=4.66
pyarrow>=15.0         # for Parquet IO
openpyxl>=3.1         # only needed if you export Excel elsewhere
swifter>=1.4          # optional; safe if present
```

> If you don’t want Parquet, you can change the readers in the script from `read_parquet` to `read_csv` and pass `parse_dates` appropriately. Parquet is much faster for large EHR tables.

---

## Converting the sample CSVs to Parquet

Run:

```bash
python scripts/csv_to_parquet.py
```

It will write .parquet files alongside the CSVs in your data/ folder.

---

## How to run on your data

1. Populate the **data/** folder in the same structure and schema.  
2. Ensure `drug_table.csv` lists your **treated drugs** and **cancer types** of interest.  
3. Make sure each treated drug appears in `valid_treat_control_drug_ingredient.csv` with `type="treat"` and that reasonable ATC-similar controls exist with `type="control"`.  
4. Convert to Parquet (or toggle the readers in the script to CSV).  
5. `python iptw_drug_analysis.py`

---

## Interpreting results

- **HR < 1** suggests lower incident cancer hazard in the treated group (vs. ATC-matched controls) within the 10-year window.  
- Trust the result only if:
  - Unbalanced covariate share is **low** after weighting (aim: ≤ 10–20%).  
  - **PH assumption** is not violated.  
  - Event counts are sufficient in both arms.  
  - Follow-up distribution is reasonable (no one-sided truncation).  

If those checks fail, treat the signal as exploratory.

---

## Troubleshooting

- **Model fails / singular matrix**  
  The script computes VIF and drops high-VIF covariates (>10) before refitting. If it still fails, the pair is skipped; look at cohort sizes and collinearity.

- **Empty or tiny cohorts**  
  The pipeline requires **≥ 500 per arm** after all filters. Relax this locally if you’re debugging, but keep in mind the variance will explode.

- **ATC control list is empty**  
  Verify the ATC name in `WHO ATC-DDD ...csv` matches the drug labels in your catalogs (case-insensitive match is used, but the tokenization matters).

- **Parquet errors**  
  Install `pyarrow` and ensure date columns are true timestamps. If you must use CSVs, switch `read_parquet` → `read_csv` and add `parse_dates`.

---

## Re-use and adaptation

- Swap `n_level` (ATC prefix depth) to tighten or loosen control similarity.  
- Adjust `N_YEAR` to other horizons (e.g., 5-year sensitivity).  
- Consider doubly-robust estimation (e.g., AIPW) if you extend the pipeline.

---
