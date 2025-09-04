# %%
"""
EHR Drug–Cancer Re-proportioning Analysis (10-year window, no test split)

Overview
--------
This script builds treated vs. control cohorts from EHR data, balances covariates
via propensity scores (logistic regression + IPTW), and estimates treatment
effects on incident cancer using Cox proportional hazards and Kaplan–Meier curves.

Key steps
---------
1) Candidate cohorts
   - Pull first and recurrent drug exposures per patient (>=30 days apart).
   - Identify treated drug and candidate controls using ATC code similarity.
2) Outcome windowing
   - Define N_YEAR risk window from drug_start_date.
   - Keep patients with sufficient follow-up or cancer within the window.
3) Covariates (confounders)
   - Demographics (age/sex/race/ethnicity) + disease phecodes prior to drug start.
4) Propensity score modeling
   - Logistic regression with GridSearchCV (AUC + weighted balance criterion).
   - Compute IPTW; assess SMD for covariate balance.
5) Survival analysis
   - Weighted CoxPH for HR; PH assumption check.
   - Weighted KM curves with CIs; checkpoint outputs to CSV.
6) Outputs
   - Per treated–control pair: HR, CI, p-value, balance metrics, missingness, follow-up stats.
   - KM curve data (treated / control) and rolling output CSVs with timestamp.

Data expectations (relative paths)
----------------------------------
- ./data/drug_exposure_by_drug/{ingredient_concept_id}.parquet
- ./data/pat_cancer.parquet, ./data/pat_record.parquet, ./data/pat_demo.parquet
- ./data/WHO ATC-DDD 2021-12-03.csv
- ./data/confounding_phecode_detail_notest.csv
- ./data/valid_treat_control_drug_ingredient.csv
- ./data/protein_tables_selected.csv

Notes
-----
- This script assumes sufficiently large cohorts (>=500 per arm) before modeling.
- We cap follow-up at N_YEAR (months) and require baseline of >=12 months.
- Sex-restricted cancers (BRCA, PRAD) are filtered to appropriate sex.
"""

import os
import re
import pandas as pd
from dateutil.relativedelta import relativedelta
import swifter  # used for opportunistic vectorized apply (not explicitly called below)
import openpyxl  # needed if reading/writing Excel elsewhere; safe to keep
import time
from datetime import datetime, date
from tqdm.notebook import tqdm

import random
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.metrics import roc_auc_score

from patsy import dmatrices
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import survival_difference_at_fixed_point_in_time_test, proportional_hazard_test


# %%
def find_similar_atc_codes(atc_name, n_level, atc_data):
    """
    Given an ATC drug name, return other ATC names sharing the same prefix up to a level.

    Parameters
    ----------
    atc_name : str
        Drug name string that appears as 'atc_name' in the WHO ATC table.
    n_level : int
        ATC hierarchy level to match on:
        1 -> first 1 char, 2 -> first 3 chars, 3 -> first 4 chars, 4 -> first 5 chars.
    atc_data : pd.DataFrame
        DataFrame with columns ['atc_code', 'atc_name'].

    Returns
    -------
    list[str]
        Unique list of similar ATC drug names (length-7 ATC codes only), excluding exact match.
    """
    # Map ATC level to prefix length
    atc_length = {1: 1, 2: 3, 3: 4, 4: 5}

    if n_level not in atc_length:
        raise ValueError("Invalid ATC level. Please choose a level between 1 and 4.")

    # Identify codes whose name matches the query (case-insensitive)
    matching_atc_codes = atc_data[atc_data['atc_name'].str.lower() == atc_name.lower()]['atc_code']

    similar_atc_names = []
    for atc_code in matching_atc_codes:
        # Take the hierarchical prefix
        first_chars = atc_code[:atc_length[n_level]]

        # Find peer codes with same prefix (restrict to 7-char codes; exclude the seed code)
        similar_atc_codes_df = atc_data[
            (atc_data['atc_code'].str.startswith(first_chars)) &
            (atc_data['atc_code'].str.len() == 7) &
            (atc_data['atc_code'] != atc_code)
        ]

        # Collect their names
        similar_atc_names.extend(similar_atc_codes_df['atc_name'].tolist())

    return list(set(similar_atc_names))


# %%
class PropensityScoreLR:
    """
    Logistic regression-based propensity score model with grid search and
    a custom scorer that rewards covariate balance (weighted SMD <= 0.1) and AUC.
    """

    def __init__(self, confounder, treatment, random_seed):
        """
        Parameters
        ----------
        confounder : np.ndarray
            Matrix of covariates (after one-hot and merging).
        treatment : np.ndarray
            Binary array (1 = treated, 0 = control).
        random_seed : int
            Seed for CV split determinism and model reproducibility.
        """
        self.confounder = confounder
        self.treatment = treatment
        self.random_seed = random_seed

    @staticmethod
    def truncate(array):
        """Trim extreme IPTW weights to avoid instability."""
        # If using percentile-based trimming, uncomment below lines:
        # lower_value = np.percentile(array, 1)
        # upper_value = np.percentile(array, 99)
        # return np.clip(array, lower_value, upper_value)
        return np.clip(array, a_min=5e-06, a_max=5e1)

    @staticmethod
    def weighted_mean(x, w):
        """Weighted mean along axis=0 for matrix x with weights w (Nx1)."""
        return np.sum(np.multiply(x, w), axis=0) / w.sum()

    @staticmethod
    def weighted_var(x, w):
        """Weighted variance with finite-sample correction."""
        m_w = PropensityScoreLR.weighted_mean(x, w)
        nw, nsw = w.sum(), (w ** 2).sum()
        var = np.multiply((x - m_w) ** 2, w)
        return np.sum(var, axis=0) * (nw / (nw ** 2 - nsw))

    @staticmethod
    def cal_IPTW(y, ps):
        """
        Compute stabilized IPTW weights.

        Parameters
        ----------
        y : np.ndarray of bool/int
            Treatment indicator (True/1 = treated).
        ps : np.ndarray of float
            Propensity scores P(T=1|X).

        Returns
        -------
        treated_w, controlled_w : (N_treated, 1), (N_control, 1)
            Truncated stabilized weights for treated and control.
        """
        ones_idx, zeros_idx = np.where(y == True), np.where(y == False)
        p_T = len(ones_idx[0]) / (len(ones_idx[0]) + len(zeros_idx[0]))
        treated_w = p_T / ps[ones_idx]
        controlled_w = (1 - p_T) / (1. - ps[zeros_idx])
        treated_w = PropensityScoreLR.truncate(treated_w)
        controlled_w = PropensityScoreLR.truncate(controlled_w)
        return np.reshape(treated_w, (len(treated_w), 1)), np.reshape(controlled_w, (len(controlled_w), 1))

    @staticmethod
    def cal_SMD(X, y, propensity_score):
        """
        Compute weighted standardized mean difference (SMD) per covariate.

        SMD = |mean_t - mean_c| / sqrt((var_t + var_c)/2), using IPTW means/vars.
        """
        ones_idx, zeros_idx = np.where(y == True), np.where(y == False)
        treated_w, controlled_w = PropensityScoreLR.cal_IPTW(y, propensity_score)
        treated_X, controlled_X = X[ones_idx], X[zeros_idx]

        treated_X_w_mu = PropensityScoreLR.weighted_mean(treated_X, treated_w)
        controlled_X_w_mu = PropensityScoreLR.weighted_mean(controlled_X, controlled_w)
        treated_X_w_var = PropensityScoreLR.weighted_var(treated_X, treated_w)
        controlled_X_w_var = PropensityScoreLR.weighted_var(controlled_X, controlled_w)

        VAR = np.sqrt((treated_X_w_var + controlled_X_w_var) / 2)
        SMD = np.divide(
            np.abs(treated_X_w_mu - controlled_X_w_mu),
            VAR,
            out=np.zeros_like(treated_X_w_mu),
            where=(VAR != 0)
        )
        return SMD

    def n_weight(self, estimator, X, y):
        """
        Custom scorer passed to GridSearchCV.

        Returns
        -------
        float
            (# of covariates balanced by SMD<=0.1) + AUC on the fold,
            i.e., larger is better (more balanced + good discrimination).
        """
        propensity_score = estimator.predict_proba(self.confounder)[:, 1]
        SMD = self.cal_SMD(self.confounder, self.treatment, propensity_score)
        ps = estimator.predict_proba(X)[:, 1]
        auc = roc_auc_score(y, ps)
        return len(np.where(SMD <= 0.1)[0]) + auc

    def fit_model(self):
        """
        Fit logistic regression PS model using 10-fold CV and the custom scorer.

        Returns
        -------
        best_estimator : sklearn.linear_model.LogisticRegression
        best_params : dict
        """
        model = LogisticRegression()
        cv = KFold(n_splits=10, shuffle=True, random_state=self.random_seed)
        param_grid = {
            'C': [0.05, 0.1, 0.5, 1],
            'penalty': ['l2'],
            'random_state': [self.random_seed],
            'solver': ['liblinear'],
            'max_iter': [1000],
        }
        grid_search = GridSearchCV(
            model,
            param_grid,
            scoring=self.n_weight,
            cv=cv,
            n_jobs=-1
        )
        grid_search.fit(self.confounder, self.treatment)
        return grid_search.best_estimator_, grid_search.best_params_


# %%
def map_race_value(x):
    """Harmonize race categories to WHITE, BLACK, ASIAN, UNKNOWN, OTHER."""
    if x in ['WHITE', 'BLACK', 'ASIAN']:
        return x
    elif x in ['UNKNOWN/TBD', 'DECLINED', 'UNABLE TO PROVIDE']:
        return 'UNKNOWN'
    else:
        return 'OTHER'


def map_ethnicity_value(x):
    """Harmonize ethnicity to HISPANIC/LATINO, NOT HISPANIC/LATINO, or UNKNOWN."""
    if x not in ['HISPANIC/LATINO', 'NOT HISPANIC/LATINO']:
        return 'UNKNOWN'
    else:
        return x


def calculate_time_period(x, column_name_1, column_name_2):
    """
    Month difference between two date columns: column_name_1 - column_name_2.

    Returns np.nan where either date is missing.
    """
    period = (x[column_name_1] - x[column_name_2])
    return period.apply(lambda v: v.days / 30.4375 if not pd.isna(v) else np.nan)


def save_output_checkpoint(output_df, output_file_path):
    """Write (or overwrite) CSV checkpoint to support long runs."""
    output_df.to_csv(output_file_path, index=False)


# %%
def retrieve_drug_record(ingredient_concept_ids):
    """
    For one or more drug ingredient_concept_ids, pull first and recurrent exposure dates.

    Rules
    -----
    - Keep exposures within [2000-01-01, 2024-02-29].
    - For each person: identify first exposure and the first recurrence >= 30 days later.
    - Drop patients without a qualifying recurrence (returns only rows where both dates exist).

    Returns
    -------
    pd.DataFrame with columns:
        person_id, drug_start_date, drug_recur_date
    """
    if not isinstance(ingredient_concept_ids, list):
        ingredient_concept_ids = [ingredient_concept_ids]

    # Read and concatenate per-ingredient exposure logs
    pat_drug = []
    for ingredient_concept_id in ingredient_concept_ids:
        pat_drug.append(
            pd.read_parquet(
                path=f'./data/drug_exposure_by_drug/{ingredient_concept_id}.parquet',
                columns=['person_id', 'drug_exposure_start_datetime']
            )
        )
    pat_drug = pd.concat(pat_drug, ignore_index=True)

    # Normalize to date and drop duplicates
    pat_drug['drug_exposure_start_datetime'] = pd.to_datetime(
        pat_drug['drug_exposure_start_datetime']
    ).dt.date
    pat_drug = pat_drug.drop_duplicates()

    # Sort by patient and exposure date
    pat_drug = pat_drug.sort_values(by=['person_id', 'drug_exposure_start_datetime'])

    # Restrict to analysis window
    pat_drug = pat_drug[
        (pat_drug['drug_exposure_start_datetime'] >= date(2000, 1, 1)) &
        (pat_drug['drug_exposure_start_datetime'] <= date(2024, 2, 29))
    ]

    # Collapse exposure dates to list per person
    pat_drug = pat_drug.groupby('person_id')['drug_exposure_start_datetime'].apply(list).reset_index()

    # Extract first date and first recurrence >= 30 days after
    def get_first_and_recur_dates(dates):
        if len(dates) < 2:
            return [None, None]
        first_date = dates[0]
        for d in dates[1:]:
            if (d - first_date).days >= 30:
                return [first_date, d]
        return [first_date, None]

    pat_drug[['drug_start_date', 'drug_recur_date']] = pd.DataFrame(
        pat_drug['drug_exposure_start_datetime'].apply(get_first_and_recur_dates).tolist(),
        index=pat_drug.index
    )

    # Drop list column and rows missing required dates
    pat_drug = pat_drug.drop(columns=['drug_exposure_start_datetime'])

    return pat_drug.dropna()


# %%
def retrieve_confounder(pat_drug_date, confounder_phecodes):
    """
    Build a person x phecode one-hot matrix using phecodes observed BEFORE drug_start_date.

    Parameters
    ----------
    pat_drug_date : pd.DataFrame
        patient-level DataFrame with ['person_id', 'drug_start_date'].
    confounder_phecodes : list[str]
        Phecode IDs to include.

    Returns
    -------
    pd.DataFrame
        Index=person_id, columns are one-hot phecode indicators (0/1).
    """
    pat_phecode = []

    for confounder_phecode in confounder_phecodes:
        # Load person-level phecode occurrences with dates
        df = pd.read_parquet(
            f'./data/confounder_by_phecode_detail/{confounder_phecode}.parquet',
            columns=['person_id', 'phecode', 'concept_date']
        )

        # Normalize date
        df['concept_date'] = pd.to_datetime(df['concept_date']).dt.date

        # Keep phecodes recorded on/before the patient's drug_start_date
        df = pd.merge(df, pat_drug_date, on='person_id', how='inner')
        df = df[df['drug_start_date'] >= df['concept_date']]

        # Only keep ID columns for one-hot
        pat_phecode.append(df[['person_id', 'phecode']])

    # Combine and one-hot the phecodes
    pat_phecode = pd.concat(pat_phecode, ignore_index=True).drop_duplicates()
    pat_phecode = pd.get_dummies(pat_phecode, columns=['phecode'])

    # Aggregate to person-level indicators
    pat_phecode = pat_phecode.groupby('person_id').agg('sum')

    return pat_phecode


# %%
def retrieve_cancer_pat(cancer_type):
    """Return patients with given cancer_type and their first cancer_start_date."""
    df = pd.read_parquet('./data/pat_cancer.parquet', columns=['person_id', 'cancer', 'cancer_start_date'])
    df['cancer_start_date'] = pd.to_datetime(df['cancer_start_date']).dt.date
    return df.loc[df['cancer'] == cancer_type, ['person_id', 'cancer_start_date']]


def retrieve_record_span(person_ids):
    """Load each person's min/max record dates, restricted to supplied person_ids."""
    df = pd.read_parquet('./data/pat_record.parquet')
    df['max_record_date'] = pd.to_datetime(df['max_record_date']).dt.date
    df['min_record_date'] = pd.to_datetime(df['min_record_date']).dt.date
    return df.loc[df['person_id'].isin(person_ids), :]


def retrieve_demo(person_ids):
    """Load demographics for supplied person_ids."""
    df = pd.read_parquet('./data/pat_demo.parquet')
    df['birth_datetime'] = pd.to_datetime(df['birth_datetime']).dt.date
    return df.loc[df['person_id'].isin(person_ids), :]


# %%
def calculate_missing_proportion(pat, pat_demo):
    """
    Compute column-wise missingness (%) for gender/race/ethnicity and
    percent of persons with any missing value.
    """
    df = pd.merge(
        pat,
        pat_demo[['person_id', 'gender', 'race', 'ethnicity']],
        on='person_id',
        how='left'
    )
    df = df[['person_id', 'gender', 'race', 'ethnicity']]

    missing_value = df.isna().mean().to_dict()
    missing_value = {key: value * 100 for key, value in missing_value.items()}

    missing_person = df.isna().any(axis=1).mean() * 100
    return missing_value, missing_person


# %%
def check_linearity(df, numeric_var, response_var):
    """
    Box–Tidwell style check for logistic-linearity of a numeric predictor.

    Adds log(numeric_var) and interaction numeric*log(numeric); if the interaction
    is significant (p<0.05), flags violation.
    """
    # Avoid log(0) by adding a tiny constant
    df[f'log_{numeric_var}'] = np.log(df[numeric_var] + 1e-6)

    formula = f'{response_var} ~ {numeric_var} + log_{numeric_var} + {numeric_var}:log_{numeric_var}'
    y, X = dmatrices(formula, data=df, return_type='dataframe')
    try:
        model = sm.Logit(y, X).fit()
        if model.pvalues.get(f'{numeric_var}:log_{numeric_var}', 1) < 0.05:
            return 'violate'
    except np.linalg.LinAlgError:
        # Singular matrix etc. — treat as inconclusive
        pass


def check_proportional(cph, cox_data):
    """
    Global PH check via lifelines' proportional_hazard_test (rank transform).

    Returns 'violate' if any covariate shows p<0.05.
    """
    results = proportional_hazard_test(cph, cox_data, time_transform='rank')
    if any(results.summary['p'] < 0.05):
        return 'violate'


# %%
# Timestamped file paths to support long jobs and partial results
date_time_suffix = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
output_file_path = f'./output_{date_time_suffix}.csv'
km_curve_treated_file_path = f"./KM_curve_treated_{date_time_suffix}.csv"
km_treated_df = pd.DataFrame()
km_curve_controlled_file_path = f"./KM_curve_controlled_{date_time_suffix}.csv"
km_controlled_df = pd.DataFrame()
output_df = pd.DataFrame()

# Metadata and crosswalk tables
atc_df = pd.read_csv('./data/WHO ATC-DDD 2021-12-03.csv')  # must include atc_code, atc_name
confounder_table = pd.read_csv('./data/confounding_phecode.csv', dtype={'Phecode': 'string'})
drug_ingredient = pd.read_csv('./data/valid_treat_control_drug_ingredient.csv')
drug_table = pd.read_csv('./data/protein_tables_selected.csv')

# ATC proximity level and random seed for reproducibility
n_level = 2
random_seed = 341

# Split drug_ingredient table into treated vs. control pools
treated_drug_ingredient = drug_ingredient.loc[drug_ingredient.type == 'treat', ['MolecularBrandName', 'ingredient_concept_id']]
controlled_drug_ingredient = drug_ingredient.loc[drug_ingredient.type == 'control', ['MolecularBrandName', 'concept_name', 'ingredient_concept_id']]

# %% [markdown]
# ## N-year
# Analysis window for incident cancer after drug_start_date

# %%
N_YEAR = 10  # follow-up cap / risk window in years

# Iterate over each treated drug row in the protein table
for index, row in tqdm(drug_table.iterrows(), total=drug_table.shape[0], desc='treated loop'):
    treated_drug = row['MolecularBrandName']

    # Ingredient(s) for the treated drug
    treated_drug_concept = treated_drug_ingredient.loc[
        treated_drug_ingredient['MolecularBrandName'] == treated_drug, 'ingredient_concept_id'
    ].to_list()

    # Patients with first exposure and qualifying recurrence for treated drug
    pat_treat = retrieve_drug_record(treated_drug_concept)

    # Identify candidate control drugs by ATC similarity to the first token
    treated_drug_first = treated_drug.split(' ')[0]
    candidate_controlled_drugs = find_similar_atc_codes(treated_drug_first, n_level, atc_df)

    # If the treated drug name has multiple tokens, search again on full name
    if treated_drug_first != treated_drug:
        add_candidate_controlled_drugs = find_similar_atc_codes(treated_drug, n_level, atc_df)
        candidate_controlled_drugs.extend(add_candidate_controlled_drugs)

    # Intersect ATC similar names with control pool available in drug_ingredient
    candidate_controlled_drugs = list(set(
        controlled_drug_ingredient.loc[
            controlled_drug_ingredient['MolecularBrandName'].isin(candidate_controlled_drugs),
            'MolecularBrandName'
        ]
    ))

    # Scan each candidate control and run the full pipeline
    for controlled_drug in tqdm(candidate_controlled_drugs, desc='controlled loop', leave=False):
        controlled_drug_concept = controlled_drug_ingredient.loc[
            controlled_drug_ingredient['MolecularBrandName'] == controlled_drug, 'ingredient_concept_id'
        ].to_list()

        # Pull control cohort analogously
        pat_control = retrieve_drug_record(controlled_drug_concept)

        # Combine treated + control; drop overlap (patients exposed to both)
        pat = pd.concat([pat_treat, pat_control], ignore_index=True)
        pat['is_treated'] = pat.person_id.isin(pat_treat.person_id)
        pat['is_controlled'] = pat.person_id.isin(pat_control.person_id)
        pat = pat[~(pat['is_treated'] & pat['is_controlled'])]

        # Ensure sufficient sample size before modeling
        if (pat['is_treated'].sum() < 500) or (pat['is_controlled'].sum() < 500):
            continue

        # Label outcome within N_YEAR window
        pat_cancer = retrieve_cancer_pat(row['Cancer'])
        pat = pd.merge(pat, pat_cancer, on='person_id', how='left')
        pat['is_cancer'] = (
            (pat['cancer_start_date'].notna()) &
            (pat['cancer_start_date'] <= pat['drug_start_date'] + relativedelta(years=N_YEAR))
        )

        # Attach record span and enforce follow-up rules
        pat_record = retrieve_record_span(pat['person_id'])
        pat = pd.merge(pat, pat_record, on='person_id', how='inner')
        # Cap max_record_date at N_YEAR after drug start
        pat['max_record_date'] = pat.apply(
            lambda x: min(x["drug_start_date"] + relativedelta(years=N_YEAR), x["max_record_date"]),
            axis=1
        )

        # Eligibility:
        # - If cancer occurs, require it to be >=3 months after drug_recur_date (avoid immediate protopathic bias).
        # - If no cancer, require follow-up through drug_recur_date (i.e., recurrence observed).
        pat = pat[
            (
                (pat['is_cancer']) &
                (pat['cancer_start_date'] >= pat['drug_recur_date'] + relativedelta(months=3))
            ) |
            (
                (~pat['is_cancer']) &
                (pat["max_record_date"] >= pat["drug_recur_date"])
            )
        ]

        # Add birth date (for age calculations)
        pat_demo = retrieve_demo(pat.person_id)
        pat = pd.merge(pat, pat_demo[['person_id', 'birth_datetime']], on='person_id', how='inner')

        # --- Time windows / ages (in months) ---
        pat['baseline_period'] = calculate_time_period(pat, 'drug_start_date', 'min_record_date')   # months before drug start
        pat['follow_up_period'] = calculate_time_period(pat, 'max_record_date', 'drug_start_date') # months at risk
        pat['time_to_event'] = calculate_time_period(pat, 'cancer_start_date', 'drug_start_date')  # months to cancer
        pat['time_to_event'] = pat['time_to_event'].fillna(pat['follow_up_period'])

        # (1) Positive time; (2) cap at N_YEAR*12 months; (3) >=12 months baseline
        pat = pat.loc[(pat['time_to_event'] > 0), :]
        pat['time_to_event'] = pat['time_to_event'].apply(lambda x: 12 * N_YEAR if x > 12 * N_YEAR else x)
        pat = pat.loc[(pat['baseline_period'] >= 12), :]

        # Ages at drug start, at cancer, and at censoring; cohort_age = age at event/censor
        pat['age'] = calculate_time_period(pat, 'drug_start_date', 'birth_datetime')
        pat['cancer_age'] = calculate_time_period(pat, 'cancer_start_date', 'birth_datetime')
        pat['record_age'] = calculate_time_period(pat, 'max_record_date', 'birth_datetime')
        pat['cohort_age'] = pat.apply(lambda x: x.cancer_age if x.is_cancer else x.record_age, axis=1)

        # Restrict to adults >=40 years at event/censor
        pat = pat[pat['cohort_age'] >= 40 * 12]

        # Missingness summary before merging additional demographic columns
        missing_value, missing_person = calculate_missing_proportion(pat, pat_demo)

        # Demographics (gender/race/ethnicity)
        pat = pd.merge(pat, pat_demo[['person_id', 'gender', 'race', 'ethnicity']], on='person_id', how='inner')

        # Encode gender: male=1, female=0
        pat['gender'] = pat['gender'].apply(lambda x: int(x == 'MALE'))

        # Sex-specific cancers
        if row['Cancer'] == 'BRCA':
            pat = pat[pat['gender'] == 0]  # females only
        elif row['Cancer'] == 'PRAD':
            pat = pat[pat['gender'] == 1]  # males only

        # Harmonize race/ethnicity and one-hot
        pat['race'] = pat.race.apply(map_race_value)
        pat['ethnicity'] = pat.ethnicity.apply(map_ethnicity_value)

        pat_race = pd.get_dummies(pat[['person_id', 'race']], columns=['race'], prefix='race')
        pat_race = pat_race[['person_id', 'race_BLACK', 'race_WHITE', 'race_OTHER']]
        pat_ethnicity = pd.get_dummies(pat[['person_id', 'ethnicity']], columns=['ethnicity'], prefix='ethnicity')
        pat_ethnicity = pat_ethnicity[['person_id', 'ethnicity_HISPANIC/LATINO', 'ethnicity_NOT HISPANIC/LATINO']]

        # Pick phecodes: cancer-specific + general (Category 'G')
        confounder_phecodes = confounder_table.loc[
            confounder_table['Category'].isin([row['Cancer'], 'G']),
            'Phecode'
        ].to_list()
        pat_phecode = retrieve_confounder(pat[['person_id', 'drug_start_date']], confounder_phecodes)

        # Build analysis table
        if (row['Cancer'] == 'BRCA') or (row['Cancer'] == 'PRAD'):
            data = pat[['person_id', 'is_treated', 'is_cancer', 'time_to_event', 'age']]
        else:
            data = pat[['person_id', 'is_treated', 'is_cancer', 'time_to_event', 'age', 'gender']]

        data = pd.merge(data, pat_race, how='inner', on='person_id')
        data = pd.merge(data, pat_ethnicity, how='inner', on='person_id')
        data = pd.merge(data, pat_phecode.reset_index(), how='left', on='person_id')

        # Missing phecodes -> 0
        data[pat_phecode.columns] = data[pat_phecode.columns].fillna(0)

        # Re-check minimum arm sizes after merges
        if (data['is_treated'].sum() < 500) or (len(data) - data['is_treated'].sum() < 500):
            continue

        # Scale age to [0,1] for stability in LR
        data['age'] = (data['age'] - data['age'].min()) / (data['age'].max() - data['age'].min())

        # Ensure boolean columns are numeric 0/1
        for col in data.columns:
            if data[col].dtype == 'bool':
                data[col] = data[col].astype(int)

        # Evaluate linearity assumption of age for PS model (advisory)
        lr_assumption = check_linearity(data, 'age', 'is_treated')

        # --- Propensity score model ---
        confounder = data.iloc[:, 4:].to_numpy().astype(float)  # all covariates after [person_id, is_treated, is_cancer, time_to_event]
        treatment = data['is_treated'].to_numpy().astype(float)

        lr = PropensityScoreLR(confounder, treatment, random_seed)
        best_model, best_params = lr.fit_model()
        propensity_score = best_model.predict_proba(confounder)[:, 1]

        # SMD and proportion of unbalanced covariates (>0.1)
        smd = PropensityScoreLR.cal_SMD(confounder, treatment, propensity_score)
        p_unbalanced = len(np.where(smd > 0.1)[0]) / len(smd)

        # Split indices and get IPTW
        ones_idx, zeros_idx = np.where(treatment == True), np.where(treatment == False)
        treated_w, controlled_w = PropensityScoreLR.cal_IPTW(treatment, propensity_score)
        combined_w = np.concatenate((treated_w, controlled_w))

        # Durations/censoring for KM/CoxPH (months)
        duration = data['time_to_event'].to_numpy()
        treated_duration, controlled_duration = duration[ones_idx], duration[zeros_idx]

        # Approximate overall censoring proportion (uses control durations as reference)
        cencor_prop = sum(1 for x in controlled_duration if x < max(controlled_duration)) / len(duration)

        cph = CoxPHFitter()

        cancer = data['is_cancer'].to_numpy()
        treated_cancer, controlled_cancer = cancer[ones_idx], cancer[zeros_idx]

        # Weighted event rates (treated vs. control) — for diagnostics
        treated_prob = np.matmul(treated_cancer, treated_w) / np.sum(treated_w)
        controlled_prob = np.matmul(controlled_cancer, controlled_w) / np.sum(controlled_w)

        # --- Cox proportional hazards model (weighted by IPTW) ---
        weight = np.zeros(len(cancer))
        weight[ones_idx] = treated_w.squeeze()
        weight[zeros_idx] = controlled_w.squeeze()

        # Start with treatment + any covariates with SMD>0.1 (i.e., adjust residual imbalance)
        cox_data = pd.DataFrame({'T': duration, 'event': cancer, 'treatment': treatment, 'weights': weight})
        cox_data = pd.concat([cox_data, data.iloc[:, 4:].iloc[:, np.where(smd > 0.1)[0]]], axis=1)

        try:
            cph.fit(cox_data, duration_col='T', event_col='event', weights_col='weights', robust=True)
            cox_assumption = check_proportional(cph, cox_data)
        except Exception as e:
            # If the model is unstable (e.g., collinearity), compute VIF and drop high-VIF covariates (>10)
            try:
                vif = pd.DataFrame()
                vif["Feature"] = cox_data.columns
                vif["VIF"] = [variance_inflation_factor(cox_data, i) for i in range(cox_data.shape[1])]

                high_vif_features = vif[vif["VIF"] > 10]["Feature"]
                high_vif_features = [ft for ft in high_vif_features if ft not in ['T', 'event', 'weights']]
                cox_data_filtered = cox_data.drop(columns=high_vif_features)

                cph.fit(cox_data_filtered, duration_col='T', event_col='event', weights_col='weights', robust=True)
                cox_assumption = check_proportional(cph, cox_data_filtered)
            except Exception as e:
                # If model still fails, skip this control and continue with next
                continue

        # Re-run PH check on full data (as in original code) — potentially redundant but kept for parity
        cox_assumption = check_proportional(cph, cox_data)

        # --- Kaplan–Meier curves (weighted) ---
        kmf_A = KaplanMeierFitter()
        kmf_B = KaplanMeierFitter()

        treated_kmf = kmf_A.fit(treated_duration, treated_cancer, label="Treated", weights=treated_w)
        controlled_kmf = kmf_B.fit(controlled_duration, controlled_cancer, label="Controlled", weights=controlled_w)

        # Extract survival functions with confidence intervals
        treated_survival_df = treated_kmf.survival_function_.rename(columns={'KM_estimate': 'KM_estimate_treated'})
        controlled_survival_df = controlled_kmf.survival_function_.rename(columns={'KM_estimate': 'KM_estimate_controlled'})

        treated_survival_df = treated_survival_df.join(treated_kmf.confidence_interval_)
        controlled_survival_df = controlled_survival_df.join(controlled_kmf.confidence_interval_)

        # Append to in-memory buffers, then checkpoint to CSV
        km_treated_df = pd.concat([km_treated_df, treated_survival_df.reset_index()], ignore_index=True)
        km_controlled_df = pd.concat([km_controlled_df, controlled_survival_df.reset_index()], ignore_index=True)

        save_output_checkpoint(km_treated_df, km_curve_treated_file_path)
        save_output_checkpoint(km_controlled_df, km_curve_controlled_file_path)

        # --- Aggregate outputs for this treated-control pair ---
        new_row = pd.Series({
            'Cancer': row['Cancer'],
            'treated drug': treated_drug,
            'controlled drug': controlled_drug,
            'unbalanced_covariate_percentage': p_unbalanced,
            'HR': cph.hazard_ratios_['treatment'],
            'CI lower': cph.summary.loc['treatment', 'exp(coef) lower 95%'],
            'CI upper': cph.summary.loc['treatment', 'exp(coef) upper 95%'],
            'P value': cph.summary.loc['treatment', 'p'],
            'treated group size': data['is_treated'].sum(),
            'controlled group size': len(data) - data['is_treated'].sum(),
            'treated cancer group size': treated_cancer.sum(),
            'controlled cancer group size': controlled_cancer.sum(),
            'mean IPTW weight': np.mean(combined_w),
            'std IPTW weight': np.std(combined_w),
            'missing proportion by column': missing_value,
            'missing proportion by person': missing_person,
            'median follow up': np.median(duration),
            'IQR follow up': (np.percentile(duration, 25), np.percentile(duration, 75)),
            'cencor proportion': cencor_prop,
            'treatment patient year': np.sum(duration) / 12,
            'proportional hazards assumption': cox_assumption,
            'linearity assumption': lr_assumption,
        })

        output_df = pd.concat([output_df, new_row.to_frame().T], ignore_index=True)
        save_output_checkpoint(output_df, output_file_path)
