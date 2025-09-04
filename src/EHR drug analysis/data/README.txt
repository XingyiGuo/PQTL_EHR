Sample EHR inputs for the drug–cancer analysis pipeline (documentation set)

Note: The production pipeline expects several files as Parquet.
These CSVs are intentionally small and are here to *illustrate schema and formatting*.
If you want to run the full script against these, convert them to Parquet first, e.g.:

    import pandas as pd
    df = pd.read_csv("pat_demo.csv", parse_dates=["birth_datetime"])
    df.to_parquet("pat_demo.parquet")

Directory layout:
- WHO ATC-DDD 2021-12-03.csv           (ATC code/name reference)
- valid_treat_control_drug_ingredient.csv
- protein_tables_selected.csv
- pat_cancer.csv                        (person_id, cancer, cancer_start_date)
- pat_record.csv                        (person_id, min_record_date, max_record_date)
- pat_demo.csv                          (person_id, birth_datetime, gender, race, ethnicity)
- drug_exposure_by_drug/*.csv           (person_id, drug_exposure_start_datetime) per ingredient_concept_id
- confounding_phecode_detail_notest.csv (Phecode, Category)
- confounder_by_phecode_detail/*.csv    (person_id, phecode, concept_date) per phecode
