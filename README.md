# Causal Inference of Drug Effects on Metabolomics and Proteomics in UK Biobank

## Description
This repository contains the downstream analysis for the paper  "﻿Causal Inference of Drug Effects on Metabolomics and Proteomics in UK Biobank". 

## Content
The repository is organized into the following directories:
- [0data_processing](0data_processing): Scripts to aggregate covariates and processing omics biomarkers in the [UK Biobank](https://www.ukbiobank.ac.uk/) for use in different blood biomarkers.
  - [rdit_data_processing](0data_processing/rdit_data_processing.py): Generates data for RDiT analysis.
  - [rdit_data_intervention_processing](0data_processing/rdit_data_intervention_processing.py): Prepares data for RDiT analysis of randomized interventions.
  - [cross_sectional_association_data_processing](0data_processing/cross_sectional_association_data_processing.R): Prepares data for cross-sectional association analysis.
  - [longitudinal_association_data_processing](0data_processing/longitudinal_association_data_processing.R): Prepares data for longitudinal association analysis.
  - [prepare_dataset_liver_muscle_damage](0data_processing/prepare_dataset_liver_muscle_damage.ipynb): Prepares dataset for associations of statin-influenced proteins with liver and muscle damage.


- ([1run_models](1run_models)): Scripts to reproduce the analyses and evaluate the causal associations between the ATC Level 4 Prescription and each omics biomarker.
  - [rdit_analysis](1run_models/rdit_analysis.py): Executes RDiT analysis.
  - [rdit_power_analysis](1run_models/rdit_power_analysis.py): Runs power analysis for the RDiT model.
  - [run_intervention_analysis](1run_models/run_intervention_analysis.py): Conducts random intervention analysis within RDiT.
  - [cross_sectional_association_analysis](1run_models/cross_sectional_association_analysis.R): Executes cross-sectional association analysis.
  - [longitudinal_association_analysis](1run_models/longitudinal_association_analysis.R): Runs longitudinal association analysis.
  - [run_liver_muscle_damage_association](1run_models/run_liver_muscle_damage_association.R): Runs associations of statin-influenced proteins with liver and muscle damage.
  - [two_sample_MR_NMR_Metabolomics.R](1run_models/two_sample_MR_NMR_Metabolomics.R): Conducts Mendelian Randomization (MR) analysis for NMR metabolites.
  - [two_sample_MR_Olink_Protein.R](1run_models/two_sample_MR_Olink_Protein.R): Conducts MR analysis for proteins measured on the Olink platform.
  - The Regression Discontinuity in Time (RDiT) model applies [rdroubst package](https://github.com/rdpackages/rdrobust/tree/master).
  - The power analysis for the RDiT model utilizes [rdpower packge](https://rdpackages.github.io/rdpower/).
  - The Mendelian Randomization model uses [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/).
