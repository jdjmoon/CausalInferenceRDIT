# CausalInferenceRDIT

## Description
Here we present the downstream analysis for the paper "﻿Causal inference of drug effects on metabolomics and proteomics in UK Biobank". 

## Content
This repository contains code to aggregate and process omics biomarkers for each method  ([0data_processing](0data_processing)) and reproduce the analyses ([1run_models](1run_models)) in the [UK Biobank](https://www.ukbiobank.ac.uk/).

- Aggregating the covariates used for the models with the blood date to compute the time between the UK Biobank's blood sampling dates and the medication prescription from the UK Biobank's GP Prescription Data.
- Aggregating the above covariates to different blood omics biomarkers to prepare for various models.
- Performing winsorization and evaluating the causal association between ATC Level 4 Prescription and each omics biomarkers.
