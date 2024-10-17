# Causal Inference of Drug Effects on Metabolomics and Proteomics in UK Biobank

## Description
This repository contains the downstream analysis for the paper  "﻿Causal Inference of Drug Effects on Metabolomics and Proteomics in UK Biobank". 

## Content
The repository is organized into the following directories:
- [0data_processing](0data_processing): Scripts to aggregate covariates and processing omics biomarkers in the [UK Biobank](https://www.ukbiobank.ac.uk/) for use in different blood biomarkers.
- ([1run_models](1run_models)): Scripts to reproduce the analyses and evaluate the causal associations between the ATC Level 4 Prescription and each omics biomarker.
  - The Regression Discontinuity in Time (RDiT) model applies [rdroubst package] (https://github.com/rdpackages/rdrobust/tree/master).
  - The Mendelian Randomization model uses [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/).
