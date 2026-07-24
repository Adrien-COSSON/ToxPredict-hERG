---
title: toxpredict-mlflow
emoji: 🏢
colorFrom: pink
colorTo: purple
sdk: docker
pinned: false
---

# ToxPredict-hERG

> Machine learning pipeline for hERG cardiotoxicity prediction from molecular structures.

---

## Overview

hERG (KCNH2) channel inhibition is a major cause of drug-induced cardiotoxicity and a critical safety endpoint in drug development. Regulatory agencies require hERG assessment for all new drug candidates.

**ToxPredict-hERG** builds a predictive model for hERG inhibitory activity (pIC50) from molecular structures, enabling early-stage virtual screening before synthesis.

---

## Dataset

- **Source:** ChEMBL (target: CHEMBL240 — KCNH2/hERG)
- **Endpoint:** IC50 (inhibition of hERG channel)
- **Volume:** ~19,000 entries
- **Processing:** filtering, unit standardization, pIC50 conversion

---

## Pipeline

```
ChEMBL extraction
      ↓
Data cleaning & curation
      ↓
Feature engineering (RDKit descriptors + Morgan fingerprints)
      ↓
ML models (Dummy → Ridge → Random Forest → XGBoost)
      ↓
DL models (MLP → GNN)
      ↓
Evaluation & interpretation
```

---

## Features

- **Physicochemical descriptors:** MW, LogP, HBD, HBA, TPSA (RDKit / Lipinski)
- **Molecular fingerprints:** Morgan/ECFP (2048-bit binary vectors)

---

## Models

| Model | Type |
|---|---|
| Dummy Regressor | Baseline |
| Ridge Regression | Linear |
| Random Forest | Ensemble |
| XGBoost | Gradient Boosting |
| MLP | Deep Learning |
| GNN (DeepChem / PyTorch Geometric) | Deep Learning |

---

## Tech Stack

- `python` — version 3.11.9
- `chembl_webresource_client` — data extraction
- `RDKit` — molecular feature engineering
- `scikit-learn`, `XGBoost` — classical ML
- `DeepChem` / `PyTorch Geometric` — deep learning
- `pandas`, `numpy`, `matplotlib`, `seaborn` — data handling & visualization

---

## Project Status

🚧 Work in progress

- [x] Target selection & feasibility analysis
- [ ] Data extraction & cleaning
- [ ] Feature engineering
- [ ] ML pipeline
- [ ] DL pipeline
- [ ] Evaluation & reporting

---

## Author

Adrien — PhD in Oncohematology | Fullstack AI (Jedha)  
Bridging pharmaceutical biology expertise with machine learning for drug discovery.