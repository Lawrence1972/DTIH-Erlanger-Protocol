[[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20635281.svg)](https://doi.org/10.5281/zenodo.20635281)](https://doi.org/10.5281/zenodo.20635817)
# DTIH — Erlanger Protocol v0.3
### Drug-Induced Liver Injury Prediction via Mechanistically-Inspired Proxy Score (MIPS)

**Author:** Lawrence Choi (최민석) · Independent Researcher · South Korea
**Status:** Active development — not peer-reviewed

---

## What this is

A CYP-metabolism-based DILI (Drug-Induced Liver Injury) risk-scoring system,
validated retrospectively against the FDA DILIrank reference list.

Built by a non-institutional researcher using public pharmacokinetic data
(DrugBank, PK-DB) and AI-assisted development.

## What this is NOT

- Not a clinical decision tool
- Not validated for biologics, IV genotoxics, or mitochondrial toxicants
- Not peer-reviewed

## Key results

| Dataset | N drugs | AUC-ROC |
|---|---|---|
| Initial benchmark (v0.2) | 38 | 0.808 |
| Extended (v0.3) | 212 | 0.846 |
| In-domain primary (v0.3-ext) | 321 | 0.796 |
| Full set incl. out-of-domain (v0.3-ext) | 392 | 0.782 |

Tier-stratified discrimination (in-domain): vMost-DILI vs vNo-DILI **AUC = 0.871**.

## Files

| File | Description |
|---|---|
| `drug_pk_database.py` | 212-drug pharmacokinetic parameter library |
| `drug_pk_database_batch2.py` | +180 drugs via therapeutic-class imputation |
| `validation_v2.py` | Full verification & validation pipeline (AUC, ROC, Monte Carlo) |
| `simulator.py` | 5-layer ODE simulation (4 drugs, proof-of-concept) |
| `validation.py` | Original v0.2 validation module |
| `main.py` | Entry point for individual-patient simulations |
| `docs/DTIH_v03_Results_Draft.docx` | Draft Results manuscript section |
| `results/*.png` | Validation figures |

## How to run

```bash
pip install -r requirements.txt
python validation_v2.py
```

This reproduces the 212-drug validation (AUC-ROC 0.846) and writes a
six-panel figure to `results/`. All results use random seed 42.

## Applicability domain (honest limitations)

1. The MIPS score is a **weighted heuristic**, not a full ODE solution. The
   5-layer ODE model (`simulator.py`) runs only for 4 reference drugs as a
   proof-of-concept.
2. **71 of 392 drugs (18%)** fall outside the applicability domain:
   - Biologics (MW > 5 kDa) — immune-mediated injury
   - IV genotoxics (platinum, alkylating) — direct DNA damage
   - Mitochondrial toxicants (NRTIs) — non-CYP mechanism
3. Pharmacokinetic parameters for **28 drugs** use therapeutic-class average
   imputation (flagged `[IMP-CLASS]` in the source).

## Context of Use

This model is optimized for screening hepatotoxicity risk of **small-molecule
drugs undergoing CYP-mediated reactive-metabolite formation**, with potential
extension to patients carrying CYP genetic polymorphisms. Non-CYP mechanisms
are explicitly out of scope.

## License

MIT License — see `LICENSE`.

## Citation

If you reference this work, please cite the Zenodo archive (DOI to be added
after upload).
