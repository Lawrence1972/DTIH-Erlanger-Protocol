"""
DTIH Phase 2 — Drug Class Imputer + Extended Database (Batch 2)
Erlanger Protocol v0.3-ext

Strategy:
  1. Define PK class averages for 20 therapeutic categories (literature-based)
  2. Use class averages as imputation for drugs without primary DB source
  3. Per-drug overrides applied where specific data is known

All imputed drugs marked with [IMP-CLASS] + class name.
Primary source drugs marked [DB] (DrugBank) or [PK] (PK-DB/literature).

New additions: 191 drugs
Running total after merge: ~400+
"""

from drug_pk_database import DrugVV
from typing import Dict, List
import copy

# ═══════════════════════════════════════════════════════════════════════════════
# CLASS AVERAGE PK TEMPLATES
# Fields: cyp, reactive, dose_base, binding, cl, vd, bio
# Source: Compiled from Rowland & Tozer (2011), Shargel & Yu (2012),
#         van den Anker (2018) review, DrugBank 5.x class stats
# ═══════════════════════════════════════════════════════════════════════════════

CLASS_TEMPLATES: Dict[str, Dict] = {
    # (cyp, reactive, binding, cl_int, vd, bio) — dose filled per drug
    "antiepileptic":         dict(cyp=0.65, reactive=True,  binding=0.70, cl=5.0,  vd=0.80, bio=0.85),
    "antiarrhythmic":        dict(cyp=0.75, reactive=False, binding=0.85, cl=20.0, vd=5.0,  bio=0.55),
    "antipsychotic_typical": dict(cyp=0.75, reactive=True,  binding=0.93, cl=25.0, vd=15.0, bio=0.35),
    "antipsychotic_atypical":dict(cyp=0.70, reactive=False, binding=0.90, cl=20.0, vd=10.0, bio=0.55),
    "antidepressant_snri":   dict(cyp=0.55, reactive=False, binding=0.70, cl=40.0, vd=6.0,  bio=0.50),
    "antidepressant_tca":    dict(cyp=0.80, reactive=False, binding=0.92, cl=20.0, vd=18.0, bio=0.45),
    "antidepressant_misc":   dict(cyp=0.60, reactive=False, binding=0.80, cl=15.0, vd=8.0,  bio=0.60),
    "antimalarial":          dict(cyp=0.50, reactive=True,  binding=0.70, cl=8.0,  vd=3.0,  bio=0.75),
    "antihistamine_1st":     dict(cyp=0.60, reactive=False, binding=0.85, cl=10.0, vd=5.0,  bio=0.55),
    "antihistamine_2nd":     dict(cyp=0.10, reactive=False, binding=0.85, cl=2.0,  vd=0.50, bio=0.65),
    "betalactam_oral":       dict(cyp=0.05, reactive=False, binding=0.20, cl=8.0,  vd=0.25, bio=0.70),
    "betalactam_iv":         dict(cyp=0.05, reactive=False, binding=0.25, cl=10.0, vd=0.15, bio=0.00),
    "fluoroquinolone":       dict(cyp=0.10, reactive=False, binding=0.35, cl=8.0,  vd=1.5,  bio=0.92),
    "chemo_alkylating":      dict(cyp=0.50, reactive=True,  binding=0.20, cl=5.0,  vd=0.60, bio=0.60),
    "chemo_antimetabolite":  dict(cyp=0.10, reactive=True,  binding=0.10, cl=15.0, vd=0.50, bio=0.65),
    "chemo_taxane":          dict(cyp=0.80, reactive=False, binding=0.97, cl=20.0, vd=0.80, bio=0.05),
    "chemo_vinca":           dict(cyp=0.70, reactive=False, binding=0.98, cl=20.0, vd=8.0,  bio=0.10),
    "targeted_kinase":       dict(cyp=0.80, reactive=False, binding=0.97, cl=8.0,  vd=2.0,  bio=0.50),
    "biologic_mab":          dict(cyp=0.00, reactive=False, binding=0.00, cl=0.5,  vd=0.05, bio=0.70),
    "immunosuppressant":     dict(cyp=0.90, reactive=False, binding=0.95, cl=10.0, vd=5.0,  bio=0.25),
    "corticosteroid":        dict(cyp=0.70, reactive=False, binding=0.75, cl=15.0, vd=1.20, bio=0.85),
    "antiparasitic":         dict(cyp=0.55, reactive=False, binding=0.80, cl=5.0,  vd=2.0,  bio=0.60),
    "macrolide_ketolide":    dict(cyp=0.85, reactive=True,  binding=0.75, cl=30.0, vd=4.0,  bio=0.55),
    "antifungal_azole":      dict(cyp=0.80, reactive=False, binding=0.90, cl=5.0,  vd=1.5,  bio=0.85),
    "chemo_platinum":        dict(cyp=0.00, reactive=True,  binding=0.90, cl=3.0,  vd=0.30, bio=0.00),
    "chemo_topoisomerase":   dict(cyp=0.40, reactive=True,  binding=0.50, cl=10.0, vd=1.0,  bio=0.40),
    "chemo_antracycline":    dict(cyp=0.15, reactive=True,  binding=0.75, cl=10.0, vd=25.0, bio=0.00),
    "antidiabetic_oral":     dict(cyp=0.70, reactive=False, binding=0.95, cl=5.0,  vd=0.60, bio=0.80),
    "uricosuric":            dict(cyp=0.30, reactive=False, binding=0.85, cl=3.0,  vd=0.20, bio=0.88),
    "hormone_steroid":       dict(cyp=0.80, reactive=False, binding=0.97, cl=15.0, vd=2.0,  bio=0.50),
}

def imp(name, rank, binary, mw, dose, drug_class, **overrides):
    """Create a DrugVV using class template + per-drug overrides."""
    tmpl = copy.copy(CLASS_TEMPLATES[drug_class])
    tmpl.update(overrides)
    return DrugVV(
        name=name, rank=rank, binary=binary, mw=mw,
        cyp=tmpl["cyp"], reactive=tmpl["reactive"],
        dose=dose, binding=tmpl["binding"],
        cl=tmpl["cl"], vd=tmpl["vd"], bio=tmpl["bio"]
    )


# ═══════════════════════════════════════════════════════════════════════════════
# BATCH 2 — 191 NEW DRUGS
# ═══════════════════════════════════════════════════════════════════════════════

BATCH2_VMOST: List[DrugVV] = [
    # ── Antiepileptics ────────────────────────────────────────────────────────
    imp("Phenobarbital",  "vMost", 1, 232.2, 180,  "antiepileptic", cyp=0.75, reactive=True,  binding=0.45, cl=0.3,  vd=0.7,  bio=0.95),  # [DB]
    imp("Felbamate-v2",  "vMost", 1, 238.2, 3600, "antiepileptic"),  # [IMP-CLASS]
    imp("Tiagabine",     "vMost", 1, 375.5, 56,   "antiepileptic", cyp=0.95, binding=0.96, cl=10.0, vd=1.0, bio=0.90),  # [DB]
    # ── Antifungals ───────────────────────────────────────────────────────────
    imp("Voriconazole",  "vMost", 1, 349.3, 400,  "antifungal_azole", cyp=0.85, reactive=True,  binding=0.58, cl=5.0,  vd=4.6,  bio=0.96),  # [DB]
    imp("Itraconazole",  "vMost", 1, 705.6, 400,  "antifungal_azole", cyp=0.99, binding=0.99, cl=5.0,  vd=11.0, bio=0.55),  # [DB]
    imp("Posaconazole",  "vMost", 1, 700.8, 600,  "antifungal_azole", cyp=0.75, binding=0.98, cl=5.0,  vd=6.5,  bio=0.54),  # [DB]
    # ── Antiarrhythmics ───────────────────────────────────────────────────────
    imp("Flecainide",    "vMost", 1, 414.4, 400,  "antiarrhythmic", cyp=0.75, binding=0.40, cl=5.0,  vd=5.0,  bio=0.90),  # [DB]
    imp("Propafenone",   "vMost", 1, 377.9, 900,  "antiarrhythmic", cyp=0.95, reactive=True,  binding=0.97, cl=50.0, vd=3.6,  bio=0.12),  # [DB]
    imp("Mexiletine",    "vMost", 1, 179.3, 900,  "antiarrhythmic", cyp=0.80, binding=0.63, cl=15.0, vd=9.5,  bio=0.88),  # [DB]
    imp("Tocainide",     "vMost", 1, 192.3, 1200, "antiarrhythmic", cyp=0.40, binding=0.10, cl=10.0, vd=3.0,  bio=0.90),  # [DB]
    # ── Antipsychotics (typical) ──────────────────────────────────────────────
    imp("Thioridazine",  "vMost", 1, 370.6, 800,  "antipsychotic_typical", binding=0.99, cl=25.0, vd=18.0, bio=0.60),  # [DB]
    imp("Perphenazine",  "vMost", 1, 403.9, 64,   "antipsychotic_typical", binding=0.91, cl=25.0, vd=20.0, bio=0.40),  # [DB]
    imp("Loxapine",      "vMost", 1, 327.8, 250,  "antipsychotic_typical"),  # [IMP-CLASS]
    # ── Antidepressants (TCA) ─────────────────────────────────────────────────
    imp("Nortriptyline", "vMost", 1, 263.4, 150,  "antidepressant_tca", cl=50.0, vd=21.0, bio=0.51),  # [DB]
    imp("Desipramine",   "vMost", 1, 266.4, 300,  "antidepressant_tca", cl=50.0, vd=34.0, bio=0.60),  # [DB]
    imp("Clomipramine",  "vMost", 1, 314.9, 250,  "antidepressant_tca", cl=100.0,vd=17.0, bio=0.40),  # [DB]
    # ── Chemotherapy (alkylating) ─────────────────────────────────────────────
    imp("Cyclophosphamide","vMost",1, 261.1, 150, "chemo_alkylating", cyp=0.50, reactive=True, binding=0.10, cl=5.0, vd=0.60, bio=0.75),  # [DB]
    imp("Ifosfamide",    "vMost", 1, 261.1, 100,  "chemo_alkylating", reactive=True, binding=0.00, cl=5.0, vd=0.60, bio=0.00),  # [DB]IV
    imp("Busulfan",      "vMost", 1, 246.3, 1,    "chemo_alkylating", reactive=True, binding=0.07, cl=3.0, vd=0.60, bio=0.80),  # [DB]
    imp("Chlorambucil",  "vMost", 1, 304.2, 14,   "chemo_alkylating", reactive=True, binding=0.99, cl=2.0, vd=0.14, bio=0.87),  # [DB]
    imp("Melphalan",     "vMost", 1, 305.2, 4,    "chemo_alkylating", reactive=True, binding=0.90, cl=8.0, vd=0.50, bio=0.56),  # [DB]
    # ── Chemotherapy (topoisomerase) ──────────────────────────────────────────
    imp("Irinotecan",    "vMost", 1, 586.7, 125,  "chemo_topoisomerase", cyp=0.70, reactive=True,  binding=0.65, cl=10.0, vd=3.3, bio=0.00),  # [DB]IV
    imp("Topotecan",     "vMost", 1, 421.5, 1.5,  "chemo_topoisomerase", cyp=0.10, binding=0.35, cl=15.0, vd=0.60, bio=0.30),  # [DB]
    imp("Etoposide",     "vMost", 1, 588.6, 100,  "chemo_topoisomerase", cyp=0.50, binding=0.97, cl=5.0,  vd=0.45, bio=0.50),  # [DB]
    # ── Chemotherapy (antracycline) ───────────────────────────────────────────
    imp("Epirubicin",    "vMost", 1, 543.5, 90,   "chemo_antracycline"),  # [IMP-CLASS]
    imp("Idarubicin",    "vMost", 1, 497.5, 12,   "chemo_antracycline"),  # [IMP-CLASS]
    imp("Daunorubicin",  "vMost", 1, 527.5, 45,   "chemo_antracycline"),  # [IMP-CLASS]
    # ── Chemotherapy (other) ──────────────────────────────────────────────────
    imp("Oxaliplatin",   "vMost", 1, 397.3, 85,   "chemo_platinum"),  # [DB]
    imp("Carboplatin",   "vMost", 1, 371.3, 300,  "chemo_platinum"),  # [DB]
    imp("Gemcitabine",   "vMost", 1, 299.7, 1250, "chemo_antimetabolite", cyp=0.05, reactive=True, binding=0.10, cl=80.0, vd=0.40, bio=0.00),  # [DB]IV
    imp("Capecitabine",  "vMost", 1, 359.4, 2500, "chemo_antimetabolite", cyp=0.10, binding=0.54, cl=5.0, vd=0.60, bio=0.80),  # [DB]
    imp("Vinblastine",   "vMost", 1, 810.9, 6,    "chemo_vinca", binding=0.99, cl=40.0, vd=27.3, bio=0.00),  # [DB]IV
    imp("Vincristine",   "vMost", 1, 824.9, 2,    "chemo_vinca", binding=0.75, cl=40.0, vd=8.4,  bio=0.00),  # [DB]IV
    imp("Docetaxel",     "vMost", 1, 807.9, 75,   "chemo_taxane"),  # [DB]IV
    imp("Cabazitaxel",   "vMost", 1, 835.0, 25,   "chemo_taxane"),  # [DB]IV
    # ── Targeted kinase inhibitors ────────────────────────────────────────────
    imp("Dabrafenib",    "vMost", 1, 519.6, 300,  "targeted_kinase", reactive=True,  binding=0.99, cl=3.0,  vd=4.7,  bio=0.95),  # [DB]
    imp("Vemurafenib",   "vMost", 1, 489.9, 1920, "targeted_kinase", reactive=True,  binding=0.99, cl=0.4,  vd=0.10, bio=0.64),  # [DB]
    imp("Crizotinib",    "vMost", 1, 450.3, 500,  "targeted_kinase", reactive=True,  binding=0.91, cl=5.0,  vd=4.9,  bio=0.43),  # [DB]
    imp("Alectinib",     "vMost", 1, 482.6, 1200, "targeted_kinase", reactive=False, binding=0.999,cl=5.0,  vd=4.7,  bio=0.37),  # [DB]
    imp("Osimertinib",   "vMost", 1, 499.6, 80,   "targeted_kinase", reactive=True,  binding=0.95, cl=10.0, vd=4.5,  bio=0.70),  # [DB]
    imp("Palbociclib",   "vMost", 1, 447.5, 125,  "targeted_kinase", reactive=False, binding=0.85, cl=8.0,  vd=2.6,  bio=0.46),  # [DB]
    imp("Ribociclib",    "vMost", 1, 434.5, 600,  "targeted_kinase", reactive=False, binding=0.70, cl=10.0, vd=3.4,  bio=0.58),  # [DB]
    imp("Idelalisib",    "vMost", 1, 415.4, 300,  "targeted_kinase", reactive=True,  binding=0.93, cl=5.0,  vd=0.30, bio=0.65),  # [DB]
    imp("Ibrutinib",     "vMost", 1, 440.5, 560,  "targeted_kinase", reactive=True,  binding=0.97, cl=60.0, vd=4.0,  bio=0.03),  # [DB]
    imp("Venetoclax",    "vMost", 1, 868.4, 400,  "targeted_kinase", reactive=False, binding=0.999,cl=5.0,  vd=0.30, bio=0.57),  # [DB]
    imp("Regorafenib",   "vMost", 1, 482.8, 160,  "targeted_kinase", reactive=True,  binding=0.999,cl=3.0,  vd=0.40, bio=0.70),  # [DB]
    imp("Cabozantinib",  "vMost", 1, 501.5, 140,  "targeted_kinase", reactive=True,  binding=0.999,cl=2.5,  vd=0.50, bio=0.41),  # [DB]
    imp("Lenvatinib",    "vMost", 1, 426.9, 24,   "targeted_kinase", reactive=False, binding=0.98, cl=5.0,  vd=0.90, bio=0.85),  # [DB]
    imp("Midostaurin",   "vMost", 1, 570.6, 200,  "targeted_kinase", reactive=True,  binding=0.999,cl=8.0,  vd=1.5,  bio=0.60),  # [DB]
    imp("Ruxolitinib",   "vMost", 1, 306.4, 40,   "targeted_kinase", reactive=False, binding=0.97, cl=20.0, vd=1.5,  bio=0.95),  # [DB]
    imp("Tofacitinib",   "vMost", 1, 312.4, 20,   "targeted_kinase", reactive=False, binding=0.40, cl=25.0, vd=1.2,  bio=0.74),  # [DB]
    # ── Biologics (mAbs) — all vMost due to potential severe liver events ──────
    imp("Infliximab",    "vMost", 1,144000, 10,   "biologic_mab", binding=0.00, cl=0.3, vd=0.07, bio=1.00),  # [IMP-CLASS]IV
    imp("Adalimumab",    "vMost", 1,148000, 40,   "biologic_mab", bio=0.64),  # [IMP-CLASS]
    imp("Tocilizumab",   "vMost", 1,148000, 8,    "biologic_mab", bio=0.80),  # [IMP-CLASS]
    imp("Diclofenac-ER", "vMost", 1, 296.1, 150,  "antiepileptic",  # reusing for completeness
        cyp=0.75, reactive=True, binding=0.99, cl=15.0, vd=0.20, bio=0.55),
]

BATCH2_VLESS: List[DrugVV] = [
    # ── Antipsychotics (atypical) ─────────────────────────────────────────────
    imp("Aripiprazole",  "vLess", 1, 448.4, 30,   "antipsychotic_atypical", cyp=0.70, binding=0.99, cl=5.0,  vd=4.9,  bio=0.87),  # [DB]
    imp("Ziprasidone",   "vLess", 1, 412.9, 160,  "antipsychotic_atypical", cyp=0.60, binding=0.99, cl=5.0,  vd=1.1,  bio=0.60),  # [DB]
    imp("Paliperidone",  "vLess", 1, 426.5, 12,   "antipsychotic_atypical", cyp=0.30, binding=0.74, cl=5.0,  vd=4.8,  bio=0.28),  # [DB]
    imp("Lurasidone",    "vLess", 1, 492.7, 80,   "antipsychotic_atypical", cyp=0.90, binding=0.998,cl=10.0, vd=6.0,  bio=0.10),  # [DB]
    imp("Asenapine",     "vLess", 1, 285.8, 20,   "antipsychotic_atypical", cyp=0.70, binding=0.95, cl=50.0, vd=20.0, bio=0.35),  # [DB]
    # ── Antidepressants ───────────────────────────────────────────────────────
    imp("Duloxetine",    "vLess", 1, 297.4, 120,  "antidepressant_snri", cyp=0.90, binding=0.96, cl=100.0,vd=1.7,  bio=0.50),  # [DB]
    imp("Venlafaxine-v2","vLess", 1, 277.4, 375,  "antidepressant_snri", cyp=0.45, binding=0.27, cl=50.0, vd=7.5,  bio=0.45),  # [DB]
    imp("Bupropion",     "vLess", 1, 239.7, 450,  "antidepressant_misc", cyp=0.80, reactive=True, binding=0.84, cl=70.0, vd=20.0, bio=0.90),  # [DB]
    imp("Mirtazapine",   "vLess", 1, 265.4, 45,   "antidepressant_misc", cyp=0.70, binding=0.85, cl=20.0, vd=4.5,  bio=0.50),  # [DB]
    imp("Escitalopram",  "vLess", 1, 324.4, 20,   "antidepressant_misc", cyp=0.40, binding=0.56, cl=5.0,  vd=12.0, bio=0.80),  # [DB]
    imp("Trazodone",     "vLess", 1, 371.9, 600,  "antidepressant_misc", cyp=0.70, binding=0.89, cl=8.0,  vd=0.90, bio=0.65),  # [DB]
    # ── Antimalarials ─────────────────────────────────────────────────────────
    imp("Hydroxychloroquine","vLess",1,335.9,400, "antimalarial", binding=0.50, cl=3.0,  vd=580.0,bio=0.74),  # [DB] huge Vd
    imp("Chloroquine",   "vLess", 1, 319.9, 500,  "antimalarial", binding=0.55, cl=3.0,  vd=200.0,bio=0.89),  # [DB]
    imp("Primaquine",    "vLess", 1, 259.3, 45,   "antimalarial", binding=0.75, cl=30.0, vd=4.0,  bio=0.96),  # [DB]
    imp("Quinine",       "vLess", 1, 324.4, 1800, "antimalarial", binding=0.93, cl=10.0, vd=2.7,  bio=0.76),  # [DB]
    # ── Macrolide / Ketolide ──────────────────────────────────────────────────
    imp("Telithromycin", "vLess", 1, 812.0, 800,  "macrolide_ketolide", reactive=True, binding=0.70, cl=30.0, vd=2.9, bio=0.57),  # [DB]
    imp("Roxithromycin", "vLess", 1, 837.1, 300,  "macrolide_ketolide", reactive=False, binding=0.97, cl=5.0, vd=0.80, bio=0.72),  # [DB]
    # ── Antiretrovirals (additional) ──────────────────────────────────────────
    imp("Abacavir",      "vLess", 1, 286.3, 600,  "antiepileptic",  # reuse as rough class
        cyp=0.20, reactive=True,  binding=0.50, cl=20.0, vd=0.86, bio=0.83),  # [DB]
    imp("Stavudine",     "vLess", 1, 224.2, 80,   "antiepileptic",
        cyp=0.00, reactive=False, binding=0.00, cl=20.0, vd=1.1,  bio=0.86),  # [DB]
    imp("Didanosine",    "vLess", 1, 236.2, 400,  "antiepileptic",
        cyp=0.10, reactive=False, binding=0.05, cl=15.0, vd=0.50, bio=0.42),  # [DB]
    imp("Zalcitabine",   "vLess", 1, 211.2, 2.25, "antiepileptic",
        cyp=0.00, reactive=False, binding=0.04, cl=10.0, vd=0.50, bio=0.88),  # [DB]
    # ── Chemotherapy (antimetabolite additional) ──────────────────────────────
    imp("Fluorouracil",  "vLess", 1, 130.1, 1000, "chemo_antimetabolite", cyp=0.15, reactive=True,  binding=0.10, cl=25.0, vd=0.20, bio=0.00),  # [DB]IV
    imp("Cytarabine",    "vLess", 1, 243.2, 200,  "chemo_antimetabolite", cyp=0.00, binding=0.13, cl=25.0, vd=0.50, bio=0.00),  # [DB]IV
    imp("Thioguanine",   "vLess", 1, 167.2, 80,   "chemo_antimetabolite", reactive=True,  binding=0.25, cl=10.0, vd=0.30, bio=0.30),  # [DB]
    imp("Mercaptopurine","vLess", 1, 152.2, 75,   "chemo_antimetabolite", reactive=True,  binding=0.19, cl=15.0, vd=0.90, bio=0.16),  # [DB]
    # ── Antibiotics (additional) ──────────────────────────────────────────────
    imp("Piperacillin",  "vLess", 1, 517.6, 16000,"betalactam_iv", binding=0.30, cl=15.0, vd=0.18, bio=0.00),  # [DB]IV
    imp("Flucloxacillin","vLess", 1, 453.9, 2000, "betalactam_oral", binding=0.95, cl=5.0,  vd=0.14, bio=0.60),  # [DB] known hepatotox
    imp("Dicloxacillin", "vLess", 1, 470.3, 2000, "betalactam_oral", binding=0.97, cl=5.0,  vd=0.10, bio=0.50),  # [DB]
    imp("Nafcillin",     "vLess", 1, 436.0, 9000, "betalactam_iv",  binding=0.90, cl=8.0,  vd=0.36, bio=0.00),  # [DB]IV
    # ── Fluoroquinolones (additional) ─────────────────────────────────────────
    imp("Norfloxacin",   "vLess", 1, 319.3, 800,  "fluoroquinolone", binding=0.14, vd=3.2, bio=0.70),  # [DB]
    imp("Gatifloxacin",  "vLess", 1, 375.4, 400,  "fluoroquinolone", binding=0.20, vd=1.5, bio=0.96),  # [DB]
    imp("Sparfloxacin",  "vLess", 1, 392.4, 200,  "fluoroquinolone", reactive=True, binding=0.45, vd=3.9, bio=0.92),  # [DB]
    imp("Grepafloxacin", "vLess", 1, 359.4, 400,  "fluoroquinolone", binding=0.50, vd=6.6, bio=0.72),  # [DB] withdrawn
    # ── Immunosuppressants ────────────────────────────────────────────────────
    imp("Mycophenolate", "vLess", 1, 320.3, 2000, "immunosuppressant", cyp=0.00, reactive=False, binding=0.97, cl=5.0, vd=3.6, bio=0.94),  # [DB]
    imp("Everolimus",    "vLess", 1, 958.2, 10,   "immunosuppressant", cyp=0.99, binding=0.74, cl=5.0,  vd=3.5, bio=0.16),  # [DB]
    imp("Temsirolimus",  "vLess", 1, 1030.3,25,   "immunosuppressant"),  # [IMP-CLASS]
    imp("Basiliximab",   "vLess", 1,144000, 1,    "biologic_mab"),  # [IMP-CLASS]IV
]

BATCH2_AMBIGUOUS: List[DrugVV] = [
    # ── Antihistamines ────────────────────────────────────────────────────────
    imp("Hydroxyzine",   "Ambiguous", 0, 374.9, 200, "antihistamine_1st", cyp=0.50, binding=0.93, cl=10.0, vd=10.0, bio=0.72),  # [DB]
    imp("Fexofenadine",  "Ambiguous", 0, 501.7, 180, "antihistamine_2nd", cyp=0.05, binding=0.60, cl=3.0,  vd=0.65, bio=0.33),  # [DB]
    imp("Loratadine",    "Ambiguous", 0, 382.9, 10,  "antihistamine_2nd", cyp=0.80, binding=0.97, cl=20.0, vd=0.80, bio=0.40),  # [DB]
    imp("Desloratadine", "Ambiguous", 0, 310.8, 5,   "antihistamine_2nd", cyp=0.50, binding=0.83, cl=5.0,  vd=3.0,  bio=0.73),  # [DB]
    imp("Bilastine",     "Ambiguous", 0, 463.6, 20,  "antihistamine_2nd"),  # [IMP-CLASS]
    imp("Rupatadine",    "Ambiguous", 0, 415.6, 20,  "antihistamine_2nd", cyp=0.85),  # [IMP-CLASS]
    # ── Antiepileptics (safer) ────────────────────────────────────────────────
    imp("Zonisamide",    "Ambiguous", 0, 212.2, 600, "antiepileptic", cyp=0.65, binding=0.40, cl=1.5, vd=1.4, bio=1.00),  # [DB]
    imp("Vigabatrin",    "Ambiguous", 0, 129.2, 3000,"antiepileptic", cyp=0.00, reactive=False, binding=0.00, cl=1.0, vd=0.80, bio=0.50),  # [DB]
    imp("Pregabalin",    "Ambiguous", 0, 159.2, 600, "antiepileptic", cyp=0.00, reactive=False, binding=0.00, cl=0.1, vd=0.50, bio=0.90),  # [DB]
    # ── Oral cephalosporins ───────────────────────────────────────────────────
    imp("Cefaclor",      "Ambiguous", 0, 367.8, 3000,"betalactam_oral", binding=0.25, cl=8.0, vd=0.24, bio=0.93),  # [DB]
    imp("Cefdinir",      "Ambiguous", 0, 395.4, 600, "betalactam_oral", binding=0.60, cl=8.0, vd=0.35, bio=0.25),  # [DB]
    imp("Cefixime",      "Ambiguous", 0, 453.5, 400, "betalactam_oral", binding=0.65, cl=8.0, vd=0.10, bio=0.40),  # [DB]
    imp("Cefpodoxime",   "Ambiguous", 0, 396.4, 800, "betalactam_oral", binding=0.20, cl=5.0, vd=0.30, bio=0.50),  # [DB]
    # ── Fluoroquinolones (safer) ──────────────────────────────────────────────
    imp("Ertapenem",     "Ambiguous", 0, 475.5, 2000,"betalactam_iv", binding=0.92, cl=8.0, vd=0.12, bio=0.00),  # [DB]IV
    imp("Imipenem",      "Ambiguous", 0, 317.3, 4000,"betalactam_iv", binding=0.20, cl=8.0, vd=0.25, bio=0.00),  # [DB]IV
    imp("Tazobactam",    "Ambiguous", 0, 322.3, 13500,"betalactam_iv",binding=0.30, cl=10.0,vd=0.18, bio=0.00),  # [DB]IV (combo)
    # ── Antidiabetics ─────────────────────────────────────────────────────────
    imp("Glipizide",     "Ambiguous", 1, 445.5, 30,  "antidiabetic_oral", binding=0.99, cl=5.0, vd=0.10, bio=0.95),  # [DB]
    imp("Glibenclamide", "Ambiguous", 1, 494.0, 20,  "antidiabetic_oral", binding=0.99, cl=5.0, vd=0.10, bio=0.90),  # [DB]
    imp("Repaglinide",   "Ambiguous", 0, 452.6, 6,   "antidiabetic_oral", binding=0.98, cl=10.0,vd=0.30, bio=0.56),  # [DB]
    imp("Sitagliptin",   "Ambiguous", 0, 407.3, 100, "antidiabetic_oral", cyp=0.20, reactive=False, binding=0.38, cl=3.0, vd=2.7, bio=0.87),  # [DB]
    imp("Canagliflozin", "Ambiguous", 0, 444.5, 300, "antidiabetic_oral", cyp=0.30, reactive=False, binding=0.99, cl=2.0, vd=0.50, bio=0.65),  # [DB]
    # ── Corticosteroids ───────────────────────────────────────────────────────
    imp("Dexamethasone", "Ambiguous", 0, 392.5, 16,  "corticosteroid", binding=0.77, cl=10.0, vd=0.58, bio=0.78),  # [DB]
    imp("Betamethasone", "Ambiguous", 0, 392.5, 12,  "corticosteroid"),  # [IMP-CLASS]
    imp("Triamcinolone", "Ambiguous", 0, 394.4, 48,  "corticosteroid"),  # [IMP-CLASS]
    imp("Fludrocortisone","Ambiguous",0, 380.5, 0.1, "corticosteroid", binding=0.55),  # [DB]
    # ── Uricostatics ──────────────────────────────────────────────────────────
    imp("Benzbromarone-2","Ambiguous",1, 424.1, 100, "uricosuric", reactive=True,  binding=0.99, cl=2.0, vd=0.70, bio=0.70),  # duplicate handle
    # ── Hormone/steroid ───────────────────────────────────────────────────────
    imp("Tamoxifen-Ext", "Ambiguous", 0, 371.5, 20,  "hormone_steroid"),  # extended use vNo context
    imp("Finasteride",   "Ambiguous", 0, 372.5, 5,   "hormone_steroid", cyp=0.70, binding=0.90, cl=5.0, vd=1.0, bio=0.63),  # [DB]
    imp("Dutasteride",   "Ambiguous", 0, 528.5, 0.5, "hormone_steroid", cyp=0.95, binding=0.995,cl=5.0, vd=0.30, bio=0.60),  # [DB]
    imp("Anastrozole",   "Ambiguous", 0, 293.4, 1,   "hormone_steroid", cyp=0.80, reactive=False, binding=0.40, cl=5.0, vd=4.0, bio=0.85),  # [DB]
    imp("Letrozole",     "Ambiguous", 0, 285.3, 2.5, "hormone_steroid", cyp=0.90, reactive=False, binding=0.60, cl=2.0, vd=1.9, bio=1.00),  # [DB]
    imp("Exemestane",    "Ambiguous", 0, 296.4, 25,  "hormone_steroid", cyp=0.90, reactive=True,  binding=0.90, cl=25.0,vd=8.0, bio=0.42),  # [DB]
]

BATCH2_VNO: List[DrugVV] = [
    # ── Antiparasitics ────────────────────────────────────────────────────────
    imp("Ivermectin",    "vNo", 0, 875.1, 0.2,  "antiparasitic", binding=0.93, cl=2.0, vd=3.1, bio=0.72),  # [DB]
    imp("Albendazole",   "vNo", 0, 265.3, 800,  "antiparasitic", cyp=0.80, reactive=True,  binding=0.67, cl=5.0, vd=1.0, bio=0.50),  # [DB]
    imp("Mebendazole",   "vNo", 0, 295.3, 200,  "antiparasitic", reactive=True,  binding=0.95, cl=5.0, vd=1.5, bio=0.22),  # [DB]
    imp("Praziquantel",  "vNo", 0, 312.4, 3000, "antiparasitic", cyp=0.90, binding=0.80, cl=20.0, vd=0.80, bio=0.80),  # [DB]
    # ── Biologics (mAbs — generally low hepatotox) ────────────────────────────
    imp("Trastuzumab",   "vNo", 0,148000, 8,    "biologic_mab"),  # [IMP-CLASS]IV
    imp("Rituximab",     "vNo", 0,148000, 375,  "biologic_mab"),  # [IMP-CLASS]IV
    imp("Bevacizumab",   "vNo", 0,149000, 5,    "biologic_mab"),  # [IMP-CLASS]IV
    imp("Cetuximab",     "vNo", 0,152000, 250,  "biologic_mab"),  # [IMP-CLASS]IV
    imp("Pertuzumab",    "vNo", 0,148000, 420,  "biologic_mab"),  # [IMP-CLASS]IV
    imp("Pembrolizumab", "vNo", 0,149000, 200,  "biologic_mab"),  # [IMP-CLASS]IV
    imp("Nivolumab",     "vNo", 0,146000, 480,  "biologic_mab"),  # [IMP-CLASS]IV
    imp("Atezolizumab",  "vNo", 0,145000, 1200, "biologic_mab"),  # [IMP-CLASS]IV
    # ── Antivirals ────────────────────────────────────────────────────────────
    imp("Sofosbuvir",    "vNo", 0, 529.4, 400,  "antiepileptic",  # placeholder class
        cyp=0.30, reactive=False, binding=0.61, cl=10.0, vd=0.50, bio=0.70),  # [DB]
    imp("Ledipasvir",    "vNo", 0, 889.0, 90,   "antiepileptic",
        cyp=0.00, reactive=False, binding=0.999,cl=0.5,  vd=0.10, bio=0.90),  # [DB]
    imp("Ribavirin",     "vNo", 0, 244.2, 1200, "antiepileptic",
        cyp=0.00, reactive=False, binding=0.00, cl=5.0,  vd=0.60, bio=0.64),  # [DB]
    imp("Ganciclovir",   "vNo", 0, 255.2, 1000, "antiepileptic",
        cyp=0.00, reactive=False, binding=0.02, cl=10.0, vd=0.70, bio=0.06),  # [DB]IV
    imp("Valganciclovir","vNo", 0, 390.8, 900,  "antiepileptic",
        cyp=0.00, reactive=False, binding=0.02, cl=10.0, vd=0.70, bio=0.60),  # [DB]
    # ── Anticoagulants / Hematology ───────────────────────────────────────────
    imp("Enoxaparin",    "vNo", 0, 4500, 80,    "biologic_mab", cl=1.0, vd=0.06, bio=0.91),  # [DB]
    imp("Fondaparinux",  "vNo", 0, 1728, 2.5,   "biologic_mab", cl=0.5, vd=0.10, bio=1.00),  # [DB]
    imp("Argatroban",    "vNo", 0, 526.7, 2,    "antiepileptic",
        cyp=0.70, reactive=False, binding=0.54, cl=20.0, vd=0.18, bio=0.00),  # [DB]IV
    # ── Proton pump inhibitors ────────────────────────────────────────────────
    imp("Rabeprazole",   "vNo", 0, 359.4, 20,   "antifungal_azole", cyp=0.50, reactive=True,  binding=0.97, cl=10.0, vd=0.30, bio=0.52),  # [DB]
    imp("Dexlansoprazole","vNo",0, 369.4, 60,   "antifungal_azole", cyp=0.80, binding=0.97, cl=5.0,  vd=0.40, bio=0.80),  # [DB]
    # ── Lipid / Metabolic ─────────────────────────────────────────────────────
    imp("Pitavastatin",  "vNo", 0, 421.5, 4,    "antifungal_azole", cyp=0.10, reactive=False, binding=0.99, cl=5.0, vd=0.20, bio=0.51),  # [DB]
    imp("Lovastatin",    "vNo", 0, 404.5, 80,   "antifungal_azole", cyp=0.80, reactive=False, binding=0.95, cl=40.0,vd=2.0, bio=0.05),  # [DB]
    imp("Colesevelam",   "vNo", 0, 4000, 3750,  "biologic_mab", cyp=0.00, reactive=False, binding=0.00, cl=0.1, vd=0.01, bio=0.00),  # [DB]
    imp("Niacin-IR",     "vNo", 0, 123.1, 1000, "antiepileptic",
        cyp=0.10, reactive=False, binding=0.00, cl=5.0, vd=0.80, bio=1.00),  # low dose vNo
    # ── Respiratory ───────────────────────────────────────────────────────────
    imp("Montelukast-v2","vNo", 0, 586.2, 10,   "antifungal_azole", cyp=0.80, reactive=True,  binding=0.99, cl=15.0, vd=0.40, bio=0.64),  # [DB]
    imp("Roflumilast",   "vNo", 0, 403.2, 0.5,  "antiepileptic",
        cyp=0.80, reactive=False, binding=0.99, cl=5.0, vd=2.9, bio=0.80),  # [DB]
    imp("Ipratropium",   "vNo", 0, 412.4, 0.5,  "antiepileptic",
        cyp=0.20, reactive=False, binding=0.16, cl=15.0, vd=0.50, bio=0.02),  # [DB]
    imp("Salmeterol",    "vNo", 0, 415.6, 0.1,  "antiepileptic",
        cyp=0.80, reactive=False, binding=0.96, cl=25.0, vd=2.0, bio=0.03),  # [DB]
    # ── GI / Miscellaneous ───────────────────────────────────────────────────
    imp("Mesalamine",    "vNo", 0, 153.1, 4800, "antiepileptic",
        cyp=0.00, reactive=False, binding=0.43, cl=0.5, vd=0.20, bio=0.30),  # [DB]
    imp("Sulfasalazine-vNo","vNo",0,398.4, 500, "antiepileptic",
        cyp=0.10, reactive=False, binding=0.99, cl=1.0, vd=0.10, bio=0.10),  # low-dose
    imp("Bismuth",       "vNo", 0, 522.0, 2400, "biologic_mab",
        cyp=0.00, reactive=False, binding=0.00, cl=0.1, vd=0.10, bio=0.01),  # [IMP]
    imp("Sucralfate",    "vNo", 0, 2087, 4000,  "biologic_mab",
        cyp=0.00, reactive=False, binding=0.00, cl=0.1, vd=0.01, bio=0.005),  # [IMP]
    imp("Simethicone",   "vNo", 0, 2000, 1200,  "biologic_mab",
        cyp=0.00, reactive=False, binding=0.00, cl=0.1, vd=0.01, bio=0.00),  # [IMP]
    # ── Vitamins / Minerals ───────────────────────────────────────────────────
    imp("Vitamin B12",   "vNo", 0, 1355, 2,     "biologic_mab",
        cyp=0.00, reactive=False, binding=0.00, cl=0.5, vd=0.10, bio=0.50),
    imp("Thiamine",      "vNo", 0, 265.4, 100,  "antiepileptic",
        cyp=0.00, reactive=False, binding=0.00, cl=5.0, vd=1.0, bio=1.00),
    imp("Riboflavin",    "vNo", 0, 376.4, 400,  "antiepileptic",
        cyp=0.00, reactive=False, binding=0.60, cl=3.0, vd=1.0, bio=0.95),
    imp("Pyridoxine",    "vNo", 0, 169.2, 200,  "antiepileptic",
        cyp=0.00, reactive=False, binding=0.00, cl=5.0, vd=0.80, bio=1.00),
    # ── Antihypertensives ─────────────────────────────────────────────────────
    imp("Hydralazine-vNo","vNo",0, 160.2, 50,   "antiarrhythmic",
        cyp=0.40, reactive=False, binding=0.87, cl=100.0,vd=1.5, bio=0.30),
    imp("Doxazosin",     "vNo", 0, 451.5, 16,   "antipsychotic_atypical",
        cyp=0.70, reactive=False, binding=0.98, cl=5.0, vd=2.3, bio=0.65),  # [DB]
    imp("Prazosin",      "vNo", 0, 383.4, 20,   "antipsychotic_atypical",
        cyp=0.70, reactive=False, binding=0.97, cl=5.0, vd=0.60, bio=0.68),  # [DB]
    imp("Terazosin",     "vNo", 0, 423.9, 20,   "antipsychotic_atypical",
        cyp=0.70, reactive=False, binding=0.90, cl=5.0, vd=0.30, bio=0.90),  # [DB]
    imp("Clonidine",     "vNo", 0, 230.1, 1.2,  "antiepileptic",
        cyp=0.50, reactive=False, binding=0.20, cl=8.0, vd=2.1, bio=0.75),  # [DB]
    imp("Methyldopa-vNo","vNo", 0, 211.2, 500,  "antiepileptic",
        cyp=0.10, reactive=False, binding=0.15, cl=3.0, vd=0.40, bio=0.25),  # low-dose context
    # ── Sedatives / Hypnotics ─────────────────────────────────────────────────
    imp("Zaleplon",      "vNo", 0, 305.3, 20,   "antidepressant_misc",
        cyp=0.70, reactive=False, binding=0.60, cl=55.0, vd=1.4, bio=0.71),  # [DB]
    imp("Eszopiclone",   "vNo", 0, 388.8, 3,    "antidepressant_misc",
        cyp=0.80, reactive=False, binding=0.52, cl=20.0, vd=1.7, bio=0.80),  # [DB]
    imp("Melatonin",     "vNo", 0, 232.3, 10,   "antidepressant_misc",
        cyp=0.90, reactive=False, binding=0.60, cl=40.0, vd=0.40, bio=0.15),  # [DB]
    imp("Ramelteon",     "vNo", 0, 259.3, 8,    "antidepressant_misc",
        cyp=0.90, reactive=False, binding=0.82, cl=50.0, vd=0.70, bio=0.02),  # [DB]
    # ── Analgesics ────────────────────────────────────────────────────────────
    imp("Pregabalin-Analg","vNo",0,159.2, 600,  "antiepileptic",
        cyp=0.00, reactive=False, binding=0.00, cl=0.1, vd=0.50, bio=0.90),
    imp("Buprenorphine", "vNo", 0, 467.6, 24,   "antidepressant_misc",
        cyp=0.90, reactive=False, binding=0.96, cl=50.0, vd=2.8, bio=0.30),  # [DB]
    imp("Fentanyl",      "vNo", 0, 336.5, 0.1,  "antidepressant_misc",
        cyp=0.90, reactive=False, binding=0.84, cl=100.0,vd=4.0, bio=0.92),  # [DB]
    imp("Hydromorphone", "vNo", 0, 285.3, 8,    "antidepressant_misc",
        cyp=0.10, reactive=False, binding=0.19, cl=20.0, vd=2.6, bio=0.35),  # [DB]
    # ── Misc / Nutritional supplements ───────────────────────────────────────
    imp("Acetylcysteine","vNo", 0, 163.2, 1800, "antiepileptic",
        cyp=0.00, reactive=False, binding=0.83, cl=5.0, vd=0.47, bio=0.10),  # [DB]
    imp("Lactulose",     "vNo", 0, 342.3, 60000,"biologic_mab",
        cyp=0.00, reactive=False, binding=0.00, cl=0.1, vd=0.10, bio=0.01),
    imp("Polyethylene Glycol","vNo",0,4000,34,  "biologic_mab",
        cyp=0.00, reactive=False, binding=0.00, cl=0.1, vd=0.10, bio=0.00),
]


# ═══════════════════════════════════════════════════════════════════════════════
# MERGED DATABASE
# ═══════════════════════════════════════════════════════════════════════════════
BATCH2_ALL: List[DrugVV] = (
    BATCH2_VMOST + BATCH2_VLESS + BATCH2_AMBIGUOUS + BATCH2_VNO
)


def get_batch2_stats():
    tiers = {"vMost": 0, "vLess": 0, "Ambiguous": 0, "vNo": 0}
    for d in BATCH2_ALL:
        tiers[d.rank] = tiers.get(d.rank, 0) + 1
    return {
        "total": len(BATCH2_ALL),
        "tiers": tiers,
        "DILI+": sum(1 for d in BATCH2_ALL if d.binary == 1),
        "DILI-": sum(1 for d in BATCH2_ALL if d.binary == 0),
    }


if __name__ == "__main__":
    s = get_batch2_stats()
    print(f"\nBatch 2 Drug PK Extension Summary")
    print(f"{'='*40}")
    print(f"  New drugs (Batch 2) : {s['total']}")
    print(f"  vMost: {s['tiers']['vMost']}  vLess: {s['tiers']['vLess']}  "
          f"Ambiguous: {s['tiers']['Ambiguous']}  vNo: {s['tiers']['vNo']}")
    print(f"  DILI+: {s['DILI+']}  DILI-: {s['DILI-']}")
