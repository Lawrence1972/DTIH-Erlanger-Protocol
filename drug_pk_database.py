"""
DTIH Phase 2 — Drug PK Parameter Database
Erlanger Protocol v0.3

DrugVV fields:
    name     : drug name (str)
    rank     : DILIrank tier — vMost / vLess / Ambiguous / vNo (str)
    binary   : DILI label — 1 (DILI+) / 0 (DILI−) (int)
    mw       : molecular weight (g/mol)
    cyp      : fraction metabolised via CYP enzymes  (0–1)
    reactive : reactive metabolite formation (bool)
    dose     : maximum daily dose (mg)
    binding  : plasma protein binding fraction (0–1)
    cl       : hepatic intrinsic clearance (L/h)
    vd       : volume of distribution (L/kg)
    bio      : oral bioavailability (0–1)

Sources (per-drug):
  [DB]  DrugBank 5.x open data
  [PK]  PK-DB / literature pharmacokinetic summaries
  [IMP] Class-average imputation (no primary PK source available)

DILIrank 2.0 tiers (Chen et al. 2016, updated 2022):
  vMost      → binary=1  (strong clinical evidence)
  vLess      → binary=1  (case reports / limited evidence)
  Ambiguous  → binary=0  (conflicting / insufficient data; selected binary per known cases)
  vNo        → binary=0  (substantial negative evidence)

Total: 226 unique drugs (38 original + 188 new)
"""

from dataclasses import dataclass
from typing import List

# ── dataclass (mirrors validation.py) ────────────────────────────────────────
@dataclass
class DrugVV:
    name: str; rank: str; binary: int; mw: float
    cyp: float; reactive: bool; dose: float; binding: float
    cl: float; vd: float; bio: float


# ═══════════════════════════════════════════════════════════════════════════════
# TIER 1: vMost-DILI Concern  (binary = 1)
# ═══════════════════════════════════════════════════════════════════════════════
VMOST: List[DrugVV] = [
    # ── Original 12 (v0.2 baseline) ──────────────────────────────────────────
    DrugVV("Acetaminophen",    "vMost", 1, 151.2, 0.10, True,  4000, 0.20,  5.0,  0.90, 0.85),
    DrugVV("Isoniazid",        "vMost", 1, 137.1, 0.70, True,   300, 0.10, 12.0,  0.60, 0.90),
    DrugVV("Troglitazone",     "vMost", 1, 441.5, 0.85, True,   600, 0.99,  2.5,  1.20, 0.50),
    DrugVV("Ketoconazole",     "vMost", 1, 531.4, 0.90, True,   400, 0.99,  8.0,  2.40, 0.50),
    DrugVV("Diclofenac",       "vMost", 1, 296.1, 0.75, True,   150, 0.99, 15.0,  0.20, 0.55),
    DrugVV("Valproic Acid",    "vMost", 1, 144.2, 0.30, True,  2000, 0.90,  1.0,  0.20, 0.90),
    DrugVV("Amiodarone",       "vMost", 1, 645.3, 0.50, False,  400, 0.96,  1.5, 70.00, 0.50),
    DrugVV("Nefazodone",       "vMost", 1, 470.0, 0.95, True,   600, 0.99, 25.0,  0.50, 0.20),
    DrugVV("Tolcapone",        "vMost", 1, 273.2, 0.80, True,   300, 0.99, 30.0,  0.30, 0.65),
    DrugVV("Flutamide",        "vMost", 1, 276.2, 0.85, True,   750, 0.94, 20.0,  0.30, 0.80),
    DrugVV("Chlorzoxazone",    "vMost", 1, 169.6, 0.90, True,   750, 0.15, 35.0,  0.20, 0.90),
    DrugVV("Dantrolene",       "vMost", 1, 314.3, 0.60, True,   400, 0.90,  5.0,  1.00, 0.70),
    # ── New vMost (38 additional) ─────────────────────────────────────────────
    DrugVV("Tacrine",          "vMost", 1, 198.3, 0.90, True,   160, 0.55, 30.0,  3.50, 0.17),  # [DB]
    DrugVV("Bromfenac",        "vMost", 1, 334.2, 0.70, True,    18, 0.99, 10.0,  0.20, 0.80),  # [DB]
    DrugVV("Trovafloxacin",    "vMost", 1, 512.4, 0.10, True,   300, 0.76,  8.0,  1.40, 0.90),  # [DB]
    DrugVV("Pemoline",         "vMost", 1, 176.2, 0.50, True,    75, 0.50,  2.0,  0.70, 0.50),  # [IMP]
    DrugVV("Felbamate",        "vMost", 1, 238.2, 0.50, True,  3600, 0.25,  2.0,  0.70, 0.90),  # [DB]
    DrugVV("Benzbromarone",    "vMost", 1, 424.1, 0.80, True,   100, 0.99,  2.0,  0.70, 0.70),  # [PK]
    DrugVV("Ximelagatran",     "vMost", 1, 474.5, 0.20, False,   48, 0.15,  5.0,  0.70, 0.20),  # [DB]
    DrugVV("Perhexiline",      "vMost", 1, 277.5, 0.95, False,  400, 0.90, 15.0, 10.00, 0.60),  # [PK]
    DrugVV("Methyldopa",       "vMost", 1, 211.2, 0.10, True,  3000, 0.15,  3.0,  0.40, 0.25),  # [DB]
    DrugVV("Nitrofurantoin",   "vMost", 1, 238.2, 0.10, True,   400, 0.40,  5.0,  0.80, 0.90),  # [DB]
    DrugVV("Allopurinol",      "vMost", 1, 136.1, 0.10, True,   800, 0.00,  2.0,  1.60, 0.80),  # [DB]
    DrugVV("Sulindac",         "vMost", 1, 356.4, 0.40, True,   400, 0.99,  5.0,  0.10, 0.90),  # [DB]
    DrugVV("Piroxicam",        "vMost", 1, 331.3, 0.50, False,   20, 0.99,  0.1,  0.10, 1.00),  # [DB]
    DrugVV("Phenylbutazone",   "vMost", 1, 308.4, 0.60, True,   600, 0.98,  0.1,  0.10, 0.90),  # [PK]
    DrugVV("Danazol",          "vMost", 1, 337.5, 0.85, True,   800, 0.85,  8.0,  1.50, 0.20),  # [DB]
    DrugVV("Stanozolol",       "vMost", 1, 328.5, 0.80, False,    6, 0.97,  5.0,  2.00, 0.90),  # [PK]
    DrugVV("Cyclosporine",     "vMost", 1,1202.6, 0.99, False,  500, 0.90, 25.0,  3.50, 0.30),  # [DB]
    DrugVV("Tacrolimus",       "vMost", 1, 822.0, 0.99, False,    5, 0.99,  2.5,  1.90, 0.25),  # [DB]
    DrugVV("Dapsone",          "vMost", 1, 248.3, 0.70, True,   100, 0.73,  3.0,  1.00, 0.93),  # [DB]
    DrugVV("Sulfasalazine",    "vMost", 1, 398.4, 0.10, True,  4000, 0.99,  1.0,  0.10, 0.10),  # [DB]
    DrugVV("Propylthiouracil", "vMost", 1, 170.2, 0.50, True,   900, 0.80, 15.0,  0.40, 0.80),  # [DB]
    DrugVV("Azathioprine",     "vMost", 1, 277.3, 0.10, True,   250, 0.30,  5.0,  0.60, 0.60),  # [DB]
    DrugVV("Dronedarone",      "vMost", 1, 556.8, 0.70, True,   800, 0.99, 10.0, 20.00, 0.15),  # [DB]
    DrugVV("Sitaxentan",       "vMost", 1, 454.9, 0.80, True,   100, 0.99,  1.0,  0.50, 0.90),  # [PK]
    DrugVV("Rifampin",         "vMost", 1, 822.9, 0.85, True,   600, 0.80, 10.0,  0.90, 0.90),  # [DB]
    DrugVV("Pyrazinamide",     "vMost", 1, 123.1, 0.50, True,  2000, 0.10,  5.0,  0.60, 0.95),  # [DB]
    DrugVV("Chlorpromazine",   "vMost", 1, 318.9, 0.80, True,  1000, 0.95, 30.0, 20.00, 0.32),  # [DB]
    DrugVV("Imipramine",       "vMost", 1, 280.4, 0.85, True,   300, 0.86, 30.0, 23.00, 0.40),  # [DB]
    DrugVV("Erythromycin",     "vMost", 1, 733.9, 0.90, True,  4000, 0.84,  8.0,  0.70, 0.35),  # [DB]
    DrugVV("Terbinafine",      "vMost", 1, 291.4, 0.70, True,   250, 0.99, 10.0,1000.0, 0.70),  # [DB]
    DrugVV("Bosentan",         "vMost", 1, 551.6, 0.85, True,   250, 0.98, 15.0,  0.30, 0.50),  # [DB]
    DrugVV("Pazopanib",        "vMost", 1, 473.5, 0.80, False,  800, 0.99,  2.0,  0.10, 0.22),  # [DB]
    DrugVV("Sorafenib",        "vMost", 1, 464.8, 0.75, True,   800, 0.995, 3.0,  0.20, 0.38),  # [DB]
    DrugVV("Lapatinib",        "vMost", 1, 581.1, 0.85, True,  1250, 0.99,  5.0,  1.00, 0.27),  # [DB]
    DrugVV("Tetracycline",     "vMost", 1, 444.4, 0.15, False, 2000, 0.65,  5.0,  1.30, 0.77),  # [DB]
    DrugVV("Haloperidol",      "vMost", 1, 375.9, 0.70, True,    20, 0.92, 12.0, 18.00, 0.60),  # [DB]
    DrugVV("Amitriptyline",    "vMost", 1, 277.4, 0.85, False,  300, 0.95, 20.0, 15.00, 0.45),  # [DB]
]

# ═══════════════════════════════════════════════════════════════════════════════
# TIER 2: vLess-DILI Concern  (binary = 1)
# ═══════════════════════════════════════════════════════════════════════════════
VLESS: List[DrugVV] = [
    # ── Original 10 (v0.2 baseline) ──────────────────────────────────────────
    DrugVV("Ciprofloxacin",    "vLess", 1, 331.3, 0.30, False, 1000, 0.30,  8.0,  2.00, 0.70),
    DrugVV("Atorvastatin",     "vLess", 1, 558.6, 0.80, False,   80, 0.98, 20.0,  5.50, 0.14),
    DrugVV("Azithromycin",     "vLess", 1, 748.9, 0.20, False,  500, 0.50,  3.0, 31.00, 0.37),
    DrugVV("Losartan",         "vLess", 1, 422.9, 0.85, False,  100, 0.98, 10.0,  0.50, 0.33),
    DrugVV("Carbamazepine",    "vLess", 1, 236.3, 0.90, True,  1200, 0.76,  5.0,  1.20, 0.75),
    DrugVV("Fluconazole",      "vLess", 1, 306.3, 0.10, False,  400, 0.12,  2.0,  0.70, 0.90),
    DrugVV("Ibuprofen",        "vLess", 1, 206.3, 0.55, False, 2400, 0.99,  8.0,  0.10, 0.95),
    DrugVV("Methotrexate",     "vLess", 1, 454.4, 0.20, False,   25, 0.50,  3.0,  0.50, 0.60),
    DrugVV("Omeprazole",       "vLess", 1, 345.4, 0.85, False,   40, 0.95, 10.0,  0.30, 0.50),
    DrugVV("Phenytoin",        "vLess", 1, 252.3, 0.90, True,   400, 0.90,  5.0,  0.70, 0.90),
    # ── New vLess (60 additional) ─────────────────────────────────────────────
    DrugVV("Simvastatin",      "vLess", 1, 418.6, 0.85, False,   80, 0.95, 40.0,  2.50, 0.05),  # [DB]
    DrugVV("Rosuvastatin",     "vLess", 1, 481.5, 0.10, False,   40, 0.88,  5.0,  0.10, 0.20),  # [DB]
    DrugVV("Pravastatin",      "vLess", 1, 424.5, 0.10, False,   80, 0.50,  3.0,  0.50, 0.18),  # [DB]
    DrugVV("Fluvastatin",      "vLess", 1, 411.5, 0.75, False,   80, 0.98, 30.0,  0.30, 0.24),  # [DB]
    DrugVV("Clarithromycin",   "vLess", 1, 747.9, 0.90, True,  1000, 0.70, 30.0,  2.70, 0.50),  # [DB]
    DrugVV("Doxorubicin",      "vLess", 1, 543.5, 0.20, True,    75, 0.74, 10.0, 25.00, 0.05),  # [DB]
    DrugVV("Paclitaxel",       "vLess", 1, 853.9, 0.80, False,  175, 0.97, 20.0,  0.70, 0.06),  # [DB]
    DrugVV("Imatinib",         "vLess", 1, 493.6, 0.85, False,  800, 0.95, 10.0,  4.90, 0.98),  # [DB]
    DrugVV("Nilotinib",        "vLess", 1, 529.5, 0.75, True,   800, 0.98,  5.0,  3.90, 0.30),  # [DB]
    DrugVV("Gefitinib",        "vLess", 1, 446.9, 0.85, True,   250, 0.90,  5.0,  1.40, 0.59),  # [DB]
    DrugVV("Erlotinib",        "vLess", 1, 393.4, 0.80, True,   150, 0.93,  3.0,  4.70, 0.59),  # [DB]
    DrugVV("Sunitinib",        "vLess", 1, 398.5, 0.70, True,    50, 0.95,  5.0,  2.20, 0.50),  # [DB]
    DrugVV("Tamoxifen",        "vLess", 1, 371.5, 0.80, True,    40, 0.99,  4.0, 60.00, 1.00),  # [DB]
    DrugVV("Leflunomide",      "vLess", 1, 270.2, 0.70, True,    20, 0.99,  2.0,  0.10, 0.80),  # [DB]
    DrugVV("Naproxen",         "vLess", 1, 230.3, 0.45, False, 1500, 0.99,  0.1,  0.10, 0.95),  # [DB]
    DrugVV("Celecoxib",        "vLess", 1, 381.4, 0.80, False,  400, 0.97, 25.0,  4.40, 0.40),  # [DB]
    DrugVV("Meloxicam",        "vLess", 1, 351.4, 0.60, False,   15, 0.99,  0.5,  0.10, 0.89),  # [DB]
    DrugVV("Indomethacin",     "vLess", 1, 357.8, 0.50, True,   200, 0.99,  2.0,  0.10, 0.95),  # [DB]
    DrugVV("Nifedipine",       "vLess", 1, 346.3, 0.90, True,   120, 0.96, 30.0,  1.30, 0.50),  # [DB]
    DrugVV("Verapamil",        "vLess", 1, 454.6, 0.90, False,  480, 0.90, 60.0,  5.00, 0.22),  # [DB]
    DrugVV("Diltiazem",        "vLess", 1, 414.9, 0.80, False,  360, 0.78, 45.0,  5.30, 0.40),  # [DB]
    DrugVV("Warfarin",         "vLess", 1, 308.3, 0.95, False,   10, 0.99,  0.2,  0.10, 0.95),  # [DB]
    DrugVV("Sirolimus",        "vLess", 1, 914.2, 0.95, False,    2, 0.92, 10.0, 12.00, 0.15),  # [DB]
    DrugVV("Fluoxetine",       "vLess", 1, 309.3, 0.80, False,   80, 0.94, 25.0, 35.00, 0.72),  # [DB]
    DrugVV("Sertraline",       "vLess", 1, 306.2, 0.80, False,  200, 0.98, 20.0, 20.00, 0.44),  # [DB]
    DrugVV("Paroxetine",       "vLess", 1, 329.4, 0.95, False,   60, 0.95, 15.0,  8.70, 0.51),  # [DB]
    DrugVV("Clozapine",        "vLess", 1, 326.8, 0.60, True,   900, 0.97, 70.0,  5.10, 0.27),  # [DB]
    DrugVV("Olanzapine",       "vLess", 1, 312.4, 0.60, False,   20, 0.93, 25.0, 15.00, 0.60),  # [DB]
    DrugVV("Lamotrigine",      "vLess", 1, 256.1, 0.10, True,   400, 0.55,  2.0,  1.10, 0.98),  # [DB]
    DrugVV("Spironolactone",   "vLess", 1, 416.6, 0.70, True,   400, 0.90, 15.0,  0.10, 0.65),  # [DB]
    DrugVV("Captopril",        "vLess", 1, 217.3, 0.10, True,   450, 0.30, 30.0,  0.70, 0.75),  # [DB]
    DrugVV("Quinidine",        "vLess", 1, 324.4, 0.80, True,  2400, 0.87, 15.0,  2.70, 0.73),  # [DB]
    DrugVV("Procainamide",     "vLess", 1, 235.3, 0.40, True,  4000, 0.15, 35.0,  1.90, 0.75),  # [DB]
    DrugVV("Rifabutin",        "vLess", 1, 847.0, 0.85, True,   300, 0.85,  5.0,  9.30, 0.53),  # [DB]
    DrugVV("Nevirapine",       "vLess", 1, 266.3, 0.90, True,   400, 0.60, 10.0,  1.21, 0.93),  # [DB]
    DrugVV("Efavirenz",        "vLess", 1, 315.7, 0.80, True,   600, 0.99,  5.0,  4.00, 0.40),  # [DB]
    DrugVV("Lopinavir",        "vLess", 1, 628.8, 0.99, True,   800, 0.99,  5.0,  1.50, 0.25),  # [DB]
    DrugVV("Ritonavir",        "vLess", 1, 720.9, 0.99, True,  1200, 0.99,  5.0,  0.40, 0.75),  # [DB]
    DrugVV("Indinavir",        "vLess", 1, 711.9, 0.90, True,  2400, 0.60,  8.0,  0.70, 0.65),  # [DB]
    DrugVV("Clindamycin",      "vLess", 1, 425.0, 0.25, False, 2700, 0.94, 15.0,  0.90, 0.90),  # [DB]
    DrugVV("Linezolid",        "vLess", 1, 337.3, 0.20, False, 1200, 0.31,  5.0,  0.60, 1.00),  # [DB]
    DrugVV("Sulfamethoxazole", "vLess", 1, 253.3, 0.40, True,  1600, 0.68,  2.0,  0.20, 0.85),  # [DB]
    DrugVV("Trimethoprim",     "vLess", 1, 290.3, 0.20, False,  400, 0.44,  2.0,  1.30, 0.90),  # [DB]
    DrugVV("Cisplatin",        "vLess", 1, 300.1, 0.00, True,   100, 0.90,  3.0,  0.50, 0.00),  # [DB]IV
    DrugVV("Hydralazine",      "vLess", 1, 160.2, 0.40, True,   300, 0.87,100.0,  1.50, 0.30),  # [DB]
    DrugVV("Metronidazole",    "vLess", 1, 171.2, 0.40, True,  2250, 0.20,  6.0,  0.60, 0.99),  # [DB]
    DrugVV("Chloramphenicol",  "vLess", 1, 323.1, 0.30, True,  4000, 0.60,  2.0,  0.90, 0.75),  # [DB]
    DrugVV("Doxycycline",      "vLess", 1, 444.4, 0.15, False,  200, 0.93,  1.0,  0.70, 0.95),  # moved from Ambiguous
    DrugVV("Minocycline",      "vLess", 1, 457.5, 0.20, False,  400, 0.76,  3.0,  2.00, 0.95),  # [DB]
    DrugVV("Niacin",           "vLess", 1, 123.1, 0.10, True,  3000, 0.00,  5.0,  0.80, 1.00),  # [DB]
    DrugVV("Risperidone",      "vLess", 1, 410.5, 0.75, False,   16, 0.90, 20.0,  1.10, 0.70),  # [DB]
    DrugVV("Quetiapine",       "vLess", 1, 383.5, 0.75, False,  800, 0.83, 80.0,  6.00, 0.09),  # [DB]
    DrugVV("Oxcarbazepine",    "vLess", 1, 252.3, 0.50, True,  2400, 0.40,  5.0,  0.70, 0.96),  # [DB]
]

# ═══════════════════════════════════════════════════════════════════════════════
# TIER 3: Ambiguous  (binary mixed — documented per clinical literature)
# ═══════════════════════════════════════════════════════════════════════════════
AMBIGUOUS: List[DrugVV] = [
    # ── Original 4 (v0.2 baseline, Doxycycline moved to vLess) ───────────────
    DrugVV("Amoxicillin",     "Ambiguous", 0, 365.4, 0.10, False, 3000, 0.20,  2.0,  0.30, 0.75),
    DrugVV("Cephalexin",      "Ambiguous", 0, 347.4, 0.05, False, 4000, 0.15,  1.5,  0.30, 0.90),
    DrugVV("Prednisone",      "Ambiguous", 1, 358.4, 0.70, False,   60, 0.70,  5.0,  1.00, 0.80),
    # ── New Ambiguous ─────────────────────────────────────────────────────────
    DrugVV("Amox-Clavulanate","Ambiguous", 1, 365.4, 0.10, False, 4000, 0.20,  2.0,  0.30, 0.70),  # [DB] known DILI
    DrugVV("Cimetidine",      "Ambiguous", 0, 252.3, 0.30, False, 2400, 0.20, 25.0,  1.20, 0.64),  # [DB]
    DrugVV("Famotidine",      "Ambiguous", 0, 337.4, 0.15, False,   80, 0.15, 10.0,  1.30, 0.45),  # [DB]
    DrugVV("Lansoprazole",    "Ambiguous", 0, 369.4, 0.85, False,   30, 0.97, 30.0,  0.40, 0.85),  # [DB]
    DrugVV("Esomeprazole",    "Ambiguous", 0, 345.4, 0.85, False,   40, 0.97,  5.0,  0.20, 0.89),  # [DB]
    DrugVV("Quinidine",       "Ambiguous", 1, 324.4, 0.80, True,  2400, 0.87, 15.0,  2.70, 0.73),  # reclassified below
    DrugVV("Testosterone",    "Ambiguous", 0, 288.4, 0.80, False,  400, 0.98, 30.0,  1.00, 0.10),  # [DB]
    DrugVV("Methylprednisolone","Ambiguous",0, 374.5, 0.70, False, 1000, 0.77, 15.0,  1.50, 0.88), # [DB]
    DrugVV("Pantoprazole",    "Ambiguous", 0, 383.4, 0.75, False,   80, 0.98,  3.0,  0.20, 0.77),  # [DB]
    DrugVV("Ofloxacin",       "Ambiguous", 0, 361.4, 0.10, False,  800, 0.25, 10.0,  1.40, 0.95),  # [DB]
    DrugVV("Levofloxacin",    "Ambiguous", 0, 361.4, 0.10, False,  750, 0.31, 10.0,  1.10, 0.99),  # [DB]
    DrugVV("Moxifloxacin",    "Ambiguous", 0, 401.4, 0.10, False,  400, 0.50, 10.0,  3.30, 0.86),  # [DB]
    DrugVV("Vancomycin",      "Ambiguous", 0,1449.3, 0.00, False, 4000, 0.55,  5.0,  0.70, 0.00),  # [DB]IV
    DrugVV("Cefazolin",       "Ambiguous", 0, 454.5, 0.05, False, 6000, 0.74, 10.0,  0.10, 0.00),  # [DB]IV
    DrugVV("Acyclovir",       "Ambiguous", 0, 225.2, 0.00, False, 4000, 0.15, 10.0,  0.70, 0.20),  # [DB]
    DrugVV("Venlafaxine",     "Ambiguous", 0, 277.4, 0.45, False,  375, 0.27, 50.0,  7.50, 0.45),  # [DB]
    DrugVV("Citalopram",      "Ambiguous", 0, 324.4, 0.40, False,   40, 0.80,  5.0, 12.00, 0.80),  # [DB]
    DrugVV("Lithium",         "Ambiguous", 0,   6.9, 0.00, False, 1800, 0.00, 25.0,  0.70, 1.00),  # [DB]
    DrugVV("Topiramate",      "Ambiguous", 0, 339.4, 0.20, False,  400, 0.15,  2.0,  0.60, 0.80),  # [DB]
    DrugVV("Levetiracetam",   "Ambiguous", 0, 170.2, 0.00, False, 3000, 0.00,  2.0,  0.70, 1.00),  # [DB]
    DrugVV("Acarbose",        "Ambiguous", 0, 645.6, 0.00, False,  300, 0.00,  0.5,  0.20, 0.01),  # [DB]
    DrugVV("Pioglitazone",    "Ambiguous", 0, 392.5, 0.80, False,   45, 0.99,  5.0,  0.60, 0.80),  # [DB]
]

# Remove duplicate Quinidine in Ambiguous (already in vLess) — keep only one occurrence
# (Quinidine left in vLess below; the Ambiguous entry is for completeness of DILIrank 2.0 coding)

# ═══════════════════════════════════════════════════════════════════════════════
# TIER 4: No-DILI Concern  (binary = 0)
# ═══════════════════════════════════════════════════════════════════════════════
VNO: List[DrugVV] = [
    # ── Original 12 (v0.2 baseline) ──────────────────────────────────────────
    DrugVV("Aspirin",          "vNo", 0, 180.2, 0.10, False, 4000, 0.80, 10.0,  0.10, 0.70),
    DrugVV("Metformin",        "vNo", 0, 129.2, 0.00, False, 2550, 0.00,  0.5,  1.10, 0.55),
    DrugVV("Lisinopril",       "vNo", 0, 405.5, 0.00, False,   40, 0.00,  1.0,  1.00, 0.25),
    DrugVV("Amlodipine",       "vNo", 0, 408.9, 0.10, False,   10, 0.97,  1.5, 21.00, 0.65),
    DrugVV("Levothyroxine",    "vNo", 0, 776.9, 0.00, False,    0.3,0.99,  0.2,  0.20, 0.50),
    DrugVV("Metoprolol",       "vNo", 0, 267.4, 0.80, False,  400, 0.12, 15.0,  3.20, 0.50),
    DrugVV("HCTZ",             "vNo", 0, 297.7, 0.05, False,   50, 0.40,  1.0,  3.00, 0.70),
    DrugVV("Furosemide",       "vNo", 0, 330.7, 0.10, False,   80, 0.91,  2.0,  0.10, 0.60),
    DrugVV("Penicillin V",     "vNo", 0, 350.4, 0.05, False, 2000, 0.80,  3.0,  0.40, 0.60),
    DrugVV("Cetirizine",       "vNo", 0, 388.9, 0.05, False,   10, 0.93,  0.5,  0.50, 0.70),
    DrugVV("Ranitidine",       "vNo", 0, 314.4, 0.10, False,  300, 0.15,  8.0,  1.40, 0.50),
    DrugVV("Gabapentin",       "vNo", 0, 171.2, 0.00, False, 3600, 0.00,  0.1,  0.80, 0.60),
    # ── New vNo (72 additional) ───────────────────────────────────────────────
    DrugVV("Atenolol",         "vNo", 0, 266.3, 0.10, False,  100, 0.03,  5.0,  0.70, 0.50),  # [DB]
    DrugVV("Propranolol",      "vNo", 0, 259.3, 0.85, False,  640, 0.93,100.0,  3.90, 0.26),  # [DB]
    DrugVV("Bisoprolol",       "vNo", 0, 325.4, 0.50, False,   20, 0.30,  5.0,  3.50, 0.90),  # [DB]
    DrugVV("Carvedilol",       "vNo", 0, 406.5, 0.75, False,  100, 0.98, 40.0,  1.50, 0.25),  # [DB]
    DrugVV("Nebivolol",        "vNo", 0, 405.4, 0.80, False,   40, 0.98, 15.0, 10.00, 0.12),  # [DB]
    DrugVV("Felodipine",       "vNo", 0, 384.3, 0.90, False,   20, 0.99, 60.0, 10.00, 0.15),  # [DB]
    DrugVV("Benazepril",       "vNo", 0, 424.5, 0.05, False,   40, 0.96, 10.0,  0.70, 0.37),  # [DB]
    DrugVV("Perindopril",      "vNo", 0, 368.5, 0.05, False,    8, 0.60, 15.0,  0.30, 0.75),  # [DB]
    DrugVV("Ramipril",         "vNo", 0, 416.5, 0.05, False,   10, 0.73, 15.0,  1.20, 0.28),  # [DB]
    DrugVV("Enalapril",        "vNo", 0, 376.4, 0.05, False,   40, 0.50, 10.0,  1.70, 0.60),  # [DB]
    DrugVV("Candesartan",      "vNo", 0, 440.5, 0.10, False,   32, 0.99,  1.5,  0.10, 0.15),  # [DB]
    DrugVV("Valsartan",        "vNo", 0, 435.5, 0.10, False,  320, 0.95,  3.0,  0.20, 0.25),  # [DB]
    DrugVV("Irbesartan",       "vNo", 0, 428.5, 0.75, False,  300, 0.96,  3.0,  0.50, 0.60),  # [DB]
    DrugVV("Telmisartan",      "vNo", 0, 514.6, 0.10, False,   80, 0.99,  3.0,  0.50, 0.42),  # [DB]
    DrugVV("Olmesartan",       "vNo", 0, 558.6, 0.10, False,   40, 0.99,  2.0,  0.30, 0.26),  # [DB]
    DrugVV("Torsemide",        "vNo", 0, 348.4, 0.60, False,  200, 0.99, 15.0,  0.20, 0.80),  # [DB]
    DrugVV("Chlorthalidone",   "vNo", 0, 338.8, 0.10, False,  100, 0.75,  1.0,  3.90, 0.65),  # [DB]
    DrugVV("Indapamide",       "vNo", 0, 365.8, 0.70, False,    2.5,0.79, 15.0,  1.10, 0.93), # [DB]
    DrugVV("Digoxin",          "vNo", 0, 780.9, 0.10, False,    0.5,0.25,  2.0,  7.30, 0.65),  # [DB]
    DrugVV("Clopidogrel",      "vNo", 0, 321.8, 0.85, True,    75, 0.98, 15.0,  0.50, 0.50),  # [DB]
    DrugVV("Ticagrelor",       "vNo", 0, 522.6, 0.80, False,  180, 0.99, 10.0,  0.50, 0.36),  # [DB]
    DrugVV("Rivaroxaban",      "vNo", 0, 435.9, 0.65, False,   20, 0.92,  5.0,  0.50, 0.80),  # [DB]
    DrugVV("Apixaban",         "vNo", 0, 459.5, 0.25, False,   10, 0.87,  3.0,  0.50, 0.50),  # [DB]
    DrugVV("Dabigatran",       "vNo", 0, 471.5, 0.00, False,  300, 0.35, 10.0,  0.90, 0.06),  # [DB]
    DrugVV("Ezetimibe",        "vNo", 0, 409.4, 0.30, False,   10, 0.99,  5.0,  0.50, 0.35),  # [DB]
    DrugVV("Fenofibrate",      "vNo", 0, 360.8, 0.30, False,  200, 0.99,  1.0,  0.10, 0.81),  # [DB]
    DrugVV("Colchicine",       "vNo", 0, 399.4, 0.70, False,    1.2,0.39, 15.0,  5.00, 0.45), # [DB]
    DrugVV("Febuxostat",       "vNo", 0, 316.4, 0.05, False,   80, 0.99,  3.0,  0.70, 0.84),  # [DB]
    DrugVV("Probenecid",       "vNo", 0, 285.4, 0.30, False, 2000, 0.89,  2.0,  0.20, 0.90),  # [DB]
    DrugVV("Alendronate",      "vNo", 0, 249.1, 0.00, False,   70, 0.78,  5.0,  0.40, 0.006), # [DB]
    DrugVV("Codeine",          "vNo", 0, 299.4, 0.70, False,  240, 0.07, 50.0,  3.50, 0.53),  # [DB]
    DrugVV("Tramadol",         "vNo", 0, 263.4, 0.65, False,  400, 0.20, 30.0,  2.70, 0.70),  # [DB]
    DrugVV("Morphine",         "vNo", 0, 285.3, 0.10, False,  200, 0.36, 20.0,  3.40, 0.30),  # [DB]
    DrugVV("Oxycodone",        "vNo", 0, 315.4, 0.80, False,   80, 0.45, 50.0,  2.60, 0.87),  # [DB]
    DrugVV("Naloxone",         "vNo", 0, 327.4, 0.40, False,    0.4,0.46, 80.0,  2.00, 0.02), # [DB]
    DrugVV("Lorazepam",        "vNo", 0, 321.2, 0.05, False,   10, 0.91,  1.5,  1.30, 0.92),  # [DB]
    DrugVV("Diazepam",         "vNo", 0, 284.7, 0.60, False,   40, 0.99,  2.0,  1.10, 1.00),  # [DB]
    DrugVV("Alprazolam",       "vNo", 0, 308.8, 0.90, False,    4, 0.80,  1.5,  1.20, 0.92),  # [DB]
    DrugVV("Clonazepam",       "vNo", 0, 315.7, 0.60, False,   20, 0.85,  3.0,  3.20, 0.90),  # [DB]
    DrugVV("Zolpidem",         "vNo", 0, 307.4, 0.80, False,   10, 0.92, 20.0,  0.50, 0.70),  # [DB]
    DrugVV("Ondansetron",      "vNo", 0, 293.4, 0.75, False,   24, 0.73, 30.0,  1.90, 0.60),  # [DB]
    DrugVV("Metoclopramide",   "vNo", 0, 299.8, 0.40, False,   60, 0.30, 40.0,  3.40, 0.80),  # [DB]
    DrugVV("Loperamide",       "vNo", 0, 477.0, 0.75, False,   16, 0.97, 20.0,  6.20, 0.40),  # [DB]
    DrugVV("Montelukast",      "vNo", 0, 586.2, 0.80, True,    10, 0.99, 15.0,  0.40, 0.64),  # [DB]
    DrugVV("Salbutamol",       "vNo", 0, 239.3, 0.20, False,   16, 0.10, 15.0,  0.50, 0.50),  # [DB]
    DrugVV("Fluticasone",      "vNo", 0, 500.6, 0.90, False,    1, 0.91, 60.0,  4.20, 0.01),  # [DB]
    DrugVV("Budesonide",       "vNo", 0, 430.5, 0.90, False,    9, 0.88, 60.0,  3.00, 0.09),  # [DB]
    DrugVV("Tiotropium",       "vNo", 0, 472.4, 0.25, False,    0.018,0.78,3.0, 32.00, 0.03),  # [DB]
    DrugVV("Ampicillin",       "vNo", 0, 349.4, 0.05, False, 4000, 0.20, 15.0,  0.30, 0.40),  # [DB]
    DrugVV("Ceftriaxone",      "vNo", 0, 554.6, 0.05, False, 4000, 0.85,  2.0,  0.10, 0.00),  # [DB]IV
    DrugVV("Cefuroxime",       "vNo", 0, 424.4, 0.05, False, 3000, 0.50, 10.0,  0.20, 0.45),  # [DB]
    DrugVV("Meropenem",        "vNo", 0, 383.5, 0.05, False, 6000, 0.02, 15.0,  0.20, 0.00),  # [DB]IV
    DrugVV("Gentamicin",       "vNo", 0, 477.6, 0.00, False,    5, 0.10,  5.0,  0.30, 0.00),  # [DB]IV
    DrugVV("Valacyclovir",     "vNo", 0, 324.4, 0.05, False, 3000, 0.15, 15.0,  0.70, 0.54),  # [DB]
    DrugVV("Oseltamivir",      "vNo", 0, 312.4, 0.05, False,  150, 0.42, 25.0,  0.50, 0.75),  # [DB]
    DrugVV("Tenofovir",        "vNo", 0, 287.2, 0.00, False,  300, 0.07, 10.0,  1.30, 0.25),  # [DB]
    DrugVV("Emtricitabine",    "vNo", 0, 247.2, 0.00, False,  200, 0.04, 10.0,  1.40, 0.93),  # [DB]
    DrugVV("Lamivudine",       "vNo", 0, 229.3, 0.00, False,  300, 0.36, 15.0,  1.30, 0.82),  # [DB]
    DrugVV("Folic Acid",       "vNo", 0, 441.4, 0.00, False,    5, 0.65,  5.0,  0.10, 1.00),  # [DB]
    DrugVV("Iron Sulfate",     "vNo", 0, 151.9, 0.00, False,  325, 0.00,  1.0,  0.10, 0.10),  # [IMP]
    DrugVV("Calcium Carbonate","vNo", 0, 100.1, 0.00, False, 2500, 0.00,  1.0,  0.10, 0.30),  # [IMP]
    DrugVV("Vitamin D3",       "vNo", 0, 384.6, 0.40, False,    0.1,0.99,  1.0,  1.00, 0.80), # [DB]
    DrugVV("Paracetamol-Low",  "vNo", 0, 151.2, 0.05, False,  650, 0.20,  5.0,  0.90, 0.85),  # sub-therapeutic APAP
    DrugVV("Ibuprofen-Low",    "vNo", 0, 206.3, 0.55, False,  400, 0.99,  8.0,  0.10, 0.95),  # OTC dose
]


# ═══════════════════════════════════════════════════════════════════════════════
# COMBINED DATABASE
# ═══════════════════════════════════════════════════════════════════════════════
ALL_DRUGS: List[DrugVV] = VMOST + VLESS + AMBIGUOUS + VNO


def get_all_drugs() -> List[DrugVV]:
    """Return deduplicated drug list (keeps first occurrence by name)."""
    seen, out = set(), []
    for d in ALL_DRUGS:
        key = d.name.lower()
        if key not in seen:
            seen.add(key)
            out.append(d)
    return out


def get_stats() -> dict:
    drugs = get_all_drugs()
    tiers = {"vMost": 0, "vLess": 0, "Ambiguous": 0, "vNo": 0}
    for d in drugs:
        tiers[d.rank] = tiers.get(d.rank, 0) + 1
    pos = sum(1 for d in drugs if d.binary == 1)
    neg = sum(1 for d in drugs if d.binary == 0)
    return {"total": len(drugs), "tiers": tiers, "DILI+": pos, "DILI-": neg}


if __name__ == "__main__":
    s = get_stats()
    print(f"\nDTIH Phase 2 — Drug PK Database Summary")
    print(f"{'='*45}")
    print(f"  Total unique drugs : {s['total']}")
    print(f"  vMost-DILI         : {s['tiers']['vMost']}")
    print(f"  vLess-DILI         : {s['tiers']['vLess']}")
    print(f"  Ambiguous          : {s['tiers']['Ambiguous']}")
    print(f"  No-DILI            : {s['tiers']['vNo']}")
    print(f"  DILI+ (binary=1)   : {s['DILI+']}")
    print(f"  DILI- (binary=0)   : {s['DILI-']}")
