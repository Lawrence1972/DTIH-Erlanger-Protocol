"""
DTIH Full Pipeline: L1→L2→L3→L4→L5 통합 시뮬레이션
Digital Twin In-silico Human — Erlanger Protocol v0.2

변경사항 (v0.1 → v0.2):
  - L4 신경-면역 조절층 코드 구현 추가
  - L5 피드백 루프 (FB1: 간손상→효소감소, FB2: 면역→혈류변화) 구현
  - V&V 파이프라인 통합
  - Windows 패키지 지원
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json, os, sys
import warnings
warnings.filterwarnings('ignore')


# ═══════════════════════════════════════════════════════════════
# 데이터 구조
# ═══════════════════════════════════════════════════════════════

@dataclass
class PatientProfile:
    name: str
    phenotype: str
    genotype_modifier: float = 1.0
    # L4 파라미터
    nerve_sensitivity: float = 1.0   # 신경 민감도 (1.0=정상)
    immune_reactivity: float = 1.0   # 면역 반응성 (1.0=정상)
    baseline_gsh: float = 1.0        # 기저 GSH 수준 (1.0=정상)

@dataclass
class DrugInput:
    name: str
    formula: str
    mw: float
    dose_mg: float
    route: str  # "oral" or "iv"
    # 대사 경로 비율
    ugt_fraction: float = 0.55
    sult_fraction: float = 0.30
    cyp_fraction: float = 0.15
    # 독성 특성
    reactive_metabolite: bool = True
    herg_ic50: float = 1000.0       # hERG IC50 (μM) — 심독성
    nerve_ec50: float = 500.0       # 신경전도 억제 EC50 (μM)
    immune_trigger: float = 0.3     # 면역 반응 유발 임계값 (정규화)
    # PK 파라미터
    ka: float = 1.5
    bioavailability: float = 0.85
    protein_binding: float = 0.20
    # PBPK
    V_blood: float = 5.0
    V_liver: float = 1.5
    V_kidney: float = 0.3
    V_heart: float = 0.3
    Q_liver: float = 90.0
    Q_kidney: float = 60.0
    Q_heart: float = 15.0
    Kp_liver: float = 2.5
    Kp_kidney: float = 1.8
    Kp_heart: float = 1.2

DRUG_LIBRARY = {
    "acetaminophen": DrugInput(
        name="Acetaminophen", formula="C₈H₉NO₂", mw=151.16,
        dose_mg=1000, route="oral",
        ugt_fraction=0.55, sult_fraction=0.30, cyp_fraction=0.15,
        reactive_metabolite=True, herg_ic50=2000, nerve_ec50=800,
        immune_trigger=0.4, protein_binding=0.20,
    ),
    "acetaminophen_od": DrugInput(
        name="Acetaminophen (과량)", formula="C₈H₉NO₂", mw=151.16,
        dose_mg=10000, route="oral",
        ugt_fraction=0.55, sult_fraction=0.30, cyp_fraction=0.15,
        reactive_metabolite=True, herg_ic50=2000, nerve_ec50=800,
        immune_trigger=0.4, protein_binding=0.20,
    ),
    "cisplatin": DrugInput(
        name="Cisplatin", formula="Cl₂H₆N₂Pt", mw=300.05,
        dose_mg=100, route="iv",
        ugt_fraction=0.10, sult_fraction=0.05, cyp_fraction=0.85,
        reactive_metabolite=True, herg_ic50=50, nerve_ec50=30,
        immune_trigger=0.15, ka=0, bioavailability=1.0, protein_binding=0.90,
        Kp_kidney=4.0,
    ),
    "doxorubicin": DrugInput(
        name="Doxorubicin", formula="C₂₇H₂₉NO₁₁", mw=543.52,
        dose_mg=60, route="iv",
        ugt_fraction=0.20, sult_fraction=0.10, cyp_fraction=0.70,
        reactive_metabolite=True, herg_ic50=10, nerve_ec50=200,
        immune_trigger=0.20, ka=0, bioavailability=1.0, protein_binding=0.75,
        Kp_heart=5.0,
    ),
}

PATIENT_LIBRARY = {
    "normal": PatientProfile("정상인 (EM)", "Normal Metabolizer", 1.0, 1.0, 1.0, 1.0),
    "cyp_inducer": PatientProfile("CYP유도자 (알코올)", "CYP2E1 Inducer", 1.0, 1.2, 1.0, 0.8),
    "ugt_poor": PatientProfile("UGT저활성 (PM)", "UGT1A Poor Metabolizer", 1.3, 1.0, 1.0, 1.0),
    "immune_high": PatientProfile("면역과반응", "Immune Hyperreactive", 1.0, 1.0, 2.0, 1.0),
    "gsh_depleted": PatientProfile("GSH저하 (영양불량)", "GSH Depleted", 1.1, 1.5, 1.2, 0.4),
    "elderly": PatientProfile("고령자 (75세)", "Elderly", 0.7, 1.3, 1.5, 0.7),
}


# ═══════════════════════════════════════════════════════════════
# L1: 유전체 → 효소 제약
# ═══════════════════════════════════════════════════════════════

def layer1_enzyme_constraints(patient: PatientProfile, drug: DrugInput) -> Dict[str, float]:
    """
    Def 1.1-1.4: 유전형 기반 효소 활성 조정
    Returns: 경로별 최대 대사 용량 비율
    """
    mod = patient.genotype_modifier
    return {
        'ugt_capacity': drug.ugt_fraction / max(mod, 0.3),
        'sult_capacity': drug.sult_fraction / max(mod, 0.3),
        'cyp_capacity': drug.cyp_fraction * max(mod, 0.5),
        'total_cl_modifier': 1.0 / max(mod, 0.3),
    }


# ═══════════════════════════════════════════════════════════════
# L2: 대사 네트워크 → 경로별 플럭스 분배
# ═══════════════════════════════════════════════════════════════

def layer2_flux_distribution(enzyme_constraints: Dict, drug: DrugInput) -> Dict[str, float]:
    """
    Def 2.1-2.6: 효소 제약 아래에서 대사 플럭스 분배 계산
    FBA 축약: 효소 용량에 비례하여 3경로로 분배
    """
    ugt = enzyme_constraints['ugt_capacity']
    sult = enzyme_constraints['sult_capacity']
    cyp = enzyme_constraints['cyp_capacity']
    total = ugt + sult + cyp

    return {
        'ugt_fraction': ugt / total,
        'sult_fraction': sult / total,
        'cyp_fraction': cyp / total,   # NAPQI 생성 경로
        'cl_modifier': enzyme_constraints['total_cl_modifier'],
        'base_cl_int': 15.0 * enzyme_constraints['total_cl_modifier'],
    }


# ═══════════════════════════════════════════════════════════════
# L3+L4+L5: 통합 PBPK-신경-면역 ODE 시스템
# ═══════════════════════════════════════════════════════════════

def run_integrated_simulation(
    drug: DrugInput,
    patient: PatientProfile,
    flux: Dict[str, float],
    t_end: float = 48.0,
) -> Dict[str, np.ndarray]:
    """
    L3 (PBPK) + L4 (신경-면역) + L5 (피드백) 통합 ODE 시스템

    상태 변수 (11개):
      y[0]  A_gut       : 장 내 약물량 (μmol)
      y[1]  C_blood     : 혈중 농도 (μM)
      y[2]  C_liver     : 간 농도 (μM)
      y[3]  C_kidney    : 신장 농도 (μM)
      y[4]  C_heart     : 심장 농도 (μM)
      y[5]  GSH         : 간 GSH 수준 (0~1, 정규화)
      y[6]  liver_damage: 간 손상 누적 (0~1)
      y[7]  theta_C     : C섬유 전도속도 비율 (0~1, Def 4.1)
      y[8]  IL6         : IL-6 수준 (정규화, Def 4.2)
      y[9]  TNFa        : TNF-α 수준 (정규화, Def 4.2)
      y[10] IL10        : IL-10 수준 (정규화, Def 4.2 — 항염증)
    """

    mw = drug.mw
    dose_umol = (drug.dose_mg / mw) * 1000
    napqi_frac = flux['cyp_fraction']
    cl_int = flux['base_cl_int']
    fu = 1.0 - drug.protein_binding

    # Well-stirred 간 청소율
    CL_hep_base = (drug.Q_liver * fu * cl_int) / (drug.Q_liver + fu * cl_int)
    # 신장 청소율 (단순화)
    CL_renal = 5.0 * fu

    def ode_system(t, y):
        A_gut, C_blood, C_liver, C_kidney, C_heart, GSH, \
            liver_dmg, theta_C, IL6, TNFa, IL10 = y

        # ─── L3: PBPK 구획 모델 (Def 3.1) ───

        # L5-FB1: 간 손상 → 효소 활성 감소 (Def 5.3 φ₁)
        enzyme_activity = max(0.1, 1.0 - liver_dmg * 0.8)
        CL_hep = CL_hep_base * enzyme_activity

        # L5-FB2: 염증 → 혈류 변화 (Def 5.3 φ₂)
        inflammation_factor = 1.0 + 0.3 * (IL6 + TNFa) - 0.1 * IL10
        Q_liver_eff = drug.Q_liver * max(0.5, min(1.5, inflammation_factor))

        # 장 흡수
        if drug.route == "oral":
            dA_gut = -drug.ka * A_gut
            absorption = drug.ka * A_gut * drug.bioavailability / drug.V_blood
        else:
            dA_gut = 0
            absorption = 0

        # 간 구획
        liver_uptake = Q_liver_eff / drug.V_liver * (C_blood - C_liver / drug.Kp_liver)
        liver_metab = CL_hep * C_liver / drug.V_liver
        dC_liver = liver_uptake - liver_metab

        # 신장 구획
        kidney_uptake = drug.Q_kidney / drug.V_kidney * (C_blood - C_kidney / drug.Kp_kidney)
        kidney_clear = CL_renal * C_kidney / drug.V_kidney
        dC_kidney = kidney_uptake - kidney_clear

        # 심장 구획
        heart_uptake = drug.Q_heart / drug.V_heart * (C_blood - C_heart / drug.Kp_heart)
        dC_heart = heart_uptake

        # 혈중 농도
        dC_blood = (absorption
                    - Q_liver_eff / drug.V_blood * C_blood
                    + Q_liver_eff / drug.V_blood * C_liver / drug.Kp_liver
                    - drug.Q_kidney / drug.V_blood * C_blood
                    + drug.Q_kidney / drug.V_blood * C_kidney / drug.Kp_kidney
                    - drug.Q_heart / drug.V_blood * C_blood
                    + drug.Q_heart / drug.V_blood * C_heart / drug.Kp_heart)

        # ─── GSH 동역학 및 간 손상 ───
        napqi_prod = liver_metab * drug.V_liver * napqi_frac
        gsh_detox_capacity = GSH * 300 * patient.baseline_gsh
        gsh_consumed = min(napqi_prod, gsh_detox_capacity)
        napqi_free = max(0, napqi_prod - gsh_consumed)

        gsh_regen = 0.08 * max(0, 1.0 - GSH)
        gsh_loss = gsh_consumed * 0.0008
        dGSH = gsh_regen - gsh_loss
        if GSH <= 0.01 and dGSH < 0:
            dGSH = 0

        # 간 손상: free NAPQI → 단백질 결합
        damage_rate = napqi_free * 0.0003 * (1.0 if drug.reactive_metabolite else 0.1)
        repair_rate = 0.005 * max(0, 1.0 - liver_dmg)
        dLiver_dmg = damage_rate - repair_rate
        if liver_dmg >= 0.99 and dLiver_dmg > 0:
            dLiver_dmg = 0

        # ─── L4a: 신경 전도 모델 (Def 4.1) ───
        # Erlanger-Gasser C섬유 전도속도 억제 (Hill 방정식)
        C_nerve = C_blood  # 약물이 혈류를 통해 신경에 도달
        n_hill = 2.0
        eta_C = (C_nerve ** n_hill) / (drug.nerve_ec50 ** n_hill + C_nerve ** n_hill)
        eta_C *= patient.nerve_sensitivity
        theta_target = max(0, 1.0 - eta_C)
        # 시간 지연 포함 (신경 반응 시간상수 ~2h)
        tau_nerve = 2.0
        dTheta_C = (theta_target - theta_C) / tau_nerve

        # ─── L4b: 면역 반응 모델 (Def 4.2) ───
        # 사이토카인 캐스케이드 ODE
        # 트리거: 간 손상 + NAPQI + 약물 직접 자극
        immune_stimulus = (
            liver_dmg * 2.0 +                              # 간세포 괴사 → DAMP 방출
            napqi_free * 0.001 +                            # NAPQI 직접 면역 자극
            max(0, C_blood / 1000 - drug.immune_trigger)   # 약물 직접 면역 반응
        ) * patient.immune_reactivity

        # IL-6 (전염증)
        dIL6 = (
            0.5 * immune_stimulus          # 외부 자극
            + 0.3 * TNFa                   # TNF-α가 IL-6 유도 (k_jk)
            - 0.2 * IL10 * IL6             # IL-10이 IL-6 억제
            - 0.15 * IL6                   # 자연 분해 (γ_k)
        )

        # TNF-α (전염증)
        dTNFa = (
            0.4 * immune_stimulus
            + 0.1 * IL6                    # 양성 피드백
            - 0.25 * IL10 * TNFa           # IL-10 억제
            - 0.2 * TNFa                   # 자연 분해
        )

        # IL-10 (항염증 — 네거티브 피드백)
        dIL10 = (
            0.1 * immune_stimulus
            + 0.2 * IL6                    # IL-6가 IL-10 유도
            + 0.15 * TNFa                  # TNF-α가 IL-10 유도
            - 0.1 * IL10                   # 자연 분해
        )

        # 음수 방지
        if IL6 <= 0 and dIL6 < 0: dIL6 = 0
        if TNFa <= 0 and dTNFa < 0: dTNFa = 0
        if IL10 <= 0 and dIL10 < 0: dIL10 = 0

        return [dA_gut, dC_blood, dC_liver, dC_kidney, dC_heart,
                dGSH, dLiver_dmg, dTheta_C, dIL6, dTNFa, dIL10]

    # 초기 조건
    if drug.route == "oral":
        y0 = [dose_umol, 0, 0, 0, 0,
              patient.baseline_gsh, 0, 1.0, 0, 0, 0]
    else:  # IV
        C0_blood = dose_umol / drug.V_blood
        y0 = [0, C0_blood, 0, 0, 0,
              patient.baseline_gsh, 0, 1.0, 0, 0, 0]

    t_eval = np.linspace(0, t_end, 800)
    sol = solve_ivp(ode_system, (0, t_end), y0, t_eval=t_eval,
                    method='RK45', max_step=0.05, rtol=1e-8, atol=1e-10)

    # ─── 파생 지표 계산 ───
    # Def 3.3: QT 위험 (hERG 억제)
    qt_risk = sol.y[4] / drug.herg_ic50  # C_heart / IC50

    # DILI 위험 종합 (간손상 + GSH 고갈 + 면역 반응)
    dili_composite = (
        sol.y[6] * 40 +                          # 간 손상 (0~40)
        np.maximum(0, 1 - sol.y[5]) * 30 +       # GSH 고갈 (0~30)
        (sol.y[8] + sol.y[9]) * 5                 # 염증 (0~30)
    )

    # 신경독성 지표 (전도속도 감소율)
    neuro_tox = (1 - sol.y[7]) * 100  # %

    return {
        't': sol.t,
        'C_blood': sol.y[1],
        'C_liver': sol.y[2],
        'C_kidney': sol.y[3],
        'C_heart': sol.y[4],
        'GSH': sol.y[5] * 100,
        'liver_damage': sol.y[6] * 100,
        'theta_C': sol.y[7] * 100,
        'IL6': sol.y[8],
        'TNFa': sol.y[9],
        'IL10': sol.y[10],
        'QT_risk': qt_risk,
        'DILI_composite': dili_composite,
        'neuro_tox': neuro_tox,
    }


# ═══════════════════════════════════════════════════════════════
# 시각화
# ═══════════════════════════════════════════════════════════════

def plot_full_results(results: Dict[str, Dict], drug_name: str, save_path: str):
    """6-panel 종합 결과 시각화"""
    fig, axes = plt.subplots(3, 2, figsize=(16, 14))
    fig.suptitle(f'DTIH L1→L2→L3→L4→L5 Integrated Simulation\n{drug_name}',
                 fontsize=15, fontweight='bold')

    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#3B8132', '#6B4C9A']

    for idx, (label, r) in enumerate(results.items()):
        c = colors[idx % len(colors)]
        lw = 2

        # (a) 혈중 농도
        axes[0,0].plot(r['t'], r['C_blood'], c=c, lw=lw, label=label)
        # (b) 간 손상 + GSH
        axes[0,1].plot(r['t'], r['liver_damage'], c=c, lw=lw, label=f'{label} 간손상', linestyle='-')
        axes[0,1].plot(r['t'], r['GSH'], c=c, lw=lw/1.5, linestyle='--', alpha=0.6)
        # (c) 사이토카인
        axes[1,0].plot(r['t'], r['IL6'], c=c, lw=lw, label=f'{label} IL-6')
        axes[1,0].plot(r['t'], r['TNFa'], c=c, lw=lw/1.5, linestyle='--', alpha=0.7)
        # (d) 신경전도
        axes[1,1].plot(r['t'], r['theta_C'], c=c, lw=lw, label=label)
        # (e) 심장 QT 위험
        axes[2,0].plot(r['t'], r['QT_risk'], c=c, lw=lw, label=label)
        # (f) DILI 종합 점수
        axes[2,1].plot(r['t'], r['DILI_composite'], c=c, lw=lw, label=label)

    # 축 설정
    for ax in axes.flat:
        ax.grid(True, alpha=0.2)
        ax.legend(fontsize=8, loc='best')

    axes[0,0].set_ylabel('혈중 농도 (μM)'); axes[0,0].set_title('(a) 혈중 약물 농도 [L3]')
    axes[0,1].set_ylabel('%'); axes[0,1].set_title('(b) 간 손상(실선) / GSH(점선) [L3]')
    axes[0,1].axhline(y=30, color='red', ls=':', alpha=0.5, label='GSH 위험선')
    axes[1,0].set_ylabel('사이토카인 수준'); axes[1,0].set_title('(c) 면역 반응: IL-6(실선), TNF-α(점선) [L4b]')
    axes[1,1].set_ylabel('C섬유 전도속도 (%)'); axes[1,1].set_title('(d) 말초신경 전도속도 [L4a — Erlanger C섬유]')
    axes[1,1].axhline(y=50, color='red', ls=':', alpha=0.5, label='신경독성 임계값')
    axes[2,0].set_ylabel('QT위험비 (C/IC50)'); axes[2,0].set_title('(e) 심장 QT 연장 위험 [L3 Def 3.3]')
    axes[2,0].axhline(y=1.0, color='red', ls=':', alpha=0.5, label='QT 위험 임계값')
    axes[2,1].set_ylabel('DILI 종합 점수'); axes[2,1].set_title('(f) DILI 종합 위험도 [L5 통합]')
    axes[2,1].axhline(y=30, color='orange', ls=':', alpha=0.5)
    axes[2,1].axhline(y=60, color='red', ls=':', alpha=0.5)

    for ax in axes[2]:
        ax.set_xlabel('시간 (hr)')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"  [저장] {save_path}")


def print_summary(results: Dict[str, Dict], drug: DrugInput):
    """콘솔 결과 요약"""
    print(f"\n{'═'*75}")
    print(f"  약물: {drug.name} ({drug.formula}) | 용량: {drug.dose_mg}mg {drug.route}")
    print(f"{'═'*75}")
    print(f"  {'환자':18s} {'Cmax':>8s} {'간손상%':>8s} {'GSH최저%':>8s} "
          f"{'C전도%':>8s} {'IL-6max':>8s} {'QTrisk':>8s} {'DILI':>8s}")
    print(f"  {'─'*68}")
    for label, r in results.items():
        cmax = np.max(r['C_blood'])
        dmg = np.max(r['liver_damage'])
        gsh = np.min(r['GSH'])
        theta = np.min(r['theta_C'])
        il6 = np.max(r['IL6'])
        qt = np.max(r['QT_risk'])
        dili = np.max(r['DILI_composite'])

        flags = ""
        if dmg > 20: flags += "🔴간 "
        if gsh < 30: flags += "🔴GSH "
        if theta < 50: flags += "🔴신경 "
        if qt > 1.0: flags += "🔴심장 "
        if il6 > 5: flags += "🟡면역 "
        if not flags: flags = "✅안전"

        print(f"  {label:18s} {cmax:8.1f} {dmg:8.1f} {gsh:8.1f} "
              f"{theta:8.1f} {il6:8.2f} {qt:8.3f} {dili:8.1f}  {flags}")


# ═══════════════════════════════════════════════════════════════
# 메인 실행
# ═══════════════════════════════════════════════════════════════

def run_scenario(drug_key: str, patient_keys: List[str] = None,
                 output_dir: str = "."):
    """하나의 약물에 대해 여러 환자 유형 시뮬레이션"""
    if drug_key not in DRUG_LIBRARY:
        print(f"  [오류] 약물 '{drug_key}' 없음. 사용 가능: {list(DRUG_LIBRARY.keys())}")
        return

    drug = DRUG_LIBRARY[drug_key]
    if patient_keys is None:
        patient_keys = list(PATIENT_LIBRARY.keys())

    print(f"\n╔{'═'*73}╗")
    print(f"║  DTIH v0.2 — L1→L2→L3→L4→L5 통합 시뮬레이션                       ║")
    print(f"║  Erlanger Protocol · Phase 1 Full Pipeline                        ║")
    print(f"╚{'═'*73}╝")

    all_results = {}
    for pk in patient_keys:
        patient = PATIENT_LIBRARY[pk]
        print(f"\n  ── {patient.name} ({patient.phenotype}) ──")

        # L1
        enzymes = layer1_enzyme_constraints(patient, drug)
        # L2
        flux = layer2_flux_distribution(enzymes, drug)
        print(f"  L2 플럭스: UGT={flux['ugt_fraction']:.1%} SULT={flux['sult_fraction']:.1%} "
              f"CYP(NAPQI)={flux['cyp_fraction']:.1%}")
        # L3+L4+L5
        result = run_integrated_simulation(drug, patient, flux)
        all_results[patient.name] = result

    print_summary(all_results, drug)

    save_path = os.path.join(output_dir, f"dtih_{drug_key}_results.png")
    plot_full_results(all_results, f"{drug.name} {drug.dose_mg}mg {drug.route}", save_path)

    # JSON 결과 저장
    json_data = {}
    for label, r in all_results.items():
        json_data[label] = {
            'Cmax_uM': float(np.max(r['C_blood'])),
            'liver_damage_pct': float(np.max(r['liver_damage'])),
            'GSH_min_pct': float(np.min(r['GSH'])),
            'C_fiber_min_pct': float(np.min(r['theta_C'])),
            'IL6_max': float(np.max(r['IL6'])),
            'TNFa_max': float(np.max(r['TNFa'])),
            'QT_risk_max': float(np.max(r['QT_risk'])),
            'DILI_composite_max': float(np.max(r['DILI_composite'])),
        }
    json_path = os.path.join(output_dir, f"dtih_{drug_key}_summary.json")
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(json_data, f, indent=2, ensure_ascii=False)
    print(f"  [저장] {json_path}")

    return all_results


def main():
    output_dir = os.path.dirname(os.path.abspath(__file__))

    # 시나리오 1: 아세트아미노펜 정상 용량
    run_scenario("acetaminophen", output_dir=output_dir)

    # 시나리오 2: 아세트아미노펜 과량
    run_scenario("acetaminophen_od", output_dir=output_dir)

    # 시나리오 3: 시스플라틴 (신독성 + 신경독성)
    run_scenario("cisplatin", output_dir=output_dir)

    # 시나리오 4: 독소루비신 (심독성)
    run_scenario("doxorubicin", output_dir=output_dir)

    print(f"\n{'═'*75}")
    print(f"  전체 시뮬레이션 완료. 결과 파일이 현재 폴더에 저장되었습니다.")
    print(f"{'═'*75}\n")


if __name__ == "__main__":
    main()
