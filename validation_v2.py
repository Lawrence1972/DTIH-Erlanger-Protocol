"""
DTIH V&V Module — Phase 2
Erlanger Protocol v0.3 | 212-drug DILIrank Expanded Validation
ASME V&V 40 Framework

Usage:
    python validation_v2.py
    python validation_v2.py --tier vMost
    python validation_v2.py --tier vLess
    python validation_v2.py --tier all
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass
import os
import sys

# ── import Phase 2 database ───────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(__file__))
from drug_pk_database import get_all_drugs, get_stats, DrugVV

# ── DILI score (same formula as v0.2 for comparability) ──────────────────────
def calc_dili_score(d: DrugVV) -> float:
    rs = min(1.0, d.cyp * (2.0 if d.reactive else 0.5)) * 40
    ds = min(1.0, np.log10(d.dose / 100) / 1.5) * 25 if d.dose > 100 else 0
    bs = d.binding * (1 - min(1, d.cl / 30)) * 15
    ml = min(1.0, d.dose * d.bio * d.cyp / 500) * 20
    return min(99, max(0, rs + ds + bs + ml))


# ── Core metrics ──────────────────────────────────────────────────────────────
def compute_metrics(drugs, threshold=30):
    scores  = [(d.name, d.rank, d.binary, calc_dili_score(d)) for d in drugs]
    actual  = np.array([s[2] for s in scores])
    sc      = np.array([s[3] for s in scores])
    pred    = (sc >= threshold).astype(int)

    TP = int(np.sum((pred == 1) & (actual == 1)))
    TN = int(np.sum((pred == 0) & (actual == 0)))
    FP = int(np.sum((pred == 1) & (actual == 0)))
    FN = int(np.sum((pred == 0) & (actual == 1)))

    sens = TP / max(TP + FN, 1)
    spec = TN / max(TN + FP, 1)
    acc  = (TP + TN) / len(scores)
    ppv  = TP / max(TP + FP, 1)
    f1   = 2 * ppv * sens / max(ppv + sens, 0.001)

    # ROC / AUC
    thresholds = np.arange(0, 100, 1)
    tpr_l, fpr_l = [], []
    for thr in thresholds:
        p  = (sc >= thr).astype(int)
        tp = np.sum((p == 1) & (actual == 1))
        fp = np.sum((p == 1) & (actual == 0))
        fn = np.sum((p == 0) & (actual == 1))
        tn = np.sum((p == 0) & (actual == 0))
        tpr_l.append(tp / max(tp + fn, 1))
        fpr_l.append(fp / max(fp + tn, 1))
    fpr_a, tpr_a = np.array(fpr_l), np.array(tpr_l)
    si  = np.argsort(fpr_a)
    auc = float(np.abs(np.trapezoid(tpr_a[si], fpr_a[si])))

    return dict(
        scores=scores, actual=actual, sc=sc, pred=pred,
        TP=TP, TN=TN, FP=FP, FN=FN,
        sens=sens, spec=spec, acc=acc, ppv=ppv, f1=f1, auc=auc,
        fpr=fpr_l, tpr=tpr_l
    )


# ── Optimal threshold search ──────────────────────────────────────────────────
def find_optimal_threshold(drugs):
    """Youden's J: maximise sensitivity + specificity - 1."""
    best_j, best_thr = -1, 30
    for thr in range(5, 95):
        m = compute_metrics(drugs, threshold=thr)
        j = m['sens'] + m['spec'] - 1
        if j > best_j:
            best_j, best_thr = j, thr
    return best_thr, best_j


# ── Tier breakdown ────────────────────────────────────────────────────────────
def tier_breakdown(all_drugs):
    tiers = {"vMost": [], "vLess": [], "Ambiguous": [], "vNo": []}
    for d in all_drugs:
        if d.rank in tiers:
            tiers[d.rank].append(d)
    return tiers


# ── Full report ───────────────────────────────────────────────────────────────
def run_vv_phase2(output_dir=".", tier_filter="all"):
    print("\n" + "=" * 65)
    print("  DTIH V&V Phase 2 — DILIrank Expanded Validation (212 drugs)")
    print("  Erlanger Protocol v0.3 | ASME V&V 40")
    print("=" * 65)

    all_drugs = get_all_drugs()

    # ── filter by tier if requested
    if tier_filter != "all":
        all_drugs = [d for d in all_drugs if d.rank == tier_filter]
        if not all_drugs:
            print(f"  [!] No drugs found for tier '{tier_filter}'")
            return None

    print(f"\n  Dataset : {len(all_drugs)} drugs  (tier_filter='{tier_filter}')")
    db_stats = get_stats()
    print(f"  Breakdown: vMost={db_stats['tiers']['vMost']}  "
          f"vLess={db_stats['tiers']['vLess']}  "
          f"Ambiguous={db_stats['tiers']['Ambiguous']}  "
          f"vNo={db_stats['tiers']['vNo']}")

    # ── default threshold = 30 (v0.2 baseline)
    m = compute_metrics(all_drugs, threshold=30)
    opt_thr, opt_j = find_optimal_threshold(all_drugs)
    m_opt = compute_metrics(all_drugs, threshold=opt_thr)

    print(f"\n  ── Fixed threshold (30) ────────────────────────────────────")
    print(f"  Sensitivity : {m['sens']:.3f}   Specificity : {m['spec']:.3f}")
    print(f"  Accuracy    : {m['acc']:.3f}   PPV         : {m['ppv']:.3f}")
    print(f"  F1 Score    : {m['f1']:.3f}   AUC-ROC     : {m['auc']:.3f}")
    print(f"\n  Confusion Matrix (thr=30):")
    print(f"               Pred DILI+  Pred DILI-")
    print(f"  Actual DILI+  {m['TP']:8d}  {m['FN']:8d}")
    print(f"  Actual DILI-  {m['FP']:8d}  {m['TN']:8d}")

    print(f"\n  ── Optimal threshold ({opt_thr}) — Youden J={opt_j:.3f} ──────────────")
    print(f"  Sensitivity : {m_opt['sens']:.3f}   Specificity : {m_opt['spec']:.3f}")
    print(f"  Accuracy    : {m_opt['acc']:.3f}   PPV         : {m_opt['ppv']:.3f}")
    print(f"  F1 Score    : {m_opt['f1']:.3f}   AUC-ROC     : {m_opt['auc']:.3f}")

    # ── Tier breakdown
    print(f"\n  ── Per-tier AUC (binary vMost/vLess=1, vNo=0) ──────────────")
    tiers = tier_breakdown(get_all_drugs())
    for pair_name, t1, t2 in [
        ("vMost vs vNo",      "vMost", "vNo"),
        ("vLess vs vNo",      "vLess", "vNo"),
        ("vMost+vLess vs vNo","combined", "vNo"),
    ]:
        if t1 == "combined":
            subset = tiers["vMost"] + tiers["vLess"] + tiers["vNo"]
        else:
            subset = tiers[t1] + tiers["vNo"]
        if len(subset) < 5:
            continue
        sub_m = compute_metrics(subset, threshold=30)
        print(f"  {pair_name:<26}: AUC={sub_m['auc']:.3f}  "
              f"Sens={sub_m['sens']:.3f}  Spec={sub_m['spec']:.3f}  "
              f"n={len(subset)}")

    # ── Monte Carlo (Acetaminophen)
    np.random.seed(42)
    apap = next((d for d in all_drugs if d.name == "Acetaminophen"), all_drugs[0])
    mc = []
    for _ in range(1000):
        cd = DrugVV(
            apap.name, apap.rank, apap.binary, apap.mw,
            max(0, apap.cyp * np.random.uniform(0.8, 1.2)), apap.reactive,
            max(1, apap.dose * np.random.uniform(0.8, 1.2)),
            min(0.99, max(0, apap.binding * np.random.uniform(0.8, 1.2))),
            max(0.1, apap.cl * np.random.uniform(0.8, 1.2)), apap.vd,
            min(1, max(0.1, apap.bio * np.random.uniform(0.8, 1.2))))
        mc.append(calc_dili_score(cd))
    mc = np.array(mc)
    print(f"\n  Monte Carlo APAP (N=1000): "
          f"Mean={np.mean(mc):.1f}  "
          f"95%CI=[{np.percentile(mc,2.5):.1f}, {np.percentile(mc,97.5):.1f}]")

    # ── Plots ─────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    fig.suptitle(
        'DTIH V&V Phase 2: Expanded DILIrank Validation (212 drugs)\n'
        'Erlanger Protocol v0.3 — ASME V&V 40 Framework',
        fontsize=13, fontweight='bold')

    COLORS = {
        'vMost': '#e74c3c', 'vLess': '#e67e22',
        'Ambiguous': '#95a5a6', 'vNo': '#27ae60'
    }

    # (a) Tier scatter
    ax = axes[0, 0]
    for i, (rank, col) in enumerate(COLORS.items()):
        sv = [s[3] for s in m['scores'] if s[1] == rank]
        ax.scatter([i] * len(sv), sv, c=col, s=40, alpha=0.6,
                   edgecolors='white', lw=0.4, label=rank)
        if sv:
            ax.plot([i - 0.25, i + 0.25], [np.mean(sv)] * 2, c=col, lw=3)
    ax.axhline(y=30, color='red', ls='--', alpha=0.5, lw=1.2, label='Thr=30')
    ax.set_xticks(range(4))
    ax.set_xticklabels(['vMost\n(n=%d)' % len(tiers['vMost']),
                         'vLess\n(n=%d)' % len(tiers['vLess']),
                         'Ambiguous\n(n=%d)' % len(tiers['Ambiguous']),
                         'vNo\n(n=%d)' % len(tiers['vNo'])])
    ax.set_ylabel('DILI Risk Score')
    ax.set_title('(a) Score by DILIrank Tier (212 drugs)')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)

    # (b) ROC curve — full dataset
    ax = axes[0, 1]
    ax.plot(m['fpr'], m['tpr'], 'b-', lw=2, label=f'All (AUC={m["auc"]:.3f})')
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='Random')
    ax.fill_between(m['fpr'], m['tpr'], alpha=0.08, color='blue')
    # Overlay vMost vs vNo
    sub = tiers['vMost'] + tiers['vNo']
    sm = compute_metrics(sub, threshold=30)
    ax.plot(sm['fpr'], sm['tpr'], 'r--', lw=1.5,
            label=f'vMost vs vNo (AUC={sm["auc"]:.3f})')
    ax.set_xlabel('FPR (1−Specificity)')
    ax.set_ylabel('TPR (Sensitivity)')
    ax.set_title('(b) ROC Curve')
    ax.legend(fontsize=9, loc='lower right'); ax.grid(True, alpha=0.2)

    # (c) Score distribution histogram
    ax = axes[0, 2]
    for rank, col in COLORS.items():
        sv = [s[3] for s in m['scores'] if s[1] == rank]
        if sv:
            ax.hist(sv, bins=20, color=col, alpha=0.55,
                    edgecolor='white', lw=0.3, label=rank)
    ax.axvline(x=30, color='red', ls='--', lw=1.5, label='Thr=30')
    ax.axvline(x=opt_thr, color='purple', ls=':', lw=1.5,
               label=f'Opt thr={opt_thr}')
    ax.set_xlabel('DILI Risk Score')
    ax.set_ylabel('Count')
    ax.set_title('(c) Score Distribution by Tier')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)

    # (d) Confusion matrix heatmap
    ax = axes[1, 0]
    cm = np.array([[m['TP'], m['FN']], [m['FP'], m['TN']]])
    im = ax.imshow(cm, interpolation='nearest', cmap='Blues')
    plt.colorbar(im, ax=ax)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(cm[i, j]), ha='center', va='center',
                    fontsize=16, fontweight='bold',
                    color='white' if cm[i, j] > cm.max() / 2 else 'black')
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(['Pred DILI+', 'Pred DILI−'])
    ax.set_yticklabels(['Actual DILI+', 'Actual DILI−'])
    ax.set_title(f'(d) Confusion Matrix (thr={30})\n'
                 f'Sens={m["sens"]:.3f}  Spec={m["spec"]:.3f}  '
                 f'F1={m["f1"]:.3f}')

    # (e) Sensitivity tornado (Acetaminophen)
    ax = axes[1, 1]
    params = [('cyp', 'CYP fraction'), ('dose', 'Daily dose'),
              ('binding', 'Protein binding'), ('cl', 'Clearance'),
              ('bio', 'Bioavailability')]
    base = calc_dili_score(apap)
    lbls, dls, dhs = [], [], []
    for attr, lbl in params:
        bv = getattr(apap, attr)
        dl2 = DrugVV(apap.name, apap.rank, apap.binary, apap.mw,
                     apap.cyp, apap.reactive, apap.dose, apap.binding,
                     apap.cl, apap.vd, apap.bio)
        setattr(dl2, attr, bv * 0.8)
        dh2 = DrugVV(apap.name, apap.rank, apap.binary, apap.mw,
                     apap.cyp, apap.reactive, apap.dose, apap.binding,
                     apap.cl, apap.vd, apap.bio)
        setattr(dh2, attr, bv * 1.2)
        lbls.append(lbl)
        dls.append(calc_dili_score(dl2) - base)
        dhs.append(calc_dili_score(dh2) - base)
    y = range(len(lbls))
    ax.barh(y, dhs, height=0.4, color='#e74c3c', alpha=0.8, label='+20%')
    ax.barh(y, dls, height=0.4, color='#3498db', alpha=0.8, label='−20%')
    ax.set_yticks(y); ax.set_yticklabels(lbls)
    ax.set_xlabel('DILI Score Change')
    ax.set_title('(e) Sensitivity Tornado (Acetaminophen)')
    ax.axvline(x=0, color='black', lw=0.5)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.2, axis='x')

    # (f) Monte Carlo
    ax = axes[1, 2]
    ax.hist(mc, bins=30, color='#9b59b6', alpha=0.7, edgecolor='white')
    ax.axvline(np.mean(mc), color='red', lw=2,
               label=f'Mean={np.mean(mc):.1f}')
    ax.axvline(np.percentile(mc, 2.5), color='orange', lw=1.5, ls='--',
               label='95% CI')
    ax.axvline(np.percentile(mc, 97.5), color='orange', lw=1.5, ls='--')
    ax.set_xlabel('DILI Risk Score')
    ax.set_ylabel('Frequency')
    ax.set_title('(f) Monte Carlo (N=1000, Acetaminophen)')
    ax.legend(fontsize=9); ax.grid(True, alpha=0.2)

    plt.tight_layout()
    path = os.path.join(output_dir, "VV_Phase2_Results.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved: {path}")

    # ── v0.2 → v0.3 comparison summary ────────────────────────────────────────
    print(f"\n  ── v0.2 → v0.3 Comparison ──────────────────────────────────")
    print(f"  {'Metric':<20} {'v0.2 (38 drugs)':>18} {'v0.3 (212 drugs)':>18}")
    print(f"  {'-'*56}")
    # v0.2 known values
    v02 = dict(AUC=0.808, Sens=0.826, Spec=0.667, F1=0.826)
    v03 = dict(AUC=m['auc'], Sens=m['sens'], Spec=m['spec'], F1=m['f1'])
    for k in ['AUC', 'Sens', 'Spec', 'F1']:
        delta = v03[k] - v02[k]
        arrow = '▲' if delta > 0.005 else ('▼' if delta < -0.005 else '≈')
        print(f"  {k:<20} {v02[k]:>18.3f} {v03[k]:>18.3f}  {arrow}{abs(delta):.3f}")

    return m['auc']


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DTIH Phase 2 V&V")
    parser.add_argument("--tier", default="all",
                        choices=["all", "vMost", "vLess", "Ambiguous", "vNo"],
                        help="Filter by DILIrank tier")
    parser.add_argument("--output", default=".",
                        help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    auc = run_vv_phase2(output_dir=args.output, tier_filter=args.tier)
    print(f"\n  Final AUC-ROC: {auc:.3f}\n")


# ══════════════════════════════════════════════════════════════
# Phase 2 Extension: calc_dili_score_v2
# Adds non-CYP mechanism terms to cover:
#   T1: Immune-mediated (biologics, allopurinol, sulfonamides)
#   T2: DNA direct damage (platinum, antimetabolite, alkylating)
#   T3: Mitochondrial toxicity (NRTIs, valproate-like)
#
# Formula (Def 2.1 extension):
#   S_v2 = S_cyp + S_dose + S_binding + S_ml   (original 4 terms, 0–75)
#         + T1_immune * 15                       (0–15)
#         + T2_genotoxic * 15                    (0–15)
#         + T3_mito * 10                         (0–10)
#   Total max = 115, clipped to 99
#
# Mechanism flags encoded via protein_binding proxy heuristic:
#   biologic_mab   : mw > 5000  (large molecules)
#   genotoxic_chemo: dose < 200 AND binding < 0.50 AND cyp < 0.30
#   mito_toxicity  : cyp == 0.0 AND reactive == False AND bio > 0.5
# ══════════════════════════════════════════════════════════════

def calc_dili_score_v2(d: DrugVV) -> float:
    """Extended DILI score with non-CYP mechanism terms."""
    # ── Original 4 terms (unchanged for backward compatibility) ──
    rs = min(1.0, d.cyp * (2.0 if d.reactive else 0.5)) * 40
    ds = min(1.0, np.log10(d.dose / 100) / 1.5) * 25 if d.dose > 100 else 0
    bs = d.binding * (1 - min(1, d.cl / 30)) * 15
    ml = min(1.0, d.dose * d.bio * d.cyp / 500) * 20
    base = rs + ds + bs + ml

    # ── T1: Immune-mediated (biologic mAb or immune sensitiser) ──
    # Proxy: very large MW (>5000 = likely biologic) OR
    #        known immune-trigger drugs (binding > 0.90, low cyp, low dose)
    is_biologic = d.mw > 5000
    is_immune_sensitiser = (d.cyp < 0.30 and d.binding > 0.80
                            and d.dose < 600 and d.reactive is False
                            and not is_biologic)
    t1 = 12.0 if is_biologic else (8.0 if is_immune_sensitiser else 0.0)

    # ── T2: Genotoxic/direct DNA damage ──
    # Proxy: very low CYP (<0.15), low-to-moderate binding,
    #        and drug is IV (bio==0) or antimetabolite pattern
    is_genotoxic = (d.cyp <= 0.15 and d.binding <= 0.50
                    and (d.bio == 0.00 or d.dose < 300))
    t2 = 10.0 if is_genotoxic else 0.0

    # ── T3: Mitochondrial (NRTI-like) ──
    # Proxy: cyp==0, not reactive, oral bioavailability >0.3, low dose
    is_mito = (d.cyp == 0.0 and d.reactive is False
               and d.bio > 0.30 and d.dose < 500)
    t3 = 8.0 if is_mito else 0.0

    return min(99, max(0, base + t1 + t2 + t3))


def compute_metrics_v2(drugs, threshold=30):
    """Run V&V using calc_dili_score_v2."""
    scores  = [(d.name, d.rank, d.binary, calc_dili_score_v2(d)) for d in drugs]
    actual  = np.array([s[2] for s in scores])
    sc      = np.array([s[3] for s in scores])
    pred    = (sc >= threshold).astype(int)

    TP = int(np.sum((pred == 1) & (actual == 1)))
    TN = int(np.sum((pred == 0) & (actual == 0)))
    FP = int(np.sum((pred == 1) & (actual == 0)))
    FN = int(np.sum((pred == 0) & (actual == 1)))

    sens = TP / max(TP + FN, 1)
    spec = TN / max(TN + FP, 1)
    acc  = (TP + TN) / len(scores)
    ppv  = TP / max(TP + FP, 1)
    f1   = 2 * ppv * sens / max(ppv + sens, 0.001)

    thresholds = np.arange(0, 100, 1)
    tpr_l, fpr_l = [], []
    for thr in thresholds:
        p  = (sc >= thr).astype(int)
        tp = np.sum((p == 1) & (actual == 1))
        fp = np.sum((p == 1) & (actual == 0))
        fn = np.sum((p == 0) & (actual == 1))
        tn = np.sum((p == 0) & (actual == 0))
        tpr_l.append(tp / max(tp + fn, 1))
        fpr_l.append(fp / max(fp + tn, 1))
    fpr_a, tpr_a = np.array(fpr_l), np.array(tpr_l)
    si  = np.argsort(fpr_a)
    auc = float(np.abs(np.trapezoid(tpr_a[si], fpr_a[si])))

    return dict(
        scores=scores, actual=actual, sc=sc, pred=pred,
        TP=TP, TN=TN, FP=FP, FN=FN,
        sens=sens, spec=spec, acc=acc, ppv=ppv, f1=f1, auc=auc,
        fpr=fpr_l, tpr=tpr_l
    )
