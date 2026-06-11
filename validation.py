"""
DTIH V&V Module: FDA DILIrank Retrospective Validation
ASME V&V 40 Framework — Erlanger Protocol
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass
import os

@dataclass
class DrugVV:
    name: str; rank: str; binary: int; mw: float
    cyp: float; reactive: bool; dose: float; binding: float
    cl: float; vd: float; bio: float

DILIRANK_38 = [
    DrugVV("Acetaminophen","vMost",1,151.2,0.10,True,4000,0.20,5.0,0.9,0.85),
    DrugVV("Isoniazid","vMost",1,137.1,0.70,True,300,0.10,12.0,0.6,0.90),
    DrugVV("Troglitazone","vMost",1,441.5,0.85,True,600,0.99,2.5,1.2,0.50),
    DrugVV("Ketoconazole","vMost",1,531.4,0.90,True,400,0.99,8.0,2.4,0.50),
    DrugVV("Diclofenac","vMost",1,296.1,0.75,True,150,0.99,15.0,0.2,0.55),
    DrugVV("Valproic Acid","vMost",1,144.2,0.30,True,2000,0.90,1.0,0.2,0.90),
    DrugVV("Amiodarone","vMost",1,645.3,0.50,False,400,0.96,1.5,70.0,0.50),
    DrugVV("Nefazodone","vMost",1,470.0,0.95,True,600,0.99,25.0,0.5,0.20),
    DrugVV("Tolcapone","vMost",1,273.2,0.80,True,300,0.99,30.0,0.3,0.65),
    DrugVV("Flutamide","vMost",1,276.2,0.85,True,750,0.94,20.0,0.3,0.80),
    DrugVV("Chlorzoxazone","vMost",1,169.6,0.90,True,750,0.15,35.0,0.2,0.90),
    DrugVV("Dantrolene","vMost",1,314.3,0.60,True,400,0.90,5.0,1.0,0.70),
    DrugVV("Ciprofloxacin","vLess",1,331.3,0.30,False,1000,0.30,8.0,2.0,0.70),
    DrugVV("Atorvastatin","vLess",1,558.6,0.80,False,80,0.98,20.0,5.5,0.14),
    DrugVV("Azithromycin","vLess",1,748.9,0.20,False,500,0.50,3.0,31.0,0.37),
    DrugVV("Losartan","vLess",1,422.9,0.85,False,100,0.98,10.0,0.5,0.33),
    DrugVV("Carbamazepine","vLess",1,236.3,0.90,True,1200,0.76,5.0,1.2,0.75),
    DrugVV("Fluconazole","vLess",1,306.3,0.10,False,400,0.12,2.0,0.7,0.90),
    DrugVV("Ibuprofen","vLess",1,206.3,0.55,False,2400,0.99,8.0,0.1,0.95),
    DrugVV("Methotrexate","vLess",1,454.4,0.20,False,25,0.50,3.0,0.5,0.60),
    DrugVV("Omeprazole","vLess",1,345.4,0.85,False,40,0.95,10.0,0.3,0.50),
    DrugVV("Phenytoin","vLess",1,252.3,0.90,True,400,0.90,5.0,0.7,0.90),
    DrugVV("Aspirin","vNo",0,180.2,0.10,False,4000,0.80,10.0,0.1,0.70),
    DrugVV("Metformin","vNo",0,129.2,0.00,False,2550,0.00,0.5,1.1,0.55),
    DrugVV("Lisinopril","vNo",0,405.5,0.00,False,40,0.00,1.0,1.0,0.25),
    DrugVV("Amlodipine","vNo",0,408.9,0.10,False,10,0.97,1.5,21.0,0.65),
    DrugVV("Levothyroxine","vNo",0,776.9,0.00,False,0.3,0.99,0.2,0.2,0.50),
    DrugVV("Metoprolol","vNo",0,267.4,0.80,False,400,0.12,15.0,3.2,0.50),
    DrugVV("HCTZ","vNo",0,297.7,0.05,False,50,0.40,1.0,3.0,0.70),
    DrugVV("Furosemide","vNo",0,330.7,0.10,False,80,0.91,2.0,0.1,0.60),
    DrugVV("Penicillin V","vNo",0,350.4,0.05,False,2000,0.80,3.0,0.4,0.60),
    DrugVV("Cetirizine","vNo",0,388.9,0.05,False,10,0.93,0.5,0.5,0.70),
    DrugVV("Ranitidine","vNo",0,314.4,0.10,False,300,0.15,8.0,1.4,0.50),
    DrugVV("Gabapentin","vNo",0,171.2,0.00,False,3600,0.00,0.1,0.8,0.60),
    DrugVV("Amoxicillin","Ambiguous",0,365.4,0.10,False,3000,0.20,2.0,0.3,0.75),
    DrugVV("Cephalexin","Ambiguous",0,347.4,0.05,False,4000,0.15,1.5,0.3,0.90),
    DrugVV("Prednisone","Ambiguous",1,358.4,0.70,False,60,0.70,5.0,1.0,0.80),
    DrugVV("Doxycycline","Ambiguous",1,444.4,0.15,False,200,0.93,1.0,0.7,0.95),
]

def calc_dili_score(d):
    rs = min(1.0, d.cyp * (2.0 if d.reactive else 0.5)) * 40
    ds = min(1.0, np.log10(d.dose/100)/1.5)*25 if d.dose > 100 else 0
    bs = d.binding * (1 - min(1, d.cl/30)) * 15
    ml = min(1.0, d.dose * d.bio * d.cyp / 500) * 20
    return min(99, max(0, rs + ds + bs + ml))

def run_vv(output_dir="."):
    print("\n" + "="*60)
    print("  V&V: FDA DILIrank Retrospective Validation")
    print("="*60)

    scores = [(d.name, d.rank, d.binary, calc_dili_score(d)) for d in DILIRANK_38]
    actual = np.array([s[2] for s in scores])
    sc = np.array([s[3] for s in scores])
    predicted = (sc >= 30).astype(int)

    TP = np.sum((predicted==1)&(actual==1))
    TN = np.sum((predicted==0)&(actual==0))
    FP = np.sum((predicted==1)&(actual==0))
    FN = np.sum((predicted==0)&(actual==1))

    sens = TP/max(TP+FN,1)
    spec = TN/max(TN+FP,1)
    acc = (TP+TN)/len(scores)
    ppv = TP/max(TP+FP,1)
    f1 = 2*ppv*sens/max(ppv+sens,0.001)

    # ROC
    thresholds = np.arange(0, 100, 2)
    tpr_l, fpr_l = [], []
    for thr in thresholds:
        p = (sc >= thr).astype(int)
        tp = np.sum((p==1)&(actual==1)); fp = np.sum((p==1)&(actual==0))
        fn = np.sum((p==0)&(actual==1)); tn = np.sum((p==0)&(actual==0))
        tpr_l.append(tp/max(tp+fn,1)); fpr_l.append(fp/max(fp+tn,1))
    fpr_a, tpr_a = np.array(fpr_l), np.array(tpr_l)
    si = np.argsort(fpr_a)
    auc = np.abs(np.trapezoid(tpr_a[si], fpr_a[si]))

    print(f"\n  Sensitivity: {sens:.3f}  |  Specificity: {spec:.3f}")
    print(f"  Accuracy:    {acc:.3f}  |  F1 Score:    {f1:.3f}")
    print(f"  AUC-ROC:     {auc:.3f}")
    print(f"\n  Confusion Matrix:")
    print(f"               Pred DILI+  Pred DILI-")
    print(f"  Actual DILI+    {TP:5d}       {FN:5d}")
    print(f"  Actual DILI-    {FP:5d}       {TN:5d}")

    # Monte Carlo
    np.random.seed(42)
    d0 = DILIRANK_38[0]
    mc = []
    for _ in range(1000):
        cd = DrugVV(d0.name,d0.rank,d0.binary,d0.mw,
            max(0,d0.cyp*np.random.uniform(0.8,1.2)),d0.reactive,
            max(1,d0.dose*np.random.uniform(0.8,1.2)),
            min(0.99,max(0,d0.binding*np.random.uniform(0.8,1.2))),
            max(0.1,d0.cl*np.random.uniform(0.8,1.2)),d0.vd,
            min(1,max(0.1,d0.bio*np.random.uniform(0.8,1.2))))
        mc.append(calc_dili_score(cd))
    mc = np.array(mc)
    print(f"\n  Monte Carlo (N=1000): Mean={np.mean(mc):.1f}, 95%CI=[{np.percentile(mc,2.5):.1f}, {np.percentile(mc,97.5):.1f}]")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle('DTIH V&V Report: Retrospective Validation against FDA DILIrank\n'
                 'ASME V&V 40 Framework - Erlanger Protocol', fontsize=14, fontweight='bold')

    ax = axes[0,0]
    for i, (rank, col, lbl) in enumerate([('vMost','#ff4444','vMost\n(High Risk)'),
        ('vLess','#ff9900','vLess\n(Moderate)'),('Ambiguous','#888888','Ambiguous\n(Uncertain)'),
        ('vNo','#00cc66','vNo\n(Safe)')]):
        sv = [s[3] for s in scores if s[1]==rank]
        ax.scatter([i]*len(sv), sv, c=col, s=60, alpha=0.7, edgecolors='white', lw=0.5)
        if sv: ax.plot([i-0.2,i+0.2],[np.mean(sv)]*2, c=col, lw=3)
    ax.axhline(y=30, color='red', ls='--', alpha=0.5, label='Threshold (30)')
    ax.set_xticks(range(4))
    ax.set_xticklabels(['vMost\n(High Risk)','vLess\n(Moderate)','Ambiguous\n(Uncertain)','vNo\n(Safe)'])
    ax.set_ylabel('DILI Risk Score'); ax.set_title('(a) DILIrank vs DTIH Prediction')
    ax.legend(fontsize=9); ax.grid(True, alpha=0.2)

    ax = axes[0,1]
    ax.plot(fpr_l, tpr_l, 'b-', lw=2, label=f'DTIH (AUC={auc:.3f})')
    ax.plot([0,1],[0,1],'k--',alpha=0.3,label='Random (0.5)')
    ax.fill_between(fpr_l, tpr_l, alpha=0.1, color='blue')
    ax.set_xlabel('FPR (1-Specificity)'); ax.set_ylabel('TPR (Sensitivity)')
    ax.set_title('(b) ROC Curve'); ax.legend(fontsize=10, loc='lower right')
    ax.grid(True, alpha=0.2)

    ax = axes[1,0]
    params = [('cyp','CYP Pathway'),('dose','Daily Dose'),('binding','Protein Binding'),
              ('cl','Clearance'),('bio','Bioavailability')]
    base = calc_dili_score(d0)
    lbls, dls, dhs = [], [], []
    for attr, lbl in params:
        bv = getattr(d0, attr)
        dl2 = DrugVV(d0.name,d0.rank,d0.binary,d0.mw,d0.cyp,d0.reactive,d0.dose,d0.binding,d0.cl,d0.vd,d0.bio)
        setattr(dl2, attr, bv*0.8); sl = calc_dili_score(dl2)
        dh2 = DrugVV(d0.name,d0.rank,d0.binary,d0.mw,d0.cyp,d0.reactive,d0.dose,d0.binding,d0.cl,d0.vd,d0.bio)
        setattr(dh2, attr, bv*1.2); sh = calc_dili_score(dh2)
        lbls.append(lbl); dls.append(sl-base); dhs.append(sh-base)
    y = range(len(lbls))
    ax.barh(y, dhs, height=0.4, color='#ff6b6b', alpha=0.8, label='+20%')
    ax.barh(y, dls, height=0.4, color='#4dabf7', alpha=0.8, label='-20%')
    ax.set_yticks(y); ax.set_yticklabels(lbls)
    ax.set_xlabel('DILI Score Change'); ax.set_title('(c) Sensitivity (Tornado)')
    ax.axvline(x=0, color='black', lw=0.5); ax.legend(); ax.grid(True, alpha=0.2, axis='x')

    ax = axes[1,1]
    ax.hist(mc, bins=30, color='#748ffc', alpha=0.7, edgecolor='white')
    ax.axvline(x=np.mean(mc), color='red', lw=2, label=f'Mean={np.mean(mc):.1f}')
    ax.axvline(x=np.percentile(mc,2.5), color='orange', lw=1.5, ls='--', label='95% CI')
    ax.axvline(x=np.percentile(mc,97.5), color='orange', lw=1.5, ls='--')
    ax.set_xlabel('DILI Risk Score'); ax.set_ylabel('Frequency')
    ax.set_title('(d) Monte Carlo (N=1000, Acetaminophen)'); ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    path = os.path.join(output_dir, "VV_Results.png")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved: {path}")
    return auc
