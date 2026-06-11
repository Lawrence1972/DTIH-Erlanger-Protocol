"""
================================================================
  DTIH Prototype v0.2 — Erlanger Protocol
  Digital Twin In-silico Human
  
  This is the COMPLETE working prototype.
  Run this single file to execute all simulations and validation.
================================================================

Output files (saved to /results folder):
  - dtih_acetaminophen_results.png      (therapeutic dose, 6 patients)
  - dtih_acetaminophen_od_results.png   (overdose 10g, 6 patients)  
  - dtih_cisplatin_results.png          (nephro+neurotoxicity)
  - dtih_doxorubicin_results.png        (cardiotoxicity)
  - VV_Results.png                      (FDA DILIrank validation, AUC-ROC)
  - *_summary.json                      (numerical data)

Requirements: Python 3.10+, numpy, scipy, matplotlib
Install:      pip install numpy scipy matplotlib
"""

import os
import sys

def main():
    print()
    print("  ╔═══════════════════════════════════════════════════╗")
    print("  ║  DTIH Prototype v0.2 — Erlanger Protocol         ║")
    print("  ║  Digital Twin In-silico Human                     ║")
    print("  ║                                                   ║")
    print("  ║  L1 Genomic → L2 Metabolic → L3 PBPK (3 organs) ║")
    print("  ║  → L4a Neural (Erlanger) → L4b Immune (Cytokine) ║")
    print("  ║  → L5 Integration (Feedback Loops)               ║")
    print("  ╚═══════════════════════════════════════════════════╝")
    print()
    
    # Output directory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(base_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    
    # ─── Part 1: Drug Simulations (L1→L2→L3→L4→L5) ───
    print("  [1/2] Running drug simulations...")
    print("  " + "─" * 50)
    
    from simulator import run_scenario
    
    run_scenario("acetaminophen", output_dir=results_dir)
    run_scenario("acetaminophen_od", output_dir=results_dir)
    run_scenario("cisplatin", output_dir=results_dir)
    run_scenario("doxorubicin", output_dir=results_dir)
    
    # ─── Part 2: V&V Validation (ASME V&V 40) ───
    print("\n  [2/2] Running V&V validation against FDA DILIrank...")
    print("  " + "─" * 50)
    
    from validation import run_vv
    auc = run_vv(output_dir=results_dir)
    
    # ─── Summary ───
    print()
    print("  ╔═══════════════════════════════════════════════════╗")
    print("  ║  ALL COMPLETE                                     ║")
    print(f"  ║  V&V AUC-ROC: {auc:.3f}                              ║")
    print("  ║                                                   ║")
    print("  ║  Results saved to: /results/                      ║")
    print("  ╚═══════════════════════════════════════════════════╝")
    print()
    print(f"  Output folder: {results_dir}")
    print()

if __name__ == "__main__":
    main()
