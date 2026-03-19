"""
MDR-TB AI Pipeline v3 — Master Pipeline Runner (Phase 14)
==========================================================
Executes all modules in sequence:
  1. Data loading (GEO, WHO, DrugBank, PubChem)
  2. Omics engine (DE, pathway, gene ranking)
  3. Epidemiology engine (trends, MDR patterns, districts)
  4. Target engine (multi-factor scoring)
  5. Compound engine (descriptors, clustering, Lipinski)
  6. ML pipeline (training + evaluation)
  7. Resistance module (gene mapping, scoring)
  8. Ranking engine (multi-objective scoring)
  9. Manuscript generator
  10. ICMR grant generator
  11. IEC protocol generator

SAFETY: Computational only. No wet-lab or synthesis steps.
"""

import sys
import time
import traceback
from pathlib import Path

# Ensure project root is on path
PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT))

from utils.helpers import load_config, log_phase, safe_save_json, ensure_dir


def run_pipeline():
    """Execute the complete MDR-TB AI pipeline."""
    start_time = time.time()
    config = load_config(str(PROJECT_ROOT / "config.yaml"))
    results = {}
    errors = []

    print("=" * 70)
    print("  MDR-TB AI Drug Discovery Pipeline v3")
    print("  SAFETY: Computational Only — No Wet-Lab Protocols")
    print("=" * 70)

    # ── PHASE 1: DATA LOADING ──────────────────────────────────
    log_phase("DATA LOADING")
    
    # GEO
    try:
        from data_connectors.geo_loader import GEOLoader
        geo = GEOLoader(config)
        expression_df = geo.load()
        results["expression"] = expression_df
        print(f"  → Expression matrix: {expression_df.shape}")
    except Exception as e:
        errors.append(f"GEO: {e}")
        traceback.print_exc()
        expression_df = None

    # WHO
    try:
        from data_connectors.who_parser import WHOParser
        who = WHOParser(config)
        who_df = who.load()
        results["who"] = who_df
        print(f"  → WHO data: {who_df.shape}")
    except Exception as e:
        errors.append(f"WHO: {e}")
        traceback.print_exc()
        who_df = None

    # CRyPTIC Genomic Atlas
    try:
        from data_connectors.cryptic_loader import CrypticLoader
        cryptic = CrypticLoader(config)
        cryptic_res = cryptic.run()
        results["cryptic"] = cryptic_res
    except Exception as e:
        errors.append(f"CRyPTIC: {e}")
        traceback.print_exc()

    # DrugBank
    try:
        from data_connectors.drugbank_loader import DrugBankLoader
        db = DrugBankLoader(config)
        drugbank_df = db.load()
        results["drugbank"] = drugbank_df
        print(f"  → DrugBank entries: {drugbank_df.shape}")
    except Exception as e:
        errors.append(f"DrugBank: {e}")
        traceback.print_exc()
        drugbank_df = None

    # PubChem
    try:
        from data_connectors.pubchem_loader import PubChemLoader
        pc = PubChemLoader(config)
        compounds_df = pc.load()
        results["compounds_raw"] = compounds_df
        print(f"  → PubChem compounds: {compounds_df.shape}")
    except Exception as e:
        errors.append(f"PubChem: {e}")
        traceback.print_exc()
        compounds_df = None

    # ── PHASE 2: OMICS ENGINE ─────────────────────────────────
    log_phase("OMICS ENGINE")
    omics_results = {}
    try:
        from engines.omics_engine import OmicsEngine
        omics = OmicsEngine(config)
        if expression_df is not None:
            omics_results = omics.run(expression_df)
            results["omics"] = omics_results
        else:
            print("  ⚠ Skipping omics — no expression data")
    except Exception as e:
        errors.append(f"Omics: {e}")
        traceback.print_exc()

    # ── PHASE 3: EPIDEMIOLOGY ENGINE ──────────────────────────
    log_phase("EPIDEMIOLOGY ENGINE")
    epi_results = {}
    try:
        from engines.epi_engine import EpidemiologyEngine
        epi = EpidemiologyEngine(config)
        if who_df is not None:
            epi_results = epi.run(who_df)
            results["epi"] = epi_results
        else:
            print("  ⚠ Skipping epi — no WHO data")
    except Exception as e:
        errors.append(f"Epi: {e}")
        traceback.print_exc()

    # ── PHASE 4: TARGET ENGINE ────────────────────────────────
    log_phase("TARGET ENGINE")
    target_results = {}
    try:
        from engines.target_engine import TargetEngine
        te = TargetEngine(config)
        de_results = omics_results.get("de_results", None)
        resistance_genes = config.get("resistance", {}).get("mdr_genes", [])
        target_results = te.run(de_results, resistance_genes)
        results["targets"] = target_results
    except Exception as e:
        errors.append(f"Target: {e}")
        traceback.print_exc()

    # ── PHASE 5: DE NOVO DRUG DESIGN ──────────────────────────
    log_phase("DE NOVO GENERATOR")
    try:
        from models.de_novo_generator import DeNovoGenerator
        import pandas as pd
        dng = DeNovoGenerator(config)
        novel_df = dng.run()
        if compounds_df is not None:
             compounds_df = pd.concat([compounds_df, novel_df], ignore_index=True)
        else:
             compounds_df = novel_df
    except Exception as e:
        errors.append(f"DeNovo: {e}")
        traceback.print_exc()

    # ── PHASE 5a: COMPOUND ENGINE ──────────────────────────────
    log_phase("COMPOUND ENGINE")
    compound_results = {}
    try:
        from engines.compound_engine import CompoundEngine
        ce = CompoundEngine(config)
        if compounds_df is not None:
            compound_results = ce.run(compounds_df)
            results["compounds"] = compound_results
            compounds_processed = compound_results.get("clustered", compounds_df)
        else:
            print("  ⚠ Skipping compound engine — no compound data")
            compounds_processed = None
    except Exception as e:
        errors.append(f"Compound: {e}")
        traceback.print_exc()
        compounds_processed = None

    # ── PHASE 5b: ALPHAFOLD / ESMFOLD STRUCTURAL PREDICTION ───
    log_phase("ALPHAFOLD / ESMFOLD PREDICTOR")
    try:
        from models.structural_predictor import ESMFoldPredictor
        from engines.alphafold_engine import AlphaFoldEngine
        
        af = AlphaFoldEngine(config)
        esm = ESMFoldPredictor(config)
        
        top_target = "DprE1"
        if target_results and "scored_targets" in target_results:
            top_target = target_results["scored_targets"].iloc[0]["Target"][:5]
            
        print(f"  [Orchestrator] Requesting structural scaffold for priority target: {top_target}")
        
        # Primary Hook: DeepMind AlphaFold EBI Database
        predicted_pdb_path = af.run(top_target)
        
        # Fallback Hook: Meta ESMFold Sequence Predictor
        if not predicted_pdb_path:
             print("  [Orchestrator] ⚠ AlphaFold structural query missed. Cascading to ESMFold inference...")
             predicted_pdb_path = esm.run(top_target)
             
    except Exception as e:
        errors.append(f"StructuralPredictors: {e}")
        traceback.print_exc()

    # ── PHASE 5c: DOCKING ENGINE ────────────────────────────
    log_phase("DOCKING ENGINE")
    docking_results = {}
    try:
        from engines.docking_engine import DockingEngine
        de = DockingEngine(config)
        top_target = "DprE1"
        if target_results and "scored_targets" in target_results:
            top_target = target_results["scored_targets"].iloc[0]["Target"][:5]
            
        if compounds_processed is not None:
            docking_results = de.run(compounds_processed, target_name=top_target)
            results["docking"] = docking_results
        else:
            print("  ⚠ Skipping docking — no compound data")
    except Exception as e:
        errors.append(f"Docking: {e}")
        traceback.print_exc()

    # ── PHASE 6: ML PIPELINE ─────────────────────────────────
    log_phase("ML PIPELINE")
    ml_results = {}
    try:
        from models.ml_pipeline import MLPipeline
        ml = MLPipeline(config)
        if compounds_processed is not None:
            ml_results = ml.run(compounds_processed)
            results["ml"] = ml_results
        else:
            print("  ⚠ Skipping ML — no compound data")
    except Exception as e:
        errors.append(f"ML: {e}")
        traceback.print_exc()

    # ── PHASE 7: RESISTANCE MODULE ────────────────────────────
    log_phase("RESISTANCE MODULE")
    resistance_results = {}
    try:
        from models.resistance_module import ResistanceModule
        rm = ResistanceModule(config)
        de_res = omics_results.get("de_results", None)
        resistance_results = rm.run(de_res)
        results["resistance"] = resistance_results
    except Exception as e:
        errors.append(f"Resistance: {e}")
        traceback.print_exc()

    # ── PHASE 8: RANKING ENGINE ───────────────────────────────
    log_phase("RANKING ENGINE")
    ranking_results = {}
    try:
        from models.ranking_engine import RankingEngine
        re = RankingEngine(config)
        target_scores = target_results.get("scored_targets", None)
        resistance_scores = resistance_results.get("scores", None)
        if compounds_processed is not None:
            ranking_results = re.run(compounds_processed, target_scores, resistance_scores)
            results["ranking"] = ranking_results
        else:
            print("  ⚠ Skipping ranking — no compound data")
    except Exception as e:
        errors.append(f"Ranking: {e}")
        traceback.print_exc()

    # ── PHASE 9: MANUSCRIPT GENERATOR ─────────────────────────
    log_phase("MANUSCRIPT GENERATOR")
    try:
        from generators.manuscript_generator import ManuscriptGenerator
        mg = ManuscriptGenerator(config)
        manuscript = mg.generate(results)
        results["manuscript_words"] = len(manuscript.split())
    except Exception as e:
        errors.append(f"Manuscript: {e}")
        traceback.print_exc()

    # ── PHASE 10: ICMR GRANT PACKAGE ─────────────────────────
    log_phase("ICMR GRANT PACKAGE")
    try:
        from generators.icmr_generator import ICMRGrantGenerator
        ig = ICMRGrantGenerator(config)
        proposal = ig.generate(results)
        results["icmr_words"] = len(proposal.split())
    except Exception as e:
        errors.append(f"ICMR: {e}")
        traceback.print_exc()

    # ── PHASE 11: IEC PROTOCOL ────────────────────────────────
    log_phase("IEC PROTOCOL")
    try:
        from generators.iec_generator import IECProtocolGenerator
        iec = IECProtocolGenerator(config)
        protocol = iec.generate(results)
        results["iec_words"] = len(protocol.split())
    except Exception as e:
        errors.append(f"IEC: {e}")
        traceback.print_exc()

    # ── SUMMARY ───────────────────────────────────────────────
    elapsed = time.time() - start_time
    print("\n" + "=" * 70)
    print("  PIPELINE COMPLETE")
    print("=" * 70)
    print(f"  ⏱ Total time: {elapsed:.1f} seconds")
    print(f"  ✓ Phases completed: {14 - len(errors)}/14")
    if errors:
        print(f"  ⚠ Errors ({len(errors)}):")
        for err in errors:
            print(f"    - {err}")
    else:
        print("  ✅ All phases completed successfully!")

    print("\n  📁 Output directories:")
    output_base = Path(config["paths"]["output_dir"])
    for subdir in ["omics", "epi", "targets", "cryptic_data", "de_novo", "compounds", "structures", "alphafold_structures", "docking", "models",
                    "resistance", "ranking", "manuscript", "icmr", "iec"]:
        d = output_base / subdir
        if d.exists():
            n_files = len(list(d.iterdir()))
            print(f"    {subdir}/  ({n_files} files)")

    # Save pipeline report
    report = {
        "elapsed_seconds": round(elapsed, 1),
        "phases_completed": 14 - len(errors),
        "phases_total": 14,
        "errors": errors,
        "status": "SUCCESS" if not errors else "PARTIAL",
    }
    ensure_dir(output_base)
    safe_save_json(report, output_base / "pipeline_report.json")

    return results, errors


if __name__ == "__main__":
    results, errors = run_pipeline()
    sys.exit(0 if not errors else 1)
