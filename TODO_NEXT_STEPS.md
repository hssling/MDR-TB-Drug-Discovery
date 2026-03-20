# TODO / Handoff

## Completed

- Ran genuine drug discovery pipeline using real ChEMBL1849 IC50 data (277 compounds, 108 active).
- Trained real QSAR models (RF, GB, LR) on Morgan ECFP4 fingerprints; best model LR test AUC 0.961.
- Computed real RDKit ADMET for all ChEMBL1849 compounds.
- Performed Tanimoto LBVS against 5 known InhA inhibitor scaffolds.
- Ranked all compounds by composite score (pIC50_norm + QSAR_prob + Tanimoto + QED).
- **Lead compound identified:** CHEMBL3125270 (IC50 = 4 nM, measured; composite score 0.8245)
  - Scaffold: pyrazole-benzofuran-pyrrolidine
  - MW 423.47 Da, cLogP 1.71, TPSA 123.46, HBD 2, HBA 5, QED 0.625
  - Lipinski pass, Veber pass, hERG low risk
  - QSAR active probability 0.957
- Wrote genuine manuscript `outputs/manuscript/manuscript_v8_genuine.md` (no fake data).
- Replaced all old manuscript files (manuscript_v7_final.md, manuscript_imrad.md) with v8 content.
- Deleted old fake manuscript files (academic_paper_v6.md, manuscript_v6_expanded.md).
- Rewrote supplementary_materials.md with real ChEMBL/RDKit/CRyPTIC data.
- Retracted fabricated output files (md_simulations, quantum_mechanics, de_novo).
- Updated ranking_summary.json to reflect real compound (CHEMBL3125270).
- Marked deprecated engines (md_engine.py, qmmm_engine.py) with DEPRECATED warnings.
- Generated 7 genuine figures using real data (no np.random for scientific values).
- Downloaded real InhA crystal structure PDB 4TZK from RCSB.
- Removed all MDR_AI_030 references from manuscript, supplementary, rankings, and code defaults.

## Genuine Drug Discovery Summary

- **Target:** InhA (Rv1484, CHEMBL1849) — enoyl-ACP reductase, cell-wall biosynthesis
  - Rationale: katG S315T mutation (freq 0.68 in MDR strains) abolishes isoniazid prodrug activation but leaves InhA intact; direct InhA inhibitors bypass this resistance
  - Crystal structure: PDB 4TZK (1.65 A resolution, NAD+ complex)
- **Dataset:** 277 compounds with measured IC50 from ChEMBL1849
- **Lead compound:** CHEMBL3125270 (IC50 = 4 nM, measured; composite score 0.8245)
  - SMILES: CCc1cc(C(=O)N[C@@H]2C[C@@H](C(N)=O)N(C(=O)c3coc4ccccc34)C2)n(CC)n1
  - All ADMET predicted computationally (not experimentally validated)
- **QSAR:** Best model LR, test AUC 0.961, F1 0.821 on real ChEMBL IC50 data
- **Scope limitations:** No docking (AutoDock Vina unavailable), no MD, no QM performed

## Current Key Files

- Canonical manuscript: `outputs/manuscript/manuscript_v8_genuine.md`
  - Also copied to: `manuscript_v7_final.md`, `manuscript_imrad.md` (same content)
- Supplementary: `outputs/manuscript/supplementary_materials.md` (real data)
- Figures: `outputs/figures/` (figure_1 through figure_7, all genuine)
- Top 10 compounds: `outputs/ranking/top_10_compounds.csv` (real ChEMBL data)
- QSAR performance: `outputs/models/qsar_performance.json` (real metrics)
- InhA structure: `outputs/structures/4TZK.pdb` (RCSB download)
- Genuine pipeline script: `scripts/genuine_drug_discovery.py`
- Figure generator: `scripts/generate_manuscript_figures.py`
- DOCX export: `scripts/export_submission_docx.py`

## Next Steps

1. Rebuild DOCX submission assets (run `python scripts/export_submission_docx.py`):
   - Ensure `outputs/submission/submission_manuscript.docx` is not open in Word first
2. Optionally generate graphical abstract or highlights for journal submission
3. Push updated package to Hugging Face / Kaggle if desired

## Notes

- Figure generator uses `/c/Python314/python.exe` on this machine (Python 3.14 + rdkit 2025.9.6).
- All computational claims are ligand-based (QSAR + Tanimoto similarity); no structure-based docking.
- All IC50 values are experimentally measured from ChEMBL1849 (not predicted).
- Manuscript explicitly acknowledges all computational limitations.
