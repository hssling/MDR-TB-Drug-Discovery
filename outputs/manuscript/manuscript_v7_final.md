# Resistance-Informed Computational Discovery of a Novel InhA-Focused Candidate for Multidrug-Resistant *Mycobacterium tuberculosis*

**Authors:** Dr Siddalingaiah H S  
**Affiliation:** Professor, Community Medicine, Shridevi Institute of Medical Sciences and Research Hospital, Tumkur  
**Correspondence:** hssling@yahoo.com | +91 8941087719  
**ORCID:** 0000-0002-4771-8285  
**Target Journal:** *Journal of Cheminformatics*  
**Article Type:** Research Article  
**Running Title:** Computational MDR-TB lead discovery workflow

---

## Abstract

**Background:** Multidrug-resistant tuberculosis (MDR-TB) remains a major obstacle to tuberculosis control and motivates rapid, reproducible strategies for prioritizing new anti-tubercular lead compounds. Computational workflows that connect resistance knowledge, target prioritization, structure-based screening, and post-docking triage may reduce the size of the experimental search space.

**Methods:** We executed a 19-phase *in silico* workflow integrating omics and epidemiological summaries, resistance-gene mapping, target scoring, AlphaFold-guided structural analysis, *de novo* compound generation, docking, molecular-dynamics proxy analysis, ADMET prediction, human off-target similarity screening, quantum-chemical analysis, retrosynthetic planning, and machine-learning-assisted ranking. The current run scored 15 validated TB targets and ranked 55 compounds.

**Results:** InhA was the highest-ranked target with a composite score of 0.8695. The top-ranked compound, `MDR_AI_030`, achieved a docking score of -9.77 kcal/mol against InhA and remained favorable across orthogonal post-docking modules. Follow-up profiling showed RMSD 0.12 nm over a 10 ns proxy molecular-dynamics run with four persistent hydrogen bonds and two water bridges; predicted gastrointestinal absorption `Yes`; blood-brain barrier permeability `No`; hERG toxicity risk `Low`; TPSA 0.0; QED 0.467; maximum Tanimoto similarity 0.146 to the human-toxin reference set; HOMO/LUMO energies of -5.5 and -2.6 eV, respectively; a band gap of 2.9 eV; and a one-step Sonogashira cross-coupling retrosynthetic route.

**Scientific Contribution:** This study presents a reproducible, resistance-aware computational manuscript package that distinguishes clearly between exported pipeline outputs and interpretive claims, and it documents an end-to-end workflow for MDR-TB lead nomination grounded in saved artefacts.

**Conclusions:** `MDR_AI_030` should be interpreted as a computationally prioritized lead rather than a validated therapeutic agent. Nonetheless, its integrated docking, stability-proxy, selectivity, quantum, and retrosynthetic profile supports continued follow-up through expanded simulations and experimental evaluation.

**Keywords:** MDR-TB, *Mycobacterium tuberculosis*, InhA, AlphaFold, generative chemistry, docking, ADMET, retrosynthesis

---

## 1. Introduction

Tuberculosis remains one of the leading infectious causes of death worldwide, and the persistence of drug-resistant *Mycobacterium tuberculosis* strains continues to erode the effectiveness of standard treatment regimens [1,2]. MDR-TB, defined by resistance to at least rifampicin and isoniazid, substantially complicates treatment by increasing regimen complexity, duration, toxicity burden, and the likelihood of treatment failure [1,2]. These clinical pressures have intensified interest in workflows that can prioritize targets and compounds more efficiently before costly experimental work begins.

The current generation of computational resources makes such workflows increasingly feasible. Resistance catalogues such as CRyPTIC organize genotype-phenotype relationships at scale and provide a structured view of clinically relevant mutations [5]. Structure prediction advances, particularly AlphaFold, provide plausible starting geometries for protein targets when experimentally determined structures are incomplete or unavailable [6]. In parallel, generative chemistry approaches can produce chemically novel candidate structures that extend beyond fixed reference libraries [4]. When these resources are integrated with ranking, docking, and post-docking profiling, they can support a coherent strategy for early-stage antimicrobial lead prioritization.

However, many computational reports overstate what their outputs demonstrate. Docking scores are not potency measurements, short or proxy molecular-dynamics summaries are not substitutes for longer physics-based campaigns, and predictive ADMET or polypharmacology signals do not establish biological safety. For publication-ready reporting, the distinction between computational nomination and validated therapeutic efficacy must therefore be explicit. A manuscript that is traceable to the exact output artefacts is more useful than one that relies on inflated language or disconnected claims.

The aim of the present study was to assemble, audit, and report an end-to-end computational discovery workflow for MDR-TB that remains tightly grounded in exported outputs. The pipeline integrates disease-context modules, resistance-informed target prioritization, compound generation, docking, machine-learning support, and post-docking triage. The resulting manuscript is intended to document a reproducible *in silico* lead-nomination process rather than to claim biological validation.

---

## 2. Materials and Methods

### 2.1 Study design

This study was conducted entirely *in silico*. The workflow combined data loading, transcriptomic and epidemiological summarization, resistance mapping, target scoring, compound generation, docking, machine-learning analysis, ADMET prediction, molecular-dynamics proxy analysis, human off-target screening, quantum-chemical estimation, retrosynthesis, and manuscript generation. No biological experiments, synthesis, animal work, or clinical interventions were performed.

### 2.2 Input data and fallback behavior

The pipeline consumes a mixture of public, derived, and simulated computational inputs. In the current project configuration, omics modules can fall back to synthetic expression matrices when live access or complete upstream data are unavailable. Epidemiological summaries were generated across 10 Indian regions covering the years 2015-2023. Resistance analysis covered 10 MDR-associated genes and 40 cataloged mutations. Compound collections combined the baseline library with *de novo* generated structures, producing 55 ranked compounds in the current run.

Because parts of the workflow can operate with fallback or simulated data, the present manuscript is framed as a computational prioritization study. Throughout the text, findings are described as generated, predicted, ranked, or prioritized rather than experimentally confirmed.

### 2.3 Omics and pathway analysis

The omics module tested 5,000 genes and identified significant differential-expression events under the configured thresholds. The current run reported 230 significantly upregulated and 196 significantly downregulated features, giving 426 significant events in total. Pathway scoring was exported to `outputs/omics/pathway_scores.csv`, where the highest enrichment score corresponded to `mtu00010_Glycolysis`.

### 2.4 Epidemiological contextualization

The epidemiology module summarized trends across 10 regions over 2015-2023. These outputs were used to contextualize disease burden and resistance trends rather than to make causal epidemiological claims. The exported summary reported an overall TB incidence trend of `Decreasing`, an MDR-TB trend of `Rising`, a highest-burden region of `Rajasthan`, and a mean MDR-TB percentage of 8.75 across the analyzed regional set.

### 2.5 Resistance-aware target prioritization

Resistance mapping catalogued mutation frequency, confidence, and associated drugs across 10 genes. In parallel, 15 validated TB targets were scored using a composite framework that incorporated druggability, essentiality, conservation, resistance association, and an expression-related term. Output tables were written to `outputs/targets/scored_targets.csv` and `outputs/resistance/resistance_scores.csv`.

### 2.6 Compound generation, ranking, and docking

The workflow processed 55 compounds, including *de novo* generated structures. The ranking module combined predicted activity, Lipinski-style compliance, target relevance, structural novelty, and resistance-related scoring into a final ranking score. Docking against InhA was exported to `outputs/docking/docking_results_InhA.csv`.

### 2.7 Post-docking triage modules

The top-ranked docked compound was carried forward into multiple post-docking modules:

- **Molecular-dynamics proxy analysis:** summarized RMSD, persistent hydrogen bonds, and water bridges over a 10 ns simulated proxy run;
- **ADMET prediction:** reported TPSA, gastrointestinal absorption, blood-brain barrier permeability, hERG risk, and QED-like scoring;
- **Polypharmacology screening:** reported fingerprint-based similarity against a human-toxin reference panel;
- **Quantum-chemical analysis:** reported HOMO, LUMO, and band-gap estimates; and
- **Retrosynthetic planning:** proposed plausible route decomposition and reaction class.

These modules were used as triage tools, not as substitutes for experimental verification.

### 2.8 Machine-learning analysis

Machine-learning outputs were derived from a 100-sample, 6-feature classification workflow using an 80/20 train-test split. Three algorithms were evaluated: random forest, gradient boosting, and logistic regression. The best-performing model in the saved output was the random forest with ROC-AUC 0.7273.

### 2.9 Reproducibility and AI-assisted drafting disclosure

All values reported in this manuscript were checked against saved CSV and JSON outputs in the workspace. AI-assisted coding and drafting tools were used under human supervision to update the manuscript and supporting generator logic. Final reported values, wording, and interpretations were reviewed against the exported artefacts by the human author, who retains full responsibility for the content.

---

## 3. Results

### 3.1 Omics, pathway, and epidemiological context

The omics module identified 426 significant differential-expression events across 5,000 genes. The top five exported pathway scores were `mtu00010_Glycolysis` (8.6509), `mtu00020_TCA_cycle` (3.7137), `mtu00190_Oxidative_phosphorylation` (2.9899), `mtu00230_Purine_metabolism` (2.8079), and `mtu00240_Pyrimidine_metabolism` (2.4112). These outputs suggest a metabolically perturbed state in the current computational run and provide a disease-context backdrop for downstream prioritization.

The epidemiology summary contextualized these results within a region-level trend analysis covering 2015-2023. The overall incidence trend was reported as decreasing, whereas the MDR-TB trend was rising. Rajasthan emerged as the highest-burden region in the current run, and the mean MDR-TB proportion across the analyzed regional set was 8.75%. Although these values are not used as causal inputs to compound-level scoring, they reinforce the rationale for emphasizing resistance-aware discovery.

### 3.2 Resistance scoring and target prioritization

The resistance module covered 10 genes, 40 cataloged mutations, and 10 drugs, with an average resistance score of 0.5085. The highest-scoring resistance-associated genes were `rpoB` (Rifampicin, 0.635), `katG` (Isoniazid, 0.620), `rrs` (Aminoglycosides, 0.575), `gyrA` (Fluoroquinolones, 0.565), and `inhA` (Isoniazid low-level resistance, 0.525). These findings are consistent with the central role of canonical MDR-associated loci in determining treatment constraints.

Target prioritization identified InhA as the top-ranked target with final score 0.8695. Other highly ranked targets included RpoB, GyrA, EmbB, and Rrs. The top target list is shown in Table 1.

**Table 1. Top-ranked targets in the current run**

| Rank | Target | Category | Final score | Priority |
| :--- | :--- | :--- | :--- | :--- |
| 1 | InhA | Cell wall | 0.8695 | High |
| 2 | RpoB | Transcription | 0.8680 | High |
| 3 | GyrA | DNA replication | 0.8650 | High |
| 4 | EmbB | Cell wall | 0.8255 | High |
| 5 | Rrs | Translation | 0.8110 | High |

The emergence of InhA as the leading target is noteworthy because it aligns a biologically established anti-TB target with the resistance-aware prioritization framework used here. This made InhA an appropriate focus for the downstream docking and post-docking workflow.

### 3.3 Compound ranking and docking

The ranking engine evaluated 55 compounds. The top-ranked compound was `MDR_AI_030` with a final ranking score of 0.736. Although the summary file reported zero compounds in a strict `n_high_priority` bucket, the ranked table still identified multiple compounds as `Promising`, with `MDR_AI_030` ranked first.

Docking results against InhA showed a leading cluster of compounds with substantially stronger predicted affinities than the repeated -4.0 kcal/mol plateau seen among many lower-ranked compounds. The top 10 docked compounds are shown in Table 2.

**Table 2. Top docking results against InhA**

| Rank | Compound ID | Binding affinity (kcal/mol) | Interacting residues | Confidence |
| :--- | :--- | :--- | :--- | :--- |
| 1 | MDR_AI_030 | -9.77 | SER-315 | High |
| 2 | MDR_AI_047 | -9.38 | HIS-132, TYR-158 | High |
| 3 | MDR_AI_026 | -8.84 | SER-315, TRP-222, TYR-158 | Medium |
| 4 | MDR_AI_010 | -8.68 | TRP-222 | Medium |
| 5 | MDR_AI_039 | -8.67 | ASP-94 | Medium |
| 6 | MDR_AI_033 | -8.53 | ASN-426 | Medium |
| 7 | MDR_AI_005 | -8.06 | TRP-222 | Medium |
| 8 | MDR_AI_019 | -8.03 | ASP-94, TYR-158 | Medium |
| 9 | MDR_AI_006 | -7.68 | SER-315, TRP-222, TYR-158 | Medium |
| 10 | MDR_AI_042 | -7.65 | TYR-158 | Medium |

Within this ranking set, `MDR_AI_030` combined the strongest docking score with a favorable downstream triage profile, making it the most coherent lead for continued analysis in the current run.

### 3.4 Post-docking triage of `MDR_AI_030`

`MDR_AI_030` was examined across all downstream triage modules. The molecular-dynamics summary reported RMSD 0.12 nm over a 10 ns simulated proxy run, together with four persistent hydrogen bonds and two water bridges. While this should not be interpreted as a substitute for a long-timescale molecular-dynamics campaign, it is consistent with a stable pose in the context of the present proxy analysis.

ADMET prediction returned TPSA 0.0, gastrointestinal absorption `Yes`, blood-brain barrier permeability `No`, hERG toxicity risk `Low`, and QED 0.467. These values suggest a lead that remains worthy of further follow-up, while also underscoring the need for proper experimental ADMET validation. The blood-brain barrier result in particular argues against describing the compound as broadly permeant or making unsupported CNS-related claims.

Polypharmacology screening produced a maximum Tanimoto similarity of 0.146 to the human-toxin reference set, which is consistent with low similarity in the context of the implemented screening panel. Again, this is a screening result rather than a direct safety claim.

Quantum-chemical analysis estimated HOMO and LUMO energies of -5.5 eV and -2.6 eV, respectively, giving a band gap of 2.9 eV. These descriptors provide an additional computational perspective on the electronic character of the lead but should be interpreted descriptively rather than mechanistically overstated.

Retrosynthesis identified a one-step Sonogashira cross-coupling route using the fragments `Brc1ccccc1` and `C#Cc1cc(C#CC(F)(F)F)cc(C(F)(F)F)c1`. This route proposal is valuable because it gives the nominated lead a plausible synthetic entry point rather than leaving it as an entirely abstract structure.

**Table 3. Integrated post-docking profile of `MDR_AI_030`**

| Module | Output | Interpretation |
| :--- | :--- | :--- |
| Docking | -9.77 kcal/mol | Strongest predicted binder in the current run |
| MD proxy | RMSD 0.12 nm; 4 persistent H-bonds; 2 water bridges | Stable pose in a short proxy simulation |
| ADMET | TPSA 0.0; GI Yes; BBB No; hERG Low; QED 0.467 | Early triage profile suitable for follow-up, not validation |
| Polypharmacology | Max Tanimoto similarity 0.146 | Low similarity to the reference human-toxin panel |
| Quantum chemistry | HOMO -5.5 eV; LUMO -2.6 eV; band gap 2.9 eV | Descriptive computational electronic profile |
| Retrosynthesis | 1 step; Sonogashira cross-coupling | Plausible synthetic entry point |

### 3.5 Machine-learning performance

Machine-learning performance was modest but adequate for prioritization support. The best model was the random forest with ROC-AUC 0.7273, accuracy 0.650, precision 0.600, recall 0.667, and F1-score 0.632. Gradient boosting showed intermediate performance, while logistic regression performed substantially worse on the current task.

**Table 4. Model performance summary**

| Algorithm | Accuracy | Precision | Recall | F1-score | ROC-AUC |
| :--- | :--- | :--- | :--- | :--- | :--- |
| random_forest | 0.650 | 0.600 | 0.667 | 0.632 | 0.727 |
| gradient_boosting | 0.600 | 0.556 | 0.556 | 0.556 | 0.626 |
| logistic_regression | 0.350 | 0.167 | 0.111 | 0.133 | 0.343 |

These results support the role of the machine-learning module as a ranking aid rather than a standalone decision engine.

---

## 4. Discussion

The central value of the current workflow lies in integration rather than in any single module. A lead that performs well only in docking can still fail triage when evaluated for stability, selectivity, or synthetic accessibility. Conversely, a compound with a modest docking score may become attractive when evaluated in a broader framework. The present pipeline addresses this by linking disease context, resistance biology, target prioritization, docking, machine-learning support, and post-docking profiling in one auditable chain.

InhA emerged as the leading target from the current scoring framework. This is biologically plausible because InhA is a well-established component of mycolic-acid biosynthesis and is directly relevant to isoniazid-related resistance biology [3,5]. The prioritization of InhA therefore provides a rational bridge between resistance-aware target scoring and structure-guided compound screening.

`MDR_AI_030` similarly stands out because it remained favorable across multiple independent outputs rather than only one. It ranked first overall, docked most strongly among the screened compounds, maintained a stable short proxy MD profile, showed low similarity to the human-toxin reference set, and remained synthetically interpretable through a simple reaction proposal. In computational lead nomination, this type of multi-module coherence is a stronger reason for prioritization than a single extreme score in isolation.

At the same time, all of these outputs remain predictive. Docking scores are not biochemical affinities. The current molecular-dynamics module reports a short, proxy-style summary and does not substitute for a full production simulation campaign. ADMET and polypharmacology findings are predictive filters rather than empirical safety or pharmacokinetic measurements. The retrosynthetic module proposes a plausible route but does not demonstrate laboratory success. These distinctions matter, especially when a manuscript is intended to be publication-ready.

Another important point is the role of fallback or simulated data. The current project architecture allows progress even when complete live datasets are unavailable, which is useful for workflow demonstration and software validation. However, it also means that biological realism depends on upstream data availability and quality. This limitation must be made explicit rather than hidden behind overstated prose.

The present manuscript therefore emphasizes traceability. All reported values were pulled from saved project artefacts, and narrative claims were constrained to what those files support. This approach is particularly important for auto-generated or AI-assisted manuscripts, where drift between the data and the prose can occur quickly if output files are not rechecked systematically.

---

## 5. Limitations

This work has several limitations. First, the workflow can operate with simulated or fallback data, so the biological realism of some upstream components is limited by data availability. Second, several modules are heuristic or proxy-based rather than full experimental or high-fidelity computational replacements. Third, the study does not include external prospective validation against independent MDR-TB screening datasets. Fourth, no wet-lab synthesis, enzymatic assay, MIC testing, or cell-based validation was performed. Finally, despite being manuscript-ready in structure and traceability, the text should still undergo domain-expert review before formal submission.

---

## 6. Conclusions

This end-to-end computational workflow prioritized InhA as the leading target and identified `MDR_AI_030` as the top-ranked lead in the current run. The integrated docking, MD-proxy, ADMET, human off-target similarity, quantum, and retrosynthetic outputs justify continued follow-up. More broadly, the project provides a reproducible, resistance-aware framework for computational MDR-TB lead nomination and a manuscript package grounded directly in exported results.

---

## Declarations

**Ethics approval and consent to participate:** Not applicable. No human participants, animals, or wet-lab experiments were involved.

**Consent for publication:** Not applicable.

**Availability of data and materials:** All reported values are traceable to saved output files in the project workspace, including `outputs/docking/docking_results_InhA.csv`, `outputs/md_simulations/md_summary.csv`, `outputs/admet/admet_toxicity_report.csv`, `outputs/polypharmacology/human_offtarget_report.csv`, `outputs/quantum_mechanics/electronic_orbitals.csv`, `outputs/retrosynthesis/synthesis_routes.csv`, `outputs/targets/scored_targets.csv`, `outputs/ranking/ranked_compounds.csv`, `outputs/omics/omics_summary.json`, and `outputs/resistance/resistance_scores.csv`.

**Code availability:** Manuscript generation and pipeline orchestration are implemented in the local project codebase, including `run_pipeline.py` and `generators/manuscript_generator.py`.

**Competing interests:** The author declares no competing interests.

**Funding:** No external funding was reported for this computational study.

**Authors' contributions:** Dr Siddalingaiah H S conceived the study, reviewed the pipeline outputs, validated the manuscript against exported artefacts, and approved the final text. AI-assisted tooling supported drafting and code editing under human supervision.

**Acknowledgements:** The project builds on public knowledge resources and open computational ecosystems for cheminformatics, structural bioinformatics, and machine learning.

---

## References

1. World Health Organization. Global tuberculosis report 2025. Geneva: World Health Organization; 2025.
2. Dheda K, Gumbo T, Maartens G, Dooley KE, McNerney R, Murray M, et al. The epidemiology, pathogenesis, transmission, diagnosis, and management of multidrug-resistant, extensively drug-resistant, and incurable tuberculosis. Lancet Respir Med. 2017;5(4):291-360.
3. Vilchèze C, Jacobs WR Jr. The mechanism of isoniazid killing: clarity through the scope of genetics. Annu Rev Microbiol. 2007;61:35-50.
4. Stokes JM, Yang K, Swanson K, Jin W, Cubillos-Ruiz A, Donghia NM, et al. A deep learning approach to antibiotic discovery. Cell. 2020;180(4):688-702.
5. The CRyPTIC Consortium. A data compendium associating the genomes of 12,289 *Mycobacterium tuberculosis* isolates with quantitative resistance phenotypes to 13 antituberculosis drugs. PLoS Biol. 2022;20(8):e3001721.
6. Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O, et al. Highly accurate protein structure prediction with AlphaFold. Nature. 2021;596(7873):583-589.
7. Koes DR, Baumgartner MP, Camacho CJ. Lessons learned in empirical scoring with smina from the CSAR 2011 benchmarking exercise. J Chem Inf Model. 2013;53(8):1893-1904.
8. Eastman P, Swails J, Chodera JD, McGibbon RT, Zhao Y, Beauchamp KA, et al. OpenMM 7: rapid development of high performance algorithms for molecular dynamics. PLoS Comput Biol. 2017;13(7):e1005659.
9. Sun Q, Berkelbach TC, Blunt NS, Booth GH, Chen S, Churchill EN, et al. PySCF: the Python-based simulations of chemistry framework. Wiley Interdiscip Rev Comput Mol Sci. 2018;8(1):e1340.
