"""RDKit-based computational tools — ADMET, fingerprints, QSAR, ranking. All local, no network."""
from __future__ import annotations

import json
from .base import BaseTool, ToolResult, Citation


_RDKIT_CITATION = Citation.now(
    source="RDKit: Open-Source Cheminformatics",
    url="https://www.rdkit.org",
    doi=""
)


class RDKitTools(BaseTool):

    @classmethod
    def get_tool_definitions(cls) -> list[dict]:
        return [
            {
                "name": "compute_admet",
                "description": (
                    "Compute physicochemical and ADMET properties for a list of SMILES strings using RDKit. "
                    "Returns MW, cLogP, TPSA, HBD, HBA, rotatable bonds, QED, "
                    "Lipinski rule-of-5 compliance, Veber oral bioavailability compliance, "
                    "GI absorption prediction, BBB penetration prediction, and hERG structural alert flag. "
                    "All values are computed locally — no external API call."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "smiles_list": {
                            "type": "array",
                            "items": {"type": "string"},
                            "description": "List of SMILES strings",
                            "minItems": 1,
                            "maxItems": 500
                        },
                        "ids": {
                            "type": "array",
                            "items": {"type": "string"},
                            "description": "Optional parallel list of compound IDs (e.g. ChEMBL IDs)"
                        }
                    },
                    "required": ["smiles_list"]
                }
            },
            {
                "name": "train_qsar_model",
                "description": (
                    "Train a QSAR binary classifier on a list of compounds with measured IC50 values. "
                    "Uses Morgan ECFP4 fingerprints and evaluates Random Forest, Gradient Boosting, "
                    "and Logistic Regression classifiers with 5-fold cross-validation. "
                    "Returns cross-validated and held-out ROC-AUC, F1, precision, recall, and accuracy. "
                    "Requires at least 30 compounds with SMILES and IC50 values."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "compounds": {
                            "type": "array",
                            "description": "List of compounds with SMILES and IC50",
                            "items": {
                                "type": "object",
                                "properties": {
                                    "id": {"type": "string"},
                                    "smiles": {"type": "string"},
                                    "ic50_nm": {"type": "number"}
                                },
                                "required": ["smiles", "ic50_nm"]
                            }
                        },
                        "activity_threshold_nm": {
                            "type": "number",
                            "default": 1000.0,
                            "description": "IC50 below this value (nM) is classified as active"
                        }
                    },
                    "required": ["compounds"]
                }
            },
            {
                "name": "rank_compounds",
                "description": (
                    "Rank compounds by composite score: "
                    "pIC50 (40%) + Tanimoto LBVS similarity (30%) + QED drug-likeness (20%) + Lipinski compliance (10%). "
                    "All IC50 values must be real measured values. "
                    "Returns a ranked table with composite scores and individual components."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "compounds": {
                            "type": "array",
                            "items": {
                                "type": "object",
                                "properties": {
                                    "id": {"type": "string"},
                                    "smiles": {"type": "string"},
                                    "ic50_nm": {"type": "number"}
                                },
                                "required": ["id", "smiles", "ic50_nm"]
                            }
                        },
                        "reference_smiles": {
                            "type": "array",
                            "items": {"type": "string"},
                            "description": "Known active reference scaffold SMILES for LBVS Tanimoto scoring"
                        },
                        "top_n": {
                            "type": "integer",
                            "default": 10,
                            "minimum": 1,
                            "maximum": 50
                        }
                    },
                    "required": ["compounds"]
                }
            }
        ]

    # ── Internal helpers ──────────────────────────────────────────────────────

    @staticmethod
    def _compute_mol_props(smiles: str) -> dict | None:
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, QED, rdMolDescriptors, AllChem

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            mw = Descriptors.ExactMolWt(mol)
            clogp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
            rings = rdMolDescriptors.CalcNumRings(mol)
            ar_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
            qed = QED.qed(mol)
            frac_csp3 = rdMolDescriptors.CalcFractionCSP3(mol)

            lipinski = (mw <= 500 and clogp <= 5 and hbd <= 5 and hba <= 10)
            veber = (rotb <= 10 and tpsa <= 140)
            gi_abs = tpsa < 140 and hbd <= 5
            bbb = tpsa < 90 and mw < 450
            herg_risk = clogp > 4.0 and ar_rings >= 2

            return {
                "mw": round(mw, 2),
                "clogp": round(clogp, 3),
                "tpsa": round(tpsa, 2),
                "hbd": hbd,
                "hba": hba,
                "rotatable_bonds": rotb,
                "rings": rings,
                "aromatic_rings": ar_rings,
                "frac_csp3": round(frac_csp3, 3),
                "qed": round(qed, 3),
                "lipinski_compliant": lipinski,
                "veber_compliant": veber,
                "gi_absorption": gi_abs,
                "bbb_penetration": bbb,
                "herg_risk": herg_risk,
                "parseable": True,
            }
        except ImportError:
            raise
        except Exception:
            return None

    @staticmethod
    def _morgan_fp(smiles: str, radius: int = 2, n_bits: int = 2048):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)

    @staticmethod
    def _tanimoto(fp1, fp2) -> float:
        from rdkit.DataStructs import TanimotoSimilarity
        return TanimotoSimilarity(fp1, fp2)

    # ── Tool implementations ──────────────────────────────────────────────────

    def compute_admet(self, smiles_list: list[str],
                      ids: list[str] | None = None) -> ToolResult:
        try:
            from rdkit import Chem  # noqa: F401
        except ImportError:
            return self._err("compute_admet", "RDKit not installed.")

        ids = ids or [f"cpd_{i+1}" for i in range(len(smiles_list))]
        results = []
        failed = []

        for cid, smi in zip(ids, smiles_list):
            props = self._compute_mol_props(smi)
            if props is None:
                failed.append(cid)
            else:
                props["id"] = cid
                props["smiles"] = smi
                results.append(props)

        if not results:
            return self._err("compute_admet", "All SMILES failed to parse with RDKit.")

        n = len(results)
        lipinski_pass = sum(1 for r in results if r.get("lipinski_compliant"))
        warnings = []
        if failed:
            warnings.append(f"{len(failed)} SMILES could not be parsed: {failed[:5]}")

        return self._ok(
            tool_name="compute_admet",
            data={
                "total_input": len(smiles_list),
                "successfully_computed": n,
                "failed_to_parse": len(failed),
                "lipinski_compliant_count": lipinski_pass,
                "compounds": results,
            },
            citations=[_RDKIT_CITATION],
            confidence=0.82,
            rationale="RDKit rule-based ADMET descriptors validated against Lipinski/Veber/Clark criteria. "
                      "Not PBPK or experimental measurement — experimental confirmation required.",
            warnings=warnings,
        )

    def train_qsar_model(self, compounds: list[dict],
                         activity_threshold_nm: float = 1000.0) -> ToolResult:
        try:
            import numpy as np
            import pandas as pd
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
            from sklearn.linear_model import LogisticRegression
            from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
            from sklearn.metrics import roc_auc_score, f1_score, accuracy_score, precision_score, recall_score
        except ImportError as e:
            return self._err("train_qsar_model", f"Missing dependency: {e}")

        if len(compounds) < 30:
            return self._err("train_qsar_model",
                             f"Only {len(compounds)} compounds provided. Need ≥30 for QSAR.")

        # Build fingerprint matrix
        fps = []
        labels = []
        valid_ids = []
        for cpd in compounds:
            smi = cpd.get("smiles", "")
            ic50 = cpd.get("ic50_nm")
            cid = cpd.get("id", "")
            if not smi or ic50 is None:
                continue
            fp = self._morgan_fp(smi)
            if fp is None:
                continue
            fps.append(list(fp))
            labels.append(1 if float(ic50) < activity_threshold_nm else 0)
            valid_ids.append(cid)

        if len(fps) < 30:
            return self._err("train_qsar_model",
                             f"Only {len(fps)} compounds had parseable SMILES. Need ≥30.")

        X = np.array(fps)
        y = np.array(labels)
        active_count = int(y.sum())
        inactive_count = len(y) - active_count

        if active_count < 5 or inactive_count < 5:
            return self._err("train_qsar_model",
                             f"Too few examples in one class: active={active_count}, inactive={inactive_count}. Need ≥5 each.")

        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.20, random_state=42, stratify=y
        )

        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        models = {
            "Random Forest": RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=42, n_jobs=-1),
            "Gradient Boosting": GradientBoostingClassifier(n_estimators=200, learning_rate=0.1, random_state=42),
            "Logistic Regression": LogisticRegression(C=0.1, max_iter=1000, class_weight="balanced", random_state=42),
        }

        model_results = {}
        best_model = None
        best_auc = -1.0

        for name, clf in models.items():
            cv_aucs = cross_val_score(clf, X_train, y_train, cv=cv, scoring="roc_auc")
            cv_f1s = cross_val_score(clf, X_train, y_train, cv=cv, scoring="f1")
            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)
            try:
                y_prob = clf.predict_proba(X_test)[:, 1]
                test_auc = roc_auc_score(y_test, y_prob)
            except Exception:
                test_auc = 0.0

            model_results[name] = {
                "cv_roc_auc_mean": round(float(cv_aucs.mean()), 4),
                "cv_roc_auc_std": round(float(cv_aucs.std()), 4),
                "cv_f1_mean": round(float(cv_f1s.mean()), 4),
                "test_roc_auc": round(test_auc, 4),
                "test_f1": round(float(f1_score(y_test, y_pred, zero_division=0)), 4),
                "test_accuracy": round(float(accuracy_score(y_test, y_pred)), 4),
                "test_precision": round(float(precision_score(y_test, y_pred, zero_division=0)), 4),
                "test_recall": round(float(recall_score(y_test, y_pred, zero_division=0)), 4),
            }
            if test_auc > best_auc:
                best_auc = test_auc
                best_model = name

        confidence = min(0.95, 0.50 + min(len(fps), 200) / 400 + best_auc * 0.20)
        warnings = []
        if len(fps) < 100:
            warnings.append(f"Small dataset ({len(fps)} compounds) — cross-validation estimates may be optimistic.")
        if best_auc < 0.70:
            warnings.append(f"Best test AUC {best_auc:.3f} is low — predictions should be treated with caution.")

        return self._ok(
            tool_name="train_qsar_model",
            data={
                "n_compounds": len(fps),
                "n_active": active_count,
                "n_inactive": inactive_count,
                "activity_threshold_nm": activity_threshold_nm,
                "fingerprint": "Morgan ECFP4 (radius=2, 2048 bits)",
                "best_model": best_model,
                "best_test_roc_auc": round(best_auc, 4),
                "models": model_results,
                "train_size": len(X_train),
                "test_size": len(X_test),
            },
            citations=[
                _RDKIT_CITATION,
                self._cite(
                    source="scikit-learn: Machine Learning in Python",
                    url="https://scikit-learn.org",
                    doi="10.5555/1953048.2078195"
                ),
                self._cite(
                    source="Rogers & Hahn (2010) Extended-Connectivity Fingerprints",
                    url="https://pubs.acs.org/doi/10.1021/ci100050t",
                    doi="10.1021/ci100050t"
                )
            ],
            confidence=confidence,
            rationale=f"QSAR trained on {len(fps)} real measured values. Best test ROC-AUC: {best_auc:.3f} ({best_model}). "
                      "Confidence scales with dataset size and model performance.",
            warnings=warnings,
        )

    def rank_compounds(self, compounds: list[dict],
                       reference_smiles: list[str] | None = None,
                       top_n: int = 10) -> ToolResult:
        try:
            import math
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            return self._err("rank_compounds", "RDKit not installed.")

        if not compounds:
            return self._err("rank_compounds", "No compounds provided.")

        # Compute ADMET + LBVS
        ref_fps = []
        if reference_smiles:
            for smi in reference_smiles:
                fp = self._morgan_fp(smi)
                if fp is not None:
                    ref_fps.append(fp)

        ic50_vals = [float(c["ic50_nm"]) for c in compounds if c.get("ic50_nm") and float(c["ic50_nm"]) > 0]
        if not ic50_vals:
            return self._err("rank_compounds", "No valid IC50 values found in compounds.")
        pic50_vals = [-math.log10(v * 1e-9) for v in ic50_vals]
        pic50_min, pic50_max = min(pic50_vals), max(pic50_vals)
        pic50_range = max(pic50_max - pic50_min, 0.001)

        scored = []
        for cpd in compounds:
            smi = cpd.get("smiles", "")
            ic50 = cpd.get("ic50_nm")
            cid = cpd.get("id", "")
            if not smi or ic50 is None:
                continue

            props = self._compute_mol_props(smi)
            if props is None:
                continue

            ic50_f = float(ic50)
            if ic50_f <= 0:
                continue

            # Normalized pIC50
            pic50 = -math.log10(ic50_f * 1e-9)
            pic50_norm = (pic50 - pic50_min) / pic50_range

            # Tanimoto LBVS
            fp = self._morgan_fp(smi)
            if fp and ref_fps:
                lbvs = max(self._tanimoto(fp, rp) for rp in ref_fps)
            else:
                lbvs = 0.0

            # Lipinski score (0 or 1)
            lip = 1.0 if props.get("lipinski_compliant") else 0.0

            # Composite
            composite = (0.40 * pic50_norm + 0.30 * lbvs +
                         0.20 * props.get("qed", 0.0) + 0.10 * lip)

            scored.append({
                "rank": 0,
                "id": cid,
                "smiles": smi,
                "ic50_nm": ic50_f,
                "pic50": round(pic50, 3),
                "pic50_norm": round(pic50_norm, 4),
                "lbvs_tanimoto": round(lbvs, 4),
                "qed": props.get("qed"),
                "lipinski_compliant": props.get("lipinski_compliant"),
                "mw": props.get("mw"),
                "clogp": props.get("clogp"),
                "tpsa": props.get("tpsa"),
                "hbd": props.get("hbd"),
                "hba": props.get("hba"),
                "herg_risk": props.get("herg_risk"),
                "gi_absorption": props.get("gi_absorption"),
                "composite_score": round(composite, 4),
            })

        if not scored:
            return self._err("rank_compounds", "No compounds could be scored (SMILES parsing failures).")

        scored.sort(key=lambda x: x["composite_score"], reverse=True)
        for i, r in enumerate(scored, 1):
            r["rank"] = i

        top = scored[:top_n]
        warnings = []
        if not ref_fps:
            warnings.append("No reference SMILES provided — LBVS Tanimoto component set to 0. Composite score reflects pIC50+QED+Lipinski only.")

        return self._ok(
            tool_name="rank_compounds",
            data={
                "total_scored": len(scored),
                "top_n": len(top),
                "composite_formula": "0.40×pIC50_norm + 0.30×LBVS_Tanimoto + 0.20×QED + 0.10×Lipinski",
                "top_compounds": top,
            },
            citations=[_RDKIT_CITATION],
            confidence=0.80,
            rationale="Composite ranking based on real measured IC50 values and RDKit-computed properties. "
                      "Weights adapted from literature (pIC50 dominant; LBVS provides scaffold similarity context).",
            warnings=warnings,
        )
