"""
MDR-TB AI Pipeline v3 — Test Suite (Phase 15)
===============================================
Basic tests for all pipeline modules.
"""

import sys
import os
import pytest
import numpy as np
import pandas as pd
from pathlib import Path

# Add project root
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from utils.helpers import (
    load_config, generate_mock_expression_matrix,
    generate_mock_who_data, generate_mock_compounds,
    generate_mock_drugbank_entries, ensure_dir
)


# ── UTILITY TESTS ─────────────────────────────────────────

class TestUtils:
    def test_load_config(self):
        config = load_config(str(PROJECT_ROOT / "config.yaml"))
        assert "project" in config
        assert config["project"]["version"] == "3.0.0"
        assert config["project"]["safety_mode"] == "computational_only"

    def test_mock_expression(self):
        df = generate_mock_expression_matrix(n_genes=100, n_samples=10)
        assert df.shape == (100, 10)
        assert df.index.name == "Gene"
        assert any(c.startswith("CTRL") for c in df.columns)
        assert any(c.startswith("MDR") for c in df.columns)

    def test_mock_who_data(self):
        df = generate_mock_who_data(regions=["A", "B"], years=[2020, 2021])
        assert len(df) == 4
        assert "Region" in df.columns
        assert "TB_Incidence_per_100k" in df.columns

    def test_mock_compounds(self):
        df = generate_mock_compounds(n=20)
        assert len(df) == 20
        assert "SMILES" in df.columns
        assert "Activity_Label" in df.columns

    def test_mock_drugbank(self):
        entries = generate_mock_drugbank_entries(n=10)
        assert len(entries) == 10
        assert "drugbank_id" in entries[0]

    def test_ensure_dir(self, tmp_path):
        d = ensure_dir(str(tmp_path / "test_sub"))
        assert d.exists()


# ── DATA CONNECTOR TESTS ──────────────────────────────────

class TestDataConnectors:
    def test_geo_loader_mock(self):
        from data_connectors.geo_loader import GEOLoader
        loader = GEOLoader()
        df = loader.load()
        assert isinstance(df, pd.DataFrame)
        assert df.shape[0] > 0
        assert df.shape[1] > 0

    def test_who_parser_mock(self):
        from data_connectors.who_parser import WHOParser
        parser = WHOParser()
        df = parser.load()
        assert isinstance(df, pd.DataFrame)
        assert "Region" in df.columns

    def test_drugbank_loader_mock(self):
        from data_connectors.drugbank_loader import DrugBankLoader
        loader = DrugBankLoader()
        df = loader.load()
        assert isinstance(df, pd.DataFrame)
        assert "drugbank_id" in df.columns

    def test_pubchem_loader_mock(self):
        from data_connectors.pubchem_loader import PubChemLoader
        loader = PubChemLoader()
        df = loader.load()
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0


# ── ENGINE TESTS ──────────────────────────────────────────

class TestOmicsEngine:
    def test_differential_expression(self):
        from engines.omics_engine import OmicsEngine
        engine = OmicsEngine()
        expr = generate_mock_expression_matrix(n_genes=200, n_samples=10)
        de = engine.differential_expression(expr)
        assert "Gene" in de.columns
        assert "log2FC" in de.columns
        assert "padj" in de.columns
        assert "significant" in de.columns
        assert len(de) == 200

    def test_pathway_scoring(self):
        from engines.omics_engine import OmicsEngine
        engine = OmicsEngine()
        expr = generate_mock_expression_matrix(n_genes=500, n_samples=10)
        de = engine.differential_expression(expr)
        pathways = engine.pathway_scoring(de, expr)
        assert "Pathway" in pathways.columns
        assert "Enrichment_Score" in pathways.columns

    def test_gene_ranking(self):
        from engines.omics_engine import OmicsEngine
        engine = OmicsEngine()
        expr = generate_mock_expression_matrix(n_genes=100, n_samples=10)
        de = engine.differential_expression(expr)
        ranked = engine.gene_ranking(de)
        assert "rank" in ranked.columns
        assert "composite_score" in ranked.columns
        assert ranked.iloc[0]["rank"] == 1


class TestEpiEngine:
    def test_compute_trends(self):
        from engines.epi_engine import EpidemiologyEngine
        engine = EpidemiologyEngine()
        who = generate_mock_who_data()
        trends = engine.compute_trends(who)
        assert "Region" in trends.columns
        assert "Slope" in trends.columns

    def test_mdr_patterns(self):
        from engines.epi_engine import EpidemiologyEngine
        engine = EpidemiologyEngine()
        who = generate_mock_who_data()
        mdr = engine.mdr_patterns(who)
        assert "Risk_Category" in mdr.columns

    def test_district_aggregation(self):
        from engines.epi_engine import EpidemiologyEngine
        engine = EpidemiologyEngine()
        who = generate_mock_who_data()
        dist = engine.district_aggregation(who)
        assert "Burden_Rank" in dist.columns


class TestTargetEngine:
    def test_score_targets(self):
        from engines.target_engine import TargetEngine
        engine = TargetEngine()
        result = engine.run()
        scored = result["scored_targets"]
        assert "Target" in scored.columns
        assert "Final_Score" in scored.columns
        assert "Priority" in scored.columns
        assert len(scored) > 0


class TestCompoundEngine:
    def test_compute_descriptors(self):
        from engines.compound_engine import CompoundEngine
        engine = CompoundEngine()
        compounds = generate_mock_compounds(n=20)
        desc = engine.compute_descriptors(compounds)
        assert "MolWt" in desc.columns

    def test_clustering(self):
        from engines.compound_engine import CompoundEngine
        engine = CompoundEngine()
        compounds = generate_mock_compounds(n=20)
        clustered = engine.cluster_compounds(compounds)
        assert "Cluster" in clustered.columns

    def test_lipinski(self):
        from engines.compound_engine import CompoundEngine
        engine = CompoundEngine()
        compounds = generate_mock_compounds(n=20)
        lip = engine.lipinski_filter(compounds)
        assert "Lipinski_Pass" in lip.columns
        assert "Lipinski_Violations" in lip.columns


# ── MODEL TESTS ───────────────────────────────────────────

class TestMLPipeline:
    def test_ml_training(self):
        from models.ml_pipeline import MLPipeline
        pipeline = MLPipeline()
        compounds = generate_mock_compounds(n=100)
        result = pipeline.run(compounds)
        assert "summary" in result
        summary = result["summary"]
        assert "best_algorithm" in summary
        assert summary.get("best_auc", 0) > 0


class TestResistanceModule:
    def test_map_genes(self):
        from models.resistance_module import ResistanceModule
        module = ResistanceModule()
        gene_map = module.map_mdr_genes()
        assert "Gene" in gene_map.columns
        assert len(gene_map) > 0

    def test_mutation_catalog(self):
        from models.resistance_module import ResistanceModule
        module = ResistanceModule()
        mutations = module.build_mutation_catalog()
        assert "Mutation" in mutations.columns
        assert len(mutations) > 0

    def test_resistance_scores(self):
        from models.resistance_module import ResistanceModule
        module = ResistanceModule()
        scores = module.compute_resistance_scores()
        assert "Resistance_Score" in scores.columns


class TestRankingEngine:
    def test_ranking(self):
        from models.ranking_engine import RankingEngine
        engine = RankingEngine()
        compounds = generate_mock_compounds(n=30)
        result = engine.run(compounds)
        ranked = result["ranked_compounds"]
        assert "Final_Rank_Score" in ranked.columns
        assert "Rank" in ranked.columns
        assert ranked.iloc[0]["Rank"] == 1


# ── GENERATOR TESTS ───────────────────────────────────────

class TestGenerators:
    def test_manuscript_generator(self):
        from generators.manuscript_generator import ManuscriptGenerator
        gen = ManuscriptGenerator()
        ms = gen.generate({})
        assert len(ms) > 100
        assert "Introduction" in ms
        assert "Methods" in ms
        assert "Results" in ms

    def test_icmr_generator(self):
        from generators.icmr_generator import ICMRGrantGenerator
        gen = ICMRGrantGenerator()
        prop = gen.generate({})
        assert len(prop) > 100
        assert "ICMR" in prop

    def test_iec_generator(self):
        from generators.iec_generator import IECProtocolGenerator
        gen = IECProtocolGenerator()
        prot = gen.generate({})
        assert len(prot) > 100
        assert "Ethics" in prot or "Ethical" in prot


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
