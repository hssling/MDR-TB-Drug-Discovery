# MDR-TB AI Drug Discovery Pipeline v5

An end-to-end, strictly computational artificial intelligence pipeline for Multi-Drug Resistant Tuberculosis (MDR-TB) drug discovery, epidemiology, and automated clinical research generation. 

Built with modular architectures spanning Transcriptomics (GEO), Epidemiological demographics (WHO), Structural Chemistry (PDB), Multi-Objective Machine Learning, and Automated Ethics/Grant reporting.

---

## 🧬 Key Innovative Capabilities

1. **Generative De Novo Drug Design**: Uses pharmacophore recombination to computationally generate totally novel molecular iterations (`models/de_novo_generator.py`).
2. **ESMFold Structural Prediction**: Automatically fetches protein 3D structures via Meta's ESMFold API for mutated targets lacking X-ray crystallography (`models/structural_predictor.py`).
3. **Advanced QSAR & LLMs**: Augments traditional molecular features (RDKit) with Deep NLP Transformers embeddings (`models/ml_pipeline.py`).
4. **Automated Documentation Generation**: Programmatically synthesizes IMRAD Manuscripts, ICMR Grant Proposals, and Institutional Ethical Committee Exemption waivers directly from runtime pipeline logs.

## 🚀 Quick Start (Local Run)

### Using Docker (Recommended)
You can launch the entire stack (Pipeline Execution + Dashboards) natively isolated.
```bash
docker-compose up --build
```
This will expose:
- `http://localhost:8501` - Epidemiology & Geo-Spatial Dashboard
- `http://localhost:8502` - Drug Discovery & Molecular Dashboard

### Local Python Execution
```bash
pip install -r requirements.txt
python run_pipeline.py
```
Outputs will gracefully generate into the `outputs/` directory.

## 🧪 CI/CD
This repository is configured with a **GitHub Actions Pipeline** (`.github/workflows/ci.yml`). Every commit automatically tests the 28 integrated Python sub-routines and executes a dry-run simulating the full biological discovery path to ensure no regressions occur.

## 🌐 Deployment (Hugging Face Spaces)
To expose your dashboards globally:
1. Create a new "Streamlit Space" on [Hugging Face](https://huggingface.co/spaces).
2. Push this repository to the HF remote. (Make sure you point `app.py` or the specific dashboard directly).
3. Example HF Sync:
```bash
git remote add space https://huggingface.co/spaces/[your-username]/mdrtb-pipeline
git push --force space main
```

## 📊 Deployment (Kaggle Integration)
To open this architecture to the global data science community:
1. Zip this entire folder and upload it as a "Kaggle Dataset".
2. Create a new Kaggle Environment Notebook and map the dataset.
3. Because the pipeline utilizes generic Python (`pandas`, `scikit-learn`, `transformers`), you can trigger the execution block straight from a single Kaggle Cell:
```python
!cp -r /kaggle/input/mdrtb-ai-pipeline-v5/* ./
!pip install -r requirements.txt
!python run_pipeline.py
```

## ⚠️ Safety Constraint Statement
This architecture is designated strictly for **computational discovery support only**. There are absolutely no wet-lab protocols, chemical synthesis pathways, or operational biological instructions anywhere within this repository. All conclusions must undergo experimental replication.
