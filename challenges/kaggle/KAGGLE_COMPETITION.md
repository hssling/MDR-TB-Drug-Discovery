# Kaggle Competition: MDR-TB Binding Affinity Prediction Challenge

**Overview**
Welcome to the MDR-TB Binding Affinity Prediction Challenge! 
Multi-Drug Resistant Tuberculosis (MDR-TB) is a global crisis. The enoyl-ACP reductase enzyme (**InhA**, UniProt: P9WGR1) is a primary target. In this competition, your objective is to predict the theoretical binding energy (kcal/mol) of novel chemical compounds against the AlphaFold-predicted structure of InhA WITHOUT running expensive AutoDock/PyRx physical simulations.

**The Dataset**
The dataset provided (`outputs/docking/docking_results_InhA.csv`) contains:
- `SMILES`: The semantic chemical representation.
- `Binding_Affinity_kcal_mol`: The target variable to predict (continuous regression).

**Evaluation Metric**
The competition is evaluated internally on **Root Mean Squared Error (RMSE)**. The goal is to build a Machine Learning model (e.g., XGBoost, ChemBERTa + Regression, GNN) that achieves an RMSE closest to 0.

**How to Enter**
1. Read the `starter_script.py` to learn how to parse the pipeline dataset and engineer features using RDKit.
2. Train your regression model on 80% of the dataset.
3. Predict the binding affinity of the testing 20%.
4. Publish your Notebook to this dataset's Kaggle page!

*Hosted by: The MDR-TB AI Pipeline v6 Initiative*
