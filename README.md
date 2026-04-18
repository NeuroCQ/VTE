"GLOBE (GLiOma platform for individualized prediction of postoperative Binary Events) is an explainable machine learning platform for preoperative prediction of postoperative venous thromboembolism (VTE) and neurological deterioration (PND) in patients undergoing glioma surgery.

🌐 Web application: https://gliomas.shinyapps.io/GLOBE/


Overview
This repository contains the complete analysis code for the study:

"GLOBE: An Explainable Machine Learning Platform for Preoperative Prediction of Thromboembolism and Neurological Deficits in Glioma"
Feiling Xiang†, Xuelian Yang†, Sijin Xiang, Mengyuan Fu, Gang Yang*
The First Affiliated Hospital of Chongqing Medical University

The repository includes:

Feature selection pipelines for VTE and PND prediction
Training and evaluation of four machine learning algorithms (GLM, RDF, GBA, SVM)
SHAP-based model interpretation
Decision curve analysis and calibration assessment
An example input data template

All models use only preoperative variables to ensure temporal validity and avoid information leakage.
