# 🧪 MSc Project Repository

This repository contains a set of C++/ROOT based tools to perform detector-driven physics studies with Monte-Carlo samples fully simulated with the International Large Detector (ILD) Geometry .

These tools were integrated into the **ILCSoft Marlin** framework and tested on simulated datasets using the NAF cluster @ Deutsches Elektronen-Synchrotron (DESY). The work aims to contribute to improved signal reconstruction strategies for future lepton collider experiments, such as the ILC.

📄 This work was carried out as part of my MSc thesis: 
**"Event Selection and Angular Reconstruction for W Boson Events at Future e+e− Colliders"**[Thesis](https://bib-pubdb1.desy.de/record/619472/)

## 📋 Features
A steering file with two main core tools for my MSc research project:
- **WWAngleCalculationProcessor**:
  Module for the reconstruction of on-shell and off-shell W bosons and their decay angles (**θ**, **φ**) for all 4f channels. 
- **WWCategorisationProcessor**:
  A classification tool that labels each 4f WW(W) event according to its decay channel topology (fully leptonic, fully hadronic, semi-leptonic), supporting systematic studies of event selection and efficiency in detector-level analyses.

## ⚙️ Supplementary Tools

- **Shell Scripts for HTCondor Automation**  
  Scripts for submitting and managing large-scale batch jobs using HTCondor. These jobs were used to generate and process additional MC samples with `WWCategorisationProcessor`, covering the full Standard Model background at √s = **250 GeV**.

- **ROOT Macros for Efficiency Optimization**  
  A set of analysis macros to visualize output, evaluate categorization performance and refine channel-based selections. These tools were used to further improve selection efficiency across different topologies.
