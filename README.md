# 🌍 SESAME Explorer

🔗 **Live app:** https://sesame.streamlit.app/  
🔗 **SESAME project:** https://github.com/A2Faisal/SESAME

---

## Overview

**SESAME Explorer** is an interactive web application for exploring and visualizing datasets from the **SESAME (Surface Earth System Analysis and Modeling Environment)** project.

The app provides a lightweight, browser-based interface to:
- Browse SESAME NetCDF datasets by Earth system sphere
- Inspect variables and metadata
- Slice and aggregate data across time and depth
- Generate publication-ready global maps
- Download figures directly for research and teaching use

The explorer is designed to support **open science**, **reproducibility**, and **research-grade data exploration** without requiring users to write code.

---

## Features

- 📂 Sphere-based dataset browsing (Atmosphere, Hydrosphere, Lithosphere, Biosphere, Cryosphere)
- 🧭 Interactive variable selection with metadata display
- ⏱ Time slicing and multi-time aggregation (mean, sum, min, max, median, std)
- 🌊 Depth slicing and multi-depth aggregation
- 🗺 Global map visualization using SESAME plotting utilities
- 🎨 Customizable colormaps and value ranges
- 💾 One-click export of figures (PNG, PDF, SVG)
- ⚡ On-demand data download with local caching

---

## Data Source

All datasets visualized in this app are part of the **SESAME Human–Earth Atlas**, described in:

> Faisal, A. A., Kaye, M., Ahmed, M., & Galbraith, E. (2025).  
> *The SESAME Human–Earth Atlas*.  
> **Scientific Data**, 12, 775.  
> https://doi.org/10.1038/s41597-025-05087-5

The atlas is also available via figshare:  
https://doi.org/10.6084/m9.figshare.28432499

---

## Relationship to the SESAME Project

This repository (**sesame-explorer**) contains **only the interactive explorer application**.

The core SESAME project, including:
- data generation pipelines
- analysis tools
- plotting utilities
- scientific workflows

is developed and maintained in the main SESAME repository:

👉 **https://github.com/A2Faisal/SESAME**

The explorer uses SESAME as a dependency (`import sesame as ssm`) and does not duplicate core functionality.

---

## Running Locally

```bash
git clone https://github.com/A2Faisal/sesame-explorer.git
cd sesame-explorer
pip install -r requirements.txt
streamlit run app.py
