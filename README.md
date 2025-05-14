# Generative Modeling of Single-cell RNA-seq Using Variational Autoencoders

This project uses **Variational Autoencoders (VAEs)** in R to explore and simulate single-cell gene expression data using the PBMC3k dataset.

## Features
- scRNA-seq preprocessing with Seurat
- Deep learning with VAE using Keras in R
- Latent space visualization and trajectory exploration
- Simulation of new cells

## How to Run
1. Install dependencies from `requirements.R`
2. Run scripts in order:
   - `01_preprocessing.R`
   - `02_train_vae.R`
   - `03_latent_analysis.R`
   - `04_generate_cells.R`

## Dataset
- [10x Genomics PBMC 3k](https://support.10xgenomics.com/single-cell-gene-expression/datasets)

## License
MIT License
