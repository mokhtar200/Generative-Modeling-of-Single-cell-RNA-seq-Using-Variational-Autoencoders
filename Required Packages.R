# requirements.R
install.packages(c("Seurat", "ggplot2", "patchwork"))
install.packages("reticulate")  # For connecting with Python Keras
reticulate::install_miniconda()
reticulate::py_install(packages = c("tensorflow", "keras", "scikit-learn"), pip = TRUE)
