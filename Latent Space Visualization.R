# scripts/03_latent_analysis.R
library(ggplot2)
library(reticulate)

encoder <- keras::load_model_hdf5("scripts/vae_encoder.h5")
latent <- encoder$predict(expr_data)

latent_df <- as.data.frame(latent)
latent_df$cell <- rownames(seurat_obj@meta.data)

# Use Seurat metadata
latent_df$cluster <- seurat_obj$seurat_clusters

ggplot(latent_df, aes(x = V1, y = V2, color = cluster)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Latent Space of Cells (VAE)")
ggsave("figures/umap_latent.png")
