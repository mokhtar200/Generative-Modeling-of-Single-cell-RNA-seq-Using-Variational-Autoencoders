# scripts/04_generate_cells.R
decoder <- keras::load_model_hdf5("scripts/vae_decoder.h5")

# Sample latent vectors
n <- 100
latent_new <- matrix(rnorm(n * latent_dim), nrow = n)
simulated_cells <- decoder$predict(latent_new)

# Save as CSV
write.csv(simulated_cells, "data/processed/simulated_cells.csv", row.names = FALSE
