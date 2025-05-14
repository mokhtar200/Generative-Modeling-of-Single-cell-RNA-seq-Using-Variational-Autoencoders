# scripts/02_train_vae.R
library(reticulate)
use_condaenv("r-reticulate")  # Ensure Keras is accessible
library(Seurat)

# Load data
seurat_obj <- readRDS("data/processed/seurat_pbmc.rds")
expr_data <- as.matrix(t(seurat_obj@assays$RNA@data[VariableFeatures(seurat_obj), ]))

# Scale to [0,1] for sigmoid decoder
expr_data <- scale(expr_data)
expr_data <- (expr_data - min(expr_data)) / (max(expr_data) - min(expr_data))

# Keras setup
keras <- import("keras")
k <- keras$backend

input_dim <- dim(expr_data)[2]
latent_dim <- 10L

# Encoder
inputs <- keras$layers$Input(shape = input_dim)
encoded <- keras$layers$Dense(128, activation = "relu")(inputs)
z_mean <- keras$layers$Dense(latent_dim)(encoded)
z_log_var <- keras$layers$Dense(latent_dim)(encoded)

sampling <- function(args) {
  z_mean <- args[[1]]
  z_log_var <- args[[2]]
  epsilon <- k$random_normal(shape = k$shape(z_mean), mean = 0.0, stddev = 1.0)
  z_mean + k$exp(0.5 * z_log_var) * epsilon
}

z <- keras$layers$Lambda(sampling)(list(z_mean, z_log_var))

# Decoder
decoder_input <- keras$layers$Input(shape = latent_dim)
decoded <- keras$layers$Dense(128, activation = "relu")(decoder_input)
outputs <- keras$layers$Dense(input_dim, activation = "sigmoid")(decoded)

# Models
encoder <- keras$Model(inputs, z)
decoder <- keras$Model(decoder_input, outputs)

vae_outputs <- decoder(z)

# Custom VAE model
vae <- keras$Model(inputs, vae_outputs)

# Loss
reconstruction_loss <- keras$losses$mean_squared_error(inputs, vae_outputs)
reconstruction_loss <- k$sum(reconstruction_loss)

kl_loss <- -0.5 * k$sum(1 + z_log_var - k$square(z_mean) - k$exp(z_log_var), axis = -1L)

vae_loss <- k$mean(reconstruction_loss + kl_loss)

vae$add_loss(vae_loss)
vae$compile(optimizer = "adam")

# Fit model
vae$fit(expr_data, epochs = 50, batch_size = 64, verbose = 1)

# Save encoder and decoder
encoder$save("scripts/vae_encoder.h5")
decoder$save("scripts/vae_decoder.h5")
