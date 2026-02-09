# Fix renv installation
# Disable sandbox temporarily
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")

# Remove renv from loaded packages if present
if ("renv" %in% loadedNamespaces()) {
  try(unloadNamespace("renv"), silent = TRUE)
}

# Install renv to user library
cat("Installing renv to user library...\n")
install.packages("renv", repos = "https://cloud.r-project.org")

# Load renv
cat("Loading renv...\n")
library(renv)

# Restore project
cat("Restoring project packages...\n")
renv::restore()

cat("\nDone! You can now open your R project normally.\n")
