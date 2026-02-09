# Fix renv::restore() issues
# This script addresses connection and repository problems

cat("Configuring renv settings...\n")

# 1. Use standard CRAN instead of RSPM to avoid connection issues
options(repos = c(
  CRAN = "https://cloud.r-project.org",
  BioCsoft = "https://bioconductor.org/packages/3.20/bioc",
  BioCann = "https://bioconductor.org/packages/3.20/data/annotation",
  BioCexp = "https://bioconductor.org/packages/3.20/data/experiment"
))

# 2. Disable problematic download methods and use more reliable ones
Sys.setenv(RENV_DOWNLOAD_METHOD = "wininet")  # Better for Windows
options(download.file.method = "wininet")

# 3. Skip GitHub packages temporarily (we'll install them separately if needed)
Sys.setenv(RENV_CONFIG_GITHUB_ENABLED = "FALSE")

# 4. Increase timeout for slow downloads
options(timeout = 300)

# 5. Disable sandbox
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")

cat("\nAttempting to restore packages with fixed settings...\n")
cat("This may take several minutes...\n\n")

# Load renv
library(renv)

# Try to restore
tryCatch({
  renv::restore(prompt = FALSE)
  cat("\n✓ Package restoration successful!\n")
}, error = function(e) {
  cat("\n✗ Some packages failed. Error:\n")
  cat(conditionMessage(e), "\n")
  cat("\nTrying to install critical packages manually...\n")
  
  # Install msigdbr manually from CRAN
  tryCatch({
    install.packages("msigdbr", repos = "https://cloud.r-project.org")
    cat("✓ msigdbr installed successfully\n")
  }, error = function(e2) {
    cat("✗ msigdbr installation failed\n")
  })
})

cat("\nDone!\n")
cat("If you still have issues, you may need to:\n")
cat("1. Check your internet connection\n")
cat("2. Check if you're behind a corporate proxy\n")
cat("3. Install problematic packages manually\n")
