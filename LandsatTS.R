devtools::install_github("logan-berner/LandsatTS", build_vignettes = TRUE)
library(rgee)
# Get the username
HOME <- Sys.getenv("HOME")

# 1. Install miniconda
reticulate::install_miniconda()

# 2. Install Google Cloud SDK
system("curl -sSL https://sdk.cloud.google.com | bash")

# 3 Set global parameters
Sys.setenv("RETICULATE_PYTHON" = sprintf("%s/.local/share/r-miniconda/bin/python3", HOME))
Sys.setenv("EARTHENGINE_GCLOUD" = sprintf("%s/google-cloud-sdk/bin/", HOME))

# 4 Instal thrid-party software
ee_install()
## IMPORTANT!! You have to <<TERMINATE R>>


# If you are a RSTUDIO user! -------------------------
# Use the <TERMINAL> to authenticate earthengine
# (why? check this -> https://twitter.com/csaybar/status/1576423665792479233)
display_terminal <- function() {
  ms1 <- "RUN THIS IN THE TERMINAL!!!\n"
  ms2 <- paste0("csaybar@csaybar-pc01:", sprintf("PATH=$PATH:%s/google-cloud-sdk/bin/", HOME), "\n")
  ms3 <- paste0("csaybar@csaybar-pc01:", "earthengine authenticate --quiet", "\n")
  cat(ms1, ms2, ms3)
}
display_terminal()

# 6. Run ee_Initialize
ee_Initialize()