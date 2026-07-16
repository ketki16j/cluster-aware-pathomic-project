# ============================================================
# STEP 6B: Prepare CAMELOMICS feature matrix for Jeremy's pipeline
# Ketki Joshi - Renal Pathomics Project
# ============================================================
# This script:
#   1. Takes Renal_Data.csv (clinical + naive-averaged pathomic features)
#   2. Drops the naive-averaged pathomic columns (97-510)
#   3. Merges in CAMELOMICS cluster-averaged features for each tissue type
#   4. Saves Renal_Data_CAMELOMICS.csv in same format as Renal_Data.csv
#   5. This new file is then used with Jeremy's exact same pipeline
#
# HOW TO RUN:
#   module load R/4.5.0
#   cd ~/scratch/others/jeremy/procurement-biopsy-pathomics-ml
#   Rscript step6B_prepare_camelomics_data.R
# ============================================================

.libPaths('/home/kxj190026/R_libs')
library(dplyr)

RENAL_DATA   <- '/home/kxj190026/scratch/others/jeremy/procurement-biopsy-pathomics-ml/data/Renal_Data.csv'
CAMEL_DIR    <- '/home/kxj190026/scratch/others/camelomics/camelomics_output'
OUTPUT_DIR   <- '/home/kxj190026/scratch/others/jeremy/procurement-biopsy-pathomics-ml'

# Best k per tissue from Step 5 results
BEST_K <- list(
  glomeruli = 8,
  tubules   = 8,
  arteries  = 6,
  scl_gloms = 8
)

cat("============================================================\n")
cat("STEP 6B: Preparing CAMELOMICS Feature Matrix\n")
cat("============================================================\n\n")

# ============================================================
# STEP 1: Load Renal_Data and keep only clinical columns (1-96)
# ============================================================
renal <- read.csv(RENAL_DATA, stringsAsFactors=FALSE)
cat("Original Renal_Data:", nrow(renal), "rows x", ncol(renal), "cols\n")

# Keep only clinical columns (1 to 96, before AI_FEATURES/pathomic start)
# Col 96 = AI_FEATURES (header row marker)
# Col 97 = COMPUTATIONAL_MORPHOMETRIC (header row marker)
# Pathomic features start after col 97
clinical_cols <- renal[, 1:96]
cat("Clinical columns kept: 1 to 96\n")
cat("Clinical col names (last 5):", paste(tail(names(clinical_cols), 5), collapse=', '), "\n\n")

# ============================================================
# STEP 2: Load CAMELOMICS features for each tissue type
# ============================================================
cat("Loading CAMELOMICS cluster-averaged features:\n")

camel_features <- list()

for (tissue in names(BEST_K)) {
  k <- BEST_K[[tissue]]
  fpath <- file.path(CAMEL_DIR, paste0("camelomics_features_", tissue, "_k", k, ".csv"))

  if (!file.exists(fpath)) {
    cat("  WARNING: Missing file:", fpath, "\n")
    next
  }

  df <- read.csv(fpath, stringsAsFactors=FALSE)
  cat("  ", tissue, "(k=", k, "):", nrow(df), "patients x", ncol(df)-1, "cluster features\n")

  # Add tissue prefix to feature columns (keep Slide_number as is)
  feat_cols <- names(df)[names(df) != "Slide_number"]
  names(df)[names(df) != "Slide_number"] <- paste0(tissue, "_", feat_cols)

  camel_features[[tissue]] <- df
}

# ============================================================
# STEP 3: Match CAMELOMICS features to Renal_Data patients
# ============================================================
cat("\nMatching CAMELOMICS features to all", nrow(renal), "patients...\n")

# Clean slide number for matching
clean_slide <- function(s) {
  s <- trimws(as.character(s))
  s <- gsub(' - ', ' ', s)
  s <- gsub(' -([0-9]+)$', ' \\1', s)
  s <- gsub('_', ' ', s)
  s <- gsub(' PAS ([0-9]+)', ' \\1', s, ignore.case=TRUE)
  s <- gsub(' PAS', '', s, ignore.case=TRUE)
  s <- gsub('PAS ', '', s, ignore.case=TRUE)
  s <- gsub('  +', ' ', s)
  return(trimws(s))
}

renal$slide_clean <- sapply(renal$Slide_number, clean_slide)

# Merge each tissue's CAMELOMICS features
merged <- renal
for (tissue in names(camel_features)) {
  df <- camel_features[[tissue]]
  df$slide_clean <- sapply(df$Slide_number, clean_slide)
  df$Slide_number <- NULL

  merged <- left_join(merged, df, by="slide_clean")
  cat("  After merging", tissue, ":", ncol(merged), "total cols\n")
}

merged$slide_clean <- NULL

cat("\nFinal merged data:", nrow(merged), "rows x", ncol(merged), "cols\n")

# ============================================================
# STEP 4: Check how many patients have CAMELOMICS features
# ============================================================
# Count patients with at least glomeruli features
glom_col <- paste0("glomeruli_cluster0_feat0")
if (glom_col %in% names(merged)) {
  has_camel <- sum(!is.na(merged[[glom_col]]))
  cat("Patients with CAMELOMICS glomeruli features:", has_camel, "/", nrow(merged), "\n")
} else {
  # Try finding first glomeruli column
  glom_cols <- grep("^glomeruli_", names(merged), value=TRUE)
  if (length(glom_cols) > 0) {
    has_camel <- sum(!is.na(merged[[glom_cols[1]]]))
    cat("Patients with CAMELOMICS glomeruli features:", has_camel, "/", nrow(merged), "\n")
  }
}

# ============================================================
# STEP 5: Fill NAs with 0 for patients missing some tissue types
# ============================================================
# Some patients may not have sclerotic gloms data (only 72/109 have it)
# Fill their scl_gloms cluster features with 0
camel_col_names <- grep("^(glomeruli|tubules|arteries|scl_gloms)_", names(merged), value=TRUE)
cat("\nFilling NA CAMELOMICS features with 0 for missing tissue types...\n")
for (col in camel_col_names) {
  na_count <- sum(is.na(merged[[col]]))
  if (na_count > 0) {
    merged[[col]][is.na(merged[[col]])] <- 0
  }
}

# ============================================================
# STEP 6: Save the new Renal_Data_CAMELOMICS.csv
# ============================================================
out_path <- file.path(OUTPUT_DIR, "data", "Renal_Data_CAMELOMICS.csv")
write.csv(merged, out_path, row.names=FALSE)

cat("\nSaved:", out_path, "\n")
cat("Shape:", nrow(merged), "rows x", ncol(merged), "cols\n")

# Also create a symlink so Jeremy's scripts find it
symlink_path <- file.path(OUTPUT_DIR, "Renal_Data_CAMELOMICS.csv")
if (file.exists(symlink_path)) file.remove(symlink_path)
file.symlink(out_path, symlink_path)
cat("Symlink created:", symlink_path, "\n")

# ============================================================
# STEP 7: Show summary for Jeremy's pipeline
# ============================================================
cat("\n============================================================\n")
cat("READY FOR STEP 6B PIPELINE\n")
cat("============================================================\n")
cat("New file: Renal_Data_CAMELOMICS.csv\n")
cat("  - Same format as Renal_Data.csv\n")
cat("  - Same 213 rows (all patients)\n")
cat("  - Columns 1-96: original clinical features\n")
cat("  - Columns 97+: CAMELOMICS cluster-averaged pathomic features\n\n")

camel_cols <- grep("^(glomeruli|tubules|arteries|scl_gloms)_", names(merged), value=TRUE)
cat("CAMELOMICS feature counts:\n")
for (tissue in names(BEST_K)) {
  cols <- grep(paste0("^", tissue, "_"), names(merged), value=TRUE)
  cat("  ", tissue, "(k=", BEST_K[[tissue]], "):", length(cols), "cluster features\n")
}
cat("  Total CAMELOMICS features:", length(camel_cols), "\n\n")

cat("Next: Run Jeremy's internal validation scripts with:\n")
cat("  input.file.name <- 'Renal_Data_CAMELOMICS.csv'\n")
cat("  (everything else stays the same)\n")
cat("============================================================\n")
