# ===========================
# MR Pipeline (Improved)
# - Restrict to intersected genes
# - Harmonise
# - Steiger directionality (filter instruments supporting exposure -> outcome)
# - MR (IVW + sensitivity)
# - Save OR / heterogeneity / pleiotropy / Steiger summary
# ===========================

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(gwasglue)
  library(TwoSampleMR)
  library(RadialMR)
  library(dplyr)
  library(tools)
})

# ---------------------------
# 0) USER SETTINGS (EDIT ME)
# ---------------------------
workDir <- "D:\\生信分析狮\\1双疾病（结肠癌、银屑疾）\\14.第一个疾病第2个数据3种细胞\\切割分第3份"
setwd(workDir)

# Input files
genes_file <- "IntersectionGenes_Combined.csv"   # first column: gene names (must match eqtl.txt$exposure)
eqtl_file  <- "eqtl.txt"                         # exposure eQTL table

# Outcome VCF (auto-pick first .vcf.gz in folder)
outcomeFiles <- list.files(pattern = "\\.vcf\\.gz$")
if (length(outcomeFiles) < 1) stop("No .vcf.gz outcome file found in the working directory.")
outcomeFile <- outcomeFiles[1]
diseaseName <- file_path_sans_ext(file_path_sans_ext(outcomeFile)) # remove .vcf.gz

# Output folder (optional)
resultDir <- "results_output1"
if (!dir.exists(resultDir)) dir.create(resultDir)

# Steiger needs case/control info for binary outcomes (strongly recommended)
# >>> Fill these using FinnGen / GWAS Catalog <<<
NCASE_OUT    <- 1396  # e.g., 12345
NCONTROL_OUT <- 174006   # e.g., 67890

# If you don't know prevalence, use sample prevalence as approximation
PREV_OUT <- if (!is.na(NCASE_OUT) && !is.na(NCONTROL_OUT) && (NCASE_OUT + NCONTROL_OUT) > 0) {
  NCASE_OUT / (NCASE_OUT + NCONTROL_OUT)
} else {
  NA_real_
}

# Filtering and QC parameters
OUTCOME_PVAL_FILTER <- 5e-06   # keep SNPs with pval.outcome > 5e-06 (your original rule)
MIN_SNP_FOR_MR      <- 2       # need >=2 SNPs to run MR
RUN_STEIGER         <- TRUE    # set FALSE if you want to skip Steiger
STEIGER_STRICT      <- TRUE    # TRUE: require steiger_dir==TRUE; FALSE: keep all if steiger fails
SAVE_PER_GENE_FILES <- FALSE   # TRUE: save per-gene harmonised datasets (bigger output)

message("========== MR Analysis Pipeline (Improved) ==========")
message(sprintf("Working directory: %s", getwd()))
message(sprintf("Outcome file: %s", outcomeFile))
message(sprintf("Disease name: %s", diseaseName))
message(sprintf("RUN_STEIGER: %s", RUN_STEIGER))

# ---------------------------
# 1) Load intersected genes
# ---------------------------
genes_data <- read.csv(genes_file, stringsAsFactors = FALSE)
target_genes <- genes_data[, 1] |> as.character() |> unique()
message(sprintf("[Step 1] Target genes: %d", length(target_genes)))
message(sprintf("Example genes: %s", paste(head(target_genes, 5), collapse = ", ")))

# ---------------------------
# 2) Load and filter eQTL
# ---------------------------
message("[Step 2] Reading eQTL data...")
eqtl_raw <- read.table(eqtl_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")

# Filter to target genes (intersection genes)
eqtl_filtered <- eqtl_raw[eqtl_raw$exposure %in% target_genes, ]
message(sprintf("eQTL rows after gene filtering: %d (genes matched: %d)",
                nrow(eqtl_filtered), length(unique(eqtl_filtered$exposure))))

unmatched_genes <- setdiff(target_genes, unique(eqtl_raw$exposure))
if (length(unmatched_genes) > 0) {
  message(sprintf("Warning: %d target genes not found in eQTL file (show first 10): %s",
                  length(unmatched_genes), paste(head(unmatched_genes, 10), collapse = ", ")))
}

# Save filtered exposure file
exposureFile <- file.path(resultDir, "eqtl_filtered.txt")
write.table(eqtl_filtered, file = exposureFile, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
message(sprintf("Saved filtered exposure: %s", exposureFile))

# Simple per-gene summary
eqtl_gene_summary <- eqtl_filtered %>%
  group_by(exposure) %>%
  summarise(
    SNP_count = n(),
    Mean_beta = mean(beta, na.rm = TRUE),
    Mean_se   = mean(se, na.rm = TRUE),
    Min_pval  = min(pval, na.rm = TRUE),
    .groups   = "drop"
  ) %>% arrange(desc(SNP_count))

write.csv(eqtl_gene_summary, file = file.path(resultDir, "eQTL_Gene_Summary.csv"), row.names = FALSE)

# ---------------------------
# 3) Read exposure data (TwoSampleMR format)
# ---------------------------
message("[Step 3] Reading exposure (TwoSampleMR format)...")
expData <- read_exposure_data(
  filename = exposureFile,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  pval_col = "pval",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  phenotype_col = "exposure",
  samplesize_col = "samplesize",  # ok if exists; if not in file, TwoSampleMR will set NA
  chr_col = "chr",
  pos_col = "pos",
  clump = FALSE
)

# If samplesize.exposure is missing/NA, set to eQTLGen max (optional)
if (!"samplesize.exposure" %in% names(expData) || all(is.na(expData$samplesize.exposure))) {
  expData$samplesize.exposure <- 31684
}

message(sprintf("Exposure loaded: %d SNPs, %d genes",
                nrow(expData), length(unique(expData$exposure))))

# ---------------------------
# 4) Read outcome VCF and convert
# ---------------------------
message("[Step 4] Reading outcome VCF and converting...")
vcf_obj <- readVcf(outcomeFile)
outcomeData <- gwasvcf_to_TwoSampleMR(vcf = vcf_obj, type = "outcome")
message(sprintf("Outcome loaded: %d SNPs", nrow(outcomeData)))

# Merge to construct a local outcome instruments file (your original style)
message("[Step 5] Merging exposure SNPs with outcome to build outcome_instruments.csv ...")
merged_raw <- merge(expData, outcomeData, by = "SNP")
message(sprintf("Merged exposure-outcome (raw): %d SNPs", nrow(merged_raw)))

outcome_instruments_file <- file.path(resultDir, "outcome_instruments.csv")
# Keep only outcome columns (exclude exposure columns)
write.csv(merged_raw[, -(2:ncol(expData))], file = outcome_instruments_file, row.names = FALSE)

outData <- read_outcome_data(
  snps = expData$SNP,
  filename = outcome_instruments_file,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "pval.outcome",
  eaf_col = "eaf.outcome"
)

# Add outcome meta needed for Steiger (binary trait)
outData$outcome <- diseaseName
outData$units.outcome <- "log odds"
outData$ncase.outcome <- NCASE_OUT
outData$ncontrol.outcome <- NCONTROL_OUT
outData$prevalence.outcome <- PREV_OUT

message(sprintf("Outcome prepared for MR: %d SNPs", nrow(outData)))

# ---------------------------
# 6) MR loop per gene
# ---------------------------
message("[Step 6] Running MR per gene...")

uniqueExp <- unique(expData$exposure)
numExp <- length(uniqueExp)

all_odds_ratios <- data.frame()
all_heterogeneity_results <- data.frame()
all_pleiotropy_results <- data.frame()
all_steiger_results <- data.frame()

clean_filename <- function(x) {
  x <- gsub("[^[:alnum:]_\\-\\.]", "_", x)
  x <- gsub("_+", "_", x)
  x <- trimws(x)
  x
}

pb <- txtProgressBar(min = 0, max = numExp, style = 3)

for (i in seq_along(uniqueExp)) {
  currentID <- uniqueExp[i]
  currentID_clean <- clean_filename(currentID)
  
  currentSubset <- expData[expData$exposure == currentID, ]
  if (nrow(currentSubset) == 0) {
    setTxtProgressBar(pb, i); next
  }
  
  # Harmonise
  harm <- tryCatch({
    harmonise_data(currentSubset, outData)
  }, error = function(e) {
    warning(sprintf("Harmonise failed for %s: %s", currentID, e$message))
    return(NULL)
  })
  if (is.null(harm) || nrow(harm) == 0) {
    setTxtProgressBar(pb, i); next
  }
  
  # TwoSampleMR keeps only good harmonised variants with mr_keep==TRUE automatically in downstream
  harm <- harm[harm$mr_keep == TRUE, ]
  if (nrow(harm) < MIN_SNP_FOR_MR) {
    setTxtProgressBar(pb, i); next
  }
  
  # Optional: outcome p-value filter (your original rule)
  harm2 <- harm[harm$pval.outcome > OUTCOME_PVAL_FILTER, ]
  if (nrow(harm2) < MIN_SNP_FOR_MR) {
    setTxtProgressBar(pb, i); next
  }
  
  # ---- Steiger directionality ----
  steiger_ok <- NA
  steiger_prop_true <- NA_real_
  steiger_min_p <- NA_real_
  nsnp_before <- nrow(harm2)
  nsnp_after <- nsnp_before
  
  harm_for_mr <- harm2
  
  if (RUN_STEIGER) {
    # Steiger needs ncase/ncontrol for binary outcome; if missing, it may fail
    st <- tryCatch({
      steiger_filtering(harm2)
    }, error = function(e) {
      warning(sprintf("Steiger failed for %s: %s", currentID, e$message))
      return(NULL)
    })
    
    if (!is.null(st) && nrow(st) > 0 && "steiger_dir" %in% names(st)) {
      steiger_prop_true <- mean(st$steiger_dir, na.rm = TRUE)
      steiger_min_p <- suppressWarnings(min(st$steiger_pval, na.rm = TRUE))
      # Filter to instruments supporting exposure -> outcome
      st2 <- st[st$steiger_dir == TRUE, ]
      nsnp_after <- nrow(st2)
      steiger_ok <- TRUE
      
      if (STEIGER_STRICT) {
        if (nrow(st2) < MIN_SNP_FOR_MR) {
          # Not enough SNPs after Steiger => skip this gene
          all_steiger_results <- rbind(
            all_steiger_results,
            data.frame(exposure = currentID, outcome = diseaseName,
                       nsnp_before = nsnp_before, nsnp_after = nsnp_after,
                       prop_steiger_dir_true = steiger_prop_true,
                       min_steiger_pval = steiger_min_p,
                       steiger_applied = TRUE, kept_for_mr = FALSE)
          )
          setTxtProgressBar(pb, i); next
        }
        harm_for_mr <- st2
      } else {
        # Non-strict: use st2 if possible, otherwise fall back to harm2
        harm_for_mr <- if (nrow(st2) >= MIN_SNP_FOR_MR) st2 else harm2
      }
    } else {
      steiger_ok <- FALSE
      # If strict and Steiger fails, skip; else continue with harm2
      if (STEIGER_STRICT) {
        all_steiger_results <- rbind(
          all_steiger_results,
          data.frame(exposure = currentID, outcome = diseaseName,
                     nsnp_before = nsnp_before, nsnp_after = NA_integer_,
                     prop_steiger_dir_true = NA_real_,
                     min_steiger_pval = NA_real_,
                     steiger_applied = FALSE, kept_for_mr = FALSE)
        )
        setTxtProgressBar(pb, i); next
      }
    }
  }
  
  # Save Steiger summary (even if not applied)
  all_steiger_results <- rbind(
    all_steiger_results,
    data.frame(exposure = currentID, outcome = diseaseName,
               nsnp_before = nsnp_before, nsnp_after = nsnp_after,
               prop_steiger_dir_true = steiger_prop_true,
               min_steiger_pval = steiger_min_p,
               steiger_applied = RUN_STEIGER, kept_for_mr = TRUE)
  )
  
  # Optionally save per-gene harmonised dataset
  if (SAVE_PER_GENE_FILES) {
    write.csv(harm_for_mr, file = file.path(resultDir, paste0("harmonised_", currentID_clean, ".csv")), row.names = FALSE)
  }
  
  # MR main
  mrResult <- tryCatch({
    mr(harm_for_mr)
  }, error = function(e) {
    warning(sprintf("MR failed for %s: %s", currentID, e$message))
    return(NULL)
  })
  if (!is.null(mrResult) && nrow(mrResult) > 0) {
    orResult <- generate_odds_ratios(mrResult)
    all_odds_ratios <- rbind(all_odds_ratios, orResult)
  }
  
  # Heterogeneity
  heteroResult <- tryCatch({
    mr_heterogeneity(harm_for_mr)
  }, error = function(e) {
    warning(sprintf("Heterogeneity failed for %s: %s", currentID, e$message))
    return(NULL)
  })
  if (!is.null(heteroResult) && nrow(heteroResult) > 0) {
    all_heterogeneity_results <- rbind(all_heterogeneity_results, heteroResult)
  }
  
  # Pleiotropy (Egger intercept)
  pleioResult <- tryCatch({
    mr_pleiotropy_test(harm_for_mr)
  }, error = function(e) {
    warning(sprintf("Pleiotropy test failed for %s: %s", currentID, e$message))
    return(NULL)
  })
  if (!is.null(pleioResult) && nrow(pleioResult) > 0) {
    all_pleiotropy_results <- rbind(all_pleiotropy_results, pleioResult)
  }
  
  setTxtProgressBar(pb, i)
}

close(pb)

# ---------------------------
# 7) Save outputs
# ---------------------------
write.csv(all_odds_ratios, file = file.path(resultDir, "all_odds_ratios.csv"), row.names = FALSE)
write.csv(all_heterogeneity_results, file = file.path(resultDir, "all_heterogeneity_results.csv"), row.names = FALSE)
write.csv(all_pleiotropy_results, file = file.path(resultDir, "all_pleiotropy_results.csv"), row.names = FALSE)
write.csv(all_steiger_results, file = file.path(resultDir, "all_steiger_results.csv"), row.names = FALSE)

cat("\n========== DONE ==========\n")
cat(sprintf("Output folder: %s\n", file.path(getwd(), resultDir)))
cat("Files:\n")
cat("  1) all_odds_ratios.csv\n")
cat("  2) all_heterogeneity_results.csv\n")
cat("  3) all_pleiotropy_results.csv\n")
cat("  4) all_steiger_results.csv\n")
cat("  5) eQTL_Gene_Summary.csv\n")
cat("  6) outcome_instruments.csv\n")
