

# ============================================================
# Genome Simulation with Artificial Motif Injection
# ============================================================
# This script simulates a synthetic genome, ensures no variants
# of defined motifs appear in the background, and then injects
# controlled motif instances into a fraction of the genes.
#
# Workflow:
#   1) Load helper functions (from functions.R)
#   2) Define parameters (edit here only!)
#   3) Compute background model from a reference FASTA
#   4) Generate synthetic genome (without motifs/variants)
#   5) Inject motifs into subset of genes
#   6) Simulate a biological vector (bio_vec)
#   7) Export outputs: FASTA and CSV files
#
# NOTE: Typically, you only need to modify Section (2).
# ============================================================


# -------------------
# 1) Load libraries and helper functions
# -------------------

library(dplyr)

# functions.R must be in the same directory:
source("src/functions.R")


# -------------------
# 2) Parameters (USER-DEFINED)
# -------------------

# === Input ===
# Input reference FASTA file:
# Used to estimate nucleotide background (base frequencies)
# and conditional probabilities (first-order Markov chain).
fasta_path  <- "data/cluster_all_filtered_longest.fa"

# === Genome size ===
# Number of synthetic genes to simulate.
n_genes     <- 7000   

# Length of each synthetic gene (in nucleotides).
gene_length <- 800   

# === Motif definition ===
# Define motifs as Position Probability Matrices (PPMs).
# Each row = position in motif, each column = probability of A/C/G/T.
# Each row MUST sum to 1.
# ⚠️ You can define MORE THAN ONE motif:
#    - Simply add more entries to the ppms list.
#    - Each motif must have a unique name (e.g., "motif_A", "motif_B").
#    - Example:
#        ppms <- list(
#          motif_A = matrix(...),
#          motif_B = matrix(...)
#        )

bases <- c("A","C","G","T")

ppms <- list(
  positive_1 = matrix(
    c(0.05,0.05,0.85,0.05,   # pos1: G-rich
      0.05,0.85,0.05,0.05,   # pos2: C-rich
      0.85,0.05,0.05,0.05,   # pos3: A-rich
      0.05,0.85,0.05,0.05,   # pos4: C-rich
      0.05,0.05,0.05,0.85,   # pos5: T-rich
      0.05,0.05,0.05,0.85),  # pos6: T-rich
    ncol = 4, byrow = TRUE,
    dimnames = list(NULL, bases)
  ),
  
  positive_2 = matrix(
    c(0.05,0.05,0.05,0.85,   # pos1: T-rich
      0.85,0.05,0.05,0.05,   # pos2: A-rich
      0.05,0.05,0.05,0.85,   # pos3: T-rich
      0.05,0.05,0.05,0.85,   # pos4: T-rich
      0.05,0.05,0.05,0.85,   # pos5: T-rich
      0.85,0.05,0.05,0.05),  # pos6: A-rich
    ncol = 4, byrow = TRUE,
    dimnames = list(NULL, bases)
  ),
  
  negative_1 = matrix(
    c(0.05,0.05,0.05,0.85,   # pos1: T-rich
      0.05,0.05,0.05,0.85,   # pos2: T-rich
      0.05,0.05,0.05,0.85,   # pos3: T-rich
      0.05,0.05,0.05,0.85,   # pos4: T-rich
      0.05,0.05,0.05,0.85,   # pos5: T-rich
      0.05,0.05,0.05,0.85),  # pos6: T-rich
    ncol = 4, byrow = TRUE,
    dimnames = list(NULL, bases)
  ),
  
  negative_2 = matrix(
    c(0.05,0.05,0.05,0.85,   # pos1: T-rich
      0.05,0.05,0.05,0.85,   # pos2: T-rich
      0.05,0.85,0.05,0.05,   # pos3: C-rich
      0.05,0.05,0.85,0.05,   # pos4: G-rich
      0.05,0.85,0.05,0.05),  # pos5: C-rich
    ncol = 4, byrow = TRUE,
    dimnames = list(NULL, bases)
  ),
  
  negative_3 = matrix(
    c(0.85,0.05,0.05,0.05,   # pos1: A-rich
      0.85,0.05,0.05,0.05,   # pos2: A-rich
      0.05,0.05,0.05,0.85,   # pos3: T-rich
      0.05,0.05,0.05,0.85,   # pos4: T-rich
      0.05,0.05,0.85,0.05),  # pos5: G-rich
    ncol = 4, byrow = TRUE,
    dimnames = list(NULL, bases)
  ),
  
  negative_4 = matrix(
    c(0.05,0.05,0.85,0.05,   # pos1: G-rich
      0.05,0.05,0.85,0.05,   # pos2: G-rich
      0.85,0.05,0.05,0.05,   # pos3: A-rich
      0.85,0.05,0.05,0.05,   # pos4: A-rich
      0.05,0.05,0.85,0.05,   # pos5: G-rich
      0.05,0.05,0.85,0.05),  # pos6: G-rich
    ncol = 4, byrow = TRUE,
    dimnames = list(NULL, bases)
  )
  
)

# Maximum Hamming distance allowed for motif variants.
# All variants within this distance from the motif consensus will be
# FORBIDDEN in the background genome.
variant_distances <- c(positive_1 = 1, positive_2 = 1,
                       negative_1 = 1, negative_2 = 1, negative_3 = 1, negative_4 = 1)

# === Fractions of injection ===
# Fractions of genes into which the motif will be injected.
# Example: fractions <- c(0.10, 0.30, 0.50)
# → Three separate simulations will be run:
#    * inject motif into 10% of genes
#    * inject motif into 30% of genes
#    * inject motif into 50% of genes
#
# For each fraction, SEPARATE output files will be created:
#   outputs/simulated_genome_frac_0.30.fasta
#   outputs/bio_vec_frac_0.30.csv
#   outputs/base_freqs_positive_1_frac_0.30.csv
fractions <- c(0.05)

# === Biological vector (bio_vec) ===
# bio_vec represents a continuous biological trait (e.g. expression).
# All genes get values from a log-normal distribution.
# Genes with injected motifs get shifted upward depending on motif score.
#
# meanlog, sdlog: baseline log-normal distribution parameters.
# multipliers: scaling factor for each motif – determines how strongly
# motif presence affects bio_vec.
#
# Formula for a motif-injected gene:
#   effective_meanlog = meanlog + multiplier * (score / threshold) * sdlog
#
# - score = log-odds of the injected motif (higher = stronger match)
# - threshold = maximum possible log-odds (S_max, "perfect match")
# - score/threshold ∈ [0,1]
# - multiplier controls how strongly motif scores shift bio_vec
#
# Example:
#   multiplier = 3, sdlog = 0.25
#   If score = threshold (perfect match), shift = 3 * 1 * 0.25 = 0.75
#   So meanlog increases from 3 → 3.75, leading to larger bio_vec values.
# ⚠️ Note: The sign of 'multiplier' determines the direction of the effect.

meanlog     <- 3     
sdlog       <- 0.25  
multipliers <- list(positive_1 = 3, positive_2 = 3,
                     negative_1 = -3, negative_2 = -3, negative_3 = -3, negative_4 = -3)

# === Output directory ===
out_dir <- "data/simulation_outputs"
dir.create(out_dir, showWarnings = FALSE)


# ============================================================
# 3) Step A: Background model from FASTA
# ============================================================
cat("Step A: Computing background model from FASTA...\n")

fasta_lines    <- readLines(fasta_path)
sequence_lines <- fasta_lines[!grepl("^>", fasta_lines)]
full_sequence  <- paste0(sequence_lines, collapse = "")
dna_vector     <- toupper(unlist(strsplit(full_sequence, "")))
dna_vector     <- dna_vector[dna_vector %in% c("A","C","G","T")]

cond_probs       <- compute_conditional_probs(dna_vector)
background_probs <- compute_background_freq(dna_vector)

cat("Background computed. Base frequencies:\n")
print(background_probs)


# ============================================================
# 4) Step B: Generate genome without motifs
# ============================================================
cat("Step B: Simulating genome without motifs...\n")

variants_list <- generate_all_variants(ppms, variant_distances)

primary_genome <- simulate_genome(
  n_genes, gene_length, cond_probs,
  background_probs, ppms, variants_list
)

cat("Genome simulation complete. Generated", n_genes, "genes.\n")


# ============================================================
# 5) Step C: Injection of motifs
# ============================================================
cat("Step C: Injecting motifs...\n")

# Precompute theoretical maximum log-odds (perfect match)
S_max <- vapply(ppms, ppm_max_logodds,
                background_probs = background_probs, numeric(1))

# Precompute thresholds: minimum log-odds among all motif variants
variant_thresholds <- vapply(names(ppms), function(name) {
  vs <- variants_list[[name]]
  if (length(vs) == 0) return(0)
  scores <- vapply(vs, function(seq_str) {
    seq_vec <- strsplit(seq_str, "", fixed = TRUE)[[1]]
    compute_match_score(seq_vec, ppms[[name]], background_probs)
  }, numeric(1))
  min(scores, na.rm = TRUE)
}, numeric(1))


# Function: inject motifs and export outputs for a given fraction
run_injection_for_fraction <- function(frac) {
  cat(sprintf(" -> Running injection for fraction %.2f\n", frac))
  
  target_n_each <- floor(n_genes * frac)
  current_genes <- primary_genome
  set_size <- length(ppms) * target_n_each
  
  # Pick which genes will get injected motifs
  pool  <- if (set_size > 0) sample(seq_len(n_genes), set_size, replace = FALSE) else integer(0)
  pools <- if (set_size > 0) split(pool, rep(names(ppms), each = target_n_each)) else
    setNames(vector("list", length(ppms)), names(ppms))
  
  injected_indices <- setNames(vector("list", length(ppms)), names(ppms))
  match_scores     <- setNames(vector("list", length(ppms)), names(ppms))
  
  bases <- c("A","C","G","T")
  
  for (motif_name in names(ppms)) {
    ppm_i <- ppms[[motif_name]]
    thr   <- variant_thresholds[[motif_name]]
    inds  <- pools[[motif_name]]
    
    scores_i <- numeric(length(inds))
    seqs_i   <- character(length(inds))  
    
    # Inject motif into selected genes
    for (j in seq_along(inds)) {
      idx <- inds[j]
      repeat {
        seq_str   <- sample_sequence_from_ppm(ppm_i)
        motif_seq <- strsplit(seq_str, "", fixed = TRUE)[[1]]
        sc        <- compute_match_score(motif_seq, ppm_i, background_probs)
        if (sc >= thr) break
      }
      scores_i[j] <- sc
      seqs_i[j]   <- seq_str
      current_genes <- inject_gene(current_genes, motif_seq, idx)
    }
    
    injected_indices[[motif_name]] <- inds
    match_scores[[motif_name]]     <- scores_i
    
    # Save per-position base frequencies for injected motifs
    L <- nchar(seqs_i[1L])
    freq_mat <- matrix(0, nrow = 4, ncol = L,
                       dimnames = list(bases, paste0("Pos", seq_len(L))))
    n <- length(seqs_i)
    for (pos in seq_len(L)) {
      chars <- substring(seqs_i, pos, pos)
      freq_mat[, pos] <- as.numeric(table(factor(chars, levels = bases))) / n
    }
    
    out_csv <- sprintf(file.path(out_dir, "base_freqs_%s_frac_%0.2f.csv"),
                       motif_name, frac)
    write.csv(round(freq_mat, 3), file = out_csv, row.names = TRUE)
  }
  
  # Simulate biological vector
  bio_vec <- simulate_bio_vector(
    n_genes          = length(current_genes),
    injected_indices = injected_indices,
    match_scores     = match_scores,
    meanlog          = meanlog,
    sdlog            = sdlog,
    multipliers      = multipliers,
    score_threshold  = S_max
  )
  
  # Export FASTA file
  fasta_lines <- unlist(mapply(function(i, gene_chars) {
    c(paste0(">gene_", i), paste0(gene_chars, collapse = ""))
  }, i = seq_along(current_genes), gene_chars = current_genes, SIMPLIFY = FALSE))
  
  fasta_path <- sprintf(file.path(out_dir, "simulated_genome_frac_%0.2f.fasta"), frac)
  writeLines(fasta_lines, con = fasta_path)
  
  # Export bio_vec CSV
  bio_df <- data.frame(
    gene_id   = paste0("gene_", seq_along(bio_vec)),
    bio_value = bio_vec
  )
  bio_path <- sprintf(file.path(out_dir, "bio_vec_frac_%0.2f.csv"), frac)
  write.csv(bio_df, file = bio_path, row.names = FALSE)
  
  cat("   Saved:", fasta_path, "and", bio_path, "\n")
  
}

# Run injection for all fractions requested
invisible(lapply(fractions, run_injection_for_fraction))


# ============================================================
# Done
# ============================================================
cat("Simulation finished. All outputs written to:", out_dir, "\n")


