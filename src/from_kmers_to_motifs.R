# =============================================================
# Packages
# =============================================================

library(tidyverse)
library(cluster)      
library(Biostrings)   
library(msa)          
library(ggseqlogo)    
library(Rtsne)


# ----------------------------
# 0) Input prep
# ----------------------------
# Load KS-test results for all k-mers.
# Compute log2 fold change = log2(mean_present / mean_absent).
# Keep only top 10% by |log2FC| and pval_fdr_log ≥ 2 (≈ p ≤ 0.01).
# Split k-mers by effect direction:
#   >0 → positive (higher when present)
#   <0 → negative (lower when present)
# ----------------------------

df <- read_csv(
  "data/ks_results/frac_0.05_ks_results.csv",
  show_col_types = FALSE
)

df <- df %>%
  mutate(log2FC = log2(mean_present / mean_absent))

threshold <- quantile(abs(df$log2FC), 0.9, na.rm = TRUE)

filtered_df <- df %>%
  filter(abs(log2FC) >= threshold & pval_fdr_log >= 2)

positive_kmers <- filtered_df$kmer[filtered_df$log2FC > 0]
negative_kmers <- filtered_df$kmer[filtered_df$log2FC < 0]

# ----------------------------
# Load genomic sequences
# ----------------------------
# Load the simulated genome (FASTA file).
# This will be used to extract actual genomic fragments that contain the significant kmers.
fasta_path <- "data/simulation_outputs/simulated_genome_frac_0.05.fasta"
seqs <- readDNAStringSet(fasta_path) 

# -------------------------------------------------------------
# Function to extract and merge overlapping k-mer hits
# -------------------------------------------------------------
# For each sequence in the genome:
# 1. Find all positions where any of the kmers occur.
# 2. Keep only true matches.
# 3. Merge overlapping or adjacent matches into larger genomic fragments.
extract_kmer_hits_union <- function(seqs, kmers) {
  if (length(kmers) == 0) return(character(0))
  
  # Regex pattern: lookahead for any of the kmers
  pattern <- paste0("(?=(", paste(kmers, collapse = "|"), "))")
  results <- list(); counter <- 1
  
  for (s in as.character(seqs)) {
    m <- gregexpr(pattern, s, perl = TRUE)[[1]]
    if (is.na(m[1]) || m[1] == -1) next  # skip if no matches
    
    starts <- m
    matched_kmers <- character(length(starts))
    
    # Confirm exact k-mer match at each start position
    for (i in seq_along(starts)) {
      found <- NA
      for (kmer in kmers) {
        subseq <- substr(s, starts[i], starts[i] + nchar(kmer) - 1)
        if (subseq == kmer) { found <- kmer; break }
      }
      matched_kmers[i] <- found
    }
    
    # Remove NA matches
    keep <- !is.na(matched_kmers) & !is.na(starts)
    starts <- starts[keep]
    matched_kmers <- matched_kmers[keep]
    if (length(starts) == 0) next
    
    # Compute end positions
    ends <- starts + nchar(matched_kmers) - 1
    keep2 <- !is.na(ends)
    starts <- starts[keep2]; ends <- ends[keep2]
    if (length(starts) == 0) next
    
    # Sort hits by start coordinate
    ord <- order(starts)
    starts <- starts[ord]; ends <- ends[ord]
    
    # Merge overlapping fragments
    merged <- list()
    cur_start <- starts[1]; cur_end <- ends[1]
    if (length(starts) > 1) {
      for (j in 2:length(starts)) {
        if (!is.na(starts[j]) && starts[j] <= cur_end + 1) {
          # Overlapping or adjacent → extend the fragment
          cur_end <- max(cur_end, ends[j], na.rm = TRUE)
        } else {
          # Save completed fragment
          merged[[length(merged) + 1]] <- substr(s, cur_start, cur_end)
          # Start a new fragment
          cur_start <- starts[j]; cur_end <- ends[j]
        }
      }
    }
    # Save final fragment
    merged[[length(merged) + 1]] <- substr(s, cur_start, cur_end)
    
    for (frag in merged) {
      results[[counter]] <- frag
      counter <- counter + 1
    }
  }
  unlist(results, use.names = FALSE)
}

# =============================================================
# Sliding q-gram distance and clustering helpers
# =============================================================

# ---- Sliding q-gram distance ----
# Compare two sequences by aligning the shorter against the longer one.
# For each offset, compare q-grams (length q substrings).
# Distance = 1 - (best fraction of matching q-grams).
sliding_q_dist <- function(s1, s2, q = 3) {
  n1 <- stringr::str_length(s1); n2 <- stringr::str_length(s2)
  if (n1 < q || n2 < q) return(1)  # too short → maximal distance
  
  # Ensure 'short' is the shorter string
  if (n1 <= n2) { short <- s1; long <- s2 } else { short <- s2; long <- s1 }
  m <- stringr::str_length(short); L <- stringr::str_length(long)
  nwin <- m - (q - 1)  # number of q-grams in short
  
  short_qgrams <- stringr::str_sub(short, 1:nwin, q:m)
  long_qgrams  <- stringr::str_sub(long, 1:(L - (q - 1)), q:L)
  
  best <- 0L
  for (off in 0:(L - m)) {
    matches <- sum(short_qgrams == long_qgrams[(1:nwin) + off])
    if (matches > best) best <- matches
    if (best == nwin) break  # perfect match → early stop
  }
  1 - best / nwin
}

# ---- Distance matrix for strings ----
# Compute full pairwise distance matrix for a set of sequences
# using sliding_q_dist.
dist_matrix_sliding_q <- function(strings, q = 3) {
  n <- length(strings)
  D <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      d_ij <- sliding_q_dist(strings[i], strings[j], q = q)
      D[i, j] <- D[j, i] <- d_ij
    }
  }
  stats::as.dist(D)
}

# ---- Pick k with standard ASW (Average Silhouette Width) ----
# For each candidate k (number of clusters), compute PAM clustering
# and evaluate mean silhouette width (cluster quality).
# The k with the maximum ASW is chosen.
pick_k_asw <- function(d, k_range = 2:8, group) {
  n <- attr(d, "Size")
  k_range <- k_range[k_range >= 2 & k_range <= max(2, n - 1)]
  asw  <- sapply(k_range, function(k) {
    fit <- pam(d, k, diss = TRUE)
    mean(silhouette(fit$clustering, d)[, "sil_width"])
  })
  
  df <- tibble(k = k_range, ASW = asw)
  print(
    ggplot(df, aes(k, ASW)) +
      geom_line() + geom_point(size = 3) +
      labs(title = paste("Average Silhouette Width across candidate k values -", group),
           x = "k", y = "Average Silhouette Width") +
      theme_minimal(base_size = 14)
  )
  df$k[which.max(df$ASW)]
}

# ---- PAM clustering with ASW-based k ----
# Run PAM clustering on a set of sequences using sliding q-gram distances.
# Number of clusters is chosen by ASW.
cluster_group_sliding_q <- function(strings, group_name, k_range = 2:8, q = 3) {
  d <- dist_matrix_sliding_q(strings, q = q)
  k <- pick_k_asw(d, k_range, group_name)
  fit <- cluster::pam(d, k, diss = TRUE)
  tibble::tibble(
    seq     = strings,
    cluster = paste0(group_name, "_", fit$clustering)
  )
}

# =============================================================
# 1) Extract windows once from genome
# =============================================================
# Extract genomic fragments containing the significant kmers.
# Store both full set (with duplicates) and unique set (for clustering).
pos_windows_all <- extract_kmer_hits_union(seqs, positive_kmers)
neg_windows_all <- extract_kmer_hits_union(seqs, negative_kmers)

pos_windows_unique <- unique(pos_windows_all)
neg_windows_unique <- unique(neg_windows_all)

# =============================================================
# 2) Cluster the windows
# =============================================================
# Perform clustering separately for positive and negative groups.
# Use sliding q-gram distance with q = 3.
k_grid <- 2:8
result_pos <- cluster_group_sliding_q(pos_windows_unique, "Positive", k_grid, q = 3)
result_neg <- cluster_group_sliding_q(neg_windows_unique, "Negative", k_grid, q = 3)

# Combine both clustering results into a single tibble
final_clusters <- bind_rows(result_pos, result_neg)

# Map cluster labels back to all (non-unique) sequences
final_clusters_expanded <- tibble(seq = c(pos_windows_all, neg_windows_all)) %>%
  left_join(final_clusters, by = "seq")

# =============================================================
# 3) t-SNE visualization (right after clustering)
# =============================================================
# Visualize sequences in 2D using t-SNE applied to the distance matrix.
plot_tsne_group <- function(final_clusters, group_name, perplexity = 30, theta = 0.5) {
  # Subset only the chosen group and keep unique sequences
  group_data <- final_clusters %>% filter(grepl(group_name, cluster))
  strings <- unique(group_data$seq)
  
  if (length(strings) < 3) {
    warning("Not enough sequences to run t-SNE for group: ", group_name)
    return(NULL)
  }
  
  # Distance matrix with sliding q-gram
  d <- dist_matrix_sliding_q(strings, q = 3)
  
  # Run t-SNE
  set.seed(42)
  tsne_out <- Rtsne::Rtsne(
    as.matrix(d),
    is_distance = TRUE,
    perplexity = min(perplexity, floor((length(strings) - 1) / 3)),
    theta = theta,
    verbose = TRUE
  )
  
  # Prepare dataframe for plotting
  plot_df <- tibble(
    seq = strings,
    tSNE1 = tsne_out$Y[, 1],
    tSNE2 = tsne_out$Y[, 2],
    cluster = group_data$cluster[match(strings, group_data$seq)]
  )
  
  # Draw scatter plot
  ggplot(plot_df, aes(tSNE1, tSNE2, color = cluster)) +
    geom_point(alpha = 0.7, size = 3) +
    theme_minimal(base_size = 14) +
    labs(title = paste("t-SNE visualization for", group_name, "group"))
}

# Run t-SNE immediately after clustering
plot_tsne_group(final_clusters_expanded, "Positive", perplexity = 10, theta = 0)
plot_tsne_group(final_clusters_expanded, "Negative", perplexity = 30, theta = 0)

# =============================================================
# 4) Motif Discovery per Cluster
# =============================================================
# For each cluster:
# 1. Collect all sequences belonging to it
# 2. Sample up to 1000 sequences (to keep alignment feasible)
# 3. Perform multiple sequence alignment
# 4. Filter out gap-dominated columns (>50% gaps)
# 5. Build a position frequency matrix (PFM)
# 6. Plot sequence logos
all_clusters <- unique(final_clusters_expanded$cluster)
motif_list <- list()

for (cluster_name in all_clusters) {
  message("--- Analyzing cluster: ", cluster_name, " ---")
  
  # Collect sequences in cluster
  cluster_seqs <- final_clusters_expanded %>%
    filter(cluster == cluster_name) %>%
    pull(seq)
  
  dna_set <- DNAStringSet(cluster_seqs)
  
  # Random sampling (max 1000 sequences)
  set.seed(42)
  subset_seqs <- sample(dna_set, size = min(length(dna_set), 1000))
  
  # Multiple Sequence Alignment
  aligned_seqs <- msaMuscle(subset_seqs)
  
  # Convert to matrix
  aligned <- as.character(unmasked(aligned_seqs))
  alignment_matrix <- do.call(rbind, strsplit(aligned, ""))
  
  # Remove gap-dominated columns (>50% gaps)
  gap_fractions <- apply(alignment_matrix, 2, function(col) mean(col == "-"))
  keep_cols <- which(gap_fractions < 0.5)
  alignment_matrix_filtered <- alignment_matrix[, keep_cols, drop = FALSE]
  
  # Position Frequency Matrix (PFM)
  pfm <- apply(alignment_matrix_filtered, 2, function(col) {
    table(factor(col, levels = c("A", "C", "G", "T")))
  })
  pfm <- as.matrix(pfm)
  pfm <- sweep(pfm, 2, colSums(pfm), FUN = "/")
  
  motif_list[[cluster_name]] <- pfm
}

# Plot all motif logos in one figure
ggseqlogo(motif_list, method = "prob") +
  labs(title = "Sequence Logos for All Clusters")




