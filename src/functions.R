
# Compute conditional probabilities (1st-order Markov model) 
# gives transition probabilities between nucleotides.
compute_conditional_probs <- function(dna_seq) {
  prev <- dna_seq[-length(dna_seq)]
  curr <- dna_seq[-1]
  bigram_table <- table(prev, curr)
  prop.table(bigram_table, margin = 1)
} 

# Compute background base frequencies (simple proportions of A/C/G/T).
compute_background_freq <- function(dna_vec) { 
  freqs <- table(dna_vec)
  freqs <- freqs / sum(freqs)
  return(freqs)
}

# Compute the maximum theoretical log-odds score (S_max) a motif can reach
# used later to normalize scores when simulating bio_vec.
ppm_max_logodds <- function(ppm, background_probs) {
  q <- as.numeric(background_probs[colnames(ppm)])
  s <- 0
  for (i in seq_len(nrow(ppm))) {
    s <- s + max(log2(ppm[i, ] / q))
  }
  s
}

# Sample one sequence directly from a PPM (row-wise).
sample_sequence_from_ppm <- function(ppm) {
  sequence <- character(nrow(ppm))
  for (i in seq_len(nrow(ppm))) {
    sequence[i] <- sample(c("A", "C", "G", "T"), size = 1, prob = ppm[i, ])
  }
  paste0(sequence, collapse = "")
}

# Generate all motif variants up to a given Hamming distance.
# ensures we know which sequences must be excluded from the background.
generate_all_variants <- function(ppms, distances) {
  bases <- c("A", "C", "G", "T")
  out <- vector("list", length(ppms))
  names(out) <- names(ppms)
  
  for (name in names(ppms)) {
    ppm <- ppms[[name]]
    motif_consensus <- apply(ppm, 1, function(p) bases[which.max(p)])
    dist_max <- distances[[name]]
    variants <- character(0)
    for (d in 0:dist_max) {
      variants_d <- generate_hamming_variants(motif_consensus, d, bases)
      variants <- unique(c(variants, variants_d))
    }
    out[[name]] <- variants
  }
  return(out)
}

# Generate all variants of a sequence at exact Hamming distance d.
generate_hamming_variants <- function(seq, d, alphabet = c("A", "C", "G", "T")) {
  if (d == 0) return(paste0(seq, collapse = ""))
  L <- length(seq)
  idx_combinations <- combn(L, d, simplify = FALSE)
  variants <- character(0)
  for (idxs in idx_combinations) {
    base_options <- lapply(idxs, function(i) setdiff(alphabet, seq[i]))
    replacements <- expand.grid(base_options, stringsAsFactors = FALSE)
    for (k in seq_len(nrow(replacements))) {
      variant <- seq
      variant[idxs] <- as.character(replacements[k, ])
      variants <- c(variants, paste0(variant, collapse = ""))
    }
  }
  unique(variants)
}

# Simulate synthetic genome:
# - Builds each gene base by base using Markov chain.
# - Rejects any sequence that matches a forbidden motif variant.
simulate_genome <- function(n_genes, gene_length, cond_probs,
                            background_probs, motifs, variants_list) {
  
  bases <- c("A", "C", "G", "T")
  variant_envs <- lapply(variants_list, function(vs) {
    e <- new.env(hash = TRUE, parent = emptyenv())
    if (length(vs)) for (s in vs) e[[s]] <- TRUE
    e
  })
  motif_lengths <- sapply(variants_list, function(vs) nchar(vs[1]))
  all_genes <- character(n_genes)
  
  for (g in seq_len(n_genes)) {
    gene_chars <- character(gene_length)
    i <- 1
    while (i <= gene_length) {
      probs <- if (i == 1) {
        as.numeric(background_probs[bases])
      } else {
        as.numeric(cond_probs[gene_chars[i - 1], bases])
      }
      rejects <- 0
      repeat {
        candidate <- sample(bases, size = 1, prob = probs)
        bad <- FALSE
        for (name in names(variants_list)) {
          L <- motif_lengths[[name]]
          if (i >= L) {
            window_str <- paste0(
              paste0(gene_chars[(i - L + 1):(i - 1)], collapse = ""),
              candidate
            )
            if (!is.null(variant_envs[[name]][[window_str]])) {
              bad <- TRUE
              break
            }
          }
        }
        if (!bad) {
          gene_chars[i] <- candidate
          i <- i + 1L
          break
        }
        rejects <- rejects + 1
        if (rejects > 10) {
          i <- max(1, i - L)
          break
        }
      }
    }
    all_genes[g] <- paste0(gene_chars, collapse = "")
  }
  return(all_genes)
}

# Compute log-odds score of a sequence vs. a PPM.
compute_match_score <- function(seq, ppm, background_probs) {
  score <- 0
  for (i in seq_along(seq)) {
    base <- seq[i]
    p <- ppm[i, base]
    q <- as.numeric(background_probs[base])
    score <- score + log2(p / q)
  }
  return(score) 
}

# Insert a motif sequence into a gene at a random position.
inject_gene <- function(all_genes, motif_seq, idx) {
  gene_str <- all_genes[[idx]]
  motif_str <- paste0(motif_seq, collapse = "")
  L <- nchar(motif_str)
  gene_len <- nchar(gene_str)
  pos <- sample(1:(gene_len - L + 1), 1)
  substr(gene_str, pos, pos + L - 1) <- motif_str
  all_genes[[idx]] <- gene_str
  return(all_genes)
}

# Simulate biological vector (bio_vec):
# - All genes start with baseline log-normal values.
# - For inserted motifs: shift meanlog upward
#   depending on motif score relative to its maximum (S_max).
simulate_bio_vector <- function(n_genes, injected_indices,
                                match_scores, meanlog, sdlog,
                                multipliers, score_threshold) {
  
  bio <- rlnorm(n_genes, meanlog = meanlog, sdlog = sdlog)
  for (motif in names(injected_indices)) {
    mul    <- multipliers[[motif]]
    idxs   <- injected_indices[[motif]]
    scores <- match_scores[[motif]]
    thr <- score_threshold[[motif]]
    for (j in seq_along(idxs)) {
      i <- idxs[j]
      s <- scores[j]
      bio[i] <- rlnorm(
        1,
        meanlog = meanlog + mul * (s / thr) * sdlog,
        sdlog   = sdlog
      )
    }
  }
  bio
}

