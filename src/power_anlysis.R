# ============================================================
# Power analysis (Normal vs Normal): KS vs Mann–Whitney
# Two single-panel plots:
#   (A) Power vs sd1 with δ = 0  (pure variance change; mu1 = mu0)
#   (B) Power vs δ   with fixed sd1 (location shift)
#
# No Bonferroni. Uses %>% pipes. Self-contained.
# ============================================================

suppressWarnings({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ---------- Cliff's δ for two independent normals ----------
# For X ~ N(mu1, sd1^2), Y ~ N(mu0, sd0^2):
#   D = X - Y ~ N(mu1 - mu0, sd0^2 + sd1^2)
#   P(X>Y) = Phi( (mu1-mu0) / sqrt(sd0^2+sd1^2) )
#   δ = P(X>Y) - P(X<Y) = 2*Phi(z) - 1,  where z = (mu1-mu0)/sqrt(sd0^2+sd1^2)
delta_two_normals <- function(mu1, mu0, sd1, sd0) {
  z <- (mu1 - mu0) / sqrt(sd0^2 + sd1^2)
  2*pnorm(z) - 1
}

# Inverse: given δ and (mu0, sd0, sd1), solve for mu1
# From δ = 2*Phi(z)-1  =>  z = qnorm((δ+1)/2)
# => mu1 = mu0 + qnorm((δ+1)/2) * sqrt(sd0^2 + sd1^2)
mu1_from_delta_two_normals <- function(delta, mu0, sd0, sd1) {
  z <- qnorm((delta + 1)/2)
  mu0 + z * sqrt(sd0^2 + sd1^2)
}

# ---------- Tests ----------
pval_ks <- function(x, y)  suppressWarnings(ks.test(x, y)$p.value)
pval_mw <- function(x, y)  suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)

# ---------- Monte Carlo power (single configuration) ----------
# Draws n_effected from X ~ N(mu1, sd1), n_not_effected from Y ~ N(mu0, sd0)
# Returns one-row data.frame with KS and MW power estimates
estimate_power_once <- function(n_effected, N_total, reps, alpha,
                                mu0, sd0, mu1, sd1) {
  n_not_effected <- N_total - n_effected
  r_y <- function(n) rnorm(n, mean = mu0, sd = sd0)  # not effected group
  r_x <- function(n) rnorm(n, mean = mu1, sd = sd1)  # effected group
  
  sig_ks <- 0
  sig_mw <- 0
  
  for (r in seq_len(reps)) {
    x <- r_x(n_effected)
    y <- r_y(n_not_effected)
    pks <- pval_ks(x, y)
    pmw <- pval_mw(x, y)
    if (is.finite(pks) && pks < alpha) sig_ks <- sig_ks + 1
    if (is.finite(pmw) && pmw < alpha) sig_mw <- sig_mw + 1
  }
  data.frame(KS = sig_ks / reps, MW = sig_mw / reps, stringsAsFactors = FALSE)
}

# ============================================================
# (A) Power vs sd1 with δ = 0  (pure variance change)
#     Here enforce mu1 = mu0 so δ = 0 by construction.
# ============================================================

run_power_vs_sd1_normal <- function(
    sd1_grid,
    mu0, sd0,
    n_effected = 50, N_total = 1000,
    alpha = 0.05, reps = 1000,
    power_ref = 0.8
) {
  
  mu1 <- mu0  # δ=0 (equal means)
  
  res <- lapply(sd1_grid, function(sd1) {
    est <- estimate_power_once(
      n_effected = n_effected, N_total = N_total, reps = reps, alpha = alpha,
      mu0 = mu0, sd0 = sd0, mu1 = mu1, sd1 = sd1
    )
    cbind(sd1 = sd1, est)
  }) %>% bind_rows()
  
  res_long <- res %>% pivot_longer(c(KS, MW), names_to = "test", values_to = "power")
  
  p <- ggplot(res_long, aes(sd1, power, color = test)) +
    geom_line() + 
    geom_point(size = 1) +
    geom_hline(yintercept = power_ref, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 13) +
    labs(
      title = "Power vs. standard deviation change (KS vs. Mann–Whitney)",
      x = "Standard deviation of affected group",
      y = "Power",
      color = "Test"
    )
  
  list(plot = p, data = res)
}

# ============================================================
# (B) Power vs δ with fixed sd1  (location shift)
#     Here compute mu1 from δ using the closed-form inversion.
# ============================================================

run_power_vs_delta_normal <- function(
    delta_grid,
    mu0, sd0,
    sd1_fixed,
    n_effected = 50, N_total = 1000,
    alpha = 0.05, reps = 1000,
    power_ref = 0.8
) {
  # No feasibility constraint like mixtures: δ ∈ (-1,1) in theory; pick sensible range (e.g., 0..0.8)
  if (any(delta_grid <= -1 | delta_grid >= 1)) {
    stop("delta values must be in (-1, 1).", call. = FALSE)
  }
  
  res <- lapply(delta_grid, function(delta) {
    mu1 <- mu1_from_delta_two_normals(delta, mu0, sd0, sd1_fixed)
    est <- estimate_power_once(
      n_effected = n_effected, N_total = N_total, reps = reps, alpha = alpha,
      mu0 = mu0, sd0 = sd0, mu1 = mu1, sd1 = sd1_fixed
    )
    cbind(delta = delta, mu1 = mu1, est)
  }) %>% bind_rows()
  
  res_long <- res %>% pivot_longer(c(KS, MW), names_to = "test", values_to = "power")
  
  p <- ggplot(res_long, aes(delta, power, color = test)) +
  geom_line() +
  geom_point(size = 1) +
  geom_hline(yintercept = power_ref, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Power vs. effect size (KS vs. Mann–Whitney)",
    x = "Effect size (Cliff’s δ)",
    y = "Power",
    color = "Test"
  )

  
  list(plot = p, data = res)
}

# ============================
# Example usage
# ============================

set.seed(123)

# Question 1: δ = 0, sweep sd1 (variance-only change)
mu0 <- 5; sd0 <- 1
sd1_grid <- seq(1, 2.5, by = 0.1)
out_A <- run_power_vs_sd1_normal(
  sd1_grid = sd1_grid,
  mu0 = mu0, sd0 = sd0,
  n_effected = 50, N_total = 1000,
  alpha = 0.05, reps = 1000, power_ref = 0.8
)
print(out_A$plot)


# Question 2: fixed sd1, sweep δ (location change)
sd1_fixed  <- 1
delta_grid <- seq(0.00, 0.60, by = 0.05)  # keep (-1,1); choose sensible range
out_B <- run_power_vs_delta_normal(
  delta_grid = delta_grid,
  mu0 = mu0, sd0 = sd0, sd1_fixed = sd1_fixed,
  n_effected = 50, N_total = 1000,
  alpha = 0.05, reps = 1000, power_ref = 0.8
)
print(out_B$plot)





# ============================================================
# KS power vs n_effected (Normal vs Normal) with Bonferroni
# - Fixed delta (default 0.10)
# - alpha_eff = alpha / K
# - Sweeps n_effected; N_total fixed
# - Returns: plot + table + n_star + alpha_eff + mu1
# ============================================================

# --------- unchanged helper ----------
# estimate power for ONE delta and ONE n_effected (re-used inside the wrapper)
estimate_power_once_KS <- function(n_effected, N_total, reps, alpha_eff,
                                   mu0, sd0, mu1, sd1) {
  n_not_effected <- N_total - n_effected
  r_y <- function(n) rnorm(n, mean = mu0, sd = sd0)
  r_x <- function(n) rnorm(n, mean = mu1, sd = sd1)
  
  sig <- 0L
  for (r in seq_len(reps)) {
    x <- r_x(n_effected)
    y <- r_y(n_not_effected)
    p <- pval_ks(x, y)
    if (is.finite(p) && p < alpha_eff) sig <- sig + 1L
  }
  sig / reps
}

# --------- NEW: Sidák + vector of deltas ----------
run_ks_power_vs_n_effected_sidak <- function(
    mu0, sd0,
    sd1,
    deltas = c(0.10, 0.20, 0.30),        # vector of Cliff’s δ to scan
    n_effected_seq = seq(5, 150, by = 5),
    N_total = 7000,
    alpha = 0.05,
    K = 21760,                           # number of hypotheses (for Sidák)
    reps = 1000,
    power_ref = 0.80
) {
  # Šidák-adjusted alpha
  alpha_eff <- 1 - (1 - alpha)^(1 / K)
  
  # For each delta: compute mu1, sweep n_effected, estimate power
  one_delta <- function(delta) {
    mu1 <- mu1_from_delta_two_normals(delta, mu0, sd0, sd1)
    df <- lapply(n_effected_seq, function(n1) {
      pow <- estimate_power_once_KS(
        n_effected = n1, N_total = N_total, reps = reps, alpha_eff = alpha_eff,
        mu0 = mu0, sd0 = sd0, mu1 = mu1, sd1 = sd1
      )
      data.frame(delta = delta, n_effected = n1, power = pow, mu1 = mu1,
                 stringsAsFactors = FALSE)
    }) %>% bind_rows()
    df
  }
  
  res <- map_dfr(deltas, one_delta) %>%
    mutate(delta = as.numeric(delta))
  
  # minimal n_effected per delta reaching the target power
  n_star_df <- res %>%
    group_by(delta) %>%
    summarize(
      n_star = {
        idx <- match(TRUE, power >= power_ref, nomatch = NA_integer_)
        if (is.na(idx)) NA_integer_ else n_effected[idx]
      },
      .groups = "drop"
    )
  
  # Plot: one colored curve per delta
  p <- ggplot(res, aes(n_effected, power, color = factor(delta), group = delta)) +
    geom_line() + 
    geom_point(size = 1) +
    geom_hline(yintercept = power_ref, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Power of KS test under Šidák correction",
      x = sprintf("Affected group size (out of N_total = %d)", N_total),
      y = "Power",
      color = "Cliff’s δ"
    ) +
    theme_minimal(base_size = 13)
  
  list(
    plot   = p,
    table  = res,        # includes columns: delta, n_effected, power, mu1
    n_star = n_star_df,  # minimal n_effected per delta to reach power_ref
    alpha_eff = alpha_eff
  )
}

# --------------- Example usage ----------------
set.seed(123)

library(purrr)

out_ks_sidak <- run_ks_power_vs_n_effected_sidak(
  mu0 = 5, sd0 = 1,
  sd1 = 1,
  deltas = c(0.30, 0.40, 0.5),    # vector of Cliff’s δ to inspect
  n_effected_seq = seq(5, 200, by = 5),
  N_total = 7000,
  alpha = 0.05,
  K = 21760,
  reps = 1000,
  power_ref = 0.80
)

print(out_ks_sidak$plot)


#deltas = c(0.10, 0.20)
#n_effected_seq = seq(50, 2000, by = 50)


# deltas = c(0.30, 0.40, 0.50),    
# n_effected_seq = seq(5, 200, by = 5)


