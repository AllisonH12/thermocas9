#!/usr/bin/env Rscript
# Canonical R limma::lmFit + limma::eBayes parity check.
#
# Usage:
#   Rscript scripts/r_limma_parity.R <input.tsv> <output.tsv>
#
# Input TSV layout:
#   col 1            : probe_id
#   cols 2..(n+1)    : sample β values (probes × samples matrix)
# Group vector is read from <input.tsv>.group (one int per line in
# column order: 1 = tumor, 0 = normal).
#
# Output TSV columns: probe_id, delta_beta, s_sq, s_tilde_sq, t_mod,
# df_tot, p_value — matching the Python `scripts/limma_ebayes.py`
# output schema so a downstream comparator can join the two by
# probe_id and compare per-probe statistics directly.
#
# This is a parity test for the pure-Python Smyth (2004)
# implementation in scripts/limma_ebayes.py. It runs the canonical
# limma::lmFit(beta_matrix, design) %>% eBayes() pipeline against
# the same (β, group) inputs the Python pipeline consumed.

suppressPackageStartupMessages({
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript scripts/r_limma_parity.R <input.tsv> <output.tsv>")
}
input_path <- args[[1]]
output_path <- args[[2]]
group_path <- paste0(input_path, ".group")

cat(sprintf("reading β matrix from %s\n", input_path))
mat <- as.matrix(read.table(
  input_path, header = TRUE, sep = "\t", row.names = 1L,
  check.names = FALSE, comment.char = ""
))
cat(sprintf("  matrix shape: %d probes x %d samples\n", nrow(mat), ncol(mat)))

group <- as.integer(scan(group_path, what = integer(0), quiet = TRUE))
stopifnot(length(group) == ncol(mat))
cat(sprintf("  group=1 (tumor): %d, group=0 (normal): %d\n",
            sum(group == 1L), sum(group == 0L)))

# Canonical limma fit. Design: intercept = normal mean, second column =
# tumor-vs-normal contrast. lmFit on the β matrix (matches the Python
# implementation's β-input choice; M-value transform is documented as a
# downstream Bioconductor convention but not required for the math).
design <- cbind(intercept = 1L, tumor_vs_normal = group)
fit <- lmFit(mat, design)
fit_eb <- eBayes(fit)

# tumor_vs_normal coefficient = β_tumor − β_normal. The Python pipeline
# emits delta_beta as (β_tumor − β_normal); use the same sign here.
delta_beta <- fit_eb$coefficients[, "tumor_vs_normal"]
s_sq <- fit$sigma^2                    # residual variance (per probe)
s_tilde_sq <- fit_eb$s2.post           # eBayes posterior variance (per probe)
t_mod <- fit_eb$t[, "tumor_vs_normal"] # moderated t-statistic
df_tot <- fit_eb$df.total              # = df.residual + df.prior, per probe
p_value <- fit_eb$p.value[, "tumor_vs_normal"]

# eBayes prior parameters (for the .meta companion).
s0_sq <- fit_eb$s2.prior
d0 <- fit_eb$df.prior

cat(sprintf("  d = %g (residual df, constant)\n", fit_eb$df.residual[1]))
cat(sprintf("  prior: s0² = %g, d0 = %g\n", s0_sq, d0))

cat(sprintf("writing %s\n", output_path))
out <- data.frame(
  probe_id = rownames(mat),
  delta_beta = delta_beta,
  s_sq = s_sq,
  s_tilde_sq = s_tilde_sq,
  t_mod = t_mod,
  df_tot = df_tot,
  p_value = p_value,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.table(out, output_path, sep = "\t", quote = FALSE, row.names = FALSE)

meta_path <- sub("\\.tsv$", ".meta.json", output_path)
cat(sprintf("writing %s\n", meta_path))
meta <- list(
  source = "R limma::lmFit + limma::eBayes",
  r_version = as.character(getRversion()),
  limma_version = as.character(packageVersion("limma")),
  n_probes = nrow(mat),
  n_samples_used = ncol(mat),
  n0 = sum(group == 0L),
  n1 = sum(group == 1L),
  d = fit_eb$df.residual[1],
  s0_sq = s0_sq,
  d0 = d0,
  d0_is_infinite = is.infinite(d0)
)
# Hand-rolled JSON to avoid pulling jsonlite as a dependency.
writeLines(c(
  "{",
  sprintf("  \"source\": \"%s\",", meta$source),
  sprintf("  \"r_version\": \"%s\",", meta$r_version),
  sprintf("  \"limma_version\": \"%s\",", meta$limma_version),
  sprintf("  \"n_probes\": %d,", meta$n_probes),
  sprintf("  \"n_samples_used\": %d,", meta$n_samples_used),
  sprintf("  \"n0\": %d,", meta$n0),
  sprintf("  \"n1\": %d,", meta$n1),
  sprintf("  \"d\": %g,", meta$d),
  sprintf("  \"s0_sq\": %g,", meta$s0_sq),
  sprintf("  \"d0\": %s,", if (meta$d0_is_infinite) "null" else format(meta$d0, digits = 10)),
  sprintf("  \"d0_is_infinite\": %s", if (meta$d0_is_infinite) "true" else "false"),
  "}"
), meta_path)

cat("done.\n")
