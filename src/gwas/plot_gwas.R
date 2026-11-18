# Manhattan + QQ plots for GEMMA results from plink.sh

setwd("/work/hs325/csgs2026")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggrepel)
  library(patchwork)
})

args      <- commandArgs(trailingOnly = TRUE)
gemma_dir <- if (length(args) >= 1) args[1] else "/work/hs325/csgs2026/gwas_out/gemma"
out_dir   <- if (length(args) >= 2) args[2] else "/work/hs325/csgs2026/gwas_out/plots"
p_col     <- if (length(args) >= 3) args[3] else "p_wald"

suggestive <- 1e-5
n_label    <- 10 # number of top SNPs labelled per run (only those above the suggestive line)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(gemma_dir, pattern = "\\.assoc\\.txt$", full.names = TRUE)
if (length(files) == 0) stop("No .assoc.txt files in ", gemma_dir)

read_gemma <- function(f) {
  fread(f) %>%
    as_tibble() %>%
    transmute(run = sub("\\.assoc\\.txt$", "", basename(f)),
              chr = as.integer(chr), rs, ps, af, beta, se,
              p = .data[[p_col]]) %>%
    filter(!is.na(p), !is.na(chr))
}

gwas <- map_df(files, read_gemma)

# load marker info
chr_info <- gwas %>%
  group_by(chr) %>%
  summarise(chr_len = max(ps), .groups = "drop") %>%
  arrange(chr) %>%
  mutate(offset = lag(cumsum(as.numeric(chr_len)), default = 0))

gwas <- gwas %>%
  left_join(chr_info, by = "chr") %>%
  mutate(pos_cum = ps + offset, logp = -log10(p))

axis_df <- chr_info %>% mutate(center = offset + chr_len / 2)

# generate p value stats
run_stats <- gwas %>%
  group_by(run) %>%
  summarise(
    n_snps = n(),
    bonf   = 0.05 / n(),
    lambda = median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1),
    n_bonf = sum(p < 0.05 / n()),
    n_sugg = sum(p < suggestive),
    .groups = "drop"
  )
print(as.data.frame(run_stats), row.names = FALSE)

# plot
manhattan_plot <- function(df, stats) {
  top <- df %>% filter(p < suggestive) %>% slice_min(p, n = n_label, with_ties = FALSE)
  ymax <- max(c(df$logp, -log10(stats$bonf))) * 1.08

  ggplot(df, aes(pos_cum, logp, colour = factor(chr %% 2))) +
    geom_point(size = 0.6, alpha = 0.8) +
    geom_hline(yintercept = -log10(stats$bonf), colour = "#C0392B", linewidth = 0.4) +
    geom_hline(yintercept = -log10(suggestive), colour = "grey50", linetype = "dashed", linewidth = 0.4) +
    geom_point(data = top, colour = "#C0392B", size = 1.2) +
    geom_text_repel(data = top, aes(label = rs), colour = "black", size = 2.5,
                    max.overlaps = Inf, min.segment.length = 0) +
    scale_colour_manual(values = c("0" = "#2C5F8A", "1" = "#8DB3D6"), guide = "none") +
    scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr, expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    labs(title = df$run[1],
         subtitle = sprintf("%s SNPs | lambda = %.3f | Bonferroni p < %.2e: %d | p < 1e-5: %d",
                            format(stats$n_snps, big.mark = ","), stats$lambda,
                            stats$bonf, stats$n_bonf, stats$n_sugg),
         x = "Chromosome", y = expression(-log[10](italic(p)))) +
    theme_classic(base_size = 11) +
    theme(panel.grid.major.y = element_line(colour = "grey92"),
          plot.subtitle = element_text(size = 9, colour = "grey30"))
}

qq_plot <- function(df, stats) {
  qq <- tibble(obs = sort(df$logp, decreasing = TRUE),
               exp = -log10(ppoints(nrow(df))))
  ggplot(qq, aes(exp, obs)) +
    geom_abline(slope = 1, intercept = 0, colour = "#C0392B", linewidth = 0.4) +
    geom_point(size = 0.6, colour = "#2C5F8A") +
    labs(title = paste(df$run[1], "QQ"),
         subtitle = sprintf("lambda = %.3f", stats$lambda),
         x = expression(Expected ~ -log[10](italic(p))),
         y = expression(Observed ~ -log[10](italic(p)))) +
    theme_classic(base_size = 11)
}

# compile plots
manhattans <- list()
for (r in unique(gwas$run)) {
  df    <- gwas %>% filter(run == r)
  stats <- run_stats %>% filter(run == r)

  p_man <- manhattan_plot(df, stats)
  manhattans[[r]] <- p_man
  ggsave(file.path(out_dir, paste0(r, "_manhattan.png")), p_man, width = 11, height = 4, dpi = 300)
  ggsave(file.path(out_dir, paste0(r, "_qq.png")), qq_plot(df, stats), width = 4.5, height = 4.5, dpi = 300)

  df %>%
    arrange(p) %>%
    slice_head(n = 50) %>%
    mutate(p_bonf = pmin(p * stats$n_snps, 1)) %>%
    select(rs, chr, ps, af, beta, se, p, p_bonf) %>%
    write_tsv(file.path(out_dir, paste0(r, "_top_hits.tsv")))
}

ggsave(file.path(out_dir, "all_manhattan.png"),
       wrap_plots(manhattans, ncol = 1),
       width = 11, height = 4 * length(manhattans), dpi = 300, limitsize = FALSE)

write_tsv(run_stats, file.path(out_dir, "run_summary.tsv"))
cat("Plots written to", out_dir, "\n")
