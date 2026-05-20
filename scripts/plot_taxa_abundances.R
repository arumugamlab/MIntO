#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# plot_taxa_abundances.R
#
# Purpose:
#   Given a phyloseq object and a set of comma-delimited taxonomic keywords,
#   this script identifies all taxa whose names match any keyword across
#   *any* taxonomic rank (Species, Genus, Family, etc.), computes their
#   relative abundances, and generates a multi-page PDF with one page per taxon.
#
#   The script supports two visualization modes:
#
#   1) Time-mode (when --time_var is provided):
#        - Shows relative abundance over time.
#        - Raw data are shown as points.
#        - LOESS-smoothed trends are shown as curves.
#        - Optional colouring by --group_var.
#        - Optional faceting by --facet_var.
#
#   2) Non-time-mode (when --time_var is NOT provided):
#        - Requires --group_var.
#        - Default (--summary_stat=median):
#        - --summary_stat=mean: mean bars with SE whiskers + jittered raw points.
#        - --summary_stat=median: boxplot + jittered raw points.
#        - Optional faceting by --facet_var.
#
#   All output is written to a single multi-page PDF, with:
#        - One taxon per page.
#        - Fixed y-axis in [0, 1] by default (relative abundance),
#          or per-taxon free scaling if --free_y is enabled.
#
# Key behaviors and design choices:
#   - Keyword matching is case-insensitive and scans all taxonomic ranks.
#   - Taxa with all-zero relative abundance after filtering are dropped and
#     reported to STDERR.
#   - Species names are assumed to already be in "Genus species" format.
#   - Canonical taxonomic ranks used here: Species, Genus, Family - will match 
#     columns in sample_data(phyloseq) case-insensitively.
#   - If only Genus/Family-level annotation is available, labels are disambiguated
#     using the taxon row name: "Genus (taxon_id)" or "Family (taxon_id)".
#
# Typical usage:
#   Rscript plot_taxa_abundances.R \
#     --phyloseq mydata.rds \
#     --keywords "Enterococcus,Lachnospiraceae" \
#     --time_var Day \
#     --group_var Treatment \
#     --facet_var Donor \
#     --out taxa_trends.pdf
#
#   Rscript plot_taxa_abundances.R \
#     --phyloseq mydata.rds \
#     --keywords "Enterococcus,Lachnospiraceae" \
#     --group_var Treatment \
#     --summary_stat mean \
#     --out taxa_bars.pdf
#
# Author: Mani Arumugam
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------
# Command-line options
# -----------------------------

option_list <- list(
  make_option(c("-p", "--phyloseq"),
              type = "character", default = NULL,
              help = "Path to RDS file containing a phyloseq object", metavar = "FILE"),
  make_option(c("-k", "--keywords"),
              type = "character", default = NULL,
              help = "Comma-separated taxon name keywords (matched against all taxonomic ranks)"),
  make_option(c("-t", "--time_var"),
              type = "character", default = NULL,
              help = "Sample metadata column to use as time (numeric). If provided, curves are drawn over time."),
  make_option(c("-g", "--group_var"),
              type = "character", default = NULL,
              help = "Sample metadata column used as main grouping variable (required if --time_var is missing)."),
  make_option(c("-f", "--facet_var"),
              type = "character", default = NULL,
              help = "Optional sample metadata column used for faceting within each taxon."),
  make_option(c("-s", "--summary_stat"),
              type = "character", default = "median",
              help = "Summary statistic when --time_var is missing: one of 'mean', 'median' [default: %default]"),
  make_option(c("-o", "--out"),
              type = "character", default = "taxa_relabund.pdf",
              help = "Output PDF filename [default: %default]"),
  make_option(c("--width"),
              type = "double", default = 8,
              help = "Figure width in inches [default: %default]"),
  make_option(c("--height"),
              type = "double", default = 6,
              help = "Figure height in inches [default: %default]"),
  make_option(c("--free_y"),
              action = "store_true", default = FALSE,
              help = "Allow y-axis to vary per taxon (default: fixed 0-1)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# -----------------------------
# Basic checks
# -----------------------------

if (is.null(opt$phyloseq)) {
  stop("Error: --phyloseq must be provided and point to a phyloseq RDS file.", call. = FALSE)
}
if (!file.exists(opt$phyloseq)) {
  stop(sprintf("Error: file not found: %s", opt$phyloseq), call. = FALSE)
}

if (is.null(opt$keywords) || nchar(trimws(opt$keywords)) == 0) {
  stop("Error: --keywords must be provided and non-empty.", call. = FALSE)
}

# normalize output filename to .pdf
out_file <- opt$out
if (tolower(tools::file_ext(out_file)) != "pdf") {
  out_file <- paste0(out_file, ".pdf")
}
opt$out <- out_file

summary_stat <- tolower(opt$summary_stat)
if (!summary_stat %in% c("mean", "median")) {
  stop("Error: --summary_stat must be one of 'mean', 'median'.", call. = FALSE)
}

# -----------------------------
# Load phyloseq object
# -----------------------------

ps <- readRDS(opt$phyloseq)
if (!inherits(ps, "phyloseq")) {
  stop("Error: the RDS file does not contain a phyloseq object.", call. = FALSE)
}

# -----------------------------
# Transform to relative abundance and drop all-zero taxa
# -----------------------------

ps_rel <- transform_sample_counts(ps, function(x) {
  s <- sum(x)
  if (s == 0) x else x / s
})

tax_sums <- taxa_sums(ps_rel)
zero_taxa <- names(tax_sums)[tax_sums == 0]

if (length(zero_taxa) > 0) {
  # message() goes to stderr by default
  message(sprintf("Dropping %d taxa with all-zero relative abundance: %s",
                  length(zero_taxa), paste(zero_taxa, collapse = ", ")))
  ps_rel <- prune_taxa(tax_sums > 0, ps_rel)
}

if (ntaxa(ps_rel) == 0) {
  stop("Error: all matched taxa had zero relative abundance; nothing to plot.", call. = FALSE)
}

# -----------------------------
# Check sample_data columns
# -----------------------------

sd <- as(sample_data(ps_rel), "data.frame")

# convenience checker
check_var_exists <- function(var_name, context) {
  if (!is.null(var_name) && !(var_name %in% colnames(sd))) {
    stop(sprintf("Error: %s '%s' not found in sample_data(). Available: %s",
                 context, var_name, paste(colnames(sd), collapse = ", ")), call. = FALSE)
  }
}

check_var_exists(opt$time_var, "time_var")
check_var_exists(opt$group_var, "group_var")
check_var_exists(opt$facet_var, "facet_var")

# -----------------------------
# Determine mode (time vs non-time)
# -----------------------------

time_mode <- !is.null(opt$time_var)

if (!time_mode) {
  # non-time-mode: group_var is required
  if (is.null(opt$group_var)) {
    stop("Error: --group_var is required because --time_var was not provided.", call. = FALSE)
  }
}

# If time_mode, ensure time_var is numeric
if (time_mode) {
  tv <- sd[[opt$time_var]]
  if (!is.numeric(tv)) {
    stop(sprintf("Error: time_var '%s' must be numeric. Please convert it before calling this script.",
                 opt$time_var), call. = FALSE)
  }
}

# -----------------------------
# Parse keywords, match taxa
# -----------------------------

keywords <- strsplit(opt$keywords, ",")[[1]]
keywords <- trimws(keywords)
keywords <- keywords[keywords != ""]
if (length(keywords) == 0) {
  stop("Error: no valid keywords parsed from --keywords.", call. = FALSE)
}

tt <- as.data.frame(tax_table(ps_rel), stringsAsFactors = FALSE)
if (nrow(tt) == 0) {
  stop("Error: tax_table(ps_rel) is empty; cannot match keywords.", call. = FALSE)
}

matched <- rep(FALSE, nrow(tt))

for (kw in keywords) {
  for (col in colnames(tt)) {
    vals <- tt[[col]]
    vals <- ifelse(is.na(vals), "", as.character(vals))
    hit <- grepl(kw, vals, ignore.case = TRUE)
    hit[is.na(hit)] <- FALSE
    matched <- matched | hit
  }
}

if (!any(matched)) {
  stop("Error: no taxa matched the provided keywords across any taxonomic ranks.", call. = FALSE)
}

taxa_keep <- rownames(tt)[matched]
message(sprintf("Matched %d taxa for keywords: %s",
                length(taxa_keep), paste(keywords, collapse = ", ")))

ps_sub <- prune_taxa(taxa_keep, ps_rel)

# -----------------------------
# Build taxon labels
# -----------------------------

tt_sub <- as.data.frame(tax_table(ps_sub), stringsAsFactors = FALSE)
tax_names <- sort(taxa_names(ps_sub))

message(sprintf("Matched taxa: %s", paste(tax_names, collapse = ", ")))

# Build a case-insensitive map from canonical ranks -> actual column names

canonical_ranks <- c("Species", "Genus", "Family")
tax_cols <- colnames(tt_sub)

rank_map <- setNames(rep(NA_character_, length(canonical_ranks)), canonical_ranks)
for (r in canonical_ranks) {
  hit <- which(tolower(tax_cols) == tolower(r))
  if (length(hit) > 0) {
    rank_map[r] <- tax_cols[hit[1]]  # take the first match if there are multiple
  }
}
message("Rank map: ", paste(names(rank_map), rank_map, sep = "->", collapse = ", "))

# Get canonical ranks for given taxon
get_rank_val <- function(row, canonical_rank) {
  colname <- rank_map[canonical_rank]
  if (is.na(colname) || !(colname %in% colnames(row))) {
    return(NA_character_)
  }
  v <- row[[colname]]
  if (is.null(v) || is.na(v) || !nzchar(v)) {
    return(NA_character_)
  }
  as.character(v)
}

# Make good taxon labels
make_taxon_label <- function(i) {
  rowname <- tax_names[i]
  row <- tt_sub[i, , drop = FALSE]

  # Species already assumed "Genus species" if present
  sp  <- get_rank_val(row, "Species")
  if (!is.na(sp)) {
    return(sp)
  }

  gn  <- get_rank_val(row, "Genus")
  if (!is.na(gn)) {
    return(gn)
  }

  fam <- get_rank_val(row, "Family")
  if (!is.na(fam)) {
    # include rowname to keep labels unique within a family
    return(paste0(fam, " (", rowname, ")"))
  }

  # fall back to taxon ID
  return(rowname)
}

taxon_labels <- vapply(seq_along(tax_names), make_taxon_label, character(1))
names(taxon_labels) <- tax_names

# -----------------------------
# Melt to long format
# -----------------------------

df <- psmelt(ps_sub)
# psmelt uses "OTU" as taxa column
if (!"OTU" %in% colnames(df)) {
  stop("Error: psmelt output does not contain 'OTU' column; unexpected phyloseq structure.", call. = FALSE)
}

df$taxon_label <- taxon_labels[as.character(df$OTU)]

# -----------------------------
# PDF device
# -----------------------------

pdf(opt$out, width = opt$width, height = opt$height)

on.exit({
  dev.off()
}, add = TRUE)

# -----------------------------
# Helper: y-scale
# -----------------------------

apply_y_scale <- function(p, df_tax) {
  if (opt$free_y) {
    # free y per taxon: ensure lower bound at 0, upper bound from data
    max_ab <- max(df_tax$Abundance, na.rm = TRUE)
    if (!is.finite(max_ab) || max_ab <= 0) {
      max_ab <- 1
    }
    p + scale_y_continuous(limits = c(0, max_ab * 1.05),
                           name = "Relative abundance")
  } else {
    # fixed y across all taxa
    p + scale_y_continuous(limits = c(0, 1),
                           name = "Relative abundance")
  }
}

# -----------------------------
# Loop over taxa and plot
# -----------------------------

taxon_order <- sort(unique(df$taxon_label))

for (tx in taxon_order) {
  df_tax <- df[df$taxon_label == tx, , drop = FALSE]
  if (nrow(df_tax) == 0) next

  if (time_mode) {
    # -------------------------
    # Time-mode: points + smooth
    # -------------------------
    aes_args <- list(
      x = as.name(opt$time_var),
      y = quote(Abundance)
    )
    if (!is.null(opt$group_var)) {
      aes_args$colour <- as.name(opt$group_var)
    }

    p <- ggplot(df_tax, do.call(aes, aes_args)) +
      geom_point(alpha = 0.4, size = 1.5) +
      geom_smooth(se = FALSE, method = "loess") +
      labs(
        title = paste0("Relative abundance over time: ", tx),
        x = opt$time_var,
        colour = if (!is.null(opt$group_var)) opt$group_var else NULL
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    if (!is.null(opt$facet_var)) {
      p <- p + facet_wrap(as.formula(paste("~", opt$facet_var)))
    }

    p <- apply_y_scale(p, df_tax)
    print(p)

  } else {
    # -------------------------
    # Non-time-mode: boxplots or barplots
    # -------------------------
    # group_var exists and is required here
    gv <- opt$group_var
    fv <- opt$facet_var

    if (summary_stat == "median") {
      # boxplot + jitter per group
      aes_args <- list(
        x = as.name(gv),
        y = quote(Abundance),
        colour = as.name(gv)
      )

      p <- ggplot(df_tax, do.call(aes, aes_args)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
        labs(
          title = paste0("Relative abundance by ", gv, ": ", tx),
          x = gv,
          colour = gv
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )

      if (!is.null(fv)) {
        p <- p + facet_wrap(as.formula(paste("~", fv)))
      }

      p <- apply_y_scale(p, df_tax)
      print(p)

    } else {
      # summary_stat == "mean": bar + whiskers + raw points (jitter)
      group_vars <- gv
      if (!is.null(fv)) {
        group_vars <- c(group_vars, fv)
      }

      df_sum <- df_tax %>%
        group_by(across(all_of(group_vars))) %>%
        summarise(
          n = n(),
          mean_abund = mean(Abundance, na.rm = TRUE),
          sd_abund = sd(Abundance, na.rm = TRUE),
          se_abund = sd_abund / sqrt(n),
          .groups = "drop"
        ) %>%
        mutate(
          y    = mean_abund,
          ymin = mean_abund - se_abund,
          ymax = mean_abund + se_abund
        )

      y_label <- "Mean relative abundance"

      aes_args <- list(
        x   = as.name(gv),
        y   = quote(y),
        fill = as.name(gv)
      )

      p <- ggplot(df_sum, do.call(aes, aes_args)) +
        geom_col(position = "dodge") +
        geom_errorbar(
          aes(ymin = ymin, ymax = ymax),
          width = 0.2,
          position = position_dodge(width = 0.9)
        ) +
        # add raw data points as jittered dots
        geom_jitter(
          data = df_tax,
          aes_string(x = gv, y = "Abundance"),
          width = 0.2,
          alpha = 0.4,
          size = 1,
          inherit.aes = FALSE
        ) +
        labs(
          title = paste0("Summary (mean) by ", gv, ": ", tx),
          x = gv,
          y = y_label,
          fill = gv
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )

      if (!is.null(fv)) {
        p <- p + facet_wrap(as.formula(paste("~", fv)))
      }

      # For summary plots, we still respect free_y vs fixed
      p <- apply_y_scale(p, df_tax)
      print(p)
    }
  }
}

message("Finished. Multi-page PDF written to: ", normalizePath(opt$out))
