# ============================================================================
# Lastovo miniCOI — community composition, DEICODE rPCA and PERMANOVA summary
# ============================================================================

required <- c("ggplot2", "dplyr", "tidyr", "vegan", "pheatmap", "ggrepel", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(ggrepel)
  library(patchwork)
})

args0 <- commandArgs(trailingOnly = FALSE)
filearg <- sub("^--file=", "", args0[grep("^--file=", args0)])
ROOT <- if (length(filearg)) {
  normalizePath(file.path(dirname(filearg[1]), ".."))
} else {
  normalizePath(getwd())
}

EXPORT <- file.path(ROOT, "exports")
OUT <- file.path(ROOT, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

read_biom_tsv <- function(path) {
  x <- read.delim(
    path, sep = "\t", header = TRUE, skip = 1,
    check.names = FALSE, quote = "", comment.char = "",
    stringsAsFactors = FALSE
  )
  names(x)[1] <- "feature_id"
  rownames(x) <- x$feature_id
  x$feature_id <- NULL
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  m
}

read_metadata <- function(path) {
  x <- read.delim(
    path, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, quote = "", comment.char = ""
  )
  names(x)[1] <- "sample_id"
  x$ITW_label <- ifelse(tolower(x$ITW) == "itw", "ITW", "noITW")
  x$mesh_label <- ifelse(x$fraction == "micro", "53 µm", "200 µm")
  x$layer_label <- ifelse(
    x$layer == "upper", "Upper",
    ifelse(x$layer == "lower", "Lower", "Whole")
  )
  x
}

# MetaZooGene taxonomy exported from this project uses a semicolon-delimited
# hierarchy in which genus/species occur at the same positions used in the
# original final-analysis script.
tax_level <- function(taxonomy, position) {
  split_tax <- strsplit(as.character(taxonomy), ";", fixed = TRUE)
  vapply(split_tax, function(v) {
    v <- trimws(v)
    if (length(v) < position) return(NA_character_)
    y <- sub("^[A-Za-z]__", "", v[position])
    if (!nzchar(y)) NA_character_ else y
  }, character(1))
}

parse_ordination <- function(path) {
  z <- readLines(path, warn = FALSE)

  parse_block <- function(label) {
    idx <- grep(paste0("^", label, "\t"), z)
    if (!length(idx)) return(NULL)
    header <- strsplit(z[idx[1]], "\t", fixed = TRUE)[[1]]
    n_rows <- as.integer(header[2])
    n_axes <- as.integer(header[3])
    block_lines <- z[(idx[1] + 1):(idx[1] + n_rows)]
    out <- read.delim(
      text = paste(block_lines, collapse = "\n"),
      sep = "\t", header = FALSE,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    names(out) <- c("ID", paste0("Axis", seq_len(n_axes)))
    out
  }

  prop_idx <- grep("^Proportion explained\t", z)
  prop <- c(NA_real_, NA_real_, NA_real_)
  if (length(prop_idx)) {
    vals <- as.numeric(strsplit(z[prop_idx[1] + 1], "\t", fixed = TRUE)[[1]])
    prop[seq_len(min(3, length(vals)))] <- vals[seq_len(min(3, length(vals)))]
  }

  feature_scores <- parse_block("Species")
  if (is.null(feature_scores)) feature_scores <- parse_block("Biplot")

  list(site = parse_block("Site"), species = feature_scores, proportion = prop)
}

read_permanova_qzv <- function(path) {
  if (!file.exists(path)) stop("Missing PERMANOVA qzv: ", path)
  listing <- utils::unzip(path, list = TRUE)
  index_path <- listing$Name[grepl("/data/index\\.html$", listing$Name)]
  if (length(index_path) != 1) stop("Could not locate data/index.html in: ", path)

  con <- unz(path, index_path)
  html <- paste(readLines(con, warn = FALSE, encoding = "UTF-8"), collapse = " ")
  close(con)
  html <- gsub("[[:space:]]+", " ", html)

  extract_row <- function(label) {
    pattern <- paste0("<th>", label, "</th>[[:space:]]*<td>([^<]+)</td>")
    hit <- regexec(pattern, html, perl = TRUE)
    value <- regmatches(html, hit)[[1]]
    if (length(value) < 2) return(NA_character_)
    trimws(value[2])
  }

  data.frame(
    sample_size = as.numeric(extract_row("sample size")),
    groups = as.numeric(extract_row("number of groups")),
    F = as.numeric(extract_row("test statistic")),
    p = as.numeric(extract_row("p-value")),
    permutations = as.numeric(extract_row("number of permutations"))
  )
}

meta <- read_metadata(file.path(EXPORT, "COI_metadata_original.tsv"))
meta19 <- meta[meta$layer %in% c("upper", "lower"), , drop = FALSE]

mat20 <- read_biom_tsv(file.path(EXPORT, "conf90_all20", "feature-table.tsv"))
tax <- read.delim(
  file.path(EXPORT, "taxonomy", "taxonomy.tsv"),
  sep = "\t", header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "", comment.char = ""
)
names(tax)[1] <- "feature_id"
tax$genus <- tax_level(tax$Taxon, 18)

common <- intersect(rownames(mat20), tax$feature_id)
mat20 <- mat20[common, , drop = FALSE]
tax_mat <- tax[match(rownames(mat20), tax$feature_id), , drop = FALSE]

valid_genus <- !is.na(tax_mat$genus) & tax_mat$genus != ""
genus_mat <- rowsum(
  mat20[valid_genus, , drop = FALSE],
  group = tax_mat$genus[valid_genus]
)
storage.mode(genus_mat) <- "numeric"

rel <- sweep(genus_mat, 2, colSums(genus_mat), "/")
taxon_order <- order(rowSums(genus_mat), decreasing = TRUE)
top10 <- rownames(genus_mat)[taxon_order[seq_len(min(10, nrow(genus_mat)))]]
top20 <- rownames(genus_mat)[taxon_order[seq_len(min(20, nrow(genus_mat)))]]

# ---------------------------------------------------------------------------
# Genus-level relative read composition
# ---------------------------------------------------------------------------
bar_long <- as.data.frame(t(rel[top10, , drop = FALSE]))
bar_long$sample_id <- rownames(bar_long)
bar_long <- bar_long %>%
  pivot_longer(cols = -sample_id, names_to = "genus", values_to = "relative") %>%
  left_join(
    meta[, c("sample_id", "ITW_label", "mesh_label", "layer_label")],
    by = "sample_id"
  ) %>%
  filter(sample_id != "N12", layer_label %in% c("Upper", "Lower"))

bar_summary <- bar_long %>%
  group_by(mesh_label, layer_label, ITW_label, genus) %>%
  summarise(
    n = n(),
    mean_relative = mean(relative),
    se_relative = ifelse(n > 1, sd(relative) / sqrt(n), NA_real_),
    .groups = "drop"
  )
write.csv(bar_summary, file.path(OUT, "top10_genus_summary.csv"), row.names = FALSE)

p_bar <- ggplot(bar_summary, aes(genus, mean_relative, fill = ITW_label)) +
  geom_col(
    position = position_dodge(0.82), width = 0.74,
    colour = "grey20", linewidth = 0.2
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_relative - se_relative),
      ymax = mean_relative + se_relative
    ),
    position = position_dodge(0.82), width = 0.24, na.rm = TRUE
  ) +
  facet_grid(layer_label ~ mesh_label) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
  labs(
    x = NULL,
    y = "Mean relative read composition",
    fill = "ITW",
    title = "Dominant genus-level composition"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ggsave(file.path(OUT, "composition_top10.png"), p_bar, width = 13, height = 8, dpi = 300)
ggsave(file.path(OUT, "composition_top10.pdf"), p_bar, width = 13, height = 8)

# ---------------------------------------------------------------------------
# Bray-Curtis clustered top-20 heatmap
# ---------------------------------------------------------------------------
hm_rel <- rel[top20, intersect(meta19$sample_id, colnames(rel)), drop = FALSE]
hm_log <- log10(hm_rel + 1e-6)
annotation <- meta19[
  match(colnames(hm_log), meta19$sample_id),
  c("sample_id", "ITW_label", "mesh_label", "layer_label")
]
annotation_col <- data.frame(
  ITW = annotation$ITW_label,
  Mesh = annotation$mesh_label,
  Layer = annotation$layer_label,
  row.names = annotation$sample_id
)
row_cluster <- hclust(vegdist(hm_rel, method = "bray"), method = "average")
col_cluster <- hclust(
  vegdist(t(rel[, colnames(hm_rel), drop = FALSE]), method = "bray"),
  method = "average"
)

png(file.path(OUT, "heatmap_top20.png"), width = 3300, height = 2400, res = 300)
pheatmap::pheatmap(
  hm_log,
  cluster_rows = row_cluster,
  cluster_cols = col_cluster,
  annotation_col = annotation_col,
  border_color = NA,
  main = "Top 20 genus-level taxa; Bray-Curtis clustering"
)
dev.off()

pdf(file.path(OUT, "heatmap_top20.pdf"), width = 11, height = 8)
pheatmap::pheatmap(
  hm_log,
  cluster_rows = row_cluster,
  cluster_cols = col_cluster,
  annotation_col = annotation_col,
  border_color = NA,
  main = "Top 20 genus-level taxa; Bray-Curtis clustering"
)
dev.off()

# ---------------------------------------------------------------------------
# DEICODE rPCA with top eight unique genus loading vectors
# ---------------------------------------------------------------------------
ord <- parse_ordination(file.path(EXPORT, "ordination", "ordination.txt"))
if (is.null(ord$site) || is.null(ord$species)) stop("Could not parse DEICODE ordination export.")

rpca <- ord$site
names(rpca)[1] <- "sample_id"
rpca <- rpca %>%
  filter(sample_id %in% meta19$sample_id) %>%
  left_join(
    meta19[, c("sample_id", "ITW_label", "mesh_label", "layer_label")],
    by = "sample_id"
  )

sp <- ord$species
names(sp)[1] <- "feature_id"
sp$loading_strength <- sqrt(
  sp$Axis1^2 + sp$Axis2^2 + if ("Axis3" %in% names(sp)) sp$Axis3^2 else 0
)
sp <- sp[order(sp$loading_strength, decreasing = TRUE), , drop = FALSE]
sp <- merge(sp, tax[, c("feature_id", "genus")], by = "feature_id", all.x = TRUE, sort = FALSE)
sp$label <- ifelse(is.na(sp$genus) | sp$genus == "", sp$feature_id, sp$genus)
sp <- sp[!duplicated(sp$label), , drop = FALSE]
sp <- head(sp, 8)

xr <- diff(range(rpca$Axis1, na.rm = TRUE))
yr <- diff(range(rpca$Axis2, na.rm = TRUE))
arrow_scale <- min(
  xr / max(abs(sp$Axis1), na.rm = TRUE),
  yr / max(abs(sp$Axis2), na.rm = TRUE)
) * 0.32
sp$xend <- sp$Axis1 * arrow_scale
sp$yend <- sp$Axis2 * arrow_scale
write.csv(sp, file.path(OUT, "rpca_top8_loadings.csv"), row.names = FALSE)

xlab <- if (!is.na(ord$proportion[1])) {
  paste0("rPCA1 (", round(100 * ord$proportion[1], 1), "%)")
} else "rPCA1"
ylab <- if (!is.na(ord$proportion[2])) {
  paste0("rPCA2 (", round(100 * ord$proportion[2], 1), "%)")
} else "rPCA2"

make_panel <- function(grouping, title) {
  ggplot(rpca, aes(Axis1, Axis2, colour = .data[[grouping]])) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_segment(
      data = sp,
      aes(x = 0, y = 0, xend = xend, yend = yend),
      inherit.aes = FALSE, colour = "black",
      arrow = grid::arrow(length = grid::unit(0.10, "inches"))
    ) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = sample_id), show.legend = FALSE, size = 2.8) +
    ggrepel::geom_text_repel(
      data = sp,
      aes(x = xend, y = yend, label = label),
      inherit.aes = FALSE, colour = "black", size = 2.8
    ) +
    labs(x = xlab, y = ylab, colour = title, title = title) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}

p_rpca <- (
  make_panel("ITW_label", "ITW") |
  make_panel("mesh_label", "Mesh size") |
  make_panel("layer_label", "Layer")
) + plot_annotation(
  title = "DEICODE robust-Aitchison rPCA; top eight genus loading vectors"
)

ggsave(file.path(OUT, "rpca_current.png"), p_rpca, width = 15, height = 5.7, dpi = 300)
ggsave(file.path(OUT, "rpca_current.pdf"), p_rpca, width = 15, height = 5.7)

# ---------------------------------------------------------------------------
# PERMANOVA statistics + Holm adjustment across the three factors
# ---------------------------------------------------------------------------
specs <- data.frame(
  Factor = c("ITW occurrence", "Mesh size", "Layer"),
  file = c("permanova_itw.qzv", "permanova_mesh.qzv", "permanova_layer.qzv"),
  stringsAsFactors = FALSE
)

rows <- lapply(seq_len(nrow(specs)), function(i) {
  z <- read_permanova_qzv(file.path(EXPORT, "beta_stats", specs$file[i]))
  df_between <- z$groups - 1
  df_within <- z$sample_size - z$groups
  r2 <- (z$F * df_between) / (z$F * df_between + df_within)

  data.frame(
    Factor = specs$Factor[i],
    N = z$sample_size,
    F = z$F,
    R2 = r2,
    p_raw = z$p,
    Permutations = z$permutations
  )
})

beta <- do.call(rbind, rows)
beta$p_Holm <- p.adjust(beta$p_raw, method = "holm")
write.csv(beta, file.path(OUT, "beta_PERMANOVA_results.csv"), row.names = FALSE)

cat("Community composition, rPCA and beta statistics complete.\n")
