#!/usr/bin/env Rscript

paf_hor <- "strand_normalized/hor_to_haplotypes.sensitive.paf"
paf_genes <- "strand_normalized/chm13_protein_coding_genes_to_haplotypes.paf"
genes_path <- "strand_normalized/chm13_chr22q_protein_coding_genes.tsv"
output_prefix <- "strand_normalized/full_chr22q_lcr_gene_layout"
windows_dir <- "strand_normalized"

suppressPackageStartupMessages(library(ggplot2))

short_name <- function(name) {
  if (startsWith(name, "CHM13#")) return("CHM13 control")
  name <- sub("^cmh", "CMH ", name)
  name <- gsub("#", " | ", name, fixed = TRUE)
  name <- sub(" \\| h[12]tg.*$", "", name)
  sub(" \\| hap", " | hap", name)
}

order_haplotypes <- function(haplotypes) {
  sample <- sub("^.*[-_]([0-9]+)(?:#|_)hap[12].*$", "\\1", haplotypes)
  haplotype <- suppressWarnings(as.integer(sub("^.*(?:#|_)hap([12]).*$", "\\1", haplotypes)))
  is_control <- startsWith(haplotypes, "CHM13#")
  sample_rank <- suppressWarnings(as.integer(sample))
  sample_rank[sample == "02"] <- -2L
  sample_rank[sample == "03"] <- -1L
  sample_rank[is.na(sample_rank)] <- .Machine$integer.max
  haplotype[is.na(haplotype)] <- .Machine$integer.max
  haplotypes[order(!is_control, sample_rank, haplotype, haplotypes)]
}

remove_redundant_105_fragments <- function(hits) {
  keep <- rep(TRUE, nrow(hits))
  for (index in which(hits$repeat_label == "105 kb")) {
    candidates <- hits[
      hits$haplotype == hits$haplotype[index] & hits$repeat_label == "41 kb",
    ]
    if (!nrow(candidates)) next
    overlap <- pmax(
      0,
      pmin(hits$target_end_0based[index], candidates$target_end_0based) -
        pmax(hits$target_start_0based[index], candidates$target_start_0based)
    )
    if (any(
      overlap / (hits$target_end_0based[index] - hits$target_start_0based[index]) >= 0.9 &
        overlap / (candidates$target_end_0based - candidates$target_start_0based) >= 0.9
    )) keep[index] <- FALSE
  }
  hits[keep, , drop = FALSE]
}

read_lcr_windows <- function(lcr) {
  read.delim(
    file.path(windows_dir, paste0("lcr", tolower(sub("LCR-", "", lcr)), "_display_windows.tsv")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

read_lcr_boxes <- function(lcr) {
  boxes <- read_lcr_windows(lcr)
  boxes$lcr <- lcr
  boxes
}

haplotype_order <- order_haplotypes(read.delim(
  file.path(windows_dir, "orientations.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)$haplotype)

read_hor_hits <- function() {
  paf <- read.delim(paf_hor, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  hits <- data.frame(
    repeat_name = paf[[1]],
    query_length_bp = paf[[2]],
    query_start_0based = paf[[3]],
    query_end_0based = paf[[4]],
    strand = paf[[5]],
    haplotype = paf[[6]],
    target_start_0based = paf[[8]],
    target_end_0based = paf[[9]],
    matches = paf[[10]],
    alignment_length_bp = paf[[11]],
    stringsAsFactors = FALSE
  )
  do.call(rbind, lapply(c("LCR-A", "LCR-B", "LCR-C", "LCR-D"), function(lcr) {
    windows <- read_lcr_windows(lcr)
    lcr_hits <- merge(hits, windows, by = "haplotype", all = FALSE)
    lcr_hits <- lcr_hits[
      lcr_hits$target_end_0based > lcr_hits$lcr_start_0based &
        lcr_hits$target_start_0based < lcr_hits$lcr_end_0based,
    ]
    lcr_hits$lcr <- rep(lcr, nrow(lcr_hits))
    lcr_hits
  }))
}

assign_label_lanes <- function(genes, minimum_spacing_bp = 80000L) {
  genes <- genes[order(genes$roi_start_0based, genes$roi_end_0based), ]
  lane_ends <- integer()
  lanes <- integer(nrow(genes))
  for (i in seq_len(nrow(genes))) {
    available <- which(lane_ends < genes$roi_start_0based[i])
    lane <- if (length(available)) available[[1]] else length(lane_ends) + 1L
    if (lane > length(lane_ends)) lane_ends[lane] <- 0L
    lane_ends[lane] <- genes$roi_end_0based[i] + minimum_spacing_bp
    lanes[i] <- lane
  }
  genes$label_lane <- lanes
  genes
}

genes <- read.delim(genes_path, check.names = FALSE, stringsAsFactors = FALSE)
genes$reference_start_mb <- 18 + genes$roi_start_0based / 1e6
genes$reference_end_mb <- 18 + genes$roi_end_0based / 1e6
genes$display_start_mb <- pmax(18, genes$reference_start_mb)
genes$display_end_mb <- pmin(23, genes$reference_end_mb)
genes <- assign_label_lanes(genes)

gene_paf <- read.delim(paf_genes, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
gene_hits <- data.frame(
  gene = sub("\\|.*$", "", gene_paf[[1]]),
  query_length_bp = gene_paf[[2]],
  query_start_0based = gene_paf[[3]],
  query_end_0based = gene_paf[[4]],
  mapping_strand = gene_paf[[5]],
  haplotype = gene_paf[[6]],
  target_start_0based = gene_paf[[8]],
  target_end_0based = gene_paf[[9]],
  matches = gene_paf[[10]],
  alignment_length_bp = gene_paf[[11]],
  mapq = gene_paf[[12]],
  stringsAsFactors = FALSE
)
gene_hits$query_coverage <- (gene_hits$query_end_0based - gene_hits$query_start_0based) /
  gene_hits$query_length_bp
gene_hits$identity <- gene_hits$matches / gene_hits$alignment_length_bp
gene_hits$score <- gene_hits$query_coverage * gene_hits$identity
gene_hits <- gene_hits[order(gene_hits$gene, gene_hits$haplotype, -gene_hits$score), ]
gene_hits <- gene_hits[!duplicated(gene_hits[c("gene", "haplotype")]), ]
gene_hits <- gene_hits[gene_hits$query_coverage >= 0.8 & gene_hits$identity >= 0.9, ]
gene_hits <- merge(gene_hits, genes, by = "gene", all.x = TRUE)

row_labels <- vapply(haplotype_order, short_name, character(1))
gene_hits$row_index <- match(gene_hits$haplotype, haplotype_order)
gene_hits$row_y <- length(haplotype_order) - gene_hits$row_index + 1
gene_hits$start_mb <- 18 + gene_hits$target_start_0based / 1e6
gene_hits$end_mb <- 18 + gene_hits$target_end_0based / 1e6
gene_hits$display_start_mb <- pmax(18, gene_hits$start_mb)
gene_hits$display_end_mb <- pmin(23, gene_hits$end_mb)
gene_hits_plot <- gene_hits[
  gene_hits$end_mb > 18 & gene_hits$start_mb < 23,
]

hor_hits <- read_hor_hits()
hor_hits$query_coverage <- (hor_hits$query_end_0based - hor_hits$query_start_0based) /
  hor_hits$query_length_bp
hor_hits$identity <- hor_hits$matches / hor_hits$alignment_length_bp
hor_hits$repeat_label <- sub("kbp$", " kb", sub("^.*consensus_", "", hor_hits$repeat_name))
hor_hits <- remove_redundant_105_fragments(hor_hits)
hor_hits$row_index <- match(hor_hits$haplotype, haplotype_order)
hor_hits$row_y <- length(haplotype_order) - hor_hits$row_index + 1
hor_hits$start_mb <- 18 + hor_hits$target_start_0based / 1e6
hor_hits$end_mb <- 18 + hor_hits$target_end_0based / 1e6
hor_hits$opacity <- pmax(0.25, pmin(1, hor_hits$query_coverage * hor_hits$identity))
hor_hits$opacity[hor_hits$repeat_label == "41 kb"] <- 0.35

lcr_boxes <- do.call(
  rbind,
  lapply(c("LCR-A", "LCR-B", "LCR-C", "LCR-D"), read_lcr_boxes)
)
boundary_sources <- read.delim(
  file.path(windows_dir, "lcr_boundary_sources.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
lcr_boxes <- merge(lcr_boxes, boundary_sources, by = c("haplotype", "lcr"), all.x = TRUE)
lcr_boxes$row_index <- match(lcr_boxes$haplotype, haplotype_order)
lcr_boxes$row_y <- length(haplotype_order) - lcr_boxes$row_index + 1
lcr_boxes$start_mb <- 18 + lcr_boxes$lcr_start_0based / 1e6
lcr_boxes$end_mb <- 18 + lcr_boxes$lcr_end_0based / 1e6
lcr_boxes$bounded <- lcr_boxes$boundary_source %in% c("both_flanks", "reference")
bounded_lcr_boxes <- lcr_boxes[lcr_boxes$bounded, ]
partial_lcr_calls <- lcr_boxes[!lcr_boxes$bounded, ]
partial_lcr_calls$known_edge_mb <- ifelse(
  partial_lcr_calls$boundary_source == "left_flank_plus_chm13_length",
  partial_lcr_calls$start_mb,
  partial_lcr_calls$end_mb
)
partial_lcr_calls$unknown_direction <- ifelse(
  partial_lcr_calls$boundary_source == "left_flank_plus_chm13_length",
  1,
  -1
)
partial_lcr_calls$open_end_mb <- partial_lcr_calls$known_edge_mb +
  partial_lcr_calls$unknown_direction * 0.10

triangle_data <- do.call(rbind, lapply(seq_len(nrow(hor_hits)), function(i) {
  hit <- hor_hits[i, ]
  if (hit$strand == "+") {
    data.frame(x = c(hit$start_mb, hit$start_mb, hit$end_mb), y = hit$row_y + c(-0.27, 0.27, 0))
  } else {
    data.frame(x = c(hit$end_mb, hit$end_mb, hit$start_mb), y = hit$row_y + c(-0.27, 0.27, 0))
  }
}))
triangle_data$hit_id <- rep(seq_len(nrow(hor_hits)), each = 3)
triangle_data$repeat_label <- rep(hor_hits$repeat_label, each = 3)
triangle_data$opacity <- rep(hor_hits$opacity, each = 3)
triangle_data$draw_order <- match(
  triangle_data$repeat_label,
  c("41 kb", "105 kb", "120 kb", "28 kb", "25 kb")
)
triangle_data <- triangle_data[order(triangle_data$draw_order, triangle_data$hit_id), ]

control_lcr_boxes <- bounded_lcr_boxes[bounded_lcr_boxes$haplotype == haplotype_order[[1]], ]
repeat_colors <- c(
  "25 kb" = "#3B6EA8",
  "28 kb" = "#E68613",
  "41 kb" = "#4C9F38",
  "105 kb" = "#CF4446",
  "120 kb" = "#8F5B8D"
)
lcr_border_colors <- c(
  "LCR-A" = "#6A3D9A",
  "LCR-B" = "#1B9E77",
  "LCR-C" = "#C58A00",
  "LCR-D" = "#D95F02"
)
gene_colors <- setNames(
  grDevices::hcl(
    h = seq(15, 375, length.out = nrow(genes) + 1)[-1],
    c = 85,
    l = 55
  ),
  genes$gene
)
genome_bands <- data.frame(ymin = c(6.5, 2.5), ymax = c(8.5, 4.5))
label_y <- 10.2 + (genes$label_lane - 1) * 0.42
max_y <- max(label_y) + 0.35

plot <- ggplot() +
  geom_rect(
    data = genome_bands,
    aes(xmin = 18, xmax = 23, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    colour = NA,
    fill = "#F5F5F5"
  ) +
  geom_rect(
    data = bounded_lcr_boxes,
    aes(
      xmin = start_mb,
      xmax = end_mb,
      ymin = row_y - 0.48,
      ymax = row_y + 0.48,
      colour = lcr
    ),
    inherit.aes = FALSE,
    fill = "#E0E0E0",
    alpha = 0.55
  ) +
  geom_segment(
    data = partial_lcr_calls,
    aes(
      x = known_edge_mb,
      xend = open_end_mb,
      y = row_y,
      yend = row_y,
      colour = lcr
    ),
    inherit.aes = FALSE,
    linewidth = 1.3,
    linetype = "dashed",
    arrow = grid::arrow(length = grid::unit(0.08, "inches"), type = "open")
  ) +
  geom_segment(
    data = partial_lcr_calls,
    aes(
      x = known_edge_mb,
      xend = known_edge_mb,
      y = row_y - 0.38,
      yend = row_y + 0.38,
      colour = lcr
    ),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_text(
    data = control_lcr_boxes,
    aes(x = (start_mb + end_mb) / 2, y = 9.45, label = lcr),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  geom_segment(
    data = gene_hits_plot,
    aes(
      x = display_start_mb,
      xend = display_end_mb,
      y = row_y + 0.36,
      yend = row_y + 0.36,
      colour = gene
    ),
    linewidth = 1.1,
    lineend = "butt"
  ) +
  geom_segment(
    data = genes,
    aes(
      x = (display_start_mb + display_end_mb) / 2,
      xend = (display_start_mb + display_end_mb) / 2,
      y = 9.62,
      yend = label_y,
      colour = gene
    ),
    inherit.aes = FALSE,
    linewidth = 0.2
  ) +
  geom_text(
    data = genes,
    aes(
      x = (display_start_mb + display_end_mb) / 2,
      y = label_y,
      label = gene,
      colour = gene
    ),
    inherit.aes = FALSE,
    size = 2.2,
    vjust = 0,
    check_overlap = FALSE
  ) +
  geom_hline(yintercept = c(8.5, 6.5, 4.5, 2.5), colour = "#BDBDBD", linewidth = 0.35) +
  scale_fill_manual(values = repeat_colors, name = "HOR consensus") +
  scale_colour_manual(values = c(repeat_colors, lcr_border_colors, gene_colors), guide = "none") +
  scale_alpha_identity(guide = "none") +
  scale_x_continuous(name = "CHM13 chr22 coordinate (Mb)", limits = c(18, 23), breaks = 18:23) +
  scale_y_continuous(
    name = NULL,
    breaks = rev(seq_along(row_labels)),
    labels = row_labels,
    limits = c(0.5, max_y),
    expand = c(0, 0)
  ) +
  labs(
    title = "chr22q11.2 LCR-A through LCR-D configurations with mapped protein-coding genes"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )

for (repeat_label in c("41 kb", "105 kb", "120 kb", "28 kb", "25 kb")) {
  plot <- plot + geom_polygon(
    data = triangle_data[triangle_data$repeat_label == repeat_label, ],
    aes(x = x, y = y, group = hit_id, fill = repeat_label, colour = repeat_label, alpha = opacity),
    linewidth = 0.3
  )
}

write.table(
  gene_hits[order(gene_hits$row_index, gene_hits$target_start_0based), ],
  paste0(output_prefix, "_gene_mappings.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
ggsave(paste0(output_prefix, ".pdf"), plot, width = 20, height = 13, units = "in", useDingbats = FALSE)
ggsave(paste0(output_prefix, ".png"), plot, width = 20, height = 13, units = "in", dpi = 300)
