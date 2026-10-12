#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
paf_path <- if (length(args) >= 1) args[[1]] else "hor_annotation/hor_to_haplotypes.sensitive.paf"
windows_dir <- if (length(args) >= 2) args[[2]] else "hor_annotation"
output_prefix <- if (length(args) >= 3) args[[3]] else "hor_annotation/full_chr22q_lcr_layout"

suppressPackageStartupMessages(library(ggplot2))

read_lcr_hits <- function(lcr) {
  windows <- read.delim(
    file.path(windows_dir, paste0("lcr", tolower(sub("LCR-", "", lcr)), "_display_windows.tsv")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  paf <- read.delim(paf_path, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
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
  hits <- merge(hits, windows, by = "haplotype", all = FALSE)
  hits <- hits[
    hits$target_end_0based > hits$lcr_start_0based &
      hits$target_start_0based < hits$lcr_end_0based,
  ]
  hits$lcr <- rep(lcr, nrow(hits))
  hits
}

read_lcr_boxes <- function(lcr) {
  boxes <- read.delim(
    file.path(windows_dir, paste0("lcr", tolower(sub("LCR-", "", lcr)), "_display_windows.tsv")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  boxes$lcr <- lcr
  boxes
}

short_name <- function(name) {
  if (startsWith(name, "CHM13#")) {
    return("CHM13 control")
  }
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

haplotype_order <- order_haplotypes(read.delim(
  file.path(windows_dir, "orientations.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)$haplotype)

lcr_names <- c("LCR-A", "LCR-B", "LCR-C", "LCR-D")
hits <- do.call(rbind, lapply(lcr_names, read_lcr_hits))
hits$query_coverage <- (hits$query_end_0based - hits$query_start_0based) / hits$query_length_bp
hits$identity <- hits$matches / hits$alignment_length_bp
hits$repeat_label <- sub("^.*consensus_", "", hits$repeat_name)
hits$repeat_label <- sub("kbp$", " kb", hits$repeat_label)
hits <- remove_redundant_105_fragments(hits)

row_labels <- vapply(haplotype_order, short_name, character(1))
hits$row_index <- match(hits$haplotype, haplotype_order)
hits$row_y <- length(haplotype_order) - hits$row_index + 1

hits$start_mb <- 18 + hits$target_start_0based / 1e6
hits$end_mb <- 18 + hits$target_end_0based / 1e6
hits$opacity <- pmax(0.25, pmin(1, hits$query_coverage * hits$identity))
hits$opacity[hits$repeat_label == "41 kb"] <- 0.35

lcr_boxes <- do.call(rbind, lapply(lcr_names, read_lcr_boxes))
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
control_lcr_boxes <- bounded_lcr_boxes[bounded_lcr_boxes$haplotype == haplotype_order[[1]], ]

triangle_data <- do.call(
  rbind,
  lapply(seq_len(nrow(hits)), function(i) {
    hit <- hits[i, ]
    if (hit$strand == "+") {
      data.frame(x = c(hit$start_mb, hit$start_mb, hit$end_mb), y = hit$row_y + c(-0.36, 0.36, 0))
    } else {
      data.frame(x = c(hit$end_mb, hit$end_mb, hit$start_mb), y = hit$row_y + c(-0.36, 0.36, 0))
    }
  })
)
triangle_data$hit_id <- rep(seq_len(nrow(hits)), each = 3)
triangle_data$repeat_label <- rep(hits$repeat_label, each = 3)
triangle_data$opacity <- rep(hits$opacity, each = 3)
triangle_data$draw_order <- match(
  triangle_data$repeat_label,
  c("41 kb", "105 kb", "120 kb", "28 kb", "25 kb")
)
triangle_data <- triangle_data[order(triangle_data$draw_order, triangle_data$hit_id), ]

repeat_colors <- c(
  "25 kb" = "#3B6EA8",
  "28 kb" = "#E68613",
  "41 kb" = "#4C9F38",
  "105 kb" = "#CF4446",
  "120 kb" = "#8F5B8D"
)
genome_bands <- data.frame(ymin = c(6.5, 2.5), ymax = c(8.5, 4.5))

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
    aes(x = (start_mb + end_mb) / 2, y = 9.35, label = lcr),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  geom_hline(yintercept = c(8.5, 6.5, 4.5, 2.5), colour = "#BDBDBD", linewidth = 0.35) +
  scale_fill_manual(values = repeat_colors, name = "HOR consensus") +
  scale_colour_manual(
    values = c(
      repeat_colors,
      "LCR-A" = "#6A3D9A",
      "LCR-B" = "#1B9E77",
      "LCR-C" = "#C58A00",
      "LCR-D" = "#D95F02"
    ),
    guide = "none"
  ) +
  scale_alpha_identity(guide = "none") +
  scale_x_continuous(
    name = "CHM13 chr22 coordinate (Mb)",
    limits = c(18, 23),
    breaks = 18:23
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = rev(seq_along(row_labels)),
    labels = row_labels,
    limits = c(0.5, 9.5),
    expand = c(0, 0)
  ) +
  labs(
    title = "chr22q11.2 LCR-A through LCR-D HOR configurations"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

for (repeat_label in c("41 kb", "105 kb", "120 kb", "28 kb", "25 kb")) {
  plot <- plot + geom_polygon(
    data = triangle_data[triangle_data$repeat_label == repeat_label, ],
    aes(x = x, y = y, group = hit_id, fill = repeat_label, colour = repeat_label, alpha = opacity),
    linewidth = 0.3
  )
}

write.table(
  hits[order(hits$row_index, hits$lcr, hits$target_start_0based), ],
  paste0(output_prefix, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
ggsave(paste0(output_prefix, ".pdf"), plot, width = 16, height = 8, units = "in", useDingbats = FALSE)
ggsave(paste0(output_prefix, ".png"), plot, width = 16, height = 8, units = "in", dpi = 300)
