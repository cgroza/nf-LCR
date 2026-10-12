#!/usr/bin/env Rscript

paf_hor <- "strand_normalized/hor_to_haplotypes.sensitive.paf"
gene_mappings_path <- "strand_normalized/full_chr22q_lcr_gene_layout_gene_mappings.tsv"
gene_reference_path <- "strand_normalized/chm13_chr22q_protein_coding_genes.tsv"
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

read_windows <- function(lcr) {
  read.delim(
    file.path(windows_dir, paste0("lcr", tolower(sub("LCR-", "", lcr)), "_display_windows.tsv")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

all_haplotype_order <- order_haplotypes(read.delim(
  file.path(windows_dir, "orientations.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)$haplotype)

plot_lcr <- function(lcr) {
  output_prefix <- file.path(
    windows_dir,
    paste0("lcr", tolower(sub("LCR-", "", lcr)), "_hor_gene_configurations_strand_normalized")
  )
  windows <- read_windows(lcr)
  haplotype_order <- all_haplotype_order
  row_labels <- vapply(haplotype_order, short_name, character(1))

  paf <- read.delim(paf_hor, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  hor_hits <- data.frame(
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
  hor_hits <- merge(hor_hits, windows, by = "haplotype", all = FALSE)
  hor_hits <- hor_hits[
    hor_hits$target_end_0based > hor_hits$lcr_start_0based &
      hor_hits$target_start_0based < hor_hits$lcr_end_0based,
  ]
  hor_hits$query_coverage <- (hor_hits$query_end_0based - hor_hits$query_start_0based) /
    hor_hits$query_length_bp
  hor_hits$identity <- hor_hits$matches / hor_hits$alignment_length_bp
  hor_hits$repeat_label <- sub("kbp$", " kb", sub("^.*consensus_", "", hor_hits$repeat_name))
  hor_hits <- remove_redundant_105_fragments(hor_hits)

  gene_hits <- read.delim(gene_mappings_path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_hits <- merge(gene_hits, windows, by = "haplotype", all = FALSE)
  gene_hits <- gene_hits[
    gene_hits$target_end_0based > gene_hits$lcr_start_0based &
      gene_hits$target_start_0based < gene_hits$lcr_end_0based,
  ]

  if (nrow(hor_hits) > 0) {
    anchors <- aggregate(target_start_0based ~ haplotype, data = hor_hits, FUN = min)
    names(anchors)[2] <- "anchor_0based"
    anchors <- merge(windows, anchors, by = "haplotype", all.x = TRUE)
    anchors$anchor_0based[is.na(anchors$anchor_0based)] <-
      anchors$lcr_start_0based[is.na(anchors$anchor_0based)]
  } else {
    anchors <- windows[c("haplotype", "lcr_start_0based")]
    names(anchors)[2] <- "anchor_0based"
  }

  hor_hits <- merge(hor_hits, anchors[c("haplotype", "anchor_0based")], by = "haplotype")
  gene_hits <- merge(gene_hits, anchors[c("haplotype", "anchor_0based")], by = "haplotype")
  hor_hits$row_index <- match(hor_hits$haplotype, haplotype_order)
  hor_hits$row_y <- length(haplotype_order) - hor_hits$row_index + 1
  gene_hits$row_index <- match(gene_hits$haplotype, haplotype_order)
  gene_hits$row_y <- length(haplotype_order) - gene_hits$row_index + 1
  hor_hits$start_kb <- (hor_hits$target_start_0based - hor_hits$anchor_0based) / 1e3
  hor_hits$end_kb <- (hor_hits$target_end_0based - hor_hits$anchor_0based) / 1e3
  gene_hits$start_kb <- (gene_hits$target_start_0based - gene_hits$anchor_0based) / 1e3
  gene_hits$end_kb <- (gene_hits$target_end_0based - gene_hits$anchor_0based) / 1e3
  gene_hits$midpoint_kb <- (gene_hits$start_kb + gene_hits$end_kb) / 2
  hor_hits$opacity <- pmax(0.25, pmin(1, hor_hits$query_coverage * hor_hits$identity))
  hor_hits$opacity[hor_hits$repeat_label == "41 kb"] <- 0.35

  genes <- read.delim(gene_reference_path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_colors <- setNames(
    grDevices::hcl(h = seq(15, 375, length.out = nrow(genes) + 1)[-1], c = 85, l = 55),
    genes$gene
  )
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
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "#F5F5F5",
      colour = NA
    ) +
    geom_hline(yintercept = c(8.5, 6.5, 4.5, 2.5), colour = "#BDBDBD", linewidth = 0.35) +
    geom_segment(
      data = gene_hits,
      aes(x = start_kb, xend = end_kb, y = row_y + 0.38, yend = row_y + 0.38, colour = gene),
      linewidth = 1.2,
      lineend = "butt"
    ) +
    geom_text(
      data = gene_hits,
      aes(x = midpoint_kb, y = row_y + 0.53, label = gene, colour = gene),
      size = 2.1,
      vjust = 0,
      check_overlap = TRUE
    )

  if (nrow(hor_hits) > 0) {
    triangles <- do.call(rbind, lapply(seq_len(nrow(hor_hits)), function(i) {
      hit <- hor_hits[i, ]
      if (hit$strand == "+") {
        data.frame(x = c(hit$start_kb, hit$start_kb, hit$end_kb), y = hit$row_y + c(-0.28, 0.28, 0))
      } else {
        data.frame(x = c(hit$end_kb, hit$end_kb, hit$start_kb), y = hit$row_y + c(-0.28, 0.28, 0))
      }
    }))
    triangles$hit_id <- rep(seq_len(nrow(hor_hits)), each = 3)
    triangles$repeat_label <- rep(hor_hits$repeat_label, each = 3)
    triangles$opacity <- rep(hor_hits$opacity, each = 3)
    triangles$draw_order <- match(
      triangles$repeat_label,
      c("41 kb", "105 kb", "120 kb", "28 kb", "25 kb")
    )
    triangles <- triangles[order(triangles$draw_order, triangles$hit_id), ]
    for (repeat_label in c("41 kb", "105 kb", "120 kb", "28 kb", "25 kb")) {
      plot <- plot + geom_polygon(
        data = triangles[triangles$repeat_label == repeat_label, ],
        aes(x = x, y = y, group = hit_id, fill = repeat_label, colour = repeat_label, alpha = opacity),
        linewidth = 0.3
      )
    }
    plot <- plot +
      scale_fill_manual(values = repeat_colors, name = "HOR consensus") +
      scale_alpha_identity(guide = "none")
  }

  plot <- plot +
    scale_colour_manual(values = c(repeat_colors, gene_colors), guide = "none") +
    scale_x_continuous(name = paste("Position relative to first mapped", lcr, "HOR (kb)")) +
    scale_y_continuous(
      name = NULL,
      breaks = rev(seq_along(row_labels)),
      labels = row_labels,
      limits = c(0.5, 9.65),
      expand = c(0, 0)
    ) +
    labs(
      title = paste("chr22q11.2", lcr, "HOR configurations with mapped protein-coding genes")
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )

  write.table(
    gene_hits[order(gene_hits$row_index, gene_hits$target_start_0based), ],
    paste0(output_prefix, "_gene_mappings.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  ggsave(paste0(output_prefix, ".pdf"), plot, width = 16, height = 8, units = "in", useDingbats = FALSE)
  ggsave(paste0(output_prefix, ".png"), plot, width = 16, height = 8, units = "in", dpi = 300)
}

for (lcr in c("LCR-A", "LCR-B", "LCR-C", "LCR-D")) plot_lcr(lcr)
