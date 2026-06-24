#!/usr/bin/env Rscript

# ==============================================================================
# Analysis script: cycling proximal tubule cells after kidney IRI
# ===============================================================================
# Purpose:
#   Reproduce cycling-status annotation, Ccnd1 summaries, cycling vs non-cycling
#   differential expression/FGSEA, and focused metabolic gene dot plots for a
#   murine kidney ischemia-reperfusion injury single-cell RNA-seq dataset.
#
# Input:
#   A processed Seurat object containing proximal tubule cells and the following
#   metadata columns:
#     - Phase: Seurat cell-cycle phase annotation; expected values include S/G2M
#     - exp.cond: Control, IRI_short, IRI_long
#     - exp.time: Control, IRI_short_1d, IRI_short_3d, IRI_short_14d,
#                 IRI_long_1d, IRI_long_3d, IRI_long_14d
#
# Data source:
#   Balzer et al., Nature Communications 2022; PMID: 35821371; GEO: GSE180420
#
# Notes:
#   - Replace INPUT_RDS with the local path to the processed Seurat object.
#   - No user-specific or institutional file paths are included.
#   - All outputs are written to OUTPUT_DIR.
# ===============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(data.table)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(viridis)
  library(ggnewscale)
  library(pheatmap)
  library(grid)
})

set.seed(123)

# -----------------------------
# User configuration
# -----------------------------
INPUT_RDS <- "path/to/processed_proximal_tubule_seurat_object.rds"
OUTPUT_DIR <- "G2M_publication_outputs"
ASSAY_USE <- "RNA"
PADJ_CUTOFF <- 0.05

# Set to TRUE to run the complete pathway-level FGSEA. This can be slow.
RUN_FGSEA <- TRUE

# Set to TRUE to run canonical metabolic program module scoring.
RUN_METABOLIC_MODULE_SCORING <- TRUE

# -----------------------------
# Helper functions
# -----------------------------
make_output_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  normalizePath(path, mustWork = FALSE)
}

check_required_metadata <- function(obj, required_cols) {
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  if (length(missing_cols) > 0) {
    stop("Missing required metadata column(s): ", paste(missing_cols, collapse = ", "))
  }
}

save_pdf <- function(plot, filename, width, height, output_dir = OUTPUT_DIR) {
  pdf(file.path(output_dir, filename), width = width, height = height)
  print(plot)
  dev.off()
}

clean_pathway_name <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^GOBP_", "", x)
  x <- gsub("^REACTOME_", "", x)
  x <- gsub("^KEGG_", "", x)
  x <- gsub("^WP_", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

assign_pathway_family <- function(x) {
  x_low <- tolower(x)
  dplyr::case_when(
    grepl("cell cycle|mitotic|m phase|dna replication|chromosome|chromatid|centromere|kinetochore|spindle|g2m|g2 m|e2f|cyclin|cdk|checkpoint|meiosis", x_low) ~
      "Cell cycle / DNA replication",
    grepl("metabolism|metabolic|glycolysis|gluconeogenesis|pyruvate|tca|citric acid|tricarboxylic|oxidative phosphorylation|respiratory chain|electron transport|atp synthesis|mitochond|oxphos|fatty acid|lipid|peroxisome|cholesterol|sterol|sphingolipid|ceramide|beta oxidation|acyl|triglyceride|amino acid|glutamine|glutathione|purine|pyrimidine|nucleotide|redox|reactive oxygen|ros", x_low) ~
      "Metabolism",
    grepl("inflamm|cytokine|interferon|tnf|nf.?kb|immune|leukocyte|chemokine|complement", x_low) ~
      "Inflammation / immune",
    grepl("apoptosis|p53|dna damage|stress|hypoxia|oxidative stress|injury|unfolded protein", x_low) ~
      "Stress / injury response",
    grepl("solute|channel|membrane potential|epithelial|sodium|potassium|transport", x_low) ~
      "Transport / epithelial function",
    TRUE ~ "Other"
  )
}

# -----------------------------
# Load data
# -----------------------------
OUTPUT_DIR <- make_output_dir(OUTPUT_DIR)
obj <- readRDS(INPUT_RDS)
DefaultAssay(obj) <- ASSAY_USE
check_required_metadata(obj, c("Phase", "exp.cond", "exp.time"))

# -----------------------------
# Experimental ordering
# -----------------------------
exp_time_order <- c(
  "Control",
  "IRI_short_1d", "IRI_short_3d", "IRI_short_14d",
  "IRI_long_1d",  "IRI_long_3d",  "IRI_long_14d"
)

focused_conditions <- c("Control", "IRI_short_3d", "IRI_long_3d")
focused_condition_labels <- c(
  Control = "Control",
  IRI_short_3d = "IRI short 3d",
  IRI_long_3d = "IRI long 3d"
)

injury_conditions <- c("Control", "IRI_short", "IRI_long")
cycling_levels <- c("Non_cycling", "Cycling")

# ==============================================================================
# 1. Define strict cycling status
# ===============================================================================
strict_cycling_genes <- c("Mki67", "Top2a", "Ube2c", "Cenpf", "Pcna", "Cdk1", "Ccnb1")
strict_cycling_genes <- intersect(strict_cycling_genes, rownames(obj))

if (length(strict_cycling_genes) == 0) {
  stop("None of the strict cycling marker genes were found in the Seurat object.")
}

strict_cycling_counts <- GetAssayData(
  obj,
  assay = ASSAY_USE,
  slot = "counts"
)[strict_cycling_genes, , drop = FALSE]

obj$Cycling_status <- "Non_cycling"
obj$Cycling_status[
  obj$Phase %in% c("S", "G2M") & Matrix::colSums(strict_cycling_counts > 0) > 0
] <- "Cycling"
obj$Cycling_status <- factor(obj$Cycling_status, levels = cycling_levels)

write.csv(
  as.data.frame(table(obj$exp.time, obj$Cycling_status)),
  file.path(OUTPUT_DIR, "cell_counts_by_time_and_cycling_status.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(obj$exp.cond, obj$Cycling_status)),
  file.path(OUTPUT_DIR, "cell_counts_by_condition_and_cycling_status.csv"),
  row.names = FALSE
)

# ==============================================================================
# 2. Cycling fraction across time points
# ===============================================================================
cycling_fraction_by_time <- obj@meta.data %>%
  filter(!is.na(Cycling_status), exp.time %in% exp_time_order) %>%
  count(exp.time, Cycling_status) %>%
  group_by(exp.time) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    exp.time = factor(exp.time, levels = exp_time_order),
    Cycling_status = factor(Cycling_status, levels = cycling_levels)
  )

write.csv(
  cycling_fraction_by_time,
  file.path(OUTPUT_DIR, "cycling_fraction_by_time.csv"),
  row.names = FALSE
)

condition_colors <- c(
  Control = "grey70",
  IRI_short_1d = "lightblue",
  IRI_short_3d = "blue",
  IRI_short_14d = "darkblue",
  IRI_long_1d = "coral",
  IRI_long_3d = "red",
  IRI_long_14d = "darkred"
)

annotation_df <- data.frame(
  exp.time = factor(exp_time_order, levels = exp_time_order),
  y = 103,
  Condition = factor(exp_time_order, levels = exp_time_order)
)

p_cycling_fraction <- ggplot(
  cycling_fraction_by_time,
  aes(x = exp.time, y = percent, fill = Cycling_status)
) +
  geom_col(position = "stack") +
  scale_fill_manual(
    name = "Cycling status",
    values = c(Non_cycling = "grey80", Cycling = "#8758AA"),
    labels = c(Non_cycling = "Non-cycling", Cycling = "Cycling")
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = annotation_df,
    aes(x = exp.time, y = y, fill = Condition),
    inherit.aes = FALSE,
    width = 0.95,
    height = 4
  ) +
  scale_fill_manual(name = "Group", values = condition_colors) +
  coord_cartesian(ylim = c(0, 108), clip = "off") +
  labs(x = NULL, y = "% cells") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

save_pdf(p_cycling_fraction, "cycling_status_percent_by_time.pdf", width = 10, height = 5)

# ==============================================================================
# 3. Ccnd1 expression by cycling status and time point
# ===============================================================================
if ("Ccnd1" %in% rownames(obj)) {
  ccnd1_summary <- FetchData(obj, vars = c("Ccnd1", "Cycling_status", "exp.time")) %>%
    mutate(
      Cycling_status = recode(as.character(Cycling_status), "Non_cycling" = "Non-cycling"),
      Cycling_status = factor(Cycling_status, levels = c("Non-cycling", "Cycling")),
      exp.time = factor(exp.time, levels = exp_time_order)
    ) %>%
    group_by(Cycling_status, exp.time) %>%
    summarise(mean_Ccnd1 = mean(Ccnd1, na.rm = TRUE), n_cells = n(), .groups = "drop")

  baseline_ccnd1 <- ccnd1_summary %>%
    filter(Cycling_status == "Non-cycling", exp.time == "Control") %>%
    pull(mean_Ccnd1)

  ccnd1_summary <- ccnd1_summary %>%
    mutate(mean_Ccnd1_fold_control = mean_Ccnd1 / baseline_ccnd1)

  write.csv(ccnd1_summary, file.path(OUTPUT_DIR, "ccnd1_summary_by_time_and_cycling_status.csv"), row.names = FALSE)

  p_ccnd1_heatmap <- ggplot(
    ccnd1_summary,
    aes(x = exp.time, y = Cycling_status, fill = mean_Ccnd1_fold_control)
  ) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(mean_Ccnd1_fold_control, 2)), size = 3) +
    scale_fill_viridis_c(option = "inferno", name = "Ccnd1 fold\nvs non-cycling\ncontrol") +
    theme_classic(base_size = 11) +
    labs(x = NULL, y = "Cell-cycle status") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_pdf(p_ccnd1_heatmap, "ccnd1_expression_heatmap_by_time_and_cycling_status.pdf", width = 7, height = 3)

  ccnd1_net_effect <- ccnd1_summary %>%
    select(exp.time, Cycling_status, mean_Ccnd1_fold_control) %>%
    pivot_wider(names_from = Cycling_status, values_from = mean_Ccnd1_fold_control) %>%
    mutate(cycling_minus_noncycling = Cycling - `Non-cycling`) %>%
    filter(!is.na(exp.time))

  write.csv(ccnd1_net_effect, file.path(OUTPUT_DIR, "ccnd1_cycling_minus_noncycling_net_effect.csv"), row.names = FALSE)

  p_ccnd1_net <- ggplot(ccnd1_net_effect, aes(x = exp.time, y = cycling_minus_noncycling)) +
    geom_col(fill = "grey70") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    theme_classic(base_size = 11) +
    labs(x = NULL, y = "Cycling - non-cycling\nCcnd1 fold change") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_pdf(p_ccnd1_net, "ccnd1_cycling_minus_noncycling_net_effect.pdf", width = 6, height = 3)
} else {
  warning("Ccnd1 was not found in rownames(obj); Ccnd1 panels were skipped.")
}

# ==============================================================================
# 4. Cycling vs non-cycling DEG and FGSEA for focused 3-day conditions
# ===============================================================================
load_mouse_pathways <- function(obj) {
  msig_cols <- msigdbr_collections(db_species = "MM")

  mouse_hallmark <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH") %>%
    mutate(collection_name = "Hallmark")
  mouse_go_bp <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M5", subcollection = "GO:BP") %>%
    mutate(collection_name = "GO_BP")
  mouse_reactome <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M2", subcollection = "CP:REACTOME") %>%
    mutate(collection_name = "Reactome")

  kegg_subcollections <- msig_cols %>%
    filter(gs_collection == "M2", grepl("KEGG", gs_subcollection)) %>%
    pull(gs_subcollection)

  mouse_kegg <- if (length(kegg_subcollections) > 0) {
    map_dfr(kegg_subcollections, ~ msigdbr(
      db_species = "MM", species = "Mus musculus", collection = "M2", subcollection = .x
    ) %>% mutate(collection_name = .x))
  } else {
    msigdbr(db_species = "MM", species = "Mus musculus", collection = "M2", subcollection = "CP:WIKIPATHWAYS") %>%
      mutate(collection_name = "WikiPathways")
  }

  pathway_df <- bind_rows(mouse_hallmark, mouse_go_bp, mouse_reactome, mouse_kegg) %>%
    mutate(pathway_collection = paste(collection_name, gs_name, sep = "__"))

  pathway_list <- split(pathway_df$gene_symbol, pathway_df$pathway_collection)
  pathway_list <- lapply(pathway_list, function(x) intersect(unique(x), rownames(obj)))
  pathway_list <- pathway_list[lengths(pathway_list) >= 10]

  pathway_annotation <- pathway_df %>%
    distinct(pathway_collection, collection_name, gs_name)

  list(pathways = pathway_list, annotation = pathway_annotation)
}

run_cycling_deg <- function(obj, exp_time_value, label, assay = ASSAY_USE, output_dir = OUTPUT_DIR) {
  meta <- obj@meta.data
  cells_cycling <- rownames(meta)[meta$exp.time == exp_time_value & meta$Cycling_status == "Cycling"]
  cells_noncycling <- rownames(meta)[meta$exp.time == exp_time_value & meta$Cycling_status == "Non_cycling"]

  if (length(cells_cycling) < 3 || length(cells_noncycling) < 3) {
    warning("Skipping ", label, ": fewer than 3 cells in at least one group.")
    return(NULL)
  }

  obj$DEG_group_tmp <- NA_character_
  obj$DEG_group_tmp[cells_cycling] <- "Cycling"
  obj$DEG_group_tmp[cells_noncycling] <- "Non_cycling"
  Idents(obj) <- "DEG_group_tmp"

  deg <- FindMarkers(
    object = obj,
    ident.1 = "Cycling",
    ident.2 = "Non_cycling",
    assay = assay,
    slot = "data",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.05,
    only.pos = FALSE
  ) %>%
    rownames_to_column("gene") %>%
    arrange(p_val_adj, desc(avg_log2FC))

  write.csv(deg, file.path(output_dir, paste0("DEG_", label, ".csv")), row.names = FALSE)
  deg
}

run_fgsea_from_deg <- function(deg, pathway_list, pathway_annotation, label, output_dir = OUTPUT_DIR) {
  if (is.null(deg)) return(NULL)

  ranks <- deg$avg_log2FC * -log10(deg$p_val + 1e-300)
  names(ranks) <- deg$gene
  ranks <- ranks[is.finite(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)

  fg <- fgsea(
    pathways = pathway_list,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    nPermSimple = 10000
  ) %>%
    as.data.frame() %>%
    left_join(pathway_annotation, by = c("pathway" = "pathway_collection")) %>%
    mutate(
      pathway_clean = clean_pathway_name(gs_name),
      pathway_family = assign_pathway_family(gs_name),
      direction = case_when(
        NES > 0 ~ "Enriched in cycling",
        NES < 0 ~ "Depleted in cycling",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(padj, desc(abs(NES)))

  if ("leadingEdge" %in% colnames(fg)) {
    fg$leadingEdge <- sapply(fg$leadingEdge, paste, collapse = ";")
  }

  data.table::fwrite(fg, file.path(output_dir, paste0("FGSEA_", label, ".csv")))
  fg
}

if (RUN_FGSEA) {
  pathway_data <- load_mouse_pathways(obj)

  deg_by_condition <- list()
  fgsea_by_condition <- list()

  for (cond in focused_conditions) {
    label <- paste0("Cycling_vs_NonCycling_", cond)
    deg_by_condition[[cond]] <- run_cycling_deg(obj, cond, label)
    fgsea_by_condition[[cond]] <- run_fgsea_from_deg(
      deg_by_condition[[cond]],
      pathway_data$pathways,
      pathway_data$annotation,
      label
    )
  }

  fgsea_long <- bind_rows(lapply(focused_conditions, function(cond) {
    fg <- fgsea_by_condition[[cond]]
    if (is.null(fg)) return(NULL)
    fg %>% mutate(condition = cond)
  }))

  write.csv(fgsea_long, file.path(OUTPUT_DIR, "FGSEA_all_results_focused_conditions_long.csv"), row.names = FALSE)

  fgsea_heatmap_df <- fgsea_long %>%
    filter(pathway_family %in% c("Cell cycle / DNA replication", "Metabolism")) %>%
    mutate(NES_sig = ifelse(!is.na(padj) & padj < PADJ_CUTOFF, NES, NA_real_)) %>%
    group_by(pathway, gs_name, pathway_clean, pathway_family) %>%
    filter(any(!is.na(NES_sig))) %>%
    ungroup()

  write.csv(fgsea_heatmap_df, file.path(OUTPUT_DIR, "FGSEA_cell_cycle_metabolism_heatmap_table.csv"), row.names = FALSE)

  if (nrow(fgsea_heatmap_df) > 0) {
    pathway_order <- fgsea_heatmap_df %>%
      group_by(pathway_clean, pathway_family) %>%
      summarise(max_abs_NES = max(abs(NES_sig), na.rm = TRUE), .groups = "drop") %>%
      arrange(pathway_family, desc(max_abs_NES)) %>%
      pull(pathway_clean)

    fgsea_heatmap_df <- fgsea_heatmap_df %>%
      mutate(
        pathway_family = factor(pathway_family, levels = c("Cell cycle / DNA replication", "Metabolism")),
        pathway_clean = factor(pathway_clean, levels = rev(unique(pathway_order))),
        condition = factor(condition, levels = focused_conditions, labels = focused_condition_labels[focused_conditions])
      )

    max_nes <- max(abs(fgsea_heatmap_df$NES_sig), na.rm = TRUE)

    p_fgsea_heatmap <- ggplot(fgsea_heatmap_df, aes(x = condition, y = pathway_clean, fill = NES_sig)) +
      geom_tile(color = "white", linewidth = 0.4) +
      facet_grid(pathway_family ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = c(-max_nes, max_nes), na.value = "grey90", name = "NES"
      ) +
      theme_classic(base_size = 13) +
      labs(x = NULL, y = NULL) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text.y = element_text(face = "bold", angle = 90)
      )

    save_pdf(
      p_fgsea_heatmap,
      "FGSEA_cell_cycle_metabolism_heatmap_focused_conditions.pdf",
      width = 7.5,
      height = max(5, 0.25 * length(unique(fgsea_heatmap_df$pathway_clean)) + 2)
    )
  }
}

# ==============================================================================
# 5. Canonical metabolic program module scoring
# ===============================================================================
canonical_metabolic_programs <- list(
  "Glycolysis / glucose metabolism" = c(
    "Hk1", "Hk2", "Gpi1", "Pfkl", "Pfkm", "Aldoa", "Aldob", "Tpi1", "Gapdh",
    "Pgk1", "Pgam1", "Eno1", "Eno2", "Pkm", "Ldha", "Slc2a1", "Slc2a2", "Pdk1", "Pdk3"
  ),
  "TCA cycle / OxPhos" = c(
    "Cs", "Aco2", "Idh3a", "Idh3b", "Ogdh", "Dlst", "Dld", "Suclg1", "Sdha", "Sdhb",
    "Fh1", "Mdh2", "Ndufa1", "Ndufa2", "Ndufb8", "Ndufs1", "Ndufs2", "Uqcrc1",
    "Uqcrc2", "Cyc1", "Cox4i1", "Cox5a", "Atp5f1a", "Atp5f1b", "Atp5mc1"
  ),
  "Fatty acid oxidation" = c(
    "Cpt1a", "Cpt2", "Slc25a20", "Acadm", "Acadvl", "Acads", "Acox1", "Acaa2", "Hadha",
    "Hadhb", "Hadh", "Ehhadh", "Eci1", "Eci2", "Decr1", "Etfa", "Etfb", "Etfdh", "Ppara", "Pdk4"
  ),
  "Fatty acid synthesis / lipogenesis" = c(
    "Acaca", "Acacb", "Fasn", "Acly", "Scd1", "Scd2", "Elovl1", "Elovl5", "Elovl6",
    "Fads1", "Fads2", "Me1", "Srebf1", "Mlxipl"
  ),
  "Phospholipid / membrane biogenesis" = c(
    "Gpam", "Gpat3", "Gpat4", "Agpat1", "Agpat2", "Agpat3", "Agpat4", "Agpat5", "Lpin1",
    "Lpin2", "Chka", "Chkb", "Pcyt1a", "Pcyt1b", "Etnk1", "Etnk2", "Pcyt2", "Cept1",
    "Chpt1", "Pemt", "Lpcat1", "Lpcat2", "Lpcat3", "Lpcat4", "Pla2g4a", "Pla2g6"
  ),
  "Cholesterol / sterol metabolism" = c(
    "Hmgcs1", "Hmgcr", "Mvk", "Pmvk", "Mvd", "Idi1", "Fdps", "Fdft1", "Sqle", "Lss",
    "Cyp51", "Msmo1", "Nsdhl", "Dhcr7", "Dhcr24", "Srebf2", "Insig1", "Ldlr"
  ),
  "Sphingolipid metabolism" = c(
    "Sptlc1", "Sptlc2", "Sptssa", "Kdsr", "Cers2", "Cers4", "Cers5", "Cers6", "Degs1",
    "Sgms1", "Sgms2", "Smpd1", "Smpd2", "Smpd3", "Asah1", "Asah2", "Sphk1", "Sphk2", "Sgpl1"
  ),
  "Amino acid metabolism" = c(
    "Gls", "Glud1", "Glul", "Got1", "Got2", "Psat1", "Phgdh", "Shmt1", "Shmt2", "Bcat1",
    "Bcat2", "Bckdha", "Bckdhb", "Ass1", "Asl", "Arg1", "Arg2", "Oat", "Prodh", "Gpt", "Gpt2"
  ),
  "Nucleotide metabolism" = c(
    "Prps1", "Prps2", "Ppat", "Gart", "Paics", "Adsl", "Atic", "Impdh1", "Impdh2", "Gmps",
    "Cad", "Dhfr", "Tyms", "Dut", "Rrm1", "Rrm2", "Tk1", "Mcm2", "Mcm3", "Mcm4"
  ),
  "Redox / glutathione metabolism" = c(
    "Gclc", "Gclm", "Gss", "Gsr", "Gsta1", "Gsta2", "Gstm1", "Gstm2", "Gstm3", "Gpx1",
    "Gpx3", "Gpx4", "Prdx1", "Prdx2", "Prdx3", "Prdx5", "Txn1", "Txn2", "Txnrd1", "Hmox1",
    "Nqo1", "Sod1", "Sod2", "Cat"
  ),
  "Peroxisomal metabolism" = c(
    "Acox1", "Acox2", "Ehhadh", "Hsd17b4", "Acaa1a", "Acaa1b", "Pex3", "Pex5", "Pex6",
    "Pex7", "Pex10", "Pex11a", "Pex11b", "Abcd1", "Abcd2", "Abcd3", "Crot", "Decr2"
  ),
  "Mitochondrial biogenesis / quality control" = c(
    "Ppargc1a", "Ppargc1b", "Nrf1", "Gabpa", "Tfam", "Tfb1m", "Tfb2m", "Polg", "Polg2",
    "Fis1", "Dnm1l", "Mfn1", "Mfn2", "Opa1", "Pink1", "Park2", "Sqstm1", "Map1lc3b"
  )
)

if (RUN_METABOLIC_MODULE_SCORING) {
  metabolic_programs_present <- lapply(canonical_metabolic_programs, function(x) intersect(unique(x), rownames(obj)))

  program_gene_counts <- tibble(
    program = names(metabolic_programs_present),
    n_genes_found = lengths(metabolic_programs_present),
    genes_found = sapply(metabolic_programs_present, paste, collapse = ";")
  )

  write.csv(program_gene_counts, file.path(OUTPUT_DIR, "canonical_metabolic_program_gene_counts.csv"), row.names = FALSE)

  metabolic_programs_present <- metabolic_programs_present[lengths(metabolic_programs_present) >= 3]

  obj <- AddModuleScore(
    object = obj,
    features = metabolic_programs_present,
    name = "canonical_met_",
    assay = ASSAY_USE
  )

  score_cols <- paste0("canonical_met_", seq_along(metabolic_programs_present))
  score_key <- tibble(score_col = score_cols, program = names(metabolic_programs_present))
  write.csv(score_key, file.path(OUTPUT_DIR, "canonical_metabolic_program_score_key.csv"), row.names = FALSE)

  score_df <- FetchData(obj, vars = c("exp.cond", "Cycling_status", score_cols)) %>%
    filter(exp.cond %in% injury_conditions, Cycling_status %in% cycling_levels) %>%
    mutate(
      exp.cond = factor(exp.cond, levels = injury_conditions),
      Cycling_status = factor(Cycling_status, levels = cycling_levels),
      group = paste(exp.cond, Cycling_status, sep = "_")
    )

  score_summary <- score_df %>%
    group_by(exp.cond, Cycling_status, group) %>%
    summarise(across(all_of(score_cols), ~ mean(.x, na.rm = TRUE)), n_cells = n(), .groups = "drop")

  score_long <- score_summary %>%
    pivot_longer(cols = all_of(score_cols), names_to = "score_col", values_to = "mean_score") %>%
    left_join(score_key, by = "score_col") %>%
    group_by(program) %>%
    mutate(z_score = as.numeric(scale(mean_score))) %>%
    ungroup()

  write.csv(score_long, file.path(OUTPUT_DIR, "canonical_metabolic_program_scores_by_condition_and_cycling_status.csv"), row.names = FALSE)

  diff_long <- score_long %>%
    select(exp.cond, Cycling_status, program, z_score) %>%
    pivot_wider(names_from = Cycling_status, values_from = z_score) %>%
    mutate(delta_z = Cycling - Non_cycling)

  diff_heatmap_mat <- diff_long %>%
    select(program, exp.cond, delta_z) %>%
    pivot_wider(names_from = exp.cond, values_from = delta_z) %>%
    column_to_rownames("program") %>%
    as.matrix()

  diff_heatmap_mat <- diff_heatmap_mat[, injury_conditions, drop = FALSE]

  pdf(file.path(OUTPUT_DIR, "metabolic_program_scores_cycling_minus_noncycling.pdf"), width = 7, height = 10)
  pheatmap(
    diff_heatmap_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    treeheight_row = 10,
    treeheight_col = 0,
    color = viridis::magma(100),
    breaks = seq(-1, 1, length.out = 101),
    border_color = NA,
    main = "Cycling - non-cycling",
    fontsize = 14,
    fontsize_row = 13,
    fontsize_col = 14,
    angle_col = 45,
    legend = TRUE,
    legend_breaks = c(-0.75, 0, 0.75),
    legend_labels = c("Lower", "0", "Higher"),
    cellwidth = 38,
    cellheight = 24
  )
  dev.off()
}

# ==============================================================================
# 6. Focused metabolic gene dot plot by condition and cycling status
# ===============================================================================
panel_genes_by_pathway <- list(
  "Glycolysis / glucose metabolism" = c("Slc2a4", "Hk2", "Pfkp", "Aldoa", "Pgam2", "Pkm", "Ldha"),
  "Amino acid metabolism" = c("Uroc1", "Acmsd", "Slc7a8", "Slc38a3", "Gls", "Gpt2", "Ass1", "Oat"),
  "Peroxisomal metabolism" = c("Acox1", "Ehhadh", "Hsd17b4", "Pex11a", "Car3", "Aldh3b2", "Aldh3b3", "Phyhd1"),
  "Fatty acid oxidation" = c("Cpt1a", "Cpt2", "Acadm", "Acadvl", "Hadha", "Hadhb", "Acaa2", "Ppara")
)

panel_genes <- unlist(panel_genes_by_pathway, use.names = FALSE)
panel_genes_present <- panel_genes[panel_genes %in% rownames(obj)]
panel_genes_missing <- setdiff(panel_genes, panel_genes_present)
writeLines(panel_genes_missing, con = file.path(OUTPUT_DIR, "focused_metabolic_genes_missing.txt"))

gene_to_pathway <- bind_rows(lapply(names(panel_genes_by_pathway), function(pathway_name) {
  data.frame(features.plot = panel_genes_by_pathway[[pathway_name]], pathway = pathway_name)
})) %>%
  filter(features.plot %in% panel_genes_present) %>%
  mutate(pathway = factor(pathway, levels = names(panel_genes_by_pathway)))

group_order_top_to_bottom <- c(
  "Control_Cycling", "IRI_short_Cycling", "IRI_long_Cycling",
  "Control_Non_cycling", "IRI_short_Non_cycling", "IRI_long_Non_cycling"
)
plot_group_order <- rev(group_order_top_to_bottom)

obj$focused_plot_group <- paste(obj$exp.cond, obj$Cycling_status, sep = "_")
obj_dot <- subset(obj, subset = exp.cond %in% injury_conditions & Cycling_status %in% cycling_levels)
obj_dot$focused_plot_group <- factor(obj_dot$focused_plot_group, levels = plot_group_order)
Idents(obj_dot) <- "focused_plot_group"

if (length(panel_genes_present) > 0) {
  dotplot_data <- DotPlot(obj_dot, features = panel_genes_present, dot.scale = 5)$data %>%
    left_join(gene_to_pathway, by = "features.plot") %>%
    mutate(
      features.plot = factor(features.plot, levels = panel_genes_present),
      id = factor(id, levels = plot_group_order),
      pathway = factor(pathway, levels = names(panel_genes_by_pathway))
    )

  write.csv(dotplot_data, file.path(OUTPUT_DIR, "focused_metabolic_gene_dotplot_data.csv"), row.names = FALSE)

  p_focused_metabolic_dotplot <- ggplot(
    dotplot_data,
    aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)
  ) +
    geom_point() +
    facet_grid(cols = vars(pathway), scales = "free_x", space = "free_x", switch = "x") +
    scale_color_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = "Scaled\navg. expr."
    ) +
    scale_size(range = c(0, 5), name = "% expr.") +
    theme_classic(base_size = 10) +
    labs(x = NULL, y = NULL) +
    guides(
      color = guide_colorbar(title.position = "top", direction = "vertical", barwidth = unit(0.35, "cm"), barheight = unit(2.2, "cm")),
      size = guide_legend(title.position = "top", direction = "vertical", ncol = 1)
    ) +
    theme(
      legend.position = "right",
      legend.box = "horizontal",
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(face = "bold", size = 8),
      panel.spacing.x = unit(0.45, "cm"),
      panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.35),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      axis.text.y = element_text(size = 8),
      plot.margin = unit(c(1, 0.5, 1, 0.5), "cm")
    )

  save_pdf(
    p_focused_metabolic_dotplot,
    "focused_metabolic_gene_dotplot_by_condition_and_cycling_status.pdf",
    width = 10,
    height = 5.5
  )
} else {
  warning("None of the focused metabolic genes were present in the object; dot plot was skipped.")
}

# -----------------------------
# Save processed object with added metadata/scores
# -----------------------------
saveRDS(obj, file.path(OUTPUT_DIR, "seurat_object_with_cycling_and_metabolic_scores.rds"))

message("Analysis complete. Outputs written to: ", OUTPUT_DIR)
