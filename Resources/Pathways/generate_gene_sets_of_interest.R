immunogenic = openxlsx::read.xlsx("Resources/Pathways/GO-BP_c5.go.bp.v2024.1.Hs.symbols.gmt.xlsx",
                                  colNames = F)
colnames(immunogenic)[1:2] = c("id", "description")

gene_sets = c("GOBP_CELL_ADHESION",
                          "GOBP_CELL_CYCLE",
                          "GOBP_DEFENSE_RESPONSE",
                          "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
                          "GOBP_ACTIVATION_OF_IMMUNE_RESPONSE",
                          "GOBP_COMPLEMENT_ACTIVATION",
                          "GOBP_CHEMOKINE_PRODUCTION",
                          "GOBP_CYTOKINE_PRODUCTION",
                          "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
                          "GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
                          "GOBP_B_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
                          "GOBP_B_CELL_MEDIATED_IMMUNITY",
                          "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                          "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
                          "GOBP_T_CELL_MEDIATED_CYTOTOXICITY",
                          "GOBP_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", 
                          "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
                          "GOBP_MACROPHAGE_ACTIVATION",
                          "GOBP_MICROGLIAL_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
                          "GOBP_NK_T_CELL_ACTIVATION",
                          "GOBP_MAST_CELL_ACTIVATION",
                          "GOBP_NEUTROPHIL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
                          "GOBP_EOSINOPHIL_ACTIVATION",
                          "GOBP_MONOCYTE_ACTIVATION",
                          "GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY")

sets_names = c("Cell adhesion",
               "Cell cycle",
               "Defense response",
               "Regulation of the immune system",
               "Activation of the immune system",
               "Complement activation",
               "Chemokine production",
               "Cytokine production",
               "Antigen processing and presentation",
               "BCR signaling pathway",
               "B-cell activation",
               "B-cell-mediated immunity",
               "TCR signaling pathway",
               "T-cell activation",
               "T-cell-mediated cytotoxicity",
               "TLR signaling pathway",
               "Leukocyte-mediated immunity",
               "Macrophage activation",
               "Microglial cell activation",
               "NK-T cell activation",
               "Mast cell activation",
               "Neutrophil activation",
               "Eosinophil activation",
               "Monocyte activation",
               "TNF-mediated signaling")

gene_sets_of_interest = immunogenic[immunogenic$id %in% gene_sets, ]
gene_sets_of_interest$id = sets_names
gene_sets_of_interest$description = "na"

# Concatenate all gene columns into a single list of genes for each set (excluding NA)
genelists = lapply(1:nrow(gene_sets_of_interest), function(i) {
  # Extract the gene columns (from the 3rd column onwards) and remove NAs
  gene_row <- na.omit(unlist(gene_sets_of_interest[i, 3:ncol(gene_sets_of_interest)]))
  # Paste the genes with tab separation
  paste(gene_row, collapse = "\t")
})

gmt = cbind(gene_sets_of_interest$id, gene_sets_of_interest$description, unlist(genelists))

# Write the result to a GMT file
write.table(gene_sets_of_interest, "Resources/Pathways/gene_sets_of_interest.gmt",
            sep = "\t", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)
