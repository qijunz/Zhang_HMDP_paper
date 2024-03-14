### Purpose: generate phenotype data of KO/taxon
### Created: 2023-08-06

source("/Users/rootqz/R/QZ_functions.R")
options(stringsAsFactors = FALSE)

setwd("/Users/rootqz/Desktop/HMDP/")

# load metadata
meta_wgs <- fread("paper/data/metadata_wgs.csv")
meta_mouse <- fread("paper/data/mouse_key.csv")

# load RSEM output, mapping to DO 1.9 million metagenes
rsem_out <- fread("JLusis_HMDP/data/Ath.HMDP.RSEM.DO1.9M.metagene.TPM.min0.in0.1samples_355mice.tsv") %>% as.data.frame()
names(rsem_out)[1] <- "gene"

# load 1.9M gene annotation, keep Bacteria gene only
gene_anno <- readRDS("JLusis_HMDP/data/DO_1.9M_NRGeneSet_KEGGanno_NCBIanno.rds") %>%
    filter(ncbi_superkingdom == "Bacteria")

# sum gene to KO level
phe_ko <- rsem_out %>%
    left_join(gene_anno %>% select(gene, KO), by = "gene", multiple = "all") %>%
    filter(KO != "") %>%
    select(-gene) %>%
    group_by(KO) %>%
    filter(n() >= 10) %>%
    summarise_all(list(~sum(.))) %>%
    column_to_rownames(var="KO") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    left_join(meta_wgs %>% select(Sample, IID, FID, PE_total), by = "Sample") %>%
    filter(PE_total >5000000) %>%
    select(-Sample, -PE_total) %>%
    select(IID, FID, everything())

fwrite(phe_ko, file = "paper/data/phe_ko.csv", sep = ",")

# sum gene to taxon level
# phylum level
phe_taxon_p <- rsem_out %>%
    left_join(gene_anno %>% select(gene, ncbi_phylum), by = "gene", multiple = "all") %>%
    filter(ncbi_phylum != "") %>%
    select(-gene) %>%
    group_by(ncbi_phylum) %>%
    filter(n() >= 1000) %>%
    summarise_all(list(~sum(.))) %>%
    column_to_rownames(var="ncbi_phylum") %>%
    t() %>%
    as.data.frame()

name_p <- as.data.frame(matrix(ncol = 2, nrow = ncol(phe_taxon_p)))
names(name_p) <- c("old", "new")
name_p$old <- names(phe_taxon_p)
name_p$new <- paste0("p_", name_p$old)
names(phe_taxon_p) <- name_p$new
freq_p <- table(gene_anno$ncbi_phylum) %>% as.data.frame()
name_p <- merge(name_p, freq_p, by.x = "old", by.y = "Var1", all.x = T, all.y = F)
name_p$level <- "phylum"

# class level
phe_taxon_c <- rsem_out %>%
    left_join(gene_anno %>% select(gene, ncbi_class), by = "gene", multiple = "all") %>%
    filter(ncbi_class != "") %>%
    select(-gene) %>%
    group_by(ncbi_class) %>%
    filter(n() >= 1000) %>%
    summarise_all(list(~sum(.))) %>%
    column_to_rownames(var="ncbi_class") %>%
    t() %>%
    as.data.frame()

name_c <- as.data.frame(matrix(ncol = 2, nrow = ncol(phe_taxon_c)))
names(name_c) <- c("old", "new")
name_c$old <- names(phe_taxon_c)
name_c$new <- paste0("c_", name_c$old)
names(phe_taxon_c) <- name_c$new
freq_c <- table(gene_anno$ncbi_class) %>% as.data.frame()
name_c <- merge(name_c, freq_c, by.x = "old", by.y = "Var1", all.x = T, all.y = F)
name_c$level <- "class"

# order level
phe_taxon_o <- rsem_out %>%
    left_join(gene_anno %>% select(gene, ncbi_order), by = "gene", multiple = "all") %>%
    filter(ncbi_order != "") %>%
    select(-gene) %>%
    group_by(ncbi_order) %>%
    filter(n() >= 1000) %>%
    summarise_all(list(~sum(.))) %>%
    column_to_rownames(var="ncbi_order") %>%
    t() %>%
    as.data.frame()

name_o <- as.data.frame(matrix(ncol = 2, nrow = ncol(phe_taxon_o)))
names(name_o) <- c("old", "new")
name_o$old <- names(phe_taxon_o)
name_o$new <- paste0("o_", name_o$old)
names(phe_taxon_o) <- name_o$new
freq_o <- table(gene_anno$ncbi_order) %>% as.data.frame()
name_o <- merge(name_o, freq_o, by.x = "old", by.y = "Var1", all.x = T, all.y = F)
name_o$level <- "order"

# family level
phe_taxon_f <- rsem_out %>%
    left_join(gene_anno %>% select(gene, ncbi_family), by = "gene", multiple = "all") %>%
    filter(ncbi_family != "") %>%
    select(-gene) %>%
    group_by(ncbi_family) %>%
    filter(n() >= 1000) %>%
    summarise_all(list(~sum(.))) %>%
    column_to_rownames(var="ncbi_family") %>%
    t() %>%
    as.data.frame()

name_f <- as.data.frame(matrix(ncol = 2, nrow = ncol(phe_taxon_f)))
names(name_f) <- c("old", "new")
name_f$old <- names(phe_taxon_f)
name_f$new <- paste0("f_", name_f$old)
names(phe_taxon_f) <- name_f$new
freq_f <- table(gene_anno$ncbi_family) %>% as.data.frame()
name_f <- merge(name_f, freq_f, by.x = "old", by.y = "Var1", all.x = T, all.y = F)
name_f$level <- "family"

# genus level
phe_taxon_g <- rsem_out %>%
    left_join(gene_anno %>% select(gene, ncbi_genus), by = "gene", multiple = "all") %>%
    filter(ncbi_genus != "") %>%
    select(-gene) %>%
    group_by(ncbi_genus) %>%
    filter(n() >= 500) %>%
    summarise_all(list(~sum(.))) %>%
    column_to_rownames(var="ncbi_genus") %>%
    t() %>%
    as.data.frame()

name_g <- as.data.frame(matrix(ncol = 2, nrow = ncol(phe_taxon_g)))
names(name_g) <- c("old", "new")
name_g$old <- names(phe_taxon_g)
name_g$new <- paste0("g_", name_g$old)
names(phe_taxon_g) <- name_g$new
freq_g <- table(gene_anno$ncbi_genus) %>% as.data.frame()
name_g <- merge(name_g, freq_g, by.x = "old", by.y = "Var1", all.x = T, all.y = F)
name_g$level <- "genus"

# species level
phe_taxon_s <- rsem_out %>%
    left_join(gene_anno %>% select(gene, ncbi_species), by = "gene", multiple = "all") %>%
    filter(ncbi_species != "") %>%
    select(-gene) %>%
    group_by(ncbi_species) %>%
    filter(n() >= 500) %>%
    summarise_all(list(~sum(.))) %>%
    column_to_rownames(var="ncbi_species") %>%
    t() %>%
    as.data.frame()

name_s <- as.data.frame(matrix(ncol = 2, nrow = ncol(phe_taxon_s)))
names(name_s) <- c("old", "new")
name_s$old <- names(phe_taxon_s)
name_s$new <- paste0("s_", name_s$old)
names(phe_taxon_s) <- name_s$new
freq_s <- table(gene_anno$ncbi_species) %>% as.data.frame()
name_s <- merge(name_s, freq_s, by.x = "old", by.y = "Var1", all.x = T, all.y = F)
name_s$level <- "species"

name_merge <- name_p %>%
    rbind(name_c) %>%
    rbind(name_o) %>%
    rbind(name_f) %>%
    rbind(name_g) %>%
    rbind(name_s)

phe_taxa <- phe_taxon_p %>%
    cbind(phe_taxon_c) %>%
    cbind(phe_taxon_o) %>%
    cbind(phe_taxon_f) %>%
    cbind(phe_taxon_g) %>%
    cbind(phe_taxon_s)

phe_taxon <- phe_taxa %>%
    rownames_to_column("Sample") %>%
    left_join(meta_wgs %>% select(Sample, IID, FID, PE_total), by = "Sample") %>%
    filter(PE_total >5000000) %>%
    select(-Sample, -PE_total) %>%
    select(IID, FID, everything())

fwrite(phe_taxon, file = "paper/data/phe_taxon.csv", sep = ",")

