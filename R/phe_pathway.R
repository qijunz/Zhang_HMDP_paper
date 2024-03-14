### Purpose: generate phenotype data of metacyc pathway, output from humann3
### Created: 2023-08-08

source("/Users/rootqz/R/QZ_functions.R")
options(stringsAsFactors = FALSE)

setwd("/Users/rootqz/Desktop/HMDP/")

# load metadata
meta_wgs <- fread("paper/data/metadata_wgs.csv")
meta_mouse <- fread("paper/data/mouse_key.csv")

# merge humann3 metacyc table
setwd("/Users/rootqz/Desktop/HMDP/JLusis_HMDP/result/humann3_out/")

for (i in 1:nrow(meta_wgs)) {
    
    wgs_id <- meta_wgs$Sample[i]
    sample_IID <-  meta_wgs$IID[i]
    
    if (file.exists(paste0("humann3_out_tbl_", wgs_id,".tar.gz"))) {
        system(paste0("tar -zxf humann3_out_tbl_", wgs_id,".tar.gz"))
        
        if (file.exists(paste0("humann3_out_tbl_", wgs_id,"/", wgs_id, "_pathabundance.tsv"))) {
            if (i == 1) {
                metacyc_tbl <- fread(paste0("humann3_out_tbl_", wgs_id,"/", wgs_id, "_pathabundance.tsv"))
                names(metacyc_tbl)[1] <- "gene_family"
                names(metacyc_tbl)[2] <- wgs_id
                
            } else {
                this_tbl <- fread(paste0("humann3_out_tbl_", wgs_id,"/", wgs_id, "_pathabundance.tsv"))
                names(this_tbl)[1] <- "gene_family"
                names(this_tbl)[2] <- wgs_id
                
                metacyc_tbl <- metacyc_tbl %>%
                    full_join(this_tbl, by = "gene_family")
            }
        }
        
        system(paste0("rm -r humann3_out_tbl_", wgs_id))
    }
}
setwd("/Users/rootqz/Desktop/HMDP/")

fwrite(metacyc_tbl, file = "paper/data/humann3_metacyc.csv", sep = ",")

# extract UNMAPPED reads, total UNINTEGRATED features, and pathways that are not broken down per organism.
metacyc_tbl_unique <- metacyc_tbl[1,]
for (i in 2:nrow(metacyc_tbl)) {
    if (!str_detect(metacyc_tbl$gene_family[i], "\\|g\\_\\_") && !str_detect(metacyc_tbl$gene_family[i], "\\|unclassified")) {
        metacyc_tbl_unique <- rbind(metacyc_tbl_unique, metacyc_tbl[i,])
    }
}

metacyc_tbl_unique[is.na(metacyc_tbl_unique)] <- 0

path_norm <- metacyc_tbl_unique[,1,drop=F]
for (i in 2:ncol(metacyc_tbl_unique)) {
    this_norm <- metacyc_tbl_unique[,..i,drop=F]
    this_norm[,1] <- this_norm[,1]/sum(this_norm[,1])*1000000
    
    path_norm <- cbind(path_norm, this_norm)
}

# create metadata for pathway definition, use shot words for pheno name
meta_path <- as.data.frame(matrix(nrow = nrow(path_norm), ncol = 2))
names(meta_path) <- c("path_id", "path_name")
meta_path$path_id <- paste0("path", c(1:nrow(path_norm)))
meta_path$path_name <- path_norm$gene_family

phe_path <- path_norm %>%
    mutate(path_id = paste0("path", c(1:nrow(meta_path)))) %>%
    select(-gene_family) %>%
    column_to_rownames("path_id") %>%
    t() %>%
    as.data.frame()

path_freq <- phe_path %>% summarise_each(list(~sum(.==0)))

path_freq_merge <- path_freq %>% 
    t() %>% 
    as.data.frame() %>%
    rownames_to_column("path_id") %>%
    mutate(absent = V1) %>%
    select(-V1) %>%
    mutate(present = 355-absent)

meta_path_out <- meta_path %>%
    left_join(path_freq_merge, by = "path_id")

fwrite(meta_path_out, file = "paper/data/metadata_pathway.csv", sep = ",")

phe_path_filtered <- phe_path[meta_wgs$Sample[which(meta_wgs$filtered_gwas == "yes")],
                              names(path_freq)[which(path_freq < length(meta_wgs$filtered_gwas == "yes")*0.8)]] %>%
    na.omit() %>%
    select(-path1, -path2) %>%
    rownames_to_column("Sample") %>%
    left_join(meta_wgs %>% select(Sample, IID, FID, PE_total), by = "Sample") %>%
    filter(PE_total >5000000) %>%
    select(-Sample, -PE_total) %>%
    select(IID, FID, everything())

fwrite(phe_path_filtered, file = "paper/data/phe_pathway.csv", sep = ",")

# PCA for 300 pathways
# load pheno data
phe_ko <- fread("paper/data/phe_ko.csv") %>% as.data.frame()
phe_taxon <- fread("paper/data/phe_taxon.csv") %>% as.data.frame()
phe_path <- fread("paper/data/phe_pathway.csv") %>% as.data.frame()

pr.out <- prcomp(phe_path %>% column_to_rownames("IID") %>% select(-FID), scale. = T)
pr.var <- pr.out$sdev^2
pve_path <- pr.var/sum(pr.var)
pca_matrix <- pr.out$x
pca_1and2 <- as.data.frame(pca_matrix[,1:2])

pca_plot <- pca_1and2 %>%
    rownames_to_column("IID") %>%
    mutate(IID = as.numeric(IID)) %>%
    left_join(meta_wgs, by = "IID") %>%
    mutate(Wave = factor(Wave)) %>%
    left_join(phe_ko %>% select(IID, K02406), by = "IID")

pca_plot %>% ggplot() + 
    geom_point(aes(x = PC1, y = PC2, color = K02406)) + 
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)) +
    xlab(paste0("PC1 (", round(pve_path[1]*100, 2), "%)")) +
    ylab(paste0("PC2 (", round(pve_path[2]*100, 2), "%)")) +
    scale_color_viridis_c(begin = 1, end = 0)
