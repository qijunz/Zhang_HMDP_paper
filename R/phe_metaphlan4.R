### Purpose: merge metaphlan4 output table
### Created: 2023-08-25

source("/Users/rootqz/R/QZ_functions.R")
options(stringsAsFactors = FALSE)

setwd("/Users/rootqz/Desktop/HMDP/")

# load wgs meta data
meta_wgs <- fread("paper/data/metadata_wgs.csv")

## SGB taxa
# merge taxon table
for (i in 1:nrow(meta_wgs)) {
    
    wgs_id <- meta_wgs$Sample[i]
    sample_IID <-  meta_wgs$IID[i]
    
    if (file.exists(paste0("/Users/rootqz/Desktop/HMDP/JLusis_HMDP/result/metaphlan4_out/", wgs_id,"_metaphlan_output.txt"))) {
        if (i == 1) {
            taxa_tbl <- fread(paste0("/Users/rootqz/Desktop/HMDP/JLusis_HMDP/result/metaphlan4_out/", wgs_id,"_metaphlan_output.txt")) %>% as.data.frame()
            names(taxa_tbl)[1] <- "clade_name"
            names(taxa_tbl)[3] <- wgs_id
            taxa_tbl$additional_species <- NULL
            
        } else {
            this_tbl <- fread(paste0("/Users/rootqz/Desktop/HMDP/JLusis_HMDP/result/metaphlan4_out/", wgs_id,"_metaphlan_output.txt")) %>% as.data.frame()
            names(this_tbl)[1] <- "clade_name"
            names(this_tbl)[3] <- wgs_id
            this_tbl$additional_species <- NULL
            
            taxa_tbl <- taxa_tbl %>%
                full_join(this_tbl, by = c("clade_name", "NCBI_tax_id"))
        }
    }
}

# get taxa metadata
meta_taxa <- taxa_tbl[,1:2]
meta_taxa$taxon_id <- paste0("taxon", c(1:nrow(meta_taxa)))

taxa_tbl[is.na(taxa_tbl)] <- 0
taxa_freq <- taxa_tbl %>% 
    column_to_rownames("clade_name") %>% 
    select(-NCBI_tax_id) %>% 
    t() %>%
    as.data.frame() %>%
    summarise_each(list(~sum(.==0))) %>%
    t() %>%
    as.data.frame() %>%
    rename_at("V1", ~"absent") %>%
    mutate(present = 355-absent) %>%
    rownames_to_column("clade_name")

taxa_tbl[taxa_tbl == 0] <- NA
taxa_abund <- taxa_tbl %>% 
    column_to_rownames("clade_name") %>% 
    select(-NCBI_tax_id) %>% 
    t() %>%
    as.data.frame() %>%
    summarise_each(list(~mean(., na.rm = T))) %>%
    t() %>%
    as.data.frame() %>%
    rename_at("V1", ~"rel_abundance") %>%
    rownames_to_column("clade_name")

taxa_split <- meta_taxa %>%
    select(clade_name) %>%
    mutate(level = NA, kingdom = NA, phylum = NA, class = NA, order = NA, family = NA, genus = NA, species = NA, sgb = NA)

for (i in 1:nrow(taxa_split)) {
    this_split <- strsplit(taxa_split$clade_name[i], "\\|")[[1]]
    
    if (length(this_split) == 1) {
        taxa_split[i,2:10] = c("kingdom", this_split[1], NA, NA, NA, NA, NA, NA, NA)
    } else if (length(this_split) == 2) {
        taxa_split[i,2:10] = c("phylum", this_split[1], this_split[2], NA, NA, NA, NA, NA, NA)
    } else if (length(this_split) == 3) {
        taxa_split[i,2:10] = c("class", this_split[1], this_split[2], this_split[3], NA, NA, NA, NA, NA)
    } else if (length(this_split) == 4) {
        taxa_split[i,2:10] = c("order", this_split[1], this_split[2], this_split[3], this_split[4], NA, NA, NA, NA)
    } else if (length(this_split) == 5) {
        taxa_split[i,2:10] = c("family", this_split[1], this_split[2], this_split[3], this_split[4], this_split[5], NA, NA, NA)
    } else if (length(this_split) == 6) {
        taxa_split[i,2:10] = c("genus", this_split[1], this_split[2], this_split[3], this_split[4], this_split[5], this_split[6], NA, NA)
    } else if (length(this_split) == 7) {
        taxa_split[i,2:10] = c("species", this_split[1], this_split[2], this_split[3], this_split[4], this_split[5], this_split[6], this_split[7], NA)
    } else if (length(this_split) == 8) {
        taxa_split[i,2:10] = c("sgb", this_split[1], this_split[2], this_split[3], this_split[4], this_split[5], this_split[6], this_split[7], this_split[8])
    }
}

meta_metaphlan4 <- meta_taxa %>%
    left_join(taxa_abund, by = "clade_name") %>%
    left_join(taxa_freq, by = "clade_name") %>%
    left_join(taxa_split, by = "clade_name") %>%
    select(taxon_id, everything())

fwrite(meta_metaphlan4, file = "paper/data/metadata_metaphlan4.csv", sep = ",")

# write out phenotype
phe_metaphlan4 <- taxa_tbl %>%
    left_join(meta_taxa %>% select(taxon_id, clade_name), by = "clade_name") %>%
    select(-clade_name, -NCBI_tax_id) %>%
    left_join(meta_metaphlan4 %>% select(taxon_id, rel_abundance, present), by = "taxon_id") %>%
    filter(rel_abundance > 0.01) %>%
    filter(present > 70) %>%
    select(-rel_abundance, -present) %>%
    column_to_rownames("taxon_id") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    left_join(meta_wgs %>% select(Sample, IID, FID, PE_total), by = "Sample") %>%
    filter(PE_total >5000000) %>%
    select(-Sample, -PE_total) %>%
    select(IID, FID, everything())

phe_metaphlan4[is.na(phe_metaphlan4)] <- 0
fwrite(phe_metaphlan4, file = "paper/data/phe_metaphlan4.csv", sep = ",")


## GTDB taxa
# merge taxa table
for (i in 1:nrow(meta_wgs)) {
    
    wgs_id <- meta_wgs$Sample[i]
    sample_IID <-  meta_wgs$IID[i]
    
    if (file.exists(paste0("/Users/rootqz/Desktop/HMDP/JLusis_HMDP/result/metaphlan4_gtdb_out/", wgs_id,"_metaphlan_gtdb.txt"))) {
        if (i == 1) {
            taxa_tbl <- fread(paste0("/Users/rootqz/Desktop/HMDP/JLusis_HMDP/result/metaphlan4_gtdb_out/", wgs_id,"_metaphlan_gtdb.txt")) %>% as.data.frame()
            names(taxa_tbl)[1] <- "clade_name"
            names(taxa_tbl)[2] <- wgs_id
            
        } else {
            this_tbl <- fread(paste0("/Users/rootqz/Desktop/HMDP/JLusis_HMDP/result/metaphlan4_gtdb_out/", wgs_id,"_metaphlan_gtdb.txt")) %>% as.data.frame()
            names(this_tbl)[1] <- "clade_name"
            names(this_tbl)[2] <- wgs_id
            
            taxa_tbl <- taxa_tbl %>%
                full_join(this_tbl, by = "clade_name")
        }
    }
}

# get taxa metadata
meta_taxa <- taxa_tbl[,1,drop=F]
meta_taxa$taxon_id <- paste0("taxon", c(1:nrow(meta_taxa)))

taxa_tbl[is.na(taxa_tbl)] <- 0
taxa_freq <- taxa_tbl %>% 
    column_to_rownames("clade_name") %>% 
    t() %>%
    as.data.frame() %>%
    summarise_each(list(~sum(.==0))) %>%
    t() %>%
    as.data.frame() %>%
    rename_at("V1", ~"absent") %>%
    mutate(present = 355-absent) %>%
    rownames_to_column("clade_name")

taxa_tbl[taxa_tbl == 0] <- NA
taxa_abund <- taxa_tbl %>% 
    column_to_rownames("clade_name") %>% 
    t() %>%
    as.data.frame() %>%
    summarise_each(list(~mean(., na.rm = T))) %>%
    t() %>%
    as.data.frame() %>%
    rename_at("V1", ~"rel_abundance") %>%
    rownames_to_column("clade_name")

taxa_split <- meta_taxa %>%
    select(clade_name) %>%
    mutate(level = NA, kingdom = NA, phylum = NA, class = NA, order = NA, family = NA, genus = NA, species = NA)

for (i in 1:nrow(taxa_split)) {
    this_split <- strsplit(taxa_split$clade_name[i], "\\;")[[1]]
    
    if (length(this_split) == 1) {
        taxa_split[i,2:9] = c("kingdom", this_split[1], NA, NA, NA, NA, NA, NA)
    } else if (length(this_split) == 2) {
        taxa_split[i,2:9] = c("phylum", this_split[1], this_split[2], NA, NA, NA, NA, NA)
    } else if (length(this_split) == 3) {
        taxa_split[i,2:9] = c("class", this_split[1], this_split[2], this_split[3], NA, NA, NA, NA)
    } else if (length(this_split) == 4) {
        taxa_split[i,2:9] = c("order", this_split[1], this_split[2], this_split[3], this_split[4], NA, NA, NA)
    } else if (length(this_split) == 5) {
        taxa_split[i,2:9] = c("family", this_split[1], this_split[2], this_split[3], this_split[4], this_split[5], NA, NA)
    } else if (length(this_split) == 6) {
        taxa_split[i,2:9] = c("genus", this_split[1], this_split[2], this_split[3], this_split[4], this_split[5], this_split[6], NA)
    } else if (length(this_split) == 7) {
        taxa_split[i,2:9] = c("species", this_split[1], this_split[2], this_split[3], this_split[4], this_split[5], this_split[6], this_split[7])
    } 
}

meta_gtdb <- meta_taxa %>%
    left_join(taxa_abund, by = "clade_name") %>%
    left_join(taxa_freq, by = "clade_name") %>%
    left_join(taxa_split, by = "clade_name") %>%
    select(taxon_id, everything())

fwrite(meta_gtdb, file = "paper/data/metadata_gtdb.csv", sep = ",")

# write out phenotype, rel abundance > 0.01% and present in > 5% of samples
phe_gtdb <- taxa_tbl %>%
    left_join(meta_taxa %>% select(taxon_id, clade_name), by = "clade_name") %>%
    select(-clade_name) %>%
    left_join(meta_gtdb %>% select(taxon_id, rel_abundance, present), by = "taxon_id") %>%
    filter(rel_abundance > 0.01) %>%
    filter(present > 17) %>%
    select(-rel_abundance, -present) %>%
    column_to_rownames("taxon_id") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    left_join(meta_wgs %>% select(Sample, IID, FID, PE_total), by = "Sample") %>%
    filter(PE_total >5000000) %>%
    select(-Sample, -PE_total) %>%
    select(IID, FID, everything())

phe_gtdb[is.na(phe_gtdb)] <- 0
fwrite(phe_gtdb, file = "paper/data/phe_gtdb.csv", sep = ",")

