### Purpose: run gwas for KOs and taxa
### Created: 2023-08-06

source("/Users/rootqz/R/QZ_functions.R")
options(stringsAsFactors = FALSE)

setwd("/Users/rootqz/Desktop/HMDP/")

# load metadata
meta_wgs <- fread("paper/data/metadata_wgs.csv")

# load pheno data
phe_ko <- fread("paper/data/phe_ko.csv") %>% as.data.frame()
phe_taxon <- fread("paper/data/phe_taxon.csv") %>% as.data.frame()
phe_path <- fread("paper/data/phe_pathway.csv") %>% as.data.frame()
phe_metaphlan4 <- fread("paper/data/phe_metaphlan4.csv") %>% as.data.frame()
phe_gtdb <- fread("paper/data/phe_gtdb.csv") %>% as.data.frame()

# load geno data for all strain
geno <- fread("GWAS/data/all_strains.tped")
src_strains <- colnames(geno)[-c(1:4)]

# gwas for trait KO/taxon/pathway
traits <- names(phe_ko)[-c(1:2)]
traits <- names(phe_taxon)[-c(1:2)]
traits <- names(phe_path)[-c(1:2)]
traits <- names(phe_metaphlan4)[-c(1:2)]
traits <- names(phe_gtdb)[-c(1:4)]

for (trait in traits) {
    cat(paste0(trait, "...\n"))
    
    pheno <- phe_gtdb[,c("FID", "IID", trait)] %>% na.omit()
    
    targets <- pheno %>% select(FID, IID)
    targets %<>% filter(FID %in% src_strains)
    
    # write .tped file
    geno[, c(colnames(geno)[1:4], targets$FID), with = F] %>%
        fwrite(paste0("paper/data/gwas_mb/tmp/", trait, ".tped"), sep = '\t', col.names = F)
    
    # write .tfam file
    tibble(FID = targets$FID %>% str_replace_all('[ /]', '.'),
           IID = targets$IID %>% str_replace_all('[ /]', '.'),
           fa = 0,
           ma = 0,
           sex = 0,
           pheno = -9) %>%
        fwrite(paste0("paper/data/gwas_mb/tmp/", trait, ".tfam"), sep = '\t', col.names = F)
    
    # write .pheno.txt file, rankZ transformed
    pheno %>% 
        mutate_at(1, ~str_replace_all(.x, '[ /]', '.')) %>%
        mutate_at(3, ~rankZ(.x)) %>%
        fwrite(paste0("paper/data/gwas_mb/tmp/", trait, ".pheno.txt"), sep = '\t', col.names = F)
    
    # write .cov.txt file
    pheno %>% 
        left_join(meta_wgs %>% select(IID, Sex, Wave), by = "IID") %>%
        select(FID, IID, Sex, Wave) %>%
        mutate_at(1, ~str_replace_all(.x, '[ /]', '.')) %>%
        mutate_at(3, ~ifelse(.x=="F", 0,1)) %>%
        fwrite(paste0("paper/data/gwas_mb/tmp/", trait, ".cov.txt"), sep = '\t', col.names = F)
    
    # run plink to get bed file
    system(paste0("/Users/rootqz/bin/plink --tfile paper/data/gwas_mb/tmp/", trait, " --allow-no-sex --maf 0.05 --geno 0.1 --make-bed --out paper/data/gwas_mb/tmp/", trait))
    
    # run GWAS Python script
    system(paste0("/Users/rootqz/anaconda3/bin/python paper/script/run_gwas.py -p paper/data/gwas_mb/tmp/", trait))
    
    # move .gwas output file
    system(paste0("mv paper/data/gwas_mb/tmp/", trait, ".gwas paper/data/gwas_mb/"))
    
    # clean space
    system(paste0("rm paper/data/gwas_mb/tmp/", trait, ".*"))
}


## merge gwas reults
# load SNP dictionary converting JAX ID to rsID
snp_dict <- fread("GWAS/data/jax_id_rsID.txt")

# KO
KEGG_gene <- fread("~/Desktop/DO/metagenomics/2020/data/KEGG_KO2gene2EC.tsv")

gwas_merge <- as.data.frame(matrix(nrow = 0, ncol = 5))
names(gwas_merge) <- c("pheno", "SNP", "Chr", "ChrPos", "PValue")

for (trait in names(phe_ko)[-c(1:2)]) {
    cat(paste0(trait, "...\n"))
    
    if (file.exists(paste0("/Volumes/ReyLab_QZ10T/AthHMDP_GWAS/result/gwas_ko/", trait,".gwas"))) {
        this_gwas <- fread(paste0("/Volumes/ReyLab_QZ10T/AthHMDP_GWAS/result/gwas_ko/", trait,".gwas"))
        this_gwas_trim <- this_gwas %>%
            mutate(pheno = trait) %>%
            select(pheno, SNP, Chr, ChrPos, PValue) %>%
            filter(PValue < 0.0001)
        
        gwas_merge <<- rbind(gwas_merge, this_gwas_trim)
    }
}

gwas_ko <- gwas_merge %>%
    left_join(snp_dict, by = c("SNP" = "snp_id")) %>%
    select(pheno, SNP, rsID, snp_chr, snp_bp_mm10, PValue) %>%
    mutate(snp_chr = factor(snp_chr, levels = c("1","2","3","4","5","6","7","8","9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"))) %>%
    left_join(KEGG_gene, by = c("pheno" = "KO"))

saveRDS(gwas_ko, file = "paper/data/gwas_ko_Pvalue10-4.rds")

# taxon
gwas_merge <- as.data.frame(matrix(nrow = 0, ncol = 5))
names(gwas_merge) <- c("pheno", "SNP", "Chr", "ChrPos", "PValue")

for (trait in names(phe_taxon)[-c(1:2)]) {
    cat(paste0(trait, "...\n"))
    
    if (file.exists(paste0("/Volumes/ReyLab_QZ10T/AthHMDP_GWAS/result/gwas_taxon/", trait,".gwas"))) {
        this_gwas <- fread(paste0("/Volumes/ReyLab_QZ10T/AthHMDP_GWAS/result/gwas_taxon/", trait,".gwas"))
        this_gwas_trim <- this_gwas %>%
            mutate(pheno = trait) %>%
            select(pheno, SNP, Chr, ChrPos, PValue) %>%
            filter(PValue < 0.0001)
        
        gwas_merge <<- rbind(gwas_merge, this_gwas_trim)
    }
}

gwas_taxon <- gwas_merge %>%
    left_join(snp_dict, by = c("SNP" = "snp_id")) %>%
    select(pheno, SNP, rsID, snp_chr, snp_bp_mm10, PValue) %>%
    mutate(snp_chr = factor(snp_chr, levels = c("1","2","3","4","5","6","7","8","9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"))) 

saveRDS(gwas_taxon, file = "paper/data/gwas_taxon_Pvalue10-4.rds")

# pathway
meta_path <- fread("paper/data/metadata_pathway.csv")

gwas_merge <- as.data.frame(matrix(nrow = 0, ncol = 5))
names(gwas_merge) <- c("pheno", "SNP", "Chr", "ChrPos", "PValue")

for (trait in names(phe_path)[-c(1:2)]) {
    cat(paste0(trait, "...\n"))
    
    if (file.exists(paste0("/Volumes/ReyLab_QZ10T/AthHMDP_GWAS/result/gwas_path/", trait,".gwas"))) {
        this_gwas <- fread(paste0("/Volumes/ReyLab_QZ10T/AthHMDP_GWAS/result/gwas_path/", trait,".gwas"))
        this_gwas_trim <- this_gwas %>%
            mutate(pheno = trait) %>%
            select(pheno, SNP, Chr, ChrPos, PValue) %>%
            filter(PValue < 0.0001)
        
        gwas_merge <<- rbind(gwas_merge, this_gwas_trim)
    }
}

gwas_path <- gwas_merge %>%
    left_join(snp_dict, by = c("SNP" = "snp_id")) %>%
    select(pheno, SNP, rsID, snp_chr, snp_bp_mm10, PValue) %>%
    mutate(snp_chr = factor(snp_chr, levels = c("1","2","3","4","5","6","7","8","9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"))) %>%
    left_join(meta_path, by = c("pheno" = "path_id"))

saveRDS(gwas_path, file = "paper/data/gwas_path_Pvalue10-4.rds")

# metaphlan4
meta_metaphlan4 <- fread("paper/data/metadata_metaphlan4.csv")

gwas_merge <- as.data.frame(matrix(nrow = 0, ncol = 5))
names(gwas_merge) <- c("pheno", "SNP", "Chr", "ChrPos", "PValue")

for (trait in names(phe_metaphlan4)[-c(1:2)]) {
    cat(paste0(trait, "...\n"))
    
    if (file.exists(paste0("paper/data/gwas_mb/", trait,".gwas"))) {
        this_gwas <- fread(paste0("paper/data/gwas_mb/", trait,".gwas"))
        this_gwas_trim <- this_gwas %>%
            mutate(pheno = trait) %>%
            select(pheno, SNP, Chr, ChrPos, PValue) %>%
            filter(PValue < 0.0001)
        
        gwas_merge <<- rbind(gwas_merge, this_gwas_trim)
    }
}

gwas_metaphlan4 <- gwas_merge %>%
    left_join(snp_dict, by = c("SNP" = "snp_id")) %>%
    select(pheno, SNP, rsID, snp_chr, snp_bp_mm10, PValue) %>%
    mutate(snp_chr = factor(snp_chr, levels = c("1","2","3","4","5","6","7","8","9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"))) %>%
    left_join(meta_metaphlan4 %>% select(clade_name, level, kingdom, phylum, class, order, family, genus, species, sgb, everything()), by = c("pheno" = "taxon_id"))

saveRDS(gwas_metaphlan4, file = "paper/data/gwas_metaphlan4_Pvalue10-4.rds")

