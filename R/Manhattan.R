### Purpose: Manhattan plot for mb gwas
### Created: 2023-08-17

source("/Users/rootqz/R/QZ_functions.R")
options(stringsAsFactors = FALSE)

setwd("/Users/rootqz/Desktop/HMDP/")

# load metadata
meta_mouse <- fread("paper/data/mouse_key.csv")
meta_wgs <- fread("paper/data/metadata_wgs.csv")
meta_path <- fread("paper/data/metadata_pathway.csv")
meta_metaphlan4 <- fread("paper/data/metadata_metaphlan4.csv")

# load GWAS results
gwas_ko <- readRDS("paper/data/gwas_ko_Pvalue10-2.rds")
gwas_taxon <- readRDS("paper/data/gwas_taxon_Pvalue10-2.rds")
gwas_path <- readRDS("paper/data/gwas_path_Pvalue10-2.rds")
gwas_metaphlan4 <- readRDS("paper/data/gwas_metaphlan4_Pvalue10-2.rds")

## Manhattan plot
color1 <- c("1", "3", "5", "7", "9", "11", "13", "15", "17", "19")
color2 <- c("2", "4", "6", "8", "10", "12", "14", "16", "18", "X")

# KO
gwas_ko <- gwas_ko %>%
    mutate(P = -log10(PValue)) %>%
    mutate(color = case_when(
        snp_chr %in% color1 ~ "color1",
        snp_chr %in% color2 ~ "color2"
    )) 

# calculate cumulative marker position
gwas_ko_cum_pos <- gwas_ko %>%
    # compute chr size
    group_by(snp_chr) %>%
    summarise(chr_len = max(snp_bp_mm10)) %>%
    
    # calculate cumulative position of each chr
    mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) %>%
    
    # add this infor to extra column
    left_join(gwas_ko, ., by=c("snp_chr" = "snp_chr")) %>%
    
    # add a cumulative position of each marker
    arrange(snp_chr, snp_bp_mm10) %>%
    mutate(pos_cum = snp_bp_mm10+total)

gwas_ko_plot <- gwas_ko_cum_pos %>%
    group_by(SNP) %>%
    summarise_at(vars(PValue), min) %>%
    ungroup() %>%
    left_join(gwas_ko_cum_pos %>% select(SNP, snp_chr, snp_bp_mm10, color, pos_cum) %>% unique())

# make x axis with chr display
axis.df <- gwas_ko_plot %>%
    group_by(snp_chr) %>%
    summarise(center=(max(pos_cum) + min(pos_cum))/2)

setEPS()
postscript("paper/figure/Manhattan_ko.eps", height = 3, width = 12)
gwas_ko_plot %>%
    ggplot(aes(x=pos_cum, y=-log10(PValue))) +
    geom_point(aes(color=color), alpha=1, size=1) +
    geom_hline(yintercept = -log10(4*10^-6), color = "#A84448", linetype = "dashed") +
    
    # custom x,y axis
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.df$center, expand = c(0,1)) +
    scale_y_continuous(limits = c(2, 9), 
                       expand = expansion(add=c(0,0))) +
    theme_classic() +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = 0, color = "white"),
          legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    scale_color_manual(values = c("color1" = "#61677A", "color2" = "#D8D9DA")) +
    ylab("-log10(P)") + 
    xlab("Mouse Chromosome")
dev.off()

# taxon
gwas_taxon <- gwas_taxon %>%
    mutate(P = -log10(PValue)) %>%
    mutate(color = case_when(
        snp_chr %in% color1 ~ "color1",
        snp_chr %in% color2 ~ "color2"
    )) 

# calculate cumulative marker position
gwas_taxon_cum_pos <- gwas_taxon %>%
    # compute chr size
    group_by(snp_chr) %>%
    summarise(chr_len = max(snp_bp_mm10)) %>%
    
    # calculate cumulative position of each chr
    mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) %>%
    
    # add this infor to extra column
    left_join(gwas_taxon, ., by=c("snp_chr" = "snp_chr")) %>%
    
    # add a cumulative position of each marker
    arrange(snp_chr, snp_bp_mm10) %>%
    mutate(pos_cum = snp_bp_mm10+total)

gwas_taxon_plot <- gwas_taxon_cum_pos %>%
    group_by(SNP) %>%
    summarise_at(vars(PValue), min) %>%
    ungroup() %>%
    left_join(gwas_taxon_cum_pos %>% select(SNP, snp_chr, snp_bp_mm10, color, pos_cum) %>% unique())

# make x axis with chr display
axis.df <- gwas_taxon_plot %>%
    group_by(snp_chr) %>%
    summarise(center=(max(pos_cum) + min(pos_cum))/2)

setEPS()
postscript("paper/figure/Manhattan_taxon.eps", height = 3, width = 12)
gwas_taxon_plot %>%
    ggplot(aes(x=pos_cum, y=-log10(PValue))) +
    geom_point(aes(color=color), alpha=1, size=1) +
    geom_hline(yintercept = -log10(4*10^-6), color = "#A84448", linetype = "dashed") +
    
    # custom x,y axis
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.df$center, expand = c(0,1)) +
    scale_y_continuous(limits = c(2, 9), 
                       expand = expansion(add=c(0,0))) +
    theme_classic() +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = 0, color = "white"),
          legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    scale_color_manual(values = c("color1" = "#61677A", "color2" = "#D8D9DA")) +
    ylab("-log10(P)") + 
    xlab("Mouse Chromosome")
dev.off()

# pathway
gwas_path <- gwas_path %>%
    mutate(P = -log10(PValue)) %>%
    mutate(color = case_when(
        snp_chr %in% color1 ~ "color1",
        snp_chr %in% color2 ~ "color2"
    )) 

# calculate cumulative marker position
gwas_path_cum_pos <- gwas_path %>%
    # compute chr size
    group_by(snp_chr) %>%
    summarise(chr_len = max(snp_bp_mm10)) %>%
    
    # calculate cumulative position of each chr
    mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) %>%
    
    # add this infor to extra column
    left_join(gwas_path, ., by=c("snp_chr" = "snp_chr")) %>%
    
    # add a cumulative position of each marker
    arrange(snp_chr, snp_bp_mm10) %>%
    mutate(pos_cum = snp_bp_mm10+total)

gwas_path_plot <- gwas_path_cum_pos %>%
    group_by(SNP) %>%
    summarise_at(vars(PValue), min) %>%
    ungroup() %>%
    left_join(gwas_path_cum_pos %>% select(SNP, snp_chr, snp_bp_mm10, color, pos_cum) %>% unique())

# make x axis with chr display
axis.df <- gwas_path_plot %>%
    group_by(snp_chr) %>%
    summarise(center=(max(pos_cum) + min(pos_cum))/2)

setEPS()
postscript("paper/figure/Manhattan_pathway.eps", height = 3, width = 12)
gwas_path_plot %>%
    ggplot(aes(x=pos_cum, y=-log10(PValue))) +
    geom_point(aes(color=color), alpha=1, size=1) +
    geom_hline(yintercept = -log10(4*10^-6), color = "#A84448", linetype = "dashed") +
    
    # custom x,y axis
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.df$center, expand = c(0,1)) +
    scale_y_continuous(limits = c(2, 9), 
                       expand = expansion(add=c(0,0))) +
    theme_classic() +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = 0, color = "white"),
          legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    scale_color_manual(values = c("color1" = "#61677A", "color2" = "#D8D9DA")) +
    ylab("-log10(P)") + 
    xlab("Mouse Chromosome")
dev.off()

# metaphlan4 taxa
gwas_metaphlan4 <- gwas_metaphlan4 %>%
    mutate(P = -log10(PValue)) %>%
    mutate(color = case_when(
        snp_chr %in% color1 ~ "color1",
        snp_chr %in% color2 ~ "color2"
    )) 

# calculate cumulative marker position
gwas_metaphlan4_cum_pos <- gwas_metaphlan4 %>%
    # compute chr size
    group_by(snp_chr) %>%
    summarise(chr_len = max(snp_bp_mm10)) %>%
    
    # calculate cumulative position of each chr
    mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) %>%
    
    # add this infor to extra column
    left_join(gwas_metaphlan4, ., by=c("snp_chr" = "snp_chr")) %>%
    
    # add a cumulative position of each marker
    arrange(snp_chr, snp_bp_mm10) %>%
    mutate(pos_cum = snp_bp_mm10+total)

gwas_metaphlan4_plot <- gwas_metaphlan4_cum_pos %>%
    group_by(SNP) %>%
    summarise_at(vars(PValue), min) %>%
    ungroup() %>%
    left_join(gwas_metaphlan4_cum_pos %>% select(SNP, snp_chr, snp_bp_mm10, color, pos_cum) %>% unique())

# make x axis with chr display
axis.df <- gwas_metaphlan4_plot %>%
    group_by(snp_chr) %>%
    summarise(center=(max(pos_cum) + min(pos_cum))/2)

setEPS()
postscript("paper/figure/Manhattan_metaphlan4.eps", height = 3, width = 12)
gwas_metaphlan4_plot %>%
    ggplot(aes(x=pos_cum, y=-log10(PValue))) +
    geom_point(aes(color=color), alpha=1, size=1) +
    geom_hline(yintercept = -log10(4*10^-6), color = "#A84448", linetype = "dashed") +
    
    # custom x,y axis
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.df$center, expand = c(0,1)) +
    scale_y_continuous(limits = c(2, 9), 
                       expand = expansion(add=c(0,0))) +
    theme_classic() +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = 0, color = "white"),
          legend.position="none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    scale_color_manual(values = c("color1" = "#61677A", "color2" = "#D8D9DA")) +
    ylab("-log10(P)") + 
    xlab("Mouse Chromosome")
dev.off()


## density of metagenomic ko GWAS, window, cnt: -2 ~ +2 Mbp 
breaks <- matrix(c(seq(0, 200, 4), seq(1, 201, 4), seq(2, 202, 4), seq(3, 203, 4)), ncol = 4)
tmp <- as.list(1:ncol(breaks)) 
for(i in 1:ncol(breaks)) {
    tmp[[i]] <- gwas_ko %>%
        mutate(snp_bp_mm10 = snp_bp_mm10/1000000) %>%
        filter(PValue < 4*10^-6) %>%
        arrange(snp_chr, snp_bp_mm10) %>%
        group_by(snp_chr) %>%
        mutate(win = cut(snp_bp_mm10, breaks = breaks[,i])) %>%
        select(pheno, snp_chr, win) %>%
        unique() %>%
        group_by(snp_chr, win) %>% 
        summarize(cnt = n()) %>%
        separate(win, into = c("tmp1", "prox", "dist", "tmp2")) %>%
        mutate(prox = as.numeric(prox), 
               dist = as.numeric(dist), 
               mid = 0.5 * (prox + dist)) %>%
        select(snp_chr, mid, cnt)
}

ko_dens <- bind_rows(tmp[[1]], tmp[[2]], tmp[[3]], tmp[[4]])
rm(tmp)

pseudo <- gwas_ko %>%
    group_by(snp_chr) %>%
    mutate(max = max(snp_bp_mm10)/1000000, min = min(snp_bp_mm10)/1000000) %>%
    select(snp_chr, max, min) %>%
    unique() %>%
    as.data.frame() %>%
    reshape(idvar = "snp_chr", v.names = "snp_bp_mm10", 
            direction = "long", varying = c("min", "max"), new.row.names = 1:40) %>%
    select(-time) %>%
    mutate(cnt = 0)

ko_dens <- ko_dens %>%
    mutate(snp_bp_mm10 = mid) %>%
    ungroup() %>%
    select(snp_chr, snp_bp_mm10, cnt) %>%
    rbind(pseudo) %>%
    mutate(snp_chr = factor(snp_chr, levels = c("1","2","3","4","5","6","7","8","9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"))) 

ko_dens_plot <- ko_dens %>%
    # compute chr size
    group_by(snp_chr) %>%
    summarise(chr_len = max(snp_bp_mm10)) %>%
    
    # calculate cumulative position of each chr
    mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) %>%
    
    # add this infor to extra column
    left_join(ko_dens, ., by=c("snp_chr" = "snp_chr")) %>%
    
    # add a cumulative position of each marker
    arrange(snp_chr, snp_bp_mm10) %>%
    mutate(pos_cum = snp_bp_mm10+total)

# make x axis with chr display
axis.df <- ko_dens_plot %>%
    group_by(snp_chr) %>%
    summarise(center=(max(pos_cum) + min(pos_cum))/2)

setEPS()
postscript("paper/figure/Manhattan_ko_density.eps", height = 1.5, width = 12)
ko_dens_plot %>% ggplot() +
    geom_line(aes(x=pos_cum, y=cnt), color = "grey30") +
    # custom x,y axis
    scale_x_continuous(labels = c(1:19, "X"), breaks = axis.df$center, expand = c(0,1)) +
    scale_y_continuous(expand = expansion(add=c(0,0))) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    ylab("# traits") + 
    xlab("Mouse Chromosome")
dev.off()

## Manhattan plot of PC1
# function to plot individual gwas result
hmdp_gwas_plot <- function(gwas_out, y_min = 0, y_max = NA) {
    
    gwas_plot <- fread(gwas_out) %>% as.data.frame()
    
    color1 <- c("1", "3", "5", "7", "9", "11", "13", "15", "17", "19")
    color2 <- c("2", "4", "6", "8", "10", "12", "14", "16", "18", "X")
    
    gwas_plot <- gwas_plot %>%
        mutate(P = -log10(PValue)) %>%
        mutate(Chr = as.character(Chr)) %>%
        mutate(Chr = case_when(Chr == "23" ~ "X",
                               TRUE ~ Chr)) %>%
        mutate(Chr = factor(Chr, levels = c("1","2","3","4","5","6","7","8","9","10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"))) %>%
        mutate(color = case_when(
            Chr %in% color1 ~ "color1",
            Chr %in% color2 ~ "color2"
        ))
    
    # calculate cumulative marker position
    gwas_plot_cum <- gwas_plot %>%
        # compute chr size
        group_by(Chr) %>%
        summarise(chr_len = max(ChrPos)) %>%
        
        # calculate cumulative position of each chr
        mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
        select(-chr_len) %>%
        
        # add this infor to extra column
        left_join(gwas_plot, ., by=c("Chr" = "Chr")) %>%
        
        # add a cumulative position of each marker
        arrange(Chr, ChrPos) %>%
        mutate(pos_cum = ChrPos+total)
    
    # make x axis with chr display
    axis.df <- gwas_plot_cum %>%
        group_by(Chr) %>%
        summarise(center=(max(pos_cum) + min(pos_cum))/2)
    
    gwas_plot_cum %>%
        ggplot(aes(x=pos_cum, y=-log10(PValue))) +
        geom_point(aes(color=color), alpha=1, size=0.7) +
        geom_hline(yintercept = -log10(4*10^-6), color = "blue", linetype = "dashed") +
        #geom_hline(yintercept = -log10(4*10^-6/2000), color = "red", linetype = "dashed") +
        
        # custom x,y axis
        scale_x_continuous(labels = c(1:19, "X"), breaks = axis.df$center, expand = c(0,1)) +
        scale_y_continuous(limits = c(y_min, y_max) ,expand = expansion(add=c(0,0))) +
        theme_classic() +
        theme(panel.background = element_blank(),
              panel.border = element_rect(fill = 0, color = "white"),
              legend.position="none",
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 15)) +
        scale_color_manual(values = c("color1" = "#13547a", "color2" = "#80d0c7")) +
        ylab("-log10(P)") + 
        xlab("Mouse Chromosome")
}

setEPS()
postscript("paper/figure/Manhattan_ko_PC1.eps", height = 2, width = 12)
hmdp_gwas_plot("paper/data/gwas_mb/PC1.gwas", y_max = 7)
dev.off()


setEPS()
postscript("paper/figure/Manhattan_metaphlan4_Akkermansia.eps", height = 2, width = 12)
hmdp_gwas_plot("/Volumes/ReyLab_QZ10T/AthHMDP_GWAS/result/gwas_metaphlan4/taxon4.gwas", y_max = 7)
dev.off()
