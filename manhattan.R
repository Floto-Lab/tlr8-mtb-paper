library(tidyverse)
pval <- read_tsv("guides_w_pval_fc_min.txt")
annot1 <- read_tsv("TKOv1-base-90k-library-91320_sequences", col_names = F) %>% mutate(library = "base")
annot2 <- read_tsv("TKOv1-supp-85k-library-85180_sequences", col_names = F) %>% mutate(library = "sup")

annot <- bind_rows(annot1, annot2)
colnames(annot) <- c("seq", "guide_annot", "gene_annot", "library")
annot <- annot %>% 
    select(guide_annot, library) %>% 
    separate(guide_annot, c("pos_chr", "gene", "strand"), sep = "_") %>% 
    separate(pos_chr, c("chr", "pos"),  sep = ":") %>% 
    separate(pos, c("pos1", "pos2"), sep = "-") %>% 
    mutate(pos1 = as.integer(pos1))
annot <- annot %>% group_by(gene, library) %>% mutate(counter = row_number(gene))
pval <- pval %>% rename(guide = gene) %>% separate(guide, c("gene", "counter"), remove = F, sep = "_") %>% mutate(counter = as.integer(counter))
pval <- pval %>% inner_join(annot)

pval <- pval %>% filter(chr != "Luciferase", chr != "Gfp", chr != "EGFP", chr != "luciferase", chr != "LacZ") 
chrs <- c(paste("chr", seq(1,22), sep = ""), "chrX", "chrY")
nCHR <- length(unique(pval$chr))
pval$BPcum <- 0 
s <- 0
nbp <- c()
for (i in chrs){
      nbp[i] <- max(pval[pval$chr == i,]$pos1)
  pval[pval$chr == i,"BPcum"] <- pval[pval$chr == i,"pos1"] + s
    s <- s + nbp[i]
}

axis.set <- pval %>% 
  group_by(chr) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- abs(floor(log10(min(pval$pval)))) + 2 
sig <-0.05 
highlight_genes <- c("TLR8","ACSL1", "SERPINB2", "CPT1B", "CD5L", "ACAD8", "PTGES2", "ENPP1", "SLC39A7", "ATG7", "ATF3", "IL32", "TRAF3", "BATF2", "CASP6") 
chrs_grey = c("chr1", "chr3", "chr5", "chr7", "chr9", "chr11", "chr13", "chr15", "chr17", "chr19", "chr21", "chrX")
#druggable proteome
prot <- read_tsv("druggable_proteome.txt")
pval <- pval %>% mutate(highlight = ifelse((gene %in% highlight_genes) & (pval < 0.05) & (gene %in% prot$hgnc_names), "infection_regulation_druggable",
                                           ifelse((gene %in% prot$hgnc_names)& (pval < 0.05), "druggable", 
                                                  ifelse(((gene %in% highlight_genes) & (pval < 0.05) ), "infection_regulation",
                                                      ifelse(chr %in% chrs_grey, "uneven_chr", "even_chr" )))))
#pval <- pval%>% mutate(dot_size = ifelse(log10(fc_min) > 2, 6, ifelse(log10(fc_min) > 1, 5, ifelse(log10(fc_min) > 1, 4, 3))))

manhplot <- ggplot(pval) +  
	#aes(x = BPcum, y = -log10(fc_min), color =  factor(chr), size = -log10(pval)) +
	aes(x = BPcum, y = -log10(pval), color =  highlight) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = c("even_chr" = "grey88", "druggable" = "red", "infection_regulation_druggable" = "green", "uneven_chr" = "grey48", "infection_regulation" = "blue")) +
  #scale_size(range = c(0.1, 10), breaks = c(5, 20, 50, 100)) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )  +
  annotate("text", x =  (pval %>% filter(guide == "TLR8_1"))$BPcum, y = -log10((pval %>% filter(guide == "TLR8_1"))$pval), label = "TLR8_1")

    ggsave("manhattan_logfc_min.png", width = 15)
