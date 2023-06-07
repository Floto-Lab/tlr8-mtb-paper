#library(CB2)
library(glue)
library(tidyverse)
library(corrplot)
library(CB2)

#load sample table
sf <- read_tsv("sample_table.txt")
#create unique identifier for each sample
sf <- sf %>% unite(sample_name, condition, experiment, remove = F) %>% 
    rename(group = condition)  
#keep only unique identifier
df_design <- sf %>% mutate(fastq_path = str_c(run_accession, ".fastq.gz"))
#base library guides
FASTA <-  "TKOv1-base-90k-library-91320_sequences.fasta"
#supplementary library guides
 cb2_count_base <- run_sgrna_quant(FASTA, df_design)
#save.image(file = "cb2_1.Rdata")
load(file = "cb2_1.Rdata")
# supplementary library guides
FASTA <-  "TKOv1-supp-85k-library-85180_sequences.fasta"

cb2_count_sup <- run_sgrna_quant(FASTA, df_design)
# save.image(file = "cb2_2.Rdata")
# load pre-calculated counts from disk
load(file = "cb2_2.Rdata")
cb2_count_base_t <- as_tibble(cb2_count_base$count, rownames = "gene")
cb2_count_sup_t <- as_tibble(cb2_count_sup$count, rownames = "gene")
# reload data sample file as it is part of the stored R objects
sf <- sf %>% select(-sample) %>% rename(sample = sample_name) %>% select(-sample_accession, -run_accession)
# convert from wide into long format
cb2_count_base_g <- cb2_count_base_t %>% 
    gather(sample, value, `1%_A`:`all_F`)
cb2_count_sup_g <- cb2_count_sup_t %>% 
    gather(sample, value, `1%_A`:`all_F`)
#count total number of guides and plot as bar plots
#supplementary
cb2_count_sup_g %>% 
    group_by(sample) %>% 
    summarize(total = sum(value)) %>% ggplot(aes(sample, total)) + geom_col() + theme_light(base_size = 20)+ 
    theme(axis.text.x = element_text(angle = 90)) + ggtitle("supplementary")
ggsave("total_count_sup.png")

#base
cb2_count_base_g %>% group_by(sample) %>% summarize(total = sum(value)) %>% ggplot(aes(sample, total)) + geom_col()+ theme_light(base_size = 20)+ 
    theme(axis.text.x = element_text(angle = 90)) +ggtitle("base")
ggsave("total_count_base.png") 

#write counts to text files
cb2_count_base_g %>% group_by(sample) %>% summarize(total = sum(value)) %>% write_tsv("total_count_base.txt")
cb2_count_sup_g %>% group_by(sample) %>% summarize(total = sum(value)) %>% write_tsv("total_count_sup.txt")

#normalization by sample size spread and then gather again
#base
cb2_count_base_g <- cb2_count_base_g %>% 
    group_by(sample) %>% 
    mutate(total_guides = sum(value)) %>% 
    mutate(value = value/(total_guides/1000000)) %>% 
    select(-total_guides) %>% 
    spread(sample,value, fill = 0) %>% 
    gather(sample, value, `1%_A`:`all_F`)

#sup
cb2_count_sup_g <- cb2_count_sup_g %>% 
    group_by(sample) %>% 
    mutate(total_guides = sum(value)) %>% 
    mutate(value = value/(total_guides/1000000)) %>% 
    select(-total_guides) %>% spread(sample,value, fill = 0) %>% 
    gather(sample, value, `1%_A`:`all_F`)

#long version of guide counts
#base
cb2_count_base_t <- cb2_count_base_g %>% 
    group_by(sample) %>%
    mutate(total_guides = sum(value)) %>% 
    mutate(value = value/(total_guides/1000000)) %>% 
    select(-total_guides) %>% 
    spread(sample,value, fill = 0) 
#supplemenatry
cb2_count_sup_t <-cb2_count_sup_g %>% group_by(sample) %>% mutate(total_guides = sum(value)) %>% mutate(value = value/(total_guides/1000000)) %>% select(-total_guides) %>% spread(sample,value, fill = 0)

#make correlation plot
png("corrplot_base.png")
corrplot(cor(as.data.frame(cb2_count_base_t)[, -1:-2], method = "spearman"), order = "hclust", title="base")
dev.off()
png("corrplot_sup.png")
corrplot(cor(as.data.frame(cb2_count_sup_t)[, -1:-2], method = "spearman"), order = "hclust", title="supplementary")
dev.off()

#count distribution density
ggplot(cb2_count_base_g, aes(value)) + geom_density(adjust = 0.5) + facet_wrap(~sample, nrow = 3)  + scale_x_continuous(trans='log2') + theme_light(base_size = 20) + 
    theme(axis.text.x = element_text(angle = 90)) + ggtitle("base")
ggsave("count_distribution_base_gg.png", width = 20, height = 10)
ggplot(cb2_count_sup_g, aes(value)) + geom_density(adjust = 0.5) + facet_wrap(~sample, nrow = 3)  + scale_x_continuous(trans='log2') + theme_light(base_size = 20)+ 
    theme(axis.text.x = element_text(angle = 90)) + ggtitle("supplementary")
ggsave("count_distribution_sup_gg.png", width = 20, height = 10)


#add extra column distinguishing the library
cb2_count_base_t <- mutate(cb2_count_base_t, exp = "base")
cb2_count_sup_t <- mutate(cb2_count_sup_t, exp = "sup")
#merge base and supplementary library
cb2_count_t <- bind_rows(cb2_count_base_t, cb2_count_sup_t)

#exp A
cb2_count_t_cor <-mutate(cb2_count_t, all_comb = (all_A +  all_F + all_D)/3, `15%` = (`15%_A` +  `15%_F` + `15%_D`)/3, `5%` = (`5%_A` +  `5%_F` + `5%_D` )/3, `1%` = (`1%_A` +  `1%_F` + `1%_D`)/3)
#combined correlation plot
png("corrplot_combined.png")
corrplot(cor(as.data.frame(cb2_count_t_cor %>% select(-exp))[, -1], method = "spearman"), order = "hclust", title="base")
dev.off()

zero_counts <- cb2_count_sup_g %>% group_by(sample) %>% mutate(drop_out = ifelse(value == 0, T , F)) %>% count(drop_out)
ggplot(zero_counts, aes(sample, n,  fill = as.factor(drop_out))) + geom_col() + ylab("Guides") + labs(fill = "is_drop_out") + theme_light(base_size = 20)+ 
    theme(axis.text.x = element_text(angle = 90)) + ggtitle("supplementary")
ggsave("drop_out_sup.png") 
zero_counts <- cb2_count_base_g %>% group_by(sample) %>% mutate(drop_out = ifelse(value == 0, T , F)) %>% count(drop_out)
ggplot(zero_counts, aes(sample, n,  fill = as.factor(drop_out))) + geom_col() + ylab("Guides") + labs(fill = "is_drop_out")+ theme_light(base_size = 20)+ 
    theme(axis.text.x = element_text(angle = 90)) + ggtitle("base")
ggsave("drop_out_base.png") 

#join base and sup libraries
cb2_count <- bind_rows(cb2_count_base_g %>% 
                       mutate(library = "base"), cb2_count_sup_g %>% 
                       mutate(library = "sup"))
 

#pooling experiments 
#bright only
bright_count_pooled <-  cb2_count %>% 
    inner_join(sf) %>% 
    filter(experiment == "A" | experiment == "D" | experiment == "F") %>% 
    select(-sample, -index, -hiseq) %>% 
    spread(group, value)  %>% 
    group_by(experiment) %>%  
    mutate(bright = `1%` + `5%` + `15%`) %>% 
    group_by(experiment) %>% mutate(all = all/sum(all) * 1000000, bright = bright/sum(bright) * 1000000 )%>% 
    group_by(gene, library) %>% 
    select(-`1%`, -`5%`, - `15%`, -all) %>% 
    spread(experiment, bright) 

#general population only
all_pooled <-cb2_count %>% 
    inner_join(sf) %>% 
    filter(experiment == "A" | experiment == "D" | experiment == "F") %>% 
    select(-sample, -index, -hiseq) %>% 
    spread(group, value)  %>% 
    group_by(experiment) %>% 
    mutate(all = all/sum(all) * 1000000 )%>% 
    group_by(gene, library)  %>% 
    summarize( all = max(all)) %>% 
    ungroup() %>% 
    mutate(all = all/sum(all) * 1000000 )

#all combined drop out 
all_pooled %>% mutate(is_dropout = ifelse(all == 0, "Yes" , "No")) %>% group_by(is_dropout, library) %>% ggplot(aes(library,  fill = is_dropout)) + geom_bar()
ggsave("drop_out_allT.png")

# density plots)
all_pooled %>% ggplot(aes(all)) + geom_density() + scale_x_continuous(trans='log2') 
ggsave("all_pooled_density.png")
bright_count_pooled %>% ggplot(aes(A)) + geom_density() + scale_x_continuous(trans='log2') 
ggsave("A_pooled_density.png")
bright_count_pooled %>% ggplot(aes(D)) + geom_density() + scale_x_continuous(trans='log2') 

ggsave("D_pooled_density.png")
bright_count_pooled %>% ggplot(aes(F)) + geom_density() + scale_x_continuous(trans='log2') 
ggsave("F_pooled_density.png")

all_pooled %>% mutate(gene_perm = sample_frac(all_pooled, size = 1)$gene)
#filter by abundance in unsorted (all)
all_pooled_filt <- filter(all_pooled, all > 0.1) %>%  unite(gene_library, c("gene", "library"))

# different pooling functions
bright_fc <- inner_join(bright_count_pooled, all_pooled) %>% 
    filter(all > .1) %>% 
    mutate(fc_A = A/all, fc_D = D/all, fc_F = F/all)  %>% 
    rowwise() %>% 
    mutate(fc_min = pmin(fc_A, fc_D, fc_F)) %>% 
    mutate(fc_max = pmax(fc_A, fc_D, fc_F)) %>% 
    mutate(fc_mean = mean(c(fc_A, fc_D, fc_F))) %>% 
    mutate(fc_median = median(c(fc_A, fc_D, fc_F))) 

#permute for every exerpiment, removing drop outs
perm1 <- all_pooled_filt %>% mutate(gene_library = sample_frac(all_pooled_filt, size = 1)$gene_library) %>% rename(all_perm1 = all) %>% 
    mutate(dropout = c(rep(FALSE, nrow(bright_fc %>% filter(A > 0))), rep(TRUE, nrow(bright_fc %>% filter(A == 0))))) %>% 
    mutate(all_perm1 = ifelse(dropout, 0, all_perm1)) %>%
    select(-dropout)
perm2 <- all_pooled_filt %>% mutate(gene_library = sample_frac(all_pooled_filt, size = 1)$gene_library) %>% rename(all_perm2 = all) %>% 
    mutate(dropout = c(rep(FALSE, nrow(bright_fc %>% filter(F > 0))), rep(TRUE, nrow(bright_fc %>% filter(F == 0))))) %>% 
    mutate(all_perm2 = ifelse(dropout, 0,  all_perm2)) %>% 
    select(-dropout)
perm3 <- all_pooled_filt %>% mutate(gene_library = sample_frac(all_pooled_filt, size = 1)$gene_library) %>% rename(all_perm3 = all) %>% 
    mutate(dropout = c(rep(FALSE, nrow(bright_fc %>% filter(D > 0))), rep(TRUE, nrow(bright_fc %>% filter(D == 0))))) %>% 
    mutate(all_perm3 = ifelse(dropout, 0, all_perm3)) %>%
    select(-dropout)

# combine permuted experiments A, F, D
perm_fc <- inner_join(all_pooled_filt, perm1) %>%  inner_join(perm2) %>% inner_join(perm3) %>%
    mutate(fc_perm1 = all_perm1/all, fc_perm2 = all_perm2/all, fc_perm3 = all_perm3/all)  %>% 
    rowwise() %>% 
    mutate(fc_min = pmin(fc_perm1, fc_perm2, fc_perm3)) %>% 
    mutate(fc_max = pmax(fc_perm1, fc_perm2, fc_perm3)) %>% 
    mutate(fc_mean = mean(c(fc_perm1, fc_perm2, fc_perm3))) %>% 
    mutate(fc_median = median(c(fc_perm1, fc_perm2, fc_perm3))) 

perm_fc <- perm_fc %>% arrange(-fc_median)
bright_fc <- bright_fc %>%   arrange(-fc_median)

rank = nrow(perm_fc)
new_rank <- vector(mode = "list", length = nrow(bright_fc))
i = nrow(bright_fc) 
while(i > 0)
{
    while(rank > 1 && perm_fc$fc_median[rank] <= bright_fc$fc_min[i])
          {rank = rank - 1}    
    new_rank[[i]] = rank
    i = i - 1
}
#ascribe pvalue based on permutation based approach
bright_fc <- bright_fc %>% 
    ungroup %>% 
    mutate(rank = unlist(new_rank)) %>% 
    mutate(pval = rank/nrow(perm_fc))

bright_fc %>% select(-rank) %>% write_tsv("guides_w_pval_fc_min.txt")

