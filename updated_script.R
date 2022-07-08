library(tidyverse)
library(glue)
library(gganimate)

processed_mapping <- read.delim("processed_Cry5B_inbred.tsv", stringsAsFactors=FALSE) %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA"))) %>%
  dplyr::select(-marker) %>%
  tidyr::unite("marker", CHROM, POS, sep = ":", remove = F) %>%
  dplyr::mutate(POS = POS/1000000) %>% 
  dplyr::select(marker:aboveBF)

genotype_matrix <- read.delim("Genotype_Matrix.tsv", stringsAsFactors=FALSE) %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA"))) %>%
  tidyr::unite("marker", CHROM, POS, sep = ":", remove = F) %>%
  dplyr::mutate(POS = POS/1000000) %>%
  tidyr::pivot_longer(cols = AB1:XZ2018,
                      names_to = "strain",
                      values_to = "allele") %>%
  dplyr::mutate(allele = dplyr::case_when(allele == "-1" ~ "REF",
                                          allele == "1" ~ "ALT",
                                          TRUE ~ "NA"),
                allele = factor(allele, levels = c("REF", "ALT")))

phenotypes <- read.delim("phenotypes.tsv", stringsAsFactors=FALSE) %>% 
  dplyr::rename(value = Cry5B) %>%
  dplyr::left_join(genotype_matrix, ., by = "strain")

## MANHATTAN PLOTS ##

BF <- processed_mapping %>% 
  dplyr::group_by(trait) %>% 
  dplyr::filter(log10p != 0) %>% 
  dplyr::distinct(marker, log10p) %>%
  dplyr::mutate(BF = -log10(0.05/sum(log10p > 0, na.rm = T))) %>%
  dplyr::ungroup() %>%
  dplyr::select(BF) %>%
  unique(.) %>%
  dplyr::slice(1) %>% # BF can be slightly different between loco and inbred... but just plot one (5.46 v 5.47...)
  as.numeric()

BF.frame <- processed_mapping %>%
  dplyr::select(trait) %>%
  dplyr::filter(!duplicated(trait)) %>%
  dplyr::mutate(BF = BF)


all_data <- processed_mapping %>%
  dplyr::inner_join(., phenotypes, by = c("marker", "CHROM", "POS")) %>%
  dplyr::mutate(sig = case_when(log10p > BF.frame$BF ~ "BF",
                                TRUE ~ "NONSIG"))

sig.colors <- c("red","black")
names(sig.colors) <- c("BF","NONSIG")
facet_scales <- "free"


all_data %>%
  dplyr::select(-c(REF:value)) %>%
  unique() %>%
  dplyr::filter(CHROM %in% c("I","II","IV")) %>%
  ggplot() + 
  theme_bw() + 
  geom_point(mapping = aes(x = POS, 
                           y = log10p,
                           colour = sig,
                           alpha = sig)) +
  scale_alpha_manual(values = c("BF" = 1,  "user" = 1, "NONSIG" = 0.25)) +
  scale_colour_manual(values = sig.colors) + 
  scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
  scale_y_continuous(limits = c(0,7.5))+
  geom_hline(yintercept = 5.47) + 
  scale_linetype_manual(values = c("BF" = 1)) +
  labs(x = "Genomic position (Mb)",
       y = expression(-log[10](italic(p)))) +
  theme(legend.position = "none", panel.grid = element_blank()) + 
  facet_grid(. ~ CHROM, scales = "free_x", space = facet_scales) + 
  #ggtitle(BF.frame$trait) +
  transition_reveal(POS) +
  shadow_trail(future = F, alpha = alpha*2, size = size/10, distance = 0.01)

animate(anim, renderer = gifski_renderer())
anim_save(anim,"anim.gif")



all_data %>%
  #dplyr::filter(CHROM %in% c("I")) %>%
  dplyr::filter(marker %in% c("I:13117")) %>%
  ggplot() +
  theme_bw() +
  aes(x = allele, y = value) +
  guides(color = "none") +
  geom_point() +
  geom_violin(alpha = 0.4, width = 0.5) +
  theme(legend.position = "bottom") +
  labs(x = "Genotype", y = "Trait Value", title = "Marker: {frame_along}") +
  transition_reveal(POS)













