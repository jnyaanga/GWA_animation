library(tidyverse)
library(glue)
library(gganimate)
library(ggbeeswarm)

######################
####  Load Data   ####
######################

# load data for manhattan plot --> ie. marker information
processed_mapping <- read.delim("processed_Cry5B_inbred.tsv", stringsAsFactors=FALSE) %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA"))) %>%
  dplyr::select(-marker) %>%
  tidyr::unite("marker", CHROM, POS, sep = ":", remove = F) %>%
  dplyr::mutate(POS = POS/1000000) %>% 
  dplyr::select(marker:aboveBF) %>%
  dplyr::distinct()

# load strain data for each marker 
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

# load phenotype data for each strain
phenotypes <- read.delim("phenotypes.tsv", stringsAsFactors=FALSE) %>% 
  dplyr::rename(value = Cry5B) %>%
  dplyr::left_join(genotype_matrix, ., by = "strain")


#####################################
####  Clean up data and combine  ####
#####################################

# clean up marker data and set significance threshold as BF
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

# combine all data
all_data <- processed_mapping %>%
  dplyr::inner_join(., phenotypes, by = c("marker", "CHROM", "POS")) %>%
  dplyr::mutate(sig = case_when(log10p > BF.frame$BF ~ "BF", TRUE ~ "NONSIG"),
                strain = as.factor(strain))

sig.colors <- c("red","black")
names(sig.colors) <- c("BF","NONSIG")
facet_scales <- "free"

################
####  PLOT  ####
################

## notes: 
## I played around with the animations a little to make sure the data was in a useable format. Feel free to update with what you found worked (transition_manual etc)

## Overall, I think the best way to animate the whole genome is probably to make individual plots for each chromosome and then merge them together later using
## https://imagemagick.org/index.php or something similar? We'll have to think about it. It just looks like it'll be too difficult to do all chromosomes in one go


# manhattan
chrI <- all_data %>%
  dplyr::select(-c(REF:value)) %>%
  unique() %>%
  dplyr::filter(CHROM %in% c("I")) %>%
  ggplot() + 
  theme_bw() + 
  aes(group = POS) +
  geom_point(mapping = aes(x = POS, y = log10p, color = sig, alpha = sig)) +
  scale_alpha_manual(values = c("BF" = 1,  "user" = 1, "NONSIG" = 0.25)) +
  scale_colour_manual(values = sig.colors) + 
  scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
  scale_y_continuous(limits = c(0,7.5)) +
  geom_hline(yintercept = 5.47) + 
  scale_linetype_manual(values = c("BF" = 1)) +
  labs(x = "", y = expression(-log[10](italic(p))), title = "Progress: {frame} of {nframes}") +
  theme(legend.position = "none", panel.grid = element_blank()) + 
  facet_grid(. ~ CHROM, scales = "free_x", space = facet_scales) + 
  transition_reveal(POS) +
  shadow_trail(future = F, alpha = alpha*2, size = size/10, distance = 0.01) 

animate(chrI, nframes=200, fps=10, duration = 20, renderer = gifski_renderer(loop = F, width = 800, height = 1200))
#anim_save("ChrI_man.gif")


# phenotype x genotype
chrIpxg <- all_data %>%
  dplyr::filter(CHROM %in% c("I")) %>%
  ggplot() +
  theme_bw() +
  aes(x = allele, y = value) +
  geom_point(aes(color = sig), position = ggbeeswarm::position_beeswarm(), shape = 21) +
  geom_boxplot(aes(color = sig), alpha = 0.4) +
  scale_color_manual(values = sig.colors) + 
  guides(color = "none") +
  theme(legend.position = "bottom") +
  labs(x = "Genotype", y = "Trait Value", title = "Progress: {frame} of {nframes}") +
  facet_grid(. ~ CHROM) +
  scale_y_continuous(limits = c(-3.5,3.5)) +
  scale_x_discrete(breaks = c("REF", "ALT")) +
  transition_manual(POS) 


animate(chrIpxg, nframes=200, duration = 20, fps=10, renderer = gifski_renderer(loop = F, width = 400, height = 400))
#settings will need to change for chrIV --> potentially split into 2 animations so that it delays at the peak markers (~frames 850-950)
#anim_save("ChrI_pxg.gif")


## MISC
# to see metadata about frames in an animation use: frame_vars.


