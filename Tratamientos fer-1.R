library(devtools)
devtools::install_github("ewallace/tidyqpcr",build_vignettes = TRUE) ## Vignettes require cowplot package
library(tidyqpcr)

library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyqpcr)
theme_set(theme_bw(base_size = 11) %+replace%
            theme(
              strip.background = element_blank(),
              panel.grid = element_blank()
            ))
# Names of target genes
target_id_levels <- c("B-actin", "CHAC1")
print(target_id_levels)

# Repeats of gene names to account for testing multiple primer sets
target_id_levels <- c(rep(gene_name_levels[1:2], each = 2))

rowkey <- tibble(well_row = LETTERS[1:4],
                 target_id = factor(target_id_levels, levels = gene_name_levels))
print(rowkey)
print(target_id_levels)

#### AQUI ES COMO LO QUIERO POR AHORA BIEN

sample_id_levels <- c("CTRL", "Erastin", "Ferrostatin-1", "Co-treatment", "NT")


colkey6 <- tibble(
  well_col = 1:5,
  sample_id = sample_id_levels,
  prep_type = " "
  
)
print(colkey6)

### vale todo good he puesto los samples 
###pone para retener la info hay que meter esto lol 

create_blank_plate(well_row = LETTERS[1:4], well_col = 1:5)

###ahora uno las cosas a ver si sirve 

plan_inhibicion_r1 <- label_plate_rowcol(plate = create_blank_plate(well_row = LETTERS[1:4], well_col = 1:5),
                                         rowkey = rowkey,
                                         colkey = colkey6)


print(plan_inhibicion_r1)

##alucinante creo que todo bien 

display_plate_qpcr(plan_inhibicion_r1)


columns <- c("well", "cq")
file_path_cq <- read.table(file = "ferroptosis con I +R copia 2.txt", header = FALSE,  sep = "", dec = ",", col.names = columns)
print(file_path_cq )
print(plan_inhibicion_r1)

plates1 = merge(x = file_path_cq, y = plan_inhibicion_r1, by = "well",
)
plates1
##parece que va pq lo lee 


display_plate_value(plates1)

ggplot(data = plates1) +
  geom_point(aes(x = sample_id, y = cq, shape = prep_type, colour = prep_type_levels),
             position = position_jitter(width = 0.2, height = 0)
  ) +
  labs(
    y = "Quantification cycle (Cq)",
    title = "All reps, unnormalized"
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


sapply(plates1, class)


plates1$cq <- as.numeric(gsub("," ,".", as.character(plates1$cq)))  # Convert one variable to numeric
plates1

plate_norm <- plates1 %>%
  calculate_deltacq_bysampleid(ref_target_ids = "B-actin") 


plate_norm

ggplot(data = plate_norm) +
  geom_point(aes(x = sample_id, y = delta_cq, shape = prep_type, colour = prep_type),
             position = position_jitter(width = 0.2, height = 0)
  ) +
  labs(y = "delta Cq") +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))             

### mirar si me deja hacer un heatmap

display_plate_value(plate_norm, value = "cq") +   # uses ggplot syntax
  scale_fill_gradient(high = "steelblue") 

### calcular abundancia relativa 

plate_med <- plate_norm %>%
  group_by(sample_id, prep_type, target_id) %>%
  summarize(
    delta_cq  = median(delta_cq, na.rm = TRUE),
    rel_abund = median(rel_abund, na.rm = TRUE)
  )

plate_med


ggplot(data = filter(plate_med, target_id != "B-actin") ) + 
  geom_point(aes(x = sample_id, y = delta_cq, shape = prep_type, colour = prep_type) ) + 
  labs(
    y = "delta Cq"
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


## ESTETICAMENTE ASI QUEDA MEJOR: 

ggplot(data = filter(plate_med, target_id != "B-actin") ) + 
  geom_point(aes(x = sample_id, y = delta_cq), color="red" ) + 
  labs(
    y = "delta Cq",
    x = "[Erastin]"
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



## abundancia RNA 

ggplot(data = filter(plate_med, target_id != "B-actin") ) + 
  geom_point(aes(x = sample_id, y = rel_abund), color="blue") + 
  labs(
    y = "RNA abundance relative to B-actin",
    x = "After 24h Treatment"
    
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



### grafica de barras pq va a quedar mejor para datos discretos

standard_dev <-  plate_norm %>% group_by(sample_id,target_id) %>%
  summarise(sd = sd(rel_abund))

chac1_sd <- subset(standard_dev, target_id  %in% "CHAC1")

p <- ggplot(filter(plate_med, target_id != "B-actin"), aes(x=sample_id, y=rel_abund)) + 
  geom_bar(stat="identity", fill = "grey", 
           position=position_dodge()) +
  geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund-chac1_sd$sd, ymax= rel_abund+chac1_sd$sd), width=.2,
                position = position_dodge(0.5), size = .5)

p +labs(title="CHAC1 RNA abundance relative to B-actin", x="Treatment", y = "Relative RNA abundance")


