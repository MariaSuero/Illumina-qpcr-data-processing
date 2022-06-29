library(devtools)
devtools::install_github("ewallace/tidyqpcr",build_vignettes = TRUE) ## Vignettes require cowplot package
library(tidyqpcr)

library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyqpcr)
install.packages("tidyverse")
library(tidyr)
theme_set(theme_bw(base_size = 11) %+replace%
            theme(
              strip.background = element_blank(),
              panel.grid = element_blank()
            ))

# Names of target genes
gene_name_values <- c("B-actin", "IL1B", "IL18")
print(gene_name_values)

# Repeats of gene names to account for testing multiple primer sets
target_id_levels <- c(rep(gene_name_values[1:3], each = 2))
print(target_id_levels)


rowkey1 <- tibble(well_row = LETTERS[1:6],
                 target_id = factor(target_id_levels, levels = gene_name_values))
print(rowkey1)

#### AQUI ES COMO LO QUIERO POR AHORA BIEN

sample_id_levels <- c("CTRL","lob1","lob2","Artem","lob+artI","lob+artI", "lob+artR", "NT")

prep_type_levels <- ""

print(sample_id_levels)
colkey4 <- tibble(
  well_col = 1:8,
  sample_id = sample_id_levels,
  prep_type = prep_type_levels
)
print(colkey4)

### vale todo good he puesto los samples 
###pone para retener la info hay que meter esto lol 

create_blank_plate(well_row = LETTERS[1:6], well_col = 1:8)

###ahora uno las cosas a ver si sirve 

plan_pyropt <- label_plate_rowcol(plate = create_blank_plate(well_row = LETTERS[1:6], well_col = 1:8),
                                       rowkey = rowkey1,
                                       colkey = colkey4)

print(plan_pyropt)

##alucinante creo que todo bien 

display_plate_qpcr(plan_pyropt)


columns <- c("well","cq")
file_path_cq <- read.table(file = "pyro primers nuevos 1.txt", header = TRUE,  sep = "", dec = ",", col.names = columns)
print(file_path_cq )

plates = merge(x = file_path_cq, y = plan_pyropt, by = "well")
plates
##parece que va pq lo lee 


cq_plate <- display_plate_value(plates)
cq_plate

## All reps, unnormalized
ggplot(data = plates) +
  geom_point(aes(x = sample_id, y = cq, shape = prep_type, colour = prep_type_levels),
             position = position_jitter(width = 0.2, height = 0)
  ) +
  labs(
    y = "Quantification cycle (Cq)",
    title = "All reps, unnormalized"
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


sapply(plates, class)


plates$cq <- as.numeric(gsub("," ,".", as.character(plates$cq)))  # Convert one variable to numeric
plates

plate_norm <- plates %>%
  calculate_deltacq_bysampleid(ref_target_ids = "B-actin") 




###grafica de todos los datos

ggplot(data = plate_norm) +
  geom_point(aes(x = sample_id, y = delta_cq, shape = prep_type, colour = prep_type),
             position = position_jitter(width = 0.2, height = 0)
  ) +
  labs(y = "delta Cq") +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  


## Normalized values heatmap 
HeatmapDELTAcq <- display_plate_value(plate_norm)
display_plate_value(plate_norm)



### calcular las medias
plate_med <- plate_norm %>%
  group_by(sample_id, prep_type, target_id) %>%
  summarize(
    delta_cq  = median(delta_cq, na.rm = TRUE),
      rel_abund = median(rel_abund, na.rm = FALSE)
  )
  
plate_med
View(plate_med)

plates_SD = merge(x = plate_med, y = standard_dev, by = "sample_id")

plates_SD
view(plates_SD)
### delta cqs de todo 

ggplot(data = filter(plate_med, target_id != "B-actin") ) + 
  geom_col(aes(x = sample_id, y = delta_cq, shape = prep_type, fill = "blue") ) + 
  labs(
    y = "delta Cq",
    x = "treatment") +
  
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



## grafica de barras

ggplot(data = filter(plate_med, target_id != "B-actin") ) + 
  geom_col(aes(x = sample_id, y = delta_cq), fill="red" ) + 
  labs(
    y = "delta Cq",
    x = "Treatment"
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


### OTRO DF 
ggplot(data = filter(plates_SD, target_id.x != "B-actin") ) + 
  geom_col(aes(x = sample_id, y = delta_cq), fill="steelblue" ) + 
  labs(
    y = "delta Cq",
    x = "Treatment"
  ) +
  
  geom_errorbar(aes(x = sample_id, y = delta_cq, ymin=delta_cq-sd, ymax=delta_cq+sd), width=.5,
                position = "identity", size = .5,   ) +
  facet_wrap(~target_id.x,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



### con sd y solo alguno datos

ggplot(data = filter(plates_SD, target_id.x != "B-actin", sample_id %in% c("Artem", "CTRL", "lob1", "lob+artI"))) + 
  geom_col(aes(x = sample_id, y = delta_cq), fill="lightblue" ) + 
  labs(
    y = "delta Cq",
    x = "Treatment"
  ) +
  
  geom_errorbar(aes(x = sample_id, y = delta_cq, ymin= delta_cq-sd, ymax= delta_cq+sd), width=.2,
                position = position_dodge(0.05), size = .5) +

  
  facet_wrap(~target_id.x,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

###solo lob y control con sd
ggplot(data = filter(plates_SD, target_id.x != "B-actin", sample_id %in% c("CTRL", "lob1"))) + 
  geom_col(aes(x = sample_id, y = delta_cq), fill="lightgreen" ) + 
  labs(
    y = "delta Cq",
    x = "Treatment"
  ) +
  
  geom_errorbar(aes(x = sample_id, y = delta_cq, ymin= delta_cq-sd, ymax= delta_cq+sd), width=.2,
                position = position_dodge(0.5), size = .5) +
  
  
  facet_wrap(~target_id.x,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


####excluir datos


ggplot(subset(plate_med, sample_id %in% c("Artem", "CTRL", "lob1", "lob+artI")) , 
       aes(x = sample_id, y = delta_cq))+
  geom_col(aes(x = sample_id, y = delta_cq, shape = prep_type, fill = "blue") ) +
  labs(
    y = "ΔCq",
    x = "Treatment"
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))





## abundancia RNA 
plate_med
ggplot(data = filter(plate_med, target_id != "B-actin") ) + 
  geom_point(aes(x = sample_id, y = rel_abund), color="red") + 
  labs(
    y = "RNA abundance relative to B-actin",
    x = "Treatment"
    
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



### barras

plate_med
ggplot(data = filter(plate_med, target_id != "B-actin") ) + 
  geom_col(aes(x = sample_id, y = rel_abund), fill="red") + 
  labs(
    y = "RNA abundance relative to B-actin",
    x = "Treatment"
    
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

##un intento tonto la vd

ggplot(data = filter(plate_med, target_id != "B-actin") ) + 
  geom_col(aes(x = sample_id, y = rel_abund), fill="red") + 
  labs(
    y = "RNA abundance relative to B-actin",
    x = "Treatment"
    
  ) +
  facet_wrap(~target_id,ncol=4) + 
  geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund-0.004, ymax= rel_abund+0.004), width=.2,
                position = position_dodge(0.5), size = .5) +
  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


###exclusion datos
 ggplot(data = filter(plate_med, target_id != "B-actin", sample_id %in% c("CTRL", "lob1" ))) + 
   scale_x_discrete(breaks = c('CTRL', 'lob1'),
                     labels = c("Control", "Lobaplatin")) +
                   geom_col(aes(x = sample_id, y = rel_abund), fill="orange") + 
                   labs(
                     y = "RNA abundance relative to B-actin",
                     x = "Treatment"
                     ) + geom_text(aes(x = sample_id, y = rel_abund, label = "*"), vjust = -5)+
                   facet_wrap(~target_id,ncol=4) + 
                   geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund-0.004, ymax= rel_abund+0.004), width=.2,
                                 position = position_dodge(0.5), size = .5) +
                   
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))                 
                 
 
 
 
 ####desviacion estandar
 standard_dev <-  plate_norm %>% group_by(sample_id,target_id) %>%
   summarise(sd = sd(rel_abund))
 standard_dev
 

pyro_filtered <- subset(standard_dev, target_id != "B-actin")

pyro_filtered
pyro_art_sd <- subset(pyro_filtered, sample_id %in% c("CTRL", "lob2", "lob+artI", "Artem"))
 
 pyro_art_sd
 
 #### FILTRO SD
 filtered_il1b_sd <- subset(pyro_art_sd, target_id != "IL18")
 filtered_il18_sd <- subset(pyro_art_sd, target_id != "IL1B")
 
 
 ###artemisina IL1B
 
 ggplot(data = filter(plate_med, target_id != "B-actin",target_id != "IL18",  sample_id %in% c("CTRL", "lob2", "lob+artI", "Artem"
 ))) + 
   geom_col(aes(x = factor(sample_id, level = c('CTRL', 'lob2', 'Artem', "lob+artI")), y = rel_abund), fill="red") +  geom_text(aes(x = sample_id, y = rel_abund, label = "***"), vjust = -2)+
   labs(
     y = "RNA abundance relative to B-actin",
     x = "Treatment"
     
   ) + scale_x_discrete(breaks = c('CTRL', 'lob2', 'Artem', "lob+artI"),
                        labels = c("Control", "Lobaplatin", "Artemisinin", "Co-treatment"))+
   facet_wrap(~target_id,ncol=4) + 
   geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund- filtered_il1b_sd$sd, ymax= rel_abund+ filtered_il1b_sd$sd), width=.2,
                 position = position_dodge(0.5), size = .5) +
   
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  
 
 ###artemisina IL18
 ggplot(data = filter(plate_med, target_id != "B-actin", target_id != "IL1B", sample_id %in% c("CTRL", "lob2", "lob+artI", "Artem"
 ))) + 
   geom_col(aes(x = factor(sample_id, level = c('CTRL', 'lob2', 'Artem', "lob+artI")), y = rel_abund), fill="red") +  geom_text(aes(x = sample_id, y = rel_abund, label = "***"), vjust = -2) +
   labs(
     y = "RNA abundance relative to B-actin",
     x = "Treatment"
     
   ) + scale_x_discrete(breaks = c('CTRL', 'lob2', 'Artem', "lob+artI"),
                        labels = c("Control", "Lobaplatin", "Artemisinin", "Co-treatment"))+
   facet_wrap(~target_id,ncol=4) + 
   geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund- filtered_il18_sd$sd, ymax= rel_abund+ filtered_il18_sd$sd), width=.2,
                 position = position_dodge(0.5), size = .5) +
   
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  
 

 
 ###artemisina con ambos biomarkers
 
 ggplot(data = filter(plate_med, target_id != "B-actin", sample_id %in% c("CTRL", "lob2", "lob+artI", "Artem"
 ))) + 
   geom_col(aes(x = factor(sample_id, level = c('CTRL', 'lob2', 'Artem', "lob+artI")), y = rel_abund), fill="red") + 
   labs(y = "RNA abundance relative to B-actin", x = "Treatment")+ scale_x_discrete(breaks = c('CTRL', 'lob2', 'Artem', "lob+artI"),
                          labels = c("Control", "Lobaplatin", "Artemisinin", "Co-treatment")) +
   facet_wrap(~target_id,ncol=4) +  geom_text(aes(x = sample_id, y = rel_abund, label = "***"), vjust = -2)+
   geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund- pyro_art_sd$sd, ymax= rel_abund+ pyro_art_sd$sd), width=.2,
                 position = position_dodge(0.5), size = .5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  
 
 
 
 ### Añadir ANOVA 
 
 analysis <- aov(data = plate_norm, rel_abund ~ sample_id)
 analysis
 summary(analysis)
 TukeyHSD(analysis)
 
 pyro_il1b_ANOVA <- subset(plate_norm, target_id  != "B-actin", target_id != "IL18",plate_norm$sample_id  %in% c("CTRL", "lob2", "lob+artI", "Artem"))
 
 pyro_il18_ANOVA <- subset(plate_norm, target_id  != "B-actin", target_id != "IL1B", plate_norm$sample_id %in% c("CTRL", "lob2", "lob+artI", "Artem"))
  
 analysisil18 <- aov(data = pyro_il18_ANOVA, rel_abund ~ sample_id)
 analysisil1b <- aov(data = pyro_il1b_ANOVA, rel_abund ~ sample_id)
 
 summary(analysisil1b)
 summary(analysisil18)
 
TukeyHSD(analysisil1b)