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
install.packages("tidyverse")
library(tidyr)

# Names of target genes
gene_name_values <- c("B-actin", "BCL2L2", "LMNB1")
print(gene_name_values)

# Repeats of gene names to account for testing multiple primer sets
target_id_levels <- c(rep(gene_name_values[1:3], each = 2))

rowkey <- tibble(well_row = LETTERS[1:6],
                 target_id = factor(target_id_levels, levels = gene_name_values))
print(rowkey)
print(target_id_levels)

#### Añadir nombres de las muestras

sample_id_levels <- c("CTRL","Rsv+querR","Querc","Rsv+QuerI1","Rsv+QuerI2","QuerRctrol","NT")
prep_type_levels <- "+Rsv"
print(sample_id_levels)
colkey4 <- tibble(
  well_col = 1:7,
  sample_id = sample_id_levels,
  prep_type = prep_type_levels
)
print(colkey4)

###Crear la placa

create_blank_plate(well_row = LETTERS[1:6], well_col = 1:7)

###Asignar valores a la tabla

plan_senescencia <- label_plate_rowcol(plate = create_blank_plate(well_row = LETTERS[1:6], well_col = 1:7),
                                          rowkey = rowkey,
                                          colkey = colkey4)

print(plan_senescencia)


display_plate_qpcr(plan_senescencia)


### leer resultados illumina 
columns <- c("well","cq")
file_path_cq <- read.table(file = "senescencia-3gen.txt", header = TRUE,  sep = "", dec = ",", col.names = columns)
print(file_path_cq )

plates = merge(x = file_path_cq, y = plan_senescencia, by = "well")
plates


display_plate_value(plates)

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


plate_norm

### visualizar los deltacq
ggplot(data = plate_norm) +
  geom_point(aes(x = sample_id, y = delta_cq, shape = prep_type, colour = prep_type),
             position = position_jitter(width = 0.2, height = 0)
  ) +
  labs(y = "delta Cq") +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))             


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




###reduciendolo a un solo marcador 

ggplot(data = filter(plate_med, target_id != "B-actin") ) + 
  geom_point(aes(x = sample_id, y = delta_cq), color="red" ) + 
  labs(
    y = "delta Cq",
    x = "Treatment"
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))




## abundancia RNA solo de bcl2l2

ggplot(data = filter(plate_med, target_id !=  "LMNB1", target_id !=  "B-actin") ) + 
  geom_col(aes(x = sample_id, y = rel_abund), fill="lightblue") + 
  labs(
    y = "RNA abundance relative to B-actin",
    x = "Treatment"
    
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

### abundancia RNA solo de bcl2l2

ggplot(data = filter(plate_med, target_id !=  "BCL2L2", target_id !=  "B-actin") ) + 
  geom_col(aes(x = sample_id, y = rel_abund), fill="lightgreen") + 
  labs(
    y = "RNA abundance relative to B-actin",
    x = "Treatment"
    
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


####desv standar 

standard_dev <-  plate_norm %>% group_by(sample_id,target_id) %>%
  summarise(sd = sd(rel_abund))



###subsets 

bcl2_sd <- subset(standard_dev, target_id  %in% "BCL2L2")

lmn_sd <- subset(standard_dev, target_id  %in% "LMNB1")

lmn_sd
###añadiendo barras en abundancia RNA solo de bcl2l2
ggplot(data = filter(plate_med, target_id !=  "LMNB1", target_id !=  "B-actin") ) + 
  geom_col(aes(x = sample_id, y = rel_abund), fill="lightblue") + 
  labs(
    y = "RNA abundance relative to B-actin",
    x = "Treatment",
    title = "BCL2L2") +
  geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund-bcl2_sd$sd, ymax= rel_abund+bcl2_sd$sd), width=.2,
                position = position_dodge(0.5), size = .5)
+ facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

###añadiendo barras en abundancia RNA solo de lmnb1

ggplot(data = filter(plate_med, target_id !=  "BCL2L2", target_id !=  "B-actin") ) + 
  geom_col(aes(x = sample_id, y = rel_abund), fill="lightblue") + 
  labs(
    y = "RNA abundance relative to B-actin",
    x = "Treatment",
    title = "LMNB1") +
  geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund-lmn_sd$sd, ymax= rel_abund+lmn_sd$sd), width=.2,
                position = position_dodge(0.5), size = .5)
+ facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

###excluir sd 

#### exclusion de datos para abund lmnb1
filtered_lmn_sd <- subset(lmn_sd, sample_id %in% c("CTRL", "Querc", "Rsv+QuerI2", "Rsv+querR"))
filtered_lmn_sd

ggplot(data = filter(plate_med, target_id !=  "BCL2L2", target_id !=  "B-actin", 
                     sample_id %in% c("CTRL", "Querc", "Rsv+QuerI2", "Rsv+querR" ))) + 
  geom_col(aes(x = sample_id, y = rel_abund), fill="lightblue") + geom_text(aes(x = sample_id, y = rel_abund, label = "***"), vjust = -2)+
  labs(
    y = "RNA abundance relative to B-actin",
    x = "Treatment",
    title = "LMNB1"
    ) +
  scale_x_discrete(breaks = c("CTRL", "Querc", "Rsv+QuerI2", "Rsv+querR"),
                   labels = c("Control", "Resveratrol","Co-treatment", "Quercetin")) +
geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund-filtered_lmn_sd$sd, ymax= rel_abund+filtered_lmn_sd$sd), width=.2,
                position = position_dodge(0.5), size = .5) + 
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))



##### exclusion de datos para abund bcl2

filtered_bcl2_sd <- subset(bcl2_sd, sample_id %in% c("CTRL", "Querc", "Rsv+QuerI2", "Rsv+querR"))
filtered_bcl2_sd


ggplot(data = filter(plate_med, target_id !=  "LMNB1", target_id !=  "B-actin", 
                     sample_id %in% c("CTRL", "Querc", "Rsv+QuerI2", "Rsv+querR" ))) + 
  geom_col(aes(x = sample_id, y = rel_abund), fill="steelblue") + geom_text(aes(x = sample_id, y = rel_abund, label = "***"), vjust = -2)+
  labs(
    y = "RNA abundance relative to B-actin",
    x = "Treatment",
    title = "BC2L2"
  ) +
  scale_x_discrete(breaks = c("CTRL", "Querc", "Rsv+QuerI2", "Rsv+querR"),
                   labels = c("Control", "Resveratrol","Co-treatment", "Quercetin"))+
  geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund-filtered_bcl2_sd$sd, ymax= rel_abund+filtered_bcl2_sd$sd), width=.2,
                position = position_dodge(0.5), size = .5) + 
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


### solo resveratrol para laminb1 -> USAR SD FILTRADA PARA CORREGIR

only_resveratrol_lmnb1<- ggplot(data = filter(plate_med, target_id !=  "BCL2L2", target_id !=  "B-actin",  sample_id %in% c("CTRL", "Querc"))) + 
  geom_col(aes(x = sample_id, y = rel_abund), fill="red") + geom_text(aes(x = sample_id, y = rel_abund, label = "***"), vjust = -0.9)+
  labs(y = "RNA abundance relative to B-actin",x = "Treatment", title = "LMNB1") +
  geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund- 0.000214, ymax= rel_abund+ 0.000214), width=.1) + facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

only_resveratrol_lmnb1 + 
  scale_x_discrete(breaks = c("CTRL", "Querc"),
                   labels = c("Control", "Resveratrol"))



### solo resveratrol para BCL2L2 -> USAR SD FILTRADA PARA CORREGIR

only_resveratrol_bcl2 <- ggplot(data = filter(plate_med, target_id !=  "LMNB1", target_id !=  "B-actin",  sample_id %in% c("CTRL", "Querc"))) + 
  geom_col(aes(x = sample_id, y = rel_abund), fill="brown") + 
  labs(y = "RNA abundance relative to B-actin",x = "Treatment",
    title = "BCL2L2")+
  geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund- 0.00015, ymax= rel_abund+ 0.00015), width=.1)+ facet_wrap(~target_id,ncol=4)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

 grafica_bcl2_rsv_complete <- only_resveratrol_bcl2+ 
  scale_x_discrete(breaks = c("CTRL", "Querc"),
                   labels = c("Control", "Resveratrol"))
 
grafica_bcl2_rsv_complete +  geom_text(aes(x = sample_id, y = rel_abund, label = "***"), vjust = -0.9)

###ANOVA y prueba de tukey
##hacer subsets

bcl2_ANOVA <- subset(plate_norm, target_id  %in% "BCL2L2")

analysis <- aov(data = bcl2_ANOVA, rel_abund ~ sample_id)

summary(analysis)
print(analysis)

TukeyHSD(analysis)

lmn_ANOVA <- subset(plate_norm, target_id  %in% "LMNB1")
analysis2 <- aov(data = lmn_ANOVA, rel_abund ~ sample_id)
summary(analysis2)
print(analysis2)
TukeyHSD(analysis2)
