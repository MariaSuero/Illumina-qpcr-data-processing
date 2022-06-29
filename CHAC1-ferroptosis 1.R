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

sample_id_levels <- c("CTRL", "30", "20", "10","5", "NT")
prep_type_levels <- "[Erastin]"

colkey3 <- tibble(
  well_col = 1:6,
  sample_id = sample_id_levels,
  prep_type = prep_type_levels
)
print(colkey3)

### vale todo good he puesto los samples 
###pone para retener la info hay que meter esto lol 

create_blank_plate(well_row = LETTERS[1:4], well_col = 1:6)

###ahora uno las cosas a ver si sirve 

plan_ferroptosis_r1 <- label_plate_rowcol(plate = create_blank_plate(well_row = LETTERS[1:4], well_col = 1:8),
                                          rowkey = rowkey,
                                          colkey = colkey3)


print(plan_ferroptosis_r1)

##alucinante creo que todo bien 

display_plate_qpcr(plan_ferroptosis_r1)


columns <- c("well", "cq")
file_path_cq <- read.table(file = "Ferropt-r1-mss copia 4.txt", header = FALSE,  sep = "", dec = ",", col.names = columns)
print(file_path_cq )
print(plan_ferroptosis_r1)

plates = merge(x = file_path_cq, y = plan_ferroptosis_r1, by = "well",
          )
plates
##parece que va pq lo lee 


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

ggplot(data = plate_norm) +
  geom_point(aes(x = sample_id, y = delta_cq, shape = prep_type, colour = prep_type),
             position = position_jitter(width = 0.2, height = 0)
  ) +
  labs(y = "delta Cq") +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 

###desv standar 

standard_dev <-  plate_norm %>% group_by(sample_id,target_id) %>%
  summarise(sd = sd(rel_abund))

view(standard_dev)
             
###  hacer un heatmap CON COLORES NEUTROS

Plate_deltacq <- display_plate_value(plate_norm, value = "cq") +   # uses ggplot syntax
  scale_fill_gradient(high = "steelblue") 
Plate_deltacq 

## no consigo meterle una leyenda :(

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

deltacq <- ggplot(data = filter(plate_med, target_id != "B-actin") ) + 
  geom_point(aes(x = sample_id, y = delta_cq), color="red" ) + 
  labs(
    y = "delta Cq",
    x = "[Erastin]"
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

deltacq



## abundancia RNA 

abundancia_rna <- ggplot(data = filter(plate_med, target_id != "B-actin")) + 
  geom_line(aes(x = sample_id, y = rel_abund , colour = prep_type_levels ), color="steelblue", group=prep_type_levels) +
  geom_point(aes(x = sample_id, y = rel_abund), color="darkblue") + 
  labs(
    y = "RNA abundance relative to B-actin",
    x = "[Erastin/µM]"
    
  ) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  


###meter sd 

 ###hago subset de la sd 

chac1_sd <- subset(standard_dev, target_id  %in% "CHAC1")
###meto la grafica

ggplot(data = filter(plate_med, target_id != "B-actin")) + 
  geom_line(aes(x = factor(sample_id, level = c('CTRL', '5', '10', "20", "30", "NT")), y = rel_abund , colour = prep_type_levels ), color="darkblue", group=prep_type_levels) +
  geom_point(aes(x = sample_id, y = rel_abund), color="steelblue") + 
  labs( y = "RNA abundance relative to B-actin", x = "[Erastin/µM]") + scale_x_discrete(breaks = c('CTRL', '5', '10', "20", "30", "NT"),
                       labels = c("Control", "5", "10", "20", "30", " ")) + geom_text(aes(x = sample_id, y = rel_abund, label = "***"), vjust = -2)+
  geom_errorbar(aes(x = sample_id, y = rel_abund, ymin= rel_abund-chac1_sd$sd, ymax= rel_abund+chac1_sd$sd), width=.1,
                position = position_dodge(0.5), size = .5) +
  facet_wrap(~target_id,ncol=4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

### regresion 
stat <- c("sample", "rna")
lm <- read.csv(file = "chac1-datos-lm.txt" ,header = FALSE,  sep = "", dec = ".", col.names = stat)
sapply(lm, class)
lm
analysis6 <- lm( sample ~ rna, data= lm)

summary(analysis6)





