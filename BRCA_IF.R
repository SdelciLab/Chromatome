#title: "BRCA1 IF"
#author: "Alisa Schmidt"
#date: "22/09/2025"

#libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)

#MDA-MB-231 IMPDH2 KO cells reconstituted with IMPDH2 (KOWT) or IMPDH2-NLS (KONLS)
#load data
data <- read.delim(
  "/Users/aschmidt/Desktop/Chromatome_paper/BRCA1_IF/Objects_Population - full nuclei.txt", 
  skip = 9)
data_cond <- data %>%
  mutate(
    Cell = case_when(
      Row == 4 ~ "KOWT",
      Row == 5 ~ "KONLS"))

table(data_cond$Cell)  
  

ggplot(data_cond, aes(x = Cell, y = full.nuclei...Intensity.NucleusBRCA.Mean, fill = Cell)) +
  geom_boxplot() + 
  labs(x = "", y = "full.nuclei...Intensity.NucleusBRCA.Mean", fill = "") +
  theme_classic() +
  scale_fill_brewer(palette = "Accent") 

median_values <- data_cond %>%
  group_by(Cell) %>%
  summarise(Median_Intensity = median(full.nuclei...Intensity.NucleusBRCA.Mean, na.rm = TRUE), .groups = "drop")
median_values

ggplot(data_cond, aes(x = Cell, y = full.nuclei...Number.of.Spots.BRCA..per.Cell, fill = Cell)) +
  geom_boxplot() + 
  labs(x = "", y = "full.nuclei...Number.of.Spots.BRCA..per.Cell", fill = "") +
  theme_classic() +
  scale_fill_brewer(palette = "Accent") 


ggplot(data_cond, aes(x = Cell, y = full.nuclei...areaum2BRCA.Area..µm....Mean.per.Cell, fill = Cell)) +
  geom_boxplot() + 
  labs(x = "", y = "full.nuclei...areaum2BRCA.Area..µm....Mean.per.Cell", fill = "") +
  theme_classic() +
  scale_fill_brewer(palette = "Accent") 
median_values <- data_cond %>%
  group_by(Cell) %>%
  summarise(Median_Intensity = median(full.nuclei...areaum2BRCA.Area..µm....Mean.per.Cell, na.rm = TRUE), .groups = "drop")
median_values



##add cell cycle information ------------------------------------------------
data_cond_Cc <- data_cond

data_cond_Cc$DAPIcorrected <- data_cond_Cc$full.nuclei...Intensity.Nucleus.DAPI.Mean - 
  data_cond_Cc$full.nuclei...Intensity.Ring.Region.DAPI.Mean

data_cond_Cc$ratio_corrected <- data_cond_Cc$DAPIcorrected * 
  data_cond_Cc$full.nuclei....Area..µm..


#check cell cycle profile 
ggplot(data_cond_Cc, aes(ratio_corrected, color = factor(Column))) + geom_density() + facet_wrap("Cell")
data_cond_Cc_KOWT <- data_cond_Cc[data_cond_Cc$Cell == "KOWT",]
a1 <- ggplot(data_cond_Cc_KOWT, aes(ratio_corrected, color = factor(Column))) + geom_density()

#Align replicates KOWT G1 peaks
density_data2 <- ggplot_build(a1)$data[[1]]
group_levels2 <- levels(factor(data_cond_Cc_KOWT$Column))
density_data2$Column <- group_levels2[density_data2$group]

#check for peak value per condtion
density2 <- density_data2 %>%
  group_by(Column) %>%
  mutate(MaxY = max(y)) %>% 
  dplyr::filter(y == MaxY) %>% 
  summarise(x = first(x), MaxY = first(MaxY), .groups = 'drop')

# First x-coordinate
First_peak <- density2$x[density2$Column == 4]

# Calculate shifts needed for each replicate
density2$shift2 = First_peak - density2$x

# Join data with shifts
density2$Column <- as.numeric(density2$Column)
 
KOWT_normalised <- data_cond_Cc_KOWT %>%
  left_join(density2, by = "Column") %>%
  mutate(adjusted_value_G1 = ratio_corrected + shift2) 
KOWT_normalised <- KOWT_normalised[, -c(47:49)]

ggplot(KOWT_normalised, aes(adjusted_value_G1, color = factor(Column))) +geom_density() 

#merge with KONLS
data_cond_Cc_KONLS <- data_cond_Cc[data_cond_Cc$Cell == "KONLS",]
data_cond_Cc_KONLS$adjusted_value_G1 <- data_cond_Cc_KONLS$ratio_corrected

CC_merged <- rbind(data_cond_Cc_KONLS,KOWT_normalised)

a2 <- ggplot(CC_merged, aes(adjusted_value_G1, color = Cell)) + geom_density() 

#align G1 -------
#Align replicates KOWT G1 peaks
density_data <- ggplot_build(a2)$data[[1]]
group_levels <- levels(factor(CC_merged$Cell))
density_data$Cell <- group_levels[density_data$group]

#check for peak value per condtion
density <- density_data %>%
  group_by(Cell) %>%
  mutate(MaxY = max(y)) %>% 
  dplyr::filter(y == MaxY) %>% 
  summarise(x = first(x), MaxY = first(MaxY), .groups = 'drop')

# KOWT x-coordinate
KOWT_peak <- density$x[density$Cell == "KOWT"]

# Calculate shifts needed for each replicate
density$shift = KOWT_peak - density$x

all_normalised <- CC_merged %>%
  left_join(density, by = "Cell") %>%
  mutate(adjusted_value_final = adjusted_value_G1 + shift) 


ggplot(all_normalised, aes(adjusted_value_final, color = Cell)) + geom_density() + xlim(0,2000000)
#remove outlier high values
all_normalised <- all_normalised[all_normalised$adjusted_value_final < 2000000, ]

#add gates ----------------
#x-coordinates
G1 <- KOWT_peak 
S <- G1 *1.61
G2 <- G1*2.22

gate1 <- G1  -189673
gate2 <- G1  + 189673
gate3 <- S + 189673
gate4 <- G2 +189673

#add G1 peak to plot
ggplot(all_normalised, aes(adjusted_value_final, color = Cell)) + geom_density() + 
   geom_vline(xintercept = G1 , col = "red", linetype = "dashed") + 
  geom_vline(xintercept = G2 , col = "orange", linetype = "dashed") + 
  geom_vline(xintercept = S , col = "yellow", linetype = "dashed") + 
   geom_vline(xintercept = gate1, col = "blue") +
   geom_vline(xintercept = gate2, col = "blue") +
  geom_vline(xintercept = gate3, col = "blue") +
  geom_vline(xintercept = gate4, col = "blue") +
 theme_bw() + labs(title= "Cell cycle profile", x = "integrated DAPI intensity", y = "", fill = "") 


data_cell_cycle_filtered <- all_normalised %>%
  mutate(Cycle = case_when(
    adjusted_value_final >= gate1 &  ratio_corrected <= gate2 ~ "G1",
    adjusted_value_final > gate2 &  adjusted_value_final <= gate3 ~ "S",
    adjusted_value_final > gate3 &  adjusted_value_final <= gate4 ~ "G2",
    TRUE ~ "Rest" 
  ))%>%
  dplyr::filter(Cycle != "Rest")

order <- c("G1","S", "G2")
data_cell_cycle_filtered$Cycle <- factor(data_cell_cycle_filtered$Cycle, levels=order)


#check integrated mean intensity
data_cell_cycle_filtered$integratedBRCA <- data_cell_cycle_filtered$full.nuclei....Area..µm.. *
  data_cell_cycle_filtered$full.nuclei...Intensity.NucleusBRCA.Mean

data_cell_cycle_filtered_outlier_removed <- data_cell_cycle_filtered %>%
  mutate(mean_value = mean(integratedBRCA, na.rm = TRUE),
         sd_value = sd(integratedBRCA, na.rm = TRUE)) %>%
  dplyr::filter(integratedBRCA >= (mean_value - 5 * sd_value) & 
                  integratedBRCA <= (mean_value + 5 * sd_value)) %>%
  dplyr::select(-mean_value, -sd_value)

ggplot(data_cell_cycle_filtered_outlier_removed, aes(x = Cycle, y = integratedBRCA, fill = Cell)) +
  geom_boxplot() + 
  labs(x = "", y = "integratedBRCA", fill = "") +
  theme_classic() +
  scale_fill_brewer(palette = "Accent") 


##Final plots -------------
##BRCA foci ----
#remove outlier, sample data 
#1)Filter outliers

df <- data_cell_cycle_filtered %>%
  mutate(
    mean_value = mean(full.nuclei...Number.of.Spots.BRCA..per.Cell, na.rm = TRUE),
    sd_value   = sd(full.nuclei...Number.of.Spots.BRCA..per.Cell, na.rm = TRUE)
  ) %>%
  dplyr::filter(
    full.nuclei...Number.of.Spots.BRCA..per.Cell >= (mean_value - 5 * sd_value),
    full.nuclei...Number.of.Spots.BRCA..per.Cell <= (mean_value + 5 * sd_value)
  ) %>%
  dplyr::select(-mean_value, -sd_value)

# 2) Find the smallest Cell×Cycle group size
min_n <- df %>% count(Cell, Cycle) %>% summarise(min_n = min(n)) %>% pull(min_n)

# 3) Sample that many rows from each group
set.seed(123) 
data_cell_cycle_filtered_outlier_removed2 <- df %>%
  group_by(Cell, Cycle) %>%
  slice_sample(n = min_n) %>%
  ungroup()


ggplot(data_cell_cycle_filtered_outlier_removed2, aes(x = Cell, y = full.nuclei...Number.of.Spots.BRCA..per.Cell, fill = Cell)) +
  geom_boxplot() +  facet_wrap(~Cycle) +
  stat_compare_means(
    method = "wilcox.test",   
    comparisons = list(c("KOWT", "KONLS"))
   #label = "p.signif"          
  ) +
  labs(x = "", y = "full.nuclei...Number.of.Spots.BRCA..per.Cell", fill = "") +
  theme_classic() +
  scale_fill_brewer(palette = "Accent") 
ggsave("Number.of.Spots.BRCA.pdf",width = 14, height = 10, units = "cm")


checkCCpics <-  data_cell_cycle_filtered_outlier_removed2[, -c(18:43,45:51)]

##Mean intensity----
#remove outlier, sample data 
#1)Filter outliers
#check integrated mean intensity
data_cell_cycle_filtered$integratedBRCA <- data_cell_cycle_filtered$full.nuclei....Area..µm.. *
  data_cell_cycle_filtered$full.nuclei...Intensity.NucleusBRCA.Mean

data_cell_cycle_filtered_outlier_removed <- data_cell_cycle_filtered %>%
  mutate(mean_value = mean(integratedBRCA, na.rm = TRUE),
         sd_value = sd(integratedBRCA, na.rm = TRUE)) %>%
  dplyr::filter(integratedBRCA >= (mean_value - 5 * sd_value) & 
                  integratedBRCA <= (mean_value + 5 * sd_value)) %>%
  dplyr::select(-mean_value, -sd_value)


# 2) Find the smallest Cell×Cycle group size
min_n <- data_cell_cycle_filtered_outlier_removed %>% count(Cell, Cycle) %>% summarise(min_n = min(n)) %>% pull(min_n)

# 3) Sample that many rows from each group
set.seed(123) 
data_cell_cycle_filtered_outlier_removed3 <- data_cell_cycle_filtered_outlier_removed %>%
  group_by(Cell, Cycle) %>%
  slice_sample(n = min_n) %>%
  ungroup()


ggplot(data_cell_cycle_filtered_outlier_removed3, aes(x = Cell, y = integratedBRCA, fill = Cell)) +
  geom_boxplot() +  facet_wrap(~Cycle) +
  stat_compare_means(
    method = "wilcox.test",   
    comparisons = list(c("KOWT", "KONLS"))
    #label = "p.signif"          
  ) + 
  labs(x = "", y = "integratedBRCA", fill = "") +
  theme_classic() +
  scale_fill_brewer(palette = "Accent") 
ggsave("integratedBRCA.pdf",width = 14, height = 10, units = "cm")
table(data_cell_cycle_filtered_outlier_removed3$Cell, data_cell_cycle_filtered_outlier_removed3$Cycle)


#yH2AX
df <- data_cell_cycle_filtered %>%
  mutate(
    mean_value = mean(full.nuclei...Number.of.Spots.H2AX..H2AX, na.rm = TRUE),
    sd_value   = sd(full.nuclei...Number.of.Spots.H2AX..H2AX, na.rm = TRUE)
  ) %>%
  dplyr::filter(
    full.nuclei...Number.of.Spots.H2AX..H2AX >= (mean_value - 5 * sd_value),
    full.nuclei...Number.of.Spots.H2AX..H2AX <= (mean_value + 5 * sd_value)
  ) %>%
  dplyr::select(-mean_value, -sd_value)

# 2) Find the smallest Cell×Cycle group size
min_n <- df %>% count(Cell, Cycle) %>% summarise(min_n = min(n)) %>% pull(min_n)

# 3) Sample that many rows from each group
set.seed(123) 
data_cell_cycle_filtered_outlier_removed2 <- df %>%
  group_by(Cell, Cycle) %>%
  slice_sample(n = min_n) %>%
  ungroup()


ggplot(data_cell_cycle_filtered_outlier_removed2, aes(x = Cell, y = full.nuclei...Number.of.Spots.H2AX..H2AX, fill = Cell)) +
  geom_boxplot() +  facet_wrap(~Cycle) +
  stat_compare_means(
    method = "wilcox.test",   
    comparisons = list(c("KOWT", "KONLS"))
    #label = "p.signif"          
  ) +
  labs(x = "", y = "full.nuclei...Number.of.Spots.H2AX..H2AX", fill = "") +
  theme_classic() +
  scale_fill_brewer(palette = "Accent") 

ggsave("Number.of.Spots.H2AX.pdf",width = 14, height = 10, units = "cm")
table(data_cell_cycle_filtered_outlier_removed2$Cell, data_cell_cycle_filtered_outlier_removed2$Cycle)

#Plot ratio of cell cycle in barplot
table(data_cell_cycle_filtered$Cell, data_cell_cycle_filtered$Cycle)

total_cond <- data_cell_cycle_filtered %>% count(Cell, name = "N_cond")
counts <- data_cell_cycle_filtered %>% count(Cell, Cycle, name = "n")
plot_df2 <- counts %>%
  left_join(total_cond, by = "Cell") %>%
  mutate(ratio_within_condition = n / N_cond)

ggplot(plot_df2, aes(x = Cell, y = ratio_within_condition, fill = Cycle)) +
  geom_col(position = "stack") +
  #scale_y_continuous(labels = percent) +
  labs(x = "", y = "Fraction of cells per condition", fill = "Cycle") +
  theme_classic() + scale_fill_brewer(palette = "Accent")

ggplot(plot_df2, aes(x = Cell, y = ratio_within_condition, fill = Cycle)) +
  geom_col(position = "dodge") +
  #scale_y_continuous(labels = percent) +
  labs(x = "", y = "Fraction of cells per condition", fill = "Cycle") +
  theme_bw() + scale_fill_brewer(palette = "Accent")