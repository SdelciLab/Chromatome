# Chromatome Analysis
pacman::p_load(piggyback, renv, here, tidyverse, targets, DEP,pheatmap,PeCorA,sva,imp4p,  vroom,
               org.Hs.eg.db,clusterProfiler,ggridges,SubCellBarCode,eulerr,scales,data.table,biomaRt,
               visNetwork,matrixStats,magick,testthat, openxlsx, janitor,seqinr)


Proteomic_Ruler <- here::here("Datasets","Processed","CCLE_prot_Ruler.txt") %>% read.delim() %>% .[-1,] %>% 
    dplyr::select(matches("Copy|Uniprot_Acc|accuracy"))%>% 
    remove_rownames() %>% 
    column_to_rownames("Uniprot_Acc") %>% 
    #mutate(across(where(is.numeric), as.numeric)) %>% 
    set_names(.,str_remove_all(names(.),"Copy\\.number\\.")) %>% 
    mutate(across(contains("_"),~log2(as.numeric(.x))),
           across(where(is.numeric), ~if_else(is.infinite(.x), NaN,.x))) %>% 
    subset(!is.nan(U2OS_BONE)) %>% 
    # subset(.,rowSums(is.na(.))<(ncol(.)/3)) %>%
    subset(Absolute.quantification.accuracy != "low") %>%
    dplyr::select(-Absolute.quantification.accuracy) %>%
    janitor::clean_names()
Proteomic_Ruler_global <- Proteomic_Ruler %>% rowMeans() %>% enframe("proteingroup","WCE") |> 
    as.data.table()
HUMAN_9606 <- fread(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"), header = F) %>%
    set_names(c("Uniprot","Type","ID")) %>% as.data.table()

Sabatini <- openxlsx::read.xlsx(here::here("Datasets","Raw","metabolism_gene_list_Sabatini.xlsx")) %>% 
    pull(Gene.Symbol)
Sabatini_Uniprot <- HUMAN_9606 %>% subset((ID %in% Sabatini) & Type == "Gene_Name") %>%
    subset(!duplicated(Uniprot)) %>% pull(Uniprot, ID)#%>% unique()


contaminant_uniprots <- seqinr::read.fasta(here::here(
    "Datasets","Raw","contaminants.fasta"),seqtype ="AA") %>% names() 

codes = openxlsx::read.xlsx(here::here("Datasets","Raw","Submission LIST.xlsx"), startRow = 3) %>% 
    .[1:64,] %>% as.data.table()
Run_names <- fread(here::here("Datasets","Raw","report.tsv"), select = "Run") %>% 
    unique() %>% 
    tidyr::separate(col = "Run",sep = "-", into = c("Experiment","run","batch","batch-2","sample","replicate")) #the dropped pieces are -new reaquisitions and the pools
Tissue <- (codes$X4 == "Sample Code") %>% which() #1  5 12 24 36 40 46 50 55 61
Tissues <- codes$X3[Tissue]
Tissue <- data.table(sample = codes$X4 %>% .[. != "Sample Code"],
                     Tissue = c(rep(Tissues[1],(5-2)),
                                rep(Tissues[2],(12-6)),
                                rep(Tissues[3],(24-13)),
                                rep(Tissues[4],(36-25)),
                                rep(Tissues[5],(40-37)),
                                rep(Tissues[6],(46-41)),
                                rep(Tissues[7],(50-47)),
                                rep(Tissues[8],(55-51)),
                                rep(Tissues[9],(61-56)),
                                rep(Tissues[10],(65-62))))
Tissue <- merge(Tissue,Run_names, on = "sample" ) %>% 
    merge(setnames(codes[X4 != "Sample Code",.(X2,X3,X4)],'X4','sample'))
setnames(Tissue, c("X2","X3"),c("Type","Cell_line"))
Tissue[,Condition:= paste(Tissue,Type,Cell_line,replicate, sep = "_")]
list_files_to_analyse <- list(main_report = "report.tsv")
Samples <-Tissue[,.(run,Condition)] 
fwrite(Samples,here::here('Datasets','Raw','sample_annoation_MS.csv'))
Load_DIA_NN_Data <- function(report_DIA_tsv_file, Samples_df){
    tmp_data_table <- fread(here::here("Datasets","Raw",report_DIA_tsv_file)) 
    tmp_data_table =tmp_data_table[Lib.Q.Value<= 0.01 & Lib.PG.Q.Value <= 0.01]
    tmp_data_table = tmp_data_table |> as.data.frame() |> janitor::clean_names() %>% 
        dplyr::select(matches("pg_max_lfq|run|protein_group|^genes$")) %>% 
        distinct() %>% 
        # mutate(old_pg =protein_group) %>%  
        filter(!str_detect(protein_group,paste0(contaminant_uniprots,collapse = "|")))
    tmp_data_table = remove_rownames(tmp_data_table) %>% 
        mutate(run_2 = run %>% str_match("-(P[:graph:]{5})-.") %>% .[,2] ) %>%  
        mutate(run = case_when(
            is.na(run_2) & str_detect(run,"M1045-2")~ "Pool_2",
            is.na(run_2) & str_detect(run,"M1045-3")~ "Pool_3",
            is.na(run_2) & str_detect(run,"M1045-4")~ "Pool_4",
            is.na(run_2) & str_detect(run,"M1045-5")~ "Pool_5",
            is.na(run_2) & str_detect(run,"M1045-7")~ "Pool_7",
            is.na(run_2) & str_detect(run,"M1045-6")~ "Pool_6",
            is.na(run_2) & str_detect(run,"M1045")~ "Pool_1",
            TRUE ~run_2
            
        )) 
    tmp_data_table = left_join(tmp_data_table,Samples_df ) %>% 
        mutate(Condition = if_else(is.na(Condition),run,Condition))
    tmp_data_table = tmp_data_table |> as.data.table() |> 
        dcast(genes+protein_group~ Condition, value.var = 'pg_max_lfq') 
    # # tmp_data_table <- pivot_wider(tmp_data_table,-run_2, names_from = Condition ,values_from = pg_max_lfq) 
    #     
    #     # group_split(protein_group) %>%
    #     # map_dfr(.x = .,~.x %>% discard_uniprot_isoforms("protein_group","genes" )) %>%
    #     # group_by(protein_group) %>%
    #     # mutate(Is_duplicated = n() > 1,
    #     #        New_Uniprot = if_else(Is_duplicated == F,protein_group,old_pg)) %>% 
    tmp_data_table[,genes:=NULL]
    tmp_data_table=    column_to_rownames(tmp_data_table |> as.data.frame(),"protein_group") %>% 
        
        # dplyr::select(any_of(tmp_data_table$Condition)) %>% 
        # discard_single_isoforms %>% 
        as.matrix() #%>% log2()
    
    tmp_data_table
    
    
}
Methods_DIA <- purrr::map(.x = list_files_to_analyse[1],
                   ~Load_DIA_NN_Data(.x,Samples))
fwrite(tmp_data_table |> as.data.frame() |> 
           rownames_to_column('PG'),here::here('Datasets','Processed','filtered_DIANN_report.csv'))

global_means = Methods_DIA$main_report |> log2() |> matrixStats::rowMeans2( na.rm = T)
Chrom_global_mean = data.table(ProteinGroups = rownames(Methods_DIA$main_report),
                               Intensity =          global_means)
Chrom_global_mean[,proteingroup:= str_remove(ProteinGroups,';[:print:]*$')]
ratios = Proteomic_Ruler_global[Chrom_global_mean, on = 'proteingroup'][,Ratio:=log2(Intensity/WCE), by = proteingroup]
geneList = ratios[is.finite(Ratio)][order(Ratio, decreasing = T)] |> pull(Ratio, proteingroup)
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              keyType = "UNIPROT",
              verbose      = FALSE)
fwrite(ego3@result,here::here('Datasets','Processed','chromatin_enrichment_FigS1A.csv'))
enrichment = ridgeplot(ego3)+scale_fill_gradient(low = 'grey50', high = 'grey90')+
    ggtitle('Chromosome enriched terms validate chromatome protocol')
enrichment
ggsave(here::here('Output','figures','chromatin_enrichment.pdf'))


tmp_data_table <- fread(here::here("Datasets","Raw",list_files_to_analyse[1])) 
tmp_data_table = janitor::clean_names(tmp_data_table[Lib.Q.Value<= 0.01 & Lib.PG.Q.Value <= 0.01] )
tmp_data_table = tmp_data_table[,.(psm_per_gene = .N), by = .(run,protein_group,genes,pg_max_lfq)] 
tmp_data_table= tmp_data_table[!str_detect(protein_group,paste0(contaminant_uniprots,collapse = "|"))]
tmp_data_table[,run_2 := run %>% str_match("-(P[:graph:]{5})-.") %>% .[,2]]
tmp_data_table[,run := case_when(
    is.na(run_2) & str_detect(run,"M1045-2")~ "Pool_2",
    is.na(run_2) & str_detect(run,"M1045-3")~ "Pool_3",
    is.na(run_2) & str_detect(run,"M1045-4")~ "Pool_4",
    is.na(run_2) & str_detect(run,"M1045-5")~ "Pool_5",
    is.na(run_2) & str_detect(run,"M1045-7")~ "Pool_7",
    is.na(run_2) & str_detect(run,"M1045-6")~ "Pool_6",
    is.na(run_2) & str_detect(run,"M1045")~ "Pool_1",
    TRUE ~run_2
    
)]
tmp_data_table[,run_2:=NULL] 
Samples = Samples |> as.data.table()
Samples[,`:=`(Tissue = str_remove_all(Condition,'_[:print:]*$') |> as.factor(),
              Cancer = str_remove(Condition,'[:print:]*?_') |> 
                  str_remove_all('_[:print:]*$') |> as.factor(),
              Sample =str_remove(Condition,'[:print:]*?_[:print:]*?_') |> str_remove('_.$') |> as.factor() )]
Samples[,rep:= 1:.N, by =.(Tissue,Cancer,Sample) ]
Samples[,id := 1:.N]
tmp_data_table = left_join(tmp_data_table,Samples ) %>% 
    mutate(Condition = if_else(is.na(Condition),run,Condition) )
tmp_data_table[,binned_psm:=case_when(
    between(psm_per_gene,1,2)~'<3',
    between(psm_per_gene,3,6)~'<6',
    between(psm_per_gene,7,10)~'<10',
    psm_per_gene>10~'10+',
    TRUE~'other')]
tmp_data_table[,binned_psm := factor(binned_psm,levels = rev(c('<3','<6','<10','10+')))]
# Set a number of 'empty bar' to add at the end of each group
data = tmp_data_table[!is.na(Tissue)][,.(prot_per_sample = .N), by = .(rep, id,binned_psm, Condition, Tissue, Cancer ,Sample)]
data$binned_psm = as.factor(data$binned_psm)
# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
nObsType <- nlevels(as.factor(data$binned_psm))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$Tissue)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$Tissue <- rep(levels(data$Tissue), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(Tissue,Cancer, Sample,id)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data <- data %>% group_by(id, Sample) %>% summarize(tot=sum(prot_per_sample))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
    group_by(Tissue) %>% 
    summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]
library(viridis)
# Make the plot
fwrite(data,here::here('Datasets','Processed','Fig1E.csv'))
p <- ggplot(data) +      
    
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=prot_per_sample, fill=binned_psm), stat="identity", alpha=0.5) +
    scale_fill_viridis(discrete=TRUE) +
    
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey80", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 500, xend = start, yend = 500), colour = "grey80", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 2000, xend = start, yend = 2000), colour = "grey80", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 4000, xend = start, yend = 4000), colour = "grey80", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 6000, xend = start, yend = 6000), colour = "grey80", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    
    # Add text showing the value of each 100/75/50/25 lines
    ggplot2::annotate("text", x = rep(max(data$id),5), 
                      y = c(0, 500, 2000, 4000, 6000), 
                      label = c("0", "500", "2000", "4000", "6000") , color="grey40", size=3 , angle=0, fontface="bold", hjust=1) +
    
    # ylim(-100,5000) +
    theme_minimal() +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm") 
    ) +
    scale_y_continuous(limits = c(-2500,8000))+
    coord_polar(start= 0.05) +
    
    # Add labels on top of each bar
    geom_text(data=label_data, aes(x=id, y=tot+50, label=Sample, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
    
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -500, label=Tissue),
              # hjust=c(1,1,0,0), colour = "black", 
              alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
ggsave(here::here('Output','figures','pept_per_prot_sample.pdf'))
# all_plots = readRDS(here::here('Output','all_plots.RDS'))
# all_plots$Prots_per_plot = p
# saveRDS(all_plots, here::here('Output','all_plots.RDS'))

ggsave(here::here('Output','Prot_per_sample.pdf'),width = 10,height = 10)
# library(clusterProfiler)
# gene = Methods_DIA$main_report[complete.cases(Methods_DIA$main_report),] |> rownames() |> str_remove(';[:print:]*$') |> unique()
# ego <- enrichGO(gene          = gene,
#                 # universe      = names(geneList),
#                 OrgDb         = org.Hs.eg.db,
#                 ont           = "CC",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05,
#                 keyType = 'UNIPROT',
#                 readable      = TRUE)
# 
# dotplot(simplify(ego), showCategory=30) + ggtitle("dotplot for ORA")
# 
# 
# Methods_DIA$main_report <-Methods_DIA$main_report %>% as.data.frame() %>% janitor::clean_names() %>% rownames_to_column("ProteinGroup")
# fwrite(Methods_DIA$main_report,here::here("Datasets","Processed","Input_report.tsv"))
# Methods_DIA <- fread(here::here("Datasets","Processed","Input_report.tsv")) %>% column_to_rownames("ProteinGroup")
# all_uniprot <- Methods_DIA %>% rownames() %>% str_split(";") %>% unlist()