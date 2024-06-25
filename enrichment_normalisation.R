library(data.table);library(lspline);library(ggplot2);library(MASS);library(stringr);library(WGCNA)
# Normalising out the enrichment
f_BIC<- function( input_data, value_name ){
    missingness_tmp <- base::crossprod(!is.na(input_data))>8 # with microproteins being very rare it's common to only have 1 point in common produce fake correlations
    #using suggested bicor
    tmp <- bicor( input_data , use = "pairwise.complete.obs", nThreads = 6)                     # Robust correlation using all pairwise complete observations
    tmp[missingness_tmp == FALSE] <- NaN
    tmp <- as.data.table( reshape2::melt( as.matrix( tmp  )))                                   # Turn distance matrix into a pair-wise data.table
    tmp <- tmp[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), value ) ]   # Rename and change to character
    tmp <- tmp[ Protein_1 > Protein_2 ]                                                         # Remove redundant pairs by keeping only A > B, removing A == B and B < A pairs
    names(tmp)[3] <- value_name                                                                 # Assign new value name
    return(tmp)
}

HUMAN_9606 <- fread(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"), header = F)  |> 
    as.data.table()
setnames(HUMAN_9606,colnames(HUMAN_9606),c("Uniprot","Type","ID"))

complexes = fread(here::here('Datasets','Raw','humap2_complexes_20200809.txt'),sep = ',') |> 
    dplyr::select(-genenames) |> subset(Confidence %in% 1:5) |> tidyr::separate_rows(Uniprot_ACCs, sep = ' ')

enrichment_plot_rep_sum_piv <-  fread(here::here("Datasets","Processed","enrichment_plot_rep_sum_piv.csv"))
data_norm <- fread(here::here("Datasets","Processed","chrom_data_norm.tsv"))
data_norm_centered = data_norm[,-1] |> as.matrix()
rownames(data_norm_centered) = data_norm$Uniprot
data_norm_centered = data_norm_centered[(is.na(data_norm_centered) |> 
                                             matrixStats::rowSums2())<30,]
data_norm_centered |> boxplot()
proteins_on_chromatin = data_norm_centered |> rownames()

proteins_on_chromatin = data.table(Proteins = data_norm$Uniprot,
           on_chromatin = data_norm$Uniprot %in%proteins_on_chromatin )
fwrite(proteins_on_chromatin,here::here('Datasets','Processed', 'proteins_on_chromatin.csv'))
proteins_on_chromatin = fread(here::here('Datasets','Processed', 'proteins_on_chromatin.csv'))
proteins_on_chromatin[,Uniprot:= str_remove_all(Proteins,';[:print:]*')]
prots_to_network = proteins_on_chromatin[on_chromatin ==T,Uniprot] |> unique()
GS_TP = fread(here::here('Datasets','Raw','ProHD2_Gold_standard_TP.gz'))
GS_TP_chrom = GS_TP[OLN_1 %in% prots_to_network & OLN_2 %in% prots_to_network]
setnames(GS_TP_chrom,c('from','to'))
nodes = data.frame(id  = prots_to_network)
nodes = left_join(nodes, complexes, by = c('id' = 'Uniprot_ACCs')) |> as.data.table()
nodes = nodes[,prots_per_complex :=.N, by = HuMAP2_ID ][order(-prots_per_complex)][,head(.SD,1),by = id]
nodes[,prots_per_complex :=.N, by = HuMAP2_ID ]
fwrite(nodes,here::here('Datasets','Processed', 'humap_complexes_chrom.csv'))

nodes[,HuMAP2_ID :=fifelse(prots_per_complex >20,HuMAP2_ID ,NA_character_)]
nodes <-  nodes %>%
    mutate(color = case_when(
        HuMAP2_ID == 'HuMAP2_02207'~"#332288",
        HuMAP2_ID == 'HuMAP2_01036' ~"#D55E00",
        HuMAP2_ID == 'HuMAP2_01573' ~"#88CCEE",
        HuMAP2_ID == 'HuMAP2_02763'~"#DDCC77",
        HuMAP2_ID == 'HuMAP2_03081' ~"#CC6677",
        HuMAP2_ID == 'HuMAP2_06754' ~"#AA4499",
        HuMAP2_ID == 'HuMAP2_00104' ~"#882255",
        HuMAP2_ID == 'HuMAP2_04055'~"white",
        HuMAP2_ID == 'HuMAP2_02572' ~"#006CD1",
               TRUE~"#CDCDCD"),
           label = fifelse(HuMAP2_ID == 'HuMAP2_02207', id,NA_character_) ) |> 
    dplyr::select(id,label,color)
complexes_to_show = nodes |> subset(color != "#CDCDCD") |> pull(id) |> unique()
GS_TP_chrom = GS_TP_chrom[from %in% complexes_to_show & to %in% complexes_to_show ]
nodes = nodes |> subset(id %in% unique(c(GS_TP_chrom$from,GS_TP_chrom$to)))

library(ggraph);library(igraph)
net <- igraph::graph_from_data_frame(d=GS_TP_chrom, vertices=nodes, directed=F) 
net <- igraph::simplify(net, remove.multiple = F, remove.loops = T) 

ggraph(net, layout = "dh") +
    geom_edge_link(color="gray30",alpha = 0.2) +
    geom_node_point(shape = 21, 
                    # fill=V(net)$color.border, 
                    size=3.5,
                    # labellabel = V(net)$label,
                    fill = V(net)$color ,stroke = 1) +
    # scale_edge_width(range = c(0.1, 2))+ # control size
    # scale_edge_size(range = c(1, 20))+ # control size
    # geom_node_text(label = V(net)$label, size=10, color="gray30", repel=T) +
    theme_void()
ggsave(here::here("Output","figures","presence_network.pdf"), width = 15,height = 10)
visnet$nodes %>% ggplot(aes(x = id, colour = -log2_FC, y = log2_FC))+
    geom_point()+
    scale_colour_gradient2(mid =c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[2] , 
                           high = c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[3],low = c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[1] )
ggsave(here::here("Output","Figures","PPI_network_high_scale.pdf"), width = 15,height = 10)



melted_data = data_norm_centered |> as.data.frame() |> 
    tibble::rownames_to_column('Uniprot') |> as.data.table() |> 
    melt(id.vars = 'Uniprot', variable.name = 'Sample',value.name = 'Abundance')

chrom_intensities_plot = copy(melted_data)
chrom_intensities_plot[,breast:=  fifelse(str_detect(Sample,'breast'),'breast','other')]
order_of_samples = chrom_intensities_plot |> 
    group_by(Sample) |> 
    summarise( median_intensity= median(Abundance,na.rm =T)) |>
    arrange(median_intensity)
# inner_join(enrichment_plot_rep_sum_piv,order_of_samples) |> 
#     ggplot(aes(x = median_intensity, y=  N))+
#     geom_point()
pre_model_histones =     ggplot(chrom_intensities_plot ,aes(y = Abundance, 
                                                            x = factor(Sample, levels = order_of_samples |> pull(Sample)),
                                                            fill = breast))+
    geom_boxplot()+
    ylab('Histone Intensity')+
    xlab('Chromatin Enriched Samples')+
    # facet_wrap('Tissue', scales = 'free_y')+
    theme_bw()+theme( axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c('darkred','grey90'))+
    ggtitle('Histone Abundance in Chromatin Enriched Samples varies even within tissues',
            subtitle = 'Incorporating this into a bayesian model could remove this effect')

fwrite(chrom_intensities_plot, here::here('Datasets','datatable','FigS2D_chromatin_average.csv'))
melted_data = melted_data[stringr::str_detect(Sample,'pbmc|astrocyte_cell_pellet_1',negate = T)]
melted_data = enrichment_plot_rep_sum_piv[melted_data, on = 'Sample']
ggplot(data =  melted_data[stringr::str_detect(Uniprot,'P13473|P84243$')],
       aes(y =Abundance, x = N))+
    geom_point()+theme_bw()+
    ggtitle('Both proteins on chromatin but to different % of their pool')+
    facet_wrap('Uniprot', scales = 'free_y')
ggsave(here::here('Output','figures','proteins_on_chrom_facets.pdf'))
fwrite( melted_data[stringr::str_detect(Uniprot,'P13473|P84243$')],here::here('Datasets','datatable','FigS2E_LAMP2_H33A.csv'))
ggplot(data =  melted_data[stringr::str_detect(Uniprot,'P01116')],
       aes(y =Abundance, x = N))+
    geom_point()+theme_bw()+
    ggtitle('Protein groups sharing peptides belonging to P01116-KRAS on chrom')+
    facet_wrap('Uniprot', scales = 'free_y')
ggsave(here::here('Output','figures','KRAS_on_chrom.pdf'),width = 9,height = 9)

residuals  = data.table()
predictors = data.table()
proteins = melted_data$Uniprot |> unique()
for(i in proteins){
    print(i)
    data_tmp  =melted_data[Uniprot == i]
    model_tmp = rlm(Abundance ~ N,  data = data_tmp)
    residuals_tmp = data.table(Samples = data_tmp$Sample,
                               residual = model_tmp$residuals,
                               Uniprot = i)
    predictor_tmp = data.table(Uniprot = i,
                               predictors =model_tmp$coefficients,
                               names_pred = names(model_tmp$coefficients))
    predictors = rbind(predictors,predictor_tmp)
    residuals = rbind(residuals_tmp,residuals)
}
residuals[str_detect(Uniprot,'P01116')]
fwrite(predictors,'predictors_onchrom.csv')
fwrite(residuals,'residuals_onchrom_unfilt.csv')
predictors = fread('predictors_onchrom.csv')
residuals_test =residuals[abs(residual)<4]
# residuals_test = residuals
normality = data.table()
for(i in unique(residuals_test$Uniprot) ){
    test_tmp =  shapiro.test(residuals_test[Uniprot == i,residual])
    normality = rbind(normality,
                      data.table(
                          Uniprot = i,
                          pval = test_tmp$p.value))
}
normality$pval |> hist()
normality[,adj_pval := p.adjust(pval,'fdr') ]
irregular_proteins = normality[adj_pval<0.1,Uniprot] |> unique()
proteins_to_remove = intersect(residuals[abs(residual)>4,Uniprot],irregular_proteins)
# residuals= fread('residuals_replicates_multiple_linear_regression.csv')
residuals =residuals[abs(residual)<4][!(Uniprot %in% proteins_to_remove)]
residuals_joined = copy(residuals)
protein_cor_piv =  residuals |> dcast(Samples ~ Uniprot, value.var = 'residual') |> 
    tibble::column_to_rownames('Samples') 

to_check = c('P10768','P05091','P11766')


residuals[str_detect(Uniprot,
                     paste(to_check, collapse = '|'))] |> View()

melted_data[str_detect(Uniprot,
                     paste(to_check, collapse = '|'))] |> View()

protein_cor = f_BIC(protein_cor_piv,'chrom_cor')
fwrite(protein_cor,here::here('Datasets','Processed','Chromatin_prot_cor.gz'))
protein_cor = fread(here::here('Datasets','Processed','Chromatin_prot_cor.gz'))
replicate_cor = copy(protein_cor)
replicate_cor[,`:=`(Protein_1  = Protein_2    ,
                    Protein_2 = Protein_1)]
corr_to_plot = rbind(protein_cor,replicate_cor) |> dcast(Protein_1~Protein_2, value.var = 'chrom_cor')|> 
    tibble::column_to_rownames('Protein_1') 
corr_to_plot = corr_to_plot[(is.na(corr_to_plot) |> matrixStats::rowSums2())<100,
                            (is.na(corr_to_plot) |> matrixStats::colSums2())<100,]

Annotation = data.frame(Uniprot_ACCs = str_remove(colnames(corr_to_plot),';[:print:]*$'),
                        ProteinGroup= colnames(corr_to_plot))
Annotation = left_join(Annotation,complexes) |> as.data.table()
Annotation[,N_prot := .N, by = HuMAP2_ID ]
Annotation[N_prot<6,HuMAP2_ID:= NA_character_]
Annotation[is.na(HuMAP2_ID),HuMAP2_ID := 'Not_complex']
Annotation = Annotation[order(HuMAP2_ID,Confidence), head(.SD,1), by =ProteinGroup]
complex_members = Annotation[!(HuMAP2_ID == 'Not_complex'),ProteinGroup]
corr_to_plot = corr_to_plot[rownames(corr_to_plot)%in% complex_members,colnames(corr_to_plot) %in% complex_members]
Annotation = Annotation[ProteinGroup %in% colnames(corr_to_plot)]
Annotation = Annotation[match(Annotation$ProteinGroup,colnames(corr_to_plot))]
library(ComplexHeatmap)
column_ha = HeatmapAnnotation(foo1 = Annotation$HuMAP2_ID)
corr_to_plot = rbind(protein_cor,replicate_cor) |> dcast(Protein_1~Protein_2, value.var = 'chrom_cor')|> 
    tibble::column_to_rownames('Protein_1') 
corr_to_plot = corr_to_plot[(is.na(corr_to_plot) |> matrixStats::rowSums2())<100,
                            (is.na(corr_to_plot) |> matrixStats::colSums2())<100,]
pdf(here::here('Output','figures','corre_heatmap.pdf'))
ht <- ComplexHeatmap::Heatmap(as.matrix(corr_to_plot)[1:200,1:200], show_row_names = F, 
                              show_column_names = F )
draw(ht)
dev.off()
fwrite(as.data.frame(corr_to_plot) |> tibble::rownames_to_column('PG'),
       here::here('Datasets','datatable','Fig3E_heatmap.csv'))

ProHD = fread(here::here('Datasets','Raw','coregulation_scores.csv'))
ProHD_limit = ProHD$coregulation_score |> quantile(probs = seq(0,1,0.005))
ProHD[,interacting_prots := coregulation_score>=ProHD_limit[200] ]
protein_cor_prohd = ProHD[protein_cor, on = c('Protein_1','Protein_2'), nomatch = NULL]
protein_cor_prohd |> 
    ggplot(aes(x = chrom_cor, fill =  as.factor(interacting_prots)))+
    geom_density(alpha = 0.3)+theme_bw()+
    scale_fill_manual(values = c('grey90','black'))+
    ggtitle('Protein_correlations reveal functional relationships')
ggsave(here::here('Output','figures','ProHD1_corre_values.pdf'))
fwrite(protein_cor_prohd,here::here('Datasets','datatable','FigS3A_proHD.csv'))
ggplot(protein_cor_piv,aes(x = Q9UQE7, y = Q14683))+
    geom_point()+theme_bw()+
    ggtitle(glue::glue('bicor {round(corr_to_plot["Q9UQE7","Q14683"],digits = 2)} SMC3 correlates with SMC1A'),
            subtitle = 'Supported by ProHD')
ggsave(here::here('Output','figures','SMC3_SMC1A_corre_values.pdf'))
fwrite(protein_cor_piv, here::here('Datasets','datatable','FigS3B_correlation.csv') )
ggplot(protein_cor_piv,aes(x = Q8IWA0, y = Q15061))+
    geom_point()+theme_bw()+
    ggtitle(glue::glue('WDR75 {round(corr_to_plot["Q8IWA0","Q15061"],digits = 2)} correlates with WDR43'),
            subtit = 'Not supported by ProHD')
ggsave(here::here('Output','figures','WDR75_WDR43_corre_values.pdf'))
ggplot(protein_cor_piv,aes(x = Q8WUM0, y = P09012))+
    geom_point()+
    ggtitle('WDR75 correlates with WDR43',
            subtit = 'Not supported by ProHD')
# LDHB = residuals[Uniprot == 'P07195']
# # LDHB = LDHB[,]
# LDHB[,`:=`(Tissue = str_remove(Samples,'_[:print:]*$'),
#            Normal = str_detect(Samples,'normal'),
#            cell_line = str_remove(Samples,'_.$'),
#            hormone_ins = str_detect(Samples,'hormone') )]   
# LDHB = LDHB[,.(mean_res = mean(residual)), by = .(Tissue,Normal,cell_line,hormone_ins)]
# LDHB = LDHB[Tissue == 'breast']
#     ggplot(LDHB,aes(y = mean_res    ,
#                x= reorder(cell_line,mean_res ), colour = Normal, fill = hormone_ins, label = cell_line))+
#     geom_col()+
#         # facet_wrap('Tissue')+
#         scale_colour_manual(values = c('black','white'))+
#         ggrepel::geom_label_repel()
residuals_joined[,Samples:= stringr::str_remove(Samples,'_.$')]
residuals_joined = residuals_joined[,.(mean_res = mean(residual)), by = .(Samples,Uniprot)]
fwrite(residuals_joined,'per_cell_line_residuals.csv')
residuals_joined = fread('per_cell_line_residuals.csv')
residuals_piv = residuals_joined[str_detect(Samples,'pool',negate = T)]|> 
    dcast(Uniprot ~ Samples, value.var =  'mean_res' )
residuals_joined[,`:=`(tissue = str_extract(Samples,'^[:print:]*?_'),
                       cell_line = str_extract(Samples,'_[:print:]*?$'))][Uniprot == 'P49588'] |> 
    ggplot(aes(x= cell_line, y = mean_res))+
    geom_col()+
    facet_wrap('tissue', scales = 'free_x') +
    scale_x_discrete(guide = guide_axis(angle = 90)) 
residuals_piv[stringr::str_detect(Uniprot,'P01116')]
for_PCA = residuals_piv[,-1]
for_PCA = impute::impute.knn(as.matrix(for_PCA))
row_sd = for_PCA$data |> matrixStats::rowSds()
annotation_col = data.table(Samples = colnames(for_PCA$data) )
annotation_col[,`:=`(normal = str_detect(Samples,'normal') |> as.numeric(),
                     tissue = str_remove(Samples,'_[:print:]*$'))]
annotation_col = annotation_col |> tibble::column_to_rownames('Samples')
ann_colors = list(
    normal  = c( `0`  ="#EFFAFD", `1` = "#38A961"),
    tissue = c(blood = "#54C8AC", brain = "#E38BC3",
                    breast = "#5EBF4A",
                    cervix = '#91CB45' , colon = '#C8A2D0', 
                    liver = '#51C27E' , lung= '#50D2F5' ,
                    pancreas = '#C6BE3B',prostate = '#A0B6DD' ,
                    skin = '#FF9E4E')
)
cols = colorRampPalette(rev(c("darkorange",'white'  ,'darkblue')))(30)
pdf(here::here("Output",'figures','tissue_heatmap.pdf'))
for_PCA$data[row_sd>(0.60),] |> pheatmap::pheatmap(annotation_col = annotation_col, color = cols,
                                                   annotation_colors =ann_colors, show_colnames = F,show_rownames = F )
dev.off()
fwrite(for_PCA$data[row_sd>(0.60),] |> as.data.frame() |> 
           tibble::rownames_to_column('PG'),here::here('Datasets','datatable','Fig2E_all_prot_heatmap.csv'))
pathways =    fread(here::here('Datasets','Raw','hsa_pathways.txt'))
pathways = pathways[Path_type == 'metabolic']
enzymes =    fread(here::here('Datasets','Raw','KEGG_genes.csv'))
enzyme_pathway_map = pathways[enzymes, on= c('Path_id' = 'pathway')
][!is.na(Path_type)][,.(Path_description,ID)] |> unique()
enzymes =  enzymes[pathway %in% pathways$Path_id][,ID]
enzymes_on_chromatin = proteins_on_chromatin[,.(Uniprot,on_chromatin)] |> unique()
enzymes_on_chromatin = HUMAN_9606[Type == 'Gene_Name'][enzymes_on_chromatin, on = 'Uniprot'
][ID %in% enzymes]
Pathways = enzyme_pathway_map[enzymes_on_chromatin, on = 'ID']

imputted = for_PCA$data
rownames(imputted) = residuals_piv$Uniprot
annotation_row = data.frame(Proteins = rownames(imputted)) |>
    dplyr::mutate(Uniprot = str_remove(Proteins,';[:print:]*$')) |> 
    dplyr::inner_join(Pathways[,.(Uniprot,Path_description)]) |> 
    dplyr::mutate(Path_description = str_remove(Path_description,' - Homo[:print:]*$')) |> 
    dplyr::select(Proteins,Path_description) |> as.data.table() |> unique()
annotation_row[,N_prots := .N, by = Path_description]
annotation_row = annotation_row[order(-N_prots),head(.SD,1), by = .(Proteins)]
annotation_row[,Path_description_simple := fifelse(N_prots>10,Path_description,'Other')]
annotation_row = annotation_row[,.(Proteins,Path_description_simple)] |> 
    tibble::column_to_rownames('Proteins') |> 
    dplyr::arrange(Path_description_simple)

ann_colors = list(
    normal  = c( `0`  ="#EFFAFD", `1` = "#38A961"),
    tissue = c(blood = "#54C8AC", brain = "#E38BC3",
               breast = "#5EBF4A",
               cervix = '#91CB45' , colon = '#C8A2D0', 
               liver = '#51C27E' , lung= '#50D2F5' ,
               pancreas = '#C6BE3B',prostate = '#A0B6DD' ,
               skin = '#FF9E4E')
)
cols = colorRampPalette(rev(c("darkorange",'white'  ,'darkblue')))(30)


pdf(here::here("Output",'figures','tissue_enzymes_heatmap.pdf'),width = 10,height = 10)
to_plot = imputted[rownames(imputted) %in% rownames(annotation_row),]
to_plot = to_plot[match(rownames(annotation_row),rownames(to_plot)),]
pheatmap::pheatmap(to_plot,
                   annotation_col = annotation_col, color = cols,
                   annotation_colors = ann_colors, show_colnames = F,
                   annotation_row = annotation_row,
                   cluster_rows = F, show_rownames = F)
dev.off()
fwrite(to_plot |> as.data.frame() |> 
           tibble::rownames_to_column('PG'),here::here('Datasets','datatable','Fig2F_all_enzy_heatmap.csv'))

tissue_comparison = dplyr::inner_join(annotation_row |> tibble::rownames_to_column('Uniprot'),
           imputted |> as.data.frame() |> tibble::rownames_to_column('Uniprot'), by = 'Uniprot') |> 
    tidyr::pivot_longer(cols = -c(1:2),names_to = 'Samples',values_to = 'Norm_abundance') |> as.data.table()
tissue_comparison[,`:=`(normal = fifelse(str_detect(Samples,'normal'),'healthy',
                                       'cancer'),
                      tissue = str_remove(Samples,'_[:print:]*$'))]


tissue_comparison[str_detect(tissue,'skin|lung|colon|breast')] |> 
    ggplot(aes(x = tissue, y = Norm_abundance  ))+
    geom_boxplot()+geom_point(aes(colour  = normal ))+
    facet_wrap('Path_description_simple')
to_test_pathways= residuals_joined[str_detect(tissue,'skin|lung|colon|breast') & str_detect(Samples,'cancer')]
pathways_pvals = data.table()
for(pathway in unique(Pathways$Path_description)){
    print(pathway)
    Uniprots_tmp = Pathways[Path_description == pathway,Uniprot]
    pathway_chrom= to_test_pathways[Uniprot %in% Uniprots_tmp]
    if(length(unique(pathway_chrom$tissue))>1){
    res_aov = aov(mean_res  ~ tissue, data = pathway_chrom
            )
    s_res_aov = summary(res_aov)
    pathways_pvals = rbind(pathways_pvals,
                           data.table(pathway = pathway,
                                      N_prots = length(unique(pathway_chrom[,Uniprot])),
                                      pval = s_res_aov[[1]][["Pr(>F)"]][1]))
    }
}
pathways_pvals[,padj := p.adjust(pval,method = 'fdr')]
pathways_plots = list()
for(pathway in pathways_pvals[order(padj)][1:15,pathway] ){
    to_plot = merge(to_test_pathways,
      Pathways[Path_description == pathway], by = 'Uniprot') 
    pathways_plots[[pathway]] = ggplot(to_plot[][, tissue:= str_remove(tissue,'_')],aes(x =tissue,
                                                   y = mean_res,
                                                   fill = tissue))+
        # add half-violin from {ggdist} package
        # stat_halfeye(
        #     # adjust bandwidth
        #     adjust = 0.5,
        #     # move to the right
        #     justification = -0.2,
        #     # remove the slub interval
        #     .width = 0,
        #     point_colour = NA
        # ) +
        
        labs(y = 'Normalised Protein Abundance', x= 'Cancer Lineage')+
        geom_boxplot(
            width = 0.2,
            position = position_nudge(x = 0.2),
            # removing outliers
            outlier.color = NA,
            alpha = 0.5
        ) +
        stat_dots(
            # ploting on left side
            side = "left",
            # adjusting position
            justification = 1.1,
            # adjust grouping (binning) of observations
            binwidth = 0.06,dotsize  = 1
        )+ ggtitle(pathway)+ theme_bw()+scale_fill_manual(values = c(skin = '#F8A04D',
                                                                     breast = '#60BB47',
                                                                     colon = '#C79DC8',
                                                                     lung = '#3EC7F3'))
}
fwrite(pathways_pvals, here::here('Datasets','Processed','pathways_pvals.csv'))
fwrite(merge(to_test_pathways,
             Pathways, by = 'Uniprot'),here::here('Datasets','datatable', 'Fig2G_pathways_enzymes.csv'))
fread(here::here('Datasets','Processed','pathways_pvals.csv')) |> View()
ggpubr::ggarrange(plotlist = pathways_plots[c(3,13,12)], ncol = 3)
ggsave(here::here('Output','Figures','KEGG_pathways_anova.pdf'), width =20, height = 6)

oxphos_III = fread('OXPHOS_complex_III_proteins_complex_mean_expression_lung.csv')
residuals_piv[str_detect(Uniprot,paste(oxphos_III$Uniprot,collapse = '|'))]
to_test_pathways= residuals_joined[str_detect(tissue,'skin|lung|colon|breast') & str_detect(Samples,'cancer')]
to_test_pathways = merge(to_test_pathways[,cell_line_name:= str_remove(Samples,'^[:print:]*_cancer_')],
      samples, by = 'cell_line_name')
mutations_pvals = data.table()
for(pathway in unique(Pathways$Path_description)){
    for(tissue_type in unique(to_test_pathways$tissue)){
    print(pathway)
    Uniprots_tmp = Pathways[Path_description == pathway,Uniprot]
    pathway_chrom= to_test_pathways[Uniprot %in% Uniprots_tmp & tissue == tissue_type]
    if(length(unique(pathway_chrom$mutations))>1){
        res_aov = aov(mean_res  ~ mutations, data = pathway_chrom
        )
        s_res_aov = summary(res_aov)
        mutations_pvals = rbind(mutations_pvals,
                               data.table(pathway = pathway,
                                          tissue = tissue_type,
                                          N_prots = length(unique(pathway_chrom[,Uniprot])),
                                          pval = s_res_aov[[1]][["Pr(>F)"]][1]))
    }
    }
}
mutations_pvals$pval |> hist()
mutations_pvals[,padj := p.adjust(pval,method = 'fdr')]
mutation_plots = list()
for(pathway in mutations_pvals[padj<0.01][,pathway] ){
    to_plot = merge(to_test_pathways,
                    Pathways[Path_description == pathway], by = 'Uniprot') 
    mutation_plots[[pathway]] = ggplot(to_plot,aes(x = factor(mutations,levels = c('no_mut','KRAS','TP53','TP53,KRAS')), y = mean_res, label = ID))+
        geom_boxplot()+theme_bw()+
        geom_point(alpha = 0.01)+
        ggrepel::geom_text_repel(data = to_plot[,head(.SD,1), by = ID], max.overlaps = 20, alpha = 0.7)+
        ggtitle(pathway)+
        facet_wrap('tissue')
}
ggpubr::ggarrange(plotlist = mutation_plots)




write.csv(imputted,here::here('Datasets','Processed','Imputted_normalised_residuals.csv'))

annot_hyperlopit <- fread(here::here("Datasets","Raw", "annot_hyperlopit.tsv"))

pca_res <- prcomp(data_norm_centered |> na.omit(), scale=T)
var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)

to_plot = pca_res$x %>% 
    as.data.frame %>%
    tibble::rownames_to_column('Uniprot') |> 
    mutate(Uniprot = Uniprot |> stringr::str_remove_all(';[:print:]*$')) |> 
    left_join(annot_hyperlopit) |> 
    mutate(compartment = fifelse(`final.assignment` %in% c('CHROMATIN','ER'),`final.assignment`,'OTHER')) 
ggplot(to_plot, aes(x=PC1,y=PC2, label = Uniprot, colour = compartment )) +
    geom_point(size=2, data = to_plot |> subset(compartment =='OTHER')) +
    geom_point(size=2, data = to_plot |> subset(compartment !='OTHER')) +
    # ggrepel::geom_text_repel(alpha = 0.7)+
    theme_bw(base_size=32) +
    labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
         y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
    # theme(legend.position="none") +
    # facet_wrap('final.assignment')+
    scale_colour_manual(values = c('darkred','darkblue','grey80'))+
    ggtitle("before_correction PCA")+
    theme(plot.title = element_text(size = 20))
ggsave(here::here("Output",'figures','before_correction_PCAProteins.pdf'))

fwrite(to_plot, here::here('Datasets','datatable','FigS2G_before_norm_PCA.csv'))

sample_res=    copy(residuals)
sample_res[,`:=`(tissue_name = stringr::str_remove(Samples,'_[:print:]*$'),
                    cell_line_name = stringr::str_remove(Samples,'_.$'),
                    normal = stringr::str_detect(Samples,'normal'))]
cancer_normal = sample_res[,.(mean_residual = mean(residual)), by = .(normal,Uniprot,tissue_name )]

cancer_normal = cancer_normal |> dcast(tissue_name + Uniprot~normal)
cancer_normal[,diff_cancer:= `FALSE` - `TRUE`]
quantiles_diff  = cancer_normal$diff_cancer |> abs() |> quantile(probs = seq(0, 1, 0.01),na.rm = T) 
cancer_normal[,change:= abs(diff_cancer)>quantiles_diff[99]]
changing_proteins = cancer_normal[change == T,Uniprot]
changing_proteins_multiple = changing_proteins |> table() |> tibble::enframe() |> dplyr::arrange(-value)
diff_prots_to_plot = changing_proteins_multiple |> head() |> pull(name)
cancer_to_plot = cancer_normal[!(tissue_name %in% c('blood','pool'))]  
    ggplot(cancer_to_plot, aes(y = tissue_name, x = diff_cancer , label = str_remove(Uniprot,';[:print:]*$')), ) +
    stat_density_ridges(aes(fill = factor(stat(quantile))),
        geom = "density_ridges_gradient",
        calc_ecdf = TRUE,
        quantiles = c(0.025, 0.975)
    ) +
    scale_fill_manual(
        name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
        labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
    )+theme_bw()+
    ggrepel::geom_label_repel(data = cancer_to_plot[change ==T & Uniprot %in% diff_prots_to_plot], aes(colour = Uniprot))+
    ggtitle('Abundance difference between cancer and normal chromatin samples')
ggsave(here::here('Output','figures','differential_proteins_cancer.pdf'))
fwrite(cancer_to_plot, here::here('Datasets','datatable','Fig3A_cancer_normal.csv'))
to_plot_cancer_normal = sample_res[Uniprot == 'P56378' & !(tissue_name %in% c('blood','pool'))]  
    ggplot(to_plot_cancer_normal[normal !=T], aes(x=   tissue_name, y = residual))+
        theme_bw()+
    geom_violin(fill = 'grey60', alpha = 0.2, adjust = 0.5)+
    geom_point(data = to_plot_cancer_normal[normal ==T][,.(residual = mean(residual)), 
                                                        by = .(cell_line_name ,normal,tissue_name,Uniprot)], 
               colour = 'red', size =2)+
        # coord_flip()+
        ggtitle('ATP5MJ ATP synthase membrane subunit j')
ggsave(here::here('Output','figures','ATP5MJ_differential_cancer.pdf')) 
fwrite(to_plot_cancer_normal, here::here('Datasets','datatable','Fig3B_ATP5MJ.csv'))
    to_plot_cancer_normal = sample_res[Uniprot == 'Q5T8D3' & !(tissue_name %in% c('blood','pool'))]  
    ggplot(to_plot_cancer_normal[normal !=T], aes(x=   tissue_name, y = residual))+
        theme_bw()+
        geom_violin(fill = 'grey60', alpha = 0.2, adjust = 0.75)+
        geom_point(data = to_plot_cancer_normal[normal ==T][,.(residual = mean(residual)), 
                                                            by = .(cell_line_name ,normal,tissue_name,Uniprot)], 
                   colour = 'red', size =2)+
        # coord_flip()+
        ggtitle('Acyl-CoA-binding domain-containing protein 5')
    ggsave(here::here('Output','figures','ACBD5_differential_cancer.pdf'))
    fwrite(to_plot_cancer_normal, here::here('Datasets','datatable','Fig3B_ACBD5.csv'))
    
    # downloaded from Epifactor 04/12/2023
    Epigenes = fread(here::here('Datasets','Raw','EpiGenes_main.csv'))
    protein_size = fread(here::here('Datasets','Raw','uniprotkb_proteome_UP000005640_AND_revi_2023_10_25.tsv.gz'))
    setnames(protein_size,'Entry','Uniprot_single')
    setnames(protein_size,'Gene Names (primary)','Gene name')
    predictors[,Uniprot_single := stringr::str_remove_all(Uniprot,';[:print:]*$')]
    protein_size = protein_size[predictors[names_pred  == 'N' & Uniprot_single  %in% 
                                               str_remove(proteins_on_chromatin[on_chromatin == T,Proteins],';[:print:]*$')], on = 'Uniprot_single']
    protein_size = NLS_predict[protein_size, on = c('Uniprot'='Uniprot_single')]
    protein_size[,is_enzyme := `Gene name` %in% enzymes]
    protein_size[, enzyme := fcase(
        `Gene name` %in% enzymes & `Gene name` %in% Epigenes$HGNC_symbol, 'Epigenetic enzyme',
        `Gene name` %in% enzymes,'enzyme',
        !(`Gene name` %in% enzymes),'non-enzyme')]
    protein_size[!is.na(NLS_found)] |> 
        ggplot(aes(x = log2(Mass), y = predictors, fill = enzyme , colour = NLS_found, label = `Gene name`  ))+
        geom_point(data= protein_size[!is.na(NLS_found) & is_enzyme ==F & NLS_found ==F], alpha = 0.2, shape = 21, size = 3)+
        geom_point(data= protein_size[!is.na(NLS_found) & is_enzyme ==F & NLS_found ==T], alpha = 0.4, shape = 21, size = 3)+
        geom_point(data= protein_size[!is.na(NLS_found) & is_enzyme ==T & NLS_found ==T], alpha = 0.6, shape = 21, size = 3)+
        geom_point(data= protein_size[!is.na(NLS_found) & is_enzyme ==T & NLS_found ==F], alpha = 0.8, shape = 21, size = 3)+
        theme_bw()+
        # facet_wrap('NLS_found')+
        ggrepel::geom_label_repel(data = protein_size[is_enzyme ==T & ( NLS_found == T)], 
                                  max.overlaps = 5, colour = 'black', fill = 'white', alpha = 0.8)+
        geom_vline(xintercept = log2(80000), linetype="dotted", linewidth = 1)+
        coord_flip()+
        scale_fill_manual(values =c('orange','darkred','grey80'))+
        scale_colour_manual(values =c('grey80','black'))+
        labs(y= 'Relative % on chromatin')+
        ggtitle('NLS containing proteins are primarily nuclear, this suggests that moonlighting enzymes would not be able to perform their canonical faction if they had an NLS \nAdditionally There seems to be a positive relationship between Enzyme (NON-N) size and % on chromatin \n suggesting that they are not limited by diffusion for entering the nucleus \ndotted line 80Kda')
    ggsave(here::here('paper_NLS_enzyme_all.pdf'), width = 9, height  =9)
    fwrite(protein_size,here::here('Datasets','datatable','FigS2L_enzyme_size_on_chrom.csv'))
    
    # KRAS_abundance = data_norm[str_detect(Uniprot,'P01116')] |>     
    #     melt(id.vars = 'Uniprot', variable.name = 'Sample',value.name = 'Abundance') 
    # KRAS_abundance[,`:=`(cell_line= str_remove(Sample,'_.$') |> str_remove('^[:print:]*?_[:print:]*?_'),
    #                      tissue = str_remove(Sample,'_[:print:]*$'))][,cell_line:=str_sub(cell_line,1,10)]
    #                      
    # KRAS_abundance_norm = KRAS_abundance[str_detect(Uniprot,'P01116'),.(mean_abundance = mean(Abundance,na.rm = T)), by = .(cell_line,tissue,Uniprot)]
    # KRAS_abundance_norm[,mean_Uniprot:= mean(mean_abundance, na.rm =T),by = Uniprot]
    # KRAS_abundance_norm[,centered_abundance:= mean_abundance - mean_Uniprot]
    # sample_order = KRAS_abundance_norm[tissue == 'lung'][Uniprot == 'P01111;P01116'][order(centered_abundance),cell_line]
    #     ggplot(KRAS_abundance_norm[tissue == 'lung'],aes(x = factor(cell_line, levels = sample_order), y = centered_abundance))+
    #     geom_col()+
    #     facet_wrap('Uniprot',scales = 'free_x')+scale_x_discrete(guide = guide_axis(angle = 90)) +
    #     ggtitle('KRAS protein groups on lung chromatin')+ 
    #         labs(x = 'samples',y = 'KRAS protein groups abundace')
    # ggsave(here::here('KRAS_PGs_on_lung_chrom_per_sample.pdf'),height = 8, width = 8)
    # 
    # ggplot(KRAS_abundance_norm,aes(x = reorder(cell_line, tissue), fill = tissue,y = centered_abundance))+
    #     geom_col()+
    #     facet_grid(Uniprot~tissue, scales = 'free_x')+scale_x_discrete(guide = guide_axis(angle = 90)) +
    #     ggtitle('KRAS protein groups on all tissues chromatin')+ 
    #     labs(x = 'samples',y = 'KRAS protein groups abundace')
    # ggsave(here::here('KRAS_PGs_on_chrom_all_tissues.pdf'),height = 8, width = 10)
    # 
    # AARS1_abundance = data_norm[str_detect(Uniprot,'P49588')] |>     
    #     melt(id.vars = 'Uniprot', variable.name = 'Sample',value.name = 'Abundance') 
    # AARS1_abundance[,`:=`(cell_line= str_remove(Sample,'_.$') |> str_remove('^[:print:]*?_[:print:]*?_'),
    #                      tissue = str_remove(Sample,'_[:print:]*$'))][,cell_line:=str_sub(cell_line,1,10)]
    # 
    # AARS1_abundance = AARS1_abundance[str_detect(Uniprot,'P49588'),.(mean_abundance = mean(Abundance,na.rm = T)), by = .(cell_line,tissue,Uniprot)]
    # AARS1_abundance[,mean_Uniprot:= mean(mean_abundance, na.rm =T),by = Uniprot]
    # AARS1_abundance[,centered_abundance:= mean_abundance - mean_Uniprot]
    # sample_order = AARS1_abundance[tissue == 'breast'][Uniprot == 'P49588'][order(centered_abundance),cell_line]
    # ggplot(AARS1_abundance[tissue == 'breast'],aes(x = factor(cell_line, levels = sample_order), y = centered_abundance))+
    #     geom_col()+
    #     facet_wrap('Uniprot',scales = 'free_x')+scale_x_discrete(guide = guide_axis(angle = 90)) +
    #     ggtitle('AARS1 protein groups on lung chromatin')+ 
    #     labs(x = 'samples',y = 'AARS1 protein groups abundace')
    # ggsave(here::here('KRAS_PGs_on_lung_chrom_per_sample.pdf'),height = 8, width = 8)
    # 
    # ggplot(AARS1_abundance,aes(x = reorder(cell_line, tissue), fill = tissue,y = centered_abundance))+
    #     geom_col()+
    #     facet_grid(Uniprot~tissue, scales = 'free_x')+scale_x_discrete(guide = guide_axis(angle = 90)) +
    #     ggtitle('KRAS protein groups on all tissues chromatin')+ 
    #     labs(x = 'samples',y = 'KRAS protein groups abundace')
    # ggsave(here::here('KRAS_PGs_on_chrom_all_tissues.pdf'),height = 8, width = 10)
    # 
    # 
    # mutations_data = read.csv(here::here('Datasets','Raw','gene_set_library_crisp.gmt'), header = F) |> as.data.table()
    # mutations_data[,cell_line:= str_extract(V1,'^[:print:]*')]
    # mutations_data[1,]
    # TP53_muts = mutations_data[str_detect(V1,'TP53'),cell_line]
    # KRAS_muts = mutations_data[str_detect(V1,'KRAS'),cell_line]
    # 
    # samples[,mutations:= fcase(
    #     !(Cell_Name %in% mutations_data$cell_line) ,'no_data',
    #     # Cell_Name %in% TP53_muts & Cell_Name %in% KRAS_muts,'TP53,KRAS',
    #     # Cell_Name %in% TP53_muts ,'TP53',
    #     Cell_Name %in% KRAS_muts,'KRAS',
    #     !(Cell_Name %in% KRAS_muts),'no_KRAS_mut')]
    # 
    # ggplot(merge(KRAS_abundance_norm, samples, all.x = T, by.x = 'cell_line',by.y = 'cell_line_name'),
    #        aes(x = reorder(cell_line,centered_abundance), fill = mutations,y = centered_abundance))+
    #     geom_col()+
    #     facet_grid(Uniprot~tissue, scales = 'free_x')+scale_x_discrete(guide = guide_axis(angle = 90)) +
    #     ggtitle('KRAS protein groups on all tissues chromatin')+ 
    #     labs(x = 'samples',y = 'KRAS protein groups abundace')
    # ggsave(here::here('KRAS_PGs_on_chrom_all_tissues_mut.pdf'),height = 8, width = 12)
    # heatmap_fatemeh = data_norm[str_detect(Uniprot, paste(to_check, collapse = '$|^'))] |> 
    #     tibble::column_to_rownames('Uniprot') |> as.matrix() |> is.na() 
    # heatmap_fatemeh_2 = heatmap_fatemeh|>  apply(2, as.numeric) 
    # heatmap_fatemeh_2 = 1-heatmap_fatemeh_2
    # rownames(heatmap_fatemeh_2) = rownames(heatmap_fatemeh)
    # # colnames(heatmap_fatemeh_2) = colnames(heatmap_fatemeh_2) |> str_remove_all('_epithelial_cell_pellet|_ephytelial_cells_|human_')
    # 
    # annotation_col_ = data.table(Samples = colnames(heatmap_fatemeh_2) )
    # annotation_col_[,`:=`(normal = fifelse(str_detect(Samples,'normal'),'healthy',
    #                                       'cancer'),
    #                      tissue = str_remove(Samples,'_[:print:]*$'))]
    # annotation_col_ = annotation_col_ |> tibble::column_to_rownames('Samples')
    # pheatmap::pheatmap(heatmap_fatemeh_2, 
    #                    annotation_col =annotation_col_, show_colnames  = F )
    
    # samples = fread(here::here('Datasets','Raw','cell_line_reference.csv'))
    # setnames(samples, 'Original_name','cell_line_name')
    # samples[,cell_line_name:=stringr::str_remove(cell_line_name,'^[:print:]*?_cancer_')]
    # samples = samples[Cell_Name != '']
    # mutations_data = fread(here::here('Datasets','Raw','gene_attribute_edges.txt.gz'), skip = 1)
    # KRAS_muts = mutations_data[GeneSym %in% c('KRAS'),CellLine]
    # TP53_muts = mutations_data[GeneSym %in% c('TP53'),CellLine]
    # samples[,mutations:= fcase(
    #     !(Cell_Name %in% mutations_data$CellLine) ,'no_data',
    #     Cell_Name %in% TP53_muts & Cell_Name %in% KRAS_muts,'TP53,KRAS',
    #     Cell_Name %in% TP53_muts ,'TP53',
    #      Cell_Name %in% KRAS_muts,'KRAS',
    #     !(Cell_Name %in% TP53_muts | Cell_Name %in% KRAS_muts),'no_mut')]
    