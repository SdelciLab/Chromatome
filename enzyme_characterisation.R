# enzyme characterisation
library(ggplot2);library(data.table);library(stringr)
imputted = read.csv(here::here('Datasets','Processed','Imputted_normalised_residuals.csv'), row.names = 1)
predictors = fread('predictors_onchrom.csv')
HUMAN_9606 <- fread(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"), header = F)  |> 
    as.data.table()
setnames(HUMAN_9606,colnames(HUMAN_9606),c("Uniprot","Type","ID"))
enrichment_plot_rep_sum_piv = fread(here::here('Datasets','Processed','enrichment_plot_rep_sum_piv.csv') )
#enzymes 
pathways =    fread(here::here('Datasets','Raw','hsa_pathways.txt'))
# pathways = pathways[Path_type == 'metabolic']
enzymes =    fread(here::here('Datasets','Raw','KEGG_genes.csv'))
enzyme_pathway_map = pathways[enzymes, on= c('Path_id' = 'pathway')
                              ][!is.na(Path_type)][,.(Path_description,ID)] |> unique()
enzymes =  enzymes[pathway %in% pathways$Path_id][,ID]

MitoGenes = readxl::read_xls(here::here('Datasets','Raw', 'Human.MitoCarta3.0.xls'), 
                             sheet = 2) |> as.data.table()
MitoGenes_oxphos = MitoGenes[str_detect(MitoCarta3.0_MitoPathways,'OXPHOS') &
                                 str_detect(MitoCarta3.0_MitoPathways,'subunit'), .(UniProt,Symbol,
                                                                               Description,
                                                                               MitoCarta3.0_MitoPathways)]
setnames(MitoGenes_oxphos,'UniProt','Uniprot')
proteins_on_chromatin = fread(here::here('Datasets','Raw', 'proteins_on_chromatin.csv'))
proteins_on_chromatin[,Uniprot:= str_remove_all(Proteins,';[:print:]*')]
enzymes_on_chromatin = proteins_on_chromatin[,.(Uniprot,on_chromatin)] |> unique()
OXPHOS_on_chromatin = enzymes_on_chromatin[MitoGenes_oxphos, on = 'Uniprot']
OXPHOS_on_chromatin=  OXPHOS_on_chromatin[!(Uniprot == 'P56385' & on_chromatin == F)]
OXPHOS_on_chromatin[,Complex_ID := str_extract(`MitoCarta3.0_MitoPathways`,'Complex[:print:]*?>') |> str_remove('( )>')]

fwrite(OXPHOS_on_chromatin,'OXPHOS_proteins_complex_annotation.csv')
fwrite(predictors[names_pred == 'N'][OXPHOS_on_chromatin[Complex_ID == 'Complex III'], on = 'Uniprot' ],here::here('Datasets','Processed','Fig_2H_OXPHOS_complex_III_proteins_complex_per_on_chrom.csv'))
residuals = fread('residuals_onchrom_unfilt.csv')
residuals_III = residuals[str_detect(Samples,'cancer')][,tissue:= str_remove(Samples,'_[:print:]*$')][,.(mean_exp = mean(residual,na.rm = T)), by = .(Uniprot,tissue)]
residuals_III = residuals_III[OXPHOS_on_chromatin[Complex_ID == 'Complex III'], on = 'Uniprot' ]
residuals_III[tissue %in% c('lung','breast','skin','colon')] |> ggplot(aes(x = tissue, y= Symbol, colour = mean_exp))+
    geom_point(size = 5)+ theme_bw()+
    scale_colour_gradient2(low = '#4B417E',high = '#F6E2A0', mid = '#E98891')
ggsave(here::here('Output','OXPHOS_III_expression.pdf'))
fwrite(residuals_III,here::here('Datasets','Processed','Fig_2H_OXPHOS_complex_III_proteins_complex_mean_expression_lung.csv'))
OXPHOS_on_chromatin_per_complex = OXPHOS_on_chromatin[,.(Perc_coverage = mean(on_chromatin,na.rm = T),
                                                         N_members = .N), by =Complex_ID ]
fwrite(OXPHOS_on_chromatin_per_complex[!is.na(Complex_ID)],here::here('Datasets','Processed','Fig_2C_OXPHOS_complexes_perc.csv'))


enzymes_on_chromatin = HUMAN_9606[Type == 'Gene_Name'][enzymes_on_chromatin, on = 'Uniprot'
                                                       ][ID %in% enzymes]
data_norm <- fread(here::here("Datasets","Processed","chrom_data_norm.tsv"))
data_norm_centered = data_norm[,-1] |> as.matrix()
rownames(data_norm_centered) = data_norm$Uniprot
data_norm_centered = data_norm_centered[(is.na(data_norm_centered) |> 
                                             matrixStats::rowSums2())<30,]


# Pathways = enzyme_pathway_map[enzymes_on_chromatin, on = 'ID']
presence_on_chrom = data_norm_centered[str_detect(rownames(data_norm_centered), 
                                                  paste(proteins_on_chromatin[on_chromatin == T,Uniprot],collapse = '|')),
                                       str_detect(colnames(data_norm_centered),'cancer')]
rownames_to_add = rownames(presence_on_chrom)
presence_on_chrom = presence_on_chrom |> is.na() |> apply(2,function (x) 1-as.numeric(x))
rownames(presence_on_chrom) = rownames_to_add
presence_on_chrom_plot = presence_on_chrom[(presence_on_chrom |> matrixStats::rowSums2())<78,] 
annotation_row = data.frame(Proteins = rownames(presence_on_chrom_plot)) |> 
    dplyr::mutate(Uniprot = str_remove(Proteins,';[:print:]*$')) |> 
    dplyr::left_join(Pathways[,.(Uniprot,Path_description)]) |> 
    dplyr::mutate(Path_description = str_remove(Path_description,' - Homo[:print:]*$')) |> 
    dplyr::select(Proteins,Path_description) |> unique() |> 
    dplyr::group_by(Proteins) |> dplyr::arrange(Proteins,Path_description) |> 
    dplyr::top_n(1) |> as.data.table()
annotation_row[,N_prots := .N, by = Path_description]
annotation_row[,Path_description_simple := fifelse(N_prots>5,Path_description,'Other')]
annotation_row = annotation_row[,.(Proteins,Path_description_simple)] |> 
    tibble::column_to_rownames('Proteins') |> 
    dplyr::arrange(Path_description_simple)
annotation_col = data.table(Samples = colnames(presence_on_chrom_plot) )
annotation_col[,`:=`(tissue = str_remove(Samples,'_[:print:]*$'))]
annotation_col[,Samples:=str_remove(Samples,'^[:print:]*_cancer_')]
annotation_col = annotation_col |> tibble::column_to_rownames('Samples')
colnames(presence_on_chrom_plot) <- str_remove(colnames(presence_on_chrom_plot),'^[:print:]*_cancer_')

pdf(here::here('Output','figures','presence_absence_all_prot.pdf'))
presence_on_chrom_plot[rownames(annotation_row),] |> 
    pheatmap::pheatmap(show_rownames = F,annotation_row =annotation_row, 
                       annotation_col = annotation_col,
                       cluster_rows = F,color = c( "white","grey40") )

dev.off()
fwrite(presence_on_chrom_plot[rownames(annotation_row),] |> as.data.frame() |> 
           tibble::rownames_to_column('PG'),here::here('Datasets','datatable','FigS2C_heatmap_pathways.csv'))
##enzymes 
pathways =    fread(here::here('Datasets','Raw','hsa_pathways.txt'))
pathways = pathways[Path_type == 'metabolic']
enzymes =    fread(here::here('Datasets','Raw','KEGG_genes.csv'))
enzyme_pathway_map = pathways[enzymes, on= c('Path_id' = 'pathway')
][!is.na(Path_type)][,.(Path_description,ID)] |> unique()
enzymes =  enzymes[pathway %in% pathways$Path_id][,ID]
proteins_on_chromatin = fread(here::here('Datasets','Raw', 'proteins_on_chromatin.csv'))
proteins_on_chromatin[,Uniprot:= str_remove_all(Proteins,';[:print:]*')]
enzymes_on_chromatin = proteins_on_chromatin[,.(Uniprot,on_chromatin)] |> unique()
enzymes_on_chromatin = HUMAN_9606[Type == 'Gene_Name'][enzymes_on_chromatin, on = 'Uniprot'
][ID %in% enzymes]
Pathways = enzyme_pathway_map[enzymes_on_chromatin, on = 'ID']

presence_on_chrom = data_norm_centered[str_detect(rownames(data_norm_centered), 
                                                  paste(enzymes_on_chromatin[on_chromatin == T,Uniprot],collapse = '|')),
                                       str_detect(colnames(data_norm_centered),'cancer')]
rownames_to_add = rownames(presence_on_chrom)
presence_on_chrom = presence_on_chrom |> is.na() |> apply(2,function (x) 1-as.numeric(x))
rownames(presence_on_chrom) = rownames_to_add
presence_on_chrom_plot = presence_on_chrom[(presence_on_chrom |> matrixStats::rowSums2())<78,] 
annotation_row = data.frame(Proteins = rownames(presence_on_chrom_plot)) |> 
    dplyr::mutate(Uniprot = str_remove(Proteins,';[:print:]*$')) |> 
    dplyr::left_join(Pathways[,.(Uniprot,Path_description)]) |> 
    dplyr::mutate(Path_description = str_remove(Path_description,' - Homo[:print:]*$')) |> 
    dplyr::select(Proteins,Path_description) |> unique() |> 
    dplyr::group_by(Proteins) |> dplyr::arrange(Proteins,Path_description) |> 
    dplyr::top_n(1) |> as.data.table()
annotation_row[,N_prots := .N, by = Path_description]
annotation_row[,Path_description_simple := fifelse(N_prots>4,Path_description,'Other')]
annotation_row = annotation_row[Path_description_simple != 'Other']
annotation_row = annotation_row[,.(Proteins,Path_description_simple)] |> 
    tibble::column_to_rownames('Proteins') |> 
    dplyr::arrange(Path_description_simple)

colour_to_plot = c('grey95','black','#82B7FF','#00D65C','#00DAE0','#D3BA00')
names(colour_to_plot)=     annotation_row$Path_description_simple |> unique()

colour_annot = list(
    Path_description_simple = colour_to_plot
)
annotation_col = data.table(Samples = colnames(presence_on_chrom_plot) )
annotation_col[,`:=`(tissue = str_remove(Samples,'_[:print:]*$'))]
annotation_col[,Samples:=str_remove(Samples,'^[:print:]*_cancer_')]
annotation_col = annotation_col |> tibble::column_to_rownames('Samples')
colnames(presence_on_chrom_plot) <- str_remove(colnames(presence_on_chrom_plot),'^[:print:]*_cancer_')

pdf(here::here('Output','figures','presence_absence_enzymes.pdf'))
presence_on_chrom_plot[rownames(annotation_row),] |> 
    pheatmap::pheatmap(show_rownames = F,annotation_row =annotation_row,  show_colnames = F,
                       annotation_col = annotation_col,
                       cluster_rows = F,color = c( "white","grey40")  )
dev.off()
fwrite(presence_on_chrom_plot |> as.data.frame() |> tibble::rownames_to_column('PG'),
       here::here('Datasets','Processed','Fig2D_enzyme_heatmap.csv'))

Perc_pathway =Pathways[,.(perc_chrom = mean(on_chromatin),
                                                                      N_prots = .N), by = Path_description]
Perc_pathway[,pathway:= str_remove(Path_description,' - Homo[:print:]*$')]
fwrite(Perc_pathway, here::here('Datasets','Processed','Fig2B_perc_on_chrom_pathway.csv'))
Perc_pathway[unique(c(6,9,14,sample(1:nrow(Perc_pathway),20)))] |> ggplot(aes(y = perc_chrom, x= reorder(pathway,perc_chrom), size = N_prots))+
    geom_col()+
    coord_flip()+theme_bw()+
    scale_size_continuous(range = c(1,5))+
    ggtitle('Perc of KEGG pathways found in 80% of samples, \ncompared to all detected enzymes ')
ggsave(here::here("Output",'figures','perc_pathway_on_chrom.pdf'),width = 10,height = 8)
KEGGs = HUMAN_9606[Type =='GeneID'][Pathways, on = 'Uniprot' ]
BP_enrichment = clusterProfiler::enrichGO(
    gene = unique(proteins_on_chromatin[on_chromatin==T,Uniprot]),
    ont = 'BP',
    OrgDb = 'org.Hs.eg.db',
    keyType = 'UNIPROT',
    universe = unique(proteins_on_chromatin[,Uniprot]))
BP_enrichment_simple = BP_enrichment |> clusterProfiler::simplify()
fwrite(BP_enrichment_simple@result,here::here('Datasets','datatable','FigS2A_BP_on_chrom.csv'))
BP_enrichment_simple@result <- BP_enrichment_simple@result |> 
    subset(str_detect(Description,'DNA|process'))
clusterProfiler::dotplot(BP_enrichment_simple, showCategory=50)+
    ggtitle('KEGG Pathways enriched on chromatin')
ggsave(here::here("Output",'figures','BP_on_chrom_stats.pdf'))
KEGG_enrichment = clusterProfiler::enrichKEGG(KEGGs[on_chromatin==T,ID],
                                              # keyType = 'UNIPROT',
                                              universe = KEGGs$ID)
fwrite(KEGG_enrichment@result,here::here('Datasets','datatable','FigS2B_KEGG_on_chrom.csv'))

clusterProfiler::dotplot(KEGG_enrichment, showCategory=30)+
    ggtitle('KEGG Pathways enriched on chromatin')
ggsave(here::here("Output",'figures','pathway_on_chrom_stats.pdf'))
annot_hyperlopit <- fread(here::here("Datasets","Raw", "annot_hyperlopit.tsv"))

pca_res <- prcomp(imputted  |>  t(), scale=F)
var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)

PCA_tissue = pca_res$x %>% 
    as.data.frame %>%
    tibble::rownames_to_column("Sample") |> 
    dplyr::mutate(Condition = str_remove(Sample,"_.$")) %>% 
    dplyr::mutate(Tissue = str_match(Condition, "^([:graph:]*?)_") %>% .[,2],
           Tissue = dplyr::if_else(is.na(Tissue),"Pool",Tissue),
           Type = str_match(Sample,"^[:graph:]*?_([:graph:]*?)_") %>% .[,2],
           Samples = str_remove_all(Sample,"^[:graph:]*?_[:graph:]*?_")) %>% 
    dplyr::left_join(enrichment_plot_rep_sum_piv) 
    ggplot(PCA_tissue,aes(x=PC1,y=PC2, label = Sample, colour = Tissue )) +
    geom_point()+
    # ggrepel::geom_text_repel(alpha = 0.7)+
    theme_bw(base_size=32) +
    labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
         y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
    # theme(legend.position="none") +
    ggtitle("original PCA")+
    # facet_wrap("Tissue")+
    theme(plot.title = element_text(size = 20))+
    stat_ellipse(na.rm = T, aes(fill = Tissue), level = 0.75)
ggsave(here::here("Output",'figures','tissue_PCA.pdf'))

fwrite(PCA_tissue,here::here('Datasets','datatable','FigS2I_tissue_PCA.csv'))

pca_res <- prcomp(imputted, scale=F)
var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)

# pca_res$x %>%
#     as.data.frame %>%
#     mutate(Uniprot = residuals_piv$Uniprot |> stringr::str_remove_all(';[:print:]*$')) |>
#     left_join(annot_hyperlopit) |>
#     ggplot(aes(x=PC1,y=PC2, label = Uniprot, colour = `final.assignment` )) +
#     geom_point(size=2) +
#     # ggrepel::geom_text_repel(alpha = 0.7)+
#     theme_bw(base_size=32) +
#     labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#          y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
#     # theme(legend.position="none") +
#     ggtitle("after_correction PCA")+
#     facet_wrap('final.assignment')+
#     theme(plot.title = element_text(size = 20))

to_plot = pca_res$x %>% 
    as.data.frame %>%
    tibble::rownames_to_column('Uniprot') |> 
    dplyr::mutate(Uniprot = Uniprot |> stringr::str_remove_all(';[:print:]*$')) |> 
    dplyr::left_join(annot_hyperlopit) |> 
    dplyr::mutate(compartment = fifelse(`final.assignment` %in% c('CHROMATIN','ER'),`final.assignment`,'OTHER')) 
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
    ggtitle("after_correction PCA")+
    theme(plot.title = element_text(size = 20))
ggsave(here::here("Output",'figures','after_correction_PCAProteins.pdf'))
fwrite(to_plot,here::here("Datasets",'datatable','FigS2H_after_norm_PCA.csv'))



predictors[,names_simplified := stringr::str_remove_all(Uniprot,';[:print:]*$')]
annot_hyperlopit[predictors, on = c('Uniprot' = 'names_simplified')][!is.na(`final.assignment`)] |> 
    ggplot(aes(x = `final.assignment`, y = predictors))+
    geom_boxplot()+
    facet_wrap('names_pred')
data_to_plot = annot_hyperlopit[predictors, on = c('Uniprot' = 'names_simplified')
] |> dplyr::left_join(HUMAN_9606[Type == 'Gene_Name'], by = 'Uniprot') |>
    dplyr::mutate(enzyme = ID %in% enzymes) 
ggplot(data_to_plot[names_pred =='N'],aes(x = `final.assignment`, y = predictors, label = ID))+
    geom_boxplot()+
    ggtitle('Perc on chromatin')+theme_bw()+ coord_flip()
ggsave(here::here('Output','figures','predictors_LOPIT.pdf'))    
fwrite(data_to_plot[names_pred =='N'],here::here('Datasets','datatable','FigS2J_slopes_lopit.csv'))

# finding proteins specifically upregulated/ downregulated in normal

# changing_cancer = residuals_piv[Uniprot %in% changing_proteins]
# changing_cancer = changing_cancer[,-1]
# changing_cancer = impute::impute.knn(as.matrix(changing_cancer))
# # row_sd = for_PCA$data |> matrixStats::rowSds()
# changing_cancer$data |> pheatmap::pheatmap()
# to_plot_pheatmap = for_PCA$data
# rownames(to_plot_pheatmap) = residuals_piv$Uniprot
# to_plot_pheatmap[changing_proteins,] |> pheatmap::pheatmap(main = 'Proteins changing between cancer and respective normal')
# changing_enzymes = changing_proteins[changing_proteins %in% metabolic_proteins$Uniprot]
# to_plot_pheatmap[changing_enzymes,] |> 
#     pheatmap::pheatmap(main = 'Enzymes changing between cancer and respective normal', cluster_cols = F)
NLS_predict  = fread(here::here("Datasets","Raw", "NLS_predict.csv"))
NLS_predict_compartments = annot_hyperlopit[NLS_predict, on = 'Uniprot'][!is.na(`final.assignment`)
][`final.assignment`!= 'unknown']  
    ggplot(final.assignment,aes(x= `final.assignment`, fill = NLS_found))+
    geom_bar()+
    coord_flip()+
    theme_bw()+
    scale_fill_manual(values = c('grey50','darkred'))+
    ggtitle('Predict_NLS')
ggsave(here::here('Output','figures','predictNLS_hyperlopit.pdf'))
fwrite(NLS_predict_compartments,here::here('Datasets','datatable','FigS1C_NLS_predict_compartments.csv'))

data_to_plot = NLS_predict[data_to_plot, on = 'Uniprot']
ggplot(data_to_plot[names_pred == 'N'][!is.na(NLS_found)],aes(x = predictors, fill = NLS_found))+
    geom_density(alpha = 0.3)+theme_bw()+
    # geom_jitter(alpha = 0.05)+
    # facet_wrap('enzyme',nrow = 2)+ 
    ggtitle('NLS sequence based signals',
            subtitle = 'agrees with Nuclear proteins have higher percentage on chromatin')
ggsave(here::here('Output','figures','predictNLS_predictor.pdf'))
fwrite(data_to_plot[names_pred == 'N'][!is.na(NLS_found)],here::here('Datasets','datatable','FigS2K_NLS_slopes.csv'))

# chromatome_proteins = residuals_piv$Uniprot |> stringr::str_remove(';[:print:]*') |> unique()
# chromatome_proteins = HUMAN_9606[Type == 'Gene_Name' & Uniprot %in% chromatome_proteins,.(Uniprot,ID)]
# opecell_interactome = fread(here::here("Datasets","Raw", "opencell-protein-interactions.csv"))
# opencell_interactome = data.table()
# for(i in chromatome_proteins$ID){
#     opencell_interactome = 
#         rbind(opencell_interactome,
#               data.table( Gene = i,
#                           interactors = setdiff(opecell_interactome[target_gene_name  == i | interactor_gene_name  == i,
#                                                                     .(target_gene_name,interactor_gene_name)] |> unlist() |> unique(),i)))
# }
# HUMAN_9606 = fread(here::here('Datasets','Raw','HUMAN_9606_idHUMAN_9606.dat'))
# setnames(HUMAN_9606,c('Uniprot','Type','ID'))
# ENSG= HUMAN_9606[Type == 'Ensembl'][,Type:=NULL]
# # HPA localisation 
HPA_location = fread(here::here('Datasets','Raw','subcellular_location.tsv'))
sort_compartments <- function(x){
    strsplit(x,';') |> unlist() |> sort() |> paste(collapse = ';')
}
HPA_location[,localication:= paste(`Main location`,`Additional location`, sep = ';') |> sort_compartments(), by = Gene]
HPA_location = HPA_location[,.(`Gene name`,localication)]
# setnames(HPA_location,c('Gene name','localication'))


# opencell_interactome = HPA_location[opencell_interactome,on = 'interactors',allow.cartesian=TRUE] |> na.omit()
# opencell_interactome[,nuclear := str_detect(localication  ,'Nucl')]
# Perc_nuclear = opencell_interactome[,.(perc_nucl = mean(nuclear),
#                                        N_int = .N), by = Gene]
lopit_genes = HUMAN_9606[Type == 'Gene_Name'][annot_hyperlopit, on = 'Uniprot'][,.(Uniprot,ID,`final.assignment`)] |> unique()
# setnames(Perc_nuclear,'Gene','ID')
# Perc_nuclear = lopit_genes[Perc_nuclear, on = 'ID'] 
# 
# ggplot(Perc_nuclear,aes(x = `final.assignment`, y = perc_nucl)) +
#     geom_jitter()
# Perc_nuclear[predictors, on = c('Uniprot' = 'names_simplified')][names_pred == 'N'][N_int>2] |> 
#     ggplot(aes(x =  predictors,y = perc_nucl, colour = log10(N_int)))+
#     geom_point()+
#     facet_wrap('final.assignment')

protein_size = fread(here::here('Datasets','Raw','uniprotkb_proteome_UP000005640_AND_revi_2023_10_25.tsv.gz'))
setnames(protein_size,'Entry','Uniprot_single')
setnames(protein_size,'Gene Names (primary)','Gene name')
predictors[,Uniprot_single := stringr::str_remove_all(Uniprot,';[:print:]*$')]
protein_size = protein_size[predictors[names_pred  == 'N'], on = 'Uniprot_single']
protein_size = lopit_genes[protein_size, on = c('Uniprot' = 'Uniprot_single')] 
protein_size = HPA_location[protein_size, on = 'Gene name']
protein_size[,N_localisation := .N, by  =localication]
protein_size[,Mass:=log10(Mass)]
ggplot(protein_size[N_localisation>10],aes(x= Mass,y = predictors))+
    geom_point()+
    facet_wrap('localication', scales = 'free')

protein_size[,is_enzyme := `Gene name` %in% enzymes ]
ggplot(protein_size,aes(x= is_enzyme ,y = predictors,label = `Gene name`))+
    # geom_boxplot()+
    geom_jitter()+
    ggrepel::geom_label_repel(data = protein_size[predictors>0.2 & is_enzyme ==T], max.overlaps = 100)+
    ggtitle('Enzymes with slope higher than 0.2')
# coord_flip()

# testing if missing proteins are less chromatin
Methods_DIA_input <- fread(here::here("Datasets","Processed","Input_report.tsv"))
detected = data.frame(ProteinGroup = Methods_DIA_input$ProteinGroup,
                         N_detected = ncol(Methods_DIA_input[,-1])-(Methods_DIA_input[,-1] |> 
                is.na()  |> apply(2,as.numeric) |> matrixStats::rowSums2()))
ggplot(detected, aes(x = N_detected))+
    geom_histogram()+
    ggtitle('Number of samples detected per protein')+theme_bw()

ggsave(here::here('Output','figures','samples_per_prot.pdf')) 
fwrite(detected,here::here('Datasets','datatable','FigS1B_detected.csv'))

detected =detected |>  dplyr::mutate(Uniprot = str_remove_all(ProteinGroup,';[:print:]*$')) |> 
    dplyr::mutate(has_NLS = Uniprot %in% NLS_predict[NLS_found ==T,Uniprot] |> as.numeric())

ggplot(detected, aes(y = N_detected, fill = as.factor(has_NLS)))+
    geom_density(alpha = 0.5)+theme_bw()+ coord_flip()+
        ggtitle('Number of samples detected per protein')+
    scale_fill_manual(values = c('grey90','black'))
ggsave(here::here('Output','figures','predictNLS_presence.pdf'))
fwrite(detected,here::here('Datasets','datatable','FigS1D_detected.csv'))

# 
# annotation_row = data.frame(ProteinGroup = Methods_DIA_input$ProteinGroup) |> 
#     dplyr::mutate(Uniprot = str_remove_all(ProteinGroup,';[:print:]*$')) |> 
#     dplyr::mutate(has_NLS = Uniprot %in% NLS_predict[NLS_found ==T,Uniprot] |> as.numeric()) |> tibble::column_to_rownames('ProteinGroup')
# annotation_row = annotation_row |> dplyr::select(has_NLS)
# # NLS_predict
# 
# annotation_row = dplyr::left_join(annotation_row,HUMAN_9606[Type == 'Gene_Name']) |> dplyr::left_join(HPA_location, by = c('ID' = 'Gene name'))
# annotation_row = annotation_row |> as.data.table()
# annotation_row[,N_location:= .N, by  = localication]
# annotation_row[,N_uniprot:= .N, by  = ProteinGroup]
# annotation_row[,nuclear := str_detect(localication,'Nucl') |> as.numeric()]
# annotation_row[,localisation := fifelse(N_location<50 | is.na(localication),'Other',localication )]
# annotation_row = annotation_row[N_uniprot==1,.(ProteinGroup,localisation,nuclear)] |> tibble::column_to_rownames('ProteinGroup')
# 
# data
# rownames(Methods_DIA_input_heatmap) = Methods_DIA_input$ProteinGroup
# 
# pheatmap_plot = Methods_DIA_input_heatmap[,enrichment_plot_rep_sum_piv[order(N),Sample]] |> 
#     pheatmap::pheatmap(cluster_cols = F, annotation_row =annotation_row, show_rownames = F,color = c("grey40", "white"), )
# pdf(here::here('Output','figures','presence_absence_heatmap.pdf'))
# pheatmap_plot
# dev.off()      

enzyme_on_chrom =  Pathways[on_chromatin ==T,Uniprot]
protein_cor = fread(here::here('Datasets','Processed','Chromatin_prot_cor.gz'))
enzyme_cor = protein_cor[str_detect(Protein_1,paste(enzyme_on_chrom,collapse = '|'))|
                str_detect(Protein_2,paste(enzyme_on_chrom,collapse = '|'))]
enzyme_cor = enzyme_cor[order(chrom_cor)]
metabolic_pathways_on_chrom = Pathways[on_chromatin == T, Path_description] |> unique()
BP_enrichment = data.table()
for(metabolic_path in metabolic_pathways_on_chrom){
    enzymes_tmp = Pathways[on_chromatin == T & Path_description ==metabolic_path, Uniprot]
    enzyme_cor_tmp =  enzyme_cor[str_detect(Protein_1,paste(enzymes_tmp,collapse = '|'))|
                                                str_detect(Protein_2,paste(enzymes_tmp,collapse = '|'))]
    enzyme_cor_tmp[,Interactor:= fifelse(str_detect(Protein_1,paste(enzymes_tmp,collapse = '|')),
                                         Protein_2,Protein_1) ]
    enzyme_cor_tmp[,Interactor:=  Interactor|> str_remove_all(';[:print:]*$')]
    enzyme_cor_tmp = enzyme_cor_tmp[,.(mean_cor = median(chrom_cor)), by = Interactor
                                    ][order(mean_cor,decreasing = T)]
    enzyme_cor_tmp = enzyme_cor_tmp |> pull(mean_cor,Interactor)
    if(length(enzyme_cor_tmp)>1){
    edo2 =  clusterProfiler::gseGO(    enzyme_cor_tmp,keyType = 'UNIPROT', 
                           ont= "BP",OrgDb = org.Hs.eg.db)
    if(nrow(edo2@result)>1){
    edo2 = as.data.table(edo2@result)
    edo2[,pathway := metabolic_path]
    BP_enrichment = rbind(BP_enrichment,edo2)}
    }
    # enrichplot::ridgeplot(edo2)
    }
fwrite(BP_enrichment,here::here('Datasets','Processed','Metabolic_corr_BP_enrichment.gz'))
BP_enrichment= fread(here::here('Datasets','Processed','Metabolic_corr_BP_enrichment.gz'))
BP_enrichment[p.adjust<0.05][,.(N_terms = .N), by = pathway] 
BP_enrichment[p.adjust<0.00001 &
              NES>0][str_detect(pathway,'folate')
                ] |> ggplot(aes(x = NES, size = setSize , y= reorder(Description,NES) ))+
    geom_point()+
    ggtitle('Pathways correlating with',
            subtitle = 'One_carbon_by_folate')
ggsave(here::here('Output','figures','BP_one_carbon_by_folate.pdf'))
# fwrite(BP_enrichment[str_detect(pathway,'folate')],here::here('Datasets','datatab'))
protein_cor_piv = imputted |> t() |> as.data.frame() |> tibble::rownames_to_column('sample') |> 
    mutate(tissue = str_remove(sample,'_[:print:]*'))

ggplot(protein_cor_piv, aes(x = Q99623, y = P35232, colour = tissue))+
    geom_point(size =3 )+theme_bw()+
    ggtitle('PHB1 - Prohibitin 1 correlated with PHB2',
            subtitle = '')
ggsave(here::here('Output','figures','PHB_cor.pdf'))
ggplot(protein_cor_piv, aes(x = P55084, y = P40939, colour = tissue))+
    geom_point(size =3 )+theme_bw()+
    ggtitle('HADHB Hydroxyacyl-CoA dehydrogenase trifunctional multienzyme complex subunit beta',
            subtitle = 'HADHA  Hydroxyacyl-CoA dehydrogenase trifunctional multienzyme complex subunit alpha')
ggsave(here::here('Output','figures','HADHB_cor.pdf'))


pathway_for_network = metabolic_pathways_on_chrom |> str_subset('folate')
BP_enrichment_network  =BP_enrichment[pathway == pathway_for_network & NES>1] 
BP_enrichment_network = BP_enrichment_network[,.(ID,Description,core_enrichment)]
BP_enrichment_network = BP_enrichment_network |> tidyr::separate_rows(core_enrichment,sep = '/') |> as.data.table()
BP_enrichment_network[,N_prot:= 1:.N, by = .(core_enrichment)]
BP_enrichment_network = BP_enrichment_network[N_prot ==1]
enzyme_net = Pathways[on_chromatin == T & Path_description ==pathway_for_network, Uniprot]
enzyme_cor_network =  enzyme_cor[str_detect(Protein_1,paste(enzyme_net,collapse = '|'))|
                                 str_detect(Protein_2,paste(enzyme_net,collapse = '|'))]
enzyme_cor_network= rbind(enzyme_cor_network,
                          enzyme_cor[str_detect(Protein_1,paste(BP_enrichment_network$core_enrichment,collapse = '|')) &
                              str_detect(Protein_2,paste(BP_enrichment_network$core_enrichment,collapse = '|'))])
enzyme_cor_network[,Interactor:= fifelse(str_detect(Protein_1,paste(enzyme_net,collapse = '|')),
                                     Protein_2,Protein_1) ]
enzyme_cor_network[,Interactor:=  Interactor|> str_remove_all(';[:print:]*$')]
enzyme_cor_network = inner_join(enzyme_cor_network,BP_enrichment_network, by = c('Interactor' = 'core_enrichment')) |> 
    as.data.table()
network_to_plot = data.table()
proteins_to_plot = unique(c(enzyme_cor_network$Protein_1,enzyme_cor_network$Protein_2))
for(prot in proteins_to_plot){
    enzyme_cor_network_tmp  =  enzyme_cor_network[str_detect(Protein_1,prot)|
                           str_detect(Protein_2,prot)][order(chrom_cor, decreasing = T)]
    enzyme_cor_network_tmp = enzyme_cor_network_tmp |> head(3)
    network_to_plot =  rbind(network_to_plot,enzyme_cor_network_tmp)
}
imputted |> t() |> as.data.frame() |> tibble::rownames_to_column('sample') |> 
    mutate(tissue = str_remove(sample,'_[:print:]*')) |> 
    ggplot(aes (x = P28340 , y= P22102, colour = tissue))+
    geom_point(size = 3)+
    theme_bw()+
    ggtitle('Corr_by tissue',
            subtitle = 'GART DNApol1')
ggsave(here::here('Output','figures','GART_DNAPOl.pdf'))
fwrite(imputted |> t() |> as.data.frame() |> tibble::rownames_to_column('sample') |> 
           mutate(tissue = str_remove(sample,'_[:print:]*')),
       here::here('Datasets','datatable','Fig3G_Cor_plot.csv'))

imputted |> t() |> as.data.frame() |> tibble::rownames_to_column('sample') |> 
    mutate(tissue = str_remove(sample,'_[:print:]*')) |> 
    ggplot(aes (x = P34897 , y= P33993, colour = tissue))+
    geom_point(size = 3)+
    theme_bw()+
    ggtitle('Corr_by tissue',
            subtitle = 'SHMT MCM7 ')
ggsave(here::here('Output','figures','SHMT_MCM7 .pdf'))
# # one carbon by folate network
# 
# visnet <- list(nodes = data.frame(
#     id = c(network_to_plot$Protein_1,network_to_plot$Protein_2) %>% unique()),
#     edges = network_to_plot |> dplyr::select(Protein_1,Protein_2,chrom_cor) %>% purrr::set_names(c("from","to",'chrom_cor')))
# visnet$edges <- visnet$edges %>% 
#     #subset(from %in% genes & to %in% genes) %>% 
#     mutate(
#         # length = round(100*(1 - chrom_cor),digits = 0),
#         # width =round(10*(chrom_cor),digits = 0),
#         folate = if_else(from %in% enzyme_net,from,to)
#     ) %>% dplyr::select(from, to,chrom_cor, folate) %>% distinct() 
# 
# # visnet$edges <-   left_join(mutate(visnet$edges ,
# #                          length = if_else((from %in% Connecting_nodes)|(to %in% Connecting_nodes),T,F)) 
# visnet$nodes <-  visnet$nodes %>% 
#     #dplyr::select(id,label, matches("Int")) %>% 
#     subset(id %in% c(visnet$edges$from,visnet$edges$to)) %>% 
#     mutate(
#         shape = if_else(id %in% enzyme_net,"square","dot"),
#            color = case_when(
#                str_detect(id,paste(enzyme_net,collapse = '|')) ~"#332288",
#                str_detect(
#                    id,paste(BP_enrichment_network[Description =='toxin transport',core_enrichment
#                                                                        ],collapse = '|')) ~"#44AA99",
#                              str_detect(id,paste(BP_enrichment_network[Description =='DNA repair',core_enrichment
#                              ],collapse = '|')) ~"#117733",
#                              str_detect(id,paste(BP_enrichment_network[Description =='chromosome organization',core_enrichment
#                              ],collapse = '|')) ~"#D55E00",
#                              str_detect(id,paste(BP_enrichment_network[Description =='regulation of DNA metabolic process',core_enrichment
#                              ],collapse = '|')) ~"#88CCEE",
#                              str_detect(id,paste(BP_enrichment_network[Description =='cellular response to DNA damage stimulus',core_enrichment
#                              ],collapse = '|')) ~"#DDCC77",
#                              str_detect(id,paste(BP_enrichment_network[Description =='establishment of protein localization',core_enrichment
#                              ],collapse = '|')) ~"#CC6677",
#                              str_detect(id,paste(BP_enrichment_network[Description =='cellular response to stress',core_enrichment
#                              ],collapse = '|')) ~"#AA4499",
#                              str_detect(id,paste(BP_enrichment_network[Description =='regulation of DNA metabolic process',core_enrichment
#                              ],collapse = '|')) ~"#882255",
#                              str_detect(id,paste(BP_enrichment_network[Description =='protein folding',core_enrichment
#                              ],collapse = '|')) ~"white",
#                              str_detect(id,paste(BP_enrichment_network[Description =='regulation of cellular component organization',core_enrichment
#                              ],collapse = '|')) ~"#006CD1",
#                              str_detect(id,paste(BP_enrichment_network[Description =='regulation of developmental process',core_enrichment
#                              ],collapse = '|')) ~"black",
#                             
#                              TRUE~"#CDCDCD"),
#            # pathway = 
#            value = if_else(id %in% enzyme_net,8,5),
#            # `font.size` = (9*value),
#            label = id) %>% 
#     mutate(label = if_else(label %in% enzyme_net,label,NA_character_ ))
# # left_join(median_enrichement, by = c("id"= "ID")) %>% 
# # mutate(color = map2color(Abundance %>% replace_na(0),mypal),
# #        color = if_else(is.na(Abundance), "grey80", color )
# # library(visNetwork)
# # visNetwork(visnet$nodes, visnet$edges, height = "1500px", width = "1500px") %>%
# #     visPhysics(solver = "forceAtlas2Based", 
# #                forceAtlas2Based = list(gravitationalConstant = -10))
# # 
# #     visIgraphLayout()
# library(ggraph);library(igraph)
# net <- igraph::graph_from_data_frame(d=visnet$edges, vertices=visnet$nodes, directed=F) 
#     net <- igraph::simplify(net, remove.multiple = F, remove.loops = T) 
#     
#     ggraph(net, layout = "dh") +
#         geom_edge_link(color="gray30",alpha = 0.2,
#                       aes(width = E(net)$chrom_cor)) +
#         geom_node_point(shape = 21, 
#                         # fill=V(net)$color.border, 
#                         size=V(net)$value,
#                         # labellabel = V(net)$label,
#                         fill = V(net)$color ,stroke = 1) +
#         scale_edge_width(range = c(0.1, 2))+ # control size
#         # scale_edge_size(range = c(1, 20))+ # control size
#         geom_node_text(label = V(net)$label, size=10, color="gray30", repel=T) +
#         theme_void()
#     ggsave(here::here("Output","figures","correlation network.pdf"), width = 15,height = 10)
#     visnet$nodes %>% ggplot(aes(x = id, colour = -log2_FC, y = log2_FC))+
#         geom_point()+
#         scale_colour_gradient2(mid =c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[2] , 
#                                high = c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[3],low = c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[1] )
#     ggsave(here::here("Output","Figures","PPI_network_high_scale.pdf"), width = 15,height = 10)
#     
# 
