library(data.table);library(ggplot2);library(stringr);library('org.Hs.eg.db');library(clusterProfiler);library(purrr);library(tibble)
library(ggalluvial)

# additional analysis chromatome
Jaccard_index_list<-function(list_of_vectors, max_jacc = 0.5,steps = 1,list_type = c("enrichement")){
    
    #removes top step most similar with other pathways and most similar until max_jacc is reached
    #samples as cols
    # list_of_vectors = Complexes_list_dataset
    length_vector <- list_of_vectors %>% lengths()
    row_max = 1
    # top_similar = 1
    # max_jacc = 0.3
    # steps = 1
    while(row_max>max_jacc){
        # list_of_vectors <- enrichment_df_BP %>% subset(`p.adjust`<0.05 ) %>% 
        #   dplyr::select(core_enrichment,Description) %>% subset(!duplicated(Description)) %>% 
        #   pull(core_enrichment,Description)  %>% purrr::map(str_split,"/") %>% flatten()
        index<-map(.x = list_of_vectors, ~.x %>% list(.) %>% rep(.,length(list_of_vectors)) %>% 
                       map2_dbl(.x = .,.y = list_of_vectors,~bayesbio::jaccardSets(.x,.y))) %>% 
            imap_dfr(.x = ., ~set_names(.x,names(list_of_vectors)) %>% enframe(name = "Pathway2",value = "JaccIndex") %>% mutate(Pathway1 = .y)) %>% 
            subset(JaccIndex != 1) %>% as.data.table() 
        setkey(index,Pathway1,Pathway2)
        index <- index[Pathway1>Pathway2]
        row_max <- index$JaccIndex %>% max()
        if(list_type == "complex"){
            index <- index[JaccIndex == row_max
            ][,top_similar:= fifelse(length_vector[Pathway1]>length_vector[Pathway2],Pathway2,Pathway1)]
            
        }else{
            index <- index[JaccIndex == row_max][,top_similar:= fifelse(length_vector[Pathway1]>length_vector[Pathway2],Pathway2,Pathway1)]
        }
        top_similar <- index$top_similar
        # pivot_wider(names_from = "Pathway2", values_from = "JaccIndex") %>% 
        # column_to_rownames("Pathway1") %>% as.matrix()
        # diag(index) <- 0
        # row_max <- index %>% matrixStats::rowMaxs() %>% set_names(row.names(index)) %>% sort(decreasing = T) %>% .[1]
        #top_similar <- index %>% matrixStats::rowSums2() %>% set_names(row.names(index)) %>% sort(decreasing = T) %>% names() %>% .[1:steps]
        # top_similar <- c(top_similar)
        list_of_vectors[which(names(list_of_vectors)%in%top_similar)]<-NULL
        print(row_max)
    }
    list_of_vectors
}
sort_compartments <- function(x){
    strsplit(x,';') |> unlist() |> sort() |> paste(collapse = ';')
}
# defining enzymes 
# downloaded from Epifactor 04/12/2023
Epigenes = fread(here::here('Datasets','Raw','EpiGenes_main.csv'))

#enzymes 
pathways =    fread(here::here('Datasets','Raw','hsa_pathways.txt'))
pathways = pathways[Path_type == 'metabolic']
enzymes =    fread(here::here('Datasets','Raw','KEGG_genes.csv'))
enzymes =  enzymes[pathway %in% pathways$Path_id][,Uniprot        ]
# HPA localisation 
HPA_location = fread(here::here('Datasets','Raw','subcellular_location.tsv'))
HUMAN_9606 <- fread(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"), header = F)  |> 
    as.data.table()
setnames(HUMAN_9606,colnames(HUMAN_9606),c("Uniprot","Type","ID"))
# reading ages of Uniprot
# from https://github.com/marcottelab/Gene-Ages/blob/master/Main/main_HUMAN.csv

ages = fread(here::here('Datasets','Raw','main_HUMAN.csv'))
setnames(ages,'V1','Uniprot')
proteins_on_chromatin = fread(here::here('Datasets','Raw', 'proteins_on_chromatin.csv'))
proteins_on_chromatin[,Uniprot:= stringr::str_remove(Proteins,';[:print:]*$')]
HPA_location_prots_location  = merge(HUMAN_9606[Type =='Gene_Name' & Uniprot %in% proteins_on_chromatin$Uniprot ],
      HPA_location[,.(`Gene name`,`Main location`,`Additional location`)], by.x = 'ID', by.y ='Gene name')
HPA_location_prots_location[,localication:= paste(`Main location`,`Additional location`, sep = ';') |> sort_compartments(), by = Uniprot]
HPA_location_prots_location[,`:=` (Nuclear = str_detect(localication,'Nucl')), by = ID]
proteins_on_chromatin  = ages[proteins_on_chromatin, on  = .(Uniprot), nomatch = NULL]
proteins_on_chromatin[,modeAge := factor(modeAge,levels =colnames(ages)[2:9]) ]
proteins_on_chromatin[,KEGG_enz := fcase(
    Uniprot %in% enzymes & Uniprot %in% Epigenes$UniProt_AC, 'Epigenetic enzyme',
    Uniprot %in% enzymes,'enzyme',
    !(Uniprot %in% enzymes),'non-enzyme')]
proteins_on_chromatin = HPA_location_prots_location[proteins_on_chromatin,on = 'Uniprot'][!is.na(ID)]
proteins_on_chromatin[,Type := fcase(Uniprot %in% enzymes & Nuclear == T,'Nuclear enzyme',
                                     Uniprot %in% enzymes & Nuclear == F,'non-Nuclear enzyme',
                                     !(Uniprot %in% enzymes),'non-enzyme')]
# unique(proteins_on_chromatin, by = "ID")[,.(Perc_prots_on_chrom = mean(on_chromatin)*100,
#                                             N_prots = .N),
#                       by = .(Type,modeAge)] |> View()
# ggplot(aes(x =Type, y = Perc_prots_on_chrom , fill =  modeAge, alpha = log2(N_prots)))+
#     geom_col(position = 'dodge')+
#     ggtitle('Percentage of proteins on chromatin divided by gene-ages',
#             subtitle = '')
ProHD = fread(here::here('Datasets','Raw','Supplementary_Table_S3.csv'))
ProHD[,Interactors := str_remove_all(Coreg_with_protein_X,';[:print:]*$') |> 
    str_remove_all('-[:print:]*$')]
ProHD_functional_enrichment= data.table()
for(i in proteins_on_chromatin[KEGG_enz =='enzyme',Uniprot] ){
    print(i)
    interactors_prohd = ProHD[str_detect(protein_X,i),Interactors] |> unique()
    if(length(interactors_prohd)>4){
    testing = enrichGO(interactors_prohd , OrgDb = 'org.Hs.eg.db', keyType = 'UNIPROT', 
             # qvalueCutoff = 0.05,
         ont = 'ALL',minGSSize = 5, universe = unique(ProHD$Interactors))
    testing = testing@result |> as.data.table()
    if(nrow(testing)>0){
        testing[,Uniprot := i]
        ProHD_functional_enrichment= rbind(ProHD_functional_enrichment,testing)
    }
}
}
fwrite(unique(ProHD_functional_enrichment), 'ProHD_enzyme_functions.csv')
ProHD_functional_enrichment = fread('ProHD_enzyme_functions.csv')
enzyme_nucle = unique(ProHD_functional_enrichment[qvalue<0.01])
enzyme_nucle[str_detect(Description,'nucleolus|nucleolar|nuclear'),.(N_terms = .N), .(Uniprot, ONTOLOGY)] |> 
    ggplot(aes(x = N_terms , y = reorder(Uniprot,N_terms) ))+
    geom_point()+
    facet_wrap('ONTOLOGY')
go_terms = enzyme_nucle[,.(Description,geneID)][,head(.SD,1),by = Description]
to_simplify = go_terms$geneID |> str_split('\\/')
names(to_simplify) = go_terms$Description
simplified = Jaccard_index_list(to_simplify, 0.3)
enzyme_nucle[Description %in% names(simplified),.N, by = .(ONTOLOGY,Description)
             ][,Description_short:= str_sub(Description,1,50)][N>2] |> 
    ggplot(aes(x  = N, y = reorder(Description_short,N)))+
    geom_point()+
    facet_wrap('ONTOLOGY',scales = 'free_y')+
    ggtitle('Processes KEGG enzymes involved in in ProHD')
# Dynamic trafficing
cytotraffic = openxlsx::read.xlsx(here::here('Datasets','Raw','1-s2.0-S0092867423005962-mmc2.xlsx'))

proteins_on_chromatin[,TurboID_cyto_nuc := Uniprot %in% 
                                                       stringr::str_remove(cytotraffic$Uniprot.accession,'-.')]
proteins_on_chromatin |> ggplot(aes(x = on_chromatin, fill= TurboID_cyto_nuc))+
    geom_bar()+
    ggtitle('Proteins which we classify as on-chromatin are enriched for proteins which \n which traffic from cytosol to nucleus (TurboID)')
    
# prot on chrom
library(ggalluvial)
library(ggplot2)
proteins_on_chromatin[on_chromatin==T][,.(Number_of_prots = .N), by = .(KEGG_enz,Nuclear,modeAge)] |> 
    ggplot(aes(y = Number_of_prots, axis2 = Nuclear, 
               axis1 = factor(modeAge,
                              levels =rev(colnames(proteins_on_chromatin)[8:15] )),fill = KEGG_enz)) +
    geom_alluvium( width = 1/12, alpha = 0.8) +
    geom_stratum(width = 1/12, fill = "grey90", color = "grey20") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), colour = 'white') +
    # scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
    scale_fill_manual(values = c('enzyme' = '#FCA331','Epigenetic enzyme'= 'darkred', `non-enzyme` = 'grey80')) +theme_bw()+
    ggtitle("ages of KEGG enzymes in nucleus and HPA localisations")
ggsave(here::here('Output','figures','ages_of_enzymes.pdf'), width = 12, height = 8)
fill_alpha <- function(fill, alpha){
    if (!is.list(fill)) {
        # Happy path for no patterns
        return(alpha(fill, alpha))
    }
    if (is_pattern(fill) || any(vapply(fill, is_pattern, logical(1)))) {
        check_device("patterns", action = "warn")
        fill <- pattern_alpha(fill, alpha)
        return(fill)
    } else {
        # We are either dealing with faulty fill specification, or we have a legend
        # key that is trying to draw a single colour. It can be given that colour
        # as a list due to patterns in other keys.
        msg <- paste0(
            "{.field fill} must be a vector of colours or list of ",
            "{.cls GridPattern} objects."
        )
        # If single colour list, try applying `alpha()`
        fill <- try_fetch(
            Map(alpha, colour = fill, alpha = alpha),
            error = function(cnd) {
                cli::cli_abort(msg, call = expr(fill_alpha()))
            }
        )
        # `length(input)` must be same as `length(output)`
        if (!all(lengths(fill) == 1)) {
            cli::cli_abort(msg)
        }
        return(unlist(fill))
    }
}
fwrite(proteins_on_chromatin[on_chromatin==T], here::here('Datasets','datatable','Fig2A_geneages.csv'))
