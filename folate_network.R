# network mapping chromatome
# read core chromatome
proteins_on_chromatin = fread(here::here('Datasets','Raw', 'proteins_on_chromatin.csv'))
proteins_on_chromatin[,Uniprot:= str_remove_all(Proteins,';[:print:]*')]
prots_to_network = proteins_on_chromatin[on_chromatin ==T,Uniprot] |> unique()

# hsa pathways
pathways =    fread(here::here('Datasets','Raw','hsa_pathways.txt'))
pathways = pathways[Path_type == 'metabolic']
enzymes =    fread(here::here('Datasets','Raw','KEGG_genes.csv'))
enzyme_pathway_map = pathways[enzymes, on= c('Path_id' = 'pathway')
][!is.na(Path_type)][,.(Path_description,ID)] |> unique()
enzymes =  enzymes[pathway %in% pathways$Path_id][,ID]
Pathways = enzyme_pathway_map[enzymes_on_chromatin, on = 'ID']

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
# fwrite(BP_enrichment,here::here('Datasets','Processed','Metabolic_corr_BP_enrichment.gz'))
BP_enrichment= fread(here::here('Datasets','Processed','Metabolic_corr_BP_enrichment.gz'))
BP_enrichment[p.adjust<0.00001 &
                  NES>0][str_detect(pathway,'folate')
                  ] |> ggplot(aes(x = NES, fill = setSize , y= reorder(Description,NES) ))+
    geom_col()+
    ggtitle('Pathways correlating with',
            subtitle = 'One_carbon_by_folate')
ggsave(here::here('Output','figures','BP_one_carbon_by_folate.pdf'))
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
BP_enrichment_network  =BP_enrichment[pathway == pathway_for_network & NES>1] |> head(7)
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
network_to_plot  = network_to_plot[chrom_cor>0.3]
imputted |> t() |> as.data.frame() |> tibble::rownames_to_column('sample') |> 
    mutate(tissue = str_remove(sample,'_[:print:]*')) |> 
    ggplot(aes (x = P28340 , y= P22102, colour = tissue))+
    geom_point(size = 3)+
    theme_bw()+
    ggtitle('Corr_by tissue',
            subtitle = 'GART DNApol1')
ggsave(here::here('Output','figures','GART_DNAPOl.pdf'))

imputted |> t() |> as.data.frame() |> tibble::rownames_to_column('sample') |> 
    mutate(tissue = str_remove(sample,'_[:print:]*')) |> 
    ggplot(aes (x = P34897 , y= P33993, colour = tissue))+
    geom_point(size = 3)+
    theme_bw()+
    ggtitle('Corr_by tissue',
            subtitle = 'SHMT MCM7 ')
ggsave(here::here('Output','figures','SHMT_MCM7 .pdf'))
# one carbon by folate network

visnet <- list(nodes = data.frame(
    id = c(network_to_plot$Protein_1,network_to_plot$Protein_2) %>% unique()),
    edges = network_to_plot |> dplyr::select(Protein_1,Protein_2,chrom_cor) %>% purrr::set_names(c("from","to",'chrom_cor')))
visnet$edges <- visnet$edges %>% 
    #subset(from %in% genes & to %in% genes) %>% 
    mutate(
        # length = round(100*(1 - chrom_cor),digits = 0),
        # width =round(10*(chrom_cor),digits = 0),
        folate = if_else(from %in% enzyme_net,from,to)
    ) %>% dplyr::select(from, to,chrom_cor, folate) %>% distinct() 

# visnet$edges <-   left_join(mutate(visnet$edges ,
#                          length = if_else((from %in% Connecting_nodes)|(to %in% Connecting_nodes),T,F)) 
visnet$nodes <-  visnet$nodes %>% 
    #dplyr::select(id,label, matches("Int")) %>% 
    subset(id %in% c(visnet$edges$from,visnet$edges$to)) %>% 
    mutate(shape = if_else(id %in% enzyme_net,"square","dot"),
           color = case_when(
               str_detect(id,paste(enzyme_net,collapse = '|')) ~"#332288",
               str_detect(id,paste(BP_enrichment_network[Description =='toxin transport',core_enrichment
               ],collapse = '|')) ~"white",
               str_detect(id,paste(BP_enrichment_network[Description =='DNA repair',core_enrichment
               ],collapse = '|')) ~"#117733",
               str_detect(id,paste(BP_enrichment_network[Description =='chromosome organization',core_enrichment
               ],collapse = '|')) ~"#D55E00",
               str_detect(id,paste(BP_enrichment_network[Description =='regulation of DNA metabolic process',core_enrichment
               ],collapse = '|')) ~"#88CCEE",
               str_detect(id,paste(BP_enrichment_network[Description =='cellular response to DNA damage stimulus',core_enrichment
               ],collapse = '|')) ~"#DDCC77",
               str_detect(id,paste(BP_enrichment_network[Description =='establishment of protein localization',core_enrichment
               ],collapse = '|')) ~"#CC6677",
               str_detect(id,paste(BP_enrichment_network[Description =='cellular response to stress',core_enrichment
               ],collapse = '|')) ~"#AA4499",
               str_detect(id,paste(BP_enrichment_network[Description =='regulation of DNA metabolic process',core_enrichment
               ],collapse = '|')) ~"#882255",
               # str_detect(id,paste(BP_enrichment_network[Description =='protein folding',core_enrichment
               # ],collapse = '|')) ~"white",
               # str_detect(id,paste(BP_enrichment_network[Description =='regulation of cellular component organization',core_enrichment
               # ],collapse = '|')) ~"#006CD1",
               # str_detect(id,paste(BP_enrichment_network[Description =='regulation of developmental process',core_enrichment
               # ],collapse = '|')) ~"black",
               
               TRUE~"#CDCDCD"),
           # pathway = 
           value = if_else(id %in% enzyme_net,8,5),
           # `font.size` = (9*value),
           label = id) %>% 
    mutate(label = if_else(label %in% enzyme_net,label,NA_character_ ))
# left_join(median_enrichement, by = c("id"= "ID")) %>% 
# mutate(color = map2color(Abundance %>% replace_na(0),mypal),
#        color = if_else(is.na(Abundance), "grey80", color )
# library(visNetwork)
# visNetwork(visnet$nodes, visnet$edges, height = "1500px", width = "1500px") %>%
#     visPhysics(solver = "forceAtlas2Based", 
#                forceAtlas2Based = list(gravitationalConstant = -10))
# 
#     visIgraphLayout()
library(ggraph);library(igraph)
net <- igraph::graph_from_data_frame(d=visnet$edges[from %in% enzyme_net| to %in% enzyme_net], vertices=visnet$nodes, directed=F) 
net <- igraph::simplify(net, remove.multiple = F, remove.loops = T) 
layout_tmp = data.table(uniprot = visnet$nodes$id,
                        hit = visnet$nodes$value)
layout_tmp[,`:=`(x = fcase(uniprot == 'P11586',1,
    uniprot == 'P22102',1,
    uniprot == 'P34897',-1,
    uniprot == 'P31939',-1,
    str_detect(uniprot,paste(BP_enrichment_network[Description =='toxin transport',core_enrichment
    ],collapse = '|')) ,rnorm(1,-0.7,0.13),
    str_detect(uniprot,paste(BP_enrichment_network[Description =='DNA repair',core_enrichment
    ],collapse = '|')),rnorm(1,0.9,0.13),
    str_detect(uniprot,paste(BP_enrichment_network[Description =='chromosome organization',core_enrichment
    ],collapse = '|')),rnorm(1,-0.01,0.13),
    str_detect(uniprot,paste(BP_enrichment_network[Description =='regulation of DNA metabolic process',core_enrichment
    ],collapse = '|')),rnorm(1,-0.3,0.13),
    str_detect(uniprot,paste(BP_enrichment_network[Description =='cellular response to DNA damage stimulus',core_enrichment
    ],collapse = '|')) ,rnorm(1,0.2,0.10)),     
    y = fcase(uniprot == 'P11586',-1,
    uniprot == 'P22102',1,
    uniprot == 'P34897',1,
    uniprot == 'P31939',-1,
    str_detect(uniprot,paste(BP_enrichment_network[Description =='toxin transport',core_enrichment
    ],collapse = '|')) ,rnorm(1,-0.1,0.13),
    str_detect(uniprot,paste(BP_enrichment_network[Description =='DNA repair',core_enrichment
    ],collapse = '|')),rnorm(1,-0.3,0.13),
    str_detect(uniprot,paste(BP_enrichment_network[Description =='chromosome organization',core_enrichment
    ],collapse = '|')),rnorm(1,+0.7,0.13),
    str_detect(uniprot,paste(BP_enrichment_network[Description =='regulation of DNA metabolic process',core_enrichment
    ],collapse = '|')),rnorm(1,-0.78,0.08),
    str_detect(uniprot,paste(BP_enrichment_network[Description =='cellular response to DNA damage stimulus',core_enrichment
    ],collapse = '|')) ,rnorm(1,-0.1,0.13))
    ), by = uniprot]
fwrite(layout_tmp, here::here('Datasets','datatable','Fig3F_network_nodes.csv'))

ggraph(net, layout_tmp[,.(x,y)]) +
    geom_edge_link(color="gray30",alpha = 0.2,
                   aes(width = E(net)$chrom_cor)) +
    geom_node_point(shape = 21, 
                    # fill=V(net)$color.border, 
                    size=V(net)$value,
                    # labellabel = V(net)$label,
                    fill = V(net)$color ,stroke = 1) +
    scale_edge_width(range = c(0.1, 2))+ # control size
    # scale_edge_size(range = c(1, 20))+ # control size
    geom_node_text(label = V(net)$label, size=10, color="gray30", repel=T) +
    theme_void()
ggsave(here::here("Output","figures","correlation network_corners.pdf"), width = 7,height = 5)7

# visnet$nodes %>% ggplot(aes(x = id, colour = -log2_FC, y = log2_FC))+
#     geom_point()+
#     scale_colour_gradient2(mid =c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[2] , 
#                            high = c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[3],low = c(scales::muted("blue"),"#f2f5f7","#9e3a3a")[1] )
# ggsave(here::here("Output","Figures","PPI_network_high_scale.pdf"), width = 15,height = 10)

