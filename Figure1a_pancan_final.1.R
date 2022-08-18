rm(list=ls())


setwd("~/Documents/ccbr1144_gameiro/analysis/8")

library(maftools)
library(ggplot2)
library(openxlsx)
# maf_file="/Users/tandonm/Documents/ccbr1144_gameiro/new_pipeline_out/somatic_paired/SNP_Indels/mutect2_out/all_somatic_variants.maf"
# filtered_maf_file="data/mutect2.filtered.maf"

source("~/Documents/helper_functions/helper_functions.tcga.R")
# file.edit("~/Documents/helper_functions/helper_functions.oncoplot.R")
# file.edit("~/Documents/helper_functions/helper_functions.tcga.R")

out_folder=file.path(getwd(),"plots","pancan_1")
tcga_data_source="~/Downloads/globus/1155/6_tcga/data"
rds_file1=file.path(tcga_data_source, "tcga_mafs.list.Rds")
rds_file2=file.path(tcga_data_source, "tcga_mafs.merged.Rds")

if (file.exists(rds_file2)) {
  mymafs <- readRDS(rds_file1)
  merged_tcga_maf <- readRDS(rds_file2)
} else {
   tcga_datasets <- c("BLCA","BRCA","CESC","COAD","ESCA",
  "GBM","HNSC","LIHC","SKCM","LUSC","OV","PAAD","PRAD","KIRC","SARC","THCA")
  #tcga_datasets <- TCGAbiolinks::getGDCprojects()$tumor
  mymafs <- lapply(tcga_datasets, function(ds) {
    returnval=NULL
    # ds="foo"
    returnval <- tryCatch( {
      get_tcga_data(ds, save_folder = tcga_data_source)
      }, 
      error=function(cond) {
        message(paste0("Download failed for ",ds))
        message(cond)
      }
    )
    return(returnval)
  })
  names(mymafs) <- tcga_datasets
  mymafs <- mymafs[!unlist(lapply(mymafs, is.null))]
  merged_tcga_maf <- merge_mafs(mymafs)
  saveRDS(mymafs, file=rds_file1)
  saveRDS(merged_tcga_maf, file=rds_file2)
}

out_folder=file.path(getwd(),"plots","TCGA_pancan_4")
if (!dir.exists(out_folder)){dir.create(out_folder, recursive = T)}
  
# pathways_df_file="data/pathwaysdf.pancan.txt"
# if (file.exists(pathways_df_file)) {
#   pathwaysdf <- read.table(pathways_df_file, sep="\t", header = T)
# } else {
  genes_file="~/Downloads/globus/1155/lab_data/selected_genes.xlsx"
  pathwaysdf <- read.xlsx(genes_file, sheet="PanCan3")
  # pathwaysdf$mgi_symbol <- pathwaysdf$Hugo_Symbol
  
  # mgi_homologs_file="data/HOM_MouseHumanSequence.rpt"
  # if (!file.exists(mgi_homologs_file)) {
  #   download.file("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", destfile = mgi_homologs_file)
  # }
  # library(data.table)
  # mouse2human <- fread(mgi_homologs_file)
  # 
  # unique_mgi <- sort(unique(pathwaysdf$mgi_symbol))
  # myids <- mouse2human$`DB Class Key`[match(unique_mgi, mouse2human$Symbol)]
  # human_ids <- unlist(lapply(myids, function(id) {
  #   returnval=NA
  #   if (!is.na(id)) {
  #     returnval=paste0(mouse2human$Symbol[mouse2human$`DB Class Key`==id & mouse2human$`Common Organism Name`=="human"],collapse=",")
  #   }
  #   return(returnval)
  # }))
  # compare_df <- data.frame(mgi_symbol=unique_mgi, human_mapped=human_ids)
  # write.xlsx(compare_df, file = "orthologs.xlsx")
  # names(mouse2human)
  
  # Basic function to convert mouse to human gene names
  # From: https://rjbioinformatics.com/2016/10/14/converting-mouse-to-human-gene-names-with-biomart-package/
  
#   require("biomaRt")
#   human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")#, host="https://uswest.ensembl.org")
#   mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")#, host="https://uswest.ensembl.org")
#   
#   mouse_ids <- biomaRt::getBM(c("mgi_symbol","ensembl_gene_id"), filters="mgi_symbol", values=pathwaysdf$mgi_symbol, mart = mouse)
#   human_ids <- biomaRt::getBM(c("mgi_symbol","ensembl_gene_id"), filters="mgi_symbol", values=pathwaysdf$mgi_symbol, mart = mouse)
#   
#   genesV2 = getLDS(attributes = c("mgi_symbol","ensembl_gene_id"), 
#                    filters = "ensembl_gene_id", 
#                    values = mouse_ids$ensembl_gene_id , 
#                    mart = mouse, 
#                    attributesL = c("hgnc_symbol","ensembl_gene_id"), 
#                    martL = human, 
#                    uniqueRows=T)
#   # biomaRt::getBM(c("mgi_symbol","ensembl_gene_id"), filters="mgi_symbol", values=pathwaysdf$mgi_symbol[54], mart = mouse)
#   # biomaRt::getBM(c("mgi_symbol","ensembl_gene_id"), filters="mgi_symbol", values=pathwaysdf$mgi_symbol[54], mart = mouse)
#   # humanx <- genesV2[, 2]
#   
#   pathwaysdf$Hugo_Symbol <- genesV2[match(pathwaysdf$mgi_symbol, genesV2[,"MGI.symbol"]),"HGNC.symbol"]
#   write.table(pathwaysdf,file=pathways_df_file,sep="\t",quote=F, row.names=F)
# }

custom_genes_df <- data.frame(Reason=stringr::str_wrap(gsub("_"," ",gsub("HALLMARK_","",pathwaysdf$Reason)), width=20),
                       Hugo_Symbol=pathwaysdf$Hugo_Symbol, 
                       stringsAsFactors = F)
custom_genes_df <- unique(custom_genes_df[!is.na(custom_genes_df$Hugo_Symbol),])




# subset_mutgenes <- pathwaysdf$Hugo_Symbol[pathwaysdf$mgi_symbol %in% c("B2m","Jak1","Psmb9")]
# subset_mafs <- lapply(mymafs, subsetMaf, genes=subset_mutgenes, mafObj=F)
# head(subset_mafs[[1]])



source("https://raw.githubusercontent.com/mtandon09/mt_helpers/main/helper_functions.oncoplot.R")

ribbon_out_folder <- file.path(out_folder, "cooccurence")
if (!dir.exists(ribbon_out_folder)){dir.create(ribbon_out_folder, recursive = T)}

merged_tcga_maf="~/Downloads/globus/1155/6_tcga/data/tcga_mafs.all.Rds"
merged_tcga_maf=readRDS(merged_tcga_maf)
merged_tcga_maf=merge_mafs(merged_tcga_maf)
plotfile=("pancan_ribbon.png")
make_single_ribbon_plot(
          merged_tcga_maf, 
          onco_genes=NULL, 
          topN=20,
          pval_high=0.01,  ## All interactions with less than this p-value will be shown
          pval_low=0.001,  ## Links with p-value less than this will be highlighted with a dashed border
          plot_frac_mut_axis=TRUE,  ## Whether or not to draw a numerical axis on the perimeter
          rotate_plot_degrees=0,   ## For custom rotation
          shrink_factor=1.3, # Higher = more shrinkage; use to control whitespace (or lack thereof) around figure 
          scale_ribbon_to_fracmut=TRUE,  ## Whether or not to scale ribbon widths to their frequency
          save_name=plotfile, 
          ribbon_color=NULL, 
          progress_func=NULL,
          sig_colors=NULL,   ## Vector of 4 colors for coloring significance
          gene_colors=NULL   ## color(s) for gene blocks
)


plotfile=file.path(ribbon_out_folder, "pancan_ribbon.all_selected_genes.png")
make_single_ribbon_plot(
          merged_tcga_maf, 
          onco_genes=unique(custom_genes_df$Hugo_Symbol), 
          topN=40,
          pval_high=0.01,  ## All interactions with less than this p-value will be shown
          pval_low=0.001,  ## Links with p-value less than this will be highlighted with a dashed border
          plot_frac_mut_axis=TRUE,  ## Whether or not to draw a numerical axis on the perimeter
          rotate_plot_degrees=0,   ## For custom rotation
          shrink_factor=1.2, # Higher = more shrinkage; use to control whitespace (or lack thereof) around figure 
          scale_ribbon_to_fracmut=TRUE,  ## Whether or not to scale ribbon widths to their frequency
          save_name=plotfile, 
          ribbon_color=NULL, 
          progress_func=NULL,
          sig_colors=NULL,   ## Vector of 4 colors for coloring significance
          gene_colors=NULL   ## color(s) for gene blocks
)

gene_lists <- split(custom_genes_df$Hugo_Symbol, custom_genes_df$Reason)

lapply(names(gene_lists), function(listname){
  currgenes <- gene_lists[[listname]]
  plotfile=file.path(paste0("pancan_ribbon.",listname,".png"))
  make_single_ribbon_plot(
    merged_tcga_maf, 
    onco_genes=currgenes, 
    topN=40,
    pval_high=0.05,  ## All interactions with less than this p-value will be shown
    pval_low=0.01,  ## Links with p-value less than this will be highlighted with a dashed border
    plot_frac_mut_axis=TRUE,  ## Whether or not to draw a numerical axis on the perimeter
    rotate_plot_degrees=0,   ## For custom rotation
    shrink_factor=1.2, # Higher = more shrinkage; use to control whitespace (or lack thereof) around figure 
    scale_ribbon_to_fracmut=TRUE,  ## Whether or not to scale ribbon widths to their frequency
    save_name=plotfile, 
    ribbon_color=NULL, 
    progress_func=NULL,
    sig_colors=NULL,   ## Vector of 4 colors for coloring significance
    gene_colors=NULL   ## color(s) for gene blocks
  )
})



if (! exists(mymafs)){
  mymafs <- readRDS(rds_file1)
}

gene_summaries <- do.call(rbind, lapply(names(mymafs), function(mafname) {
  # mafname="BLCA"
  mafobj <- mymafs[[mafname]]
  returnval <- mafobj@gene.summary
  returnval <- returnval[returnval$Hugo_Symbol %in% custom_genes_df$Hugo_Symbol, c("Hugo_Symbol","AlteredSamples")]
  returnval$cohort_pct <- returnval$AlteredSamples/as.numeric(mafobj@summary$summary[3])
  returnval$dataset <- paste0(stringr::str_pad(mafname, width = 5,side = "left"),
                              # "[",
                              stringr::str_pad(mafobj@summary[["summary"]][[3]], width = 4,side = "left")
                              # "]"
                              )
  returnval$plot_value <- round(returnval$cohort_pct*100, 1)
  return(returnval)
  })
)
gene_summaries$dataset <- gsub("([[:digit:]]+)", "[\\1]", gene_summaries$dataset)

# mhc_genes <- custom_genes_df$Hugo_Symbol[custom_genes_df$Reason=="MHC"]
# this <- gene_summaries[gene_summaries$Hugo_Symbol]


library(dplyr)
gene_summaries <- left_join(gene_summaries, custom_genes_df)


library(ggplot2)
library(reshape2)
# plotdata <- melt(gene_summaries, id.vars=c("dataset","Reason","cohort_pct"), measure.vars="Hugo_Symbol")
plotdata <- gene_summaries
plotdata$dataset <- factor(plotdata$dataset, levels=sort(unique(plotdata$dataset), decreasing = T))

myplot <- ggplot(plotdata) +
  facet_wrap(~Reason) +
  geom_col(aes(x=plot_value, y=dataset, fill=Hugo_Symbol))

ggsave(file.path(out_folder,"pancan_test.1.pdf"),plot=myplot, width=16, height=16)


gene_text_data <- split(plotdata$Hugo_Symbol, plotdata$Reason)
# prim_colors <- setNames(rainbow(length(gene_text_data)), gene_text_data)
library(ggsci)
# prim_colors <- setNames(pal_rickandmorty((palette = c("schwifty")))(length(gene_text_data)), names(gene_text_data))
# prim_colors <- setNames(pal_aaas((palette = c("default")))(length(gene_text_data)), names(gene_text_data))
prim_colors <- setNames(pal_d3((palette = c("category10")))(length(gene_text_data)), names(gene_text_data))
# prim_colors <- setNames(pal_futurama((palette = c("planetexpress")))(length(gene_text_data)), names(gene_text_data))
# prim_colors <- setNames(pal_jco((palette = c("default")))(length(gene_text_data)), names(gene_text_data))
# prim_colors <- setNames(pal_npg((palette = c("nrc")))(length(gene_text_data)), names(gene_text_data))
# library(ggthemes)
# prim_colors <- setNames(gdocs_pal()(length(gene_text_data)), names(gene_text_data))
# prim_colors <- setNames(calc_pal()(length(gene_text_data)), names(gene_text_data))


xmax=round(max(plotdata$plot_value) * 3, 2)
gene_text <- do.call(rbind, lapply(names(gene_text_data), function(x){
  # x=names(gene_text_data)[1]
  mygenes=sort(unique(gene_text_data[[x]]))
  start_color=prim_colors[x]
  end_color="grey90"
  
  colfunc <- colorRampPalette(c(start_color, end_color))
  gene_colors <- colfunc(length(mygenes)+1)[1:length(mygenes)]
  names(gene_colors) <- mygenes
  data.frame(gene=mygenes, Reason=x, 
             # xpos=0.15, 
             xpos=xmax*0.8,
             # ypos=1:length(mygenes) + 1, gene_color=gene_colors)
             ypos=(1:length(mygenes) + 1)*1.5, gene_color=gene_colors)
  })
)

custom_colors <- setNames(gene_text$gene_color, gene_text$gene)
myplot <- ggplot(plotdata) +
  facet_wrap(~Reason) +
  geom_col(aes(x=plot_value, y=dataset, fill=Hugo_Symbol), color="black", size=0.3, show.legend=F) +
  # geom_text(data=gene_text, aes(x=xpos, y=ypos, label=gene, color=gene),fontface=2,inherit.aes = F, show.legend=F) +
  geom_label(data=gene_text, aes(x=xpos, y=ypos, label=gene, color=gene),fill="white", size=4,fontface=2,inherit.aes = F, show.legend=F) +
  scale_fill_manual(values=custom_colors) +
  scale_color_manual(values=custom_colors) +
  xlim(c(0, xmax)) +
  xlab("Patients with mutations (% of cohort)") + ylab("TCGA Cancer Type") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold.italic", size=16)
  )
# myplot
ggsave(file.path(out_folder,"pancan_test.2.pdf"),plot=myplot, width=10, height=10)



gene_text_data <- split(plotdata, plotdata$Reason)
newplotdata <- do.call(rbind, lapply(gene_text_data, function(currdf) {
  # currdf <- gene_text_data[[1]]
  currgenes <- sort(unique(currdf$Hugo_Symbol))
  geneidx <- setNames(paste0("gene_",1:length(currgenes)), currgenes)
  
  currdf$gene_idx <- geneidx[currdf$Hugo_Symbol]
  return(currdf)
}))

gene_text <- do.call(rbind, lapply(split(newplotdata, newplotdata$Reason), function(currdf) {
  mycols <- c("Hugo_Symbol","Reason","gene_idx")
  unique(currdf[,..mycols])
}))


uniq_genes <- sort(unique(newplotdata$gene_idx))

# library(RColorBrewer)
# custom_colors <- setNames(brewer.pal("Set1", n=length(uniq_genes)), uniq_genes)
library(ggsci)
# custom_colors <- setNames(pal_jco((palette = c("default")))(length(uniq_genes)), uniq_genes)
custom_colors <- setNames(pal_npg((palette = c("nrc")))(length(uniq_genes)), uniq_genes)

# gene_text$xpos=round(max(newplotdata$plot_value) * 1.1, 2)
# gene_text$ypos=as.numeric(unlist(lapply(strsplit(gene_text$gene_idx,"_"),"[[",2)))+1
gene_text$xpos=xmax*0.8
gene_text$ypos=(as.numeric(unlist(lapply(strsplit(gene_text$gene_idx,"_"),"[[",2)))+1)*1.5

myplot <- ggplot(newplotdata) +
  facet_wrap(~Reason) +
  geom_col(aes(x=plot_value, y=dataset, fill=gene_idx), color="black", size=0.4, show.legend=F) +
  # geom_text(data=gene_text, aes(x=xpos, y=ypos, label=Hugo_Symbol, color=gene_idx),fontface=2,inherit.aes = F, show.legend=F) +
  geom_label(data=gene_text, aes(x=xpos, y=ypos, label=Hugo_Symbol, color=gene_idx),fill="white", size=4,fontface=2,inherit.aes = F, show.legend=F) +
  scale_fill_manual(values=custom_colors) +
  scale_color_manual(values=custom_colors) +
  xlim(c(0, xmax)) +
  xlab("Patients with mutations (% of cohort)") + ylab("TCGA Cancer Type") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold.italic", size=16)
  )
  # theme_linedraw(base_size = 14)
# myplot
ggsave(file.path(out_folder,"pancan_test.3.pdf"),plot=myplot, width=10, height=10)

# sort(table(gene_summaries$))




