rm(list=ls())

## Path to current working directory
# wd="/data/CCBR/projects/ccbr1144/manuscript_data/code_for_figures/exome"
wd="~/Documents/ccbr1144_gameiro/analysis/minnar_manuscript/prep_for_biowulf"
setwd(wd)

library(maftools)
library(dplyr)
## Contains custom functions for filtering and plotting MAF file data
source("scripts/helper_functions.R")

output_dir="figures/2"
if (!dir.exists(output_dir)){ dir.create(output_dir, recursive = T)}

## Path to the unfiltered variants from Dragen, converted to MAF using vcf2maf, and then gzipped
raw_maf_file="data/minnar_variants_maftools.maf.gz"
mafobj <- filter_maf_chunked(raw_maf_file, chunk_lines=100000, grep_vcf_filter_col="PASS|.",
                   # n_callers=ncallers, variant_caller="consensus",
                   #savename = filtered_maf_file,
                   non_silent_only=F,
                   # t_alt_min=1, t_alt_max=1e12, t_depth_min=2,t_depth_max=1e12, 
                   t_alt_min=10, t_alt_max=1e12, t_depth_min=50,t_depth_max=1e12, 
                   tumor_freq_min=0.05,tumor_freq_max=1,
                   n_alt_min=0,n_alt_max=2,n_depth_min=0,n_depth_max=1e12,
                   norm_freq_min=0,norm_freq_max=1
                   )

targets_coverage=49569003
## This value can be computed from the Agilent targets bed file using this function
## found in the helper_functions.R script
# targets_file="~/Documents/exome_bed_files/SureSelect_mm10.bed"
# targets_coverage <- compute_exome_coverage(targets_file)

library(ggplot2)
custom_burdenplot <- file.path(output_dir,"supp_fig_1c.tiff")
mysamples <- c(TC1a9="TC1A9",CMT64="CMT64",CMT64="RVP3")
burdenplot <- make_burden_plot(mafobj, mb_covered = targets_coverage/1e6, sample_order = mysamples,
                               plotType = "Barplot", add_median = 10)
burdenplot <- burdenplot + theme_linedraw(base_size = 14)
ggsave(custom_burdenplot, width=6, height=4, plot = burdenplot, dpi=320)





## Get genes of interest from msigdb
library(msigdbr)
m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol, gs_subcat)
## Select just ones matching "interferon gamma"
pathwaysdf <- data.frame(m_t2g)
mypaths <- pathwaysdf[grepl("INTERFERON_GAMMA",pathwaysdf$gs_name, ignore.case = T),]

## Set for make_oncoplot function (and format stuff nicely)
genes_df <- data.frame(Reason=stringr::str_wrap(gsub("_"," ",gsub("HALLMARK_","",mypaths$gs_name)), width=10),
                       Hugo_Symbol=mypaths$gene_symbol, 
                       stringsAsFactors = F)

maf.filtered <- mafobj
genes_for_oncoplot <- genes_df
include_all = T
show_sample_names=T

## Make matrix of mutations per sample and put samples in specified order
oncomat <- createOncoMatrix(maf.filtered, g=genes_for_oncoplot$Hugo_Symbol, add_missing = F)$oncoMatrix
custom_order <- match(mysamples, colnames(oncomat), nomatch=0)
oncomat <- oncomat[,custom_order]
oncomat.plot <- oncomat

### Set the height of the plot based on number of genes
onco_height=NULL
if (is.null(onco_height)) {
  onco_height=max(round(0.2*nrow(oncomat.plot),0),6)
}

oncomat.plot <- gsub("_"," ",oncomat.plot)

onco_height=max(round(0.15*nrow(oncomat.plot),0),6)
# onco_width=max(c(onco_height*0.75,max(round(0.4*ncol(oncomat.plot),0))))
onco_width=5

variant_type_data <- data.frame(maf.filtered@variant.classification.summary)
rownames(variant_type_data) <- variant_type_data$Tumor_Sample_Barcode
colnames(variant_type_data) <- gsub("_"," ",colnames(variant_type_data))
# variant_type_data <- variant_type_data[,c(-1,-ncol(variant_type_data)), drop=F]
variant_type_data <- variant_type_data[,ncol(variant_type_data), drop=F]
variant_type_data <- variant_type_data[match(colnames(oncomat.plot), rownames(variant_type_data)),
                                       rev(order(colSums(variant_type_data))), drop=F]
# browser()
# var_anno_colors <- mutation_colors[match(colnames(variant_type_data), names(mutation_colors))]
library(ComplexHeatmap)
mb_covered <- targets_coverage/1e6
variant_type_data <- variant_type_data/mb_covered
burden_plotarea <- unit(1,"inches")
# ha = HeatmapAnnotation("Mutations\n/MB\n\n\n" = anno_barplot(variant_type_data,
burden_ha = HeatmapAnnotation("Mutations/MB" = anno_barplot(variant_type_data,
                                             # gp = gpar(fill = var_anno_colors),
                                             gp = gpar(fill = "grey40", color="black",lwd=2),
                                             axis_param = list(gp=gpar(fontsize=11)),
                                             bar_width = 0.9,
                                             border = F,
                                             height=burden_plotarea),
                       annotation_name_side = "left",
                       annotation_name_gp = gpar(fontsize=12),
                       annotation_name_offset = unit(0.2, "npc"),
                       annotation_name_rot = 90
                       
)

# pctanno_data <- paste0(prettyNum(rowSums(oncomat.plot!="")/ncol(oncomat.plot)*100, digits=0),"%")
# pctanno_data <- paste0("  ",pctanno_data,"  ")
pctanno_data <- paste0(rowSums(oncomat.plot!=""),"/", ncol(oncomat.plot))
pctanno_plotarea <- unit(0.5,"inches")
pctanno_ha = HeatmapAnnotation("simplename" = anno_text(pctanno_data,
# pctanno_ha = rowAnnotation("Frequency\nof each gene\nmutation" = anno_text(pctanno_data,
                                             gp = gpar(fontsize = 11),
                                             show_name = T
                                             # gp = gpar(fill = "grey40", color="black",lwd=2),
                                             # axis_param = list(gp=gpar(fontsize=9)),
                                             # bar_width = 0.9,
                                             # border = F,
                                        #      width=pctanno_plotarea
                                        ),
                           which="row",
                           width=pctanno_plotarea,
                           annotation_label=c("Frequency\nof each gene\nmutation"),
                           show_annotation_name = F,
                           annotation_name_rot = 0,
                           annotation_name_side = "top"
                       
)


# draw(pctanno_ha)
# oncomat.plot <- oncomat.plot[rev(order(oncomat.plot[,mysamples[1]], oncomat.plot[,mysamples[2]], oncomat.plot[,mysamples[3]])),]
# split_idx=rep("IFN-gamma Response", nrow(oncomat.plot))
# pathway_label=expression(paste0("IFN-",gamma," Response"))
pathway_label="IFN-γ Response"
split_idx=rep(pathway_label, nrow(oncomat.plot))
names(split_idx) <- rownames(oncomat.plot)
# split_idx=data.frame(rep("IFN-gamma Response", nrow(oncomat.plot)))
# colnames(split_idx) <- "Pathway"
# rownames(split_idx) <- rownames(oncomat.plot)
split_colors=list(Pathway=c("red3"))
names(split_colors$Pathway) <- pathway_label
oncomat.plot <- oncomat.plot[order(oncomat.plot[,1],oncomat.plot[,2], rownames(oncomat.plot)),]
onco_base_default <- oncoPrint(oncomat.plot, alter_fun = alter_fun, 
                               col=mutation_colors, 
                               # row_order=1:nrow(oncomat.plot),
                               column_order = 1:ncol(oncomat.plot),
                               name="oncoplot",
                               row_title="",
                               row_names_gp = gpar(fontsize = 14),
                               column_names_gp = gpar(fontsize = 12),
                               show_pct = F,
                               # pct_side = "right",
                               row_split=split_idx,
                               # bottom_annotation = myanno,
                               left_annotation = rowAnnotation(Pathway = split_idx, col=split_colors, annotation_width = unit(0.3, "mm"), show_annotation_name=F),
                               # right_annotation = NULL,
                               # right_annotation = pctanno_ha,
                               top_annotation = burden_ha,
                               # row_names_side = "left",
                               show_column_names = show_sample_names,
                               column_labels = names(mysamples),
                               column_names_rot = 0,
                               column_names_centered = T
                               )#,
# onco_base_default
savename <- file.path(output_dir,"fig_1b.tiff")
# pdf(file = savename,height=onco_height,width=onco_width)
tiff(file = savename,height=onco_height,width=onco_width, units="in",res=320)
# onco_base_default
# plot(1:10)
draw(onco_base_default, 
     heatmap_legend_side = "bottom",
     padding = unit(c(1, 1, 1, 1), "mm"))
# draw(onco_base_default)#, use_raster = T, raster_device="tiff")
# decorate_annotation("Mutations/MB", {
#   # grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
#   # grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(4, "mm"), just = "left")
#   grid.lines(x = c(0+0.5, ncol(oncomat.plot)+0.5), y=c(10,10),default.units = "native",
#              gp = gpar(lty=2, col="grey30",lwd=3))
# })
# 
dev.off()



