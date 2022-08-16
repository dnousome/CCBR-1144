rm(list=ls())


setwd("~/Documents/ccbr1144_gameiro/analysis/8")

library(maftools)
library(ggplot2)
# maf_file="/Users/tandonm/Documents/ccbr1144_gameiro/new_pipeline_out/somatic_paired/SNP_Indels/mutect2_out/all_somatic_variants.maf"
# filtered_maf_file="data/mutect2.filtered.maf"
# out_folder=file.path(getwd(),"plots","mutect2")

# maf_file="/Users/tandonm/Documents/ccbr1144_gameiro/new_pipeline_out/somatic_paired/SNP_Indels/merged_somatic_variants/all_somatic_variants.maf"
# filtered_maf_file="data/merged.filtered.maf"
# out_folder=file.path(getwd(),"plots","merged")

# maf_file="~/Documents/ccbr1144_gameiro/SF_delivery/maf/somatic_combined.maf.csv"
maf_file="~/Documents/ccbr1144_gameiro/SF_delivery/maf/dragen.all.maf"
filtered_maf_file="data/dragen.filtered.maf"
out_folder=file.path(getwd(),"plots","dragen_2")

if (! dir.exists(out_folder)){dir.create(out_folder, recursive = T)}
# mafobj <- read.maf(maf_file)
source("~/Documents/helper_functions/helper_functions.oncoplot.R")
if (! file.exists(filtered_maf_file)) {
  mafobj_full <- filter_maf_chunked(maf_file, chunk_lines=100000, 
                                    grep_vcf_filter_col="PASS|\\.|panel_of_normals",
                                    # grep_vcf_filter_col="PASS",
                     # n_callers=ncallers, variant_caller="consensus",
                     savename = filtered_maf_file,
                     non_silent_only=F,
                     # variant_caller="consensus",n_callers=1,
                     # t_alt_min=1, t_alt_max=1e12, t_depth_min=2,t_depth_max=1e12, 
                     t_alt_min=5, t_alt_max=1e12, t_depth_min=20,t_depth_max=1e12, 
                     tumor_freq_min=0.05,tumor_freq_max=1,
                     n_alt_min=0,n_alt_max=5,n_depth_min=0,n_depth_max=1e12,
                     norm_freq_min=0,norm_freq_max=0.01
  )
  
} else {
  mafobj_full <- read.maf(filtered_maf_file)
}

# mafobj_full@variants.per.sample
output_tbl <- make_variant_table(mafobj_full, extra_cols = c("Variant Callers"="set","Filters"="FILTER"))
library(openxlsx)
# write.xlsx(output_tbl, file=file.path(out_folder,"somatic_variants.consensus.xlsx"), overwrite = T)
write.xlsx(output_tbl, file=file.path(out_folder,"somatic_variants.dragen.xlsx"), overwrite = T)

genes_file="~/Documents/ccbr1144_gameiro/lab_data/selected_genes.xlsx"
pathwaysdf <- read.xlsx(genes_file, sheet="oncoplot")
custom_genes_df <- data.frame(Reason=stringr::str_wrap(gsub("_"," ",gsub("HALLMARK_","",pathwaysdf$Reason)), width=10),
                       Hugo_Symbol=pathwaysdf$Hugo_Symbol, 
                       stringsAsFactors = F)

# library(msigdbr)
# ## Get human gene sets from msigdb
# # m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
# m_t2g <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME") %>% 
#   dplyr::select(gs_name, gene_symbol, gs_subcat)
# ## Select just ones matching "interferon gamma"
# pathwaysdf <- data.frame(m_t2g)
# table(pathwaysdf$gs_subcat)
# mypaths <- pathwaysdf[grepl("INTERFERON_GAMMA|ANTIGEN|JAK",pathwaysdf$gs_name, ignore.case = T),]
# # mypaths <- pathwaysdf[grepl("JAK",pathwaysdf$gs_name, ignore.case = T),]
# mypaths$gs_name[grepl("INTERFERON_GAMMA", mypaths$gs_name)] <- "IFN-gamma"
# mypaths$gs_name[grepl("ANTIGEN",mypaths$gs_name)] <- "Antigen Presentation"
# mypaths$gs_name[grepl("JAK",mypaths$gs_name)] <- "JAK/STAT/IL12"
# db_genes_df <- data.frame(Reason=stringr::str_wrap(gsub("_"," ",gsub("REACTOME_","",mypaths$gs_name)), width=10),
#                        Hugo_Symbol=mypaths$gene_symbol, 
#                        stringsAsFactors = F)
# db_genes_df <- db_genes_df[!duplicated(db_genes_df$Hugo_Symbol),]

source("~/Documents/helper_functions/helper_functions.oncoplot.R")
targets_file="~/Documents/exome_bed_files/SureSelect_mm10.bed"
targets_coverage <- compute_exome_coverage(targets_file)/1e6



# WT_samples <- c("8_MC38_wt","1_EMT6_wt","21_Balb_c_tail","22_C57BL_6_tail")
# cell_line_samples1 <- c("20_TC1","14_TC1A9","12_CMT64","13_RVP3")
# cell_line_samples2 <- c("15_RM1","18_LLC","20_TC1","23_Kp1",
#                         "16_CT26","17_4T1","19_TS_A")
# cell_line_samples3 <- c("8_MC38_wt","1_EMT6_wt",
#                         "15_RM1","18_LLC","20_TC1","23_Kp1",
#                         "16_CT26","17_4T1","19_TS_A")
balbc_samples <- c("21_Balb_c_tail",
                   "16_CT26",
                   "1_EMT6_wt",
                   "17_4T1",
                   "19_TS_A"
                   )
b6_samples <- c("22_C57BL_6_tail",
                "8_MC38_wt",
                "12_CMT64",
                "13_RVP3",
                "15_RM1",
                "18_LLC",
                "20_TC1",
                "14_TC1A9",
                "23_Kp1",
                "21_Balb_c_tail"
                )
MC38_samples <- c("8_MC38_wt",
                  "11_MC38_B2m_3",
                  "9_MC38_Jak1_3",
                  "10_MC38_Jak1_11"
                  )

EMT6_samples <- c("1_EMT6_wt",
                  "2_EMT6_442_5_B2m_KO",
                  "3_EMT6_442_6_B2m_KO",
                  "4_EMT6_450_4_Jak1_KO",
                  "5_EMT6_450_6_Jak1_KO",
                  "6_EMT6_454_5_LMP2_KO",
                  "7_EMT6_454_6_LMP2_KO"
)

allgrps <- list(#WT=WT_samples,
                MC38_wt=MC38_samples,
                EMT6_wt=EMT6_samples,
                BALBC_tail=balbc_samples,
                C57B6_tail=b6_samples
                )
names(allgrps) <- unlist(lapply(allgrps, "[[", 1))

# rawmaf <- read.maf(maf_file)
# this <- subsetMaf(rawmaf, tsb=MC38_samples, genes = c("Jak1","B2m") )

if (!dir.exists(out_folder)){dir.create(out_folder, recursive = T)}

source("~/Documents/helper_functions/helper_functions.oncoplot.R")
# file.edit("~/Documents/helper_functions/helper_functions.oncoplot.R")
library(dplyr)

# custom_filtered_mafdata <- list()

# for (grpname in names(allgrps)) {
custom_filtered_mafdata <- lapply(names(allgrps), function(grpname){
  # grpname=names(allgrps)[1]
  mysamples <- allgrps[[grpname]]
  mafobj <- subsetMaf(mafobj_full, tsb=mysamples, dropLevels=F)
  # mafobj <- subsetMaf(mafobj_full, tsb=mysamples)
  # this <- subsetMaf(mafobj_full, tsb=mysamples, genes = c("Jak1","B2m") )
  mafdata <- rbind(mafobj@data, mafobj@maf.silent)
  # mafobj <- mafobj_full
  mafdata_by_sample <- split(mafdata, mafdata$Tumor_Sample_Barcode)
  id_cols <- c("Chromosome","Start_Position","End_Position","Hugo_Symbol","HGVSp_Short")

  mutations_list <- lapply(mafdata_by_sample, function(currmafdata, mycols) {
    id_string <- apply(currmafdata[,..mycols],1,paste0, collapse="_")
    return(id_string)
  }, id_cols)

  wt_muts=mutations_list[[mysamples[1]]]
  wt_removed <- do.call(rbind, lapply(names(mutations_list), function(listname, mut_list, exclude_muts, all_muts){
    # listname=names(mutations_list)[1]
    # mut_list=mutations_list
    # exclude_muts=wt_muts
    # all_muts=mafdata_by_sample
    mut_strs <- mut_list[[listname]]
    idx <- which(!mut_strs %in% exclude_muts)
    all_muts[[listname]][idx,]
  }, mutations_list, wt_muts, mafdata_by_sample)
  )
  return(wt_removed)
})

unfiltered_varcount <- mafobj_full@variants.per.sample
colnames(unfiltered_varcount) <- c("Tumor_Sample_Barcode","UnsubtractedVariants")

pairs_data <- unique(mafobj_full@data[,c("Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode")])
wt_subtract_df <- do.call(rbind,lapply(names(allgrps), function(x){data.frame(Tumor_Sample_Barcode=allgrps[[x]], WT_subtracted=x)}))
wt_subtract_df <- wt_subtract_df[!wt_subtract_df$Tumor_Sample_Barcode==wt_subtract_df$WT_subtracted,]
wt_subtract_df$WT_subtracted[! wt_subtract_df$WT_subtracted %in% unfiltered_varcount$Tumor_Sample_Barcode] <- NA
pairs_data <- dplyr::left_join(pairs_data, wt_subtract_df)

wt_filtered_maf <- read.maf(do.call(rbind, custom_filtered_mafdata))
burden_data <- wt_filtered_maf@variants.per.sample
burden_data$burden <- burden_data$Variants/targets_coverage

# unfiltered_varcount <- mafobj_full@variants.per.sample
burden_data <- dplyr::left_join(unfiltered_varcount, burden_data)

burden_data <- dplyr::left_join(burden_data, pairs_data)


library(openxlsx)
write.xlsx(burden_data, file = file.path(out_folder,paste0("all_burden_data.xlsx")), rowNames = F)

for (grpname in names(allgrps)) {
  # mafobj@variants.per.sample
  mysamples <- allgrps[[grpname]]
  mysamples <- mysamples
  mafobj <- subsetMaf(wt_filtered_maf, tsb=mysamples, dropLevels=F)
  burden_data <- mafobj@variants.per.sample
  
  curr_out=file.path(out_folder, grpname)
  if (!dir.exists(curr_out)){dir.create(curr_out, recursive = T)}
  source("~/Documents/helper_functions/helper_functions.oncoplot.R")
  custom_oncoplot_file1 <- file.path(curr_out,"customonco.pdf")
  make_oncoplot2(mafobj,cohort_freq_thresh = NULL, genes_to_plot = custom_genes_df, 
                 use_clinvar_anno = T, show_sample_names = T, total_mut = T,
                 custom_column_order = mysamples,
                 include_all = T, include_all_genes = T,
                 title_text = NULL,
                 savename=custom_oncoplot_file1)

  burden_data$burden <- burden_data$Variants/targets_coverage
  burden_data$Tumor_Sample_Barcode <- factor(burden_data$Tumor_Sample_Barcode, levels=mysamples)
  custom_burdenplot <- file.path(out_folder,grpname, paste0(grpname,".count.pdf"))
  xaxis_text <- element_text(angle=30, hjust=1)
  burdenplot <- ggplot(burden_data) + 
    geom_col(aes(x=Tumor_Sample_Barcode, y=Variants), fill="cornflowerblue") +
    theme_linedraw(base_size = 12) +
    xlab("") + ylab("Number of somatic mutations") +
    theme(
      axis.text.x = xaxis_text,
      axis.ticks.x = element_blank(),
      legend.position = "right",
      legend.text = element_text(size=8),
      legend.title = element_blank(),
      legend.key.height = unit(0.01,"npc"),
      legend.key.width =  unit(0.02,"npc"),
      plot.title = element_text(hjust = 0.5),
      panel.border = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank())
  ggsave(custom_burdenplot, width=5, height=4, plot = burdenplot)
  
  custom_burdenplot <- file.path(out_folder,grpname, paste0(grpname,".burden.pdf"))
  xaxis_text <- element_text(angle=30, hjust=1)
  burdenplot <- ggplot(burden_data) + 
    geom_col(aes(x=Tumor_Sample_Barcode, y=burden), fill="cornflowerblue") +
    theme_linedraw(base_size = 12) +
    xlab("") + ylab("Mutations per Mb") +
    theme(
      axis.text.x = xaxis_text,
      axis.ticks.x = element_blank(),
      legend.position = "right",
      legend.text = element_text(size=8),
      legend.title = element_blank(),
      legend.key.height = unit(0.01,"npc"),
      legend.key.width =  unit(0.02,"npc"),
      plot.title = element_text(hjust = 0.5),
      panel.border = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank())
  ggsave(custom_burdenplot, width=5, height=4, plot = burdenplot)
  write.table(burden_data, file = file.path(out_folder,grpname,paste0(grpname,".burden_data.txt")), quote = F, row.names = F, sep="\t")
}  

source("~/Documents/helper_functions/helper_functions.oncoplot.R")

grpname="EMT6"
mysamples <- allgrps[[grpname]]
mafobj <- subsetMaf(mafobj_full, tsb=mysamples, dropLevels=F)
burden_data <- mafobj@variants.per.sample
burden_data$burden <- burden_data2$Variants/targets_coverage

custom_burdenplot <- file.path(out_folder,grpname, paste0(grpname,".count.pdf"))
xaxis_text <- element_text(angle=30, hjust=1)
burdenplot <- ggplot(burden_data) + 
  geom_col(aes(x=Tumor_Sample_Barcode, y=Variants), fill="cornflowerblue") +
  theme_linedraw(base_size = 12) +
  xlab("") + ylab("Number of somatic mutations") +
  theme(
    axis.text.x = xaxis_text,
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size=8),
    legend.title = element_blank(),
    legend.key.height = unit(0.01,"npc"),
    legend.key.width =  unit(0.02,"npc"),
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank())
ggsave(custom_burdenplot, width=5, height=4, plot = burdenplot)

custom_burdenplot <- file.path(out_folder,grpname, paste0(grpname,".burden.pdf"))
xaxis_text <- element_text(angle=30, hjust=1)
burden_data2 <- burden_data
burdenplot <- ggplot(burden_data) + 
  geom_col(aes(x=Tumor_Sample_Barcode, y=burden), fill="cornflowerblue") +
  theme_linedraw(base_size = 12) +
  xlab("") + ylab("Mutations per Mb") +
  theme(
    axis.text.x = xaxis_text,
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size=8),
    legend.title = element_blank(),
    legend.key.height = unit(0.01,"npc"),
    legend.key.width =  unit(0.02,"npc"),
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank())
ggsave(custom_burdenplot, width=5, height=4, plot = burdenplot)

source("~/Documents/helper_functions/helper_functions.oncoplot.R")
mysamples <- EMT6_samples[!grepl("wt",EMT6_samples,ignore.case = T)]
mafobj <- subsetMaf(mafobj_full, tsb=mysamples)

custom_burdenplot1 <- file.path(out_folder,"customburden.count.pdf")
burdenplot <- make_burden_plot(mafobj,sample_order = mysamples,
                             plotType = "Barplot", add_median = T)
ggsave(custom_burdenplot1, width=5, height=4, plot = burdenplot)




source("~/Documents/helper_functions/helper_functions.oncoplot.R")
mysamples <- cell_line_samples1[grepl("wt",cell_line_samples1,ignore.case = T)]
mafobj <- subsetMaf(mafobj_full, tsb=mysamples)
custom_burdenplot2 <- file.path(out_folder,"customburden.count.2.pdf")
burdenplot <- make_burden_plot(mafobj_full,sample_order = mysamples,
                             plotType = "Barplot", add_median = T)
ggsave(custom_burdenplot2, width=5, height=4, plot = burdenplot)


library(MAFDash)
plotlist <- list("summary_plot"= T,
                 # "Silent vs. Non-Silent Mutations"= myplot,
                 # "TCGA Compare"= tcgacomp_image,
                 "Burden"= custom_burdenplot,
                 "Oncoplot - Select Genes" = custom_oncoplot_file1,
                 "Oncoplot - Select Pathways" = custom_oncoplot_file2
                 # "Selected Oncoplot" = custom_oncoplot_file,
                 # "Ti/Tv" = titv_plot_file
                 # "heatmap" = F,
                 # "Sample Signatures" = myplot2,
                 # "COSMIC Signatures" = sighm_image
                 # "Sample Signatures (with silent)" = myplot3,
                 # "COSMIC Signatures (with silent)" = sighm_image2
)


## Render dashboard
# dest_lib_path=file.path("~/Documents/my_tools/MAFdash/Rpackage/build/lib_build")
# library(MAFDashRPackage,lib.loc = dest_lib_path)
library(MAFDash)
html_filename=file.path("ccbr1144.summary.html")
getMAFDashboard(MAFfilePath = mafobj,
                plotList = plotlist,
                outputFileName = html_filename, 
                outputFileTitle = paste0("CCBR-1144: WES Murine Cell Lines")
)

file.copy("/private/var/folders/cf/n051mjpx001cwsqsrky6cvl1jtm0m6/T/RtmpIjs64R/ccbr1144.summary.html",html_filename, overwrite = T)



source("~/Documents/helper_functions/helper_functions.mutsig.R")
library(MutationalPatterns)
# maf.filtered <- read.maf(maf_file)
mut_mat_mouse <- make_tnm(mafobj, use_silent_mutations = T)
make_signature_plot(mut_mat_mouse, savepath=file.path("sbs.signatures.pdf"))

# source("scripts/helper_functions.mutsig.R")
cosmic_data <- load_cosmic_data()
cosmic_sigdata <- cosmic_data$signatures
cosmic_sigs <- as.matrix(cosmic_sigdata[match(rownames(mut_mat_mouse),cosmic_sigdata$Somatic.Mutation.Type),
                                        grepl("^SBS",colnames(cosmic_sigdata))])

fit_res_cosmic <- fit_to_signatures(mut_mat_mouse, cosmic_sigs)
contributions_raw_cosmic <- t(fit_res_cosmic$contribution)
contributions_cosmic <- contributions_raw_cosmic/rowSums(contributions_raw_cosmic)
contributions_cosmic <- contributions_cosmic[,match(rownames(cosmic_data$etio), colnames(contributions_cosmic))]
# contributions_cosmic <- contributions_cosmic[,colSums(contributions_cosmic) > 0.1]
contributions_cosmic <- contributions_cosmic[,colSums(contributions_cosmic) > 0]

sigplotdata <- reshape2::melt(contributions_cosmic)
colnames(sigplotdata) <- c("Sample","Signature","Relative_Contribution")
# sigplotdata$Etiology <- cosmic_data$etio$Etiology[match(sigplotdata$Signature, rownames(cosmic_data$etio))]
sigplotdata$Etiology <- paste0("[", cosmic_data$etio$Etiology_category[match(sigplotdata$Signature, rownames(cosmic_data$etio))], "] ",
                               cosmic_data$etio$Etiology[match(sigplotdata$Signature, rownames(cosmic_data$etio))])
sigplotdata$text <- paste0("Sample: ", sigplotdata$Sample,"\n",
                           "Signature: ", sigplotdata$Signature,"\n",
                           "Etiology: ", sigplotdata$Etiology,"\n")

sbs_colors <- cosmic_data$etio$annotation_color
names(sbs_colors) <- rownames(cosmic_data$etio)
sbs_colors <- sbs_colors[match(colnames(contributions_cosmic), names(sbs_colors), nomatch=0)]

etio_colors <- cosmic_data$etio$category_color
names(etio_colors) <- paste0("[", cosmic_data$etio$Etiology_category,"] ",cosmic_data$etio$Etiology)
etio_colors <- etio_colors[names(etio_colors) %in% sigplotdata$Etiology]

sig_contrib_plot <- ggplot(sigplotdata) + 
  geom_col(aes(x=Sample, y=Relative_Contribution, fill=Etiology, text=text), color="grey20", size=0.5) +
  scale_fill_manual(values=etio_colors) +
  xlab("") + ylab("Relative Contribution of Signature") +
  theme_linedraw()
ggsave("signature_contributions.pdf", width=6, height=4, plot = sig_contrib_plot)

plotly::ggplotly(sig_contrib_plot, tooltip=c("text"))













titv_plot_file <- file.path(getwd(),"titv.pdf")
pdf(titv_plot_file,width=8, height=8)
plotTiTv(titv(mafobj),showBarcodes = T)
dev.off()



library(msigdbr)
# library(stringr)
## Get human gene sets from msigdb
# m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol, gs_subcat)


## For oncoplot and other MAF plotting functions
# source("https://raw.githubusercontent.com/mtandon09/mt_helpers/main/helper_functions.oncoplot.R")


## Select just ones matching "interferon gamma"
pathwaysdf <- data.frame(m_t2g)
mypaths <- pathwaysdf[grepl("INTERFERON_GAMMA",pathwaysdf$gs_name, ignore.case = T),]

## Set up a data frame with columns "Reason" and "Hugo_Symbol"
## The "Reason" value is used to label the plot, so here I'm replacing _ with spaces and adding text wrapping with stringr
genes_df <- data.frame(Reason=stringr::str_wrap(gsub("_"," ",gsub("HALLMARK_","",mypaths$gs_name)), width=10),
                       Hugo_Symbol=mypaths$gene_symbol, 
                       stringsAsFactors = F)


source("~/Documents/helper_functions/helper_functions.oncoplot.R")
## Sizing can get tricky with larger plots; it will be adjusted automatically for the number of samples and genes if saving to pdf
make_oncoplot2(mafobj,cohort_freq_thresh = 0.7, genes_to_plot = genes_df, use_clinvar_anno = T, show_sample_names = T, savename="customonco.pdf")




