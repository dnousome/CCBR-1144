
filter_maf_tbl <- function(maftbl, # A tibble for filtering
                           flag_genes="default",  # Character vector of gene symbols to exclude; or set to 'default' to remove preset list of 50 genes from the FLAG paper
                           #save_name=NULL,
                           grep_vcf_filter_col="PASS|\\.",  # Regex to select from the FILTER column
                           non_silent_only=F,     # If TRUE, only non-silent mutations are returned (same definition as maftools)
                           t_alt_min=2,           # Minimum count for ALT allele in tumor sample
                           t_alt_max=1e12,        # Maximum count for ALT allele in tumor sample
                           t_depth_min=5,         # Minimum total depth (REF+ALT) in tumor sample
                           t_depth_max=1e12,      # Maximum total depth (REF+ALT) in tumor sample
                           tumor_freq_min=0.01,
                           tumor_freq_max=1,
                           n_alt_min=0,
                           n_alt_max=1,
                           n_depth_min=0,
                           n_depth_max=1e12,
                           norm_freq_min=0,
                           norm_freq_max=0.01,
                           gnomAD_AF_min=0,
                           gnomAD_AF_max=0.001,
                           AF_min=0,
                           AF_max=0.001,
                           ExAC_AF_min=0, 
                           ExAC_AF_max=0.001, 
                           n_callers=2,
                           variant_caller=NULL) {
  
  if (length(flag_genes)==0) {
    flag_genes <- c()
  } else if (flag_genes[1]=="default") {
    flag_genes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
  }
  # browser()
  require(tibble)
  require(dplyr)
  df <- as_tibble(maftbl)
  maf_df.raw <- df[df$Hugo_Symbol != "Hugo_Symbol",]
  
  mafgenome <- detect_maf_genome(maf_df.raw)[[1]]
  if (mafgenome %in% c("hg19","hg38")) {
    maf_df.raw <- maf_df.raw[!df$Hugo_Symbol %in% flag_genes,]
  }
  
  if ("FILTER" %in% colnames(maf_df.raw)) {
    if (!is.null(grep_vcf_filter_col)) {
      maf_df.raw <- maf_df.raw[grepl(grep_vcf_filter_col,maf_df.raw$FILTER),]
    }
  } else {
    message("FILTER column not found; skipping...")
  }
  
  #Variant Classification with High/Moderate variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                   "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                   "In_Frame_Ins", "Missense_Mutation")
  if ("Variant_Classification" %in% colnames(maf_df.raw)) {
    if (non_silent_only) {
      # maf_df.raw <- maf_df.raw[grepl(grep_vcf_filter_col,pull(maf_df.raw,FILTER)),]
      maf_df.raw <- maf_df.raw[maf_df.raw$Variant_Classification %in% vc.nonSilent,]
    }
  } else {
    message("Variant_Classification column not found; skipping...")
  }
  
  caller_set_column="set"
  filter_caller=rep(TRUE,nrow(maf_df.raw))
  if (! is.null(variant_caller)) {       ### Set 'variant_caller' to NULL to skip any filtering based on caller
    if (caller_set_column %in% colnames(maf_df.raw)) {
      maf_df.raw$set[maf_df.raw$set=="" & maf_df.raw$Hugo_Symbol=="Hugo_Symbol"] <- caller_set_column  ## This can happen if MAFs are lazily cat-ed together (mutiple header lines in the file)
      maf_df.raw$set[maf_df.raw$set==""] <- "N.A."
      if (variant_caller == "consensus") {   ### Set 'variant_caller' to 'consensus' to keep variants by two or more callers
        # filter_caller <- grepl("-|Intersection", maf_df.raw$set)
        filter_caller <- unlist(lapply(strsplit(maf_df.raw$set,"-"), function(x) {length(x)>=n_callers | "Intersection" %in% x}))
      } else {                             ### Set 'variant_caller' to one of the callers (mutect, mutect2, vardict, or strelka) to get only that caller
        # filter_caller <- grepl(paste0(variant_caller,"[|-]|Intersection"), maf_df.raw$set)
        filter_caller <- unlist(lapply(strsplit(maf_df.raw$set,"-"), function(x) {any(unique(c(variant_caller,"Intersection")) %in% x)}))
      }
    }
  }
  maf_df.raw <- maf_df.raw[filter_caller,]
  # browser()

  if (!"tumor_freq" %in% colnames(maf_df.raw)) {
    # if (! all(c("t_alt_count","t_depth") %in% colnames(maf_df.raw))) {
      # stop("Can't find t_alt_count or t_depth columns")
    # }
    if (! "t_alt_count" %in% colnames(maf_df.raw)) {
      maf_df.raw$t_alt_count <- NA
    }
    if (! "t_depth" %in% colnames(maf_df.raw)) {
      maf_df.raw$t_depth <- NA
    }
    maf_df.raw$tumor_freq <- as.numeric(maf_df.raw$t_alt_count)/as.numeric(maf_df.raw$t_depth)
  }
  if (!"norm_freq" %in% colnames(maf_df.raw)) {
    # if (! all(c("n_alt_count","n_depth") %in% colnames(maf_df.raw))) {
    #   maf_df.raw$norm_freq <- rep(0, nrow(maf_df.raw))
    # } else {
    #   maf_df.raw$norm_freq <- as.numeric(maf_df.raw$n_alt_count)/as.numeric(maf_df.raw$n_depth)
    # }
    
    if (! "n_alt_count" %in% colnames(maf_df.raw) ) {
      maf_df.raw$n_alt_count <- NA
    }
    if (! "n_depth" %in% colnames(maf_df.raw)) {
      maf_df.raw$n_depth <- NA
    }
    maf_df.raw$norm_freq <- as.numeric(maf_df.raw$n_alt_count)/as.numeric(maf_df.raw$n_depth)
  }
  
  
  maf_num_filter_columns <- list("t_alt_count"=c(min=t_alt_min, max=t_alt_max),
                                 "t_depth"=c(min=t_depth_min, max=t_depth_max),
                                 "tumor_freq"=c(min=tumor_freq_min, max=tumor_freq_max),
                                 "n_alt_count"=c(min=n_alt_min, max=n_alt_max),
                                 "n_depth"=c(min=n_depth_min, max=n_depth_max),
                                 "norm_freq"=c(min=norm_freq_min, max=norm_freq_max),
                                 "gnomAD_AF"=c(min=gnomAD_AF_min, max=gnomAD_AF_max),
                                 "AF"=c(min=AF_min, max=AF_max),
                                 "ExAC_AF"=c(min=ExAC_AF_min, max=ExAC_AF_max)
  )
  # browser()
  numfilter_columns <- names(maf_num_filter_columns)[names(maf_num_filter_columns) %in% colnames(maf_df.raw)]
  notfound <- setdiff(names(maf_num_filter_columns), numfilter_columns)
  if (length(notfound) > 0 ) {
    message(paste0("Couldn't find these columns; skipping filtering for these: ", paste0(notfound, collapse=",")))
  }
  # maf_df.raw$Matched_Norm_Sample_Barcode[is.na(maf_df.raw$Matched_Norm_Sample_Barcode)] <- 0
  if (all(maf_df.raw$Tumor_Sample_Barcode==maf_df.raw$Matched_Norm_Sample_Barcode, na.rm=T)) {
    warning("Normal IDs match tumor IDs, skipping n_alt_count, n_depth, and norm_freq filters...")
    numfilter_columns <- setdiff(numfilter_columns, c("n_alt_count","n_depth","norm_freq"))
  }
  
  return_df <- maf_df.raw
  if (length(numfilter_columns)>0) {
    all_num_filters <- lapply(numfilter_columns, function(col_name) {
      currdata <- as.numeric(pull(maf_df.raw,col_name))
      currdata[is.na(currdata)] <- 0  ## Keeps NAs
      filter_vec <- currdata >= maf_num_filter_columns[[col_name]]["min"] & currdata <= maf_num_filter_columns[[col_name]]["max"]
      return(filter_vec)
    })
    names(all_num_filters) <- numfilter_columns
    # browser()
    # print(lapply(all_num_filters, sum))
    final_num_filters <- Reduce("&", all_num_filters)
    
    return_df <- maf_df.raw[final_num_filters,]
  }
  return(return_df)
}



filter_maf_chunked <- function(maf, chunk_lines=10000, savename=NULL,...) {
  require(maftools)
  clindata <- NULL
  if ("MAF" %in% class(maf)) {
    filtered_df <- filter_maf_tbl(
                            rbind(maf@data, maf@maf.silent),
                            ...)
    clindata <- maf@clinical.data
    
  } else if (file.exists(maf)) {
    require(readr)
    readr_filterfunc <- function(df, pos) {
      filter_maf_tbl(df,...)
    }
    filtered_df <- read_tsv_chunked(maf,chunk_size = chunk_lines, col_types = cols(), callback = DataFrameCallback$new(readr_filterfunc), comment="#")
  } else {
    # stop(paste0("Don't know what to do with input type '",class(maf),"'"))
    stop(paste0("MAF file not found: '",maf,"'"))
  }
  
  if ( ! is.null(savename) ) {
    if (!dir.exists(dirname(savename))) {dir.create(dirname(savename), recursive = T)}
    write.table(filtered_df, file=savename,sep="\t", quote=F, row.names=F)
  }
  maf.filtered <- read.maf(filtered_df, clinicalData = clindata)
  invisible(maf.filtered)
}



make_burden_plot <- function(maf.filtered, plotType=NULL, mb_covered=NULL, 
                             sample_order=NULL,
                             save_data_to_file=NULL, 
                             add_median=T) {
  
  require(dplyr)
  num_var_data <- maf.filtered@variants.per.sample
  colnames(num_var_data) <- c("Tumor_Sample_Barcode","Variants_filtered")
  num_var_data$mut_burden_count <- num_var_data$Variants_filtered
  num_var_data$mut_burden <- num_var_data$mut_burden_count
  y_label_text="Mutation Count"
  if (is.numeric(mb_covered)) {
    print("Normalizing mutation count by covered bases...")
    num_var_data$mut_burden <- num_var_data$mut_burden/mb_covered
    y_label_text="Mutation Burden (mutations/Mb)"
  }
  
  nsamples=nrow(num_var_data)
  if (is.null(plotType)) {
    plotType <- ifelse(nrow(num_var_data) > 15, "Dotplot", "Barplot")
    print(paste0("Using plot type: ", plotType))
  }
  # browser()
  ## Re-jigger the factor levels so they're ordered by decreasing mutation burden (for the plotting)
  num_var_data$Tumor_Sample_Barcode <- factor(num_var_data$Tumor_Sample_Barcode,
                                              levels=num_var_data$Tumor_Sample_Barcode[order(num_var_data$mut_burden, decreasing = T)])
  
  ########################################################
  #### 5. Generate plots for mutation burden
  
  ## Pick colors
  median_mut_burdens <- num_var_data %>% summarise(median=median(mut_burden))
  
  num_var_data$xlabel <- factor(num_var_data$xlabel,
                                levels=num_var_data$xlabel[order(num_var_data$mut_burden, decreasing = T)])
  num_var_data$hoverlabel <- paste0("Sample: ",num_var_data$Tumor_Sample_Barcode,"\nMutations: ", num_var_data$mut_burden)
  
  ### Mutation burden stacked with variant classification counts
  ### Works better for smaller cohorts, say < 20
  variant_type_per_sample <- as.data.frame(maf.filtered@variant.classification.summary)
  var_type.melt <- reshape2::melt(variant_type_per_sample, id.vars="Tumor_Sample_Barcode",variable.name="Classification",value.name="mutation_count")
  var_type.melt$mut_burden <- var_type.melt$mutation_count
  if (is.numeric(mb_covered)) {
    var_type.melt$mut_burden <- var_type.melt$mut_burden/mb_covered
  }
  
  
  plotdata <- var_type.melt[var_type.melt$Classification != "total",]
  tsb_order <- variant_type_per_sample$Tumor_Sample_Barcode[order(variant_type_per_sample$total, decreasing = T)]
  if (!is.null(sample_order)) {
    tsb_order <- sample_order
  }
  plotdata$Tumor_Sample_Barcode <- factor(as.character(plotdata$Tumor_Sample_Barcode),
                                          levels=tsb_order)
  plotdata$Classification <- gsub("_"," ",plotdata$Classification)
  
  class_means <- plotdata %>% group_by(Classification) %>% summarise(mean=mean(mut_burden))
  plotdata$Classification <- factor(as.character(plotdata$Classification),
                                    levels=class_means$Classification[order(class_means$mean, decreasing = F)])
  
  my_class_colors <- mutation_colors[names(mutation_colors) %in% levels(plotdata$Classification)]
  
  plotdata$hoverlabel <- paste0("Sample: ",plotdata$Tumor_Sample_Barcode,"\nMutations: ", plotdata$mut_burden)
  
  if (plotType=="Barplot") {
    if (length(unique(plotdata$Tumor_Sample_Barcode)) <= 20) {
      xaxis_text <- element_text(angle=30, hjust=1)
    } else {
      xaxis_text <- element_blank()
    }
    burden_plot <- ggplot(plotdata, aes(x=Tumor_Sample_Barcode, y=mut_burden, text=hoverlabel)) +
      geom_bar(aes(fill=Classification), stat="identity",width=1,size=0.3, color="black") +
      scale_fill_manual(values=my_class_colors) +
      theme_linedraw(base_size = 12) +
      xlab("") + ylab(y_label_text) +
      # geom_hline(data = median_mut_burdens, aes(yintercept=median),linetype="dashed", color="grey60") +  ### Is screwed up with ggplotly
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
  } else {
    
    require(ggbeeswarm)
    ### Mutation Burden - Scatter/Dot plot
    ### Works better for larger cohorts
    alpha_val=1
    point_cex=2
    if (nrow(num_var_data) > 200) {
      alpha_val=0.5
    } else if (nrow(num_var_data) > 20) {
      alpha_val=0.8
    }
    burden_plot <- ggplot(num_var_data, aes(x=1, y=mut_burden, text=hoverlabel)) +
      # geom_beeswarm(color="blue3",cex=2,size=5,dodge.width=0.2,priority='density', alpha=alpha_val) +
      geom_quasirandom(color="blue3",width=0.3,size=5,alpha=alpha_val, method="quasirandom", bandwidth=0.1) +
      # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
      #              geom = "crossbar", width = 0.7, color="gray70", size = 0.2) +
      scale_y_log10()+
      theme_linedraw(base_size = 12) +
      ggtitle("Mutation Burden") +
      ylab(y_label_text) + xlab("") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  median_val=median(var_type.melt[var_type.melt$Classification== "total","mut_burden"])
  if (is.numeric(add_median)) {
    median_val=add_median
    add_median=T
  }
  
  if (add_median) {
    median_mut_burdens <- data.frame(median=median_val)
    burden_plot <- burden_plot + geom_hline(data = median_mut_burdens, aes(yintercept=median),linetype="dashed", color="grey60")
  }
  
  
  ## Write data to a file for external plotting if desired
  if (!is.null(save_data_to_file)) {
    if (dir.exists(dirname(save_data_to_file))) {
      outdata <- as.data.frame(maf.filtered@variant.classification.summary)
      outdata$total_per_mb <- outdata$total/mb_covered
      outdata$mb_covered <- mb_covered
      print(paste0("Saving plot data to ", save_data_to_file))
      write.table(outdata, file = save_data_to_file, sep="\t", quote=F,row.names = F)
    } else {
      warning("Path for output data file not found! Skipping...")
    }
  }
  
  
  return(burden_plot)
  
}



detect_maf_genome <- function(maf) {
  if ("data.frame" %in% class(maf)) {
    maf_tbl <- maf
  } else if ("MAF" %in% class(maf)) {
    maf_tbl <- maf@data
  } else {
    stop("Argument 'maf' must be a data.frame, data.table or MAF object")
  }
  if (!"NCBI_Build" %in% colnames(maf_tbl)) {
    warning("No genome information in MAF obj.")
    return(NA)
  }
  my_genome = unique(maf_tbl$NCBI_Build)
  if (length(my_genome) > 1) {
    warning("Multiple genomes listed in MAF obj. Trying the first one")
    my_genome <- my_genome[1]
  }
  return_genome <- switch(my_genome, GRCh38 = "hg38", GRCh37 = "hg19", 
                          GRCm38 = "mm10", NA)
  my_chrs <- unique(maf_tbl$Chromosome)
  add_chr = sum(grepl("^chr", my_chrs)) < length(my_chrs)
  pkg_prefix = ifelse(return_genome == "mm10", "BSgenome.Mmusculus.UCSC.", 
                      "BSgenome.Hsapiens.UCSC.")
  genome_package = paste0(pkg_prefix, return_genome)
  return(list(genome = return_genome, add_chr = add_chr, bsgenome_pkg = genome_package))
}



compute_exome_coverage <- function(targets_bed_file, out_file=NULL) {
  ##### This function will read the target regions BED file and
  #####  compute the sum of the lengths of the regions
  require(GenomicRanges)
  
  ## This bit will only read in the first three columns
  num_fields <- max(count.fields(targets_bed_file, sep = "\t"))
  my_classes <- c("character","integer","integer", rep("NULL", num_fields-3))
  
  ## Read the BED file as a table
  bed_data <- read.table(targets_bed_file, sep="\t",colClasses = my_classes, 
                         stringsAsFactors = F)
  colnames(bed_data) <- c("chr","start","end")
  
  ## Convert to a GenomicRanges object
  bed.gr <- makeGRangesFromDataFrame(bed_data)
  
  ## Collapse any overlapping features
  bed.gr.collapse <- reduce(bed.gr)
  
  ## Sum up the widths of each range
  total_exome_coverage = sum(width(bed.gr.collapse))
  
  if (! is.null(out_file)) {
    ## Write to a file
    write.table(total_exome_coverage,file = out_file, col.names =F, row.names = F)
  } else {
    return(total_exome_coverage)
  }
}


### Define colors for mutation types
mutation_colors <- c(Nonsense_Mutation="#ad7aff",Missense_Mutation="#377EB8",Frame_Shift_Del="#4DAF4A",
                     In_Frame_Ins="#ff008c",Splice_Site="#FF7F00",Multi_Hit="#FFFF33",Frame_Shift_Ins="#A65628",
                     In_Frame_Del="#f781bf",Translation_Start_Site="#400085",Nonstop_Mutation="#b68dfc",
                     Amp="green2",Del="darkred",
                     no_variants="#d6d6d6", Pathogenic="black",VUS="grey50")
names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
### List defining functions for color and shape of cells in oncoplot
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # "0" = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
  #             gp = gpar(fill = "#CCCCCC", col = NA))
  # },
  "Nonsense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Nonsense Mutation"], col = NA))
  },
  "Missense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Missense Mutation"], col = NA))
  },
  "Frame Shift Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Frame Shift Del"], col = NA))
  },
  "In Frame Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["In Frame Ins"], col = NA))
  },
  "Splice Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Splice Site"], col = NA))
  },
  "Multi Hit" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Multi Hit"], col = NA))
  },
  "Frame Shift Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Frame Shift Ins"], col = NA))
  },
  "In Frame Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["In Frame Del"], col = NA))
  },
  "Nonstop Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Nonstop Mutation"], col = NA))
  },
  "Translation Start Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Translation Start Site"], col = NA))
  },
  "Amp" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.25, "mm"), h-unit(0.25, "mm"),
              gp = gpar(fill = mutation_colors["Amp"], col = NA))
  },
  "Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.25, "mm"), h-unit(0.25, "mm"),
              gp = gpar(fill = mutation_colors["Del"], col = NA))
  },
  "no variants" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              # gp = gpar(fill = "#e0e0e0", col = NA))
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  "Pathogenic" = function(x, y, w, h) {
    # grid.points(x, y, pch = 18, size=w, gp=gpar(col=col["pathogenic"]))
    # grid.rect(x, y, w*0.7, h*0.2,
    #           gp = gpar(fill = col["pathogenic"], col = NA))
    # grid.rect(x, y, w*0.1, h*0.7,
    #           gp = gpar(fill = col["pathogenic"], col = NA))
    grid.rect(x, y, w*0.8, h*0.8,
              gp = gpar(col = mutation_colors["Pathogenic"], fill = NA, lwd=5))
  },
  "VUS" = function(x, y, w, h) {
    # grid.points(x, y, pch = 3, size=w,gp=gpar(col=col["VUS"], lwd=3))
    # grid.rect(x, y, w*0.2, h-unit(0.5, "mm"),
    #           gp = gpar(fill = col["VUS"], col = NA))
    grid.rect(x, y, w*0.8, h*0.8,
              gp = gpar(col = mutation_colors["VUS"], fill = NA, lwd=5))
  }
)




### Cretaes matrix for oncoplot from maf file
### Adapted from maftools: https://github.com/PoisonAlien/maftools/blob/master/R/oncomatrix.R
createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE){
  require(maftools)
  require(dplyr)
  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }
  
  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)
  
  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      vc = c("")
      names(vc) = 0
      
      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }
  
  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = unique(g))]
  }
  
  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x){
                                x = unique(as.character(x))
                                xad = x[x %in% c('Amp', 'Del')]
                                xvc = x[!x %in% c('Amp', 'Del')]
                                
                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                  # xvc = paste0(xvc, collapse="|")
                                }
                                
                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)
  
  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])
  
  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)
  
  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)
  
  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }
  
  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
  
  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)
  
  
  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId
    
    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]
    
    mdf = mdf[, -ncol(mdf)]
    
    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix
    
    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort
    
    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy
    
    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}

