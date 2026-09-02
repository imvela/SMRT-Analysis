################################################################################
# Prepare the workspace
rm(list = ls()); graphics.off(); cat("\014")

################################################################################
# User editable variables

# filenames
name_data_sample = "20260831-SMRT005C01-P2.csv"   # Experimental group data
name_plate_sample = "SMRT005C01_P2_plate_map.csv" # Data plate layout

# experimental parameters
exp_code = "SMRT005C01"
plate_number = 2
experimental_date = "20260831" #Date experiment was run on LCMS
VIAL = TRUE

# dilutions
CPDil = 1 # calibration plate dilution
SPDil = 1 # sample plate dilution

# Do not edit below
################################################################################
# Part 0 - Libraries and parameters
################################################################################

###
# load libraries
library(ggplot2); library(tidyr); library(tidyverse)
library(dplyr); library(ggforce); 
library(scales); library(stringr);
library(ggpubr);library(DescTools);
library(gridExtra);library(rlang);
library(htmltools);
#library(htmlwidgets);

###
# Define purity levels
ILApurity = 0.99 # purity of ILA used in calibration
IAldpurity = 0.97 # purity of IAld used in calibration
###

# Create max measurable values, accounting for dilutions
IPA_max = (SPDil/CPDil)*100
ILA_max = (SPDil/CPDil)*100
IAld_max = (SPDil/CPDil)*100
IAcr_max = (SPDil/CPDil)*75
Trp_max = (SPDil/CPDil)*550
# Set directory
base_dir = Sys.glob("~/Library/CloudStorage/GoogleDrive-*@kingdomsuperculture.com/Shared drives/Ingredient-Development/1. Active Projects/SMRT/Notebook/SMRT-QQQ")
setwd(base_dir)
###
# Define date of data analysis and folder for results
analysis_date <- as.character(format(Sys.Date(), "%Y%m%d"))
folder <- paste0("Indole-ProcessedData/",analysis_date,"-",exp_code,"-P", plate_number,"-Data",experimental_date,"/")
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
name_data_sample = paste0("Indole-QQQdata/", name_data_sample)   # Experimental group data
name_plate_sample = paste0("Indole-platemaps/", name_plate_sample) # Data plate layout

################################################################################
# Part 1 - Pre-processed data and files
################################################################################

###
# read plate layout
platemap <- read.csv(name_plate_sample, stringsAsFactors = TRUE)

###
# read data
DataQQQ <- read.csv(name_data_sample, stringsAsFactors = TRUE)

# rename variables
names(DataQQQ) <- c("Sample","DataFile","SampleType","Level", "day-time", "well_96",
                    "Trp.area", "Trp.conc","Trp.ret","Trp.SN","Trp.trpd-d5",
                    "ILA.area", "ILA.conc","ILA.ret","ILA.SN","ILA.IAAd7",
                    "IAld.area", "IAld.conc","IAld.ret","IAld.SN","IAld.IAAd7",
                    "IAcr.area", "IAcr.conc","IAcr.ret","IAcr.SN","IAcr.IAAd7",
                    "IPA.area", "IPA.conc","IPA.ret","IPA.SN","IPA.IAAd7",
                    )

# reformat the data
DataQQQ <- DataQQQ[c(2:length(DataQQQ$Sample)),]
DataQQQ$well_96 <- str_extract(DataQQQ$well_96, '\\b\\w+$')
DataQQQ <- DataQQQ[!is.na(DataQQQ$well_96),]

DataQQQ$Sequence <- gsub("^[^_S]*_S", "", DataQQQ$DataFile)
DataQQQ$Sequence<-substr(DataQQQ$Sequence,nchar(DataQQQ$Sequence)-4,nchar(DataQQQ$Sequence))
DataQQQ$Sequence <- as.numeric(gsub("*.d", "", DataQQQ$Sequence))

DataQQQ[,c("Trp.area", "Trp.conc","Trp.ret","Trp.SN","Trp.trpd-d5",
           "ILA.area", "ILA.conc","ILA.ret","ILA.SN","ILA.IAAd7",
           "IAld.area", "IAld.conc","IAld.ret","IAld.SN","IAld.IAAd7",
           "IAcr.area", "IAcr.conc","IAcr.ret","IAcr.SN","IAcr.IAAd7",
           "IPA.area", "IPA.conc","IPA.ret","IPA.SN","IPA.IAAd7")] <- lapply(DataQQQ[,c("Trp.area", "Trp.conc","Trp.ret","Trp.SN","Trp.trpd-d5",
                                                                                        "ILA.area", "ILA.conc","ILA.ret","ILA.SN","ILA.IAAd7",
                                                                                        "IAld.area", "IAld.conc","IAld.ret","IAld.SN","IAld.IAAd7",
                                                                                        "IAcr.area", "IAcr.conc","IAcr.ret","IAcr.SN","IAcr.IAAd7",
                                                                                        "IPA.area", "IPA.conc","IPA.ret","IPA.SN","IPA.IAAd7")], function(x) as.numeric(as.character(x)))

#DataQQQ <- DataQQQ[grep("Q2", DataQQQ$DataFile),]
#DataQQQ <- DataQQQ[DataQQQ$Sample >=14,]
#DataQQQ <- DataQQQ[DataQQQ$Sample <=300,]

# separate by sample type for processing
DataQQQ_ExtStd <- DataQQQ[DataQQQ$Sample == "241025_R8",]
#DataQQQ <- DataQQQ[DataQQQ$Sample != "241025_R8",]
DataQQQ_QC <- DataQQQ[DataQQQ$SampleType == "QC",]
DataQQQ_Blank <- DataQQQ[DataQQQ$SampleType == "Blank",]
DataQQQ_Cal <- DataQQQ[DataQQQ$SampleType == "Cal",]
DataQQQ <- DataQQQ[DataQQQ$SampleType == "Sample",]
## Only change if bad injection occurs
# DataQQQ<-DataQQQ[DataQQQ$ILA.IAAd7>4000,]
if(nchar(DataQQQ$well_96[1])>3){
  DataQQQ$well_96<-substr(DataQQQ$well_96, 3, nchar(DataQQQ$well_96))
}
DataQQQ <- merge(platemap,DataQQQ, by = "well_96")
# DataQQQ<-DataQQQ[DataQQQ$sample_id!="CIHP206X_D4_3",]
rm(platemap, name_data_sample, name_plate_sample)

################################################################################
# Part 2 - Internal and External QC
################################################################################

###
# essential conversion for plotting heat map
convertL <- function(ALP){
  NN <- which(LETTERS == ALP)
  return(NN)
}

### Internal control assessment
DataQQQ$row_96 <- as.numeric(lapply(DataQQQ$row_96, convertL))
QC_Int = as.character("Not evaluated")
if(!is_empty(DataQQQ)){
  graph_id = "Internal control assessment"
  deg <- 3;  cut_sd <- 2 # polynomial degree & residual cutoff
  
  # First fit → flag outliers, add InternalQC column ------------------- ##
  m1  <- lm(ILA.IAAd7 ~ poly(Sequence, deg, raw = TRUE), DataQQQ)
  sig <- sd(resid(m1))
  DataQQQ <- DataQQQ %>%
    mutate(InternalQC = if_else(abs(resid(m1)) > cut_sd * sig, "fail", "pass"))
  
  # Re-fit on “pass” points ------------------------------------------- ##
  m2 <- lm(ILA.IAAd7 ~ poly(Sequence, deg, raw = TRUE),
           data = filter(DataQQQ, InternalQC == "pass"))
  
  # Prediction grid for smooth curve ---------------------------------- ##
  pred <- tibble(Sequence = seq(min(DataQQQ$Sequence, na.rm = TRUE),
                                max(DataQQQ$Sequence, na.rm = TRUE), length.out = 200)) |>
    mutate(fit = predict(m2, newdata = cur_data()))
  
  # Plot --------------------------------------------------------------- ##
  p<- ggplot(DataQQQ, aes(Sequence, ILA.IAAd7)) +
    geom_point(aes(colour = InternalQC == "fail")) +
    scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "black"), guide = "none") +
    geom_line(data = pred, aes(y = fit),
              colour = "grey40", linetype = "dashed", linewidth = 0.8) +
    labs(title = graph_id,
         x = "Sequence", y = "IAAd7 area") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_blank(),
          axis.line        = element_line())
  
  save_name <- paste0(folder,"1-IndoleQQQ-InternalQC.png")
  ggsave(plot = p, width = 7, height = 5, dpi = 400, filename = save_name)
  
  fail_points <- DataQQQ$DataFile[DataQQQ$InternalQC == "fail"]
  Nfail = length(fail_points)
  Ffail = length(fail_points)/length(DataQQQ$DataFile)*100
  QC_Int = paste0(Nfail," samples failed, corresponding to ",Ffail, "% of samples")
}

### External control assessment
QC_Ext = as.character("Not evaluated")
if(length(DataQQQ_ExtStd$Sample)>0){
  ## Adjust for dilutions
  
  DataQQQ_ExtStd$ILA.conc <- DataQQQ_ExtStd$ILA.conc*30/(CPDil*ILApurity)
  DataQQQ_ExtStd$IAld.conc <- DataQQQ_ExtStd$IAld.conc*30/(CPDil*IAldpurity)
  DataQQQ_ExtStd$IPA.conc <- DataQQQ_ExtStd$IPA.conc*30/(CPDil)
  DataQQQ_ExtStd$Trp.conc <- DataQQQ_ExtStd$Trp.conc*30/(CPDil)
  DataQQQ_ExtStd$IAcr.conc <- DataQQQ_ExtStd$IAcr.conc*30/(CPDil)
  
  
  DataQQQ_ExtStd <- group_by(DataQQQ_ExtStd, Sample)
  DataQQQ_ExtStd_sum <- summarise(DataQQQ_ExtStd, 
                                  meanILA = mean(ILA.conc), sdILA= sd(ILA.conc),
                                  meanIAld = mean(IAld.conc), sdIAld= sd(IAld.conc))
  
  
  DataQQQ_ExtStd_sum <-unique(DataQQQ_ExtStd_sum)
  
  DataQQQ_ExtStd_sum <- as.data.frame(lapply(DataQQQ_ExtStd_sum, function(x) {
    if(is.numeric(x)) round(x, 1) else x
  }))
  
  # Specify the width and height
  fixed_width <- 6  # Adjust the width as needed
  fixed_height <- nrow(DataQQQ_ExtStd_sum) * 0.45  # Adjust the height based on the number of rows
  
  # Save the dataframe as a PNG image
  savename <- paste0(folder, "2-IndoleQQQ-ExternalQC.png")
  png(savename, width = fixed_width, height = fixed_height, units = "in", res = 450)
  grid.table(DataQQQ_ExtStd_sum)
  dev.off()
  DataQQQ_ExtStd_sum$date<-experimental_date
  DataQQQ_ExtStd_sum$experiment<-exp_code
  DataQQQ_ExtStd_sum$Sample<-"241025_R8"
  Data_QC_all<-as.data.frame(read.csv("Indole-QCdata/Experiment_indole_standard_SMRT.csv"))
  ## check if data has already been added to data frame and if not add to data frame and rewrite to file
  if(!(DataQQQ_ExtStd_sum$date%in%Data_QC_all$date &&DataQQQ_ExtStd_sum$experiment%in%Data_QC_all$experiment) ){
    # DataQQQ_ExtStd_sum$X<-length(Data_QC_all$X)+1
    Data_QC_all<-Data_QC_all[,c("Sample","meanILA","sdILA","meanIAld","sdIAld","date","experiment")]
    Data_QC_all<-rbind(Data_QC_all,DataQQQ_ExtStd_sum)
    write.csv(Data_QC_all,"Indole-QCdata/Experiment_indole_standard_SMRT.csv")
  }
  
  # ── MAIN WRAPPER ──────────────────────────────────────────────────────────────
  analyse_indole_experiment <- function(before   = NA,      # e.g. "20250201"
                                        after    = NA,      # e.g. "20250601"
                                        thr_ila  = 0.25,    # ILA threshold
                                        thr_iald = 0.25,    # IAld threshold
                                        today    = Sys.Date(),
                                        out_dir  = ".") {
    
    file <- "Indole-QCdata/Experiment_indole_standard_SMRT.csv"
    
    # ── load & date-filter rows ────────────────────────────────────────────────
    d <- read.csv(file) %>%
      mutate(date      = as.Date(as.character(date), "%Y%m%d"),
             keep_date = (is.na(before) | date >= as.Date(before, "%Y%m%d")) &
               (is.na(after)  | date <= as.Date(after , "%Y%m%d")))
    
    m0 <- d %>% filter(keep_date) %>%                      # first-pass means
      summarise(ila  = mean(meanILA , na.rm = TRUE),
                iald = mean(meanIAld, na.rm = TRUE))
    
    d <- d %>%
      mutate(outlier = keep_date &
               (abs(meanILA  - m0$ila )/m0$ila  > thr_ila  |
                  abs(meanIAld - m0$iald)/m0$iald > thr_iald),
             excl     = !keep_date | outlier,
             date_num = as.numeric(date))                  # numeric for LM
    
    # ── plotting & stats helper ────────────────────────────────────────────────
    plot_series <- function(y, ttl, file_name) {
      fit      <- lm(reformulate("date_num", y), data = d, subset = !excl)
      pred_val <- round(predict(fit,
                                tibble(date_num = as.numeric(today))), 2)
      avg_val  <- round(mean(d[[y]][!d$excl], na.rm = TRUE), 2)
      
      g <- ggplot(d, aes(date, .data[[y]], colour = excl)) +
        geom_point() +
        geom_smooth(data = d %>% filter(!excl),
                    method = "lm", se = FALSE, colour = "blue") +
        geom_hline(yintercept = avg_val, linetype = "dashed") +
        scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "red"),
                            guide = "none") +
        scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +
        labs(title = ttl,
             subtitle = sprintf("Predicted on %s: %.2f   |   Mean: %.2f",
                                format(today, "%Y-%m-%d"),
                                pred_val, avg_val)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
      
      ggsave(file.path(out_dir, file_name), g,
             width = 12, height = 5, dpi = 400)
      list(predicted = pred_val, average = avg_val)
    }
    
    # ── run & return stats ─────────────────────────────────────────────────────
    ila_stats  <- plot_series("meanILA",  "ILA experiment-to-experiment",
                              "3-IndoleQQQ-ILA-Exp2Exp.png")
    iald_stats <- plot_series("meanIAld", "IAld experiment-to-experiment",
                              "4-IndoleQQQ-IAld-Exp2Exp.png")
    
    invisible(list(ila_stats = ila_stats, iald_stats = iald_stats))
  }
  
  ## EXAMPLE CALL
  res <- analyse_indole_experiment(before  = "20250201",
                                   after   = "20260320",
                                   thr_ila = 0.25,
                                   thr_iald= 0.25,
                                   today = as.Date(as.character(DataQQQ_ExtStd_sum$date), "%Y%m%d"),
                                   out_dir = folder)
  Stat1 = round(DataQQQ_ExtStd_sum$meanILA[1]/res$ila_stats$predicted*100,0)
  Stat2 = round(DataQQQ_ExtStd_sum$meanILA[1]/res$ila_stats$average*100,0)
  Stat3 = round(DataQQQ_ExtStd_sum$meanIAld[1]/res$iald_stats$predicted*100,0)
  Stat4 = round(DataQQQ_ExtStd_sum$meanIAld[1]/res$iald_stats$average*100,0)
  
  QC_Ext_ILA = as.character(paste0("To predicted ILA value: ", Stat1,"%; To expected ILA value: ", Stat2,"%"))
  QC_Ext_IAld = as.character(paste0("To predicted IAld value: ", Stat3,"%; To expected IAld value: ", Stat4,"%"))
  QC_Ext_result <- ifelse(abs(Stat2-100) < 25 & abs(Stat4-100) < 25, "External QC values are within expected range", "External QC values are NOT within expected range")
  
  
  ### QC control assessment
  if(!is_empty(DataQQQ_ExtStd)){
    # Graph view
    graph_id = "Reference standard assessment"
    p1 <- ggplot(data = DataQQQ_QC, aes(Sequence, ILA.area)) + 
      geom_point() +
      labs(title = paste0(as.character(round(100*(max(DataQQQ_QC$ILA.area)-min(DataQQQ_QC$ILA.area))/max(DataQQQ_QC$ILA.area),0)),"%")) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), 
            axis.line = element_line())
    p2 <- ggplot(data = DataQQQ_QC, aes(Sequence, IAld.area)) + 
      geom_point() +
      labs(title = paste0(as.character(round(100*(max(DataQQQ_QC$ILA.conc)-min(DataQQQ_QC$ILA.conc))/max(DataQQQ_QC$ILA.conc),0)),"%")) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), 
            axis.line = element_line())
    
    save_name <- paste0(folder,"5-IndoleQQQ.png")
    
    figure <- ggarrange(p1,p2,
                        labels = c("ILA", "IAld"),
                        ncol = 2, nrow = 1)
    figure <- annotate_figure(figure, top = text_grob(graph_id,
                                                      color = "black", face = "bold", size = 16))
    
    ggsave(plot = figure, width = 8, height = 6, dpi = 400, filename = save_name)
  }
  
  ###
  # Summary of QC
  p3 <- ggplot() +
    annotate("text", x = 0.5, y = 0.95,  label = "Summary of QC", fontface = "bold") +
    annotate("text", x = 0.5, y = 0.90, label = paste0("Sample type: Vial is", VIAL)) +
    annotate("text", x = 0.5, y = 0.80,  label = "Internal standard", fontface = "bold") +
    annotate("text", x = 0.5, y = 0.75, label = QC_Int) +
    annotate("text", x = 0.5, y = 0.65, label = "External standard",fontface = "bold") +
    annotate("text", x = 0.5, y = 0.60,  label = QC_Ext_ILA) +
    annotate("text", x = 0.5, y = 0.55,  label = QC_Ext_IAld) +
    annotate("text", x = 0.5, y = 0.50,  label = QC_Ext_result) +
    annotate("text", x = 0.5, y = 0.1, label = "") +
    theme(
      plot.background  = element_rect(fill = "white", colour = "white"), 
      panel.background = element_rect(fill = "white", colour = "white"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), axis.text.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks = element_blank()
    )
  
  save_name <- paste0(folder, "6-IndoleQQQ-Summary.png")
  ggsave(filename = save_name, plot = p3, width = 5, height = 5, dp = 400, bg = "white")  
  rm(p1,p2,p3)
  rm(DataQQQ_Blank, DataQQQ_Cal, DataQQQ_QC, DataQQQ_ExtStd)
  
  }
 

################################################################################
# Part 3 - Data processing
################################################################################
###
# Remove wells to omits

#DataQQQ <- DataQQQ[!(DataQQQ$DataFile %in% fail_points),]

DataQQQ$group<-as.factor(DataQQQ$sample_id)
#### Remove Calibration from plates
DataQQQ<-DataQQQ[DataQQQ$Level=="",]

## Adjust for dilutions
DataQQQ$ILA.conc <- DataQQQ$ILA.conc*SPDil/(CPDil*ILApurity)
DataQQQ$IAld.conc <- DataQQQ$IAld.conc*SPDil/(CPDil*IAldpurity)
DataQQQ$IPA.conc <- DataQQQ$IPA.conc*SPDil/(CPDil)
DataQQQ$Trp.conc <- DataQQQ$Trp.conc*SPDil/(CPDil)
DataQQQ$IAcr.conc <- DataQQQ$IAcr.conc*SPDil/(CPDil)

###
# Plate view
if (VIAL == FALSE) {
  
  makeStdGraph <- function(p1){
    # Graph all data on a per dilution basis
    p1 <- p1 + coord_equal()
    p1 <- p1 + 
      scale_x_continuous(breaks = 1:24, expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_continuous(breaks = 1:16, labels = LETTERS[1:16], expand = expansion(mult = c(0.01, 0.01)), trans = reverse_trans())
    p1 <- p1 + 
      scale_fill_gradient(low = "white", high = "purple") + 
      labs( x = "Col", y = "Row") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),)
    return(p1)
  }
  
  graph_id <- "Plate view of individual compounds"
  subtitle_name <- "ILA"
  p1 <- ggplot(data = DataQQQ) + 
    geom_circle(aes(x0 = column_96, y0 = row_96, r = 0.45, fill = ILA.conc))
  p1 <- makeStdGraph(p1)
  p1
  
  subtitle_name <- "IPA"
  p2 <- ggplot(data = DataQQQ) + 
    geom_circle(aes(x0 = column_96, y0 = row_96, r = 0.45, fill = IPA.conc))
  p2 <- makeStdGraph(p2)
  p2
  
  save_name <- paste0(folder,"7-IndoleQQQ.png")
  
  figure <- ggarrange(p1,p2,
                      labels = c("ILA", "IPA"),
                      ncol = 1, nrow = 2) + bgcolor("White")
  
  figure <- annotate_figure(figure, top = text_grob(graph_id, 
                                                    color = "black", face = "bold", size = 14))
  
  ggsave(plot = figure, width = 10, height = 10, dpi = 400, filename = save_name)
  
  rm(p1,p2)
}

### 
# Graph view
makeStdGraph <- function(p1){
  # Graph all data on a per dilution basis
  p1 <- p1 +  
    theme(axis.text = element_text(color = "black", size = 18, face = "bold"),
          axis.line = element_line(colour = "black"),
          axis.ticks.length = unit(-0.25, "cm"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 0.5,vjust =0.5),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(color="black", size = 18, face = "bold"),
    )
  return(p1)
}

DataQQQsub <- DataQQQ

p1 <- ggplot(DataQQQsub, aes(group, ILA.conc)) +geom_boxplot()+
  #stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = .25,colour = "gray") +
  #stat_summary(fun = mean, geom = "errorbar", aes(ymin = ..y.., ymax = ..y..), size = 1.2, width = .6, colour = "gray") +
  geom_jitter(aes(colour = InternalQC), width = .35, height = 0, size = 2, alpha = .9) +
  scale_colour_manual(values = c("fail" = "red", "pass" = "black"), guide  = "none") +
  geom_hline(yintercept = ILA_max, linetype = "dotted") +
  labs(y = "ILA titer [mg/L]", x = "Group", title = "ILA titers for all groups")
p1 <- makeStdGraph(p1)

##— Plot with colours driven by InternalQC (red = fail, black = pass) —##
p2 <- ggplot(DataQQQsub, aes(group, IAld.conc)) +geom_boxplot()+
  #stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = .25,colour = "gray") +
  #stat_summary(fun = mean, geom = "errorbar", aes(ymin = ..y.., ymax = ..y..), size = 1.2, width = .6, colour = "gray") +
  geom_jitter(aes(colour = InternalQC), width = .35, height = 0, size = 2, alpha = .9) +
  scale_colour_manual(values = c("fail" = "red", "pass" = "black"), guide  = "none") +
  geom_hline(yintercept = IAld_max, linetype = "dotted") +
  labs( y = "IAld titer [mg/L]", x = "Group", title = "IAld titers for all groups")
p2 <- makeStdGraph(p2)

p5 <- ggplot(data = DataQQQsub, aes(x = group, y = Trp.conc))
p5 <- p5 + geom_boxplot() + geom_point() +geom_hline(yintercept = Trp_max,linetype = 'dotted')+
  labs(y= "Tryptophan titer [mg/L]", x = "Group", title = "Tryptophan titers for all groups") 
p5 <- makeStdGraph(p5)

p6 <- ggplot(data = DataQQQsub, aes(x = group, y = IPA.conc))
p6 <- p6 + geom_boxplot() + geom_point() +geom_hline(yintercept = IPA_max,linetype = 'dotted')+
  labs(y= "IPA titer [mg/L]", x = "Group", title = "IPA titers for all groups") 
p6 <- makeStdGraph(p6)

p7 <- ggplot(data = DataQQQsub, aes(x = group, y = IAcr.conc))
p7 <- p7 + geom_boxplot() + geom_point() +geom_hline(yintercept = IAcr_max,linetype = 'dotted')+
  labs(y= "Indole-Acrylic Acid titer [mg/L]", x = "Group", title = "Indole-Acrylic Acid titers for all groups") 
p7 <- makeStdGraph(p7)

figure <- ggarrange(p1,p5,p2, ncol = 3, nrow = 1) + bgcolor("White")
graph_id <- paste0("Titers of non main compounds per group")
figure <- annotate_figure(figure, top = text_grob(graph_id, 
                                                  color = "black", face = "bold", size = 14))

save_name <- paste0(folder,"8-IndoleQQQ.png")
ggsave(plot = figure, width = 14, height = 8, dpi = 400, filename = save_name)

figure <- ggarrange(p6 ,p7, ncol = 2, nrow = 1) + bgcolor("White")
graph_id <- paste0("Titers of all main compounds per group")
figure <- annotate_figure(figure, top = text_grob(graph_id, 
                                                  color = "black", face = "bold", size = 14))

save_name <- paste0(folder, "9-IndoleQQQ.png")

ggsave(plot = figure, width = 14, height = 8, dpi = 400, filename = save_name)
rm(p1,p2,p5,p6,p7)

write.csv(DataQQQsub, paste0(folder, exp_code,"-processed-P",plate_number, ".csv"))

################################################################################
# Part 4 - Further Data processing
################################################################################
DataQQQsub <- group_by(DataQQQsub, group)
DATAEYE_sum <- summarise(DataQQQsub,
                         meanILA = mean(ILA.conc), sdILA= sd(ILA.conc),
                         meanIAld = mean(IAld.conc), sdIAld= sd(IAld.conc),
                         meanIPA = mean(IPA.conc), sdIPA= sd(IPA.conc),
                         meanTrp = mean(Trp.conc), sdTrp= sd(Trp.conc),
                         meanIAcr = mean(IAcr.conc), sdIAcr= sd(IAcr.conc)
                         )
DATAEYE_sum <-unique(DATAEYE_sum)

DATAEYE_sum <- as.data.frame(lapply(DATAEYE_sum, function(x) {
  if(is.numeric(x)) round(x, 1) else x
}))

DATAEYE_sum_corr <- summarise(DataQQQsub[DataQQQsub$InternalQC == "pass",],
                         meanILA_QC = mean(ILA.conc), sdILA_QC = sd(ILA.conc),
                         meanIAld_QC = mean(IAld.conc), sdIAld_QC = sd(IAld.conc),
                         meanIPA_QC = mean(IPA.conc), sdIPA_QC= sd(IPA.conc),
                         meanTrp_QC = mean(Trp.conc), sdTrp_QC= sd(Trp.conc),
                         meanIAcr_QC = mean(IAcr.conc), sdIAcr_QC= sd(IAcr.conc))
DATAEYE_sum_corr <-unique(DATAEYE_sum_corr)

DATAEYE_sum_corr <- as.data.frame(lapply(DATAEYE_sum_corr, function(x) {
  if(is.numeric(x)) round(x, 1) else x
}))

DATAEYE_sum <- merge(DATAEYE_sum, DATAEYE_sum_corr, by = "group", all = TRUE)

# Specify the width and height
fixed_width <- 25  # Adjust the width as needed
fixed_height <- nrow(DATAEYE_sum) * 0.35  # Adjust the height based on the number of rows

# Save the dataframe as a PNG image
savename <- paste0(folder, exp_code,"_summary-P",plate_number, ".png")
png(savename, width = fixed_width, height = fixed_height, units = "in", res = 450)
grid.table(DATAEYE_sum)
dev.off()

write.csv(DATAEYE_sum, paste0(folder, exp_code,"-summary-P",plate_number, ".csv"))

################################################################################
# Part 5 - Further Data processing for powder or solid samples
################################################################################
### Calculate weight
DataQQQsub$Trp.g<-(DataQQQsub$Trp.conc)*(1/DataQQQsub$concentration.mg.mg)
DataQQQsub$IAcr.g<-(DataQQQsub$IAcr.conc)*(1/DataQQQsub$concentration.mg.mg)
DataQQQsub$ILA.g<-(DataQQQsub$ILA.conc)*(1/DataQQQsub$concentration.mg.mg)
DataQQQsub$IAld.g<-(DataQQQsub$IAld.conc)*(1/DataQQQsub$concentration.mg.mg)
DataQQQsub$IPA.g<-(DataQQQsub$IPA.conc)*(1/DataQQQsub$concentration.mg.mg)

p1 <- ggplot(DataQQQsub, aes(group, IPA.g)) +geom_boxplot()+
  #stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = .25,colour = "gray") +
  #stat_summary(fun = mean, geom = "errorbar", aes(ymin = ..y.., ymax = ..y..), size = 1.2, width = .6, colour = "gray") +
  geom_jitter(aes(colour = InternalQC), width = .35, height = 0, size = 2, alpha = .9) +
  scale_colour_manual(values = c("fail" = "red", "pass" = "black"), guide  = "none") +
  #geom_hline(yintercept = ILA_max, linetype = "dotted") +
  labs(y = "IPA titer [mg/g]", x = "Group", title = "IPA titers for all groups")
p1 <- makeStdGraph(p1)

##— Plot with colours driven by InternalQC (red = fail, black = pass) —##
p2 <- ggplot(DataQQQsub, aes(group, IAcr.g)) +geom_boxplot()+
  #stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = .25,colour = "gray") +
  #stat_summary(fun = mean, geom = "errorbar", aes(ymin = ..y.., ymax = ..y..), size = 1.2, width = .6, colour = "gray") +
  geom_jitter(aes(colour = InternalQC), width = .35, height = 0, size = 2, alpha = .9) +
  scale_colour_manual(values = c("fail" = "red", "pass" = "black"), guide  = "none") +
  #geom_hline(yintercept = IAld_max, linetype = "dotted") +
  labs( y = "Indole Acrylic Acid titer[mg/g]", x = "Group", title = "Indole Acrylic Acid titers for all groups")
p2 <- makeStdGraph(p2)

figure <- ggarrange(p1 ,p2, ncol = 2, nrow = 1) + bgcolor("White")
graph_id <- paste0("Titers of all main compounds per group")
figure <- annotate_figure(figure, top = text_grob(graph_id, 
                                                  color = "black", face = "bold", size = 14))

save_name <- paste0(folder, "9b-IndoleQQQ.png")

ggsave(plot = figure, width = 14, height = 8, dpi = 400, filename = save_name)
rm(p1,p2)
write.csv(DataQQQsub, paste0(folder, exp_code,"-processed-weight-P",plate_number, ".csv"))

DataQQQsub <- group_by(DataQQQsub, group)
DATAEYE_sum <- summarise(DataQQQsub,
                         meanIPA_g = mean(IPA.g), sdIPA_g= sd(IPA.g),
                         meanIAcr_g = mean(IAcr.g), sdSkatole_g= sd(IAcr.g))
DATAEYE_sum <-unique(DATAEYE_sum)

DATAEYE_sum <- as.data.frame(lapply(DATAEYE_sum, function(x) {
  if(is.numeric(x)) round(x, 8) else x
}))

DATAEYE_sum_corr <- summarise(DataQQQsub[DataQQQsub$InternalQC == "pass",],
                              meanIPA_QC = mean(IPA.g), sdIPA_QC = sd(IPA.g),
                              meanIAcr_QC = mean(IAcr.g), sdIAcr_QC = sd(IAcr.g))
DATAEYE_sum_corr <-unique(DATAEYE_sum_corr)

DATAEYE_sum_corr <- as.data.frame(lapply(DATAEYE_sum_corr, function(x) {
  if(is.numeric(x)) round(x, 8) else x
}))

DATAEYE_sum <- merge(DATAEYE_sum, DATAEYE_sum_corr, by = "group", all = TRUE)

# Specify the width and height
fixed_width <- 15  # Adjust the width as needed
fixed_height <- nrow(DATAEYE_sum) * 0.35  # Adjust the height based on the number of rows

# Save the dataframe as a PNG image
savename <- paste0(folder, exp_code,"_summary_gg-P",plate_number, ".png")
png(savename, width = fixed_width, height = fixed_height, units = "in", res = 450)
grid.table(DATAEYE_sum)
dev.off()

write.csv(DATAEYE_sum, paste0(folder, exp_code,"-summary_gg-P",plate_number, ".csv"))

