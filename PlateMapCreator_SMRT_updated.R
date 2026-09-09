## Created: 20230424
# Last edited: 20260806 by AF
# Author: Anthony Frederick
# Department: Microscopy & Optics(M&O)
# Purpose: Accept csv file of a platemap and add an experimental code,to quickly
# adjust the number of replicates, and code number to match experimental card
# file for further processing.
### Set filename
###############################################################################
# Inputs to change
### 
# The plate file should have the following columns
# sample_id, well, row,column, plate_id
plate_file = "https://docs.google.com/spreadsheets/d/1VdTedIb0KNoxs9qdP_AxM-9WGffc0zGrxH6wQQpPHpk/edit?usp=sharing"

#CHANGE to experiment name
expCode <-"Sample_test"
nrep<-1# Used for determining number of replicates
# Used for randomization of sample, change only when doing a DIFFERENT plate submission
# choose whatever number you want
seed=sample.int(99999,1)

###############################################################################

# load libraries
library(ggplot2); library(tidyr); library(dplyr); library(ggforce); library(scales);
library(DescTools); library(PMCMRplus); library(ggpubr); library(tidyverse);

####
# Google Sheet to publish tabs into ("nobo-val-plate-map"; do not change).
link = "https://docs.google.com/spreadsheets/d/16B28cVPxArJWAeVUA8W80AfhYGk1OmWuKiBl3PGNfUQ/edit?usp=sharing"
plate_map_dir = Sys.glob("~/Library/CloudStorage/GoogleDrive-*@kingdomsuperculture.com/Shared drives/Science/Platform/AssayManagerData/IPA-LCMS-QQQ/IPA-LCMS-QQQ-platemap/")
biomek_file_dir = Sys.glob("~/Library/CloudStorage/GoogleDrive-*@kingdomsuperculture.com/Shared drives/Ingredient-Development/1. Active Projects/SMRT/Notebook/SMRT-QQQ/Biomek-files/")
#' Read a submitted plate from a Google Sheet URL, a Drive link, or a CSV path.
read_plate_file <- function(plate_file) {
  if (grepl("docs.google.com/spreadsheets", plate_file, fixed = TRUE)) {
    if (!requireNamespace("googlesheets4", quietly = TRUE)) {
      stop("googlesheets4 is required to read a Google Sheets URL.")
    }
    as.data.frame(googlesheets4::read_sheet(ss = plate_file))
    
  } else if (grepl("drive.google.com", plate_file, fixed = TRUE)) {
    if (!requireNamespace("googledrive", quietly = TRUE)) {
      stop("googledrive is required to read a Google Drive link.")
    }
    file_id <- regmatches(
      plate_file,
      regexpr("(?<=id=)[^&]+|(?<=/d/)[^/?]+", plate_file, perl = TRUE)
    )
    tmp <- tempfile(fileext = ".csv")
    googledrive::drive_download(googledrive::as_id(file_id), path = tmp, overwrite = TRUE)
    read.csv(tmp)
    
  } else if (grepl("\\.csv$", plate_file, ignore.case = TRUE)) {
    read.csv(plate_file)
    
  } else {
    stop("Unsupported input: provide a Google Sheets URL, Google Drive link, or a .csv path/URL.")
  }
}
push_to_sheet <- function(df, link, sheet_name, dry_run) {
  if (dry_run) {
    message(sprintf("[dry-run] would write Google Sheet tab: %s", sheet_name))
    return(invisible(NULL))
  }
  if (!requireNamespace("googlesheets4", quietly = TRUE)) {
    stop("googlesheets4 is required to write Google Sheet tabs.")
  }
  googlesheets4::sheet_add(ss = link, sheet = sheet_name)
  googlesheets4::sheet_write(data = df, ss = link, sheet = sheet_name)
  invisible(NULL)
}


# read source plate map
plate_map <- read_plate_file(plate_file)

plate_map$source_plate<-"source_plate"
all_wells<-plate_map$well
only_sample<-plate_map[plate_map$sample_id!="",]
### Reactor Optimization
# Remove T-1 
only_sample<-only_sample[!grepl(x=only_sample$sample_id,pattern = "_T-1"),]


num_sample<-length(only_sample$sample_id)
num_samp_rep<-num_sample*nrep

num_plates<-ceiling(num_samp_rep/96)
num_of_rep_plate<-floor(96/num_sample)
### Create empty data frame
test<-plate_map[plate_map$column==30,]
## Fill dataframe with the designated number of replicates
# replicate count
for (r in 1:nrep) {
  if(r<=nrep){
    test<-rbind(test,only_sample)
  }
}

### 
test$target_plate<-"target_plate"
test$well_96<-"A1"
test$volume<-100
random_file<-test[test$column==30,]
random_file<-random_file[,c("target_plate","source_plate","well_96","well")]
colnames(random_file)<-c("target_plate","source_plate","target_well","source_well")
name_plate<-test$plate_id[1]
set.seed(seed) # Set the seed for the randomization
for (c in 1:num_plates) {
  # Set the seed for so that randomization can be repeated within an experiment 
  # seed changes for each plate so randomization is different
  set.seed(seed+c-1)
  temp<-head(test, 96)
  temp$well_96<-sample(all_wells, size = length(temp$well))
  temp$row_96<-substr(temp$well_96,1,1)
  temp$column_96<-substr(temp$well_96,2,3)
  temp$target_plate<-paste0("BPO_plate_",c)
  write.csv(temp,paste(plate_map_dir,expCode,"_plate_map_AssayFormat_P",c,".csv",sep = ""))
  temp_sheet <- temp
  temp_sheet$plate_id<-paste0(name_plate,"_plate_map_AssayFormat_P",c)
  push_to_sheet(temp_sheet, link, paste0(name_plate,"_plate_map_AssayFormat_P",c), FALSE)
  
  temp<-temp[temp$sample_id!="",]
  #Re-label
  temp<-temp[,c("target_plate","source_plate","well_96","well")]
  colnames(temp)<-c("target_plate","source_plate","target_well","source_well")
  random_file<-rbind(random_file,temp)
  #remove 96 from test
  if(length(test$sample_id)>=96){
    test<-test[-(1:96), ]
  } else {
    test<-test[-(1:length(test$sample_id)), ]
  }
}

write.csv(random_file,paste0(biomek_file_dir,expCode,"_plate_map_Biomek_file_ALL.csv"))

msg <- sprintf('display dialog "✅ Script completed successfully! There are %d target plates needed for randomization on biomek"   buttons {"OK"} with title "All Done" with icon note', num_plates)
system(paste("osascript -e", shQuote(msg)))



