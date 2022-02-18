#Load the packages
library(flowCore)
library(flowAI)
library(ggcyto)
library(tidyverse)
# library(cytoverse)

#How to get help
??flowCore

#Load a single file
fcs1 <- read.FCS("FlowRepository_FR-FCM-ZZZU_files/0001.FCS")
fcs1
summary(fcs1)
names(fcs1)
exprs(fcs1)
each_col(fcs1, median)
keyword(fcs1)
keyword(fcs1)$FILENAME
keyword(fcs1)$FCSversion
keyword(fcs1)$SPILL

#Compensation
spillover(fcs1) # compensation matrix calculated from single stain controls (stored in fcs files)

fcs1_comp <- fcs1 %>% 
  spillover() %>% 
  .$SPILL %>% 
  compensate(fcs1, .)

fcs1_comp@exprs

#Cleaning: clean intermittent events on time scale
# For a set of FCS files, flow_auto_qc performs a complete and automatic quality control. 
# It consists in the detection and removal of anomalies by checking three properties of flow cytometry:
# 1) flow rate, 2) signal acquisition, 3) dynamic range.


fcs1_comp_clean <- fcs1_comp %>% 
  flow_auto_qc()

fcs1_comp_clean

keyword(fcs1_comp_clean) <- keyword(fcs1)
fcs1_comp_clean
??flowAI

#Transformation: FCS data are on linear scale, we need to transform the data on log scale for visualisation except time, FSC, SSC
colnames(fcs1_comp_clean)
exprs(fcs1_comp_clean)
trans <- fcs1_comp_clean[, 3:10] %>% 
  colnames() %>% 
  estimateLogicle(fcs1_comp_clean, .)

fcs1_comp_clean_log <- fcs1_comp_clean[, 3:10] %>% 
  colnames() %>% 
  estimateLogicle(fcs1_comp_clean, .) %>% 
  flowCore::transform(fcs1_comp_clean, .)

#Visualise the results
??ggcyto
autoplot(fcs1_comp_clean) # without log transformation it's difficult to visualise fluor channels
autoplot(fcs1_comp_clean_log)

fcs1_comp_clean_log
fcs1_comp_clean_log@parameters@data

fcs1_comp %>% 
  autoplot(x = "Time", y = "FSC-A", bins = 256)

fcs1_comp_clean_log %>% 
  autoplot(x = "Time", y = "FSC-A", bins = 256)

fcs1_comp_clean_log %>% 
  autoplot(x = "FSC-A", y = "SSC-A", bins = 256)

fcs1_comp_clean_log %>% 
  autoplot(x = "MHCII", y = "CD14", bins = 256)

fcs1_comp_clean_log %>% 
  autoplot(x = "MHCII", y = "CD14", bins = 300)

fcs1_comp_clean_log %>% 
  autoplot(x = "CD14", y = "CD11c", bins = 128)

fcs1_comp_clean_log %>% 
  autoplot(x = c("CD123", "IL12"))

autoplot(fcs1_comp_clean_log, x="PE-Cy7-A", y="PerCP-Cy5-5-A", bins = 256)
autoplot(fcs1_comp_clean_log, x="Time", y="FSC-A", bins = 128)
autoplot(transform(fcs1_comp,trans), x="Time", y="FSC-A", bins = 128)

#In a flowSet
myfiles <- list.files(path="data_fcs", pattern=".FCS$")
fcs.all <- flowCore::read.flowSet(myfiles, path="data_fcs/")
fcs.all
fcs.all[[1]]
spillover(fcs.all[[1]])

fcs_clean <- fcs.all %>% 
  # obtain compenstation matrix/spill matrix
  fsApply(
    function(.x){.x %>% spillover() %>% .$SPILL},
    simplify = FALSE
  ) %>% 
  # do compensation
  mapply(compensate, fcs.all, ., SIMPLIFY = FALSE) %>% 
  as(., "flowSet") %>% 
  # do quality control: remove low quality data
  flow_auto_qc()

fcs_clean_log <- fcs_clean %>% 
  # create log transformation object
  fsApply(estimateLogicle, 
          colnames(fcs_clean)[!colnames(fcs_clean) %in% c("FSC-A", "SSC-A", "Time")], 
          simplify = FALSE) %>% 
  # apply the transformation on the clean data
  flowCore::transform(fcs_clean, .)

fs_comp <-compensate(fs, spillover(fs[[1]])$SPILL)
fs_comp_clean <- flow_auto_qc(fs_comp)
trans <- estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[[1]][,3:10]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)

#fsApply
??fsApply
fsApply(fcs.all, each_col, median)
