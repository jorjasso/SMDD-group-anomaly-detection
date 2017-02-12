#Box plots for the experiments
rm(list=ls())
setwd("~/Documents/GITProjects/SMDD-group-anomaly-detection")
library(dplyr,warn.conflicts = FALSE)
library(ggplot2)
library(lazyeval)
library(geometry)
library(openxlsx)
library(ggvis)
library(caret)
library(randomForest)
library(mlbench)
library(caretEnsemble)
library(rpart)
library(parallel)
library(doParallel)
library(glmnet)
require(gridExtra)
library(hydroGOF)
library(reshape2)

#DATASET SLOAN
#-------------
path_1<-"output/workspaceSloanRandomFirstTypeAnomalies1.mat.csv"
path_2<-"output/workspaceSloanSecondTypeAnomalies2.mat.csv"
path_3<-"output/workspaceSloanThirdTypeAnomalies3.mat.csv"

data_set_1<-tbl_df(read.csv(path_1,head=FALSE,sep=","))
data_set_2<-tbl_df(read.csv(path_2,head=FALSE,sep=","))
data_set_3<-tbl_df(read.csv(path_3,head=FALSE,sep=","))

names(data_set_1)<-c("SMDD.C.k.1", "SMDD", "SMDD.N", "OCSMM", "SVDD", 
                     "SMDD.C.k.9","SMDD.C.k.92","SMDD.C.k.94","SMDD.C.k.96","SMDD.C.k.98",
                     "SMDD.C.N.k.1","SMDD.C.N.k.9","SMDD.C.N.k.92","SMDD.C.N.k.94","SMDD.CP.N.k.96","SMDD.C.N.k98"
                     )

names(data_set_2)<-names(data_set_1)
names(data_set_3)<-names(data_set_1)

#data_set_1
#data_set_2
df<-data_set_3 %>%
  dplyr::select(SVDD,OCSMM,SMDD,SMDD.C.k.1,SMDD.N,SMDD.C.N.k.1)

data<-melt(as.data.frame(df))

data<-data %>%
  mutate(G=ifelse(variable%in% list("SVDD","OCSMM","SMDD","SMDD.C.k.1"),1,0))

data %>% ggvis(~factor(variable), ~value,stroke=~G) %>% layer_boxplots() %>%
scale_numeric("y", domain = c(0.4,1)) %>%
add_axis("y", 
         title = "AUC",
         title_offset = 50,
         subdivide = 1,
         properties = axis_props(labels = list(fontSize = 15),title = list(fontSize = 15)))  %>%
  add_axis("x",
           title = "Models",
           subdivide = 0,
           title_offset = 90,
           properties = axis_props(labels = list(fontSize = 15,angle = 45,align = "left", baseline = "middle"),
                                   title = list(fontSize = 15))) %>%
  add_legend("stroke", title = "Kernel",
             properties = legend_props(
               title = list(fontSize = 16),
               labels = list(fontSize = 14, dx = 5),
               symbol = list(stroke = "black", strokeWidth = 2,
                             shape = "square", size = 200),
               legend = list(
                 x = scaled_value("x", 4.5),
                 y = scaled_value("y", 30)
               )
             ))
