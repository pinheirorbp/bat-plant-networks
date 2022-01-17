# define the working directory
# in R Studio you can do it using the following line of code
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# what do run?
run_vaznull=F
run_restnull=F
PLOT=T
# packages and data files
library(bipartite)
library(igraph)
matrices= list(
  dryfor_dry=read.table("dryforest.dryseason.txt"),
  dryfor_wet=read.table("dryforest.wetseason.txt"),
  wetfor_dry=read.table("laselva.dryseason.txt"),
  wetfor_wet=read.table("laselva.wetseason.txt"),
  dryfor_year=read.table("dryforest.wholeyear.txt"),
  wetfor_year=read.table("laselva.wholeyear.txt"),
  dryfor_normal_wet=read.table("dryforest.normalwet.txt")
)
# Observed network structure
NET_INDICES=data.frame(row.names = names(matrices))
MOD=list()
binMOD=list()
for (i in 1:7){
  wmatrix=matrices[[i]]
  binmatrix=wmatrix
  binmatrix[binmatrix>0]=1
  # 1- modularity
  # weighted
  MOD[[i]]=DIRT_LPA_wb_plus(wmatrix)
  NET_INDICES$WMOD[i]=MOD[[i]]$modularity
  # binary
  binMOD[[i]]=DIRT_LPA_wb_plus(binmatrix)
  NET_INDICES$binMOD[i]=binMOD[[i]]$modularity
  # 2- overall nestedness
  # weighted
  WNODA=nest.smdm(wmatrix, weighted = T,decreasing = "abund", constraints = c(MOD[[i]]$Row_labels,MOD[[i]]$Col_labels))
  WNODF=nest.smdm(wmatrix, weighted = T,decreasing = "fill", constraints = c(MOD[[i]]$Row_labels,MOD[[i]]$Col_labels))
  NET_INDICES$WNODA[i]=WNODA$WNODAmatrix
  NET_INDICES$WNODF[i]=WNODF$WNODFmatrix
  # binary
  NODF=nest.smdm(binmatrix, weighted = F, constraints = c(binMOD[[i]]$Row_labels,binMOD[[i]]$Col_labels))
  NET_INDICES$NODF[i]=NODF$NODFmatrix
  #3- connectance
  NET_INDICES$connectance[i]=networklevel(wmatrix,index = "connectance")
  #4- dimensions
  NET_INDICES$nrow[i]=nrow(wmatrix)
  NET_INDICES$ncol[i]=ncol(wmatrix)
  NET_INDICES$dim[i]=NET_INDICES$ncol[i]+NET_INDICES$nrow[i]
  #5- compartments
  NET_INDICES$compartments[i]=networklevel(wmatrix,index = "number of compartments")
  #6- robustness
  NET_INDICES$robustnessHL[i]=networklevel(wmatrix,index ="robustness")[1]
  NET_INDICES$robustnessLL[i]=networklevel(wmatrix,index ="robustness")[2]
  #9- Nestedness SM and DM
  # weighted
  NET_INDICES$WNODA_SM[i]=WNODA$WNODA_SM_matrix
  NET_INDICES$WNODA_DM[i]=WNODA$WNODA_DM_matrix
  NET_INDICES$WNODF_SM[i]=WNODF$WNODF_SM_matrix
  NET_INDICES$WNODF_DM[i]=WNODF$WNODF_DM_matrix
  # binary
  NET_INDICES$NODF_SM[i]=NODF$NODF_SM_matrix
  NET_INDICES$NODF_DM[i]=NODF$NODF_DM_matrix
}
### Lista de espécies ####
species=c(sort(unique(c(rownames(matrices[[5]]),rownames(matrices[[6]])))), sort(unique(c(colnames(matrices[[5]]),colnames(matrices[[6]])))))
##### PLOTS #####
### Dry Forest Whole Year ####
net_dryfor_year=graph_from_incidence_matrix(matrices[[5]], weighted = T)
# L_dryfor=layout_with_fr(net_dryfor_year,weights = NULL)
# L_dryfor=norm_coords(L_dryfor)
#save(L_dryfor, file="L_dryfor.RData")
# The position of nodes vary each time one run the algorithm
# Thus we saved the organization applied in the figure at the paper
load("L_dryfor.RData")
V(net_dryfor_year)$size=c(rowSums(matrices[[5]]),colSums(matrices[[5]]))^(1/2)*2.5
V(net_dryfor_year)$color<-c("#08995B","#DF480C")[1+V(net_dryfor_year)$type]
V(net_dryfor_year)$shape<-c("square","circle")[1+V(net_dryfor_year)$type]
E(net_dryfor_year)$width=E(net_dryfor_year)$weight
## Plots are saved in the folder "plots"
if(PLOT){
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/dry_forest_year.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_dryfor_year,vertex.label=NA, layout=L_dryfor, rescale=F)
  text(0,-1.2, "Dry Forest - Whole Year")
  dev.off()
  dev.off()
  ##
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/lab_dry_forest_year.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_dryfor_year,vertex.label.cex=.8, vertex.label.dist=.6, vertex.label.color="black", vertex.label.rotation= pi, layout=L_dryfor, rescale=F)
  text(0,-1.2, "Dry Forest - Whole Year")
  dev.off()
  dev.off()
}
##
##
net_dryfor_year2=net_dryfor_year
V(net_dryfor_year2)$name=as.character(match(names(V(net_dryfor_year)),species))
par(mar=c(1,1,1,1), bg=NA)
V(net_dryfor_year2)$color<-c("#CFEBDF","#F7DCD0")[1+V(net_dryfor_year)$type]
##
if(PLOT){
  png(filename = "plots/networks/labn_dry_forest_year.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_dryfor_year2,vertex.label.cex=1, vertex.label.dist=.7
       ,vertex.label.color="black", layout=L_dryfor, rescale=F, vertex.frame.color= "gray60",edge.color="gray90", vertex.label.font=2)
  text(0,-1.2, "Dry Forest - Whole Year")
  dev.off()
  dev.off()
}
## Dry Forest Dry Season
net_dryfor_dry=graph_from_incidence_matrix(matrices[[1]], weighted = T)
V(net_dryfor_dry)$size=c(rowSums(matrices[[1]]),colSums(matrices[[1]]))^(1/2)*2.5
V(net_dryfor_dry)$color<-c("#08995B","#DF480C")[1+V(net_dryfor_dry)$type]
V(net_dryfor_dry)$shape<-c("square","circle")[1+V(net_dryfor_dry)$type]
E(net_dryfor_dry)$width=E(net_dryfor_dry)$weight
L_dryfor_dry=L_dryfor[match(names(V(net_dryfor_dry)),names(V(net_dryfor_year))),]
##
if(PLOT){
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/dry_forest_dry.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_dryfor_dry,vertex.label=NA, layout=L_dryfor_dry, rescale=F)
  text(0,-1.2, "Dry Forest - Dry Season")
  dev.off()
  dev.off()
  ##
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/lab_dry_forest_dry.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_dryfor_dry,vertex.label.cex=.8, vertex.label.dist=.6, vertex.label.color="black", vertex.label.rotation= pi, layout=L_dryfor_dry, rescale=F)
  text(0,-1.2, "Dry Forest - Dry Season")
  dev.off()
  dev.off()
}
## Dry Forest Wet Season
net_dryfor_wet=graph_from_incidence_matrix(matrices[[2]], weighted = T)
V(net_dryfor_wet)$size=c(rowSums(matrices[[2]]),colSums(matrices[[2]]))^(1/2)*2.5
V(net_dryfor_wet)$color<-c("#08995B","#DF480C")[1+V(net_dryfor_wet)$type]
V(net_dryfor_wet)$shape<-c("square","circle")[1+V(net_dryfor_wet)$type]
E(net_dryfor_wet)$width=E(net_dryfor_wet)$weight
L_dryfor_wet=L_dryfor[match(names(V(net_dryfor_wet)),names(V(net_dryfor_year))),]
##
if(PLOT){
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/dry_forest_wet.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_dryfor_wet,vertex.label=NA, layout=L_dryfor_wet, rescale=F)
  text(0,-1.2, "Dry Forest - Wet Season")
  dev.off()
  dev.off()
  ##
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/lab_dry_forest_wet.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_dryfor_wet,vertex.label.cex=.8, vertex.label.dist=.6, vertex.label.color="black", vertex.label.rotation= pi, layout=L_dryfor_wet, rescale=F)
  text(0,-1.2, "Dry Forest - Wet Season")
  dev.off()
  dev.off()
}
## Dry Forest Wet Season - Normal Year
net_dryfor_normal_wet=graph_from_incidence_matrix(matrices[[7]], weighted = T)
V(net_dryfor_normal_wet)$size=c(rowSums(matrices[[7]]),colSums(matrices[[7]]))^(1/2)*2.5
V(net_dryfor_normal_wet)$color<-c("#08995B","#DF480C")[1+V(net_dryfor_normal_wet)$type]
V(net_dryfor_normal_wet)$shape<-c("square","circle")[1+V(net_dryfor_normal_wet)$type]
E(net_dryfor_normal_wet)$width=E(net_dryfor_normal_wet)$weight

#L_dryfor_normal=layout_with_fr(net_dryfor_normal_wet,weights = NULL)
#L_dryfor_normal=norm_coords(L_dryfor_normal)
#save(L_dryfor_normal, file="L_dryfor_normal.RData")
load("L_dryfor_normal.RData")

##
if(PLOT){
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/dry_forest_normal_wet.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_dryfor_normal_wet,vertex.label=NA, layout=L_dryfor_normal, rescale=F)
  text(0,-1.2, "Dry Forest - Normal Year - Wet Season")
  dev.off()
  dev.off()
  ##
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/lab_dry_forest_normal_wet.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_dryfor_normal_wet,vertex.label.cex=.8, vertex.label.dist=.6, vertex.label.color="black", vertex.label.rotation= pi, layout=L_dryfor_normal, rescale=F)
  text(0,-1.2, "Dry Forest - Normal Year - Wet Season")
  dev.off()
  dev.off()
}
### Wet Forest Whole Year ####
net_wetfor_year=graph_from_incidence_matrix(matrices[[6]], weighted = T)
# L_wetfor=layout_with_fr(net_wetfor_year,weights = NULL)
# L_wetfor=norm_coords(L_wetfor)
# save(L_wetfor, file="L_wetfor.RData")
load("L_wetfor.RData")
V(net_wetfor_year)$size=c(rowSums(matrices[[6]]),colSums(matrices[[6]]))^(1/2)*2.5
V(net_wetfor_year)$color<-c("#08995B","#DF480C")[1+V(net_wetfor_year)$type]
V(net_wetfor_year)$shape<-c("square","circle")[1+V(net_wetfor_year)$type]
E(net_wetfor_year)$width=E(net_wetfor_year)$weight
##
if(PLOT){
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/wet_forest_year.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_wetfor_year,vertex.label=NA, layout=L_wetfor, rescale=F)
  text(0,-1.2, "Wet Forest - Whole Year")
  dev.off()
  dev.off()
  ##
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/lab_wet_forest_year.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_wetfor_year,vertex.label.cex=.8, vertex.label.dist=.6, vertex.label.color="black", vertex.label.rotation= pi, layout=L_wetfor, rescale=F)
  text(0,-1.2, "Wet Forest - Whole Year")
  dev.off()
  dev.off()
}
##
net_wetfor_year2=net_wetfor_year
V(net_wetfor_year2)$name=as.character(match(names(V(net_wetfor_year)),species))
V(net_wetfor_year2)$color<-c("#CFEBDF","#F7DCD0")[1+V(net_wetfor_year)$type]
##
if(PLOT){
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/labn_wet_forest_year.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_wetfor_year2,vertex.label.cex=1, vertex.label.dist=.7
       ,vertex.label.color="black", layout=L_wetfor, rescale=F, vertex.frame.color= "gray60",edge.color="gray90", vertex.label.font=2)
  text(0,-1.2, "Wet Forest - Whole Year")
  dev.off()
  dev.off()
}
## Wet Forest Dry Season
net_wetfor_dry=graph_from_incidence_matrix(matrices[[3]], weighted = T)
V(net_wetfor_dry)$size=c(rowSums(matrices[[3]]),colSums(matrices[[3]]))^(1/2)*2.5
V(net_wetfor_dry)$color<-c("#08995B","#DF480C")[1+V(net_wetfor_dry)$type]
V(net_wetfor_dry)$shape<-c("square","circle")[1+V(net_wetfor_dry)$type]
E(net_wetfor_dry)$width=E(net_wetfor_dry)$weight
L_wetfor_dry=L_wetfor[match(names(V(net_wetfor_dry)),names(V(net_wetfor_year))),]
##
if(PLOT){
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/wet_forest_dry.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_wetfor_dry,vertex.label=NA, layout=L_wetfor_dry, rescale=F)
  text(0,-1.2, "Wet Forest - Dry Season")
  dev.off()
  dev.off()
  ##
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/lab_wet_forest_dry.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_wetfor_dry,vertex.label.cex=.8, vertex.label.dist=.6, vertex.label.color="black", vertex.label.rotation= pi, layout=L_wetfor_dry, rescale=F)
  text(0,-1.2, "Wet Forest - Dry Season")
  dev.off()
  dev.off()
}
## Wet Forest Wet Season
net_wetfor_wet=graph_from_incidence_matrix(matrices[[4]], weighted = T)
V(net_wetfor_wet)$size=c(rowSums(matrices[[4]]),colSums(matrices[[4]]))^(1/2)*2.5
V(net_wetfor_wet)$color<-c("#08995B","#DF480C")[1+V(net_wetfor_wet)$type]
V(net_wetfor_wet)$shape<-c("square","circle")[1+V(net_wetfor_wet)$type]
E(net_wetfor_wet)$width=E(net_wetfor_wet)$weight
L_wetfor_wet=L_wetfor[match(names(V(net_wetfor_wet)),names(V(net_wetfor_year))),]
##
if(PLOT){
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/wet_forest_wet.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_wetfor_wet,vertex.label=NA, layout=L_wetfor_wet, rescale=F)
  text(0,-1.2, "Wet Forest - Wet Season")
  dev.off()
  dev.off()
  ##
  par(mar=c(1,1,1,1), bg=NA)
  png(filename = "plots/networks/lab_wet_forest_wet.png",width = 15,height = 15,units = "cm",res = 600, bg = NA)
  plot(net_wetfor_wet,vertex.label.cex=.8, vertex.label.dist=.6, vertex.label.color="black", vertex.label.rotation= pi, layout=L_wetfor_wet, rescale=F)
  text(0,-1.2, "Wet Forest - Wet Season")
  dev.off()
  dev.off()
}
### FIGURA 2 ####
if(PLOT){
  svg(filename = "plots/networks/fig2.svg",width = 50/2.5,height = 35/2.5,bg = NA)
  layout(matrix(c(1,6,0,7,0,8,
                  0,4,0,3,0,5,
                  2,9,0,10,0,11),byrow = T,nrow = 3,ncol =6),widths =c(.8,3,.2,3,.2,3),heights = c(3,.6,3))
  # 1
  par(mar=c(0,0,0,0), bg=NA)
  plot.new()
  text(x=0.5,y=0.5,cex=3,"DRY FOREST",srt=90,)
  # 2
  plot.new()
  text(x=0.5,y=0.5,cex=3,"RAINFOREST",srt=90,)
  # 3
  plot.new()
  rect(xleft = .3,xright = .7,ybottom = .35,ytop = .65, lwd=3, border ="#C55A11")
  text(x=0.5,y=0.5,cex=2.5,"DRY SEASON", col ="#C55A11" )
  # 4
  plot.new()
  rect(xleft = .3,xright = .7,ybottom = .35,ytop = .65, lwd=3)
  text(x=0.5,y=0.5,cex=2.5,"WHOLE YEAR")
  # 5
  plot.new()
  rect(xleft = .3,xright = .7,ybottom = .35,ytop = .65, lwd=3, border ="#2E75B6")
  text(x=0.5,y=0.5,cex=2.5,"WET SEASON", col ="#2E75B6" )
  #6
  plot(net_dryfor_year2,vertex.label.cex=2.5, vertex.label.dist=.7
       ,vertex.label.color="black", layout=L_dryfor, rescale=F, vertex.frame.color= "gray60",edge.color="gray90", vertex.label.font=2)
  #7
  plot(net_dryfor_dry,vertex.label=NA, layout=L_dryfor_dry, rescale=F)
  #8
  plot(net_dryfor_wet,vertex.label=NA, layout=L_dryfor_wet, rescale=F)
  #9
  plot(net_wetfor_year2,vertex.label.cex=2.5, vertex.label.dist=.7
       ,vertex.label.color="black", layout=L_wetfor, rescale=F, vertex.frame.color= "gray60",edge.color="gray90", vertex.label.font=2)
  #10
  plot(net_wetfor_dry,vertex.label=NA, layout=L_wetfor_dry, rescale=F)
  #11
  plot(net_wetfor_wet,vertex.label=NA, layout=L_wetfor_wet, rescale=F)
  dev.off()
}
### Dissimilarity of interactions ####
library(betalink)
networks=list(net_dryfor_dry, net_dryfor_wet, net_wetfor_dry, net_wetfor_wet, net_dryfor_year, net_wetfor_year, net_dryfor_normal_wet)
names(networks)=names(matrices)
BETA0=network_betadiversity(networks)
BETA=BETA0[c(1,11,12,19),]
BETA
### Null models ####
if(run_vaznull){
  VAZNULLMATRICES=list()
  for (i in 1:7){
    VAZNULLMATRICES[[i]]=vaznull(matrices[[i]],N=1000)
  }
  VAZNULLRESULTS=list()
  VAZNULL_INDICES=data.frame(row.names = names(matrices))
  for (i in 1:7){
    VAZNULLRESULTS[[i]]=list()
    names(VAZNULLRESULTS)[i]=names(matrices[i])
    print(paste("Starting null model:", i, names(matrices[i])))
    for(j in 1:1000){
      if(j%%10==0){print(paste("null matrix n:",j))}
      NULLMAT=VAZNULLMATRICES[[i]][[j]]
      # 1- modularity
      VAZNULLRESULTS[[i]]$MOD[[j]]=DIRT_LPA_wb_plus(NULLMAT)
      VAZNULLRESULTS[[i]]$WMOD[j]=VAZNULLRESULTS[[i]]$MOD[[j]]$modularity
      # 2- overall nestedness
      # weighted
      WNODA=nest.smdm(NULLMAT, weighted = T,decreasing = "abund", constraints = c(MOD[[i]]$Row_labels,MOD[[i]]$Col_labels))
      WNODF=nest.smdm(NULLMAT, weighted = T,decreasing = "fill", constraints = c(MOD[[i]]$Row_labels,MOD[[i]]$Col_labels))
      VAZNULLRESULTS[[i]]$WNODA[j]=WNODA$WNODAmatrix
      VAZNULLRESULTS[[i]]$WNODF[j]=WNODF$WNODFmatrix
      #3- connectance
      VAZNULLRESULTS[[i]]$connectance[j]=networklevel(NULLMAT,index = "connectance")
      #4- dimensions
      VAZNULLRESULTS[[i]]$nrow[j]=nrow(NULLMAT)
      VAZNULLRESULTS[[i]]$ncol[j]=ncol(NULLMAT)
      VAZNULLRESULTS[[i]]$dim[j]=sum(dim(NULLMAT))
      #5- compartments
      VAZNULLRESULTS[[i]]$compartments[j]=networklevel(NULLMAT,index = "number of compartments")
    }
    # Summary of null models
    VAZNULL_INDICES$WMOD_mean[i]=mean(VAZNULLRESULTS[[i]]$WMOD)
    VAZNULL_INDICES$WMOD_sd[i]=sd(VAZNULLRESULTS[[i]]$WMOD)
    VAZNULL_INDICES$WNODA_mean[i]=mean(VAZNULLRESULTS[[i]]$WNODA)
    VAZNULL_INDICES$WNODA_sd[i]=sd(VAZNULLRESULTS[[i]]$WNODA)
    VAZNULL_INDICES$WNODF_mean[i]=mean(VAZNULLRESULTS[[i]]$WNODF)
    VAZNULL_INDICES$WNODF_sd[i]=sd(VAZNULLRESULTS[[i]]$WNODF)
    VAZNULL_INDICES$connectance_mean[i]=mean(VAZNULLRESULTS[[i]]$connectance)
    VAZNULL_INDICES$connectance_sd[i]=sd(VAZNULLRESULTS[[i]]$connectance)
    VAZNULL_INDICES$compartments_mean[i]=mean(VAZNULLRESULTS[[i]]$compartments)
    VAZNULL_INDICES$compartments_sd[i]=sd(VAZNULLRESULTS[[i]]$compartments)
  }
  VAZNULL_INDICES$WMOD_z=(NET_INDICES$WMOD-VAZNULL_INDICES$WMOD_mean)/VAZNULL_INDICES$WMOD_sd
  VAZNULL_INDICES$WNODF_z=(NET_INDICES$WNODF-VAZNULL_INDICES$WNODF_mean)/VAZNULL_INDICES$WNODF_sd
  VAZNULL_INDICES$WNODA_z=(NET_INDICES$WNODA-VAZNULL_INDICES$WNODA_mean)/VAZNULL_INDICES$WNODA_sd
  save(VAZNULL_INDICES,VAZNULLMATRICES,VAZNULLRESULTS, file="vaznull.RData")
}
if(!run_vaznull){load("vaznull.RData")}
# rm so far #
rm(BETA0, binmatrix, NODF,wmatrix,WNODA,WNODF,i)
## Restricted null ##
if(run_restnull){
  source("Restricted-Null-Model/PosteriorProb.R")
  source("Restricted-Null-Model/RestNullModel.R")
  RESTNULLMATRICES=list()
  for (i in 1:7){
    MProb=PosteriorProb(M=as.matrix(matrices[[i]]), R.partitions =MOD[[i]]$Row_labels, C.partitions = MOD[[i]]$Col_labels,Prior.Pij = "degreeprob", Conditional.level ="areas")
    RESTNULLMATRICES[[i]]=RestNullModel(M=as.matrix(matrices[[i]]),Pij.Prob = MProb,allow.degeneration = F,Numbernulls = 1000,connectance = T,byarea = T,R.partitions = MOD[[i]]$Row_labels,C.partitions = MOD[[i]]$Col_labels)
    print(i)
  }
  RESTNULLRESULTS=list()
  RESTNULL_INDICES=data.frame(row.names = names(matrices))
  for (i in 1:7){
    RESTNULLRESULTS[[i]]=list()
    names(RESTNULLRESULTS)[i]=names(matrices[i])
    print(paste("Starting null model:", i, names(matrices[i])))
    for(j in 1:1000){
      if(j%%10==0){print(paste("null matrix n:",j))}
      NULLMAT=RESTNULLMATRICES[[i]][[j]]$NullMatrix
      # 1- nestedness SM
      # weighted
      WNODA=nest.smdm(NULLMAT, weighted = T,decreasing = "abund", constraints = c(MOD[[i]]$Row_labels,MOD[[i]]$Col_labels))
      WNODF=nest.smdm(NULLMAT, weighted = T,decreasing = "fill", constraints = c(MOD[[i]]$Row_labels,MOD[[i]]$Col_labels))
      RESTNULLRESULTS[[i]]$WNODA_SM[j]=WNODA$WNODA_SM_matrix
      RESTNULLRESULTS[[i]]$WNODF_SM[j]=WNODF$WNODF_SM_matrix
      RESTNULLRESULTS[[i]]$WNODA_DM[j]=WNODA$WNODA_DM_matrix
      RESTNULLRESULTS[[i]]$WNODF_DM[j]=WNODF$WNODF_DM_matrix
    }
    # Summary of null models
    RESTNULL_INDICES$WNODA_SM_mean[i]=mean(RESTNULLRESULTS[[i]]$WNODA_SM)
    RESTNULL_INDICES$WNODA_SM_sd[i]=sd(RESTNULLRESULTS[[i]]$WNODA_SM)
    RESTNULL_INDICES$WNODF_SM_mean[i]=mean(RESTNULLRESULTS[[i]]$WNODF_SM)
    RESTNULL_INDICES$WNODF_SM_sd[i]=sd(RESTNULLRESULTS[[i]]$WNODF_SM)
    RESTNULL_INDICES$WNODA_DM_mean[i]=mean(RESTNULLRESULTS[[i]]$WNODA_DM)
    RESTNULL_INDICES$WNODA_DM_sd[i]=sd(RESTNULLRESULTS[[i]]$WNODA_DM)
    RESTNULL_INDICES$WNODF_DM_mean[i]=mean(RESTNULLRESULTS[[i]]$WNODF_DM)
    RESTNULL_INDICES$WNODF_DM_sd[i]=sd(RESTNULLRESULTS[[i]]$WNODF_DM)
  }
  RESTNULL_INDICES$WNODA_SM_z=(NET_INDICES$WNODA_SM-RESTNULL_INDICES$WNODA_SM_mean)/RESTNULL_INDICES$WNODA_SM_sd
  RESTNULL_INDICES$WNODF_SM_z=(NET_INDICES$WNODF_SM-RESTNULL_INDICES$WNODF_SM_mean)/RESTNULL_INDICES$WNODF_SM_sd
  RESTNULL_INDICES$WNODA_DM_z=(NET_INDICES$WNODA_DM-RESTNULL_INDICES$WNODA_DM_mean)/RESTNULL_INDICES$WNODA_DM_sd
  RESTNULL_INDICES$WNODF_DM_z=(NET_INDICES$WNODF_DM-RESTNULL_INDICES$WNODF_DM_mean)/RESTNULL_INDICES$WNODF_DM_sd
  save(RESTNULL_INDICES,RESTNULLMATRICES,RESTNULLRESULTS, file="restnull.RData")
}
if(!run_restnull){load("restnull.RData")}