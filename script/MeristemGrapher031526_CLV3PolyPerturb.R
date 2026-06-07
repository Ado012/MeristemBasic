library(Biostrings)



#test to make sure all points are shown on the graph

LimitsTest<-function(datavalues, ymin, ymax, dataname)
{
  
  limitsexceeded = 0
  
  dataframesize = nrow(datavalues)
  
  for (i in 1:dataframesize)
  {
    if ((datavalues[,'value'][i] > ymax) || (datavalues[,'value'][i] < ymin))
    {
      limitsexceeded = dataname
      break
    }
  }
  
  limitsexceeded
}





MeristemGrapher<-function(graphName,graphTitle,graphXtitle,graphYtitle, graphYmin, graphYlimit, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
{

mergedDynamicData<-do.call("cbind", list(fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data))




coltitles<-c("FM1L1","FM1L2","FM1L3","FM1L4","FM1L5","FM1L6","FM1L7",
             "FM2L1","FM2L2","FM2L3","FM2L4","FM2L5","FM2L6","FM2L7",
             "FM3L1","FM3L2","FM3L3","FM3L4","FM3L5","FM3L6","FM3L7",
             "FM4L1","FM4L2","FM4L3","FM4L4","FM4L5","FM4L6","FM4L7",
             "FM5L1","FM5L2","FM5L3","FM5L4","FM5L5","FM5L6","FM5L7",
             "FM6L1","FM6L2","FM6L3","FM6L4","FM6L5","FM6L6","FM6L7",
             "FM7L1","FM7L2","FM7L3","FM7L4","FM7L5","FM7L6","FM7L7",
             "FM8L1","FM8L2","FM8L3","FM8L4","FM8L5","FM8L6","FM8L7",
             "FM9L1","FM9L2","FM9L3","FM9L4","FM9L5","FM9L6","FM9L7",
             "FM10L1","FM10L2","FM10L3","FM10L4","FM10L5","FM10L6","FM10L7"
)


colnames(mergedDynamicData)<-coltitles



mergedDynamicData[is.na(mergedDynamicData)] = 0





markerexpgenotype<-c(rep("1",5),rep("2",5),rep("3",5),rep("4",5),rep("5",5),rep("6",5),rep("7",5),rep("8",5),rep("9",5),rep("10",5))
simgenotype<-c(rep("1",140),rep("2",140),rep("3",140),rep("4",140),rep("5",140),rep("6",140),rep("7",140),rep("8",140),rep("9",140),rep("10",140))

simgenotype<-melt(simgenotype)
colnames(simgenotype)<-c("simgenotype")
mergedDynamicData<-melt(mergedDynamicData, id.vars = NULL) #had to add id.vars due to dimer measure not working
mergedDynamicData<-cbind(mergedDynamicData,simgenotype)


#test to make sure graph does not have points off screen
limittest <- LimitsTest(mergedDynamicData, 0, graphYlimit, graphName)

#replace 0 with graphYmin 031526
resultGraph <- ggplot(mergedDynamicData, aes(x=variable, y=value, colour=simgenotype)) + ylab(graphYtitle) + xlab("genotype")+ ggtitle(graphTitle)  +
  theme(plot.title = element_text(hjust = 0.5))+ geom_jitter(shape=16, position=position_jitter(0.2))+coord_cartesian(ylim=c(graphYmin,graphYlimit)) + 
  scale_color_manual(values=c("purple","blue",'green4','orange','slateblue','maroon','coral' ,'tomato','deeppink','tan3')) +
  scale_fill_manual(values=c("purple","blue",'green4','orange','slateblue','maroon','coral','tomato','deeppink','tan3'))



ggsave(filename = graphName,
       width = 3700,
       height = 600,
       plot = resultGraph,
       dpi = 100,
       units = "px",
       device = 'png')


limittest
}

MeristemGrapherCompactGraph<-function(graphName,graphTitle,graphXtitle,graphYtitle, graphYmin, graphYlimit, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
{
  
  fWT_Data = fM6_Data
  fM4_Data = fM7_Data
  f9M_Data = fM8_Data
  fTM_Data = fM9_Data
  f9I_Data = fM10_Data
  
  
  
  
  #added 12/10/20
  mergedDynamicData<-do.call("rbind", list(fWT_Data,fM4_Data,f9M_Data,fTM_Data,f9I_Data))
  #WT 970M CLV3-2 970M4 970I
  
  
  
  #added 12/10/20
  colnames(mergedDynamicData)<-c("L1","L2","L3","L4","L5","L6","L7")
  
  
  
  mergedDynamicData[is.na(mergedDynamicData)] = 0
  
  
  
  
  
  
  simgenotype<-c(rep("FWTL1",20),
                 rep("FWTL2",20),
                 rep("FWTL3",20),
                 rep("FWTL4",20),
                 rep("FWTL5",20),
                 rep("FWTL6",20),
                 rep("FWTL7",20))
  
  simgenotype2<-c(rep("1. wildtype",20),
                  rep("1. wildtype",20),
                  rep("1. wildtype",20),
                  rep("1. wildtype",20),
                  rep("1. wildtype",20),
                  rep("1. wildtype",20),
                  rep("1. wildtype",20))
  
  
  
  simgenotype<-melt(simgenotype)
  colnames(simgenotype)<-c("layer")
  mergedDynamicData<-melt(mergedDynamicData)
  mergedDynamicData<-cbind(mergedDynamicData,simgenotype)
  
  
  #addition 121020 
  simgenotype2<-melt(simgenotype2)
  colnames(simgenotype2)<-c("simgen")
  mergedDynamicData<-cbind(mergedDynamicData,simgenotype2)
  
  
  fWTMean_Data<-colMeans(na.rm=TRUE,fWT_Data)

  
  fWTSD_Data<-sapply(fWT_Data,sd, na.rm= TRUE)


  fWTMean_Data[is.nan(fWTMean_Data)]<-0


  fWTSD_Data[is.na(fWTSD_Data)]<-0


  fWTSD_Data <- data.frame(matrix(unlist(fWTSD_Data), nrow=length(fWTSD_Data), byrow=TRUE))


  colnames(fWTSD_Data) <- 'standdev'


  
  fWTMean_Data <- data.frame(matrix(unlist(fWTMean_Data), nrow=length(fWTMean_Data), byrow=TRUE))


  
  
  colnames(fWTMean_Data) <- 'mean'


  
  fWTMean_Data<-cbind(fWTMean_Data,fWTSD_Data)


  
  
  #added 12/10/20
  mergedDynamicDataMean<-do.call("rbind", list(fWTMean_Data))
  
  
  #added 12/10/20
  #colnames(mergedDynamicDataMean)<-c("L1","L2","L3","L4","L5","L6","L7")
  
  layerIDMean<-c(rep(c("L1","L2","L3","L4","L5","L6","L7"),1))
  
  
  
  layerGenotypeIDMean<-c(rep("FWTL1",1),rep("FWTL2",1),rep("FWTL3",1),rep("FWTL4",1),rep("FWTL5",1),rep("FWTL6",1),rep("FWTL7",1))
  
  
  
  genotypeIDMean<-c(rep("1. wildtype",1),rep("1. wildtype",1),rep("1. wildtype",1),rep("1. wildtype",1),rep("1. wildtype",1),rep("1. wildtype",1),rep("1. wildtype",1))
  
  
  
  
  layerIDMean<-melt(layerIDMean)
  colnames(layerIDMean)<-c("layer")
  #mergedDynamicDataMean<-melt(mergedDynamicDataMean)
  mergedDynamicDataMean<-cbind(mergedDynamicDataMean,layerIDMean)
  
  layerGenotypeIDMean<-melt(layerGenotypeIDMean)
  colnames(layerGenotypeIDMean)<-c("layerGenotype")
  mergedDynamicDataMean<-cbind(mergedDynamicDataMean,layerGenotypeIDMean)
  
  #addition 121020 
  genotypeIDMean<-melt(genotypeIDMean)
  colnames(genotypeIDMean)<-c("SimGenotype")
  mergedDynamicDataMean<-cbind(mergedDynamicDataMean,genotypeIDMean)
  #mergedDynamicDataMean$value=as.numeric(mergedDynamicDataMean$value)
  
  #lff<-as.numeric(mergedDynamicDataMean$value)
  
  
  
  
  
  resultGraph <-ggplot(data = mergedDynamicDataMean) + theme(plot.title = element_text(hjust = 0.5)) +
    geom_line(aes(x=layer, y=mean,group=SimGenotype, colour=SimGenotype, linewidth = 2)) +
    geom_pointrange(aes(x=layer, y=mean,ymin=mean-standdev, ymax=mean+standdev, colour=SimGenotype)) +
    ylim(graphYmin, graphYlimit)+
    labs(y= graphYtitle, x = "Layer")+ ggtitle(graphTitle) +
    theme(plot.title = element_text(size=18),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.x = element_text(size=11),
          axis.text.y = element_text(size=11),                         
          legend.title=element_text(size=16),
          legend.text=element_text(size=12)
    ) +
    scale_color_manual(values=c("purple","blue",'green','orange','slateblue','maroon','coral')) +
    scale_fill_manual(values=c("purple","blue",'green','orange','slateblue','maroon','coral'))
  
  
  ggsave(filename = graphName,
         width = 1000,
         height = 600,
         plot = resultGraph,
         dpi = 100,
         units = "px",
         device = 'png')
  
  #code for line or ribbon standard dev
  #geom_pointrange(aes(ymin=collatedMean-collatedStdDev, ymax=collatedMean+collatedStdDev))
  #geom_ribbon(data = mergedDynamicDataMean, aes(x=layer, y=mean, group= SimGenotype, fill = SimGenotype, ymin=mean+standdev, ymax=mean-standdev),  alpha=0.25)
  
}





MeristemGrapherCollated<-function(filetitle, graphtitle,DM1,DM2,DM3,DM4,DM5, DM6,DM7,DM8,DM9, DM0)
{
  limittestrecord <-c()
  
  
  fM1_Data<-csvSheller(DM1,"filteredfiles", targetColumn = "WUSNuc_WS")
  fM2_Data<-csvSheller(DM2,"filteredfiles", targetColumn = "WUSNuc_WS")
  fM3_Data<-csvSheller(DM3,"filteredfiles", targetColumn = "WUSNuc_WS")
  fM4_Data<-csvSheller(DM4,"filteredfiles", targetColumn = "WUSNuc_WS")
  fM5_Data<-csvSheller(DM5,"filteredfiles", targetColumn = "WUSNuc_WS")
  fM6_Data<-csvSheller(DM6,"filteredfiles", targetColumn = "WUSNuc_WS")
  fM7_Data<-csvSheller(DM7,"filteredfiles", targetColumn = "WUSNuc_WS")
  fM8_Data<-csvSheller(DM8,"filteredfiles", targetColumn = "WUSNuc_WS")
  fM9_Data<-csvSheller(DM9,"filteredfiles", targetColumn = "WUSNuc_WS")
  fM10_Data<-csvSheller(DM0,"filteredfiles", targetColumn = "WUSNuc_WS")
  
  
  filetitlefinal <-paste("WUSNuc_",filetitle)
  graphtitlefinal <-paste("WUSNuc_",graphtitle)
  
  limittest <-MeristemGrapher(filetitlefinal,graphtitlefinal,"placeholder","wusNuc level", -5, 150, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  limittestrecord <-append(limittestrecord,limittest)
  
  
  
  fM1_Data<-csvSheller(DM1,"filteredfiles", targetColumn = "WUSCyto")
  fM2_Data<-csvSheller(DM2,"filteredfiles", targetColumn = "WUSCyto")
  fM3_Data<-csvSheller(DM3,"filteredfiles", targetColumn = "WUSCyto")
  fM4_Data<-csvSheller(DM4,"filteredfiles", targetColumn = "WUSCyto")
  fM5_Data<-csvSheller(DM5,"filteredfiles", targetColumn = "WUSCyto")
  fM6_Data<-csvSheller(DM6,"filteredfiles", targetColumn = "WUSCyto")
  fM7_Data<-csvSheller(DM7,"filteredfiles", targetColumn = "WUSCyto")
  fM8_Data<-csvSheller(DM8,"filteredfiles", targetColumn = "WUSCyto")
  fM9_Data<-csvSheller(DM9,"filteredfiles", targetColumn = "WUSCyto")
  fM10_Data<-csvSheller(DM0,"filteredfiles", targetColumn = "WUSCyto")
  
  
  filetitlefinal <-paste("WUSCyto_",filetitle)
  graphtitlefinal <-paste("WUSCyto_",graphtitle)
  
  limittest <-MeristemGrapher(filetitlefinal,graphtitlefinal,"placeholder","wusCyto level", 0, 65, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  limittestrecord <-append(limittestrecord,limittest)
  
  
  
  
  
  
  
  fM1_Data<-csvSheller(DM1,"filteredfiles", targetColumn = "WUSRNA")
  fM2_Data<-csvSheller(DM2,"filteredfiles", targetColumn = "WUSRNA")
  fM3_Data<-csvSheller(DM3,"filteredfiles", targetColumn = "WUSRNA")
  fM4_Data<-csvSheller(DM4,"filteredfiles", targetColumn = "WUSRNA")
  fM5_Data<-csvSheller(DM5,"filteredfiles", targetColumn = "WUSRNA")
  fM6_Data<-csvSheller(DM6,"filteredfiles", targetColumn = "WUSRNA")
  fM7_Data<-csvSheller(DM7,"filteredfiles", targetColumn = "WUSRNA")
  fM8_Data<-csvSheller(DM8,"filteredfiles", targetColumn = "WUSRNA")
  fM9_Data<-csvSheller(DM9,"filteredfiles", targetColumn = "WUSRNA")
  fM10_Data<-csvSheller(DM0,"filteredfiles", targetColumn = "WUSRNA")
  
  
  filetitlefinal <-paste("WUSRNA_",filetitle)
  graphtitlefinal <-paste("WUSRNA_",graphtitle)
  
  limittest <- MeristemGrapher(filetitlefinal,graphtitlefinal,"placeholder","wusRNA level", -20, 400, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  limittestrecord <-append(limittestrecord,limittest)
  
  
  
  
  
  fM1_Data<-csvSheller(DM1,"filteredfiles", targetColumn = "CLV3Sig1")
  fM2_Data<-csvSheller(DM2,"filteredfiles", targetColumn = "CLV3Sig1")
  fM3_Data<-csvSheller(DM3,"filteredfiles", targetColumn = "CLV3Sig1")
  fM4_Data<-csvSheller(DM4,"filteredfiles", targetColumn = "CLV3Sig1")
  fM5_Data<-csvSheller(DM5,"filteredfiles", targetColumn = "CLV3Sig1")
  fM6_Data<-csvSheller(DM6,"filteredfiles", targetColumn = "CLV3Sig1")
  fM7_Data<-csvSheller(DM7,"filteredfiles", targetColumn = "CLV3Sig1")
  fM8_Data<-csvSheller(DM8,"filteredfiles", targetColumn = "CLV3Sig1")
  fM9_Data<-csvSheller(DM9,"filteredfiles", targetColumn = "CLV3Sig1")
  fM10_Data<-csvSheller(DM0,"filteredfiles", targetColumn = "CLV3Sig1")
  
  
  filetitlefinal <-paste("CLV3Sig_",filetitle)
  graphtitlefinal <-paste("CLV3Sig_",graphtitle)
  
  limittest <- MeristemGrapher(filetitlefinal,graphtitlefinal,"placeholder","CLV3 RNA level", 0, 17, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  limittestrecord <-append(limittestrecord,limittest)
  
  
  
  fM1_Data<-csvSheller(DM1,"filteredfiles", targetColumn = "CLV3_Peptide")
  fM2_Data<-csvSheller(DM2,"filteredfiles", targetColumn = "CLV3_Peptide")
  fM3_Data<-csvSheller(DM3,"filteredfiles", targetColumn = "CLV3_Peptide")
  fM4_Data<-csvSheller(DM4,"filteredfiles", targetColumn = "CLV3_Peptide")
  fM5_Data<-csvSheller(DM5,"filteredfiles", targetColumn = "CLV3_Peptide")
  fM6_Data<-csvSheller(DM6,"filteredfiles", targetColumn = "CLV3_Peptide")
  fM7_Data<-csvSheller(DM7,"filteredfiles", targetColumn = "CLV3_Peptide")
  fM8_Data<-csvSheller(DM8,"filteredfiles", targetColumn = "CLV3_Peptide")
  fM9_Data<-csvSheller(DM9,"filteredfiles", targetColumn = "CLV3_Peptide")
  fM10_Data<-csvSheller(DM0,"filteredfiles", targetColumn = "CLV3_Peptide")
  
  filetitlefinal <-paste("CLV3Peptide_",filetitle)
  graphtitlefinal <-paste("CLV3Peptide_",graphtitle)
  
  limittest <- MeristemGrapher(filetitlefinal,graphtitlefinal,"placeholder","CLV3Peptide level", 0, 250, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  limittestrecord <-append(limittestrecord,limittest)
  
  
  
  fM1_Data<-csvSheller(DM1,"filteredfiles", targetColumn = "Monomer")
  fM2_Data<-csvSheller(DM2,"filteredfiles", targetColumn = "Monomer")
  fM3_Data<-csvSheller(DM3,"filteredfiles", targetColumn = "Monomer")
  fM4_Data<-csvSheller(DM4,"filteredfiles", targetColumn = "Monomer")
  fM5_Data<-csvSheller(DM5,"filteredfiles", targetColumn = "Monomer")
  fM6_Data<-csvSheller(DM6,"filteredfiles", targetColumn = "Monomer")
  fM7_Data<-csvSheller(DM7,"filteredfiles", targetColumn = "Monomer")
  fM8_Data<-csvSheller(DM8,"filteredfiles", targetColumn = "Monomer")
  fM9_Data<-csvSheller(DM9,"filteredfiles", targetColumn = "Monomer")
  fM10_Data<-csvSheller(DM0,"filteredfiles", targetColumn = "Monomer")
  
  
  filetitlefinal <-paste("Monomer_",filetitle)
  graphtitlefinal <-paste("Monomer_",graphtitle)
  
  limittest <- MeristemGrapher(filetitlefinal,graphtitlefinal,"placeholder","Monomer level", 0, 15, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  limittestrecord <-append(limittestrecord,limittest)
  
  
  
  fM1_Data<-csvSheller(DM1,"filteredfiles", targetColumn = "Dimer")
  fM2_Data<-csvSheller(DM2,"filteredfiles", targetColumn = "Dimer")
  fM3_Data<-csvSheller(DM3,"filteredfiles", targetColumn = "Dimer")
  fM4_Data<-csvSheller(DM4,"filteredfiles", targetColumn = "Dimer")
  fM5_Data<-csvSheller(DM5,"filteredfiles", targetColumn = "Dimer")
  fM6_Data<-csvSheller(DM6,"filteredfiles", targetColumn = "Dimer")
  fM7_Data<-csvSheller(DM7,"filteredfiles", targetColumn = "Dimer")
  fM8_Data<-csvSheller(DM8,"filteredfiles", targetColumn = "Dimer")
  fM9_Data<-csvSheller(DM9,"filteredfiles", targetColumn = "Dimer")
  fM10_Data<-csvSheller(DM0,"filteredfiles", targetColumn = "Dimer")
  
  
  
  filetitlefinal <-paste("Dimer_",filetitle)
  graphtitlefinal <-paste("Dimer_",graphtitle)
  
  limittest <- MeristemGrapher(filetitlefinal,graphtitlefinal,"placeholder","Dimer level", 0, 15, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  limittestrecord <-append(limittestrecord,limittest)
  
  
  print(limittestrecord)
}






