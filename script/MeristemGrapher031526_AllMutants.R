library(Biostrings)

#adapt LimitsTest to AllMutants
LimitsTest<-function(datavalues, ymin, ymax, dataname)
{
  
  limitsexceeded = 0
  
  dataframesize = nrow(datavalues)
  
  for (i in 1:dataframesize)
  {
    valuemax = datavalues[,'mean'][i] + datavalues[,'standdev'][i] 
    valuemin = datavalues[,'mean'][i] - datavalues[,'standdev'][i] 
    
    
    if (( valuemax > ymax) || (valuemin < ymin))
    {
      limitsexceeded = dataname
      break
    }
  }
  
  limitsexceeded
}




MutantComponentGraph<-function(graphName,graphTitle,graphXtitle,graphYtitle, graphYmin, graphYlimit, fMut_Data, layerprefix, mutantname)
{
  
  #added 12/10/20
  #mergedDynamicData<-do.call("rbind", list(fWT_Data,fM4_Data,f9M_Data,fTM_Data,f9I_Data))
  #WT 970M CLV3-2 970M4 970I
  
  fMutMean_Data<-colMeans(na.rm=TRUE,fMut_Data)
  
  
  fMutSD_Data<-sapply(fMut_Data,sd, na.rm= TRUE)
  
  
  if (any(is.na(fMutMean_Data)))
  {print("NAs in 2nd dataframe")
    print(fMutMean_Data)}
  
  
  fMutMean_Data[is.nan(fMutMean_Data)]<-0
  
  
  if (any(is.na(fMutSD_Data)))
  {print("NAs in 3nd dataframe")
    print(fMutSD_Data)}
  
  
  fMutSD_Data[is.na(fMutSD_Data)]<-0
  
  
  fMutSD_Data <- data.frame(matrix(unlist(fMutSD_Data), nrow=length(fMutSD_Data), byrow=TRUE))
  
  
  colnames(fMutSD_Data) <- 'standdev'
  
  
  
  fMutMean_Data <- data.frame(matrix(unlist(fMutMean_Data), nrow=length(fMutMean_Data), byrow=TRUE))
  
  
  
  
  colnames(fMutMean_Data) <- 'mean'
  
  
  
  fMutMean_Data<-cbind(fMutMean_Data,fMutSD_Data)
  
  
  
  
  #added 12/10/20
  mergedDynamicDataMean<-do.call("rbind", list(fMutMean_Data))
  
  
  #added 12/10/20
  #colnames(mergedDynamicDataMean)<-c("L1","L2","L3","L4","L5","L6","L7")
  
  layerIDMean<-c(rep(c("L1","L2","L3","L4","L5","L6","L7"),1))
  
  
  
  layerGenotypeIDMean<-c(rep(paste0(layerprefix,"L1"),1),rep(paste0(layerprefix,"L2"),1),rep(paste0(layerprefix,"L3"),1),rep(paste0(layerprefix,"L4"),1),
                         rep(paste0(layerprefix,"L5"),1),rep(paste0(layerprefix,"L6"),1),rep(paste0(layerprefix,"L7"),1))
  
  
  
  genotypeIDMean<-c(rep(mutantname,1),rep(mutantname,1),rep(mutantname,1),rep(mutantname,1),rep(mutantname,1),rep(mutantname,1),rep(mutantname,1))
  
  
  
  
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
  
  
  print('graphdata')
  print(mergedDynamicDataMean)
  
  limittest <- LimitsTest(mergedDynamicDataMean, graphYmin, graphYlimit, graphName)
  
  
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
  
  
  limittest
}





MeristemGrapherCompactGraph<-function(graphName,graphTitle,graphXtitle,graphYtitle, graphYmin, graphYlimit, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
{
  
  fWT_Data = fM6_Data
  fM4_Data = fM7_Data
  f9M_Data = fM8_Data
  fTM_Data = fM9_Data
  f9I_Data = fM10_Data
  
  limittestrecord <-c()
  
  
  graphNameIndividual = paste0("WT",graphName)
  graphTitleIndividual = paste0("WT",graphTitle)
  
  limittest <-MutantComponentGraph(graphNameIndividual,graphTitle,graphXtitle,graphYtitle, graphYmin, graphYlimit, fWT_Data, "FWT", "wildtype")
  
  limittestrecord <-append(limittestrecord,limittest)
  
  
  graphNameIndividual = paste0("M4",graphName)
  graphTitleIndividual = paste0("M4",graphTitle)
  
  limittest <-MutantComponentGraph(graphNameIndividual,graphTitleIndividual,graphXtitle,graphYtitle, graphYmin, graphYlimit, fM4_Data, "FM4", "970M4")
  
  limittestrecord <-append(limittestrecord,limittest)
  
  
  graphNameIndividual = paste0("9M",graphName)
  graphTitleIndividual = paste0("9M",graphTitle)
  
  limittest <-MutantComponentGraph(graphNameIndividual,graphTitleIndividual,graphXtitle,graphYtitle, graphYmin, graphYlimit, f9M_Data, "F9M", "970m")
  
  limittestrecord <-append(limittestrecord,limittest)
  
  
  
  graphNameIndividual = paste0("9I",graphName)
  graphTitleIndividual = paste0("9I",graphTitle)
  
  limittest <-MutantComponentGraph(graphNameIndividual,graphTitleIndividual,graphXtitle,graphYtitle, graphYmin, graphYlimit, f9I_Data, "F9I", "970i")
  
  limittestrecord <-append(limittestrecord,limittest)
  
  print(limittestrecord)
  
}





MeristemGrapherCollated<-function(filetitle, graphtitle,DM1,DM2,DM3,DM4,DM5, DM6,DM7,DM8,DM9, DM0)
{
  
  
  
  
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
  
  MeristemGrapherCompactGraph(filetitlefinal,graphtitlefinal,"placeholder","wusNuc level", -5, 150, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  
  

  
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
  
  MeristemGrapherCompactGraph(filetitlefinal,graphtitlefinal,"placeholder","wusCyto level", 0, 75, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  
  

  
  
  
  
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
  
  MeristemGrapherCompactGraph(filetitlefinal,graphtitlefinal,"placeholder","wusRNA level", -30, 400, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  

  
  
  
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
  
  MeristemGrapherCompactGraph(filetitlefinal,graphtitlefinal,"placeholder","CLV3 RNA level", -0.5, 6, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  
  

  
  
  
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
  
  MeristemGrapherCompactGraph(filetitlefinal,graphtitlefinal,"placeholder","CLV3Peptide level", 0, 400, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  
  

  
  
  
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
  
  MeristemGrapherCompactGraph(filetitlefinal,graphtitlefinal,"placeholder","Monomer level", -2, 15, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  

  
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
  
  MeristemGrapherCompactGraph(filetitlefinal,graphtitlefinal,"placeholder","Dimer level", -2, 15, fM1_Data,fM2_Data,fM3_Data,fM4_Data,fM5_Data,fM6_Data,fM7_Data,fM8_Data,fM9_Data,fM10_Data)
  

  
}





