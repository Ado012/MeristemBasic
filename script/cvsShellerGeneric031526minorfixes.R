library(reshape2)
library(ggplot2)
library(dplyr)

setwd("/home/ace/Documents/Software_Projects/MeristemBasic/Results/Processed_DataFiles")


#cbind.fill replacement...should check 022021
cbindPad <- function(...){
  args <- list(...) #take all the arguments received and put them in a list
  n <- sapply(args,nrow) #find the number of rows for each argument
  mx <- max(n) #get the max number of rows
  pad <- function(x, mx){
    if (nrow(x) < mx){
      nms <- colnames(x)
      padTemp <- matrix(NA, mx - nrow(x), ncol(x))
      colnames(padTemp) <- nms
      if (ncol(x)==0) {
        return(padTemp)
      } else {
        return(rbind(x,padTemp))
      }
    }
    else{
      return(x)
    }
  }
  rs <- lapply(args,pad,mx)
  return(do.call(cbind,rs))
}




distanceFinder<-function(point1,point2)
{
  
  
  
  distancex=as.numeric(point1[1])-as.numeric(point2[1])
  
  
  distancey=as.numeric(point1[2])-as.numeric(point2[2])
  distancez=as.numeric(point1[3])-as.numeric(point2[3])
  
  distance=sqrt(distancex*distancex+distancey*distancey+distancez*distancez)
  
  
}

#112525drop the species concentration > 0 requirement. Also rename clv3value to specievalue to make it clearer. Hopefully will fix issue with Nas used to fill 
#missing values by increasing available data points
sheller<-function(basePoint, result, boundaries, targetColumn)
{
  
  filteredset<-c()
  
  for (i in 1:nrow(result))
  {
    #get coordinates of results
    xpoint<-as.double(result[i,][1])
    ypoint<-as.double(result[i,][2])
    zpoint<-as.double(result[i,][3])
    
    #calculate layer and distance from center
    testPoint<-c(xpoint,ypoint,zpoint)
    axialPoint<-c(0,0,zpoint)
    
    distanceCalc<-distanceFinder(point1=testPoint,point2=basePoint)
    axialDistance<-distanceFinder(point1=testPoint,point2=axialPoint)
  
      speciesValue=result[,targetColumn][i]
    
      #filter points to retrieve those in the proper boundaries
    if (distanceCalc > boundaries[1] && distanceCalc < boundaries[2] &&  axialDistance < boundaries[3])
    {
      filteredset<-rbind(filteredset,result[i,])
      
    }
    
    
    
    
    
    
  }
  
  #filtereddf<-data.frame(filteredset)
  #cut results to a consistent size of the 20 highest values. Supposedly will still work if there are less than 20 rows
  #filteredsetcuttosize<-head(sort(filteredset[,targetColumn],decreasing=TRUE),n=20)
  #per auto suggest replaced targetcolumn with .data[[targetcolumn]] to properly arrange even though it doesn't matter 031526
  filteredsetcuttosize<-arrange(filteredset,desc(.data[[targetColumn]]))
  filteredsetcuttosize <-head(filteredsetcuttosize, n= 20)
  
  #colnames(filteredset)<-c("x","y","z","radius","WUSRNA","WUSNuc_WS","WUSCyto","CLV3Sig1","CLV3_Peptide", "StochasticTimeOverFlow","CLV3Sig2", "MarkerOverFlow", "Monomer", "Dimer","MonomerMarker","DimerMarker", "neighbors?")
  #pass results
  filteredsetcuttosize
  
}


csvSheller<-function(inputfile,resultspath, targetColumn)
{
  
  meristemdata<-read.csv(inputfile,sep=' ',header=T)
  
  
  #Layers are reckoned from the bottom center of the meristem
  base<-c(0,0,0)
  filteredset<-c()
  
  
  
  
  #L1
  layerboundaries<-c(8.5,100,2.5)
  filteredL1<-sheller(basePoint=base, result=meristemdata, boundaries=layerboundaries, targetColumn=targetColumn)
  filteredL1Values<-data.frame(filteredL1[,targetColumn])
  
  if (nrow(filteredL1Values) < 1)
    filteredL1Values<-data.frame(matrix(NA, nrow = 20, ncol = 1))
  
  #change 8 to 7.5 for sim 
  #L2
  layerboundaries<-c(8,8.5,2.5)
  filteredL2<-sheller(basePoint=base,result=meristemdata, boundaries=layerboundaries, targetColumn=targetColumn)
  filteredL2Values<-data.frame(filteredL2[,targetColumn])
  
  if (nrow(filteredL2Values) < 1)
    filteredL2Values<-data.frame(matrix(NA, nrow = 20, ncol = 1))
  
  #L3
  layerboundaries<-c(6.5,8,2.5)
  filteredL3<-sheller(basePoint=base,result=meristemdata, boundaries=layerboundaries, targetColumn=targetColumn)
  filteredL3Values<-data.frame(filteredL3[,targetColumn])
  
  if (nrow(filteredL3Values) < 1)
    filteredL3Values<-data.frame(matrix(NA, nrow = 20, ncol = 1))
  
  #L4
  layerboundaries<-c(6,7,2.5)
  filteredL4<-sheller(basePoint=base,result=meristemdata, boundaries=layerboundaries, targetColumn=targetColumn)
  filteredL4Values<-data.frame(filteredL4[,targetColumn])
  
  if (nrow(filteredL4Values) < 1)
    filteredL4Values<-data.frame(matrix(NA, nrow = 20, ncol = 1))
  
  #L5
  layerboundaries<-c(5,6,2.5)
  filteredL5<-sheller(basePoint=base,result=meristemdata, boundaries=layerboundaries, targetColumn=targetColumn)
  filteredL5Values<-data.frame(filteredL5[,targetColumn])
  
  if (nrow(filteredL5Values) < 1)
    filteredL5Values<-data.frame(matrix(NA, nrow = 20, ncol = 1))
  
  #L6
  layerboundaries<-c(4,5,2.5)
  filteredL6<-sheller(basePoint=base,result=meristemdata, boundaries=layerboundaries, targetColumn=targetColumn)
  filteredL6Values<-data.frame(filteredL6[,targetColumn])
  
  if (nrow(filteredL6Values) < 1)
    filteredL6Values<-data.frame(matrix(NA, nrow = 20, ncol = 1))
  
  
  #L7
  layerboundaries<-c(3,4,2.5)
  filteredL7<-sheller(basePoint=base,result=meristemdata, boundaries=layerboundaries, targetColumn=targetColumn)
  filteredL7Values<-data.frame(filteredL7[,targetColumn])
  
  if (nrow(filteredL7Values) < 1)
    filteredL7Values<-data.frame(matrix(NA, nrow = 20, ncol = 1))
  
  
  #mergedLayerData<-cbind.fill(mergedLayerData,c(filteredL1CLV3Values,filteredL2CLV3Values,filteredL3CLV3Values,filteredL4CLV3Values,
                                                #filteredL5CLV3Values,filteredL6CLV3Values,filteredL7CLV3Values), fill =NA)

  
  mergedLayerData<-cbindPad(filteredL1Values,filteredL2Values)
  mergedLayerData<-cbindPad(mergedLayerData,filteredL3Values)
  mergedLayerData<-cbindPad(mergedLayerData,filteredL4Values)
  mergedLayerData<-cbindPad(mergedLayerData,filteredL5Values)
  mergedLayerData<-cbindPad(mergedLayerData,filteredL6Values)
  mergedLayerData<-cbindPad(mergedLayerData,filteredL7Values)
  
  #mergedLayerData<-cbind.fill(mergedLayerData,c(filteredL1CLV3Values,filteredL2CLV3Values,filteredL3CLV3Values,filteredL4CLV3Values,
  #filteredL5CLV3Values,filteredL6CLV3Values,filteredL7CLV3Values), fill =NA)
  
  df <- data.frame(matrix(ncol = 2, nrow = 20))
  #added in extra columns if short may want to revise later 022321
  if (ncol(mergedLayerData)==5)
    mergedLayerData<-cbind(mergedLayerData,df)
  
  if (nrow(mergedLayerData) < 20) #addition from cvsshellergeneric2
  {
    rowstobeadded = 20-nrow(mergedLayerData)
    
    df <- data.frame(matrix(ncol = 7, nrow = rowstobeadded))
    
    colnames(df) <- colnames(mergedLayerData)
    
    mergedLayerData<-rbind(mergedLayerData,df)
    
    
  }
  
  constructedname<-paste(inputfile,"filtered",sep="")
  constructednamefinal<-paste(constructedname,".csv",sep="")
  
  destination<-paste(resultspath,"/",constructednamefinal,sep='')
  
 columnnames <-colnames(mergedLayerData)
 
 #added on 041124 to systemize column namessince a zeroed out dataframe will not merge
 #changed on 031426 to multiple strings per auto sugesstion even though it doesn't matter
 names(mergedLayerData) <- c("layer1", "layer2", "layer3", "layer4", "layer5", "layer6", "layer7")
  
  write.table(mergedLayerData, file=destination, row.names=FALSE, quote=F)

  mergedLayerData
  
  
}




