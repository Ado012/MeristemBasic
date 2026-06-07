#csvsplitter
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

csvSplitter<-function(inputfile,resultspath)
{

result<-read.csv(inputfile,sep=' ',header=F)

resultdfInitial<-data.frame(result)

rowsizeInitial<-nrow(resultdfInitial)

resultdf<-result[3:rowsizeInitial,] #subset out data

dir.create(resultspath)

k=1
g=1

while (g < 100)
{  
  
    
    print(k+1)
    print(k+1365)
  
resultpiece<-resultdf[(k):(k+1365),]

iteratorchar<-as.character(g)

if(g<10)
  iteratorchar<-paste("00",iteratorchar,sep="")

else if(g==10 || g<100)
  iteratorchar<-paste("0",iteratorchar,sep="")
  
if (g==10)
print(resultpiece)
  
  
constructedname<-paste("pararesults",iteratorchar,sep="")
constructednamefinal<-paste(constructedname,".csv",sep="")

destination<-paste(resultspath,"/",constructednamefinal,sep='')

write.table(resultpiece, file=destination, row.names=FALSE, quote=F, col.names = c("x","y","z","radius","WUSRNA","WUSNuc_WS","WUSCyto","CLV3Sig1","CLV3_Peptide","StochasticTimeOverFlow","CLV3Sig2", "MarkerOverFlow", "Monomer", "Dimer","MonomerMarker","DimerMarker", 
                                                                                   "ckReceptor","ckLigand","ckComplex","WExport", "neighbors?"))

#henrik annotations
#write.table(resultpiece, file=constructednamefinal, row.names=FALSE, quote=F, col.names = c("x","y","z","WUSRNA","WUSNuc_WS","CLV3MRNA","Y","Kanadi","KanadiSig","A","B","AG1","AS1","AG2","AS2a","AS2b","var1","var2"))


k=k+1366+1
g=g+1
}

}

csvSplitter(args[1],args[2])
