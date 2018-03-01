load('CoGAPS.nP5.Rda')

library('org.Hs.eg.db')
library('NCBI2R')



simplicityGENES <- function(As, Ps) {
  # rescale p's to have max 1
  pscale <- apply(Ps,1,max)
  
  # rescale A in accordance with p's having max 1
  As <- sweep(As, 2, pscale, FUN="*")
  
  # find the A with the highest magnitude
  Arowmax <- t(apply(As, 1, function(x) x/max(x)))
  pmax <- apply(As, 1, max)
  
  # determine which genes are most associated with each pattern
  ssl <- matrix(NA, nrow=nrow(As), ncol=ncol(As),
                dimnames=dimnames(As))
  for (i in 1:ncol(As)) {
    lp <- rep(0, ncol(As))
    lp[i] <- 1
    ssl.stat <- apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp)))
    ssl[order(ssl.stat),i] <- 1:length(ssl.stat)
  }
  
  return(ssl)
  
}

z <- nP5$Amean/(nP5$Asd+.Machine$double.eps)
simGenes5 <- simplicityGENES(As = nP5$Amean, Ps = nP5$Pmean)

CG <- nP5

library('Hmisc')

for (p in 1:ncol(simGenes5)) {
  pdf(sprintf('SimplicityGenes%d.pdf',p))
  sortSim <- names(sort(simGenes5[,p],decreasing=F))
  
  geneThresh <- min(which(simGenes5[sortSim,p] > 
                            apply(simGenes5[sortSim,],1,min)))
  
  markerGenes <- names(sort(simGenes5[,p],
                            decreasing=F)[1:geneThresh])
  
  markerGenes <- markerGenes[apply(normLumiDat.Gene[markerGenes,],1,max) > 5]
  
  if (length(markerGenes)>1) {
    markerGenes <- markerGenes[simGenes5[markerGenes,p] <=
                                 apply(simGenes5[markerGenes,],1,min)]
  }
  
  if (!file.exists(sprintf('MarkerGenes%d.csv',p))) {
    
    markerGeneInfo <- select(org.Hs.eg.db,keys=markerGenes,
                             columns=c('SYMBOL','ENTREZID'),keytype='SYMBOL')
    markerGeneInfo$info <- ""
    i <- 0
    for (g in markerGeneInfo$ENTREZID) {
      i <- i+1
      if (is.na(g)) next
      markerGeneInfo[i,'info'] <- GetGeneInfo(g)$genesummary
    }
    
    write.table(markerGeneInfo,file=sprintf('MarkerGenes%d.csv',p),row.names=F,sep=",")
  }
  
  for (g in markerGenes) {
    errbar(x=as.numeric(substr(grep('S',
                                    colnames(CG$Pmean),value=T),1,2)),
           y=normLumiDat.Gene[g,grep('S',colnames(CG$Pmean))],
           yplus=normLumiDat.Gene[g,grep('S',colnames(CG$Pmean))] +
             normLumiDat.Gene.SD[g,grep('S',colnames(CG$Pmean))],
           yminus=normLumiDat.Gene[g,grep('S',colnames(CG$Pmean))] -
             normLumiDat.Gene.SD[g,grep('S',colnames(CG$Pmean))],
           col='black',lty=2,lwd=2,type='o',
           ylim=c(0,max(normLumiDat.Gene[g,]+normLumiDat.Gene.SD[g,])),
           xlab='passage', ylab=g)
    errbar(x=as.numeric(substr(grep('R',colnames(CG$Pmean),
                                    value=T),1,2)),
           y=normLumiDat.Gene[g,grep('R',colnames(CG$Pmean))],
           yplus=normLumiDat.Gene[g,grep('R',colnames(CG$Pmean))] +
             normLumiDat.Gene.SD[g,grep('R',colnames(CG$Pmean))],
           yminus=normLumiDat.Gene[g,grep('R',colnames(CG$Pmean))] -
             normLumiDat.Gene.SD[g,grep('R',colnames(CG$Pmean))],
           col='red',lty=1,lwd=2,type='o',add=T,errbar.col='red')
    legend('bottomleft',lty=c(2,1),pch=c(19,19),col=c('black','red'),
           legend=c('Sensitive','Resistant'))
    
    
  }
  
  
  
  dev.off()
}

## gene set analysis with pathways downloaded from MSigDB

library('GSA')
h <- GSA.read.gmt('h.all.v5.0.symbols.gmt')
names(h$genesets) <- h$geneset.names
h <- h$genesets
h <- sapply(h,intersect,row.names(normLumiDat.Gene))

c2 <- GSA.read.gmt('c2.cp.v5.0.symbols.gmt')
names(c2$genesets) <- c2$geneset.names
c2 <- c2$genesets
c2 <- sapply(c2,intersect,row.names(normLumiDat.Gene))


c2.PID <- c2[substr(names(c2),1,4)=="PID_"]
c2.PID <- c2.PID[sapply(c2.PID,length)>=5]

library('CoGAPS')
h.GSStat <- calcCoGAPSStat(Amean = nP5$Amean,
                           Asd = nP5$Asd,
                           GStoGenes = h)

c2.GSStat <- calcCoGAPSStat(Amean=nP5$Amean,
                            Asd=nP5$Asd,
                            GStoGenes=c2.PID)

# up in sensitive
names(which(h.GSStat$GSUpreg[2,] < 0.05))
#[1] "HALLMARK_HYPOXIA"                   
#[2] "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
#[3] "HALLMARK_APICAL_SURFACE"            
#[4] "HALLMARK_UV_RESPONSE_UP"            
#[5] "HALLMARK_PANCREAS_BETA_CELLS"    

# up in resistant
names(which(h.GSStat$GSUpreg[3,] < 0.05))
#[1] "HALLMARK_TNFA_SIGNALING_VIA_NFKB"          
#[2] "HALLMARK_IL6_JAK_STAT3_SIGNALING"          
#[3] "HALLMARK_ESTROGEN_RESPONSE_EARLY"          
#[4] "HALLMARK_ESTROGEN_RESPONSE_LATE"           
#[5] "HALLMARK_INTERFERON_ALPHA_RESPONSE"        
#[6] "HALLMARK_INTERFERON_GAMMA_RESPONSE"        
#[7] "HALLMARK_COMPLEMENT"                       
#[8] "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
#[9] "HALLMARK_INFLAMMATORY_RESPONSE"            
#[10] "HALLMARK_COAGULATION"                      
#[11] "HALLMARK_ALLOGRAFT_REJECTION"              
#[12] "HALLMARK_KRAS_SIGNALING_UP"                
#[13] "HALLMARK_KRAS_SIGNALING_DN"

h.GSStat.Vis <- t(apply(h.GSStat$GSUpreg,2,function(x){
  apply(sweep(nP5$Pmean,1,-log(x),FUN="*"),2,sum)
}))

h.GSStat.Vis.SD <- t(apply(h.GSStat$GSUpreg,2,function(x){
  apply(sweep(nP5$Psd,1,-log(x),FUN="*"),2,sum)
}))

h.GSStat.Vis.norm <- apply(h.GSStat.Vis,1,sum)
h.GSStat.Vis <- sweep(h.GSStat.Vis,1,h.GSStat.Vis.norm,FUN="/")
h.GSStat.Vis.SD <- sweep(h.GSStat.Vis.SD,1,h.GSStat.Vis.norm,
                         FUN="/")

pdf(sprintf('HallmarkGS.pdf',p))
for (p in names(sort(apply(h.GSStat$GSUpreg[3:5,],2,min)))) {
  
  if (min(h.GSStat$GSUpreg[,p])>0.1) next
  
  P.Exprs <- nP5$Amean[h[[p]],]%*%nP5$Pmean
  P.Exprs <- sweep(P.Exprs,1,apply(P.Exprs,1,sum),FUN="/")
  
  matplot(x=as.numeric(substr(grep('S',
                                   colnames(h.GSStat.Vis),
                                   value=T),1,2)),
          t(P.Exprs[,grep('S',colnames(P.Exprs))]),col='grey',
          ylim=c(0,pmax(max(h.GSStat.Vis[p,]+h.GSStat.Vis.SD[p,]),
                        max(P.Exprs ))),
          xlab='passage', ylab=p,type='l')
  
  matlines(x=as.numeric(substr(grep('R',
                                    colnames(h.GSStat.Vis),
                                    value=T),1,2)),
           t(P.Exprs[,grep('R',colnames(P.Exprs))]),col='#FF000055',
           ylim=c(0,max(h.GSStat.Vis[p,]+h.GSStat.Vis.SD[p,])),
           xlab='passage', ylab=p,type='l')
  
  errbar(x=as.numeric(substr(grep('S',
                                  colnames(h.GSStat.Vis),
                                  value=T),1,2)),
         y=h.GSStat.Vis[p,grep('S',colnames(h.GSStat.Vis))],
         yplus=h.GSStat.Vis[p,grep('S',colnames(h.GSStat.Vis))] +
           h.GSStat.Vis.SD[p,grep('S',colnames(h.GSStat.Vis))],
         yminus=h.GSStat.Vis[p,grep('S',colnames(h.GSStat.Vis))] -
           h.GSStat.Vis.SD[p,grep('S',colnames(h.GSStat.Vis))],
         col='black',lty=2,lwd=2,type='o',
         ylim=c(0,max(h.GSStat.Vis[p,]+h.GSStat.Vis.SD[p,])),
         add=T)
  
  
  errbar(x=as.numeric(substr(grep('R',
                                  colnames(h.GSStat.Vis),
                                  value=T),1,2)),
         y=h.GSStat.Vis[p,grep('R',colnames(h.GSStat.Vis))],
         yplus=h.GSStat.Vis[p,grep('R',colnames(h.GSStat.Vis))] +
           h.GSStat.Vis.SD[p,grep('R',colnames(h.GSStat.Vis))],
         yminus=h.GSStat.Vis[p,grep('R',colnames(h.GSStat.Vis))] -
           h.GSStat.Vis.SD[p,grep('R',colnames(h.GSStat.Vis))],
         col='red',lty=1,lwd=2,type='o',errbar.col='red', 
         add=T)
  
  legend('bottomleft',lty=c(2,1),pch=c(19,19),col=c('black','red'),
         legend=c('Sensitive','Resistant'))
  
  
}
dev.off()


h.GSSig <- names(which(apply(h.GSStat$GSUpreg,2,min) > 0.1))



## Gene set analysis of TRANSFAC targets
load('TRANSFAC_Genes_2014.Rda')
TF2Gene <- sapply(TF2Gene,intersect,row.names(normLumiDat.Gene))
TF2Gene <- TF2Gene[sapply(TF2Gene,length)>=5]
TF2Gene <- TF2Gene[sapply(TF2Gene,length)<200]
calcCoGAPSStat(Amean=nP5$Amean,
               Asd=nP5$Asd,
               GStoGenes=TF2Gene) -> TFStats

# upreg
names(which(TFStats$GSUpreg[2,]<0.05))

# downreg
names(which(TFStats$GSUpreg[3,]<0.05))

TF.GSStat.Vis <- t(apply(TFStats$GSUpreg,2,function(x){
  apply(sweep(nP5$Pmean,1,-log(x),FUN="*"),2,sum)
}))

TF.GSStat.Vis.SD <- t(apply(TFStats$GSUpreg,2,function(x){
  apply(sweep(nP5$Psd,1,-log(x),FUN="*"),2,sum)
}))

TF.GSStat.Vis.norm <- apply(TF.GSStat.Vis,1,sum)
TF.GSStat.Vis <- sweep(TF.GSStat.Vis,1,TF.GSStat.Vis.norm,FUN="/")
TF.GSStat.Vis.SD <- sweep(TF.GSStat.Vis.SD,1,TF.GSStat.Vis.norm,
                          FUN="/")

pdf(sprintf('TFGS.pdf',p))
for (p in names(sort(apply(TFStats$GSUpreg[3:5,],2,min)))) {
  
  if (min(TFStats$GSUpreg[,p])>0.1) next
  
  P.Exprs <- nP5$Amean[TF2Gene[[p]],]%*%nP5$Pmean
  P.Exprs <- sweep(P.Exprs,1,apply(P.Exprs,1,sum),FUN="/")
  
  matplot(x=as.numeric(substr(grep('S',
                                   colnames(TF.GSStat.Vis),
                                   value=T),1,2)),
          t(P.Exprs[,grep('S',colnames(P.Exprs))]),col='grey',
          ylim=c(0,pmax(max(TF.GSStat.Vis[p,]+TF.GSStat.Vis.SD[p,]),
                        max(P.Exprs ))),
          xlab='passage', ylab=p,type='l')
  
  matlines(x=as.numeric(substr(grep('R',
                                    colnames(TF.GSStat.Vis),
                                    value=T),1,2)),
           t(P.Exprs[,grep('R',colnames(P.Exprs))]),col='#FF000055',
           ylim=c(0,max(TF.GSStat.Vis[p,]+TF.GSStat.Vis.SD[p,])),
           xlab='passage', ylab=p,type='l')
  
  errbar(x=as.numeric(substr(grep('S',
                                  colnames(TF.GSStat.Vis),
                                  value=T),1,2)),
         y=TF.GSStat.Vis[p,grep('S',colnames(TF.GSStat.Vis))],
         yplus=TF.GSStat.Vis[p,grep('S',colnames(TF.GSStat.Vis))] +
           TF.GSStat.Vis.SD[p,grep('S',colnames(TF.GSStat.Vis))],
         yminus=TF.GSStat.Vis[p,grep('S',colnames(TF.GSStat.Vis))] -
           TF.GSStat.Vis.SD[p,grep('S',colnames(TF.GSStat.Vis))],
         col='black',lty=2,lwd=2,type='o',
         ylim=c(0,max(TF.GSStat.Vis[p,]+TF.GSStat.Vis.SD[p,])),
         add=T)
  
  
  errbar(x=as.numeric(substr(grep('R',
                                  colnames(TF.GSStat.Vis),
                                  value=T),1,2)),
         y=TF.GSStat.Vis[p,grep('R',colnames(TF.GSStat.Vis))],
         yplus=TF.GSStat.Vis[p,grep('R',colnames(TF.GSStat.Vis))] +
           TF.GSStat.Vis.SD[p,grep('R',colnames(TF.GSStat.Vis))],
         yminus=TF.GSStat.Vis[p,grep('R',colnames(TF.GSStat.Vis))] -
           TF.GSStat.Vis.SD[p,grep('R',colnames(TF.GSStat.Vis))],
         col='red',lty=1,lwd=2,type='o',errbar.col='red', 
         add=T)
  
  legend('bottomleft',lty=c(2,1),pch=c(19,19),col=c('black','red'),
         legend=c('Sensitive','Resistant'))
  
  
}
dev.off()

library('network')
markerGenes3Graph <- read.table('tabdelimited.EXmXd_Oi8Rwa.txt',
                                header=F,sep="\t",stringsAsFactors=F)

markerGene3Orig <- read.table('MarkerGenes3.csv',header=T,sep=",",
                              stringsAsFactors=F)

for (i in 1:nrow(markerGenes3Graph)) {
  if (markerGenes3Graph[i,1] %in% row.names(normLumiDat.Gene) & 
        markerGenes3Graph[i,2] %in% row.names(normLumiDat.Gene)) next
  
  if (!(markerGenes3Graph[i,1] %in% row.names(normLumiDat.Gene))) {
    newGene <- intersect(select(org.Hs.eg.db,
                                keys = markerGenes3Graph[i,1],
                                keytype='SYMBOL',
                                columns=c('ALIAS'))$ALIAS,
                         markerGene3Orig$SYMBOL)
    if (length(newGene)>1) {stop('more than one gene found')}
    if (length(newGene)==0) {
      warning('no alias found for ', markerGenes3Graph[i,1])
      newGene <- markerGenes3Graph[i,1]
    }
    message(markerGenes3Graph[i,1], " to ", newGene)
    markerGenes3Graph[i,1] <- newGene
  }
  
  if (!(markerGenes3Graph[i,2] %in% row.names(normLumiDat.Gene))) {
    
    newGene <- intersect(select(org.Hs.eg.db,
                                keys = markerGenes3Graph[i,2],
                                keytype='SYMBOL',
                                columns=c('ALIAS'))$ALIAS,
                         markerGene3Orig$SYMBOL)
    if (length(newGene)>1) {stop('more than one gene found')}
    if (length(newGene)==0) {
      warning('no alias found for ', markerGenes3Graph[i,2])
      newGene <- markerGenes3Graph[i,2]
    }
    message(markerGenes3Graph[i,2], " to ", newGene)
    markerGenes3Graph[i,2] <- newGene
    
  }
}

markerGenes3inNetwork <- intersect(unique(c(markerGenes3Graph$V1,
                                            markerGenes3Graph$V2)),
                                   markerGene3Orig$SYMBOL)

markerGenes3Network <- matrix(0,nrow=length(markerGenes3inNetwork),
                              ncol=length(markerGenes3inNetwork),
                              dimnames=list(markerGenes3inNetwork,
                                            markerGenes3inNetwork))

for (i in 1:nrow(markerGenes3Graph)) {
  
  if (!(markerGenes3Graph$V1[i] %in% markerGenes3inNetwork)) next
  if (!(markerGenes3Graph$V2[i] %in% markerGenes3inNetwork)) next
  
  
  markerGenes3Network[markerGenes3Graph$V1[i],
                      markerGenes3Graph$V2[i]] <- markerGenes3Graph$V15[i]
}
network(markerGenes3Network) -> markerGenes3Network

normLumiDat.Gene.marker3.Z <- sweep(normLumiDat.Gene[markerGenes3inNetwork,],
                                    1,apply(normLumiDat.Gene[markerGenes3inNetwork,],
                                            1, mean))

normLumiDat.Gene.marker3.Z <- normLumiDat.Gene.marker3.Z / 
  normLumiDat.Gene.SD[markerGenes3inNetwork,]

pdf('NetworkHeatmap.pdf',width=20,height=10)
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}



layout(matrix(c(1,2,3),nrow=1),widths=c(0.1,.45,.45))


for (tidx in unique(substr(colnames(normLumiDat.Gene),1,2))) {
  color.bar(greenred(15)[c(2:7,9:14)],min=-3,nticks=13)
  
  set.seed(1234)
  plot(markerGenes3Network,usearrows = F,edge.col='grey',
       vertex.col=as.character(cut(normLumiDat.Gene.marker3.Z[,paste(tidx,'S',sep='-')],breaks=c(-Inf,seq(from=-3,to=3,by=0.5),Inf),labels=greenred(15)[c(1:7,9:15)])),
       displaylabels=T,label.cex=0.5,
       vertex.border=0)
  title(sprintf('Sensitive (%s)',tidx))
  set.seed(1234)
  plot(markerGenes3Network,usearrows = F,edge.col='grey',
       vertex.col=as.character(cut(normLumiDat.Gene.marker3.Z[,paste(tidx,'R',sep='-')],breaks=c(-Inf,seq(from=-3,to=3,by=0.5),Inf),labels=greenred(15)[c(1:7,9:15)])),
       displaylabels=T,label.cex=0.5,vertex.border=0)
  title(sprintf('Resistant (%s)',tidx))
}
dev.off()




