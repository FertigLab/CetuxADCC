load('NormGeneData.Rda')

library('CoGAPS')
normLumiDat.Gene.SD[is.na(normLumiDat.Gene.SD)] <- 
  0.1*normLumiDat.Gene[is.na(normLumiDat.Gene.SD)]

normLumiDat.Gene.SD <- pmax(normLumiDat.Gene.SD,
                            0.1*normLumiDat.Gene)

nP4 <- gapsRun(D=normLumiDat.Gene,S = normLumiDat.Gene.SD,
        nFactor="4", nEquil="50000", nSample="50000")

save(list=ls(), file='CoGAPS.nP4')
