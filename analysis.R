library(Matrix)
##load("conserveAAhas.Rdata")
#load("check_var_out2.Rdata")
setwd("/home/eliot/Dropbox/Projects/Conserved_snps")
# outData = get(load("conservedcommon_marmoset_9483.Rdata"))
# randMat = get(load("permutationscommon_marmoset_9483.Rdata"))

outDataFs = list.files(path=".", pattern="conserved")
randMatFs = list.files(path=".", pattern="permutations")

for (fileInd in 1:1){##length(outDataFs)){
  
  load(outDataFs[fileInd])
  load(randMatFs[fileInd])

  outData = outDataMat
  rm(outDataMat)
  
  myColNames =
    c("rsID", "chr", "pos", "qryID", "mmAAsSTR", "mmAlsSTR", "mmFramesSTR", "mmFxnSTR", "hsRSID",
      "hsAAsSTR", "hsAlsSTR", "hsFramesSTR", "hsFxnSTR", "constatus", "alignIdentity")
  
  ##load("mm_hs_dbsnp.Rdata")
  ##outData = get(load("check_var_out_one_multi.Rdata"))
  colnames(outData) = myColNames
  
  randConMat = randConMat[!is.na(outData[,"alignIdentity"]),]
  outData = outData[ !is.na(outData[,"alignIdentity"]),  ]
  randMat = randConMat
  rm(randConMat)
  
  # goodInds = grep("NP_", outData[,"qryID"])
  # outDataBk = outData
  # outData = outData[goodInds,]
  # 
  # randMatBk = randMat
  # randMat = randMat[goodInds,]
  
  ## load("/Users/ezh/projects/conserved_SNPs/sequences_data/global_aln_pers2.Rdata")
  
  ## ## mmDBSNPalleles =
  ## ##     paste(mmDBSNP[,"allele_ref"], mmDBSNP[,"allele"], sep=";")
  ## ## hsDBSNPalleles =
  ## ##     paste(hsDBSNP[,"allele_ref"], hsDBSNP[,"allele"], sep=";")
  
  ## mmAAs = strsplit(outData[,"mmAAsSTR"], ";")
  ## mmRSids = sapply(strsplit(outData[,1], ";"),'[[',1)
  
  ## mmID = paste(mmRSids, outData[,"mmProtRef"],sep="|")
  ## matchInd = match(mmID, outPer[,1])
  ## removeInds = which(is.na(matchInd))
  
  ## matchIndFinal = matchInd[-removeInds]
  ## ##outDataBK = outData
  ## outData = outData[-removeInds,]
      
      ## missenseInds = grep("missense",outData[,"mmFxnSTR"]) 
  ## outData = outData[m    issenseInds, ]
  
  ## make histogram 
  mStart = 0; mEnd = 100
  int = floor((mEnd - mStart)/7)
  intervals = seq(mStart, mEnd, int)
  
  consCount = matrix(NA, nrow=4, ncol=length(intervals)-1)
  alignIdent = as.numeric(outData[,"alignIdentity"])
  
  nShuf = 5000
  
  randPercentiles = c()
  for ( ind in 1:(length(intervals)-1) ){
    
    lBound = intervals[ind]
    uBound = intervals[ind+1]
    currRandMat = randMat[which((alignIdent>lBound & alignIdent<=uBound)), ] 
    consStats = outData[ which((alignIdent>lBound & alignIdent<=uBound)), "constatus"]
    consCount[1,ind] = length(which(consStats==0))
    consCount[2,ind] = length(which(consStats==1))
    consCount[3,ind] = length(which(consStats==2))
    consCount[4,ind] = length(which(consStats==3))
    ##print(length(which(currRandMat==0))/nrow(currRandMat))
    qq = sum(apply(currRandMat, 1,function(x) all(x==0)))/nrow(currRandMat)
    print(qq)
    
    currRandVec = rep(0,nShuf)
    
    for (i in 1:nShuf){
      cc1 = length(which(currRandMat[,i]==0))
      cc2 = length(which(currRandMat[,i]==1))
      cc3 = length(which(currRandMat[,i]==2))
      cc4 = length(which(currRandMat[,i]==3))
      currRandVec[i] = (cc2+cc3++cc4)/(cc1+cc2+cc3+cc4)
    }
    ##currRandVec = c(currRandVec,  (sum(consCount[2:3,ind]))/sum(consCount[,ind]) )
    
    randPercentiles = rbind(randPercentiles, quantile(currRandVec, c(.025, .975)) )
  }
  
  x = c()
  for (ind in 1:(length(intervals)-1) ){
    xpoint = (intervals[ind] + intervals[ind+1])/2
    x = c(x, xpoint)
  }
  
  totCounts = apply(consCount, 2, sum)
  ## generate percentages
  pers = apply(consCount[2:4, ], 2, sum)/apply(consCount[c(1,3,4),], 2, sum)

}
##pers = apply(consCount[2:3,],2 ,sum )/totCounts*100

##library(ggplot2)

setwd("/Users/eliot/Desktop/memes/r_analysis")
##load("alleleDist.Rdata")

outPers = cbind(id=floor(x),pers=pers)
outPers = as.data.frame(outPers)

#persEx = c(28.42012, 25.81620, 29.19390, 30.38013, 33.35870, 34.15842)
#pers = c(24.03860 24.76292 27.54479 29.79003 32.63492 34.29092)
#totCounts = c(6839  9385 19314 15955 17659 24715)

Counts = totCounts
pp = ggplot(outPers, aes(x=id, y=pers)) + labs(x="Peptide Identity (%)", y="Conserved SNPs (%)")
pp = pp + geom_point(aes(size=Counts)) + xlim(30, 100) + ylim(15, 40)
##pp = pp + geom_point(aes(x=x, y=persEx, size=totCountsEx)) 

##png(file="mygraphic.png")
png(file="mygraphic.png",width=2000,height=1626,res=300)
pp
dev.off()
##pp + theme(legend.text = element_text(as.character(counts*1000)))

mmcondist = c(1544,15396,2811,2106,6251,3113)
mmintdist = c(3538,14557,3648,3563,14466,2741)

Alleles = c("A/T", "A/G", "A/C", "T/G", "T/C", "G/C")

names(mmcondist) = Alleles
names(mmintdist) = Alleles

##newall = c("A/C","A/G","A/T","G/C","T/C","T/G")

##newmmcondist = mmcondist[newall]
##newmmintdist = mmintdist[newall]

mmdata = 
  data.frame(al=Alleles, mc=mmConDist/sum(mmConDist), mce=mmConDistEs/sum(mmConDistEs),
             mnc=mmNonConDist/sum(mmNonConDist), mnce=mmNonConDistEs/sum(mmNonConDistEs), 
             mi=mmIntDist/sum(mmIntDist), ms=mmSynDist/sum(mmSynDist))
hsdata = 
  data.frame(al=Alleles, hc=hsConDist/sum(hsConDist), 
             hnc=hsNonConDist/sum(hsNonConDist), hi=hsIntDist/sum(hsIntDist),
             hs=hsSynDist/sum(hsSynDist))

mmdata = mmdata[order(mmdata$mc),]
mmdata$al=factor(mmdata$al,levels=mmdata$al)
colnames(mmdata) = c("Alleles", "Conserved", "Conserved Ess", 
                     "Non-conserved", "Non-conserved Ess", "Intronic", "Synonymous")

hsdata = hsdata[order(hsdata$hc),]
hsdata$al=factor(hsdata$al,levels=hsdata$al)
colnames(hsdata) = c("Alleles", "Conserved", "Non-conserved", "Intronic", "Synonymous")


library(reshape2)

mmdata = melt(mmdata, id.var=c("Alleles"))
colnames(mmdata)[2:3] = c("Type", "Frequency")

pp = ggplot(mmdata, aes(x=Alleles, y=Frequency, fill=Type) )
pp = pp + geom_bar( stat="identity", position="dodge") + 
  scale_fill_manual(values=c("#b772e2", "#9a1ae3", "#5b94e9", "#0948e8", "#35dcd2", "#7bc945")) +
  labs(y="Frequency (%)")

png(file="all_freq.png",width=2000,height=1219,res=300)
pp
dev.off()

### hs allele freq plot
hsdata = melt(hsdata, id.var=c("Alleles"))
colnames(hsdata)[2:3] = c("Type", "Frequency")

pp = ggplot(hsdata, aes(x=Alleles, y=Frequency, fill=Type) )
pp = pp + geom_bar( stat="identity", position="dodge") + 
  scale_fill_manual(values=c("#b772e2", "#9a1ae3", "#5b94e9", "#7bc945")) +
  labs(y="Frequency (%)")

png(file="all_freq_hs.png",width=2000,height=1219,res=300)
pp
dev.off()


###### sift plots
load("siftdata.Rdata")



persNOMINCON = apply(consCount[2:4, ], 2, sum)/totCounts
persExact = consCount[2, ]/totCounts
persStrong = consCount[3, ]/totCounts
persWeak = consCount[4, ]/totCounts 

persExToRest = consCount[2, ]/apply(consCount[c(3,4),], 2, sum) 

############ compare mutation distributions ##################
##############################################################

#alleleList = head(outData[,"hsAlsSTR"],n=50)

changeClasses = combn(c("A","T","G","C"),2)

getMutDist = function(alleleElement){
  if(length(alleleElement)>2){
    return(rep(0,6))
  }
  apply(
    changeClasses,2, function(x){ length(setdiff(x,alleleElement))==0 }
  )
}

applyMutDist = function(alleleList, mysep){
  ## need to change punctuation
  alleleList = strsplit(alleleList, mysep, fixed=T)
  counts = t(sapply(alleleList, getMutDist))
  ##counts = counts[!is.na(counts)]
  return(apply(counts,2,sum) )
  
}

mmConAlleles = outData[outData[,"constatus"]!="0","mmAlsSTR"]
mmNonConAlleles = outData[outData[,"constatus"]=="0","mmAlsSTR"]

hsConAlleles = outData[outData[,"constatus"]!="0","hsAlsSTR"]
hsNonConAlleles = outData[outData[,"constatus"]=="0","hsAlsSTR"]

#mmDBSNPdist = applyMutDist(mmDBSNPalleles)
#hsDBSNPdist = applyMutDist(hsDBSNPalleles)

mmConDist = applyMutDist(mmConAlleles, ";")
mmNonConDist = applyMutDist(mmNonConAlleles, ";")

hsConDist = applyMutDist(hsConAlleles, ";")
hsNonConDist = applyMutDist(hsNonConAlleles, ";")

unconDelID = unconVarsID[unconVarsID %in% rsIDs] 
conDelID = conVarsID[conVarsID %in% rsIDs]

hsDELCONall = outData[match(conDelID, outData[,"hsRSID"]), "hsAlsSTR"]
hsDELCONall = hsDELCONall[!is.na(hsDELCONall)]

hsDELunCONall = outData[match(unconDelID, outData[,"hsRSID"]), "hsAlsSTR"]
hsDELunCONall = hsDELunCONall[!is.na(hsDELunCONall)] 

hsDELconDist = applyMutDist(hsDELCONall, ";")
hsDELunconDist = applyMutDist(hsDELunCONall, ";")

mmIntFs =
  list.files(path="/Users/ezh/projects/conserved_SNPs/parseDBSNP_data/mouse/intron/29July16",
             pattern="*.Rdata", full.names=T)
hsIntFs =
  list.files(path="/Users/ezh/projects/conserved_SNPs/parseDBSNP_data/human/intron",
             pattern="*.Rdata", full.names=T)

mmIntFs = mmIntFs[1:2]
hsIntFs = hsIntFs[1:2]

getAlleles = function(fileList){
  outAll = c()
  for (file in fileList){
    currData = get(load(file))
    outAll = rbind(outAll, currData)
  }
  outAll = as.matrix(outAll)
  return(outAll)
}

mmIntAll = getAlleles(mmIntFs)
hsIntAll = getAlleles(hsIntFs)

mmSampInd = sample(1:length(mmIntAll), 50000)
hsSampInd = sample(1:length(hsIntAll), 50000)    

mmIntAllsub = mmIntAll[mmSampInd] 
hsIntAllsub = hsIntAll[hsSampInd]

mmIntDist = applyMutDist(mmIntAllsub)
hsIntDist = applyMutDist(hsIntAllsub)

## start here
## dist for all the mm dbSNP and hs dbSNP SNPs load("dbsnpALLDIST.Rdata")
##load("all_alleleDist.Rdata")

###############
mmConMat = rbind(mmConDist, mmIntDist)
mmNonConMat = rbind(mmNonConDist, mmIntDist)

hsConMat = rbind(hsConDist, hsIntDist)
hsNonConMat = rbind(hsNonConDist, hsIntDist)
###############

## chi-square test
getCramer = function(chisquareMat){
  xVal = chisq.test(chisquareMat)[["statistic"]]
  cramerV = sqrt( xVal/sum(chisquareMat) )
  return(cramerV)
}

mmConCram = getCramer(mmConMat)
mmNonCram = getCramer(mmNonConMat)

hsConCram = getCramer(hsConMat)
hsNonCram = getCramer(hsNonConMat)

mmCram = getCramer(rbind(mmConDist, mmNonConDist)) 
hsCram = getCramer(rbind(hsConDist, hsNonConDist))

mmtohsCram = getCramer(rbind(mmConDist, hsConDist))

## qq1 = getCramer(rbind(mmConDist, mmDBSNPdist))
## qq2 = getCramer(rbind(hsConDist, hsDBSNPdist))

delconMat = rbind(hsDELconDist, hsDELunconDist)

chisq.test(delconMat)
## chisq.test(mmNonConMat)

## chisq.test(hsConMat)
## chisq.test(hsNonConMat) 

## chisq.test(rbind(mmConDist, mmNonConDist))
## chisq.test(rbind(hsConDist, hsNonConDist))

############### counts of transversions and inversion statusus ####
###################################################################

## conMat = outData[outData[,"constatus"]!=0,]
## nonConMat = outData[outData[,"constatus"]==0,]

## conMatIn = conMat[,c("mmAlsSTR", "mmFramesSTR", "mmFxnSTR",
##                      "hsAlsSTR", "hsFramesSTR", "hsFxnSTR")]
## nonConMatIn = nonConMat[,c("mmAlsSTR", "mmFramesSTR", "mmFxnSTR",
"hsAlsSTR", "hsFramesSTR", "hsFxnSTR")]

########## percentage of conserved in deleterious SNPs
## find the percentage in all orthologous human variants

library(VariantAnnotation)

hsFlatFs =
  list.files(path="/Users/ezh/projects/conserved_SNPs/parseDBSNP_data/human/31July16", pattern="*.Rdata",
             full.names=T)

#hsFlatFs = hsFlatFs[1]

hsFlat = c()
for (hsFlatF in hsFlatFs){
  currFlat = get(load(hsFlatF))
  currFlat = currFlat[grep("missense", currFlat[,"fxn"]), ]
  ##NMid = sub("\\..*","", currFlat[,"prot_acc_ref"])
  NMid = sub("\\..*","", currFlat[,"mrna_acc_ref"])
  NPid = currFlat[,"prot_acc_ref"]
  currFlat = cbind(currFlat[,"rsID"], NMid, NPid)
  hsFlat = rbind(hsFlat, currFlat)
}

mmFlatFs =
  list.files(path="/Users/ezh/projects/conserved_SNPs/parseDBSNP_data/mouse/31July16", patter="*.Rdata",
             full.names=T)

mmFlat = c()
for (mmFlatF in mmFlatFs){
  currFlat = get(load(mmFlatF))
  currFlat = currFlat[grep("missense", currFlat[,"fxn"]), ]
  NMid = sub("\\..*","", currFlat[,"mrna_acc_ref"])
  currFlat = cbind(currFlat[,"rsID"], NMid)
  mmFlat = rbind(mmFlat, currFlat)
}


hsDBSNP = hsFlat
mmDBSNP = mmFlat

## load homgene human genes
load("/Users/ezh/reference/homgenes.Rdata")
outData = get(load("check_var_con_multi_all.Rdata"))
colnames(outData) = myColNames 
outData = outData[grep("missense", outData[,"mmFxnSTR"]),]

unconVars = outData[outData[,"constatus"]=="0",]
conVars = outData[outData[,"constatus"]!="0",]

horthVars = hsFlat[hsFlat[,"NMid"]%in%humanGenes, ]
morthVars = mmFlat[mmFlat[,"NMid"]%in%mouseGenes, ]

horthVars = unique(horthVars)
morthVars = unique(morthVars)

##### get out NPIDs
library(Biostrings)
protSeqHSF = "/Users/ezh/reference/refseq/protein/human/human_protein_merged.faa"
protSeqshs = readAAStringSet(protSeqHSF)
names(protSeqshs) = sapply(strsplit(names(protSeqshs),"|", fixed=T),'[[',4)

exactConVars = conVars[conVars[,"constatus"]!="0",]
mhsRSIDS = unlist(strsplit(exactConVars[,"hsRSID"], ";"))
mmatchInd = match(mhsRSIDS, horthVars[,1])
outNPIDS = horthVars[mmatchInd,3]
outNPIDS = outNPIDS[-which(is.na(outNPIDS))]

npid = as.data.frame(table(outNPIDS))
mind = match(npid$outNPIDS, names(protSeqshs))

removeInds = which(is.na(mind))

if (length(removeInds)>0){
  npid = npid[-removeInds,]
  mind = mind[-removeInds]
}

mlengths = width(protSeqshs[mind])
npid = cbind(npid,norm=npid$Freq/mlengths)

outNPIDs = npid[order(npid$norm, decreasing=T),]

## take top 2000
outNPIDs = as.character(outNPIDs[1:2000,1])

write(outNPIDs, file="outNPIDs_conserved2.txt")
#####

## find deleterious SNPS
hsVCF = readVcf("/Users/ezh/reference/dbSNP/human/VCF/sift_annotate/orth_filtered.recode_SIFTpredictions.vcf", "hs")

hsVCFInfo = unlist(info(hsVCF)[,48])
hsVCFInfoFilt = hsVCFInfo[-grep("WARNING", hsVCFInfo)]
##hsVCFsplit = strsplit(hsVCFInfo,'|')

hsDel = hsVCFInfoFilt[grep("DELETERIOUS", hsVCFInfoFilt)]
hsTol = hsVCFInfoFilt[-grep("DELETERIOUS", hsVCFInfoFilt)] 

hsDel = cbind( sapply(strsplit(hsDel,'|',fixed=T), '[[', 12), "DEL" )
hsTol = cbind( sapply(strsplit(hsTol,'|',fixed=T), '[[', 12), "TOL" )

## siftDB =
##     read.csv("/Users/ezh/reference/sift_ensemble/result.txt", sep='\t', header=F)

##siftDB = siftDB[,c("V1","V6")]
##siftDB = as.matrix(siftDB)

siftDBhs = rbind(hsDel,hsTol)
siftDBhs = unique(siftDBhs)

mmVCF = readVcf("/Users/ezh/reference/dbSNP/mouse/vcf/sift_annotate/orth_filtered.recode_SIFTpredictions.vcf", "hs")

mmVCFInfo = unlist(info(mmVCF)[,9])
mmVCFInfoFilt = mmVCFInfo[-grep("WARNING", mmVCFInfo)]
##hsVCFsplit = strsplit(hsVCFInfo,'|')

mmDel = mmVCFInfoFilt[grep("DELETERIOUS", mmVCFInfoFilt)]
mmTol = mmVCFInfoFilt[-grep("DELETERIOUS", mmVCFInfoFilt)]

mmDel = cbind( sapply(strsplit(mmDel,'|',fixed=T), '[[', 12), "DEL" )
mmTol = cbind( sapply(strsplit(mmTol,'|',fixed=T), '[[', 12), "TOL" )

siftDBmm = rbind(mmDel, mmTol)
siftDBmm = unique(siftDBmm)
##tolStat = siftDB[match(orthVars[,1], siftDB[,1]), 2]

getPercentage = function(orthVars, conVarsID, unconVarsID, status, siftDB){
  
  ## orthVars = horthVars
  ## siftDB = siftDBhs
  ## status="TOL"
  ## conVarsID = hsconVarsID; unconVarsID = hsunconVarsID
  
  tolStat = siftDB[match(orthVars[,1], siftDB[,1]), 2]    
  rsIDs = unique(orthVars[which(tolStat==status),1])
  
  rsIDsCon = unlist(strsplit(conVarsID, ";"))
  rsIDsUnCon = unlist(strsplit(unconVarsID, ";")) 
  
  ## sum(conVarsID%in%rsIDs) / length(conVarsID)
  ## sum(unconVarsID%in%rsIDs) / length(unconVarsID)
  
  ## sum(conVarsID%in%rsIDs) / (sum(conVarsID%in%rsIDs) + sum(unconVarsID%in%rsIDs))
  ## sum(unconVarsID%in%rsIDs) / (sum(conVarsID%in%rsIDs) + sum(unconVarsID%in%rsIDs))
  
  ## perCon = sum(conVarsID%in%rsIDs)/length(conVarsID)*100
  ## perUncon = sum(unconVarsID%in%rsIDs)/length(unconVarsID)*100
  
  ## perCon = sum(rsIDs%in%conVarsID)/length(rsIDs)*100
  ## perUncon = sum(rsIDs%in%unconVarsID)/length(rsIDs)*100
  perCon = sum(rsIDs%in%conVarsID)
  perUncon = sum(rsIDs%in%unconVarsID)
  
  ##print(head(tolStat))
  ##print(c(perCon, perUncon))
  
  return(c(perCon, perUncon))
}

hsconVarsID = unique(unlist(strsplit(conVars[,"hsRSID"],";")))
hsunconVarsID = unlist(strsplit(unconVars[,"hsRSID"],";"))
hsunconVarsID = unique(hsunconVarsID[hsunconVarsID!="NA"])

mmconVarsID = unique(unlist(strsplit(conVars[,"rsID"],";")))
mmunconVarsID = unlist(strsplit(unconVars[,"rsID"],";"))
mmunconVarsID = unique(mmunconVarsID[mmunconVarsID!="NA"])

hsconVarsID = hsconVarsID[hsconVarsID%in%siftDBhs[,1]]
hsunconVarsID = hsunconVarsID[hsunconVarsID%in%siftDBhs[,1]]

mmconVarsID = mmconVarsID[mmconVarsID%in%siftDBmm[,1]]
mmunconVarsID = mmunconVarsID[mmunconVarsID%in%siftDBmm[,1]]

## huam to mouse per of tol conserved
pers = getPercentage(horthVars, hsconVarsID, hsunconVarsID,"TOL",siftDBhs)
## human to mouse per of untol conserved
pers2 = getPercentage(horthVars, hsconVarsID, hsunconVarsID,"DEL",siftDBhs) 

## mouse to hum
pers3 = getPercentage(morthVars, mmconVarsID, mmunconVarsID,"TOL",siftDBmm) 
pers4 = getPercentage(morthVars, mmconVarsID, mmunconVarsID,"DEL",siftDBmm) 

#################################### check MAF's in alfred DB ############

alfredF = "/Users/ezh/reference/Alfred/alfred_short2.txt"
##alfredDB = read.csv(alfredF, sep="\t", header=F)

hsrsidCon = outData[outData[,14]!="0", 9]
hsrsidCon = hsrsidCon[!is.na(hsrsidCon)]
hsrsidCon = unlist(strsplit(hsrsidCon,";"))

hsrsidUnCon = outData[outData[,14]=="0", 9]
hsrsidUnCon = hsrsidUnCon[which(hsrsidUnCon!="NA")]
hsrsidUnCon = unlist(strsplit(hsrsidUnCon, ";"))

## filter for missense variants


## getMAF = function(crsid){

##     cMat = alfDB[alfDB[,1]==crsid,2:3]
##     ## As = median(cMat[which(cMat[,1]=="A"),2])
##     ## Ts = cMat[which(cMat[,1]=="T"),2]
##     ## Gs = cMat[which(cMat[,1]=="G"),2]
##     ## Cs = cMat[which(cMat[,1]=="C"),2]

## }

conMatchInd = match(hsrsidCon, alfDB[,1])



## humRSIDS =
##     read.csv("/Users/ezh/projects/conserved_SNPs/hsRSids.txt")

## create histogram of conserved SNPs
####################################################
### compare alignment tools

genOne = get(load("check_var_out2.Rdata"))
genTwo = get(load("check_var_out2.Rdata"))
genThree = get(load("check_var_con_multi_all.Rdata"))


colnames(genThree)[1:14] = colnames(genOne)

g1Con = genOne[,14]
g2Con = genTwo[,14]
g3Con = genThree[,14]

##matchInds = match()

getQuery = function(mat){
  return(paste(mat[,1], mat[,4], sep='|')) 
}

g1q = getQuery(genOne)
g2q = getQuery(genTwo)
g3q = getQuery(genThree)

getMatchedMat = function(q1, q2, q3, mat1, mat2, mat3){
  
  ## q1 = g1q
  ## q2 = g2q
  ## q3 = g3q
  ## mat1 = genOne
  ## mat2 = genTwo
  ## mat3 = genThree
  
  
  matchInds = match(q1, q2)
  removeInds = which(is.na(matchInds))
  
  if (length(removeInds)!=0){
    matchInds = matchInds[-removeInds]
    mat1 = mat1[-removeInds,]
    q1 = q1[-removeInds]
    
  }
  
  mat2 = mat2[matchInds,]
  
  matchInds = match(q1, q3)
  removeInds = which(is.na(matchInds))
  matchInds = matchInds[-removeInds]
  
  mat1 = mat1[-removeInds,]
  mat2 = mat2[-removeInds,]
  mat3 = mat3[matchInds, ]
  
  return( list(mat1, mat2, mat3) )
}

out = getMatchedMat(g1q, g2q, g3q, genOne, genTwo, genThree)

genOne2 = out[[1]]
genTwo2 = out[[2]]
genThree2 = out[[3]]

genOneConID = which(genOne2[,14] =="1" & genThree2[,14]=="0" )

ind = 128
qq = rbind(genOne2[ind,], genThree2[ind,-15])