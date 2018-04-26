#Hi Ben or whomever ends up working on this code! Assuming you already have the random forest model built,
#this code will take the new data, predict the subzones, use the edatopic overlap to determine likley site
#series, and produce summary tables. The first part of the code up to the function and loop section as line 175
#hasn't really been changed, but almost everything after that has been. For creating the summary statistics,
#there are some cutoffs/constants that might be changed or could even be adjustable on the web tool. I've 
#tried to comment these. Please don't hesitate to get in touch if anything doesn't makes sense. Thanks!

.libPaths("E:/R packages")
#.libPaths("F:/R/Packages")
require (RGtk2)
require(plyr)
require (rChoiceDialogs)
require (data.table)
require(doBy)
require (utils)
require(labdsv)
require(tools )
require(svDialogs)
require(tcltk)
require(randomForest)
require(foreach)
require(dplyr)
require(reshape2)
require(reshape)
library(doParallel)
require(data.table)
#===============================================================================
# Clear environment and add in previous rF model 
#===============================================================================

rm(list=ls())

###Function to find the proportion edatopic overlap###
edaOverlap <- function(dat1, dat2){
  return(length(dat1[dat1 %in% dat2])/length(dat2))
}

######Functions used in suitability rules####

bifurcTrend <- function(x){ ##percent improve, stable, and same - bifurcated model in summary data
  if(x[1] >= 25 & x[3] >= 25){
    return("TRUE")
  }else{
    return("FALSE")
  }
}

newSuit <- function(x){ ##Calculate new suitability accounting for current suitability
  suitLookup <- data.frame(SuitCurr = c(1,2,3,5), Value = c(-0.5, -0.3, 0.3, 0.6)) ##VALUES FOR CURRENT SUIT COULD BE ADJUSTED
  val <- suitLookup$Value[match(x[5], suitLookup$SuitCurr)]
  return(val+(1*x[1]+2*x[2]+3*x[3]+5*x[4]))
}

newSuitnoCurrent <- function(x){ ##New suitability without current suitability
  return(1*x[1]+2*x[2]+3*x[3]+5*x[4])
}

###stepSum calculates the difference between suitability in each time period####
stepSum <- function(x){##x is dataframe with colums Period, Spp, CurrentSuit, NewSuit
  if(x[1] == 2025){
    return(as.numeric(x[3])-as.numeric(x[4]))
  }else if(x[1] == 2055){
    y <- as.numeric(rawDat[rawDat$Spp == x[2] & rawDat$FuturePeriod == 2025, "NewSuit"])
    return(y - as.numeric(x[4]))
  }else{
    y <- as.numeric(rawDat[rawDat$Spp == x[2] & rawDat$FuturePeriod == 2055, "NewSuit"])
    return(y - as.numeric(x[4]))
  }
}

###Calculate percent of models improving, same, or declining####
modDir2 <- function(x, direction){#columns 1,2,3,X,Currentsuit
  curr <- x[5]
  if(curr == 5 | is.na(curr)){
    curr <- 4
  }
  x[5] <- 0
  x <- c(0,x)
  improve <- sum(x[1:curr])
  same <- x[curr + 1]
  decline <- sum(x[(curr+2):length(x)])
  if(direction == "Improve"){
    return(round(improve*100, digits = 0))
  }else if(direction == "Stable"){
    return(round(same*100, digits = 0))
  }else if(direction == "Decline"){
    return(round(decline*100, digits = 0))
  }else{
    return("error")
  }
}

###Flag model agreement in mid rotation period
modAgree <- function(x){
  m <- max(x)
  if(m >= 80){
    return("High")
  }else if(m >= 60){
    return("Moderate")
  }else{
    return("Low")
  }
}

combineList <- function(...) {##Combine multiple dataframe in foreach loop
  mapply(FUN = rbind, ..., SIMPLIFY=FALSE)
}

replaceZone <- function(x){##x is dataframe with MergedBGC, Table, SS_NoSpace
  pattern <- x[1]
  replacement <- x[2]
  SS_NoSpace <- x[3]
  if(!is.na(replacement)){ ##|replacement != "NS"
    return(sub(pattern,replacement,SS_NoSpace))
  }else{
    return(SS_NoSpace)
  }
}

wd=tk_choose.dir()
setwd(wd)

###enter model file name
fname="BGCv10_2000Pt_Rnd_Normal_1961_1990MSY_RFmodelKiriFinal.Rdata"
#fname = (file.choose())
load(fname)

#choose data file: ClimateBC output including all variables and all future periods
fplot=(file.choose())

Columns = c("GCM", "ID1", "ID2", "Latitude", "Longitude", "Elevation", "AHM", "bFFP",
            "CMD07","DD5_sp","EMT","Eref_sm","EXT","FFP","MCMT","MSP",
            "PPT07","PPT08", "PPT05","PPT06","PPT09", "SHM","TD","Tmax_sp","Tmin_at",
            "Tmin_sm","Tmin_wt", "PPT_at","PPT_wt", "PAS","eFFP",
            "Eref09","MAT","Tmin_sp","CMD")

Y1 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
fplot=basename(file_path_sans_ext(fplot))
Y1 <- Y1[!is.na(Y1$ID1),]

gc()


colnames (Y1) [1:3] = c("Model", "SiteNo", "BGC")

Y1$Model <- as.factor (Y1$Model)
Y1$SiteNo <- as.factor(Y1$SiteNo)
Y1$BGC <-gsub(" ", "", Y1$BGC, fixed = TRUE)
Y1$BGC <- as.factor(Y1$BGC)

#####generate some addtional variables
####
Y1$PPT_MJ <- Y1$PPT05 + Y1$PPT06 # MaY/June precip
Y1$PPT_JAS <- Y1$PPT07 + Y1$PPT08 + Y1$PPT09 # July/Aug/Sept precip
Y1$PPT.dormant <- Y1$PPT_at + Y1$PPT_wt # for calculating spring deficit
Y1$CMD.def <- 500 - (Y1$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
Y1$CMD.def [Y1$CMD.def < 0] <- 0 #negative values set to zero = no deficit
Y1$CMDMax <- Y1$CMD07
Y1$CMD.total <- Y1$CMD.def + Y1$CMD

##Specify variables in model
model = "5.1"
VarList = c("AHM", "bFFP","CMD.total","DD5_sp","EMT","Eref_sm","EXT","FFP","MCMT","MSP",
            "PPT_JAS","PPT_MJ","PPT06","SHM","TD","Tmax_sp","Tmin_at","Tmin_sm","Tmin_wt",
            "PAS","CMD.def","CMDMax","eFFP","Eref09","MAT","PPT07","Tmin_sp")
List = c("Model", "SiteNo", "BGC")
Y1save = Y1
Y1$BGC  <- as.factor(Y1$BGC)
Y1.sub=Y1[,names(Y1) %in% c(List,VarList)]

Y1.sub$BGC <- gsub(" ", "", Y1.sub$BGC, fixed = TRUE)

Y1.sub$BGC <- factor(Y1.sub$BGC, levels=unique(Y1.sub$BGC))

##Predict future subzones######
Y1.sub$BGC.pred <- predict(BGCmodel, Y1.sub[,-c(1:2)]) 
gc()

Y1.sub$PlotNo <- attr(Y1.sub, "row.names")
S1 <- subset(Y1, select=c("Model", "SiteNo", "Latitude", "Longitude", "Elevation"))
S2 <- subset(Y1.sub, select=c("BGC", "BGC.pred"))   

Y2.sub <- cbind.data.frame(S1, S2)

Y2.sub <- Y2.sub[,c("Model", "SiteNo", "BGC", "BGC.pred", "Latitude", "Longitude", "Elevation")]

Y2.sub$BGC.pred <- gsub(" ", "", Y2.sub$BGC.pred, fixed = TRUE)

Y2.sub$BGC.pred <- factor(Y2.sub$BGC.pred, levels=unique(Y2.sub$BGC.pred))

Y2.sub$BGC <- gsub(" ", "", Y2.sub$BGC, fixed = TRUE)

Y2x <- as.character(Y2.sub$Model)
Ystr <- strsplit(Y2x, "_")
Y4 <- matrix(unlist(Ystr), ncol=3, byrow=TRUE)
Y2.sub <- cbind(Y4, Y2.sub)
colnames(Y2.sub)[1:3]=c("GCM","Scenario", "FuturePeriod" )

Y3.sub = Y2.sub

Y3.sub$Num <- 1

Y3.sub$BGC.pred <- gsub(" ", "", Y3.sub$BGC.pred, fixed = TRUE)

#===============================================================================
#
#          Import suitability and edatopic tables
#           
#===============================================================================

edatopename="EdatopicWithSpecialCurrent"
edatopename2=paste(wd,"/",edatopename,".csv",sep="")
E1 <-read.csv(edatopename2,stringsAsFactors=FALSE,na.strings=".")

E1 <- E1[,-c(5)]
E1 <- unique(E1)

###create list of focal BGCs & edatopic space
e1 <- as.list(unique(Y3.sub$BGC), all.names=FALSE)
edatopic1 <- E1[E1$MergedBGC %in% e1,]
edatopic1$Codes[edatopic1$Codes == ""] <- NA

###create list of predicted BGCs and edatopic space
e2 <- as.list (unique(Y3.sub$BGC.pred), all.names=FALSE)
edatopic2 <- E1[E1$MergedBGC %in% e2,]
edatopic2$Codes[edatopic2$Codes == ""] <- NA

#################Import and Build Tree species suitability##############
treesuit="TreeSppSuit_v10.6"
treesuit2=paste(wd,"/",treesuit,".csv",sep="")
S1 <- read.csv(treesuit2,stringsAsFactors=F,na.strings=".")
S1 <- unique(S1)

#===============================================================================
# Builds list of all BGCs, Future BGCs, and Site Series
#===============================================================================


Y3.sub1 <-Y3.sub
Y3.sub1$BGC <- gsub("[[:space:]]","",Y3.sub1$BGC)
Y3.sub1$BGC.pred <- gsub("[[:space:]]","",Y3.sub1$BGC.pred)
Y3.sub1 <- Y3.sub1[,c("SiteNo","FuturePeriod","BGC","BGC.pred")]
Y3.sub1 <- Y3.sub1[order(Y3.sub1$SiteNo, Y3.sub1$FuturePeriod, Y3.sub1$BGC,Y3.sub1$BGC.pred),]

BGCminProp <- 0 ##to exclude BGCs with low prediction rates
average <- "No"##"Yes" ##

#####Average points (if specified) and remove BGCs with low predictions####

if(average == "Yes"){
  Y3.sub1$BGC.len <- ave(Y3.sub1$BGC, Y3.sub1$FuturePeriod, Y3.sub1$BGC, FUN = length)
  Y3.sub1$Pred.len <- ave(Y3.sub1$BGC, Y3.sub1$FuturePeriod, Y3.sub1$BGC, Y3.sub1$BGC.pred, FUN = length)
  Y3.sub1$BGC.len <- as.numeric(Y3.sub1$BGC.len); Y3.sub1$Pred.len <- as.numeric(Y3.sub1$Pred.len)
  Y3.sub1$BGC.prop <- Y3.sub1$Pred.len/Y3.sub1$BGC.len
  Y3.sub1 <- Y3.sub1[Y3.sub1$BGC.prop > BGCminProp,] ###Change this to a proportion to exclude BGCs with low prediction
  Y3.sub1 <- within(Y3.sub1, SiteNo <- match(BGC, unique(BGC)))
}else{
  Y3.sub1$BGC.len <- ave(Y3.sub1$BGC, Y3.sub1$SiteNo, Y3.sub1$FuturePeriod, Y3.sub1$BGC, FUN = length)
  Y3.sub1$Pred.len <- ave(Y3.sub1$BGC, Y3.sub1$SiteNo, Y3.sub1$FuturePeriod, Y3.sub1$BGC, Y3.sub1$BGC.pred, FUN = length)
  Y3.sub1$BGC.len <- as.numeric(Y3.sub1$BGC.len); Y3.sub1$Pred.len <- as.numeric(Y3.sub1$Pred.len)
  Y3.sub1$BGC.prop <- Y3.sub1$Pred.len/Y3.sub1$BGC.len
  Y3.sub1 <- Y3.sub1[Y3.sub1$BGC.prop > BGCminProp,] ###Change this to a proportion to exclude BGCs with low prediction
}

Y3.sub1 <- within(Y3.sub1,
                  BGC.len <- ave(BGC, SiteNo, FuturePeriod, BGC, FUN = length),
                  Pred.len <- ave(BGC, SiteNo, FuturePeriod, BGC, BGC.pred, FUN = length),
                  BGC.len <- as.numeric(BGC.len), Pred.len <- as.numeric(Pred.len),
                  BGC.prop <- Pred.len/BGC.len)

Y3.sub1 <- Y3.sub1[order(Y3.sub1$SiteNo,Y3.sub1$FuturePeriod,Y3.sub1$BGC,Y3.sub1$BGC.pred),]
BGC.pred.out <- unique(Y3.sub1[,c(1:4,7)])

BGClist = (unique(Y3.sub1$BGC))
FuturePeriod.list <- as.list(unique(Y3.sub1$FuturePeriod))
BGCfutures.list <- as.list(unique(Y3.sub1$BGC.pred)) ### to use later on to limit the site series
BGCfocalE <- edatopic1[edatopic1$MergedBGC %in% Y3.sub1$BGC  ,] ### extracts edatopic space for BGC focal of SiteNo
BGCfutureE <- edatopic2[edatopic2$MergedBGC %in% Y3.sub1$BGC.pred  ,] #extracts edatopic info only for predicted BGCs
Y3.sub1$SiteNo <- as.character(Y3.sub1$SiteNo)
SiteNo.list = as.list(unique(Y3.sub1$SiteNo))

gc()



###OPTIONAL: Set up to run loops in parallel###
require(doParallel)
set.seed(123321)
coreNo <- makeCluster(detectCores() - 1)
registerDoParallel(coreNo, cores = detectCores() - 1)
Cores <- as.numeric(detectCores()-1)
clusterEvalQ(coreNo, .libPaths("E:/R packages"))
#cl = makePSOCKcluster(6)
#registerDoParallel(cl)
#clusterEvalQ(cl, .libPaths("E:/R packages"))

#=======================================================================
#             LOAD FUNCTIONS USED IN LOOP AND FOR SUMMARY STATS-- NEW!!!
#=======================================================================


#=======================================================================================  
#####Nested foreach loops to caclulate site series - now matches special SS together#####
#    To run loops in parallel, change the %do% in the first loop to %dopar%
#    Would also need to load all packages used inside loop onto workers using .packages command
#========================================================================================
SiteNo.suit <-  foreach(SNL = SiteNo.list, .combine = rbind, .packages = c("doBy","foreach")) %do% {##  for each SiteNo in the data
  
  options(warn=0)
  
  cat("===========================================","\n")
  cat(SNL, "\n", "===========================================","\n")
  
  SiteFuture.suit <- foreach(i = FuturePeriod.list, .combine = rbind)  %do% {
    
    Y3.each <- Y3.sub1[Y3.sub1$SiteNo %in% SNL ,] ## extracts data for each site
    Y3.each <- Y3.each[Y3.each$FuturePeriod %in% i,] ##extracts data for each time period
    Y3.siteno <- as.list(unique (Y3.each$SiteNo))
    Y3.BGC <- as.list(unique(Y3.each$BGC))
    Y3.BGC.pred<- unique(Y3.each$BGC.pred)
    BGCfocalE <- edatopic1[edatopic1$MergedBGC %in% Y3.BGC  , ] ### extracts edatopic space for BGC focal of SiteNo
    BGCfutureE <- edatopic2[edatopic2$MergedBGC %in% Y3.BGC.pred  , ] #extracts edatopic info only for predicted BGCs
    
    Y3.SSlist = as.list(unique(BGCfocalE$SS_NoSpace))
    
    FTS2 <-  foreach(SS = Y3.SSlist, .combine =rbind, .packages = c("doBy","foreach")) %dopar% {      ##  for each site series for a SiteNo BGC
  
      options(warn=0)
      
      cat("===========================================","\n")
      
      cat(SS, "\n", "===========================================","\n")
      
      SSfocal <- BGCfocalE[BGCfocalE$SS_NoSpace %in% SS ,] ###find focal site series cells
      SSfocalE <- as.list(unique(SSfocal$Edatopic))
      
      ##select site series only with some edatopic overlap with SSfocal
      
      SSfutureE <- BGCfutureE[BGCfutureE$Edatopic %in% SSfocalE,]
      futureZones <- unique(SSfutureE$MergedBGC)
      futureSS.names <- unique(SSfutureE$SS_NoSpace)
      if(length(SSfutureE$Edatopic) > 0){
        
        SSfocal <- BGCfocalE[BGCfocalE$SS_NoSpace %in% SS,]
        SSfuture <- BGCfutureE[BGCfutureE$SS_NoSpace %in% futureSS.names,]
        
        ##match site series within each projected subzone
        futureSS <- foreach(futSS = futureZones, .combine = rbind) %do% {
          dat <- SSfuture[SSfuture$MergedBGC == futSS,]
          fut <- dat[dat$Edatopic %in% SSfocal$Edatopic,]
          if(any(fut$Codes[!is.na(fut$Codes)] %in% SSfocal$Codes[!is.na(SSfocal$Codes)]) & (length(fut$Codes[!is.na(fut$Codes)])/length(fut$Edatopic)) > 0.1){###Are there matchin special edatopic cells?
            oldSp <- unique(SSfocal$Codes[!is.na(SSfocal$Codes)])
            newSp <- unique(fut$SS_NoSpace[match(oldSp, fut$Codes)])
            dat <- unique(dat[dat$SS_NoSpace %in% newSp, -c(4:5)]) #which ones have the new special edatope
            dat$alloverlap <- 1/length(dat$SS_NoSpace)
            
          }else{
            if(all(is.na(SSfocal$Codes))){
              dat <- dat[is.na(dat$Codes),] ####remove special edatopes
            }
              dat <- foreach(x = unique(as.character(dat$SS_NoSpace)), .combine = rbind) %do% {
                dat1 <- dat[dat$SS_NoSpace == x,]
                Overlap <- edaOverlap(SSfocal$Edatopic, dat1$Edatopic)
                Revoverlap <- edaOverlap(dat1$Edatopic, SSfocal$Edatopic)
                dat1 <- unique(dat1[-c(4:5)])
                dat1$alloverlap <- Overlap*Revoverlap
                dat1 <- as.data.frame(dat1)
            }
            ##if not special, loop through each SS to calculate overlap
            
          }
          dat
        }
        
        if(length(futureSS) > 0){
        futureSS <- futureSS[futureSS$alloverlap > 0,] ###Adjust this to remove low overlap SS
        Y3.each <- unique(Y3.each)
        Y3.each <- Y3.each[Y3.each$BGC.pred %in% futureSS$MergedBGC,] ###Remove predictions with no overlap and recalculate
        Y3.each$BGC.prop <- Y3.each$BGC.prop/sum(Y3.each$BGC.prop)
        futureSS <- merge(futureSS, Y3.each[,c("BGC.pred","BGC.prop")], by.x = "MergedBGC", by.y = "BGC.pred", all.x = TRUE)
        
        
        ####Calculate the SS ratio
        SSoverlap <- summaryBy(alloverlap~MergedBGC, data=futureSS, id = 'SS_NoSpace', FUN=c(sum))
        futureSS$overlaptot<- SSoverlap$alloverlap.sum[match(futureSS$MergedBGC, SSoverlap$MergedBGC )]
        futureSS$SSratio <- futureSS$alloverlap/futureSS$overlaptot
        #summaryBy(SSratio ~ MergedBGC, data=futureSS, FUN=c(sum)) #### for checking that SSratio sums to 100%
        #summaryBy(SSratio ~ SS_NoSpace, data=futureSS, FUN=c(sum)) #Ummm?
        
        
        ####Calculated the overall site series probability
        futureSS$SSprob <- (futureSS$BGC.prop * futureSS$SSratio)
        sum(futureSS$SSprob)
        
        futureSS$SSCurrent <- rep(SS,length(futureSS$SSprob))
        futureSS <- futureSS[,-c(4:7)]
        
        futureSS$FuturePeriod <- as.character(i)  
        futureSS$SiteNo <- as.character(SNL)
        
        futureSS <- as.data.frame(futureSS)
        }else{
          warning("No matching edatopes - some predictions have been removed")
        }
      }
      
    } #For each Site
  } #For each year
  
} # for all


#################End of FOREACH LOOP###################
#######################################################
#####SiteNo.suit now contains projected BGC units as SS_NoSpace######
SiteNo.suit <- SiteNo.suit[!is.na(SiteNo.suit$SSprob),]

Crosswalk <- read.csv("Crosswalk.csv")
###Import Data for stocking Standards####
StandDat <- read.csv("StockStands_v11.csv")

StandInfo <- read.csv("StockingInfo.csv")
StandInfo <- unique(StandInfo)

StandHeight <- read.csv("StockingHeight.csv", stringsAsFactors = FALSE)
StandHeight$SiteSeries <- gsub("[[:space:]]","",StandHeight$SiteSeries)

StandDat$Footnotes <- apply(StandDat[,8:12], 1, function(x) paste(x[!is.na(x)], collapse = ","))
StandDat <- StandDat[,-c(8:12)]
StandDat$SppFN <- paste(StandDat$Species,"(",StandDat$Footnotes,")", " ", sep = "")
StandDat$Species <- as.character(StandDat$Species)
StandDat$SppFN <- ifelse(StandDat$Footnotes != "", paste(StandDat$Species,"(",StandDat$Footnotes,")", " ", sep = ""),
                         StandDat$Species)

FeasTrajLookup <- data.frame(SuitDiff = c(0,1,2,-1,-2,9,8,7), Flag = c("S","I1","I2","D1","D2","A1","A2","A3"))

####Merge with Spp suitability and create output tables#####
############################################################


StockingNames <- c("Unit","Standard","SSName","Primary","Preferred","Secondary","Acceptable","Tertiary",
"StockingTarget","StockingMINpa","StockingMINp","StockingDelay","AssessmentEarliest","AssessmentLatest","H1","H2","H3","H4","H5")

##Need to select region for creating reference guide####
region <- "Pr George"


#===================================================================================
#####Foreach loops to calculate summary statistics within each current BGC unit#####
######allOutput is a list with 3 dataframes#########################################
#===================================================================================
allOutput <- foreach(Site = unique(SiteNo.suit$SiteNo), .combine =  combineList, .multicombine = TRUE) %:%
  foreach(SS = unique(SiteNo.suit$SSCurrent[SiteNo.suit$SiteNo == Site]), .combine =  combineList, .multicombine = TRUE, .packages = "reshape") %dopar%{
    
    #==========================================================================================
    ###Prepare Data
    #==========================================================================================
    
    test <- SiteNo.suit[SiteNo.suit$SSCurrent == SS & SiteNo.suit$SiteNo == Site,-2]
    SuitSx <- S1[,c(2:4)]
    SuitSx <- SuitSx[order(SuitSx$Unit, SuitSx$Spp, -SuitSx$Suitability),] ##Use most conservative suitability if multiple exist
    SuitSx <- SuitSx[!duplicated(SuitSx[,1:2]),]
    
    ###Check if any subzones are in the crosswalk table####
    if(any(test$MergedBGC %in% Crosswalk$Modeled)){
      ####Merge with crosswalk table and update SS names####
      test <- merge(test, Crosswalk[!Crosswalk$Modeled %in% S1$BGC,], by.x = "MergedBGC", by.y = "Modeled", all.x = TRUE)
      test$NewSS <- apply(test[,c("MergedBGC","Tables","SS_NoSpace")],1,replaceZone)
      comb <- merge(test, SuitSx, by.x = "NewSS", by.y = "Unit", all.x = TRUE) ##Merge SS with suitability tabl
    }else{
      comb <- merge(test, SuitSx, by.x = "SS_NoSpace", by.y = "Unit", all.x = TRUE) ##Merge SS with suitability table
    }
    
    comb$Spp[is.na(comb$Spp)] <- "NS"
    comb <- comb[order(comb$FuturePeriod, comb$SS_NoSpace),]
    comb$Int <- interaction(comb$SS_NoSpace,comb$FuturePeriod)
    index <- which(comb$Int[-1] != comb$Int[-length(comb$Int)]) #index at each factor level
    index <- c(0,index, length(comb$SSprob))
    allSp <- unique(comb$Spp)
    allSp <- allSp[!is.na(allSp)]
    
    ####following loop creates rows for not suitable species in each projected Site Series###
    for(i in 1:(length(index) - 1)){
      Sp <- comb$Spp[(index[i]+1):index[i+1]]
      missing <- allSp[!allSp %in% Sp]
      if(length(missing) > 0){
        new <- comb[rep(index[i+1], length(missing)),]
        new$Spp <- missing
        new$Suitability <- NA
        comb <- rbind(comb, new)
      }
    }
    
    comb <- comb[order(comb$FuturePeriod, comb$SS_NoSpace, comb$Spp),]
    
    ####Create dummy matrix with all suitabilities incase some are missing###
    comb$Suitability[is.na(comb$Suitability)] <- "X"
    temp <- comb[c(1,1,1,1),]
    temp$Spp <- "T_Rex"
    temp$Suitability <- c("1","2","3","X")
    temp$SSprob <- 0
    
    comb <- rbind(comb,temp)
    
    ####new data frame with proportion votes for each suit class##
    numVotes <- cast(comb, Spp + FuturePeriod + SSCurrent ~ Suitability, value = "SSprob", fun.aggregate = sum) ###votes for each suitability
    ##numVotes$Sum <- rowSums(numVotes[,c(4:7)]) ##check votes sum to 1 (ignoring rounding errors)
    
    ####Remove all dinosaurs###
    numVotes <- numVotes[numVotes$Spp != "T_Rex",]
    
    ###merge in current suitability
    Zone <- sub("/.*","",SS)
    NewSS <- SS
    #####Assign NewSS the subzone from crosswalk table if predicted subzone not in suit table
    if(Zone %in% Crosswalk$Modeled & !Zone %in%  S1$BGC){
      NewSS <- sub(Zone, as.character(Crosswalk$Tables[Crosswalk$Modeled == Zone]),SS)
    }
    SuitSub <- SuitSx[SuitSx$Unit == NewSS,2:3]
    numVotes <- merge(numVotes, SuitSub, by = "Spp", all.x = TRUE)
    numVotes$Suitability[is.na(numVotes$Suitability)] <- 5
  
    #==========================================================================================
    ###Apply suitability rules to create summary output and raw data output
    #============================================================================================
    
    ###Create statistics for raw data####
    rawDat <- as.data.frame(numVotes)
    rawDat$ModAgree <- apply(rawDat[,4:7],1,FUN = max)
    rawDat$NewSuit <- apply(rawDat[,4:7],1,FUN = newSuitnoCurrent)
    rawDat <- rawDat[order(rawDat$Spp, rawDat$FuturePeriod),]
    rawDat$Suitability <- ifelse(rawDat$Suitability > 3.5,4,rawDat$Suitability) #
    rawDat$NewSuit <- ifelse(rawDat$NewSuit > 3.5,4,rawDat$NewSuit) #
    rawDat$SuitDiff <- apply(rawDat[,c("FuturePeriod","Spp","Suitability","NewSuit")],1,stepSum)
    
    ##THE BELOW CUTOFFS COULD BE ADJUSTED
    rawDat$PeriodTraj <- ifelse(rawDat$SuitDiff > 1.5, "Strongly Improving", 
                                ifelse(rawDat$SuitDiff > 0.2, "Improving",
                                       ifelse(rawDat$SuitDiff > -0.2, "No Change",
                                              ifelse(rawDat$SuitDiff <= -0.2, "Declining", "Strongly Declining"))))
    rawDat$Risk <- ifelse(rawDat$X >= 0.5, "Very High", ##CUTOFFS COULD BE ADJUSTED
                            ifelse(rawDat$X >= 0.35, "High",
                                   ifelse(rawDat$X >= 0.2, "Moderate","Low")))
    rawDat$modAgrClass <- ifelse(rawDat$ModAgree >= 0.8, "High",
                                 ifelse(rawDat$ModAgree >= 0.6, "Moderate","Low"))
    
    ##Feasibility of establishment using 2025 time period####
    Feas <- numVotes[numVotes$FuturePeriod == 2025,]
    Feas$FeasEstab <- Feas$`1`+Feas$`2`+(Feas$`3`*0.75) ###CONSTANT HERE COULD BE ADJUSTED
    Feas$NewSuit <- apply(Feas[,c(4:8)],1,newSuit)
    Feas$NewSuit <- round(Feas$NewSuit, digits = 0)
    Feas$NewSuit <- ifelse(Feas$NewSuit == 0, 1,
                           ifelse(Feas$NewSuit >= 4, 10, Feas$NewSuit))
    Feas$Estab.Risk <- ifelse(Feas$FeasEstab >= 0.8, "Low CC Risk", ##CUTOFFS COULD BE ADJUSTED
                             ifelse(Feas$FeasEstab >= 0.65, "Moderate CC Risk",
                                    ifelse(Feas$FeasEstab >= 0.5, "Considerable CC Risk","High Risk")))
    ###Feas$Estab.Risk <- ifelse(Feas$Estab.Risk == "High Risk" & Feas$NewSuit < 3.5, "Considerable CC Risk", Feas$Estab.Risk) ##Incase we get some weird output where it is suitable but says High Risk
    Feas$Suitability[Feas$Suitability == 5] <- 10
    Feas$SuitDiff <- Feas$Suitability - Feas$NewSuit
    Feas <- merge(Feas, FeasTrajLookup, by = "SuitDiff", all.x = TRUE)
    Feas$Flag <- as.character(Feas$Flag)
    Feas$Flag <- ifelse(Feas$Suitability == 10 & Feas$NewSuit == 10, "Not In",
                        ifelse(Feas$NewSuit == 10, "Removing", Feas$Flag))
    
    ####2055 data used for summary statistics####
    dat55 <- numVotes[numVotes$FuturePeriod == 2055, c(1,4:8)]
    dat55$Improve <- apply(dat55[,2:6],1,modDir2,direction = "Improve")
    dat55$Stable <- apply(dat55[,2:6],1,modDir2,direction = "Stable")
    dat55$Decline <- apply(dat55[,2:6],1,modDir2,direction = "Decline")
    dat55$Bifurc <- apply(dat55[,7:9],1,bifurcTrend)###test for bifurcation based on percent improve/decline
    dat55$NewSuit <- apply(dat55[,2:6],1,FUN = newSuitnoCurrent)
    dat55$Suit25 <- apply(Feas[,c(5:8)],1,newSuitnoCurrent) ###Get 2025 suitability
    dat55$NewSuit[dat55$NewSuit > 3.5] <- 4
    dat55$Suit25[dat55$Suit25 > 3.5] <- 4
    dat55$SuitDiff <- dat55$Suit25 - dat55$NewSuit
    dat55$OverallTraj <- ifelse(dat55$SuitDiff >= 2, "Strongly Improving",
                                ifelse(dat55$SuitDiff >= 0.5, "Improving",
                                       ifelse(dat55$SuitDiff >=-0.5, "No Change",
                                              ifelse(dat55$SuitDiff >= -2, "Declining", "Strongly Declining"))))
    dat55$OverallTraj <- ifelse(dat55$Bifurc, "Bifurcating", dat55$OverallTraj)
    dat55$ModAgree <- apply(dat55[,7:9],1,FUN = modAgree)
                                
    dat55 <- dat55[order(dat55$Spp),]
    Feas <- Feas[order(Feas$Spp),]
    
    ###create summary output
    outputSum <- data.frame(CurrentSS = Feas$SSCurrent, Species = Feas$Spp, CurrentSuit = Feas$Suitability, NewSuit = Feas$NewSuit, Estab.Trajectory = Feas$Flag, Estab.Risk = Feas$Estab.Risk,  
                            TrendMidRot = dat55$OverallTraj, ModelAgree = dat55$ModAgree, Improve = dat55$Improve, Stable = dat55$Stable, Decline = dat55$Decline)
    outputSum$SiteNo <- Site
    outputSum$CurrentSuit[outputSum$CurrentSuit == 10] <- "X"
    outputSum$NewSuit[outputSum$NewSuit == 10] <- "X"
    outputSum <- outputSum[,c(length(outputSum),1:length(outputSum)-1)]
    #############Create raw output##################
    outputRaw <- rawDat[,c(1:3,4:8,10,12:14)]
    outputRaw[,4:9] <- round(outputRaw[,4:9], digits = 2)
    outputRaw[,4:7] <- outputRaw[,4:7]*100
    outputRaw$SiteNo <- Site
    outputRaw <- outputRaw[,c(length(outputRaw),1:length(outputRaw)-1)]
    
    #==========================================================================================
    ###Create Stocking Standards Table
    #==========================================================================================
    dat <- outputSum
    dat <- dat[,2:5]
    NewSS <- SS
    Zone <- sub("/.*","",SS)
    if(!any(StandDat$SS_NoSpace %in% SS) & any(Crosswalk$Modeled %in% Zone)){ ##NewSS from crosswalk
      NewSS <- sub(Zone, as.character(Crosswalk$Tables[Crosswalk$Modeled == Zone]),SS)
    }else if(!any(StandDat$SS_NoSpace %in% SS) & !any(Crosswalk$Modeled %in% Zone)){###If not in either crosswalk or data, assign No SS
      NewSS <- "No SS"
    }
    
    ###use StockingStandards, Assesment, and height data to recreate reference guide with climate change rows
    if(NewSS == "No SS" | !any(StandDat$SS_NoSpace %in% NewSS)){ ##If No SS or no matching enteries, create NA entery
      StockStands <- as.data.frame(matrix(nrow = 1, ncol = 19))
      colnames(StockStands) <- StockingNames
      StockStands$Unit <- SS
      StockStands$Standard <- "No Data Available"
    }else{
      StockSub <- StandDat[StandDat$SS_NoSpace == NewSS & StandDat$Region == region,]
      Standard <- as.character(StockSub$Standard[1])
      if(length(StockSub$Region) == 0){
        region2 <- as.character(StandDat$Region[StandDat$SS_NoSpace == NewSS][1])
        StockSub <- StandDat[StandDat$SS_NoSpace == NewSS & StandDat$Region == region2,]
        Standard <- as.character(StockSub$Standard[1])
      }
      dat <- merge(dat,StockSub[,c(4,9)], by = "Species", all.x = TRUE)
      AssessSub <- StandInfo[StandInfo$Standard == Standard,][1,]
      if(length(AssessSub$Region) == 0){###If a subzone isn't found in specified region, choose from different region (theoretically shouldn't be an issue, but there's still some gaps in the data)
        region2 <- as.character(StandInfo$Region[StandInfo$SS_NoSpace == NewSS][1])
        AssessSub <- StandInfo[StandInfo$Region == region2 & StandInfo$SS_NoSpace == NewSS,]
        if(length(AssessSub$Region) == 0){
          temp <- matrix(nrow = 1, ncol = length(AssessSub))
          colnames(temp) <- colnames(AssessSub)
          AssessSub <- rbind(AssessSub, temp)
        }
      }
      HeightSub <- StandHeight[StandHeight$Standard == Standard,]
    
      if(length(HeightSub$Region) > 0){
        HeightSub$Comb <- paste(HeightSub$Species,": ", HeightSub$Height, sep = "") ##paste species to their heights
        HeightSub <- t(HeightSub$Comb)
        if(length(HeightSub) >= 5){
          HeightSub <- HeightSub[1:5]
          HeightSub <- as.data.frame(HeightSub)
          HeightSub <- t(HeightSub)
        }else{
          fill <- matrix(nrow = 1, ncol = (5 - length(HeightSub))) ##to account for different numbers of species
          HeightSub <- cbind(HeightSub,fill)
        }
      }else{
        HeightSub <- matrix(nrow = 1, ncol = 5)
      }
      
      colnames(HeightSub) <- c("H1","H2","H3","H4","H5")
      ###current reference guide
      StockStands <- data.frame(Unit = SS, 
                                Standard = StockSub$Standard[1],
                                SSName = AssessSub$SiteSeriesName,
                                Primary = paste(dat$SppFN[dat$CurrentSuit == 1 & dat$Species != "NS" & !is.na(dat$SppFN)], collapse = ","),
                                Preferred = paste(StockSub$SppFN[StockSub$PreferredAcceptable == "P" & StockSub$Species != "NS" & !is.na(StockSub$SppFN)], collapse = ","),
                                Secondary = paste(dat$SppFN[dat$CurrentSuit == 2 & dat$Species != "NS" & !is.na(dat$SppFN)], collapse = ","),
                                Acceptable = paste(StockSub$SppFN[StockSub$PreferredAcceptable == "A" & StockSub$Species != "NS" & !is.na(StockSub$SppFN)], collapse = ","),
                                Tertiary = paste(dat$SppFN[dat$CurrentSuit == 3 & dat$Species != "NS" & !is.na(dat$SppFN)], collapse = ","))
      StockStands <- cbind(StockStands, AssessSub[,7:12], HeightSub)
      
    }
    
    ####climate change informed
    StockStands2 <- data.frame(Unit = unique(outputSum$CurrentSS), 
                               Standard = "ClimateChange",
                               SSName = NA,
                               Primary = paste(dat$Species[dat$NewSuit == 1], collapse = ","),
                               Preferred = NA,
                               Secondary = paste(dat$Species[dat$NewSuit == 2], collapse = ","),
                               Acceptable = NA,
                               Tertiary = paste(dat$Species[dat$NewSuit == 3], collapse = ","))
    fill <- matrix(nrow = 1, ncol = (length(StockStands) - length(StockStands2)))
    StockStands2 <- cbind(StockStands2, fill)
    colnames(StockStands2) <- colnames(StockStands)
    
    #combine
    StockStands <- rbind(StockStands, StockStands2)
    StockStands <- cbind(c(Site, Site), StockStands)
    colnames(StockStands)[1] <- "SiteNo"
    
    list(outputRaw,outputSum, StockStands)
  }###end of loop


write.csv(allOutput[[2]], "SummaryExample.csv")
write.csv(allOutput[[1]], "RawDataExample.csv")
write.csv(allOutput[[3]], "ReferenceGuide.csv")
