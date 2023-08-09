###Loading libraries required for the PoW model
library("stringr")
library("ape")
library("nleqslv")
##########################################

#######
#Insert path to log file for the inferred short-term substitution rate (*.log output file from BEAST v1.10.4)
filepath_rates <- "/Path/to/directory/*.log"

#Insert path to the ultrametric distance trees (*.trees output file from BEAST v1.10.4) 
#and log file for the standard HKY substitution model used to construct the ultrametric distance trees  
filepath_distances <- "/Path/to/directory/*.trees"
filepath_HKYsubstitutionModel <- "/Path/to/directory/*.log"

#######

#######
#Import the log file for the inferred short-term rate
con = readLines(paste(filepath_rates, sep = ""))
#Import the log file for the standard HKY substitution model used to construct the ultrametric distance trees
con3 = readLines(paste(filepath_HKYsubstitutionModel, sep = ""))
#find the column number for meanRate or manually enter the column number for the posterior rate distribution
meanrate_column <- which(strsplit(con[3],"\t")[[1]]=="meanRate")
#find the position of the kappa and ACTG base frequencies or manually enter their median values below
kappa_column <- which(strsplit(con3[3],"\t")[[1]]=="kappa")
frequencies1_column <- which(strsplit(con3[3],"\t")[[1]]=="frequencies1")
#remove the first 10% of log file as burn-in
con <- con[c(5:length(con))]
con <- con[c((round(length(con)/10)+1):length(con))]
con3 <- con3[c(5:length(con3))]
con3 <- con3[c((round(length(con3)/10)+1):length(con3))]
#######

#sample_size sets the number of samples taken at random from the post-burn-in posterior distributions
sample_size=100

#######
#selecting states from the log file at random
sampled_states <- sample(con,sample_size)
sampled_states2 <- sample(con3,sample_size)
#######

#######
#read meanRate and other parameters of the substitution model

rates <- c()
freq1 <- c()
freq2 <- c()
freq3 <- c()
freq4 <- c()
K <- c()
for(i in 1:length(sampled_states)){
  rates <- c(rates,as.numeric(strsplit(sampled_states[[i]],"\t")[[1]][[meanrate_column]]))
  freq1 <- c(freq1,as.numeric(strsplit(sampled_states2[[i]],"\t")[[1]][[frequencies1_column]]))
  freq2 <- c(freq2,as.numeric(strsplit(sampled_states2[[i]],"\t")[[1]][[frequencies1_column+1]]))
  freq3 <- c(freq3,as.numeric(strsplit(sampled_states2[[i]],"\t")[[1]][[frequencies1_column+2]]))
  freq4 <- c(freq4,as.numeric(strsplit(sampled_states2[[i]],"\t")[[1]][[frequencies1_column+3]]))
  K <- c(K,as.numeric(strsplit(sampled_states2[[i]],"\t")[[1]][[kappa_column]]))
}

#######

#use the best fit value of muMax for the PoW transformation of every sampled distance tree
#for RNA viruses, muMax = 3.65*10**(-2)
#for single stranded DNA viruses, muMax = 2.0*10**(-2)
#for double stranded DNA viruses, muMax = 3.0*10**(-3)
muMax = 3.65*10**(-2)
muMaxs=rep(muMax,sample_size)

#uncomment below if you want to add variation around the value of muMax
#muMax = 3.65*10**(-2)
#muMax_stdev=10.0*10**-3
#muMaxs=rnorm(sample_size,muMax,muMax_stdev)

#con2 reads the distance tree file
con2 = readLines(paste(filepath_distances, sep = ""))

for(i in 1:length(con2)){
  if(paste(strsplit(con2[[i]],'')[[1]][c(1:4)],collapse = '')=='tree'){
    break
  }
}

#tree_list includes a list of constructed distance trees from BEAST, removing the first 10% as burn-in
tree_list <- con2[c(i:(length(con2)-1))] 
tree_list <- tree_list[c((round(length(tree_list)/10)+1):length(tree_list))]
sampled_trees <- sample(tree_list,sample_size)

#create a new tree file based on the subsampling in the previous step
#Insert the directory you would like <sampled_trees.trees> to be exported
recreated_treefile <- c(con2[c(1:(i-1))],sampled_trees,'End;')
write.table(recreated_treefile, "/Path/to/directory/sampled_trees.trees", sep = "\n", quote = FALSE, row.names = FALSE)

#import <sampled_trees.trees> created in the previous step from the directory
filename <- "/Path/to/directory/sampled_trees.trees"
importedtree <- read.nexus(filename, tree.names = NULL, force.multi = FALSE)

#transform each of the sampled distance trees using the PoW method 
for(el in 1:sample_size){
  print(el)
  muMax = muMaxs[[el]]
  meanRate = rates[[el]]
  pA = freq1[[el]]
  pC = freq2[[el]]
  pG = freq3[[el]]
  pT = freq4[[el]]
  kap = K[[el]]
  
  pY = pT + pC
  pR = pA + pG
  delta = 0.1
  betaMax = muMax/(2*((pT*pC + pA*pG)*kap + pY*pR))
  steps = seq(0,log10(betaMax)+9,delta)
  steps2 = seq(1,length(steps),1)
  lambda <- function(l){
    final <- sum(exp(l*(1/delta*steps+1))*(2*((pT*pC + pA*pG)*kap + pY*pR))*10**(-9+steps))/sum(exp(l*steps2)) - meanRate
    return (final)
  }
  
  LAMBDA = nleqslv(0,lambda)$x
  proportions = exp(LAMBDA*(1/delta*steps+1))/sum(exp(LAMBDA*steps2))
  ratespread = 10**(-9+steps)
  
  e2 <- function(t) exp(-ratespread*t)
  e3 <- function(t) exp(-(pR*kap+pY)*ratespread*t)
  e4 <- function(t) exp(-(pY*kap+pR)*ratespread*t)
  pTC <- function(t) pC+(pC*pR)/pY*e2(t)-pC/pY*e4(t)
  pAG <- function(t) pG+(pG*pY)/pR*e2(t)-pG/pR*e3(t)
  pTA <- function(t) pA*(1-e2(t))
  pTG <- function(t) pG*(1-e2(t))
  pCA <- function(t) pA*(1-e2(t))
  pCG <- function(t) pG*(1-e2(t))
  s1 <- function(t) sum(proportions*2*pT*pTC(t))
  s2 <- function(t) sum(proportions*2*pA*pAG(t))
  v <- function(t) sum(proportions*(2*pT*pTA(t)+2*pT*pTG(t)+2*pC*pCA(t)+2*pC*pCG(t)))
  a1 <- function(t) -log(1-(pY*s1(t))/(2*pT*pC)-v(t)/(2*pY))
  a2 <- function(t) -log(1-(pR*s2(t))/(2*pA*pG)-v(t)/(2*pR))
  b <- function(t) -log(1-v(t)/(2*pY*pR))
  k1 <- function(t) (a1(t)-pR*b(t))/(pY*b(t))
  k2 <- function(t) (a2(t)-pY*b(t))/(pR*b(t))
  func <- function(t){
    final <- 2*pT*pC/pY*(a1(t)-pR*b(t))+(2*pA*pG/pR)*(a2(t)-pY*b(t))+2*pY*pR*b(t) - totalBranchLength
    return(final)
  }
  
  
  allBranches=mapply(c,importedtree[[el]]$edge[,1],importedtree[[el]]$edge[,2],(importedtree[[el]]$edge.length),SIMPLIFY = FALSE)
  
  Allpaths <- list(list())
  for(j in 0:importedtree[[el]]$Nnode+1){
    branches <- list()
    for(i in 1:(length(nodepath(importedtree[[el]],j,importedtree[[el]]$Nnode+2))-1)){
      branches[[i]] <- c(nodepath(importedtree[[el]],j,importedtree[[el]]$Nnode+2)[i+1],nodepath(importedtree[[el]],j,importedtree[[el]]$Nnode+2)[i])
    }
    Allpaths[[j]] <- c(branches)
  }
  
  trackBranches <- list()
  globcount=0
  for(k in 1:length(Allpaths)){
    totalBranchLength=0
    count=0
    convertedHeight<-list()
    for(i in 1:length(Allpaths[[k]])){
      for(j in 1:length(allBranches)){
        if(sum(allBranches[[j]] %in% Allpaths[[k]][[i]],na.rm = TRUE)==2){
          matches=0
          count = count + 1
          globcount = globcount + 1
          trackBranches[[globcount]] <- c(allBranches[[j]][c(1:2)])
          for(m in 1:length(trackBranches)){
            if(sum(trackBranches[[m]] %in% Allpaths[[k]][[i]],na.rm = TRUE)==2){
              matches = matches + 1
            }
          }
          if(count==1 && matches==1){
            totalBranchLength = totalBranchLength + allBranches[[j]][[3]]
            convertedHeight[[count]] <- c(nleqslv(0,func)$x)
            importedtree[[el]]$edge.length[[j]]=convertedHeight[[count]]
          }
          if(count>1 && matches==1){
            totalBranchLength = totalBranchLength + allBranches[[j]][[3]]
            convertedHeight[[count]] <- c(nleqslv(0,func)$x)
            importedtree[[el]]$edge.length[[j]]=convertedHeight[[count]]-convertedHeight[[count-1]]
          }
          if(count>1 && matches>1){
            TRUE
          }
        }
      }
    }
  }
}
#save the PoW-transformed trees in a directory
#this file can then be directly imported to TreeAnnotator to construct a PoW-transformed maximum clade credibility tree
ape::write.nexus(importedtree, file="/Path/to/directory/PoWtransformed.trees")


