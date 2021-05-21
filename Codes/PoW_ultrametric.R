###Libraries required for the analyses
library("ape")
library("stats")
library("nleqslv")
library("scales")
library("treeio")
library("ggtree")
require(tidytree)
##########################################

#include the path directory to your tree file here: path/to/directory/filename.tree
directory <- "path/to/directory/filename.tree"
d <- read.beast(directory)
t <- as_tibble(d)

###uncomment below to visualise the imported tree
#ggtree(d) + geom_text(aes(label=node), hjust=-.3) 

###uncomment below to set the PoW parameters for HCV
#muMax = 3.98*10**(-2)
#meanRate = 1.2*10**(-3)

###uncomment below to set the PoW parameters for sarbecovirus
#muMax = 3.98*10**(-2)
#meanRate = 5.5*10**(-4)

###uncomment below to set the PoW parameters for sarbecovirus
#muMax = 4.0*10**(-2)
#meanRate = 6.11*10**(-6)

delta = 0.2
steps = seq(0,log10(muMax)+9,delta)
steps2 = seq(1,length(steps),1)
lambda <- function(l){
  final <- sum(exp(l*(1/delta*steps+1))*10**(-9+steps))/sum(exp(l*steps2)) - meanRate
  return (final)
}

LAMBDA = nleqslv(0,lambda)$x
proportions = exp(LAMBDA*(1/delta*steps+1))/sum(exp(LAMBDA*steps2))
ratespread = 10**(-9+steps)
alpha = 3/4

func <- function(t){
  final <- sum(alpha*proportions*(1-exp(-ratespread*t/alpha))) - totalBranchLength
  return(final)
}

nTips=(length(t$node)+1)/2
tip_dates=c()
for(i in c(1:nTips)){
  temp=strsplit(t$label[[i]], "")
  count=length(temp[[1]])
  while((temp[[1]][[count]]!='/' && temp[[1]][[count]]!='|') && count > 1){
    count = count - 1
  }
  if(count>1){
    tip_dates <- c(tip_dates,as.double(paste(temp[[1]][c((count+1):(count+1+3))],collapse = '')))
  }
}

if(length(tip_dates)>0){
  for(i in c(1:nTips)){
    t$height[[i]] =tip_dates[[i]]
  }
}

ultimate_count=0
check_parent=c()
if(length(tip_dates)>0){
  while(ultimate_count!=(length(t$node)-length(tip_dates))){
    for(i in c(1:length(t$node))){
      node_1=parent(d,i)
      for(j in c(1:length(t$node))){
        node_2=parent(d,j)
        if(node_1==node_2 && i!=j && !(node_1 %in% check_parent) && t$height[[i]]>1 && t$height[[j]]>1){
          ultimate_count = ultimate_count + 1
          check_parent <- c(check_parent,node_1)
          t$height[[node_1]]=(t$height[[i]]+t$height[[j]])/2
        }
      }
    }
  }
}

if(length(tip_dates)>0){
  for(i in c(1:length(t$node))){
    t$height[[i]] = max(tip_dates) - t$height[[i]]
  }
}

node_track <- c()
height_track <- c()
height_track_low <- c()
height_track_high <- c()
internal_height <- c()
internal_node <- c()
internal_height_low <- c()
internal_height_high<- c()
for(i in c(1:nTips)){
  tracker=i
  while(tracker!=rootnode(d)){
    temp_1=tracker
    tracker = parent(d,tracker)
    temp_2=tracker
    if(temp_1 %in% c(1:nTips) && !(is.na(t$height_median[[temp_1]]))){
      height_1 <- t$height_median[[temp_1]]
      height_1_low <- t$height_0.95_HPD[[temp_1]][[1]]
      height_1_high <- t$height_0.95_HPD[[temp_1]][[2]]
      totalBranchLength = t$height_median[[temp_2]]
      height_2 <- c(nleqslv(0,func)$x)
      node_track <- c(node_track,temp_1)
      height_track <- c(height_track,height_1)
    } else if(temp_1 %in% c(1:nTips) && is.na(t$height_median[[temp_1]])){
      height_1 <- 0
      totalBranchLength = t$height_median[[temp_2]]
      height_2 <- c(nleqslv(0,func)$x)
      node_track <- c(node_track,temp_1)
      height_track <- c(height_track,height_1)
    } else if(!(temp_1 %in% c(1:nTips))){
      totalBranchLength = t$height_median[[temp_1]]
      height_1 <- c(nleqslv(0,func)$x)
      totalBranchLength = t$height_0.95_HPD[[temp_1]][[1]]
      height_1_low <- c(nleqslv(0,func)$x)
      totalBranchLength = t$height_0.95_HPD[[temp_1]][[2]]
      height_1_high <- c(nleqslv(0,func)$x)
      totalBranchLength = t$height_median[[temp_2]]
      height_2 <- c(nleqslv(0,func)$x)
      if(!(temp_1 %in% node_track)){
        internal_node <- c(internal_node,temp_1)
        internal_height <- c(internal_height,(t$height[[temp_1]] + height_1))
        internal_height_low <- c(internal_height_low,(t$height[[temp_1]] + height_1_low))
        internal_height_high <- c(internal_height_high,(t$height[[temp_1]] + height_1_high))
      }
      node_track <- c(node_track,temp_1)
      height_track <- c(height_track,height_1)
      height_track_low <- c(height_track_low,height_1_low)
      height_track_high <- c(height_track_high,height_1_high)
    }
    if((height_2 - height_1)<0){
      print('Warning: negative branch length detected here:')
      print(c(temp_1,temp_2,height_1,height_2))
    }
    t$branch.length[[temp_1]] <- (height_2 - height_1)
  }
}

for(i in c(1:length(internal_height))){
  t$height_median[[internal_node[[i]]]] <- as.double(internal_height[[i]])
  t$height_0.95_HPD[[internal_node[[i]]]][[1]] <- as.double(internal_height_low[[i]])
  t$height_0.95_HPD[[internal_node[[i]]]][[2]] <- as.double(internal_height_high[[i]])
}

totalBranchLength = t$height_median[[rootnode(d)]]
height_1 <- c(nleqslv(0,func)$x)
t$height_median[[rootnode(d)]] = as.double(height_1)

totalBranchLength = t$height_0.95_HPD[[rootnode(d)]][[1]]
height_1 <- c(nleqslv(0,func)$x)
t$height_0.95_HPD[[rootnode(d)]][[1]] = as.double(height_1)

totalBranchLength = t$height_0.95_HPD[[rootnode(d)]][[2]]
height_1 <- c(nleqslv(0,func)$x)
t$height_0.95_HPD[[rootnode(d)]][[2]] = as.double(height_1)

d=as.phylo(t)

#uncomment below to see the PoW-transformed tree
#ggtree(d) + geom_text(aes(label=label), hjust=-.3) + geom_text(aes(x=branch, label=scientific(branch.length)), vjust=-1.0, color='black')


### uncomment below to see the divergence time for the following pairs of sarbecovirus sequences:
### SARS-CoV-1 and the closest bat relative:
# t$label[[23]]
# t$height_median[[parent(d,23)]]
# t$height_0.95_HPD[[parent(d,23)]]
### SARS-CoV-2 and the RaTG13:
# t$label[[58]]
# t$height_median[[parent(d,58)]]
# t$height_0.95_HPD[[parent(d,58)]]
### SARS-CoV-2 and Pangoline:
# t$label[[1]]
# t$height_median[[parent(d,1)]]
# t$height_0.95_HPD[[parent(d,1)]]


### uncomment below to see the Root height of the tree:
# t$height_median[[rootnode(d)]]
# t$height_0.95_HPD[[rootnode(d)]]

###uncomment below to save the PoW transformed tree in your local directory
#write.beast(as.treedata(t), file = 'path/to/directory/transformed_tree.tree', translate = TRUE, tree.name = "tree")
