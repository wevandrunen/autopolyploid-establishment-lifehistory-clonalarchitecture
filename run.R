############################################
#### Running the Simulation in Parallel ####
############################################

# install the following packages
install.packages(c("tidyverse", "foreach", "doParallel", "data.table"))


########## Running Full Parameter Combination Sets ##########
# Use this code to run the simulation for a full parameter combination set ("combos"), replicated nreps times
# See below to run the simulation for a single parameter set, replicated nreps times

rm(list=ls())

library(tidyverse)
library(foreach)
library(doParallel)
library(data.table)
options(scipen = 999)

ncores <- detectCores() - 4 # If running on personal computer, recommend using n-1 cores to avoid crashes
registerDoParallel(ncores)
source("main.R", local=TRUE) # main.R (and reproduction.R, recruitment.R) must be located in the same directory

## set overall parameters
nreps <- 20 # number of replicates (NOTE: to maximize runtime:output, nreps should be a multiple of no_cores)
maxgen <- 500 # maximum number of generations
D <- 30 # population lattice dimensions
N <- 900 # population size

## parameter values per simulation set
# Example runs the core dispersal-inviability parameter set for all equal survivorship probabilites between cytotypes, no clonal reproduction, and no unreduced gamete production
ug_values <- c(0) # unreduced gamete production 
disp_values <- c(0.5,2,5) # pollen (dp) and seed (ds) average dispersal distances
ks2 <- c(0.1,0.5,0.9) # diploid selfed-seed inviability 
ks4 <- ks2 # tetraploid selfed-seed inviability 
sv2 <- seq(0,0.9,by=0.1) # diploid per generation survival probability
sv4 <- sv2 # tetraploid per generation survival probability
c2 <- c(0) # diploid ramet production
c4 <- c2 # tetraploid ramet production
dc2 <- c(0) # diploid ramet dispersal distances
dc4 <- dc2 # tetraploid ramet dispersal distances

## Obtain all parameter combinations
# NOTE: ensure parameter terms in "combos" are correct 
# adjust expand.grid terms as needed (e.g., remove sv4 if sv2=sv4, and retain both terms if sv2!=sv4)
# For the given parameter values above (90 combinations):
combos <- as.matrix(expand.grid(ug_values,disp_values,ks2,sv2,c2,dc2))
# All possible terms:
# combos <- as.matrix(expand.grid(ug_values,disp_values,ks2,ks4,sv2,sv4,c2,c4,dc2,dc4)) 
ncomb <- length(combos[,1])

Sys.time()
## Run Simulation
paraVect<-rep(0,times=13) # initialize vector for holding parameter values
for(i in 1:ncomb){
  
  # NOTE: adjust indexing as needed to pull correct parameters out of combos matrix
  paraVect[1] <- 1 # starting number of polyploids
  paraVect[2] <- 5 # ovules produced per individual
  paraVect[3] <- combos[i,1] # ug
  paraVect[4] <- combos[i,2] # dp
  paraVect[5] <- combos[i,2] # ds
  paraVect[6] <- combos[i,3] # ks2
  paraVect[7] <- combos[i,3] # ks4
  paraVect[8] <- combos[i,4] # sv2
  paraVect[9] <- combos[i,4] # sv4
  paraVect[10] <- combos[i,5] # c2
  paraVect[11] <- combos[i,5] # c4
  paraVect[12] <- combos[i,6] # dc2
  paraVect[13] <- combos[i,6] # dc4
  
  # custom combine function for parallel loop output
  custcom <- function(List1,List2) {
    masters <- c(List1[[1]], List2[[1]]) #append master arrays to list
    geninfos <- rbind(List1[[2]],List2[[2]]) #rbind geninfo matrixes together
    return(list(master=masters,geninfo=geninfos))
  }
  
  # run main simulation function
  store<-foreach(nr = 1:nreps,.combine="custcom",.init=list(list(),matrix(NA,ncol=60,nrow=0)),.packages='tidyverse',.inorder=FALSE) %dopar% {
    modeloutput<-main(D,N,paraVect,maxgen)
    list(list(modeloutput[[1]]),modeloutput[[2]])
  }
  
  # write generational data to file
  # insert description file name
  fwrite(store[[2]],file=paste(paste(c(i,D,maxgen,nreps),collapse="_"),c("file_name_gendata"),c(".csv"),sep=""),na="NA")
  
}
Sys.time()
stopImplicitCluster()

## Explanation of output:
# The simulation returns a 3-dimensional array and a matrix for each replicate
# In the code above, only the matrices (the per generation stats) are being written to output, one file per parameter combination [i]

# Arrays:
# The arrays are stored in a list (store$master, store[[1]]), where each entry (store$master[[1]],store$master[[2]], etc) contains the population information from one simulation replicate
# The first dimension goes from 1:maxgen (total generations), the second from 1:N (population size)
# The third dimension stores 5 variables, the per individual information for each generation in the simulation
# [x location, y location, genetID of individual, ploidy of individual, within generation plantID of individual]
# These are used to produce the "snapshots" of 2D population structure below (for a single parameter set)

# Matrices:
# The matrices are appended together as each replicate finishes (into store$geninfo, store[[2]]), so that the final output is a large matrix with maxgen*nreps rows
# Each row contains the "generation stats" for a single generation, within a particular simulation replicate
# There are 60 columns: 
# [generation #, number of 4x, number of 2x surviving this generation, number of 4x surviving this generation, 
# number of 2x genets, number of 4x genets, avg 4x genet size (number ramets), size of initital 4x genet, 
# number of 22 sexual propagule types produced by last generation (see Table S1 for types), number of 2x ramets produced, number of 4x ramets produced,
# number of 22 sexual propagule types recruited from last generation (see Table S1 for types), number of 2x ramets recruited, number of 4x ramets recruited, 
# avg empty spaces surrounding each surviving 2x individual, avg empty spaces surrounding each surviving 4x individual]
# Because of the order of steps taken per generation, remember that these stats contain a mix of information about the current *and previous* generation
# Notably, the recorded seeds made and recruited are those produced by individuals in the previous generation 





########## Running Single Parameter Set ##########

rm(list=ls())

library(tidyverse)
library(foreach)
library(doParallel)
library(data.table)
options(scipen = 999)

no_cores <- detectCores() - 1 # If running on personal computer, recommend using n-1 cores to avoid crashes
registerDoParallel(no_cores)
source("main.R", local=TRUE) # main.R (and reproduction.R, recruitment.R) must be located in the same directory

## set overall parameters
nreps <- 20 # number of replicates (NOTE: to maximize runtime:output, nreps should be a multiple of no_cores)
maxgen <- 500 # maximum number of generations
D <- 30 # population lattice dimensions
N <- 900 # population size

## single parameter set 
Vect<-c(1,5,0) # pstart, ovnumb, ug
dpsVect<-c(0.5,0.5) # dp, ds; pollen (dp) and seed (ds) average dispersal distances
ksVect<-c(0.1,0.1) # ks2, ks4; diploid and tetraploid selfed-seed inviability
svVect<-c(0.5,0.9) # sv2, sv4; diploid and tetraploid per generation survival
cVect<-c(0,0) # c2, c4; diploid and tetraploid ramet production per generation
dcVect<-c(0,0) # dc2, dc4; diploid and tetraploid ramet dispersal distances
paraVect<-c(Vect,dpsVect,ksVect,svVect,cVect,dcVect) 

# custom combine function for parallel loop output
custcomb <- function(List1,List2) {
  masters <- c(List1[[1]], List2[[1]]) #append master arrays to list
  geninfos <- rbind(List1[[2]],List2[[2]]) #rbind geninfo matrixes together
  return(list(master=masters,geninfo=geninfos))
}

## Run Simulation
Sys.time()
store<-foreach(nr = 1:nreps,.combine="custcomb",.init=list(list(),matrix(NA,ncol=60,nrow=0)),.packages='tidyverse',.inorder=FALSE) %dopar% {
  modeloutput<-main(D,N,paraVect,maxgen)
  list(list(modeloutput[[1]]),modeloutput[[2]])
}
Sys.time()
stopImplicitCluster()

## write generational data to file
fwrite(store[[2]],file=paste(paste(c(D,N,paraVect),collapse="_"),c("test"),c(".csv"),sep=""),na="NA")



########## Basic Visualizations for Single Parameter Set ##########

## number of tetraploids per generation. Each replicate is a separate line.
store2<-as.data.frame(store[[2]])
store3<-bind_cols(store2,Gen=as.vector(rep(c(1:maxgen),times=nreps)),RunNo=as.vector(rep(c(1:nreps),each=maxgen)))
ggplot(store3,aes(x=Gen,y=V2,colour=RunNo,group=RunNo))+
  geom_line(size=1.5) + 
  scale_y_continuous(limits=c(0,N)) +
  scale_x_continuous(limits=c(0,maxgen))+
  labs(x="Generation", y="Number of Tetraploids")
  theme_bw()+theme(legend.position="none") 

## 2D expansion of tetraploids into diploid population in 10 generation snapshots
# for single replicate (here, replicate #4)
snapshots<-c(1,seq(10,50,by=10))
df<-data.frame()
for(i in snapshots){
  dfsnap<-data.frame(Snap=rep(i,times=N),X = store$master[[4]][i,,1],Y = store$master[[4]][i,,2],Ploidy = store$master[[4]][i,,4],GenetID = store$master[[4]][i,,3])
  df<-rbind(df,dfsnap)}
ggplot(df,aes(X, Y)) + 
  geom_tile(data=filter(df,Ploidy==4),aes(fill=as.factor(GenetID))) + 
  geom_tile(data=filter(df,Ploidy==2),fill="lightgrey")+
  coord_fixed() + facet_wrap(~Snap,nrow=round(length(snapshots)/5))+
  theme_void() + scale_y_continuous(breaks=(seq(0,D,by=5)+0.5))+ scale_x_continuous(breaks=(seq(0,D,by=5)+0.5))+
  theme(legend.position="none",panel.grid.major=element_line(linetype=2,colour="slategrey"),strip.text=element_text(size=14,face="bold"),panel.ontop=TRUE,axis.ticks=element_line(colour="white")) 

