##################################
#### Main Simulation Function ####
##################################

# Description:
# This function performs 4 primary actions: 
# 1) initialize the simulations 
# 2) oversee the reproduction/recruitment steps per generation 
# 3) store per generation information 
# 4) output the final data when simulations are complete


main<- function(D,N,paraVect,maxgen) 
{
  ########## Parameters ##########
  
  # D - individuals are situated on a D x D grid, where each grid unit can contain one individual
  # N - population size, constant
  G <- D*D # grid size (total # sites in population)
  
  pstart <- paraVect[1] # starting number of polyploids
  ovnumb <- paraVect[2] # ovules produced per individual 
  ug <- paraVect[3] # rate of unreduced gamete production 
  dp <- paraVect[4] # average pollen dispersal distance, in grid units
  ds <- paraVect[5] # average seed dispersal distance, in grid units
  ks2 <- paraVect[6] # probabilty diploid selfed seed is inviable 
  ks4 <- paraVect[7]# probabilty tetraploid selfed seed is inviable 
  sv2 <- paraVect[8] # diploid probability of survival to the next generation
  sv4 <- paraVect[9] # tetraploid probability of survival to the next generation
  c2 <- paraVect[10] # diploid ramet production 
  c4 <- paraVect[11] # tetraploid ramet production
  dc2 <- paraVect[12] # diploid average ramet dispersal distance, in grid units
  dc4 <- paraVect[13] # tetraploid average ramet dispersal distance, in grid units
  # maxgen - maximum number of generations  
  
  
  ########## Main Storage Objects ##########
  
  ## array of locations and traits of each individual for all generations
  master<-array(NA,dim=c(maxgen,N,5)) # gen, N rows, columns = [x,y,genetID,ploidy,plantID], 
  gendata<-matrix(NA,nrow=maxgen,ncol=60) # 60 columns = [gen,# polyploids, #2xsurvive, #4xsurvive, #2xgenets, #4xgenets, (neighbour info), 24 propagule types produced last gen, 24 propagule types recruited from last gen], free 2x spots, free 4x spots
  
  ## per generation population information
  initialinfo <- matrix(NA,nrow=N,ncol=5) # x, y, genetID, ploidy, plantID
  previnfo <- matrix(NA,nrow=N,ncol=5) # x, y, genetID, ploidy, plantID
  
  
  
  ########## Initial Generation ##########
  
  ## randomly assign initial individuals to cells
  initialsites<-arrayInd(sample(c(1:G),size=N,replace=FALSE),.dim=c(D,D))
  initialinfo[,c(1:2)]<-initialsites
  # randomized global genet IDs
  initialinfo[,3]<-sample(c(1:N),size=N,replace=FALSE) 
  
  ## randomly assign pstart individual(s) to be tetraploid
  ploidyvect<-rep(2,times=N)
  if(pstart == 0){
    ploidyvect<-ploidyvect
  } else {
    ploidyvect[sample(1:N,size=pstart,replace=FALSE)]<-4
  }
  initialinfo[,4]<-ploidyvect 
  pzero<-initialinfo[initialinfo[,4]==4,3] # record genet ID of original polyploid
  
  ## assign within-generation plant IDs
  initialinfo[,5]<-c(1:N)
  
  ## copy initial population to master array
  master[1,,] <- initialinfo
  gendata[1,] <- c(1,pstart,0,0,(N-pstart),pstart,rep(0,times=4),rep(0,times=48),0,0) # columns = generation, starting polyploids, 2xsurvival, 4xsurvival,#2xgenets, #4xgenets, (neighbourinfo), 24 propagule types produced, 24 propagule types germinated, free 2x spots, free 4x spots
  
  
  
  ########## Run Simulation ##########
  
  for (w in 2:maxgen)
  {
    previnfo <- master[(w-1),,] # previous generation population info [x,y,genetID,ploidy,plantID]
    
    ## Reproduction step
    source("reproduction.R", local=TRUE) # reproduction.R must be in the same directory
    reproductionlist <- reproduction(D,N,paraVect,previnfo) 
    # Reproduction output: 
    prevgenprops <- reproductionlist[[1]] # number of each propagule type made (sexual and clonal)
    seedpool <- reproductionlist[[2]] # pool of viable seeds 
    rametpool <- reproductionlist[[3]] # pool of ramets
    neighbourtotal4x <- reproductionlist[[4]] # average number of neighbours per tetraploid
    neighboursame4x <- reproductionlist[[5]] # average number of neighbours per tetraploid of the same cytotype
    neighbourtotal2x <- reproductionlist[[6]] # average number of neighbours per diploid
    neighboursame2x <- reproductionlist[[7]] # average number of neighbours per diploid of the same cytotype
    
    ## Survival & Recruitment steps
    source("recruitment.R",local=TRUE) # recruitment.R must be in the same directory
    recruitmentlist <- recruitment(D,N,paraVect,previnfo,seedpool,rametpool,pzero)
    # Survival & Recruitment output:
    prevgenseedgerm <- recruitmentlist[[1]] # number of each propagule type recruited (sexual and clonal)
    n2xsurvivors<- recruitmentlist[[3]] # number of surviving diploids from previous generation
    n4xsurvivors<- recruitmentlist[[4]] # number of surviving tetraploids from previous generation
    n2xgenet <- recruitmentlist[[5]] # number of diploid genets in population
    n4xgenet <- recruitmentlist[[6]] # number of tetraploid genets in population
    largest2xgenet <- recruitmentlist[[7]] # size of largest diploid genet
    largest4xgenet <- recruitmentlist[[8]] # size of largest tetraploid genet
    average4xgenet <- recruitmentlist[[9]] # average number of ramets per tetraploid genet
    sizepzero <- recruitmentlist[[10]] # current size of initial tetraploid genet
    average2xfreehood <- recruitmentlist[[11]] # average number of local vacant sites per diploid
    average4xfreehood <- recruitmentlist[[12]] # average number of local vacant sites per tetraploid
    
    npolyploids <- length(recruitmentlist[[2]][which(recruitmentlist[[2]][,4] == 4),4]) # total number of tetraploid individuals this generation
    
    ## population information for current generation [x,y,genetID,ploidy,plantID]
    currentinfo <- recruitmentlist[[2]]
    # carry current population info into next loop iteration
    master[w,,] <- currentinfo
    
    ## output results per generation
    gendata[w,1] <- w # generation number 
    gendata[w,2] <- npolyploids # number of tetraploids
    gendata[w,3] <- n2xsurvivors # number of diploid survivors from last generation
    gendata[w,4] <- n4xsurvivors # number of tetraploid survivors from last generation
    gendata[w,5] <- n2xgenet # number of diploid genets
    gendata[w,6] <- n4xgenet # number of tetraploid genets
    gendata[w,7] <- average4xgenet # tetraploid genet average size 
    gendata[w,8] <- sizepzero # current size of initial tetraploid genet
    gendata[w,9] <- neighboursame2x # average number of neighbours per diploid of the same cytotype
    gendata[w,10] <- neighboursame4x # average number of neighbours per tetraploid of the same cytotype
    gendata[w,11:34] <- prevgenprops # number of all propagule types produced by the previous generation
    gendata[w,35:58] <- prevgenseedgerm # number of those propagule types recruited this generation
    gendata[w,59] <- average2xfreehood # average number of local vacant sites per diploid
    gendata[w,60] <- average4xfreehood # average number of local vacant sites per tetraploid
    
    
    ###### Stop Conditions ###### 
    # simulation will terminate if conditions are met
    
    ## if ug = 0 and all polyploids have been excluded
    if(pstart > 0 && npolyploids == 0 && ug == 0) {
      gendata[(w+1):maxgen,2] <- 0
      break
    }
    
    ## if 4x = N (full tetraploid establishment)
    # disable to continue monitoring population dynamics
    if(npolyploids == N){
      gendata[(w+1):maxgen,2] <- N
      break
    }
    
  }
  
  mainlist<-list(master,gendata)
  
  return(mainlist)	
  
}