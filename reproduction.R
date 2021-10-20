###############################
#### Reproduction Function ####
###############################

# Description:
# This function carries out the reproduction step in the simulations, as called by main.R. There are five primary actions: 
# 1) calculate neighbourhood cytotype compositions (per cytotype)
# 2) disperse pollen and choose pollen donors for each ovule produced
# 3) determine sexual offspring crosstypes and viability (see table S1 in the Supplementary Materials)
# 4) produce clonal offspring
# 5) store all clonal and viable sexual offspring into a common pool for use in the recruitment step


reproduction<- function(D,N,paraVect,previnfo)
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
  
  
  ########## Storage Objects ##########
  
  ## previous generation
  prevgen <- matrix(NA,nrow=N,ncol=5) 
  prevgen <- previnfo # Matrix N x 4 [x, y, genetID, ploidy, plantID]
  
  
  
  ########## Average Neighbourhood Composition Prior to Mating ##########
  
  ## calculate the average number of neighbours, and total same cytotype neighbours, for diploids and tetraploids
  area_neighbours<-function(r){
    
    d=2# 5x5 square around focal plant
    areainfo<-rep(0,times=2)
    
    xrange<-c((prevgen[r,1]-d):(prevgen[r,1]+d))
    xrange[which(xrange < 1)] <- xrange[which(xrange < 1)] + D
    xrange[which(xrange > D)] <- xrange[which(xrange > D)] - D
    yrange<-c((prevgen[r,2]-d):(prevgen[r,2]+d))
    yrange[which(yrange < 1)] <- yrange[which(yrange < 1)] + D
    yrange[which(yrange > D)] <- yrange[which(yrange > D)] - D
    arearange<-expand.grid(xrange,yrange)
    
    # which plants in the population are in focal plants area? (prevgen row#)
    sitematch<-function(x){which(prevgen[,1]==arearange[x,1] & prevgen[,2]==arearange[x,2])}
    neighbours<-unlist(sapply(c(1:length(arearange[,1])),FUN=sitematch))
    neighbours<-neighbours[neighbours != r] # remove focal cell
    
    total<-length(neighbours)
    sameploidy<-length(which(prevgen[neighbours,4] == prevgen[r,4]))
    areainfo<-c(total,sameploidy)
  }
  
  ## apply neighbourhood function to all individuals in population
  allneighbours<-t(sapply(c(1:N),FUN=area_neighbours))
  avg4xtotalarea<-mean(allneighbours[which(prevgen[,4] == 4),1],na.rm=TRUE) # total neighbours per tetraploid
  avg4xsamearea<-mean(allneighbours[which(prevgen[,4] == 4),2],na.rm=TRUE) # same cytotype neighbours per tetraploid
  avg2xtotalarea<-mean(allneighbours[which(prevgen[,4] == 2),1],na.rm=TRUE) # total neighbours per diploid
  avg2xsamearea<-mean(allneighbours[which(prevgen[,4] == 2),2],na.rm=TRUE) # same cytotype neighbours per diploid
  
  
  
  ########## Pollen Dispersal and Ovule Fertilization ##########
  
  ## object to store ovule information
  ovuleprod <- array(NA,dim=c(N,ovnumb,2)) # seeds (col) per individual (row), [pollen donor ID, seed type]
  
  ## calculate distances between all individuals, wrap around 2D toroidal boundaries
  distarray<-array(NA,dim=c(N,N,2))
  distx<-function(x){
    nowrap<-abs(prevgen[x,1]-prevgen[,1])
    wrapcond<-which(nowrap > D/2,arr.ind=TRUE)
    nowrap[wrapcond]<-(D-nowrap)[wrapcond]
    nowrap}
  disty<-function(x){
    nowrap<-abs(prevgen[x,2]-prevgen[,2])
    wrapcond<-which(nowrap > D/2,arr.ind=TRUE)
    nowrap[wrapcond]<-(D-nowrap)[wrapcond]
    nowrap}
  distarray[,,1]<-t(sapply(c(1:length(prevgen[,1])),FUN=distx))
  distarray[,,2]<-t(sapply(c(1:length(prevgen[,1])),FUN=disty))
  # euclidean distances between all individuals, note that the max possible distance after wrapping is sqrt(2)*(D/2)
  distmat<-matrix(sqrt((distarray[,,1]^2)+(distarray[,,2]^2)) ,nrow=length(prevgen[,1]),ncol=length(prevgen[,1]))
  
  ## calculate pollen density from all plants at site of each individual (including self)
  beta = 1 # shape parameter. 1 = exponential, such that D_G(0) = 0 is highest density.
  alpha = beta/dp # rate parameter = beta/avg
  # gamma probability distribution
  D_G <- function(x){((alpha^beta)*(x^(beta-1))*exp(-alpha*x))/gamma(beta)}
  # pollen density from every other individual in pop
  probmat<-apply(distmat,MARGIN=1,FUN=D_G) 
  
  ## choose pollen donors (by within generation plantID) for all ovules weighted by probability density
  # pollen produced per individual is assumed to be unlimited
  donor<-function(x){sample(prevgen[,5],size=ovnumb,prob=probmat[x,],replace=TRUE)}
  # successful pollen donors for each ovule per individual
  ovuleprod[,,1]<-t(sapply(c(1:N),FUN=donor)) 
  
  
  
  ########## Cross-type and Seed Viability ##########
  
  ## Determine cross-type of all fertilized ovules
  
  # Diploid Mothers (dam + sire) (probabilty)
  # 1 - SHOOT inviable self, 1n + 1n, (1-ug)(1-ug)k2
  # 2 - SHOOT viable self, 1n + 1n, (1-ug)(1-ug)(1-k2)
  # 3 - SHOOT 3x inviable self, 1n + 2n OR 2n + 1n, 2*ug(1-ug)
  # 4 - SHOOT inviable self, 2n + 2n, ug*ug*k2
  # 5 - SHOOT viable self, 2n + 2n, ug*ug*(1-k2)
  # 6 - RAMET inviable self, 1n + 1n, (1-ug)(1-ug)k2
  # 7 - RAMET viable self, 1n + 1n, (1-ug)(1-ug)(1-k2)
  # 8 - RAMET 3x inviable self, 1n + 2n OR 2n + 1n, 2*ug(1-ug)
  # 9 - RAMET inviable self, 2n + 2n, ug*ug*k2
  # 10 - RAMET viable self, 2n + 2n, ug*ug*(1-k2)
  # 11 - 3x between, 1n + 2n, (1-ug) 
  # 12 - viable 4x between, 2n + 2n, ug
  # 13 - viable within, 1n + 1n, (1-ug)(1-ug)
  # 14 - 3x within, 1n + 2n OR 2n + 1n, 2*ug(1-ug)
  # 15 - viable within, 2n + 2n, ug*ug
  
  # Tetraploid Mothers (dam + sire) (probabilty)
  # 16 - SHOOT inviable self, 2n + 2n, k4
  # 17 - SHOOT viable self, 2n + 2n, (1-k4)
  # 18 - RAMET inviable self, 2n + 2n, k4
  # 19 - RAMET viable self, 2n + 2n, (1-k4)
  # 20 - 3x between,, 2n + 1n, (1-ug)
  # 21 - viable between, 2n + 2n, ug
  # 22 - viable within, 2n + 2n, 1
  
  ## function takes info on ovules and pollen donors per individual, 
  # and determines the crosstype of each seedling according to the probabilities above (see table S1 in Supplementary Materials)
  viability<-function(x){
    
    # temporary crosstype storage
    fert_type<-rep(0,times=ovnumb)
    
    ## conditions for selfing:
    # within shoot, which pollen donors have the same genetID & plantID as the ovule producer
    selfshootcond<-which(prevgen[ovuleprod[x,,1],3] == prevgen[x,3] & prevgen[ovuleprod[x,,1],5] == prevgen[x,5]) 
    # between ramets, which pollen donors have the same genetID but NOT plant ID as the ovule producer
    selframetcond<-which(prevgen[ovuleprod[x,,1],3] == prevgen[x,3] & prevgen[ovuleprod[x,,1],5] != prevgen[x,5]) 
    
    ## conditions for outcrossing:
    donorinds<-sapply(c(1:ovnumb), function(r) which(prevgen[,5] == ovuleprod[x,r,1])) # row number of pollen donor in pinfo (plantID)
    # between-cytotype mating, and NOT the same genet ID
    betweencond<-which(prevgen[ovuleprod[x,,1],3] != prevgen[x,3] & prevgen[x,4] != prevgen[donorinds,4])
    # within-cytotype mating, and NOT the same genet ID 
    withincond<-which(prevgen[ovuleprod[x,,1],3] != prevgen[x,3] & prevgen[x,4] == prevgen[donorinds,4])
    
    # diploid mothers
    if(prevgen[x,4] == 2){
      # within shoot selfing, between ramet selfing, between cytotype outcrossing, within cytotype outcrossing
      fert_type[selfshootcond]<-sample(c(1,2,3,4,5),size=length(selfshootcond),prob=c(ks2*(1-ug)^2,(1-ks2)*(1-ug)^2,2*ug*(1-ug),ks2*ug*ug,(1-ks2)*ug*ug),replace=TRUE)
      fert_type[selframetcond]<-sample(c(6,7,8,9,10),size=length(selframetcond),prob=c(ks2*(1-ug)^2,(1-ks2)*(1-ug)^2,2*ug*(1-ug),ks2*ug*ug,(1-ks2)*ug*ug),replace=TRUE)
      fert_type[betweencond] <- sample(c(11,12),size=length(betweencond),prob=c(1-ug,ug),replace=TRUE) 
      fert_type[withincond] <- sample(c(13,14,15),size=length(withincond),prob=c((1-ug)^2,2*ug*(1-ug),ug*ug),replace=TRUE)
    } 
    
    # polyploid mothers  
    if(prevgen[x,4] == 4){
      # within shoot selfing, between ramet selfing, between cytotype outcrossing, within cytotype outcrossing
      fert_type[selfshootcond]<-sample(c(16,17),size=length(selfshootcond),prob=c(ks4,1-ks4),replace=TRUE) 
      fert_type[selframetcond]<-sample(c(18,19),size=length(selframetcond),prob=c(ks4,1-ks4),replace=TRUE) 
      fert_type[betweencond] <- sample(c(20,21),size=length(betweencond),prob=c(1-ug,ug),replace=TRUE)  
      fert_type[withincond] <- 22 
    }
    
    fert_type
    
  }
  
  # get cross-types for each ovule per individual
  ovuleprod[,,2]<-t(sapply(c(1:N),FUN=viability)) 
  
  
  
  ########## Summarize Cross-Types ##########
  
  ## types of seeds made
  seedtypes <- rep(0,times=24) # see above for seed type codes, 23 + 24 = diploid and tetraploid ramets
  names(seedtypes)<-as.character(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"))
  seedsum<-table(ovuleprod[,,2])
  seedsum2<-as.vector(seedsum)
  names(seedsum2)<-names(seedsum)
  seedtypes[names(seedsum2)]<-seedsum2
  
  ## list of viable seeds only for recruitment
  viableseeds<-c(2,5,7,10,12,13,15,17,19,21,22)
  viablecond<-which(ovuleprod[,,2] %in% viableseeds,arr.ind=TRUE)
  viablecondinds<-which(ovuleprod[,,2] == 2 | ovuleprod[,,2] == 5 | ovuleprod[,,2] == 7 | ovuleprod[,,2] == 10 | ovuleprod[,,2] == 12 | ovuleprod[,,2] == 13 | ovuleprod[,,2] == 15 | ovuleprod[,,2] == 17 | ovuleprod[,,2] == 19 | ovuleprod[,,2] == 21 | ovuleprod[,,2] == 22, arr.ind=TRUE)
  
  ## info on viable seedling recruits (source locations, parents, ids, etc)
  # goes in to seedling recruitment pool
  seedrecruits<-matrix(NA,nrow=length(viablecondinds[,1]),ncol=8) # source x, source y, ploidy, sv, mother ID, pollen donor ID, seed type
  seedrecruits[,1]<-prevgen[,1][viablecondinds[,1]] # x
  seedrecruits[,2]<-prevgen[,2][viablecondinds[,1]] # y
  seedrecruits[,3]<-prevgen[,4][viablecondinds[,1]] # ploidy
  seedrecruits[,4]<-prevgen[,3][viablecondinds[,1]] # mother (genetID)
  seedrecruits[,5]<-prevgen[ovuleprod[,,1][viablecond],3] # father (genetID)
  seedrecruits[,6]<-viablecondinds[,1] # mother (plantID)
  seedrecruits[,7]<-ovuleprod[,,1][viablecond] # pollen donor (plantID)
  seedrecruits[,8]<-ovuleprod[,,2][viablecond] # seed type
  
  ## fix ploidy from viable crosses involving a diploid UG ovule: 5,10,12,15
  seedrecruits[which(seedrecruits[,8] == 5 | seedrecruits[,8] == 10 | seedrecruits[,8] == 12 | seedrecruits[,8] == 15),3] <- 4
  
  
  
  ########## Clonal Reproduction ##########
  
  ## each generation, each 2x can make c2 ramets, each 4x can make c4 ramets
  # these go into a new ramet recruit pool 
  
  # make a list of individuals making the clones (ie, prevgen), format like seed recruits
  clones<-matrix(NA,nrow=length(prevgen[,1]),ncol=9)
  clones[,1:2]<-prevgen[,1:2] # xy coordinates of ramet sources
  clones[,3]<-prevgen[,4] # ploidy
  clones[,4]<-prevgen[,3] # mother genetID
  clones[,5]<-prevgen[,3] # "father" genetID
  clones[,6]<-prevgen[,5] # mother plantID
  clones[,7]<-prevgen[,5] # "father" plantID
  clones[,8]<-ifelse(prevgen[,4] == 2, 23, 24) # 23,24 is "seedtype=ramet" for 2x,4x
  clones[,9]<-ifelse(prevgen[,4] == 2, c2, c4) # number of ramets
  
  # repeat rows by last column (number of ramets being made per individual)
  clonesdf<-as.data.frame(clones)
  rametrecruits<-as.matrix(clonesdf %>% uncount(V9))
  
  
  ## add ramets to seed counts to get all propagule types produced this generation
  propaguletypes<-seedtypes
  propaguletypes[23]<-length(which(rametrecruits[,8] == 23))
  propaguletypes[24]<-length(which(rametrecruits[,8] == 24))
  
  
  
  ########## Pass Information To Main ##########
  
  reproduction_list<-list(propaguletypes,seedrecruits,rametrecruits,avg4xtotalarea,avg4xsamearea,avg2xtotalarea,avg2xsamearea)
  
  return(reproduction_list)  
  
} 