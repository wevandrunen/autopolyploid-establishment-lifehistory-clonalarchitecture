##############################
#### Recruitment Function ####
##############################

# Description:
# This function carries out the survival and recruitment steps, as called by main.R. There are five primary actions:
# 1) determine which individuals survive or die, and where in the population empty sites are left for offspring recruitment
# 2) calculate the average number of empty sites around each individual (per cytotype)
# 3) disperse available seeds and ramets into the population
# 4) determine which offspring are recruited to the empty sites
# 5) create new generation, composed of surviving individuals and newly recruited offspring


recruitment<-function(D,N,paraVect,previnfo,seedpool,rametpool,pzero)
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
  prevgen <- previnfo # Matrix N x 4 [x, y, genetID, ploidy,plantID]
  
  
  
  ######## Survival ##########
  
  ## individuals survive to next year with probability sv
  # determine ramet death (1 = survive, 0 = die)
  survival<-function(x){
    if(prevgen[x,4] == 2){
      # diploid survival
      sample(c(0,1),size=1,prob=c(1-sv2,sv2)) 
    } else {
      # tetraploid survival
      sample(c(0,1),size=1,prob=c(1-sv4,sv4))
    }
  }
  survivors<-unlist(lapply(c(1:N),FUN=survival)) # vector of individuals that live (1) or die (0)
  # which individuals live?
  prevgensurvive<-matrix(prevgen[which(survivors == 1),],nrow=length(which(survivors == 1)),ncol=5)
  # how many of each ploidy survive?
  numb2xsurvive <- length(which(prevgensurvive[,4] == 2))
  numb4xsurvive <- length(which(prevgensurvive[,4] == 4))
  
  
  
  ########## Identify New Recruit Sites ##########
  
  ## seed germination sites, opened up by the dead
  nsites<-N-length(prevgensurvive[,1])
  germsites<-matrix(prevgen[which(survivors == 0), c(1:2)],nrow=nsites,ncol=2)
  # store recruitment sites for current generation population information
  currentgen <- matrix(NA,nrow=length(germsites[,1]),ncol=5) 
  currentgen[,c(1:2)]<-germsites
  
  
  
  ########## Average Free Sites in Neighbourhood ##########
  
  ## calculate the average number of vacant sites around diploid vs. tetraploid individuals
  freearea<-function(r){
    
    d=2# 5x5 square around focal plant
    areainfo<-rep(0,times=2)
    
    xrange<-c((prevgen[r,1]-d):(prevgen[r,1]+d))
    xrange[which(xrange < 1)] <- xrange[which(xrange < 1)] + D
    xrange[which(xrange > D)] <- xrange[which(xrange > D)] - D
    yrange<-c((prevgen[r,2]-d):(prevgen[r,2]+d))
    yrange[which(yrange < 1)] <- yrange[which(yrange < 1)] + D
    yrange[which(yrange > D)] <- yrange[which(yrange > D)] - D
    arearange<-expand.grid(xrange,yrange)
    
    # which spots around each plant emptied? (currentgen row#) NOTE: this CAN include the focal plant
    sitematch<-function(x){which(currentgen[,1]==arearange[x,1] & currentgen[,2]==arearange[x,2])}
    freeneighbours<-unlist(sapply(c(1:length(arearange[,1])),FUN=sitematch))
    
    total<-length(freeneighbours)
    total
  }
  
  allfreeneighbours<-t(sapply(c(1:N),FUN=freearea))
  avg4xfreearea<-mean(allfreeneighbours[which(prevgen[,4] == 4)],na.rm=TRUE) # average vacant sites around each tetraploid
  avg2xfreearea<-mean(allfreeneighbours[which(prevgen[,4] == 2)],na.rm=TRUE) # average vacant sites around each diploid
  
  
  
  ########## Calculate Seed Probability Densities at New Recruit Sites ##########
  
  ## calculate distances between recruitment sites and seed pool source locations, wrap around 2D toroidal boundaries,
  # and calculate the probability density of each seed in the pool at each recruitment site
  
  if(length(seedpool[,1]) == 0){
    # if there are no viable seeds
    germprobmat<-matrix(nrow=length(currentgen[,1]),ncol=0)
  } else {
    if(length(currentgen[,1]) == 0){
      # if there are no sites to fill
      germprobmat<-matrix(nrow=0,ncol=length(seedpool[,1]))
    } else {
      germdistarray<-array(NA,dim=c(length(currentgen[,1]),length(seedpool[,1]),2)) # rows = new germ sites, cols = seed , [x,y]
      germdistx<-function(x){
        nowrap<-abs(currentgen[x,1]-seedpool[,1])
        wrapcond<-which(nowrap > D/2,arr.ind=TRUE)
        nowrap[wrapcond]<-(D-nowrap)[wrapcond]
        nowrap}
      germdisty<-function(x){
        nowrap<-abs(currentgen[x,2]-seedpool[,2])
        wrapcond<-which(nowrap > D/2,arr.ind=TRUE)
        nowrap[wrapcond]<-(D-nowrap)[wrapcond]
        nowrap}
      germdistarray[,,1]<-t(sapply(c(1:length(currentgen[,1])),FUN=germdistx))
      germdistarray[,,2]<-t(sapply(c(1:length(currentgen[,1])),FUN=germdisty))
      # each row is a new germ site, each column is a potential seed:
      germdistmat<-matrix(sqrt((germdistarray[,,1]^2)+(germdistarray[,,2]^2)) ,nrow=length(currentgen[,1]),ncol=length(seedpool[,1]))
      
      ## calculate seed probability at recruitment site for each potential seed
      beta=1 # shape parameter. 1 = exponential such that D_G(0) = 0 is highest density.
      alpha=beta/ds # rate parameter = beta/average
      # gamma prob dist
      D_G <- function(x){((alpha^beta)*(x^(beta-1))*exp(-alpha*x))/gamma(beta)}
      
      germprobmat<-matrix(apply(germdistmat,MARGIN=2,FUN=D_G),nrow=length(currentgen[,1]),ncol=length(seedpool[,1]))
      
    }
  }
  
  
  
  ########## Calculate Ramet Probability Densities at New Recruit Sites ##########
  
  ## calculate distances between recruitment sites and ramet pool source locations, wrap around 2D toroidal boundaries,
  # and calculate the probability density of each ramet in the pool at each recruitment site
  
  if(length(rametpool[,1]) == 0){
    # if no ramets are produced
    ramprobmat<-matrix(nrow=length(currentgen[,1]),ncol=0)
  } else {
    if(length(currentgen[,1]) == 0){
      # if there are no sites to fill
      germprobmat<-matrix(nrow=0,ncol=length(rametpool[,1]))
    } else {
      ramdistarray<-array(NA,dim=c(length(currentgen[,1]),length(rametpool[,1]),2)) # rows = new sites, cols = ramets , [x,y]
      ramdistx<-function(x){
        nowrap<-abs(currentgen[x,1]-rametpool[,1])
        wrapcond<-which(nowrap > D/2,arr.ind=TRUE)
        nowrap[wrapcond]<-(D-nowrap)[wrapcond]
        nowrap}
      ramdisty<-function(x){
        nowrap<-abs(currentgen[x,2]-rametpool[,2])
        wrapcond<-which(nowrap > D/2,arr.ind=TRUE)
        nowrap[wrapcond]<-(D-nowrap)[wrapcond]
        nowrap}
      ramdistarray[,,1]<-t(sapply(c(1:length(currentgen[,1])),FUN=ramdistx))
      ramdistarray[,,2]<-t(sapply(c(1:length(currentgen[,1])),FUN=ramdisty))
      # each row is a new site, each column is a potential ramet
      ramdistmat<-matrix(sqrt((ramdistarray[,,1]^2)+(ramdistarray[,,2]^2)) ,nrow=length(currentgen[,1]),ncol=length(rametpool[,1]))
      
      ## calculate ramet probability density at recruitment site for each potential ramet
      # normal density function:
      sigmasq <- 0.15 # variance for distribution
      D_N <- function(x){
        if(rametpool[x,3] == 2){avg=dc2}
        if(rametpool[x,3] == 4){avg=dc4}
        exp((-(ramdistmat[,x]-avg)^2)/2*(sigmasq))/sqrt(2*pi*(sigmasq))}
      
      # exponential density function:
      # beta=1
      # D_N <- function(x){
      #   if(rametpool[x,3] == 2){alpha=beta/dc2}
      #   if(rametpool[x,3] == 4){alpha=beta/dc4}
      #   ((alpha^beta)*(x^(beta-1))*exp(-alpha*x))/gamma(beta)}
      
      ramprobmat<-matrix(sapply(c(1:length(rametpool[,1])),FUN=D_N),nrow=length(currentgen[,1]),ncol=length(rametpool[,1]))
    }
  }
  
  
  
  ########## Recruit Seeds and Ramets ##########
  
  ## add pools and probability matrixes for seeds and ramets together
  totalpool<-rbind(seedpool,rametpool) # columns = [x,y,ploidy,mother,father,seed type]
  totalprobmat<-cbind(germprobmat,ramprobmat)
  
  ## choose recruits for new sites
  if(length(currentgen[,1]) == 0){
    # if there are no sites to fill
    successprops<-matrix(nrow=0,ncol=8) 
  } else {
    totalpoolID<-c(1:length(totalpool[,1]))
    successpropnumber<-c()
    # once a seed or ramet is chosen for a site, it is removed from the recruitment pool
    for(i in 1:length(currentgen[,1])){
      if(i == 1){
        successpropnumber[i]<-sample(totalpoolID,size=1,prob=totalprobmat[i,])
      } else {
        successpropnumber[i]<-sample(totalpoolID2,size=1,prob=totalprobmat[i,-c(successpropnumber)])
      }
      totalpoolID2<-totalpoolID[-c(successpropnumber)]
    }
    successprops<-matrix(totalpool[successpropnumber,],nrow=length(currentgen[,1]),ncol=8) 
    # columns = [x,y,ploidy,mother genet,father genet,mother plantid, father plantid,propagule type]
  }
  
  
  
  ########## New Recruit Information ##########
  
  ## add ploidy of recruits to current generation population information
  currentgen[,4]<-successprops[,3] 
  
  ## seeds get a new genetid that is currently unused, ramets retain the genetid of their parent
  # ramet ids
  currentgen[successprops[,8] == 23 | successprops[,8] == 24,3] <- successprops[successprops[,8] == 23 | successprops[,8] == 24,4]
  # seed ids (+ 1 to account for pzero)
  totalids<-c(1:(N+1))
  # taken ids = pzero + survivor ids + ids of any ramet recruits
  freeids<-setdiff(totalids,unique(c(pzero,prevgensurvive[,3],successprops[successprops[,8] == 23 | successprops[,8] == 24,4])))
  currentgen[successprops[,8] != 23 & successprops[,8] != 24,3] <- sample(freeids,size=length(which(successprops[,8] != 23 & successprops[,8] != 24)),replace=FALSE)
  
  
  
  ########## Add Recruits to Survivors ##########
  
  ## append current generation to previous generation survivors
  currentgen_all <- rbind(currentgen,prevgensurvive) # columns = [x,y,genetID,ploidy,plantID]
  # assign new within generation plantids
  currentgen_all[,5]<-c(1:N)
  
  ## how many genets of each ploidy are there?
  numb2xgenet <- length(unique(currentgen_all[currentgen_all[,4] == 2,3]))
  numb4xgenet <- length(unique(currentgen_all[currentgen_all[,4] == 4,3]))
  
  ## for each ploidy, how big is the largest genet?
  largest2xgenet <- sort(table(currentgen_all[currentgen_all[,4] == 2,3]),decreasing=TRUE)[1]
  largest4xgenet <- sort(table(currentgen_all[currentgen_all[,4] == 4,3]),decreasing=TRUE)[1]
  # average tetraploid genet size
  average4xgenet <- mean(table(currentgen_all[currentgen_all[,4] == 4,3])) 
  # current size of initial tetraploid genet
  sizepzero <- length(currentgen_all[which(currentgen_all[,3] == pzero),1])
  
  
  
  ########## Summarize Recruit Cross-types/Ramets ##########
  
  ## record recruit types
  proptypes <- rep(0,times=24) # see Reproduction function for propagule type codes
  names(proptypes)<-as.character(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"))
  propsum<-table(successprops[,8])
  propsum2<-as.vector(propsum)
  names(propsum2)<-names(propsum)
  proptypes[names(propsum2)]<-propsum2
  
  
  
  ########## Pass Information Back to Main ##########
  
  recruitment_list<-list(proptypes,currentgen_all,numb2xsurvive,numb4xsurvive,numb2xgenet,numb4xgenet,largest2xgenet,largest4xgenet,average4xgenet,sizepzero,avg2xfreearea,avg4xfreearea)
  
  return(recruitment_list)  
  
}
