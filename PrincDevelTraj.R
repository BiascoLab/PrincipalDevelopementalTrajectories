library("RANN")
library("parallel")

options(bitmapType='cairo')

## functions definition
norm_vec <- function(x) sqrt(sum(x^2))

linMap <- function(x, from, to){
  (x - min(x)) / max(x - min(x)) * (to - from) + from}

average_term=function(index){
  i=index
  v = samples[i,];
  v_neighbors_idx = samples_original_neighbor_idx[i,!is.na(samples_original_neighbor_idx[i,])]
  v_neighbors_dist = samples_original_neighbor_dist[i,!is.na(samples_original_neighbor_idx[i,])]

  neighbor_size = length(v_neighbors_idx);
  if (neighbor_size >= 1){
    t_= as.matrix(copy_original[v_neighbors_idx, ])
    average_weight_ =exp(v_neighbors_dist^2 * iradius16)
    average_weight_sum_=sum(average_weight_)
    if (average_weight_sum_ > 1e-20){
      average_weight_M=matrix((average_weight_/average_weight_sum_),ncol=neighbor_size)
      average_term_ = average_weight_M%*%t_
    }
  }
  return(average_term_)
}

repulsion_term=function(index){
  i=index
  v = samples[i,];
  v_neighbors_idx = samples_samples_neighbor_idx[i,!is.na(samples_samples_neighbor_idx[i,])]
  v_neighbors_dist = samples_samples_neighbor_dist[i,!is.na(samples_samples_neighbor_idx[i,])]
  
  neighbor_size = length(v_neighbors_idx);
  if(neighbor_size>neighbor_size_limit){
    t_= as.matrix(samples[v_neighbors_idx, ])
    v_=matrix(as.numeric(v),ncol = inputDimen,nrow=neighbor_size,byrow = T)
    diff_ = v_ - t_
    fff=prcomp( diff_,scale. = F,center =F)
    rep_weights=diff_%*%fff$rotation[,1]
    
    rep_weights=apply(rep_weights,1, function(x) norm_vec(x))
    repulsion_weight_ = exp((rep_weights^2) * iradiusS16) * (1.0/rep_weights)^1.0;
    repulsion_weight_M=matrix(c(0,(repulsion_weight_[-1]/sum(repulsion_weight_[-1]))),ncol=neighbor_size)
    AtDist=(fff$rotation[,1])%*%t(fff$rotation[,1])%*%t(diff_)
    
    repulsion_term_ =repulsion_weight_M %*% t(AtDist)
  }else{
    repulsion_term_=c(0,0)
  }
  return(repulsion_term_)
}

#############################################################
# Data points input file is expected to have one observation point per line.
fileInput="input.txt" 
coord=read.table(file=fileInput,sep=",",header = F)

inputDimen=dim(coord)[2]
# Rescale input space 
lowBoundRescInp=0
uppBoundRescInp=1

copy_original=apply(coord,2,function(x) {linMap(x,lowBoundRescInp,uppBoundRescInp)} )

#Parameters:
radius= 0.02        # neighborhood size
mu= 0.3             # regularization power
neighbor_size_limit=2 
n.cores=4 # cpu cores to be used 


samples=copy_original

#Compute Average Term
radius2 = radius^2;
iradius16 = (-1)/(2*radius^2);
radiusS =radius;
radiusS2 = radius^2;
iradiusS16 = (-1)/(2*radius^2);
  
iter=1
Delta=1
Eps=0.001
while(Delta>Eps){
  iter=iter+1
  
  samples_original_neighbor_idx = c();
  samples_original_neighbor_dist = c();
  
  
  original_samples=nn2(copy_original, query = samples, searchtype="radius", radius = radius*3,k =dim(copy_original)[1])
  
  samples_original_neighbor_idx=original_samples$nn.idx
  samples_original_neighbor_idx[which(samples_original_neighbor_idx==0)]=NA
  samples_original_neighbor_dist=original_samples$nn.dists
  samples_original_neighbor_dist[which(samples_original_neighbor_idx==0)]=NA
  
  AV=mcmapply(average_term,c(1:dim(samples)[1]),mc.cores=n.cores,SIMPLIFY = F)
  average_term_=do.call(rbind,AV)
  
  ################################  
  samples_samples_neighbor_idx = c();
  samples_samples_neighbor_dist = c();
  
  samples_samples=nn2(samples, query = samples, searchtype="radius", radius = radius*3,k =dim(samples)[1])
  
  samples_samples_neighbor_idx=samples_samples$nn.idx
  samples_samples_neighbor_idx[which(samples_samples_neighbor_idx==0)]=NA
  samples_samples_neighbor_dist=samples_samples$nn.dists
  samples_samples_neighbor_dist[which(samples_samples_neighbor_idx==0)]=NA
  
  RV=mcmapply(repulsion_term,c(1:dim(samples)[1]),mc.cores=n.cores,SIMPLIFY = F)
  repulsion_term_=do.call(rbind,RV)
  samplesPrev=samples
  samples = (average_term_ +(mu * repulsion_term_));

  Delta=(sum((samples-samplesPrev)^2))
  print(Delta)
  if (inputDimen==2){
    plot(copy_original[,1],copy_original[,2],type='p',col=rainbow(12,alpha = 0.2)[2],xlim=c(0,1),ylim=c(0,1),sub=Delta)
    points(samples,type="p",pch=20,col=rainbow(12,alpha = 0.7)[1]);
  }
}

