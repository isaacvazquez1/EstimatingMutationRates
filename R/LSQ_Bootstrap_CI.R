# LSQ_Bootstrap_CI
#
#
#
#
# You can learn more about package authoring with RStudio at:
#
#
LSQ_Bootstrap_CI<-function(Replicates, alpha, FreeCores, File, population_size, GC,
                           dilution,discretization, WT2R_bounds, WT2M_bounds, M2R_bounds){
  library(parallel)

  cl <- makeCluster(detectCores() - FreeCores)
  cl

  library(doParallel)
  registerDoParallel(cl)

  rep=Replicates
  output<-foreach(k=1:rep)%dopar%{
    EstimatingMutationRates::LSQ_Estimation(File, population_size, GC, dilution, discretization, WT2R_bounds, WT2M_bounds, M2R_bounds)
  }

  result_1=c()
  for(j in 1:rep){
    result_1=rbind(result_1,output[[j]])
  }

  #Estimation
  estimate_WT2R=mean(result_1[,1])
  estimate_WT2M=mean(result_1[,2])
  estimate_M2R=mean(result_1[,3])

  confidence_intervals=matrix(0,3,3)
  confidence_intervals[1,1]=estimate_WT2R
  confidence_intervals[2,1]=estimate_WT2M
  confidence_intervals[3,1]=estimate_M2R

  #95% confidence region
  alpha
  for(i in 1:ncol(result_1)){
    aux=result_1[,i]
    aux=sort(aux)
    lower=round(alpha/2*(nrow(result_1)+1))
    upper=round((1-alpha/2)*(nrow(result_1)+1))
    if(lower==0){
      lower=1
    }
    if(upper>nrow(result_1)){
      upper=nrow(result_1)
    }
    confidence_intervals[i,2]=aux[lower]
    confidence_intervals[i,3]=aux[upper]
  }


  colnames(confidence_intervals)=c('Estimate', 'Lower bound', 'Upper bound')
  rownames(confidence_intervals)=c(paste('Wildtype to Mutant | GC',GC), paste('Wildtype to Mutator | CG',GC), paste('Mutator to Mutant mutator  | GC',GC))
  return(confidence_intervals)
}
