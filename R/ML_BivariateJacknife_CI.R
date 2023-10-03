# ML_BivariateJackknife_CI
#
#
#
#
# You can learn more about package authoring with RStudio at:
#
#
ML_BivariateJackknife_CI<-function(Replicates, alpha, FreeCores, File, population_size, resample_size, GC, dilution_lowerbound, dilution_upperbound, discretization,
                          WT2R_bounds, WT2M_bounds, strength){

  ## Bootstrap procedure with replacement
  library(parallel)

  cl <- makeCluster(detectCores() - FreeCores)
  cl

  library(doParallel)
  registerDoParallel(cl)
  rep=Replicates
  output<-foreach(k=1:rep)%dopar%{
    EstimatingMutationRates::ML_Bivariate(File, population_size, resampling = T, resample_size, GC, dilution=runif(1,dilution_lowerbound,dilution_upperbound),
                                           discretization, WT2R_bounds, WT2M_bounds, strength, alpha=0.05)
  }

  result_1=c()
  for(j in 1:rep){
    result_1=rbind(result_1,output[[j]][[1]])
  }


  #Estimation


  estimate_WT2R=mean(result_1[,1])
  estimate_WT2M=mean(result_1[,2])
  estimate_M2R=mean(result_1[,3])



  confidence_intervals=matrix(0,3,3)
  confidence_intervals[1,1]=estimate_WT2R
  confidence_intervals[2,1]=estimate_WT2M
  confidence_intervals[3,1]=estimate_M2R


  #(1-alpha)100% confidence region


  aux_WT2R=qt(p= alpha/2, df=rep-1, lower.tail = F)*var(result_1[,1])/sqrt(rep)
  aux_WT2M=qt(p= alpha/2, df=rep-1, lower.tail = F)*var(result_1[,2])/sqrt(rep)
  aux_M2R=qt(p= alpha/2, df=rep-1, lower.tail = F)*var(result_1[,3])/sqrt(rep)


  confidence_intervals[1,2]=estimate_WT2R-aux_WT2R
  confidence_intervals[1,3]=estimate_WT2R+aux_WT2R
  confidence_intervals[2,2]=estimate_WT2M-aux_WT2M
  confidence_intervals[2,3]=estimate_WT2M+aux_WT2M
  confidence_intervals[3,2]=estimate_M2R-aux_M2R
  confidence_intervals[3,3]=estimate_M2R+aux_M2R


  colnames(confidence_intervals)=c('Estimate', 'Lower bound', 'Upper bound')
  rownames(confidence_intervals)=c(paste('Wildtype to Mutant | GC',GC), paste('Wildtype to Mutator | CG',GC), paste('Mutator to Mutant mutator  | GC',GC))
  return(confidence_intervals)
}


#FILENAME="QM-Mutants.txt"
#ML_BivariateJackknife_CI(Replicates=2, alpha=0.05, FreeCores=5, File=FILENAME, population_size=2^33, resample_size=50, GC=1, dilution_lowerbound=0.00005, dilution_upperbound=2^20/2^33, discretization=25, WT2R_bounds=c(-8, -5), WT2M_bounds=c(-7, -3), strength=500)
