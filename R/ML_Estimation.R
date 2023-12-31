# ML_Estimation
#
# Estimations by the maximum likelihood methods
#
#
# You can learn more about package authoring with RStudio at:
#
#


#Definition of probability-generating functions

h <- function(z) {
  (1 - z)^((1 - z) / z)
}

px <- function(z,Am,m1,Amm) {
  (h(z)^{Am}) * (h(h(z)^m1)^{Amm})
}

psb <- function(z,pb) {
  (1 - pb) + pb * z
}

ps <- function(z,p) {
  (1 - p) + p * z
}

pgb <- function(z,Am,m1,pb,Amm) {
  px(psb(z,pb),Am,m1,Amm)
}

pgs <- function(z,p,Am,m1,Amm){
  px(ps(z,p),Am,m1,Amm)
}

ML_Estimation=function(File, population_size, resampling, resample_size, GC, dilution, discretization,
                       WT2R_bounds, WT2M_bounds, M2R_bounds, alpha){

  #Dilution process: Binomial sampling, with probability ps1, of resistant cell in the j-th growth cycle
  binomial_sampling_GC_j<-function(td,j,ps1){
    a=c()
    for (i in 1:nrow(td)){
      a=c(a,rbinom(n = 1, size = (td[ , 2*j-1]+td[ , 2*j])[i], prob = ps1))
    }
    return(a)
  }

  #Recovery of the k-th probability state of the pmf
  #WT2R  Mutation rate from wildtype to resistant cells
  #WT2M  Mutation rate from wildtype to mutator cells
  #pop   Final population size
  #N     Total number of probability states
  pmf_coeff<-function(k,WT2R,WT2M,M2R,pop,N,ps1){
    sum=0
    for(n in 0:(N-1)){
      #Auxiliary variable to transform the pgf into a characteristic function
      aux_char=exp((2i*pi*n)/N)
      #Auxiliary variable to apply the FFT algorithm
      aux=(-1*2i*pi*n*k)/N
      #FFT
      sum=sum + pgs(aux_char,p=ps1,Am=10^{WT2R}*pop,m1=10^{WT2M},Amm=10^{M2R}*pop)*exp(aux)
    }
    #Real part
    return(Re(sum/N))
  }

  discretization_fun<-function(range, n){
    mu_min=range[1]
    mu_max=range[2]
    a=seq(0,n)
    result=mu_min*(mu_max/mu_min)^{a/n}
    return(result)
  }

  #log-likelihood function
  ll_fun<-function(data,WT2R,WT2M,M2R,pop_size,N,ps1){
    sum(log(pmf_coeff(data+1,WT2R,WT2M,M2R,pop_size,N,ps1)))
  }

  #Reading the file
  td <- read.table(File,header = T,sep="")
  #Binomial sampling of resistant cells in the j-th growth cycle, j=1
  data<-binomial_sampling_GC_j(td, j=GC,p=dilution)

  if(resampling==T){
    data=sample(data, resample_size, replace = T)
  }

  #Maximum value in data
  N=max(data)
  N


  #Sample space discretization
  paramters_WT2R=WT2R_bounds
  WT2R=discretization_fun(paramters_WT2R,discretization)
  paramters_WT2M=WT2M_bounds
  WT2M=discretization_fun(paramters_WT2M,discretization)
  paramters_M2R=M2R_bounds
  M2R=discretization_fun(paramters_M2R,discretization)

  #Initial maximum
  max=ll_fun(data,WT2R = WT2R[1], WT2M = WT2M[1], M2R=M2R[1], population_size, N, dilution)

  powers=c()
  output=c()
  for(i in WT2R){
    for(j in WT2M){
      for(k in M2R){
        aux=ll_fun(data, WT2R = i, WT2M = j, M2R = k,population_size, N, dilution)
        output=rbind(output, c(10^{i}, 10^{j}, 10^{k}, aux))
        if(max<aux){
          max=aux
          powers=c(i,j,k)
        }
      }
    }
  }

  #Maximum value
  max
  #Estimations on log_10 scale
  names(powers)=c('Wildtype to Mutant', 'Wildtype to Mutator', 'Mutator to Mutant mutator')


  df=as.data.frame(output)


  df_confidence<-df[df[,4] >= max-qchisq(1-alpha,3)/2, ]

  CI=matrix(0,3,2)
  row.names(CI)=c('Wildtype to Mutant', 'Wildtype to Mutator', 'Mutator to Mutant mutator')
  colnames(CI)=c('Lower bound','Upper bound')
  CI[1,1]=min(df_confidence$V1)
  CI[1,2]=max(df_confidence$V1)
  CI[2,1]=min(df_confidence$V2)
  CI[2,2]=max(df_confidence$V2)
  CI[3,1]=min(df_confidence$V3)
  CI[3,2]=max(df_confidence$V3)

  out=list(10^{powers}, CI, output)
  names(out)<-c('Estimations', 'ML_Confidence_Interval', 'Log_likelihood')
  return(out)
}

#ML_Estimation(File=FILENAME, population_size=2^33, resampling=F, resample_size=50, GC=2, dilution=0.0001,discretization=25, WT2R_bounds=c(-8, -5), WT2M_bounds=c(-7, -3), M2R_bounds=c(-7, -3), alpha=0.05)
