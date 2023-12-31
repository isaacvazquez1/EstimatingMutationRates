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


ML_Bivariate=function(File, population_size, resampling, resample_size, GC, dilution, discretization,
                       WT2R_bounds, WT2M_bounds, strength, alpha){
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
  pmf_coeff<-function(k,WT2R,WT2M,M2R,pop,N,ps1,strength){
    sum=0
    for(n in 0:(N-1)){
      #Auxiliary variable to transform the pgf into a characteristic function
      aux_char=exp((2i*pi*n)/N)
      #Auxiliary variable to apply the FFT algorithm
      aux=(-1*2i*pi*n*k)/N
      #FFT
      sum=sum + pgs(aux_char,p=ps1,Am=10^{WT2R}*pop,m1=10^{WT2M},Amm=strength*10^{WT2R}*pop)*exp(aux)
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
  ll_fun<-function(data,WT2R,WT2M,M2R,pop_size,N,ps1,strength){
    sum(log(pmf_coeff(data+1,WT2R,WT2M,M2R,pop_size,N,ps1,strength)))
  }
  #Reading the file
  td <- read.table(File,header = F,sep="")
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

  #Initial maximum
  max=ll_fun(data, WT2R = WT2R[1], WT2M = WT2M[1], M2R=WT2R[1], population_size, N, dilution, strength)

  powers=c()
  output=c()
  for(i in WT2R){
    for(j in WT2M){
      aux=ll_fun(data, WT2R = i, WT2M = j, M2R = i,population_size, N, dilution, strength)
      output=rbind(output, c(10^{i}, 10^{j}, strength*10^{i}, aux))
      if(max<aux){
        max=aux
        powers=c(10^{i}, 10^{j}, strength*10^{i})
      }
    }
  }

  #Maximum value
  max

  names(powers)=c('Wildtype to Mutant', 'Wildtype to Mutator', 'Mutator to Mutant mutator')
  df=as.data.frame(output)


  df_confidence<-df[df[,4] >= max-qchisq(1-alpha,2)/2, ]

  CI=matrix(0,3,2)
  row.names(CI)=c('Wildtype to Mutant', 'Wildtype to Mutator', 'Mutator to Mutant mutator')
  colnames(CI)=c('Lower bound','Upper bound')
  CI[1,1]=min(df_confidence$V1)
  CI[1,2]=max(df_confidence$V1)
  CI[2,1]=min(df_confidence$V2)
  CI[2,2]=max(df_confidence$V2)
  CI[3,1]=min(df_confidence$V3)
  CI[3,2]=max(df_confidence$V3)

  out=list(powers, CI, output)
  names(out)<-c('Estimations', 'ML_Confidence_Interval', 'Log_likelihood')
  return(out)
}


#(File, population_size, resampling, resample_size, GC, dilution, discretization, WT2R_bounds, WT2M_bounds, strength)
#FILENAME="QM-Mutants.txt"
#ML_Bivariate(File=FILENAME, population_size=2^33, resampling = F, resample_size=0, GC=1, dilution=0.0001, discretization=25, WT2R_bounds=c(-8, -5), WT2M_bounds=c(-7, -3), strength= 1000, alpha=0.05)
