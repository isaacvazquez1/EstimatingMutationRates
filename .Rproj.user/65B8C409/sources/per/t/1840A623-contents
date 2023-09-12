# LSQ_Estimation
#
# Estimations by direct adjustment of the probability-generating function using
# a least squares method
#
#
# You can learn more about package authoring with RStudio at:
#
#


#Definition of probability-generating functions

h <- function(z){
  (1 - z)^((1 - z) / z)
}

px <- function(z,Am,m1,Amm) {
  (h(z)^Am) * (h(h(z)^m1)^Amm)
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


#Empirical probability function at the j-th growth cycle
epgf_gc_j<-function(td, j, ps1){
  a=c()
  for (i in 1:nrow(td)){
    a=c(a,rbinom(n = 1, size = (td[ , 2*j-1]+td[ , 2*j])[i], prob = ps1))
  }
  epgf<-function(x){
    return(sum(x^a)/nrow(td))
  }
  return(epgf)
}



LSQ_Estimation<-function(File, population_size, GC, dilution, discretization,
                         WT2R_bounds, WT2M_bounds, M2R_bounds){

  #Reading the file
  td <- read.table(File,header = F,sep="")
  n_cycles<-GC
  nn1 <- population_size
  lbm <- 10^WT2R_bounds[1]
  ubm <- 10^WT2R_bounds[2]
  lbm1 <- 10^WT2M_bounds[1]
  ubm1 <- 10^WT2M_bounds[2]
  lbmm <- 10^M2R_bounds[1]
  ubmm <- 10^M2R_bounds[2]
  ps1 <- dilution
  pb1 <- 2^(20)/(2^(33))
  epgf_1=epgf_gc_j(td, 1, ps1)
  for (j in 2:n_cycles){
    assign(paste0('epgf_', j), epgf_gc_j(td, j, ps1))
  }

  zr <- seq(0.4, 0.9, by = 0.1)

  mr <- (lbm * (ubm / lbm)^(0:discretization / discretization)) * nn1
  m1r <- (lbm1 * (ubm1 / lbm1)^(0:discretization / discretization))
  mmr <- (lbmm * (ubmm / lbmm)^(0:discretization / discretization)) * nn1



  tls1 <- list()
  lls <- list()

  tls <- list()
  upper_mne_j <- 1000000

  for (j in 1:n_cycles){
    mne_aux=paste0('mne', j)
    assign(mne_aux, upper_mne_j)
  }

  for (mi in mr) {
    for (m1i in m1r) {
      for (mmi in mmr) {
        for (i in 1:n_cycles){
          assign(paste0('tp', i), 0)
        }
        for(z in zr){
          pg1 <- pgs(z,ps1,mi,m1i,mmi)
          tp1=tp1+(pg1 - epgf_1(z))^2
          for (i in 2:n_cycles){
            pg_aux=paste0('pg', i-1)
            pg.var=get(pg_aux)
            pg_aux_actual=paste0('pg', i)
            assign(pg_aux_actual, pgb(pg.var,mi,m1i,pb1,mmi))
            pg_actual.var=get(pg_aux_actual)
            epgf_aux=paste0('epgf_', i)
            epgf.var=get(epgf_aux)
            tp_aux=paste0('tp',i)
            tp.var=get(tp_aux)
            s=tp.var+(pg_actual.var-epgf.var(z))^2
            assign(tp_aux, s)
          }
        }

        tpt1 <- tp1
        for (i in 2:n_cycles){
          tpt_aux=paste0('tpt', i-1)
          tpt.var=get(tpt_aux)
          tp_aux_actual=paste0('tp',i)
          tp_actual.var=get(tp_aux_actual)
          s=tpt.var+tp_actual.var
          tpt_aux_actual=paste0('tpt', i)
          assign(tpt_aux_actual, s)
        }

        for (i in 1:n_cycles){
          tpt_aux=paste0('tpt', i)
          tpt.var=get(tpt_aux)
          mne_aux=paste0('mne', i)
          mne.var=get(mne_aux)
          ls_aux=paste0('ls',i)
          if (tpt.var < mne.var) {
            assign(mne_aux,tpt.var)
            assign(ls_aux, list(log(mi / nn1), log(m1i), log(mmi / nn1), log(tpt.var)))
          }
        }

        #tls <- c(tls, list(list(log(mi / nn1), log(m1i), log(mmi / nn1)),list(log(tpt1), log(tpt2), log(tpt3))))
      }
    }
  }


  for(i in 1:n_cycles){
    aux=paste0('ls',i)
    aux.var=get(aux)
    esti=paste0('Estimations_GC_',i)
    assign(esti,exp(sapply(aux.var, `[[`, 1)[1:3]))
    esti.var=get(esti)
    names(esti.var)=c(paste0('Wildtype to Mutant | GC', i), paste0('Wildtype to Mutator | GC', i),
                      paste0('Mutator to Mutant | GC', i))
    print(esti.var)
  }

}





