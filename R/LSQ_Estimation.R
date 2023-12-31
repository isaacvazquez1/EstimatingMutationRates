# LSQ_Estimation
#
# Estimations by direct adjustment of the probability-generating function using
# a least squares method
#
#
# You can learn more about package authoring with RStudio at:
#
#


LSQ_Estimation=function(File, population_size, GC, dilution, discretization, WT2R_bounds, WT2M_bounds, M2R_bounds){

  #Reading the file
  td <- read.table(File,header = F,sep="")
  n_cycles<-GC
  nn1 <- population_size
  ni <- discretization
  lbm <- 10^WT2R_bounds[1]
  ubm <- 10^WT2R_bounds[2]
  lbm1 <- 10^WT2M_bounds[1]
  ubm1 <- 10^WT2M_bounds[2]
  lbmm <- 10^M2R_bounds[1]
  ubmm <- 10^M2R_bounds[2]
  ps1 <- dilution
  pb1 <- 2^(20)/(2^(33))

  h <- function(z) {
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




  #resample_size=50
  #td=sample(td, resample_size, replace = T)
  #head(td)


  epgf_gc_j<-function(j, td, ps1){
    a=c()
    for (i in 1:nrow(td)){
      a=c(a,rbinom(n = 1, size = (td[ , 2*j-1]+td[ , 2*j])[i], prob = ps1))
      #a=sample(a, size=50, replace = T)
    }

    epgf<-function(x){
      return(sum(x^a)/nrow(td))
    }
    return(epgf)
  }
  zr <- seq(0.4, 0.9, by = 0.1)

  mr <- (lbm * (ubm / lbm)^(0:ni / ni)) * nn1

  m1r <- (lbm1 * (ubm1 / lbm1)^(0:ni / ni))

  mmr <- (lbmm * (ubmm / lbmm)^(0:ni / ni)) * nn1



  for (i in 1:GC){
    aux=paste0('tls', i)
    assign(aux, c())
  }
  #tls1 <- c()
  #tls2 <- c()
  #tls3 <- c()


  lls <- c()
  tls <- c()

  for (i in 1:GC){
    aux=paste0('epgf', i)
    assign(aux, epgf_gc_j(i, td, ps1))
  }
  #epgf1 <- epgf_gc_j(1)
  #epgf2 <- epgf_gc_j(2)
  #epgf3 <- epgf_gc_j(3)

  for (i in 1:GC){
    aux=paste0('mne', i)
    assign(aux, 1000000)
  }
  #mne1 <- 1000000
  #mne2 <- 1000000
  #mne3 <- 1000000


  for (mi in mr) {
    for (m1i in m1r) {
      for (mmi in mmr) {
        #tp1=0
        #tp2=0
        #tp3=0
        for (i in 1:GC){
          aux=paste0('tp', i)
          assign(aux, 0)
        }
        for(z in zr){
          #pgs(z,p,Am,m1,Amm)             #Parameters to be used
          pg1 <- pgs(z,ps1,mi,m1i,mmi)
          #pgb(z,Am,m1,pb,Amm)            #Parameters to be used
          #pg2 <- pgb(pg1,mi,m1i,pb1,mmi)
          #pg3 <- pgb(pg2,mi,m1i,pb1,mmi)
          if(GC>1){
            for (i in 2:GC){
              aux=paste0('pg', i)
              past=paste0('pg', i-1)
              past.var=get(past)
              assign(aux, pgb(past.var,mi,m1i,pb1,mmi))
            }
          }
          #tp1=tp1+(pg1 - epgf1(z))^2
          #tp2=tp2+(pg2 - epgf2(z))^2
          #tp3=tp3+(pg3 - epgf3(z))^2
          for (i in 1:GC){
            theoretical=paste0('pg', i)
            theoretical.var=get(theoretical)
            empirical=paste0('epgf', i)
            empirical.var=get(empirical)
            aux=paste0('tp',i)
            assign(aux, get(aux) + (empirical.var(z) - theoretical.var)^2 )
          }
        }
        tpt1 <- tp1
        #tpt2 <- tp1 + tp2
        #tpt3 <- tp1 + tp2 + tp3
        if(GC>1){
          for (i in 2:GC){
            aux=paste0('tpt', i)
            tp_now=paste0('tp', i)
            tp_now.var=get(tp_now)
            past=paste0('tpt', i-1)
            past.var=get(past)
            assign(aux, past.var + tp_now.var)
          }
        }
        #if (tpt1 < mne1) {
        #  mne1 <- tpt1
        #  ls1 <-  c(log(mi / nn1), log(m1i), log(mmi / nn1), log(tpt1))
        #}
        #if (tpt2 < mne2) {
        #  mne2 <- tpt2
        #  ls2 <- c(log(mi / nn1), log(m1i), log(mmi / nn1), log(tpt2))
        #}
        #if (tpt3 < mne3) {
        #  mne3 <- tpt3
        #  ls3 <-c(log(mi / nn1), log(m1i), log(mmi / nn1), log(tpt3))
        #}
        for(i in 1:GC){
          tpt.name=paste0('tpt',i)
          tpt.var=get(tpt.name)
          bound.name=paste0('mne',i)
          bound.var=get(bound.name)
          if(tpt.var < bound.var){
            assign(bound.name, tpt.var)
            assign(paste0('ls',i), c(log(mi / nn1), log(m1i), log(mmi / nn1), log(tpt.var))  )
          }
        }

        #tls <- cbind(tls, c(log(mi / nn1), log(m1i), log(mmi / nn1)),
        #             c(log(tpt1), log(tpt2), log(tpt3))
        #)
      }
    }
  }

  #tls1 <- cbind(tls1, ls1)
  #tls2 <- cbind(tls2,ls2)
  #tls3 <- cbind(tls3, ls3)
  for(i in 1:GC){
    assign(paste0('tls',i), get(paste0('ls',i)))
  }
  #lls <- cbind(lls, tls)

  aux=c(tls1[1:3])
  if(GC>1){
    for(i in 2:GC){
      assign(paste0('aux'), cbind(aux, get(paste0('tls',i))[1:3]))
    }
  }
  #aux<-cbind(tls1[1:3,],tls2[1:3,], tls3[1:3,])
  aux=t(exp(aux))
  naming=c('Growth cycle 1:')
  if(GC>1){
    for(i in 2:GC){
      a=paste0('Growth cycle ',i,':')
      naming=c(naming,a)
    }
  }
  row.names(aux)=naming
  colnames(aux)=c(' WT to Mutant', ' WT to Mutator', ' Mutator to Mutant')

  return(aux)
}
