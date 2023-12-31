# QuantileM_GrowthCycles
#
# Phenotypic dynamic simulations by the quantile function model along n growth cycles
#
#
# You can learn more about package authoring with RStudio at:
#
#

#Bottleneck process
sampling2<-function(population_size,mutants,mutators,mutator_mutants,sample_size){
  s=sample.int(n=4,size=sample_size,replace = T, prob = c(1-(mutants/population_size+mutators/population_size+                                                                mutator_mutants/population_size),
                                                          mutants/population_size,mutators/population_size,
                                                          mutator_mutants/population_size))
  total_wildtype=length(subset(s, s==1))
  total_mutants=length(subset(s, s==2))
  total_mutators=length(subset(s, s==3))
  total_mutator_mutants=length(subset(s, s==4))
  result=matrix(,1,4)
  result[1,1]=total_wildtype
  result[1,2]=total_mutants
  result[1,3]=total_mutators
  result[1,4]=total_mutator_mutants
  colnames(result)=c('Wildtype','Mutants','Mutators','Mutator mutants')
  return(result)
}

#Second layer of the quantile function model

f1<-function(initial_pop,Final_pop, final_time, mutation_rate){
  if(is.na(Final_pop) || Final_pop==0 || is.na(final_time) || final_time==0){
    result<-c(0,0)
  }
  else{
    #growth rate
    #if Final_pop=1, then gr=0 and 1/gr is undefined
    if(Final_pop==1){
      result<-c(0,0)
    }
    else{
      #Suppose the populattion growths as initial_pop*exp (gr* t), t in [0,final_time].
      #Therefore, initial_pop*exp (gr*final_time)= Final_pop. Implying
      #gr=log(Final_pop/initial_pop) /final_time
      gr<- log(Final_pop/initial_pop)/final_time
      #Number of mutations
      mutations<-rpois(1,Final_pop*mutation_rate)
      if(mutations==0){
        result<-c(0,0)
      }
      else{
        #Sample of size 'mutations' from a uniform variable on [0,1]

        uniform<-runif(mutations)
        #Time when mutants appeared
        t<-1/gr* log(uniform*(Final_pop-initial_pop)/initial_pop+1)
        #Time between the final time and the time when mutants appeared
        time_diff<- final_time - t
        mutants_offspring<-rep(-1,mutations)
        for(i in 1:mutations){
          mutants_offspring[i]<-1+rnbinom(1,1,prob=exp(-gr*time_diff[i]))
          if(is.na(rnbinom(1,1,prob=exp(-gr*time_diff[i])))){
            print(time_diff[i])
          }
        }
        totalmutants<-sum(mutants_offspring)
        totalmutants
        #Original mutants and final raised mutants
        result<-c(mutations, totalmutants)
      }
    }
  }
  return(result)
}

#First layer of the quantile function

f2<-function(initial_pop, Final_pop, final_time, mutation_rate, mutation_rate_mutators){
  #growth factor
  r<-log(Final_pop/initial_pop)/final_time
  #Number of mutations
  mutations<-rpois(1,Final_pop*mutation_rate)
  if(mutations==0){
    result<-c(0,0,0,0)
  }
  else{
    #Sample of size 'mutations' from a uniform variable on [0,1]
    uniform<-runif(mutations)
    #Time when mutants appeared
    t<-1/r*  log(uniform*(Final_pop-initial_pop)/initial_pop+1)
    #Time between the final time and the time when mutants appeared
    time_diff<- final_time - t
    #Number of mutants raised by each mutation
    #Create an auxiliary vector with -1's
    mutants<-rep(-1,length(time_diff))
    #Number of mutants as a sample from a r.v. with negative binomial distribution and parameters (1,     exp(-r*time_diff[i]) )
    for(i in 1:mutations){
      #A pure-birth process
      mutants[i]<-1+rnbinom(1,1,prob=exp(-r*time_diff[i]))
    }
    totalmutants<-sum(mutants)
    totalmutants
    #Create an auxiliary vector with -1's for the mutator mutants


    mutator_offspring<-rep(-1,mutations)
    mutators<-rep(-1,mutations)
    for(i in 1:mutations){
      a<-f1(initial_pop=1, Final_pop = mutants[i], time_diff[i], mutation_rate_mutators)
      mutator_offspring[i]<-a[2]
      mutators[i]<-a[1]
    }
    total_mutator_mutants<-sum(mutator_offspring)
    total_mutators<-sum(mutators)
    #Original mutants, mutants at the end, mutators, mutator-mutants
    result<-c(mutations, totalmutants, total_mutators, total_mutator_mutants)
  }
  return(result)
}

#Quantile function model along n growth cycles

QuantileM_GrowthCycles<-function(initial_pop,final_pop,mu_WR,mu_WM,mu_MR,n_cycles,bottle_neck){
  #Variables to storage the output
  for (i in 1:n_cycles){
    assign(paste0('Mutations_to_mutants_GC_', i), -1)
  }
  for (i in 1:n_cycles){
    assign(paste0('Mutants_GC_', i), -1)
  }
  for (i in 1:n_cycles){
    assign(paste0('Mutations_to_mutators_GC_', i), -1)
  }
  for (i in 1:n_cycles){
    assign(paste0('Mutators_GC_', i), -1)
  }
  for (i in 1:n_cycles){
    assign(paste0('Mutations_to_mutant_mutators_GC_', i), -1)
  }
  for (i in 1:n_cycles){
    assign(paste0('Mutant_mutators_GC_', i), -1)
  }
  #Assume the population grows as f(t)= 1*e^{t*log(2)}, t\in [0, t_f], such that f(t_f)=final_pop.
  #Therefore, t_f= log(final_pop/1)/log(2)
  t=log(final_pop/1)/log(2)
  #Resistant cells
  a1<-f2(initial_pop, final_pop, t, mu_WR, 0)
  Mutations_to_mutants_GC_1=a1[1]
  Mutants_GC_1=a1[2]
  #Mutator and resistant mutator cells
  a2<-f2(initial_pop, final_pop, t, mu_WM, mu_MR)
  Mutations_to_mutators_GC_1=a2[1]
  Mutators_GC_1=a2[2]
  Mutations_to_mutant_mutators_GC_1=a2[3]
  Mutant_mutators_GC_1=a2[4]
  cycle=2
  while (n_cycles-cycle >=0) {
    #An unbiased random sample from the immediate previous growth cycle
    s=sampling2(population_size=final_pop,
                mutants=get(paste0('Mutants_GC_', cycle-1)),
                mutators=get(paste0('Mutators_GC_', cycle-1)),
                mutator_mutants=get(paste0('Mutant_mutators_GC_', cycle-1)),
                sample_size=bottle_neck)
    #Initializing counter variables
    old_mutants=0
    old_mutators=0
    old_mutator_mutants=0
    #Assume the population grows as f(t)=bottle_neck*e^{t*log(2)}, t\in [0, t_f], such that f(t_f)=final_pop.
    #Therefore, t_f= log(final_pop/bottle_neck)/log(2)
    t1=log(final_pop/bottle_neck)/log(2)
    #Resistant, mutator, and resistant mutator cells from the sample 'appeared' on the new culture at time t=0.
    time_difference=t1-0
    if(s[,'Mutants']>0){
      for (k in seq(1,s[,'Mutants'])) {
        #Each resistant cell has from t=0 to t1 to give rise to its offspring.
        #Let's note the new growth rate is log(final_pop/bottle_neck)/t1. Thus the offspring size of each
        #resistant cell can be seen as a sample from a r.v. with negative binomial distribution and parameters
        #(1, eˆ{-r(time_difference)}) plus one. Since the support of the distribution is {0,1,2,...}, but we are
        #modelling a pure-birth process
        old_mutants=old_mutants+ 1+rnbinom(1,1,prob=exp(-(log(final_pop/bottle_neck)/t1)*time_difference))
      }
    }

    if(s[,'Mutators']>0){
      for (k in seq(1,s[,'Mutators'])){
        old_mutators=old_mutators+ 1+rnbinom(1,1,prob=exp(-log(final_pop/bottle_neck)/t1*time_difference))
      }
    }

    if(s[,'Mutator mutants']>0){
      for (k in seq(1,s[,'Mutator mutants'])){
        old_mutator_mutants=old_mutator_mutants+ 1+rnbinom(1,1,prob=exp(-log(final_pop/bottle_neck)/t1*time_difference))
      }
    }
    a1<-f2(bottle_neck,final_pop,t1,mu_WR,0)
    #New mutations leading to resistant cells
    dummy_v=paste0('Mutations_to_mutants_GC_', cycle)
    var=get(dummy_v)
    var=a1[1]
    assign(dummy_v,var)
    #New resistant cells
    dummy_v=paste0('Mutants_GC_', cycle)
    var=get(dummy_v)
    #Total number of resistant cells in the culture
    var=a1[2]+old_mutants
    assign(dummy_v,var)
    a2<-f2(bottle_neck,final_pop,t1,mu_WM,mu_MR)
    #New mutations leading to mutator cells
    dummy_v=paste0('Mutations_to_mutators_GC_', cycle)
    var=get(dummy_v)
    var=a2[1]
    assign(dummy_v,var)
    #New mutator cells
    dummy_v=paste0('Mutators_GC_', cycle)
    var=get(dummy_v)
    #Total number of mutator cells in the culture
    var=a2[2]+old_mutators
    assign(dummy_v,var)
    #New mutations leading to resistant mutator cells
    dummy_v=paste0('Mutations_to_mutant_mutators_GC_', cycle)
    var=get(dummy_v)
    var=a2[3]
    #New resistant mutator cells
    assign(dummy_v,var)
    dummy_v=paste0('Mutant_mutators_GC_', cycle)
    var=get(dummy_v)
    #Total number of resistant mutator cells in the culture
    var=a2[4]+old_mutator_mutants
    assign(dummy_v,var)

    #Go to the next growth cycle
    cycle=cycle+1
  }
  #Initialization of the output variable
  output=rep(0,6*n_cycles)
  #Defining the output entry by entry with the information of each growth cycle
  naming=c()
  for (j in 1:n_cycles) {
    h=6*(j-1)
    dummy_v=paste0('Mutations_to_mutants_GC_', j)
    output[h+1]=get(dummy_v)
    naming=c(naming, dummy_v)
    dummy_v=paste0('Mutants_GC_', j)
    output[h+2]=get(dummy_v)
    naming=c(naming, dummy_v)
    dummy_v=paste0('Mutations_to_mutators_GC_', j)
    output[h+3]=get(dummy_v)
    naming=c(naming, dummy_v)
    dummy_v=paste0('Mutators_GC_', j)
    output[h+4]=get(dummy_v)
    naming=c(naming, dummy_v)
    dummy_v=paste0('Mutations_to_mutant_mutators_GC_', j)
    output[h+5]=get(dummy_v)
    naming=c(naming, dummy_v)
    dummy_v=paste0('Mutant_mutators_GC_', j)
    output[h+6]=get(dummy_v)
    naming=c(naming, dummy_v)
  }
  names(output)=naming
  return(output)
}
