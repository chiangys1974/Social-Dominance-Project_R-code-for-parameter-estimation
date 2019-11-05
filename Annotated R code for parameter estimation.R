#### R code for parameter estimation 
#### Acknowledgement: Special thanks to Tsung-Tang Lee for assistance in refining the code
#### Warning: the following code takes about 2 days to run in the following computational environment:
####          Processor: Intel(R) Core(TM) i7-8700 CPU @ 3.2GHZ; Memory: 64GB

## Import supporting libraries
  library(stringr)
  library(data.table)
  library(pipeR)
  library(doParallel)
  library(foreach)

## Import data and set up constants
  N=10 # Number of contestants
  data=read.table("temp/data  v1.txt",header=TRUE) # data is available upon request made to the authors
  data=as.matrix(data)

# Match-up matrix (Table 1)
  Records=c(2,0,0,0,1,0,0,1,0,1, 
          1,2,0,1,1,0,0,0,0,0,
          1,1,2,1,1,0,0,0,0,0,
          1,0,0,2,1,0,0,1,0,1,
          0,0,0,0,2,0,1,1,1,1,
          1,1,1,1,1,2,0,0,0,0,
          1,1,1,1,0,1,2,0,0,0,
          0,1,1,0,0,1,1,2,1,0,
          1,1,1,1,0,1,1,0,2,0,
          0,1,1,0,0,1,1,1,1,2) # we code diagonal (self-self) as "2"
  RM=matrix(Records,byrow=TRUE,nrow=10,ncol=10)

## Customized function for parameter optimization (minimization of the difference between predicted and observed values)
func1 = function(x ){
  # x[1] alpha; x[2] beta; x[3] gamma; x[4] k  
  Score = NULL
  for (i in 1:10) {
    
    # Beta measure
    Si=which(RM[i,]==1) # i's losers
    Normalized_Beta=1/(((N-1)-(rowSums(RM[Si,])-2))^x[2]) # -2 here to substract the diagonal value in the match-up matrix
    Beta_Measure=sum(Normalized_Beta)
    
    # Gamma measure
    Pi=which(RM[i,]==0) # i's winners
    Normalized_gamma=1/((rowSums(RM[Pi,])-2)^x[3])
    Gamma_Measure=sum(Normalized_gamma)
    
    # Integrated measure
    Synthesis=x[1]*Beta_Measure-(1-x[1])*Gamma_Measure
    
    # converted to the Sigmoidd function
    Sy=1/(1+exp(-1*x[4]*Synthesis))
    Score=c(Score,100*Sy)
  }
    Delta=sum((data[person,1:10]-(Score))^2 )
}

## Set the boundaries for parameters alpha, beta, gamma and k
## Note that the first three parameters are bound between 0 and 1
  lower = c(0,0,0,-10)
  upper = c(1,1,1,10)


## Code for parallel computation 
  myCoreNums <- detectCores()
  cl <- makeCluster(myCoreNums-1)
  registerDoParallel(cl)

## Set up data files
  allresult = data.table()
  result.list1 = list()
  result.list2 = list()
  result.list3 = list()

## Estimation of fitted parameter values for each participant
for (person in 1:nrow(data)) {
  
  # Random Initialization
  templist = NULL
  ptm <- proc.time() # record the computation time
  templist = foreach(i=1:10000 , .combine = rbind)  %dopar% { # a specific function employed in the library "foreach"
    set.seed(10000+i) # set up random seeds so that readers can replicate the simulation result
    in.alpha = sample(0:1000000 , 1)/1000000 # initial values for alpha
    in.beta = sample(0:1000000 , 1)/1000000  #                for beta
    in.gamma = sample(0:1000000 , 1)/1000000 #                for gamma
    in.k = sample(-10000000:10000000 , 1)/1000000 #           for k
    
    # Use the optimal function in R to get the best-fit values 
    temp.result = optim(par=c(in.alpha,in.beta,in.gamma,in.k) ,func1 , 
                        lower=lower
                        ,upper=upper ,method = 'L-BFGS-B') 
                  # "func1" (defined above) is the target for optimization 
    
    # Pursuit of accuracy to the 6th decimal points
    # save the estimation results
    return.mat = c('person' = person , setseed = 10000+i ,'ini.a' = in.alpha , 'ini.b' = in.beta  ,'ini.g' = in.gamma , 'ini.k' = in.k , 
      'alpha' = round(temp.result$par[[1]],5) ,'beta' = round(temp.result$par[[2]],5) , 'gamma' = round(temp.result$par[[3]],5) ,
      'k' = round(temp.result$par[[4]],5) , 'delta' = temp.result$value, 'convergence' = temp.result$convergence )

  } # end loop for 10,000 trials
  
  # the following code works to save memory space
    templist = as.data.table(templist) # combine data 
    if(nrow(allresult) < 1 ){ # 1st row special case
      allresult = templist
    }else{
      allresult = rbindlist(list(allresult , templist)  ,use.names = T ,fill = T)
    }
    parLapplyTime <- proc.time() - ptm # record running time for estimation for each participant  
      cat('person:' , person , '\n' )
     
}# end loop for person

stopImplicitCluster() # shut down the library "doParallel"

allresult[,min.del := ifelse(delta == min(delta) , T , F ) , by = 'person'] 
# For each person, among the 10,000 initial values, identify the fittest values (see the min function)
mindel.res = allresult[min.del== T,] 
# Discard the remaining 99999 inferior cases
mindel.res = mindel.res[!duplicated(person)]
# If more than one fittest value, selects the first one
