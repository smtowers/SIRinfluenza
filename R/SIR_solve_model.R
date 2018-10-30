#' A function to calculate the deterministic estimates of the number
#' infected versus time and the average number of descendant infections (ANDI)
#' for a deterministic Susceptible, Infected, Recovered (SIR) model
#'
#' @param N the population size
#' @param R0 the reproduction number
#' @param gamma the recovery rate, in units one over the units of time
#' @param fsusc the initial fraction susceptible
#' @param time_end the end time of the simulation
#' @param delta_t time step to solve model
#' @return list model compartment values vs time and ANDI
#'
#' @export
SIR_solve_model = function(N,
                           R0,
                           gamma,
                           fsusc,
                           time_end,
                           delta_t=0.01){

  I_0 = 1      
  S_0 = fsusc*N-I_0
  R_0 = (1-fsusc)*N
 
  # note that time is in units of 1/gamma
  tbeg  = 0           
  beta = R0*gamma     
  tend = time_end

  vparameters = c(gamma=gamma,beta=beta)
  inits       = c(S=S_0,I=I_0,R=R_0)

  iend = 1
  while (iend>0.001){
    vt = seq(0,tend,delta_t)
    sirmodel = as.data.frame(deSolve::lsoda(inits, vt, SIRfunc, vparameters))
    sirmodel = subset(sirmodel,!is.na(S+I+R))
    iend = sirmodel$I[nrow(sirmodel)]
    tend = tend*2
  }
  minS = min(sirmodel$S,na.rm=T)
  final = (max(sirmodel$S,na.rm=T)-min(sirmodel$S,na.rm=T))

  incidence = beta*sirmodel$S*sirmodel$I/N
  ################################################################################## 
  ################################################################################## 
  # find the time at which the incidence was at a maximum
  ################################################################################## 
  iind = which.max(incidence)
  time_of_max = sirmodel$time[iind]

  return(list(results    = sirmodel,
              final_size = final,
              incidence  = incidence,
              time_of_max = time_of_max
             ))
} # end solve_model function definition



