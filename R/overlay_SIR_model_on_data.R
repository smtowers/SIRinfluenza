#' A function to calculate the deterministic estimates of the number
#' infected versus time and the average number of descendant infections (ANDI)
#' for a deterministic Susceptible, Infected, Recovered (SIR) model
#' and overlay the results on Midwest 2008 influenza outbreak data
#'
#' @param R0 the reproduction number
#' @return NULL
#'
#' @export
overlay_SIR_model_on_data = function(R0){

  cat("R0=",R0,"\n")
  ########################################################################
  # N is approximately the population of IL IN MI MN OH WI (CDC region 5)
  ########################################################################
  N = 52000000
  
  ########################################################################
  # the recovery rate of influenza is approximate 1/3 days^{-1}
  ########################################################################
  gamma = 1/3

  ########################################################################
  # solve the SIR model with these parameters
  ########################################################################
  fsusc = 1.0
  delta_t = 0.1
  time_end = 365

  a = SIR_solve_model(N,R0,gamma,fsusc,time_end,delta_t)

  ########################################################################
  # now load the data and plot it
  ########################################################################
  #load("SIRinfluenza/data/midwest_flu_2008.rda")

  ########################################################################
  # find the time at which the weekly incidence was at a maximum
  # and multiply by 7 to get days
  ########################################################################
  iind = which.max(midwest_flu_2008$num_influenza_cases)
  peak_time = (midwest_flu_2008$week[iind]-0.5)*7

  ########################################################################
  # obtain the model incidence prediction and normalise the results to
  # the data
  # also calculate the time shift needed such that the peak of the 
  # model will align with the peak of the data
  ########################################################################
  delt = peak_time-a$time_of_max
  xtime = a$results$time+delt
  ypred = a$incidence
  i = which(xtime>=0&xtime<=max(midwest_flu_2008$week*7))
  xtime = xtime[i]
  ypred = ypred[i]
  ypred = ypred*sum(midwest_flu_2008$num_influenza_cases)/sum(ypred)/delta_t

  ########################################################################
  # aggregate the model prediction within week, and then compare the 
  # results to the data and calculate goodness-of-fit statistics
  ########################################################################
  xtime_week = as.integer(xtime/7)
  g = aggregate(ypred,by=list(xtime_week),FUN="sum")
  y = g[[2]] 
  y = y[1:nrow(midwest_flu_2008)]
  Pearson = sum((y-midwest_flu_2008$num_influenza_cases)^2/y)
  neg_Pois_loglike = sum(y-midwest_flu_2008$num_influenza_cases*log(y))
  
  ########################################################################
  # plot the data with the model overlaid
  ########################################################################
  ymax = 1.4*max(midwest_flu_2008$num_influenza_cases)

  plot((midwest_flu_2008$week-0.5)*7,midwest_flu_2008$num_influenza_cases,xlab="Days from Jan 1, 2008",ylab="Weekly \043 identified influenza cases",pch=20,cex=4,ylim=c(0,ymax),main="Midwest influenza cases, early 2008",col="black",font.lab=2)
  lines(xtime,ypred*7,col="red",lwd=8)

  ########################################################################
  # overlay text with the Pearson chi-squared and
  # Poisson negative log likelihood statistics
  ########################################################################
  par(new=T)
  plot(c(0,1),c(0,1),col="white",axes=F,xlab="",ylab="",xaxt="n",yaxt="n")
  text(0.50,0.80,paste("Pearson chi^2=",round(Pearson,0),sep=""),adj=c(0,1),cex=0.8,font=2)
  text(0.50,0.75,paste("Poisson neg log likelihood=",round(neg_Pois_loglike,0),sep=""),adj=c(0,1),cex=0.8,font=2)
  legend("topleft",legend=c("Data",paste("SIR model w/ R0=",R0,sep="")),col=c("black","red"),lwd=5,bty="n",cex=1.0)
  
  return()
} # end overlay_SIR_model_on_data function definition



