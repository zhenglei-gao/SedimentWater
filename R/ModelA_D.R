#library(ccSolve)

### Parameters to be estimated in Kenner
# CO2@0: 0.625822
# D_m: 0.48847
# D_p: 0.0409886
# Koc: 983.724
# Ms@0: 0.712317
# Mw@0: 0.644232
# NER@0: 1.38462
# Ps@0: 1.08405
# Pw@0: 95.8947
# Zs: 1.97313
# Zwc: 7.26934
# dKoc: 0.197997
# dkaer: 0.124355
# focsed: 4.99524
# ksmc: 1.19434
# ksmn: 2.81008
# kspm: 0.000651739
# kspn: 0.0016952
# kwmc: 1.49565
# kwpm: 0.00149592
# sigma: 4.93111
# theta: 0.705821



#' Define the toxswa functions in C
#' 
#' Defind the ODE functions
#'
#' @param n number of sediment layers
#' @param language the language used
#'
#' @return a list of functions in language
#' @export
gen.TOXSWA <- function(n=4,language=c("C","F95", "Fortran")){
  language <- match.arg(language)
  if(language=="C"){
    if(n==1){
      odeCode <- " 
  double fdiss= 1/(1+Kd*rhob/theta);
  dPw = -kpwa*Pw-2*Dp/(Zs0)*(Pw/Zw-fdiss*Ps0/Zs0);
  dPs0 = -kpse*Ps0+2*Dp/(Zs0)*(Pw/Zw-fdiss*Ps0/Zs0);
"
    }else{
      odeCode <- " 
  double fdiss= 1/(1+Kd*rhob/theta);
  dPw = -kpwa*Pw-2*Dp/(Zs0)*(Pw/Zw-fdiss*Ps0/Zs0);
  dPs0 = -kpse*Ps0+2*Dp/(Zs0)*(Pw/Zw-fdiss*Ps0/Zs0)-2*Dp*fdiss/(Zs0+Zs1)*(Ps0/Zs0-Ps1/Zs1);
"
      if(n>2){
        for(i in 1:(n-2)){
          addCode <- paste0("dPs",i, "= -delta*kpse*Ps",i,"+2*Dp*fdiss/(Zs",i-1,"+Zs",i,")*(Ps",i-1,
                            "/Zs",i-1,"-Ps",i,"/Zs",i,")-2*Dp*fdiss/(Zs",i,"+Zs",i+1,")*(Ps",i,"/Zs",i,"-Ps",i+1,
                            "/Zs",i+1,");")
          odeCode <- paste0(odeCode,addCode,"\n")
        }
        i <- n-1
        addn <- paste0("dPs",i, "= -delta*kpse*Ps",i,"+2*Dp*fdiss/(Zs",i-1,"+Zs",i,")*(Ps",i-1,
                       "/Zs",i-1,"-Ps",i,"/Zs",i,")")
        odeCode <- paste0(odeCode,addn,";\n")
      }
    }
    return(odeCode)
  }
}
# n <- 10
# yini <- c(100,rep(0,n))
# SWModel <- gen.TOXSWA(n=n)
# names(yini) <- c("Pw",paste0("Ps",0:(n-1)))
# Zs <- 2
# Zsed <- rep(Zs/n,n)
# names(Zsed) <- paste0("Zs",0:(n-1))
# parms <- c(kpwa = log(2)/50, kpse = log(2)/50, Kd = 50,
#            rhob=2.05,theta=0.5,Zw=5,Dp=0.05,delta=1,Zsed)
# parms2 <- c(kpwa = log(2)/50, kpse = log(2)/1, Kd = 50,
#            rhob=2.05,theta=0.1,Zw=5,Dp=0.5,delta=1,Zsed)
# cToxswaA <- compile.ode(SWModel, language = "C", 
#                         parms = parms,y=yini) 
# code(cToxswaA)
# checkTime <- FALSE
# if(checkTime) system.time(replicate(1000,{outcA <- ode (func = cToxswaA, y = yini, parms = parms, 
#                           times = 0:100)}))
# plotSW <- function(dta){
#   dta <- data.frame(dta)
#   dta <- data.frame(dta[,1:2],Ps=apply(dta[,3:ncol(dta)],1,sum))
#   dta <- melt(dta,id.vars=1)
#   p <- ggplot(dta,aes(x=time,y=value,color=variable))+geom_line()
#   p
# }
# outcA <- ode (func = cToxswaA, y = yini, parms = parms, 
#               times = 0:100)
# plotSW(outcA)
# yini2 <- yini
# yini2[1] <- 1
# outcA2 <- ode (func = cToxswaA, y = yini2, parms = parms, 
#               times = 0:100)
# plotSW(outcA2)
# ---------------------
# The F95 version 
# ---------------------                        

# ---------------------
# + parameters, C code
# ---------------------

checkRelations <- function(AllParVec){
  Kd <- Koc*focsed/100
  rhob <- 2.5*(1-theta)
  fdiss <- 1/(1+Kd*rhob/theta)
  
}
mkin_wide_to_long <- function (wide_data, time = "t") 
{
  colnames <- names(wide_data)
  if (!(time %in% colnames)) 
    stop("The data in wide format have to contain a variable named ", 
         time, ".")
  vars <- subset(colnames, colnames != time)
  n <- length(colnames) - 1
  long_data <- data.frame(name = rep(vars, each = length(wide_data[[time]])), 
                          time = as.numeric(rep(wide_data[[time]], n)), value = as.numeric(unlist(wide_data[vars])), 
                          row.names = NULL)
  return(long_data)
}

#' Forward calculation
#' 
#' Forward calculation
#'
#' @param cModel A lists including the differential equations, the parameters, the states, the initials. 
#' @param optParms optimazation ODE parameters
#' @param fixParms fixed parameters 
#' @param optInis optimazation initial parameters
#' @param fixInis fixed parameters 
#' @param obsTimes observation times
#' @param transParms transformation functions
#' 
#' @return a data frame
#' 
#' @export
SW_Model_calc <- function (cModel, optParms,fixParms,
                           optInis,fixInis,obsTimes=1:100,
                           transParms = NULL, ...) 
{
  ## model is the compiled code?
  ## data observed.
  
  if(!is.null(transParms)){
    
  }else{
    ## replace OptParms to parms and yini
    parms <- c(optParms,fixParms)
    yini <- c(optInis,fixInis)
    outcA <- ode (func = cModel, y = yini, parms = parms, 
                  times = obsTimes)
  }
  dta <- data.frame(outcA)
  dta <- data.frame(dta[,1:2],Ps=apply(dta[,3:ncol(dta)],1,sum))
  dta <- melt(dta,id.vars=1)
  names(dta) <- c("Time","Media","pred")
  dta <- dta %>%
    mutate(Media = plyr::mapvalues(Media, c("Pw","Ps"), c("Water", "Sediment")))
  levels(dta$Media) <- c("Water","Sediment")
  return(dta)
}
#' Forward calculation
#' 
#' Forward calculation
#'
#' @param data data frame
#' @param cModel A lists including the differential equations, the parameters, the states, the initials. 
#' @param optParms optimazation ODE parameters
#' @param fixParms fixed parameters 
#' @param optInis optimazation initial parameters
#' @param fixInis fixed parameters 
#' @param transParms transformation functions
#' @param ... observation times
#' 
#' @return a data frame
#' 
#' @export
costSW <- function(data,cModel, optParms,fixParms,
                   optInis,fixInis,
                   transParms = NULL, ...){
  obsTimes <- unique(data[,1])
  #ldata <- melt(data,id.vars=1)
  pred <- SW_Model_calc(cModel, optParms,fixParms,optInis,fixInis,obsTimes,
                                    transParms = NULL, ...)
  ##browser()
  ldata <- left_join(data,pred,by=c("Time", "Media"))
  ldata$resid <- ldata$pred- ldata$Concentration
  return(ldata)
}
#' Forward calculation for residual only
#' 
#' Forward calculation
#'
#' @param par parameters
#' @param data data frame
#' @param cModel model definition
#' @param fixParms fixed parameters
#' @param fixInis fixed initial values
#' @param transParms tranformation functions
#' @param sepn seperater index
#' @param ... other unused yet inputs
#' 
#' @return a vector of residuals
#' @export
calcResid <- function(par,data,cModel,fixParms,fixInis,transParms,sepn,...){
#   if(names(data)[2]=="Media"){
#     ## in case there are multiple time points. No problem!
#     wdata <- cast(data,Time~Media,value=Concentration)
#   }else{
#     data <- melt(data,id.vars = 1)
#   }
  optParms <- par[1:sepn]
  npar <- length(par)
  if(npar > sepn) optInis <- par[(sepn+1):npar] else optInis <- NULL
  obsTimes <- unique(data[,1])
  #ldata <- melt(data,id.vars=1)
  pred <- SW_Model_calc(cModel, optParms,fixParms,optInis,fixInis,obsTimes,
                        transParms = NULL, ...)
  ldata <- left_join(data,pred,by=c("Time", "Media"))
  resid <- ldata$pred- ldata$Concentration
  return(resid)
}
#' Get Y
#'
#' Get Y
#' 
#' @param logpar logtransformed parameters
#' @param data data frame
#' @param cModel model definition
#' @param fixParms fixed parameters
#' @param fixInis fixed initial values
#' @param transParms tranformation functions
#' @param sepn seperater index
#' @param ... other unused yet inputs
#'
#' @export
getY <- function(logpar,data,cModel,fixParms,fixInis,transParms,sepn,...){
  #   if(names(data)[2]=="Media"){
  #     ## in case there are multiple time points. No problem!
  #     wdata <- cast(data,Time~Media,value=Concentration)
  #   }else{
  #     data <- melt(data,id.vars = 1)
  #   }
  names(logpar) <- c("kpwa", "kpse", "Kd")
  par <- exp(logpar)
  optParms <- exp(logpar[1:sepn])
  npar <- length(logpar)
  if(npar > sepn) optInis <- par[(sepn+1):npar] else optInis <- NULL
  obsTimes <- unique(data[,1])
  #ldata <- melt(data,id.vars=1)
  pred <- SW_Model_calc(cModel, optParms,fixParms,optInis,fixInis,obsTimes,
                        transParms = NULL, ...)
  ldata <- left_join(data,pred,by=c("Time", "Media"))

  return(ldata$pred)
}
#' Plot a data frame
#'
#' Plot a SW dta 
#' 
#' @param dta a sw data frame
#' 
#' @return a plot
#' @export
plotSW <- function(dta,origDta=NULL){
  dta <- data.frame(dta)
  dta <- data.frame(dta[,1:2],Ps=apply(dta[,3:ncol(dta)],1,sum))
  dta <- melt(dta,id.vars=1)
  names(dta) <- c("Time", "Media","Concentration")
  dta$Media <- mapvalues(dta$Media,from=c("Pw","Ps"),to=c("Water","Sediment"))
  
  p <- ggplot(dta,aes(x=Time,y=Concentration,color=Media))+geom_line()
  if(!is.null(origDta)) p <- p+ geom_point(data=origDta,aes(x=Time,y=Concentration,color=Media))
  p
}


#' Simulate a data frame
#'
#' Simu a SW dta 
#' 
#' @param dta a sw data frame
#' 
#' @return a simulated data frame with 
#' @export
simSW <- function(func,yini,parms,times,sigmaS,plotting=FALSE){
  dta <- ode (func = cToxswaA, y = yini, parms = parms, 
              times = times)
  nt <- length(times)
  if(length(sigmaS)==1) sigmaS <- rep(sigmaS,2)
  dta <- data.frame(dta[,1:2],Ps=apply(dta[,3:ncol(dta)],1,sum))
  dta$Pw <- keeppos(rnorm(nt,dta$Pw,sigmaS[1]))
  dta$Ps <- keeppos(rnorm(nt,dta$Ps,sigmaS[2]))
  dta <- melt(dta,id.vars=1)

  names(dta) <- c("Time", "Media","Concentration")
  dta$Media <- mapvalues(dta$Media,from=c("Pw","Ps"),to=c("Water","Sediment"))
  if(plotting) {
    p <- ggplot(dta,aes(x=Time,y=Concentration,color=Media))+geom_point()+geom_smooth()
    print(p)
  }
  return(dta)
}

keeppos <- function(x){
  ind <- (x<0)
  x[ind] <- -x[ind]
  return(x)
}