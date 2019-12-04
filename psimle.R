library(bbmle)

psimle.simulate.power.data = function(alpha=0.902,beta=0.797,sigma=0.131,Npoints=40) {
  targetvalue = sample(seq(5,25),size=Npoints,replace=T)
  sd = sigma * targetvalue
  predicted.y = alpha*targetvalue^beta
  response = NA
  for(i in seq(1,length(targetvalue))) {
    response[i] = rnorm(1,mean=predicted.y[i],sd=sd[i])
  }
  return(data.frame(targetvalue,response))
}


psimle.powermodel <- function(targetvalue, response, remove.outliers = FALSE, outliers.sd = 3.0)
{
  #Variable Initialization
  number.of.outliers<-0 
  
  #This function is only used if outliers are being removed. If outliers are not being removed,
  #this function will be skipped. 
  if(remove.outliers==TRUE)
  {
    #Step One: Estimate the Parameters over all data (i.e., data that might have outliers)
    temp.model<-mle2(response~dnorm(mean=alpha*targetvalue^beta,sd=(alpha*targetvalue^beta)*sigma),
                     optimizer="nlminb",
                     start=list(alpha=1, beta=1, sigma=1),
                     data=list(response=response,targetvalue=targetvalue))
    temp.alpha <- coef(temp.model)[1]
    temp.beta <- coef(temp.model)[2]
    temp.sigma <- coef(temp.model)[3]
    
    #Step Two: Create upper and lower boundaries that match each targetvalue's mean and SD
    temp.model.mean <- temp.alpha*targetvalue^temp.beta
    temp.model.sd <- temp.model.mean * temp.sigma
    upperboundary <- (temp.model.mean) + (outliers.sd*(temp.model.sd))
    lowerboundary <- (temp.model.mean) - (outliers.sd*(temp.model.sd))
    
    #Step Three: Construct data array that removes any values above upper or below lower boundary
    outlier.data <- data.frame(targetvalue=targetvalue,response=response,lowerlimit=lowerboundary,upperlimit=upperboundary)
    outlier.data <- outlier.data[response < upperboundary & response > lowerboundary,]
    
    #Step Four: Put data back into original names for second estimation (now with no outlier data)
    number.of.outliers <- length(response)-length(outlier.data$targetvalue)
    targetvalue <- outlier.data$targetvalue
    response <- outlier.data$response
  }
  
  
  #Estimate Parameters
  powermodel<-mle2(response~dnorm(mean=alpha*targetvalue^beta,sd=(alpha*targetvalue^beta)*sigma),
                   optimizer="nlminb",
                   start=list(alpha=1, beta=1, sigma=1),
                   data=list(response=response,targetvalue=targetvalue))
  power.alpha <- coef(powermodel)[1]
  power.beta <- coef(powermodel)[2]
  power.sigma <- coef(powermodel)[3]
  
  power.alpha.se <- stdEr(powermodel)[1]
  power.beta.se <- stdEr(powermodel)[2]
  power.sigma.se <- stdEr(powermodel)[3]
  
  out = list("name"="power",
             "trials"=length(response),
             "model"=powermodel,
             "logLik"=as.numeric(logLik(powermodel)),
             "alpha"=power.alpha,
             "beta"=power.beta,
             "sigma"=power.sigma,
             "alpha.se"=power.alpha.se,
             "beta.se"=power.beta.se,
             "sigma.se"=power.sigma.se,
             "outliers"=remove.outliers,
             "num.outliers"=number.of.outliers)
  
  if(remove.outliers) {
    out$alpha.pre=temp.alpha
    out$beta.pre=temp.beta
    out$sigma.pre=temp.sigma
  }
  
  return(out)
}

psimle.linearmodel <- function(targetvalue, response, remove.outliers = FALSE, outliers.sd = 3.0)
{
  #Variable Initialization
  number.of.outliers<-0 
  
  #This function is only used if outliers are being removed. If outliers are not being removed,
  #this function will be skipped. 
  if(remove.outliers==TRUE)
  {
    #Step One: Estimate the Parameters over all data (i.e., data that might have outliers)
    temp.model<-mle2(response~dnorm(mean=alpha+targetvalue*beta,sd=(alpha+targetvalue*beta)*sigma),
                     optimizer="nlminb",
                     start=list(alpha=1, beta=1, sigma=1),
                     data=list(response=response,targetvalue=targetvalue))
    temp.alpha <- coef(temp.model)[1]
    temp.beta <- coef(temp.model)[2]
    temp.sigma <- coef(temp.model)[3]
    
    #Step Two: Create upper and lower boundaries that match each targetvalue's mean and SD
    temp.model.mean <- temp.alpha+targetvalue*temp.beta
    temp.model.sd <- temp.model.mean * temp.sigma
    upperboundary <- (temp.model.mean) + (outliers.sd*(temp.model.sd))
    lowerboundary <- (temp.model.mean) - (outliers.sd*(temp.model.sd))
    
    #Step Three: Construct data array that removes any values above upper or below lower boundary
    outlier.data <- data.frame(targetvalue=targetvalue,response=response,lowerlimit=lowerboundary,upperlimit=upperboundary)
    outlier.data <- outlier.data[response < upperboundary & response > lowerboundary,]
    
    #Step Four: Put data back into original names for second estimation (now with no outlier data)
    number.of.outliers <- length(response)-length(outlier.data$targetvalue)
    targetvalue <- outlier.data$targetvalue
    response <- outlier.data$response
  }
  
  
  #Estimate Parameters
  linearmodel<-mle2(response~dnorm(mean=alpha+targetvalue*beta,sd=(alpha+targetvalue*beta)*sigma),
                    optimizer="nlminb",
                    start=list(alpha=1, beta=1, sigma=1),
                    data=list(response=response,targetvalue=targetvalue))
  linear.alpha <- coef(linearmodel)[1]
  linear.beta <- coef(linearmodel)[2]
  linear.sigma <- coef(linearmodel)[3]
  
  linear.alpha.se <- stdEr(linearmodel)[1]
  linear.beta.se <- stdEr(linearmodel)[2]
  linear.sigma.se <- stdEr(linearmodel)[3]
  
  out = list("name"="linear",
             "trials"=length(response),
             "model"=linearmodel,
             "logLik"=as.numeric(logLik(linearmodel)),
             "alpha"=linear.alpha,
             "beta"=linear.beta,
             "sigma"=linear.sigma,
             "alpha.se"=linear.alpha.se,
             "beta.se"=linear.beta.se,
             "sigma.se"=linear.sigma.se,
             "outliers"=remove.outliers,
             "num.outliers"=number.of.outliers)
  
  if(remove.outliers) {
    out$alpha.pre=temp.alpha
    out$beta.pre=temp.beta
    out$sigma.pre=temp.sigma
  }
  
  return(out)
}

psimle.logmodel <- function(targetvalue, response, remove.outliers = FALSE, outliers.sd = 3.0)
{
  #Variable Initialization
  number.of.outliers<-0 
  
  #This function is only used if outliers are being removed. If outliers are not being removed,
  #this function will be skipped. 
  if(remove.outliers==TRUE)
  {
    #Step One: Estimate the Parameters over all data (i.e., data that might have outliers)
    temp.model<-mle2(response~dnorm(mean=alpha+log(targetvalue)/beta,sd=(alpha+log(targetvalue)/beta)*sigma),
                     optimizer="nlminb",
                     start=list(alpha=1, beta=1, sigma=1),
                     data=list(response=response,targetvalue=targetvalue))
    temp.alpha <- coef(temp.model)[1]
    temp.beta <- coef(temp.model)[2]
    temp.sigma <- coef(temp.model)[3]
    
    #Step Two: Create upper and lower boundaries that match each targetvalue's mean and SD
    temp.model.mean <-temp.alpha+log(targetvalue)/temp.beta
    temp.model.sd <- temp.model.mean * temp.sigma
    upperboundary <- (temp.model.mean) + (outliers.sd*(temp.model.sd))
    lowerboundary <- (temp.model.mean) - (outliers.sd*(temp.model.sd))
    
    #Step Three: Construct data array that removes any values above upper or below lower boundary
    outlier.data <- data.frame(targetvalue=targetvalue,response=response,lowerlimit=lowerboundary,upperlimit=upperboundary)
    outlier.data <- outlier.data[response < upperboundary & response > lowerboundary,]
    
    #Step Four: Put data back into original names for second estimation (now with no outlier data)
    number.of.outliers <- length(response)-length(outlier.data$targetvalue)
    targetvalue <- outlier.data$targetvalue
    response <- outlier.data$response
  }
  
  
  #Estimate Parameters
  logmodel<-mle2(response~dnorm(mean=alpha+log(targetvalue)/beta,sd=(alpha+log(targetvalue)/beta)*sigma),
                 optimizer="nlminb",
                 start=list(alpha=1, beta=1, sigma=1),
                 data=list(response=response,targetvalue=targetvalue))
  log.alpha <- coef(logmodel)[1]
  log.beta <- coef(logmodel)[2]
  log.sigma <- coef(logmodel)[3]
  
  log.alpha.se <- stdEr(logmodel)[1]
  log.beta.se <- stdEr(logmodel)[2]
  log.sigma.se <- stdEr(logmodel)[3]
  
  out = list("name"="log",
             "trials"=length(response),
             "model"=logmodel,
             "logLik"=as.numeric(logLik(logmodel)),
             "alpha"=log.alpha,
             "beta"=log.beta,
             "sigma"=log.sigma,
             "alpha.se"=log.alpha.se,
             "beta.se"=log.beta.se,
             "sigma.se"=log.sigma.se,
             "outliers"=remove.outliers,
             "num.outliers"=number.of.outliers)
  
  if(remove.outliers) {
    out$alpha.pre=temp.alpha
    out$beta.pre=temp.beta
    out$sigma.pre=temp.sigma
  }
  
  return(out) 
}

psimle.countmodel <- function(targetvalue, response, remove.outliers = FALSE, outliers.sd = 3.0)
{
  #Variable Initialization
  number.of.outliers<-0 
  
  #This function is only used if outliers are being removed. If outliers are not being removed,
  #this function will be skipped. 
  if(remove.outliers==TRUE)
  {
    #Step One: Estimate the Parameters over all data (i.e., data that might have outliers)
    temp.model<-mle2(response~dnorm(mean=alpha+targetvalue*beta,sd=(alpha+targetvalue*beta)^(1/sigma)),
                     optimizer="nlminb",
                     start=list(alpha=1, beta=1, sigma=1),
                     data=list(response=response,targetvalue=targetvalue))
    temp.alpha <- coef(temp.model)[1]
    temp.beta <- coef(temp.model)[2]
    temp.sigma <- coef(temp.model)[3]
    
    #Step Two: Create upper and lower boundaries that match each targetvalue's mean and SD
    temp.model.mean <- temp.alpha+targetvalue*temp.beta
    temp.model.sd <- temp.model.mean * temp.sigma
    upperboundary <- (temp.model.mean) + (outliers.sd*(temp.model.sd))
    lowerboundary <- (temp.model.mean) - (outliers.sd*(temp.model.sd))
    
    #Step Three: Construct data array that removes any values above upper or below lower boundary
    outlier.data <- data.frame(targetvalue=targetvalue,response=response,lowerlimit=lowerboundary,upperlimit=upperboundary)
    outlier.data <- outlier.data[response < upperboundary & response > lowerboundary,]
    
    #Step Four: Put data back into original names for second estimation (now with no outlier data)
    number.of.outliers <- length(response)-length(outlier.data$targetvalue)
    targetvalue <- outlier.data$targetvalue
    response <- outlier.data$response
  }
  
  
  #Estimate Parameters
  countmodel<-mle2(response~dnorm(mean=alpha+targetvalue*beta,sd=(alpha+targetvalue*beta)^(1/sigma)),
                   optimizer="nlminb",
                   start=list(alpha=1, beta=1, sigma=1),
                   data=list(response=response,targetvalue=targetvalue))
  count.alpha <- coef(countmodel)[1]
  count.beta <- coef(countmodel)[2]
  count.sigma <- coef(countmodel)[3]
  
  count.alpha.se <- stdEr(countmodel)[1]
  count.beta.se <- stdEr(countmodel)[2]
  count.sigma.se <- stdEr(countmodel)[3]
  
  out = list("name"="count",
             "trials"=length(response),
             "model"=countmodel,
             "logLik"=as.numeric(logLik(countmodel)),
             "alpha"=count.alpha,
             "beta"=count.beta,
             "sigma"=count.sigma,
             "alpha.se"=count.alpha.se,
             "beta.se"=count.beta.se,
             "sigma.se"=count.sigma.se,
             "outliers"=remove.outliers,
             "num.outliers"=number.of.outliers)
  
  if(remove.outliers) {
    out$alpha.pre=temp.alpha
    out$beta.pre=temp.beta
    out$sigma.pre=temp.sigma
  }
  
  return(out) 
}

psimle.discontmodel <- function(targetvalue, response, discontPoint = 0)
{
  #Estimate Parameters
  dismodel<-mle2(psimle.discontNLL,
                 optimizer="nlminb",
                 start=list(alpha=1, beta.a=1, beta.b=1, sigma.a=1, sigma.b=1),
                 data=list(discontPoint = discontPoint, response=response,targetvalue=targetvalue))
  dis.alpha <- coef(dismodel)[1]
  dis.betaA <- coef(dismodel)[2]
  dis.betaB <- coef(dismodel)[3]
  dis.sigmaA <- coef(dismodel)[4]
  dis.sigmaB <- coef(dismodel)[5]
  
  dis.alpha.se <- stdEr(dismodel)[1]
  dis.betaA.se <- stdEr(dismodel)[2]
  dis.betaB.se <- stdEr(dismodel)[3]
  dis.sigmaA.se <- stdEr(dismodel)[4]
  dis.sigmaB.se <- stdEr(dismodel)[5]
  
  return(list("name"="discont",
              "trials"=length(response),
              "model"=dismodel,
              "logLik"=as.numeric(logLik(dismodel)),
              "alpha"=dis.alpha,
              "betaA"=dis.betaA,
              "betaB"=dis.betaB,
              "sigmaA"=dis.sigmaA,
              "sigmaB"=dis.sigmaB,
              "alpha.se"=dis.alpha.se,
              "betaA.se"=dis.betaA.se,
              "betaB.se"=dis.betaB.se,
              "sigmaA.se"=dis.sigmaA.se,
              "sigmaB.se"=dis.sigmaB.se,
              "discontPoint"=discontPoint))
}


psimle.discontNLL <- function(alpha, beta.a, beta.b, sigma.a, sigma.b, discontPoint, response, targetvalue)
{
  res <- rep(0,length(targetvalue))
  for(i in 1:length(targetvalue)) 
  {
    if(targetvalue[i]<=discontPoint)
    {
      res[i]<- dnorm(response[i], (alpha*targetvalue[i]^beta.a), sigma.a*(alpha*targetvalue[i]^beta.a), log= T)
    }
    else
    {
      res[i] <- dnorm(response[i], (alpha*targetvalue[i]^beta.b), sigma.b*(alpha*targetvalue[i]^beta.b), log= T)
    }
  }
  nll.value <- -sum(res)
}

psimle.AIC<-function(psimle.model1, psimle.model2, critical = 3)
{ 
  if(min(psimle.model1$trials,psimle.model2$trials)<160)
  {
    aic.table<-AICctab(psimle.model1$model,psimle.model2$model,weights=TRUE,sort=FALSE,nobs=min(psimle.model1$trials,psimle.model2$trials))
    m1.dAIC <- aic.table$dAIC[1]
    m2.dAIC <- aic.table$dAIC[2]
  }
  else
  {
    aic.table<-AICtab(psimle.model1$model,psimle.model2$model,weights=TRUE,sort=FALSE)
    m1.dAIC <- aic.table$dAIC[1]
    m2.dAIC <- aic.table$dAIC[2]
  }
  
  dAIC <- max(m1.dAIC,m2.dAIC)
  
  better.model <- "neither"  
  if((m1.dAIC==0)&&(dAIC >= critical))
    better.model <- psimle.model1$name
  if((m2.dAIC==0)&&(dAIC >= critical))
    better.model <- psimle.model2$name
  
  m1.weight <- round(aic.table$weight[1],3)
  m2.weight <- round(aic.table$weight[2],3)
  
  weight.ratio = max(m1.weight/m2.weight, m2.weight/m1.weight)
  
  return(list("betterModel"=better.model,                 
              "dAIC"=dAIC,
              "m1weight"=m1.weight,
              "m2weight"=m2.weight,
              "wratio"=weight.ratio))                
}


psimle.likelihoodratio<-function(psimle.model, testAlpha = coef(psimle.model$model)[1], testBeta = coef(psimle.model$model)[2], testSigma = coef(psimle.model$model)[3])
{  
  model.ll<- -2*c(logLik(psimle.model$model))
  test.ll<- 2*psimle.model$model@minuslogl(testAlpha,testBeta,testSigma)
  degfree <- nargs()-1
  chisq.value <- abs(model.ll - test.ll)
  p.value<- pchisq(chisq.value,degfree,lower.tail=FALSE)
  
  return(list("modelLogLik"=model.ll,
              "testLogLik"=test.ll,
              "chiSq"=chisq.value,
              "degfree"=degfree,
              "pValue"=p.value))
}

psimle.waldtest<-function(estValue, testValue, estSE)
{
  degfree<- 1
  chisq.value<- ((estValue-testValue)^2)/(estSE^2)
  p.value <- pchisq(chisq.value,degfree,lower.tail=FALSE)
  return(list("chiSq"=chisq.value,
              "degfree"=degfree,
              "pValue"=p.value))
}

psimle.generatePlots<-function(targetvalue, response, power=NA, linear=NA, log=NA, counting=NA, discontinuous=NA, outlier=NA, filename=NA) {
  
  if(is.na(outlier)) {
    outlier = 2.5
  }
  
  value = 300 #the unit for the 2:3 ratio of the width to height in the output image
  #png(filename=filename,width=2*value,height=3*value,units="px",res=120) #used for png
  if(!is.na(filename)) {
    pdf(file=filename,width=6,height=9)
  }
  par(mfrow=c(3,2),mgp = c(1.8, 0.5, 0),mar=c(3,3,3,1))
  #mgp=(the distance of the axis labels from the axis, the distance from the axis tick labels to the axis tick marks, the distance away from the graph for the axes)
  #mar=(bottom margin, the left margin, the top margin, the right margin)
  
  range = max(targetvalue)-min(targetvalue)
  padding = 1 #a percentage indicating how much padding to add to the x-axis
  
  x = seq(min(targetvalue)-padding*range,max(targetvalue)+padding*range,length.out=100) #the x values used to plot model lines and fan lines
  #If log is being plotted, automatically only include greater than 0
  if(!is.na(log)[1]) {
    x = x[x > 0]
  }
  
  allYValues = c()
  if(!is.na(power)[1]) {
    y = power$alpha*x^power$beta
    thesd = y*power$sigma
    power.y1 = y+thesd*outlier
    power.y2 = y-thesd*outlier
    allYValues = c(allYValues,power.y1,power.y2)
  }
  if(!is.na(linear)[1]) {
    y = linear$beta*x+linear$alpha
    thesd = y*linear$beta
    linear.y1 = y+thesd*outlier
    linear.y2 = y-thesd*outlier
    print(linear.y2)
    print(y)
    allYValues = c(allYValues,linear.y1,linear.y2)
  }
  if(!is.na(log)[1]) {
    y = log$alpha+log(x)/log$beta
    thesd = y*log$sigma
    log.y1 = y+thesd*outlier
    log.y2 = y-thesd*outlier
    allYValues = c(allYValues,log.y1,log.y2)
  }
  if(!is.na(counting)[1]) {
    y = counting$beta*x+counting$alpha
    thesd = y^(1/counting$sigma)
    counting.y1 = y+thesd*outlier
    counting.y2 = y-thesd*outlier
    allYValues = c(allYValues,counting.y1,counting.y2)
  }
  
  
  xlab = "Target Value"
  ylab = "Response"
  ylower = min(allYValues,na.rm=T)
  yupper = max(allYValues,na.rm=T)
  ylim = c(ylower,yupper) #used if you don't want the fan lines cut off in the figure 
  ylim = c(min(response),max(response))
  
  plot(targetvalue,response,ylim=ylim,main="Raw Data",xlab=xlab,ylab=ylab,lwd=1,pch=21,bg="gray",panel.first=grid(lty=6))
  if(!is.na(power)[1]) {
    plot(targetvalue,response,ylim=ylim,main="Power",bg="gray",lwd=1,pch=21,xlab=xlab,ylab=ylab,panel.first=grid(lty=6))
    lines(x,power$alpha*x^power$beta,col="blue",lwd=2)
    lines(x,power.y1,lty=2)
    lines(x,power.y2,lty=2)
  }
  if(!is.na(linear)[1]) {
    plot(targetvalue,response,ylim=ylim,main="Linear",bg="grey",lwd=1,pch=21,xlab=xlab,ylab=ylab,panel.first=grid(lty=6))
    lines(x,linear$beta*x+linear$alpha,col="blue",lwd=2)
    lines(x,linear.y1,lty=2)
    lines(x,linear.y2,lty=2)
  }
  if(!is.na(log)[1]) {
    plot(targetvalue,response,ylim=ylim,main="Log",bg="grey",lwd=1,pch=21,xlab=xlab,ylab=ylab,panel.first=grid(lty=6))
    lines(x,log$alpha+log(x)/log$beta,col="blue",lwd=2)
    lines(x,log.y1,lty=2)
    lines(x,log.y2,lty=2)
  }
  if(!is.na(counting)[1]) {
    plot(targetvalue,response,ylim=ylim,main="Counting",bg="grey",lwd=1,pch=21,xlab=xlab,ylab=ylab,panel.first=grid(lty=6))
    lines(x,counting$beta*x + counting$alpha,col="blue",lwd=2)
    lines(x,counting.y1,lty=2)
    lines(x,counting.y2,lty=2)
  }
  if(!is.na(discontinuous)[1]) {
    splitAt = discontinuous$discontPoint
    plot(targetvalue,response,ylim=ylim,main="Discontinuous",bg="grey",lwd=1,pch=21,xlab=xlab,ylab=ylab,panel.first=grid(lty=6))
    lines(x[x<=splitAt],discontinuous$alpha*x[x<=splitAt]^discontinuous$betaA,col="cadetblue3",lwd=2)
    lines(x[x>=splitAt],discontinuous$alpha*x[x>=splitAt]^discontinuous$betaB,col="cornflowerblue",lwd=2)
  }
  
  if(!is.na(filename)) {
    dev.off()
  }
}
