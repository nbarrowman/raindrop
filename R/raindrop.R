#'
#' plotraindrops: plot a list of raindrops
#'
#' @description
#' For plotting raindrops
#'
#' @author Nick Barrowman <nbarrowman@cheo.on.ca>
#'
#' @param logliks    Log likelihoods
#' @param theta      Corresponding parameter values.
#' @param cutoff     decrease from maximum in the log likelihood
#' @param type       type of curve
#' @param tfactor    not sure
#' @param col        fill color
#' @param lwd        line width
#' @param xlim       x-limits for plot
#' @param K          number of log likelihoods (not necessary to specify)
#' @param ylim       y-limits for plot
#' @param xaxis      draw the x-axis?
#'
#' @examples
#' studies <- data.frame(trial=c(1,1),x=c(1,15),n=c(11,25),rx=c(1,0))
#' theta <- seq(-8,1,length=200)
#' # create quadratic log likelihood
#' study1fit <- glm(cbind(x,n-x)~rx,family=binomial,data=studies[studies$trial==1,])
#' table <- summary(study1fit)$coef
#' mlelogOR <- table["rx","Estimate"] 
#' ASElogOR <- table["rx","Std. Error"]
#' quadloglik <- -0.5*(theta-mlelogOR)^2/ASElogOR^2
#' # create conditional log likelihood
#' condloglik <- log(clik(theta,1,10,15,10))
#' maxcondloglik <- max(condloglik)
#' plotraindrops(list("Wald"=quadloglik,"Conditional"=condloglik),theta)
#'
#' @export

plotraindrops <- function(logliks,theta,cutoff=-2,type="loglik",tfactor=1.0,
  col=2,lwd=1.0,xlim=range(theta),K=length(logliks),ylim=c(1,3*K+2),xaxis=TRUE) {

  if (length(cutoff)==1) { cutoff <- rep(cutoff,K) }
  if (length(col)==1) { col <- rep(col,K) }
  if (length(tfactor)==1) { tfactor <- rep(tfactor,K) }
  labels <-names(logliks)
  plot(0,0,type="n",xlim=xlim,ylim=ylim,axes=F,xlab="",ylab="")
  if (xaxis) axis(1)
  for (i in 1:K) {
    y <- 3*(K-i+1)
    plotraindrop(logliks[[i]],theta,y=y,cutoff=cutoff[i],type=type,tfactor=tfactor[i],col=col[i],lwd=lwd)
    mtext(labels[i],at=y,side=2,adj=1,las=1)
  }
}



#'
#' plotraindrop: plot a raindrop
#'
#' @description
#' For plotting raindrops
#'
#' @author Nick Barrowman <nbarrowman@cheo.on.ca>
#'
#' @param loglik     Log likelihod
#' @param theta      Corresponding parameter values.
#' @param y          y-value
#' @param cutoff     decrease from maximum in the log likelihood
#' @param type       type of curve
#' @param tfactor    not sure
#'
#' @export
plotraindrop <- function(loglik,theta,y,cutoff=-2,type="loglik",tfactor=1.0,
  col=2,lwd=1.0,rot=FALSE) {
  l <- loglik - max(loglik)
  low <- l < cutoff
  yval <- rep(y,length(theta))
  yval[!low] <- NA

  pvalue <- (1-pchisq(cutoff*l,1))/2  # one-sided p-value

  L <- exp(l)

  par(err=-1)
  if (type=="loglik") {
    select <- l >= cutoff
    THETA <- theta[select]
    ll <- l[select]
    thickness <- tfactor*(1-ll/cutoff)
    if (rot) {
      polygon(c(y-thickness,rev(y+thickness)),c(THETA,rev(THETA)),col=col,border=F)
      lines(y+thickness,THETA,lwd=lwd)
      lines(y-thickness,THETA,lwd=lwd)
    } else {
      polygon(c(THETA,rev(THETA)),c(y-thickness,rev(y+thickness)),col=col,border=F)
      lines(THETA,y+thickness,lwd=lwd)
      lines(THETA,y-thickness,lwd=lwd)
    }
  } else
  if (type=="pvalue") {
    thickness <- tfactor*pvalue
    thickness[l< cutoff] <- 0
    polygon(c(theta,rev(theta)),c(y-thickness,rev(y+thickness)),col=col)
    lines(theta,y+thickness)
    lines(theta,y-thickness)
    lines(theta,yval,col=0)
    lines(theta,yval,lty=8)
  } else
  if (type=="lik") {
    thickness <- tfactor*L/2
    thickness[l< cutoff] <- 0
    polygon(c(theta,rev(theta)),c(y-thickness,rev(y+thickness)),col=col)
    lines(theta,y+thickness)
    lines(theta,y-thickness)
    lines(theta,yval,col=0)
    lines(theta,yval,lty=8)
  } else {
    stop(paste("Unknown type:",type))
  }
  par(err=0)
}



#'
#' clik: conditional likelihood function for a 2x2 table
#'
#' @description
#' conditional likelihood function for a 2x2 table
#'
#' @author Nick Barrowman <nbarrowman@cheo.on.ca>
#'
#' @param theta      log odds ratio
#' @param aa         number in top left
#' @param bb         number in top right
#' @param cc         number in bottom left
#' @param dd         number in bottom right
#'
#' @export
clik <- function(theta,aa,bb,cc,dd) {
  lo <- max(0,aa-dd)
  hi <- min(aa+bb,aa+cc)
  num <- nCm(aa+bb,aa)*nCm(cc+dd,cc)*exp(theta*aa)
  den <- 0 
  for (z in lo:hi) {
    den <- den+nCm(aa+bb,z)*nCm(cc+dd,aa+cc-z)*exp(theta*z)
  }
  return(num/den)
}



"nCm"<-
function(n, m, tol = 1e-08)
{
#  DATE WRITTEN:  7 June 1995               LAST REVISED:  10 July 1995
#  AUTHOR:  Scott Chasalow
#
#  DESCRIPTION: 
#        Compute the binomial coefficient ("n choose m"),  where n is any 
#        real number and m is any integer.  Arguments n and m may be vectors;
#        they will be replicated as necessary to have the same length.
#
#        Argument tol controls rounding of results to integers.  If the
#        difference between a value and its nearest integer is less than tol,  
#        the value returned will be rounded to its nearest integer.  To turn
#        off rounding, use tol = 0.  Values of tol greater than the default
#        should be used only with great caution, unless you are certain only
#        integer values should be returned.
#
#  REFERENCE: 
#        Feller (1968) An Introduction to Probability Theory and Its 
#        Applications, Volume I, 3rd Edition, pp 50, 63.
#
	len <- max(length(n), length(m))
	out <- numeric(len)
	n <- rep(n, length = len)
	m <- rep(m, length = len)
	mint <- (trunc(m) == m)
	out[!mint] <- NA
	out[m == 0] <- 1	# out[mint & (m < 0 | (m > 0 & n == 0))] <-  0
	whichm <- (mint & m > 0)
	whichn <- (n < 0)
	which <- (whichm & whichn)
	if(any(which)) {
		nnow <- n[which]
		mnow <- m[which]
		out[which] <- ((-1)^mnow) * Recall(mnow - nnow - 1, mnow)
	}
	whichn <- (n > 0)
	nint <- (trunc(n) == n)
	which <- (whichm & whichn & !nint & n < m)
	if(any(which)) {
		nnow <- n[which]
		mnow <- m[which]
		foo <- function(j, nn, mm)
		{
			n <- nn[j]
			m <- mm[j]
			iseq <- seq(n - m + 1, n)
			negs <- sum(iseq < 0)
			((-1)^negs) * exp(sum(log(abs(iseq))) - lgamma(m + 1))
		}
		out[which] <- unlist(lapply(seq(along = nnow), foo, nn = nnow, 
			mm = mnow))
	}
	which <- (whichm & whichn & n >= m)
	nnow <- n[which]
	mnow <- m[which]
	out[which] <- exp(lgamma(nnow + 1) - lgamma(mnow + 1) - lgamma(nnow - 
		mnow + 1))
	nna <- !is.na(out)
	outnow <- out[nna]
	rout <- round(outnow)
	smalldif <- abs(rout - outnow) < tol
	outnow[smalldif] <- rout[smalldif]
	out[nna] <- outnow
	out
}
