summary.jrfit <-
function (object,int=!object$fitint,df=NULL,...) 
{
    x <- as.matrix(object$x)
    ublock <- unique(object$block)
    m <- length(ublock)
    est <- object$coef[1:object$P]
    ses <- sqrt(diag(object$varhat))[1:object$P]
    tstat <- est/ses
	if( is.null(df) ) df<-object$DF
	pval<-2*pt(abs(tstat),df=df,lower.tail=FALSE)
		
    coef <- cbind(est, ses, tstat,pval)
#	if( !int ) coef<-coef[-1,,drop=FALSE]
#	if( !is.matrix(coef) ) coef<-matrix(coef,ncol=4)
    colnames(coef) <- c("Estimate", "Std. Error", "t-value","p.value")
    ans <- list(coefficients = coef)
	ans$call<-object$call
    class(ans) <- "summary.jrfit"
    ans
}
