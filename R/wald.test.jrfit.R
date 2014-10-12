wald.test.jrfit<-function(fit,K,asymptotic=FALSE) {

# assumes K is full row rank

	est<-K%*%fit$coef
	var.est<-K%*%fit$varhat%*%t(K)

	statistic<-t(est)%*%solve(var.est)%*%est

	q<-nrow(K)

	if( asymptotic ) {
		p.value<-pchisq(statistic,q,lower.tail=FALSE)
		df<-q
	} else {
		statistic<-statistic/q
		p.value<-pf(statistic,q,fit$DF,lower.tail=FALSE)
		df<-c(q,fit$DF)
	}

	list(statistic=statistic,p.value=p.value,asymptotic=asymptotic,df=df)


}

