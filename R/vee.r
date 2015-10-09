vee<-function(ehat,center,method='dhl',scores=wscores) {

	disp0<-function(e,scores=wscores) {
		r <- rank(e, ties.method = "first")/(length(e) + 1)
		getScores(scores, r) %*% e
	}


	estf<-switch(method,
		'dhl'=signedrank,
		'mm'=median
	)

	m<-tapply(ehat,center,estf)	
	
	sigb2<-switch(method,
		'dhl'=(pi/3)*((disp0(m,scores)/length(unique(center)))^2),
		'mm'=mad(m)^2
	)

	ehatp<-ehat-model.matrix(~as.factor(center)-1)%*%m

	sige2<-switch(method,
		'dhl'=(pi/3)*((disp0(ehatp,scores)/length(center))^2),
		'mm'=mad(ehatp)^2
	)

	list(sigb2=sigb2,sige2=sige2)

}

