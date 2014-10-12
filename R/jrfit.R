jrfit <-
function (x, y, block, yhat0 = NULL, scores = wscores, fitint = NULL, 
    var.type = 'sandwich',fitblock=FALSE,tuser=NULL,...) 
{

  call<-match.call()

    if( var.type == 'sandwich' ) {
		v1<-tsand
	}
    if( var.type == 'cs' ) {
		v1<-tcs
	}
	if( var.type == 'ind' ) {
		v1<-tind
	}
    if( var.type == 'user' ) {
		if( !exists('tuser') ) stop('tuser not defined') 
		v1<-tuser
    }

    x <- as.matrix(x)
    x1 <- as.matrix(cbind(rep(1, nrow(x)), x))
    qrx1 <- qr(x1)

    if (is.null(fitint)) {
        if (qrx1$rank == ncol(x1)) {
            x <- x1
			fitint<-TRUE
		} else {
			fitint<-FALSE
		}
    } else {
        if (fitint) {
            x <- x1
		}
    }

	P<-ncol(x) # number of non-nuisance parameters

    Q <- as.matrix(qr.Q(qrx1))
    q1 <- Q[, 1]
    xq <- as.matrix(Q[, 2:qrx1$rank])

	if( fitblock ) {

		z<-model.matrix(~as.factor(block)-1)
		z<-z[,2:ncol(z)]

		# get basis matrix for z perp x1
		QZ<-cbind(Q,z)
		qrxz<-qr(QZ)

		Qxz<-qr.Q(qrxz)

		zq<-Qxz[,(qrx1$rank+1):qrxz$rank]

		xq<-cbind(xq,zq)

#		xz<-cbind(x,z)
#		qrxz<-qr(xz)
#		x<-xz[,1:qrxz$rank]

		x<-cbind(x,zq)

#		x1<-cbind(x1,z)
#    	qrx1 <- qr(x1)

	}


    if (is.null(yhat0)) {
        beta0 <- suppressWarnings(rq(y ~ xq - 1)$coef)
    } else {
        beta0 <- lsfit(xq, yhat0, intercept = FALSE)$coef
    }
    fit <- jaeckel(xq, y, beta0, scores = scores)
    if (fit$convergence != 0) {
    	fit <- jaeckel(xq, y, fit$par, scores = scores)
    	if (fit$convergence != 0) warning("Convergence status not zero in jaeckel")
	}
    betahat <- fit$par
    yhat <- xq %*% betahat
    ehat <- y - yhat
    alphahat <- median(ehat)
    ehat <- ehat - alphahat
    yhat <- yhat + alphahat
    bhat <- lsfit(x, yhat, intercept = FALSE)$coefficients
    tauhat <- gettauF0(ehat, ncol(xq), scores)
    xxpxi <- x %*% chol2inv(chol(crossprod(x)))
    A1 <- crossprod(xxpxi, q1)
    A2 <- crossprod(xxpxi, xq)
    sigma0 <- sigmastar(ehat, block, ncol(xq) + 1)
    taus <- taustar(ehat, ncol(xq) + 1)
    V1 <- v1(ehat, xq, block, scores = scores)
    varhat <- sigma0 * taus * taus * tcrossprod(A1) + tauhat * 
        tauhat * tcrossprod(A2 %*% V1, A2)

#	DF<-ifelse(fitblock,length(y)-ncol(xq)-1,length(unique(block)))
	DF<-switch(var.type,
		cs = length(y)-ncol(xq)-1-1,
		ind = length(y)-ncol(xq)-1,
		sandwich = length(unique(block)),
		user = length(unique(block)) # need to figure out how to specify this.
	)

    res <- list(coefficients = bhat, residuals = ehat, fitted.values = yhat, 
        varhat = varhat, x = x, y = y, block=block, tauhat = tauhat, tauhats=taus, qrx1 = qrx1, 
        disp = fit$value, scores = scores, v1 = v1, fitint=fitint, P=P, var.type=var.type,DF=DF)
	res$call<-call
    class(res) <- "jrfit"
    res
}
