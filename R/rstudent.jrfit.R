rstudent.jrfit <-
function (model,...) 
{
	fit<-model
    N <- nrow(fit$x)
    p <- ncol(fit$x)
    block <- fit$block
    ublock <- unique(block)
    m <- length(ublock)
#    nvec <- vector(m, mode = "numeric")
#    for (k in 1:m) {
#        nvec[k] <- sum(block == ublock[k])
#    }
	nvec<-tapply(block,block,length)
    K <- sum(nvec * (nvec - 1)) - p
    ehat <- fit$residual
    r <- rank(e, ties.method = "first")/(length(e) + 1)
    a <- getScores(fit$scores, r)
    d11 <- crossprod(ehat, a)/(N - p)
    ds11 <- crossprod(ehat, sign(ehat))/(N - p)
    g11 <- crossprod(sign(ehat), a)/(N - p)
    d12 <- ds12 <- g12 <- 0
    for (k in 1:m) {
        ehatk <- ehat[block == ublock[k]]
        ak <- a[block == ublock[k]]
        D12 <- outer(ehatk, ak)
        diag(D12) <- 0
        d12 <- d12 + sum(D12)
        Ds12 <- outer(ehatk, sign(ehatk))
        diag(Ds12) <- 0
        ds12 <- ds12 + sum(Ds12)
        G12 <- outer(sign(ehatk), ak)
        diag(G12) <- 0
        g12 <- g12 + sum(G12)
    }
    d12 <- d12/K
    ds12 <- ds12/K
    g12 <- g12/K
    bhat <- vector(m, mode = "numeric")
    epsilonhat <- vector(m, mode = "numeric")
    for (k in 1:m) {
        ehatk <- ehat[block == ublock[k]]
        bhat[k] <- median(ehatk)
        epsilonhat[block == ublock[k]] <- ehatk - bhat[k]
    }
    sigmab <- mad(bhat)^2
    sigmaepsilon <- mad(epsilonhat)^2
    sigmahat <- sigmab + sigmaepsilon
	rhohatR<-sigmab/sigmahat

##############################
### Get rhohat stuff ###
### added 10/15/12 ###
        totals<-total<-0
        for(k in 1:m) {
                ak<-a[block==ublock[k]]
                ek<-ehat[block==ublock[k]]
                nvec[k]<-length(ak)
				if( nvec[k] > 1 ) {
                	for( i in 1:(nvec[k]-1) ) {
                        	for( j in (i+1):nvec[k] ) {
                                	total<-total+ak[i]*ak[j]
                                	totals<-total+sign(ek[i])*sign(ek[j])
                        	}
                	}
				}
        }

        M<-sum(choose(nvec,2))-p
        rho<-total/M
        rhos<-totals/M
##############################


    bo <- order(block)
    ehat <- fit$resid[bo]
    X <- fit$x[bo, ]
    Q <- qr.Q(qr(X))
    H <- Q %*% t(Q)
    JN <- matrix(1, nrow = N, ncol = N)
    C1 <- sigmahat
    C2 <- fit$tauhats^2/N^2 * sum(nvec * (1 + (nvec - 1) * rhos))
    C3 <- fit$tauhat^2 * diag(H %*% Bmat(1, rho, block[bo]) %*% H)
    C4 <- vector(N, mode = "numeric")
    ind <- 1
    for (k in 1:m) {
        C4[ind:(ind + nvec[k] - 1)] <- (nvec[k] - 1) * ds12
        ind <- ind + nvec[k]
    }
    C4 <- -1 * fit$tauhats/N * (C4 + ds11)
    C5 <- -fit$tauhat * diag(Bmat(d11, d12, block[bo]) %*% H)
    C6 <- C4
    C7 <- fit$tauhat * fit$tauhats/N * diag(JN %*% Bmat(g11, g12, block[bo]) %*% H)
    C8 <- -1 * fit$tauhat * diag(H %*% Bmat(d11, d12, block[bo]))
    C9 <- fit$tauhats * fit$tauhat/N * diag(H %*% Bmat(d11, d12, block[bo]) %*% JN)
    C <- C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9

### fix for est ses which are neg ###
### correspondence with McKean indicates we could use ###
# diag( (I-H) sigmahat^2 B(rhohat) (I-H) )

	if( sum( C <= 0 ) > 0 ) {
		H2<-diag(1,nrow(H)) - H		
		CH2<-diag(H2%*%Bmat(1,rhohatR)%*%H2)
		C[C<=0]<-CH2[C<=0]
	}

    rs <- ehat/sqrt(C)
    rs[bo] <- rs
    rs


}
