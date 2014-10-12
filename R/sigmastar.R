sigmastar <-
function (ehat, block, p) 
{
    ublock <- unique(block)
    m <- length(ublock)
    nvec <- vector(m, mode = "numeric")
    signe <- sign(ehat)
    total <- 0
    for (k in 1:m) {
        signek <- signe[block == ublock[k]]
        nvec[k] <- length(signek)
		if( nvec[k] > 1 ) {
        	for (i in 1:(nvec[k] - 1)) {
            	for (j in (i + 1):nvec[k]) {
                	total <- total + signek[i] * signek[j]
            	}
        	}
		}
    }
    M <- sum(choose(nvec, 2)) - p
    rhos <- total/M
    nstar <- sum(nvec * (nvec - 1))/sum(nvec)
    1 + nstar * rhos
}
