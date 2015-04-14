tcs <-
function (ehat, X, block, ahat = NULL, scores = wscores, eps = 1e-04) 
{
    if (is.null(ahat)) {
       r <- rank(ehat, ties.method = "first")/(length(ehat) + 1)
       ahat <- getScores(scores, r)
    }
    nvec <- tapply(ehat, block, length)
    M <- sum(choose(nvec, 2)) - ncol(X)
    myfunc1 <- function(ak) {
        Ak <- tcrossprod(ak)
        Ak[upper.tri(Ak)]
    }
    rhohat <- sum(unlist(tapply(ahat, block, myfunc1)))/M
    K0 <- -1/(max(nvec) - 1)
    if (rhohat < K0) 
        rhohat <- K0 + eps
    if (rhohat > 1) 
        rhohat <- 1 - eps
    myfunc <- function(i, x, rhohat) {
        xi <- as.matrix(x[i, ])
        ni <- nrow(xi)
        sk <- matrix(rhohat, nrow = ni, ncol = ni)
        diag(sk) <- 1
        t(xi) %*% sk %*% xi
    }
    Reduce("+", tapply(1:nrow(X), block, myfunc, X, rhohat))
}
