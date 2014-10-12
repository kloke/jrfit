tE <-
function (ehat, X, block, ahat = NULL, scores = wscores, fs.correct = TRUE) 
{
    if (is.null(ahat)) 
        ahat <- getScores(scores, ehat)
    A <- matrix(ahat, nrow = m, byrow = TRUE)
    sphi <- cor(A)
    myfunc <- function(i, x, sk) {
        xi <- x[i, ]
        t(xi) %*% sk %*% xi
    }
    vhat <- Reduce("+", tapply(1:nrow(X), block, myfunc, X, sphi))
    if (fs.correct == TRUE) {
        m <- length(unique(block))
        vhat <- m/(m - ncol(X)) * vhat
    }
    vhat
}
