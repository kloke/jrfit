tsand <-
function (ehat, x, block, a = NULL, scores = wscores, fs.correct = TRUE) 
{
    x <- as.matrix(x)
    if (is.null(a)) 
        a <- getScores(scores, ehat)
    Xa <- cbind(x, a)
    myfunc <- function(i, Xa) {
        Xa <- Xa[i, ,drop=FALSE]
        p <- ncol(Xa) - 1
        tcrossprod(crossprod(Xa[, 1:p,drop=FALSE], Xa[, p + 1,drop=FALSE]))
    }
    vhat <- Reduce("+", tapply(1:nrow(Xa), block, myfunc, Xa))
    if (fs.correct == TRUE) {
        m <- length(unique(block))
        if (m > ncol(x)) {
            vhat <- m/(m - ncol(x)) * vhat
        }
    }
    vhat
}
