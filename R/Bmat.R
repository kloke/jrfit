Bmat <-
function (v1, v2, block) 
{
    ublock <- unique(block)
    m <- length(ublock)
    N <- length(block)
    B <- matrix(0, nrow = N, ncol = N)
    ind <- 1
    for (k in 1:m) {
        nk <- sum(block == ublock[k])
        Bk <- diag(rep(v1 - v2, nk),nk) + v2
        B[ind:(ind + nk - 1), ind:(ind + nk - 1)] <- Bk
        ind <- ind + nk
    }
    B
}
