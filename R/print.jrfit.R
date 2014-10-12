print.jrfit <-
function (x, ...) 
{
    coef <- coef(x)
    cat("\nCoefficients:\n")
    print(coef, ...)
}
