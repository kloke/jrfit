print.summary.jrfit <-
function (x, digits = max(5, .Options$digits - 2), ...) 
{
    cat("\nCoefficients:\n")
#    est <- x$coef
#    print(format(est, digits = digits), quote = F)
	printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)
#	printCoefmat(x$coefficents,P.values=TRUE,has.Pvalue=TRUE)
}
