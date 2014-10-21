calculateCN <- function(P,Q){

  if (length(P)!=length(Q)) stop

  myCN <- rep(NA,length(P))

  for (r in 1:length(P)){

    # f <- function (CN) ((0.371-0.2*(1000/CN - 10))^2) / (0.371+0.8*(1000/CN - 10)) - 0.03

    f <- function (CN) ((P[r]-0.2*(1000/CN - 10))^2) / (P[r]+0.8*(1000/CN - 10)) - Q[r]

    # uniroot(f, c(0, 100), myP = P[r], myQ = Q[r], tol = 0.1)
    # temp <- try(uniroot(f, c(0, 100), myP = P[r], myQ = Q[r], tol = 0.1)$root, TRUE)
    # curve(f, 0, 100)
    x <- c()

    for (c in 1:101) x[c] <- f(c-1)

    minRange <- min(which(!is.na(x) & x>0)-1)

    maxRange <- max(which(!is.na(x) & x<0)-1)

    print(paste("CN for event n.",r,"is calculated in the range from",minRange, "to", maxRange))

    temp <- round(uniroot(f, c(minRange, maxRange), tol = 0.10)$root,0)

    if (is.numeric(temp)) myCN[r] <- temp

  }

  return(myCN)

}
