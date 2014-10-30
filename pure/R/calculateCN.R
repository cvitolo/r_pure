#' Calculate the Curve Number from P and Q
#'
#' @param P maximum precipitation for an event (it can also be a vector, for multiple events)
#' @param Q maximum discharge for an event (it can also be a vector, for multiple events)
#'
#' @author Claudia Vitolo
#'
#' @return Curve Number, in the range [0,100]
#'

calculateCN <- function(P,Q){

  if (length(P)!=length(Q)) stop

  S <- 5*P + 10*Q - sqrt( (5*P+10*Q)^2 -25*(P^2-P*Q) )

  CN <- 1000/(S + 10)

  return(CN)

}
