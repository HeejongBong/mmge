#' Single edge inference
#'
#' This function performs statistical tests on edge-wise null hypotheses. 
#'
#' @param W.hat.S List with length `m` of 2-dimensional array (`q*q`) or 3-dimensional array (`q*q*num.lambda`) containing the spatial precision matrix estimates over `m` sessions.
#' @param S.hat.T List with length `m` of 2-dimensional array (`p*p`) containing the temporal covariance matrix estimates over `m` sessions.
#' @param ns Number of trials. A numeric for all sessions or a sequence of a different order for each session can be given.
#' @return List output containing the inference result
#'
#'    `R.hat.S` List of the estimated partial correlation matrix at each session.
#'
#'    `p.value` The p-values for the entries of the precision matrix to be zero at each session.
#'
#'    `T.hat` Estimated test statistics \
#'
#'    `Var.T.hat` Estimated variances for the test statistics
#' @export
inf.single.edge <- function(W.hat.S, S.hat.T, ns){
    # No. of sessions
    if( length(W.hat.S) != length(S.hat.T) ){
        stop("The lengths of W.hat.S and S.hat.T do not match")
    }
    m = length(W.hat.S)
    if( length(ns) == 1 ){
        ns = rep(ns,m)
    }
    if( length(ns) != m ){
        stop("The length of ns does not match the number of sessions")
    }
    
    # Dimension of W.hat.S: no.'s of space points and penalty parameters
    if( dim(W.hat.S[[1]])[1] != dim(W.hat.S[[1]])[2] ){
        stop("The precision matrices should be square")
    }
    q = dim(W.hat.S[[1]])[1]
    if( length(dim(W.hat.S[[1]])) == 2 ){
        W.hat.S = lapply(W.hat.S, function(W){array(W, c(q,q,1))})
    } else if( length(dim(W.hat.S[[1]])) != 3 ){
        stop("The precision matrices should be given 2- or 3-dimensional")
    }
    num.lambda = dim(W.hat.S[[1]])[3] # No. of penalty parameters

    # Dimension of S.hat.S: no. of time points
    if( dim(S.hat.T[[1]])[1] != dim(S.hat.T[[1]])[2] ){
        stop("The covariance matrices should be square")
    }
    if( length(dim(S.hat.T[[1]])) != 2 ){
        stop("The covariance matrices should be given 2- dimensional")
    }
    p = dim(S.hat.T[[1]])[1] 
    
    # Rho.hat.S, spatial partial correlation matrix
    R.hat.S = lapply(W.hat.S, function(x){
        dgWh.S = array(rep(apply(x, 3, diag),q), c(q,num.lambda,q))

        return(x / aperm(sqrt(dgWh.S),c(1,3,2)) / aperm(sqrt(dgWh.S),c(3,1,2)))
    })
    
    # Test statistics
    T.hat = Reduce(f = "+", lapply(1:m, function(l){
        sqrt(ns[l]*p/m) * R.hat.S[[l]]
    }), accumulate=FALSE)
    
    # Standard error of test statistics
    Var.T.hat = Reduce(f = "+", lapply(1:m, function(l){
        (sum(S.hat.T[[l]]**2)/p/m) * (1 - R.hat.S[[l]]**2)**2
    }), accumulate=FALSE)
    
    # P-values
    p.value = 2 * (1 - pnorm(abs(T.hat) / sqrt(Var.T.hat)))

    return(list(R.hat.S=R.hat.S, p.value=p.value, T.hat=T.hat, Var.T.hat=Var.T.hat))
}

#' Multiple edge inference
#'
#' This function performs a statistical test on a global null hypothesis over a given
#' edge set.
#'
#' @param W.hat.S List with length `m` of 2-dimensional array (`q*q`) or 3-dimensional array (`q*q*num.lambda`) containing the spatial precision matrix estimates over `m` sessions.
#' @param S.hat.T List with length `m` of 2-dimensional array (`p*p`) containing the temporal covariance matrix estimates over `m` sessions.
#' @param ns Number of trials. A numeric for all sessions or a sequence of a different order for each session can be given.
#' @param E.set 2-dimensional integer array of size `num.E*2` containing the pairs of index for the nodes the target edges connect. If not provided, it is postulated from argument `E.mat`. 
#' @param E.mat 2-dimensional boolen array of size `q*q` containing the adjacency matrix of the target edges. if not provided, it is postulated from arguemnt `E.set`.
#' @param num.bst Integer value for the number of bootstrap iterations. If it is 0, the bootstrap test is not performed.
#' @return List of output containing the inference result.
#'
#'    `R.hat.S` List of the estimated partial correlation matrix at each session.
#'
#'    `T.hat.E` Estimated test statistics 
#'
#'    `Cov.T.E` Estimated covariance matrix for the test statistics
#'
#'    `l.inf.Z` Bootstrap samples of l-infinity norm of zetas
#'
#'    `p.value` p-values for the global null hypothesis. Not returned, if `num.bst` is not positive
#' @export
inf.multiple.edge <- function(W.hat.S, S.hat.T, ns, E.set=NULL, E.mat=NULL, num.bst=0){
    # No. of sessions
    if( length(W.hat.S) != length(S.hat.T) ){
        stop("The lengths of W.hat.S and S.hat.T do not match")
    }
    m = length(W.hat.S)
    if( length(ns) == 1 ){
        ns = rep(ns,m)
    }
    if( length(ns) != m ){
        stop("The length of ns does not match the number of sessions")
    }
    
    # Dimension of W.hat.S: no.'s of space points and penalty parameters
    if( dim(W.hat.S[[1]])[1] != dim(W.hat.S[[1]])[2] ){
        stop("The precision matrices should be square")
    }
    q = dim(W.hat.S[[1]])[1]
    if( length(dim(W.hat.S[[1]])) == 2 ){
        W.hat.S = lapply(W.hat.S, function(W){array(W, c(q,q,1))})
    } else if( length(dim(W.hat.S[[1]])) != 3 ){
        stop("The precision matrices should be given 2- or 3-dimensional")
    }
    num.lambda = dim(W.hat.S[[1]])[3] # No. of penalty parameters

    # Dimension of S.hat.S: no. of time points
    if( dim(S.hat.T[[1]])[1] != dim(S.hat.T[[1]])[2] ){
        stop("The covariance matrices should be square")
    }
    if( length(dim(S.hat.T[[1]])) != 2 ){
        stop("The covariance matrices should be given 2- dimensional")
    }
    p = dim(S.hat.T[[1]])[1] 
    
    # Setting E.mat and E.set
    if( is.null(E.mat) & is.null(E.set) ){
        stop("One of E.mat and E.set should be provided")
    }
    if( is.null(E.mat) ){
        if( dim(E.set)[2] != 2 ){
            stop("Coordinates should be given in E.set")
        }
        
        E.mat = array(FALSE, c(q,q))
        E.mat[E.set] = TRUE
    } else{
        if( dim(E.mat)[1] != dim(E.mat)[2] ){
            stop("E.mat should be a boolean square matrix")
        }
    }
    E.mat = E.mat + t(E.mat)
    E.mat[lower.tri(E.mat, diag = TRUE)] = 0
    E.mat = E.mat != 0
    E.set = which(E.mat, arr.ind=TRUE)

    num.E = sum(1*E.mat)
    
    # Rho.hat.S, spatial partial correlation matrix
    R.hat.S = lapply(W.hat.S, function(x){
        dgWh.S = array(rep(apply(x, 3, diag),q), c(q,num.lambda,q))
        return(x / aperm(sqrt(dgWh.S),c(1,3,2)) / aperm(sqrt(dgWh.S),c(3,1,2)))
    })
    
    Cov.T.E = cov.multiple.edge(R.hat.S, S.hat.T, ns, E.set, E.mat)
    
    # Test statistics
    T.hat = Reduce(f = "+", lapply(1:m, function(l){
        sqrt(ns[l] * p / m) * R.hat.S[[l]]
    }), accumulate=FALSE)
    T.hat.E = array(T.hat[E.mat], c(num.E,num.lambda))
    l.inf.T.E = apply(abs(T.hat.E), 2, max)
    
    # Bootstrap sample
    if( num.bst <= 0 ){
        return(list(R.hat.S=R.hat.S, T.hat.E=T.hat.E, Cov.T.E=Cov.T.E))
    } else {    
        zetas = array(apply(Cov.T.E, 3, function(A){
            MASS::mvrnorm(num.bst, rep(0, num.E), A)}), c(num.bst,num.E,num.lambda))
        l.inf.Z = apply(abs(zetas), c(1, 3), max)
        p.value = apply(t(l.inf.Z) > l.inf.T.E, 1, mean)

        return(list(R.hat.S=R.hat.S, T.hat.E=T.hat.E, Cov.T.E=Cov.T.E,
                    l.inf.Z=l.inf.Z, p.value=p.value))
    }
}

#' Covariance of multiple edge test statistics
#'
#' This function calculates the covariance matrix of test statistics for multiple edge
#' inference.
#'
#' @param W.S List with length `m` of 2-dimensional array (`q*q`) or 3-dimensional array (`q*q*num.lambda`) containing the spatial precision or partial-correlation matrices over `m` sessions.
#' @param S.T List with length `m` of 2-dimensional array (`p*p`) containing the temporal covariance matrices over `m` sessions.
#' @param ns Number of trials. A numeric for all sessions or a sequence of a different order for each session can be given.
#' @param E.set 2-dimensional integer array of size `num.E*2` containing the pairs of index for the nodes the target edges connect. If not provided, it is postulated from argument `E.mat`. 
#' @param E.mat 2-dimensional boolen array of size `q*q` containing the adjacency matrix of the target edges. if not provided, it is postulated from arguemnt `E.set`.
#' @return Estimated covariance matrix for the test statistics
#' @export
cov.multiple.edge <- function(W.S, S.T, ns, E.set=NULL, E.mat=NULL){
    # No. of sessions
    if( length(W.S) != length(S.T) ){
        stop("The lengths of W.hat.S and S.hat.T do not match")
    }
    m = length(W.S)
    if( length(ns) == 1 ){
        ns = rep(ns,m)
    }
    if( length(ns) != m ){
        stop("The length of ns does not match the number of sessions")
    }
    
    # Dimension of W.hat.S: no.'s of space points and penalty parameters
    if( dim(W.S[[1]])[1] != dim(W.S[[1]])[2] ){
        stop("The precision matrices should be square")
    }
    q = dim(W.S[[1]])[1]
    if( length(dim(W.S[[1]])) == 2 ){
        W.S = lapply(W.S, function(W){array(W, c(q,q,1))})
    } else if( length(dim(W.S[[1]])) != 3 ){
        stop("The precision matrices should be given 2- or 3-dimensional")
    }
    num.lambda = dim(W.S[[1]])[3] # No. of penalty parameters

    # Dimension of S.hat.S: no. of time points
    if( dim(S.T[[1]])[1] != dim(S.T[[1]])[2] ){
        stop("The covariance matrices should be square")
    }
    if( length(dim(S.T[[1]])) != 2 ){
        stop("The covariance matrices should be given 2- dimensional")
    }
    p = dim(S.T[[1]])[1] 
    
    # Setting E.mat and E.set
    if( is.null(E.mat) & is.null(E.set) ){
        stop("One of E.mat and E.set should be provided")
    }
    if( is.null(E.mat) ){
        if( dim(E.set)[2] != 2 ){
            stop("Coordinates should be given in E.set")
        }
        
        E.mat = array(FALSE, c(q,q))
        E.mat[E.set] = TRUE
    } else{
        if( dim(E.mat)[1] != dim(E.mat)[2] ){
            stop("E.mat should be a boolean square matrix")
        }
    }
    E.mat = E.mat + t(E.mat)
    E.mat[lower.tri(E.mat, diag = TRUE)] = 0
    E.mat = E.mat != 0
    E.set = which(E.mat, arr.ind=TRUE)
    
    num.E = sum(1*E.mat)
    ExE = array(c(aperm(array(rep(E.set, num.E), c(num.E,2,num.E)), c(1,3,2)),
                  aperm(array(rep(E.set, num.E), c(num.E,2,num.E)), c(3,1,2))),
                c(num.E, num.E, 4))
    
    # Rho.hat.S, spatial partial correlation matrix
    R.S = lapply(W.S, function(x){
        dgW.S = array(rep(apply(x, 3, diag),q), c(q,num.lambda,q))
        return(x / aperm(sqrt(dgW.S),c(1,3,2)) / aperm(sqrt(dgW.S),c(3,1,2)))
    })
    
    # Covariance of test statistics
    Cov.T.E = Reduce(f='+', lapply(1:m, function(l){
        (sum(S.T[[l]]**2)/p/m) * array(apply(ExE, c(1,2), function(x){
            i1 = x[1]; j1 = x[2]; i2 = x[3]; j2 = x[4]
            R = R.S[[l]]
            return(
                R[i1,i2,]*R[j1,j2,] + R[i1,j2,]*R[i2,j1,] + 0.5*R[i1,j1,]*R[i2,j2,]
                  *(R[i1,i2,]^2 + R[j1,j2,]^2 + R[i1,j2,]^2 + R[i2,j1,]^2)
                - R[i1,i2,]*R[i2,j2,]*R[i2,j1,] - R[i1,i2,]*R[i1,j1,]*R[i1,j2,]
                - R[j1,j2,]*R[i2,j2,]*R[i1,j2,] - R[j1,j2,]*R[i2,j1,]*R[i1,j1,])
        }), c(num.E,num.E,num.lambda))
    }), accumulate=FALSE)
    
    return(array(apply(Cov.T.E, 3, function(A){ as.matrix(Matrix::nearPD(A)$mat) }),
                 c(num.E, num.E, num.lambda) ))
}