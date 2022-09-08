#' Estimate spatial covariance matrices
#'
#' This function estimates the spatial covariance parameter from  
#' matrix-variate observation data over multiple sessions. 
#'
#' @param data List with length `m` of 3-dimensional array (`p*q*n.l`) containing `n.l` matrix observations of dimension `p*q` over sessions `l` = 1, ..., m.
#' @param lambdas Sequence of group lasso penalty hyperparameters. Default is determined by `p`, `q`, `m`, and `n.l`'s. 
#' @param eps Threshold value for convergence. Default is 1e-4.
#' @param maxit Maximum number of iteration. Default is 10,000.
#' @param verbose Boolean with default `FALSE`.
#' @return List output containing the spatial estimation result
#' 
#'   `W.hat.S` List of Omega.hat.S, the estimate of the spatial precision matrix, at each session.
#'
#'   `S.hat.S` List of Sigma.hat.S, the estimate of the spatial covariance matrix, at each session.
#'
#'   `b.hat.S` List of beta.hat.S, the estimate of the spatial regression parameters, at each session.
#'
#'   `P.hat.S` List of Phi.hat.S, the estimate of the spatial residual covariance, at eachs ession.
#' @export
est.spatial<- function(data, lambdas=NULL, eps=1e-4, maxit=1e+4, verbose=FALSE){
    # Dimension of data
    m = length(data) # No. of sessions
    p = dim(data[[1]])[1] # No. of timesteps
    q = dim(data[[1]])[2] # No. of spatial channels
    
    # No. of trials in each session
    ns = unlist(lapply(data, function(x){dim(x)[3]}))
    sess.id = unlist(lapply(1:m, function(l){rep(l,ns[l])}))

    # Default lambda
    if(is.null(lambdas)){
        lambdas = c(sqrt((m+log(q*m*min(ns)*p))/(min(ns)*p)))
    }
    num.lambda = length(lambdas)

    e.hat.S = array(0, dim=c(q,p,sum(ns),num.lambda))    # epsilon.hat.S
    b.hat.S = array(-1, dim=c(q,q,m,num.lambda))         # beta.hat.S
    for(j in 1:q){
        if(verbose){ 
            cat('spatial group lasso at i=', j,'\n') 
            flush.console()
        }

        # Stack spatial observations across trials and times
        sd.X.group = lapply(data, function(x){
            sqrt(apply(x[,-j,]**2, 2, mean))})
        X.tilde.group = matrix(bdiag(lapply(1:m, function(l){
            matrix(aperm(data[[l]][,-j,], c(2,1,3)), 
                   nrow=q-1, ncol=p*ns[l]) / sd.X.group[[l]]})), 
            nrow=(q-1)*m, ncol=p*sum(ns))
        Y.group = matrix(unlist(lapply(data, function(x){c(x[,j,])})),
                         nrow=1, ncol=p*sum(ns))
        
        # Group Lasso
        ret.gglasso = gglasso::gglasso(
            x=t(X.tilde.group), y=t(Y.group), group=rep(1:(q-1), m), lambda=lambdas,
            eps = eps, maxit = maxit, intercept=FALSE)
        
        # epsilon.hat.S
        e.hat.S[j,,,] = c(Y.group) - predict(ret.gglasso, t(X.tilde.group))
        # beta.hat.S
        b.hat.S[-j,j,,] = ret.gglasso$beta / unlist(sd.X.group)
    }    
    
    # List epsilon.hat.S and beta.hat.S across sessions
    e.hat.S = lapply(1:m, function(l){
        array(e.hat.S[,,sess.id==l,], c(q,p,ns[l],num.lambda))
    })
    b.hat.S = lapply(1:m, function(l){
        array(b.hat.S[,,l,], c(q,q,num.lambda))
    })
    
    #Phi.bar.S
    P.bar.S = lapply(1:m, function(l){
        array(apply(array(e.hat.S[[l]], c(q,p*ns[l],num.lambda)), 3,
            function(x){x %*% t(x) / (p*ns[l])}), c(q,q,num.lambda))
    })
    #Phi.hat.S
    P.hat.S = lapply(1:m, function(l){
        dgPb.S = array(rep(apply(P.bar.S[[l]],3,diag),q), c(q,num.lambda,q))
        bhPb.S = b.hat.S[[l]] * aperm(dgPb.S, c(1,3,2))

        return(- P.bar.S[[l]] - bhPb.S - aperm(bhPb.S, c(2,1,3)))
    })
    
    #Omega.hat.S
    W.hat.S = lapply(P.hat.S, function(x){
        dgPh.S = array(rep(apply(x,3,diag),q), c(q,num.lambda,q))

        return(x / aperm(dgPh.S,c(1,3,2)) / aperm(dgPh.S,c(3,1,2)))
    })
    #Sigma.hat.S
    S.hat.S = lapply(W.hat.S, function(x){
        array(apply(x,3,solve), c(q,q,num.lambda))
    })
    
    return(list(W.hat.S=W.hat.S, S.hat.S=S.hat.S, b.hat.S=b.hat.S, P.hat.S=P.hat.S))
}

#' Estimate temporal covariance matrices
#'
#' This function estimates the temporal covariance parameter from  
#' matrix-variate observation data over multiple sessions. 
#'
#' @param data List with length `m` of 3-dimensional array (`p*q*n`) containing `n` matrix observations of dimension `p*q` over `m` sessions
#' @param order.T Order for the temporal AR models. A numeric for all sessions or a sequence of a different order for each session can be given. If not provided, the default value is determined using `decay.T
#' @param decay.T Decay rate of AR parameter. A numeric for all sessions or a sequence of a different rate for each session can be given. Default is 0.5.
#' @param verbose Boolean with default `FALSE`
#' @return List output containing the temporal estimation result
#'
#'   `W.hat.T` List of Omega.hat.T, the estimate of the temporal precision matrix, at each session.
#'
#'   `S.hat.T` List of Sigma.hat.T, the estimate of the temporal covariance matrix, at each session.
#'
#'   `b.hat.T` List of beta.hat.T, the estimate of the temporal regression parameters, at each session.
#' 
#'   `P.hat.T` List of Phi.hat.T, the estimate of the temporal residual variances, at each session
#' @export
est.temporal <- function(data, order.T=NULL, decay.T=NULL, verbose=FALSE){
    # Dimension of data
    m = length(data) # No. of sessions
    p = dim(data[[1]])[1] # No. of time points
    q = dim(data[[1]])[2] # No. of space points
    
    # No. of trials in each session
    ns = unlist(lapply(data, function(x){dim(x)[3]}))
    sess.id = unlist(lapply(1:m, function(l){rep(l,ns[l])}))
    
    # Default values for decay.T and order.T
    if(is.null(decay.T)){
        decay.T = rep(0.5, m)
    }
    if(is.null(order.T)){
        order.T =  floor((ns*q)^(1/(2*decay.T+2)))
    }
    if(length(order.T)==1){
        order.T = rep(order.T, m)
    } 
    if(length(order.T) != m){
        stop("The length of order.T should match that of data")
    }
  
    e.hat.T = array(0, dim=c(p,q,sum(ns)))
    b.hat.T = array(0, dim=c(p,p,m))
    
    # Trivial calculation on t=1
    if(verbose){ 
        cat('temporal AR fit at t=', 1,'\n') 
        flush.console()
    }
    e.hat.T[1,,] = unlist(lapply(data, function(x){x[1,,]}))
    for(t in 2:p){
        if(verbose){ 
            cat('temporal AR fit at t=', t,'\n') 
            flush.console()
        }
        # Stack spatial observations across trials and times
        X.AR = lapply(1:m, function(l){
            matrix(data[[l]][max(1,t-order.T[l]):(t-1),,],
            nrow=min(t-1,order.T[l]), ncol=q*ns[l])})
        Y.AR = lapply(1:m, function(l){matrix(data[[l]][t,,],
            nrow=1, ncol=q*ns[l])})
        
        # Fit autoregressive model 
        XX.AR = lapply(1:m, function(l){
            X.AR[[l]] %*% t(X.AR[[l]])})
        XY.AR = lapply(1:m, function(l){
            X.AR[[l]] %*% t(Y.AR[[l]])})

        b.hat.AR = lapply(1:m, function(l){
            MASS::ginv(XX.AR[[l]]) %*% XY.AR[[l]]})
        e.hat.AR = lapply(1:m, function(l){
            Y.AR[[l]] - t(b.hat.AR[[l]]) %*% X.AR[[l]]})

        e.hat.T[t,,] = unlist(e.hat.AR)
        b.hat.T[1:(t-1),t,] = unlist(lapply(1:m, function(l){
            b.hat.l = array(0, c(t-1))
            b.hat.l[max(1,t-order.T[l]):(t-1)] = b.hat.AR[[l]]
            return(b.hat.l)
        }))
    }
    
    # List epsilon.hat.T and beta.hat.T across sessions
    e.hat.T = lapply(1:m, function(l){
        e.hat.T[,,sess.id==l]
    })
    b.hat.T = lapply(1:m, function(l){
        b.hat.T[,,l]
    })
    
    # Phi.hat.T
    P.hat.T = lapply(1:m, function(l){
        apply(e.hat.T[[l]]**2, 1, sum) / (ns[l]*q - pmin((1:p)-1,order.T[l]))  
    })
    # Sigma.hat.T
    S.hat.T = lapply(1:m, function(l){
        invIA = MASS::ginv(diag(p) - b.hat.T[[l]])
        S.hat.l = t(invIA) %*% diag(P.hat.T[[l]]) %*% invIA
        return(S.hat.l / sum(diag(S.hat.l)) * p)
    })
    # Omega.hat.T
    W.hat.T = lapply(1:m, function(l){
        MASS::ginv(S.hat.T[[l]])
    })
  
    return(list(W.hat.T=W.hat.T, S.hat.T=S.hat.T, b.hat.T=b.hat.T, P.hat.T=P.hat.T))
}