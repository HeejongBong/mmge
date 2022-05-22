##
##  vector_variate_methods.R
##
##  This code implements baseline methods, which are compared 
##  in IPyhton Notebooks at `/EXAMPLE/`. These implementations
##  are not included in package `mmge`. You should manually
##  import the functions by 
##       `source('/EXTERNAL/vector_variate_methods.R')`.
##
##  Created by Heejong Bong and Zongge Liu on 05/21/22.
##  


################################################################
# Spatial covariance estimate by Ren, Kang, Fan, & Lv. (2019). #
################################################################
est.spatial.rkfl<- function(data, lambdas=NULL, 
                            verbose=FALSE, whiten=FALSE){
    # Heterogeneous Group Square-root Lasso is done by package `HGSL`
    require(HGSL)
    
    # Dimension of data
    m = length(data) # No. of sessions
    p = dim(data[[1]])[1] # No. of timesteps
    q = dim(data[[1]])[2] # No. of spatial channels
    
    # No. of trials in each session
    ns = sapply(data, function(x){dim(x)[3]})
    sess.id = c(sapply(1:m, function(l){rep(l,ns[l])}))

    # Default lambda
    if(is.null(lambdas)){
        lambdas = c(1)
    }
    num.lambda = length(lambdas)
    
    # Whitening temporal correlation for vector-variate methods
    if(whiten){
        data = whiten.per.session(data)
    }
    
    e.hat.S = array(0, dim=c(q,p,sum(ns),num.lambda))
    b.hat.S = array(-1, dim=c(q,q,m,num.lambda))
    for(j in 1:q){
        if(verbose){ 
            cat('spatial HGSL at i=', j,'\n') 
            flush.console()
        }
        
        # Stack spatial observations across trials and times
        X.group = matrix(bdiag(lapply(1:m, function(l){
            matrix(aperm(data[[l]][,-j,], c(2,1,3)), 
                   nrow=q-1, ncol=p*ns[l])})), 
            nrow=(q-1)*m, ncol=p*sum(ns))
        Y.group = matrix(unlist(lapply(data, function(x){c(x[,j,])})),
                         nrow=1, ncol=p*sum(ns))

        # Apply Heterogeneous Group Square-root Lasso
        index = Reduce(f='+', rbind(rep(1,m),ns*p-1), accumulate=TRUE)
        b.hat.HGSL = HGSL::S_TISP_Path(
            X=t(X.group), y=t(Y.group), grps=rep(1:(q-1), m), 
            k=m, index=index, lambda=lambdas)
        
        # epsilon.hat
        e.hat.S[j,,,] = c(Y.group) - t(X.group) %*% b.hat.HGSL
        # beta.hat
        b.hat.S[-j,j,,] = b.hat.HGSL
    }    
    
    # List epsilon.hat.S and beta.hat.S across sessions
    e.hat.S = lapply(1:m, function(l){
        array(e.hat.S[,,sess.id==l,], c(q,p,ns[l],num.lambda))
    })
    b.hat.S = lapply(1:m, function(l){
        array(b.hat.S[,,l,], c(q,q,num.lambda))
    })
    
    # Phi.bar
    P.bar.S = lapply(1:m, function(l){
        array(apply(array(e.hat.S[[l]], c(q,p*ns[l],num.lambda)), 3,
            function(x){x %*% t(x) / (p*ns[l])}), c(q,q,num.lambda))
    })
    # Phi.hat
    P.hat.S = lapply(1:m, function(l){
        dgPb.S = array(rep(apply(P.bar.S[[l]],3,diag),q), c(q,num.lambda,q))
        bhPb.S = b.hat.S[[l]] * aperm(dgPb.S, c(1,3,2))

        return(- P.bar.S[[l]] - bhPb.S - aperm(bhPb.S, c(2,1,3)))
    })
    
    # Omega.hat
    W.hat.S = lapply(P.hat.S, function(x){
        dgPh.S = array(rep(apply(x, 3, diag),q), c(q,num.lambda,q))

        return(x / aperm(dgPh.S,c(1,3,2)) / aperm(dgPh.S,c(3,1,2)))
    })
    # Sigma.hat
    S.hat.S = lapply(W.hat.S, function(x){
        array(apply(x, 3, solve), c(q,q,num.lambda))
    })
    return(list(W.hat.S=W.hat.S, S.hat.S=S.hat.S, b.hat.S=b.hat.S, P.hat.S=P.hat.S))
}

inf.single.rkfl = function(W.hat.S, ns){
    # No. of sessions
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
        W.hat.S = array(W.hat.S, c(q,q,1))
    }
    if( length(dim(W.hat.S[[1]])) != 3 ){
        stop("The precision matrices should be given 2- or 3-dimensional")
    }
    num.lambda = dim(W.hat.S[[1]])[3] # No. of penalty parameters
    
    # Test statistics
    T.hat = Reduce(f = "+", lapply(1:m, function(l){
        sqrt(ns[l]/m) * W.hat.S[[l]]
    }), accumulate=FALSE)
    
    # Standard error of test statistics
    Var.T.hat = Reduce(f = "+", lapply(1:m, function(l){
        W = W.hat.S[[l]]
        aperm(apply(array(which(matrix(1,q,q) == 1, arr.ind=TRUE), c(q,q,2)),
              c(1,2), function(e){
            i = e[1]; j = e[2]
            return((1/m) * (W[i,j,]^2 + W[i,i,]*W[j,j,]))
        }), c(2,3,1))
    }), accumulate=FALSE)
    
    # P-values
    p.value = 1 - pnorm(T.hat / sqrt(Var.T.hat))

    return(list(p.value=p.value, T.hat=T.hat, Var.T.hat=Var.T.hat))
}


###############################################################
# Spatial covariance estimate by Cai, Li, Liu, & Xie. (2016). #
###############################################################
est.spatial.cllx<- function(data, lambdas=NULL, 
                            verbose=FALSE, solver="ECOS", whiten=FALSE){
    require(CVXR)
    
    # Dimension of data
    m = length(data) # No. of sessions
    p = dim(data[[1]])[1] # No. of timesteps
    q = dim(data[[1]])[2] # No. of spatial channels
    
    # No. of trials in each session
    ns = sapply(data, function(x){dim(x)[3]})
    sess.id = c(sapply(1:m, function(l){rep(l,ns[l])}))

    # Default lambda
    if(is.null(lambdas)){
        lambdas = c(1)
    }
    num.lambda = length(lambdas)
    
    # Whitening temporal correlation for vector-variate methods
    if(whiten){
        data = whiten.per.session(data)
    }
  
    # Spatial sample covariance matrices
    S.bar.S = lapply(1:m, function(l){
        cor(array(aperm(data[[l]], c(1,3,2)), c(p*ns[l], q)))
    })
    
    W1.hat.S = array(0, dim=c(q,q,m,num.lambda))
    for (i in 1:length(lambdas)){
        # Setting lambda
        lambda.i = lambdas[i]
        if(verbose){ cat('lambda:', lambda.i,'\n') }
        
        for (j in 1:q){
            if(verbose){ 
                cat('spatial estimation using CVXR at i=', j,'\n') 
                flush.console()
            }
            
            # Variable Setting for CVXR
            ej = array(0, dim=c(q,1)); ej[j] = 1
            beta.CVXR = Variable(q, m)
            
            # Objective function
            obj = do.call(max_elemwise, sapply(1:m, function(l){
                cvxr_norm(beta.CVXR[,l], p=1)
            }))
            # Constraints
            constraints = list(
                sum(ns) * lambda.i^2 - Reduce(f='+', lapply(1:m, function(l){
                    ns[l] * (S.bar.S[[l]] %*% beta.CVXR[,l] - ej)^2
                }), accumulate=FALSE) >= 0
            )
            # Minimiation
            ret.CVXR = solve(Problem(Minimize(obj), constraints), solver=solver)
            W1.hat.S[,j,,i] = ret.CVXR$getValue(beta.CVXR)
        }
    }
    
    # beta.hat.S, note that it is not the same beta as in regression based methods
    W1.hat.S = lapply(1:m, function(l){
        array(W1.hat.S[,,l,], c(q,q,num.lambda))
    })
    
    # Omega.hat.S
    W.hat.S = lapply(1:m, function(l){
        pmin(abs(W1.hat.S[[l]]), aperm(abs(W1.hat.S[[l]]), c(2,1,3)),
             c(q,q,num.lambda)) * sign(W1.hat.S[[l]])
    })
    # Sigma.hat
    S.hat.S = lapply(W.hat.S, function(x){
        array(apply(x, 3, solve), c(q,q,num.lambda))
    })
    
    return(list(W1.hat.S=W1.hat.S, W.hat.S=W.hat.S, S.hat.S=S.hat.S))
}


#####################################################
# Spatial covariance estimate by Lee & Liu. (2015). #
#####################################################
est.spatial.ll<- function(data, lambdas1=NULL, lambdas2=NULL, nu=NULL, 
                          verbose=FALSE, solver="ECOS", whiten=FALSE){
    require(CVXR)
    
    # Dimension of data
    m = length(data) # No. of sessions
    p = dim(data[[1]])[1] # No. of timesteps
    q = dim(data[[1]])[2] # No. of spatial channels
    
    # No. of trials in each session
    ns = sapply(data, function(x){dim(x)[3]})
    sess.id = c(sapply(1:m, function(l){rep(l,ns[l])}))

    # Default lambda
    if(is.null(lambdas1)){
        lambdas1 = c(1)
    }
    if(is.null(lambdas2)){
        lambdas2 = c(1)
    }
    num.l1 = length(lambdas1); num.l2 = length(lambdas2)
    
    # Default nu
    if (is.null(nu)){
        nu = 1/m
    }
    
    # Whitening temporal correlation for vector-variate methods
    if(whiten){
        data = whiten.per.session(data)
    }
    
    # Spatial sample covariance matrices
    S.bar.S = lapply(1:m, function(l){
        cor(array(aperm(data[[l]], c(1,3,2)), c(p*ns[l], q)))
    })
    
    R.hat = array(0, dim = c(q,q,m, num.l1,num.l2))
    M.hat = array(0, dim = c(q,q,   num.l1,num.l2))
    for(i in 1:num.l1){ for(j in 1:num.l2){
        l1.i = lambdas1[i]; l2.j = lambdas2[j]
        if(verbose){
            cat('regress lambda1:', l1.i, 'lambda2:', l2.j,'\n')
            flush.console()
        }
        for(k in 1:q){
            ek = array(0, dim=c(q,1)); ek[k] = 1

            beta.CVXR = Variable(q,m)
            m.CVXR = Variable(q)

            obj = (
                cvxr_norm(m.CVXR, p=1) 
                + nu * Reduce(f='+', lapply(1:m, function(l){
                    cvxr_norm(beta.CVXR[,l], p=1)
                }), accumulate=FALSE)
            )
            constraints = c(
                list(m^2 * l1.i^2 - (Reduce(f='+', lapply(1:m, function(l){
                    S.bar.S[[l]] %*% (beta.CVXR[,l] + m.CVXR) - ek
                })))^2 >= 0),
                lapply(1:m, function(l){
                    l2.j^2 - (S.bar.S[[l]] 
                    %*% (beta.CVXR[,l]+m.CVXR)-ek)^2 >= 0
                }),
                list(Reduce(f='+', lapply(1:m, function(l){
                    beta.CVXR[,l] })) == 0
                )
            )
            ret.CVXR = psolve(Problem(Minimize(obj), constraints),
                              solver="ECOS", ignore_dcp=TRUE)

            R.hat[,k,,i,j] = ret.CVXR$getValue(beta.CVXR)
            M.hat[,k,i,j] = ret.CVXR$getValue(m.CVXR)
        }
    }}
    
    R.hat = lapply(1:m, function(l){
        array(R.hat[,,l,,], c(q,q,num.l1,num.l2))
    })
    
    W0.hat.S = lapply(1:m, function(l){
        R.hat[[l]] + M.hat
    })
    
    W.hat.S = lapply(1:m, function(l){
        pmin(abs(W0.hat.S[[l]]), aperm(abs(W0.hat.S[[l]]), c(2,1,3,4)),
             c(q,q,num.l1,num.l2)) * sign(W0.hat.S[[l]])
    })
    
    S.hat.S = lapply(W.hat.S, function(x){
        array(apply(x, c(3,4), solve), c(q,q,num.l1,num.l2))
    })
    
    return(list(W0.hat.S=W0.hat.S, W.hat.S=W.hat.S, S.hat.S=S.hat.S))
}

############################
# data whitening procedure #
############################
whiten.per.session = function(data){
    # Use existing package `whitening`
    require(whitening)
    
    return(lapply(1:m, function(l){
        aperm(array(whitening::whiten(
                array(aperm(data[[l]], c(3,2,1)), c(ns[l]*q,p))),
            c(ns[l], q, p)), c(3,2,1))
    }))
}