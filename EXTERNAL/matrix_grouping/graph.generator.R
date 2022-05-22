dyn.load("lib/SFGen.so")
graph.generator <- function(d = 50, graph = "random", g = NULL, prob = NULL, vis = FALSE, verbose = FALSE){  
  gcinfo(FALSE)
  if(verbose) cat("Generating data from the multivariate normal distribution with the", graph,"graph structure....")
  if(is.null(g)){
    g = 1
    if(graph == "hub" || graph == "cluster"){
      if(d > 40)	g = ceiling(d/20)
      if(d <= 40) g = 2
    }
  }
  
  if(graph == "random"){
    if(is.null(prob))	prob = min(1, 3/d)
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }
  
  if(graph == "cluster"){
    if(is.null(prob)){
      if(d/g > 30)	prob = 0.3
      if(d/g <= 30)	prob = min(1,6*g/d)
    }
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }  
  
  
  # parition variables into groups
  g.large = d%%g
  g.small = g - g.large
  n.small = floor(d/g)
  n.large = n.small+1
  g.list = c(rep(n.small,g.small),rep(n.large,g.large))
  g.ind = rep(c(1:g),g.list)
  rm(g.large,g.small,n.small,n.large,g.list)
  gc()
  
  # build the graph structure
  theta = matrix(0,d,d);
  if(graph == "band"){
    for(i in 1:g){
      diag(theta[1:(d-i),(1+i):d]) = 1
      diag(theta[(1+i):d,1:(d-1)]) = 1
    }	
  }
  if(graph == "cluster"){
    for(i in 1:g){
      tmp = which(g.ind==i)
      tmp2 = matrix(runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
      tmp2 = tmp2 + t(tmp2)		 	
      theta[tmp,tmp][tmp2<prob] = 1
      rm(tmp,tmp2)
      gc()
    }
  }
  if(graph == "hub"){
    for(i in 1:g){
      tmp = which(g.ind==i)
      theta[tmp[1],tmp] = 1
      theta[tmp,tmp[1]] = 1
      rm(tmp)
      gc()
    }
  }
  if(graph == "random"){
    tmp = matrix(runif(d^2,0,0.5),d,d)
    tmp = tmp + t(tmp)
    theta[tmp < prob] = 1
    #theta[tmp >= tprob] = 0
    rm(tmp)
    gc()
  }
  
  if(graph == "scale-free"){
    out = .C("SFGen",dd0=as.integer(2),dd=as.integer(d),G=as.integer(theta),package="huge")
    theta = matrix(as.numeric(out$G),d,d)
  }
  diag(theta) = 0
  
  if(verbose) cat("done.\n")
  rm(verbose)
  gc()
  
  sim_g = list(theta = Matrix(theta,sparse = TRUE), sparsity= sum(theta)/(d*(d-1)), graph.type=graph)
  class(sim_g) = "sim_g" 
  return(sim_g)
}


omega.generator.simple = function(theta,v = NULL, u = NULL,vis=FALSE){
  d = dim(theta)[1]
  if(is.null(u)) u = 0.1
  if(is.null(v)) v = 0.3
  omega = theta*v
  # make omega positive definite and standardized
  diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
  sigma = cov2cor(solve(omega))
  omega = matrix(solve(sigma))
  omega[abs(omega)<1e-10] = 0
  omega = matrix(omega,d,d)
  sigma = solve(omega)
  return(list(omega=omega,sigma=sigma))
}

omega.generator = function(theta,v = NULL, u = NULL,vis=FALSE){
  d = dim(theta)[1]
  if(is.null(u)) u = 0.1
  if(is.null(v)) v = 0.3
  omega = theta*v
  # make omega positive definite and standardized
  diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
  sigma = cov2cor(solve(omega))
  omega = matrix(solve(sigma))
  omega[abs(omega)<1e-10] = 0
  omega <- Matrix(matrix(omega,d,d))
  
  # generate multivariate normal data
  x = mvrnorm(n,rep(0,d),sigma)
  sigmahat = cor(x)
  
  # graph and covariance visulization
  if(vis == TRUE){
    fullfig = par(mfrow = c(1, 3), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
    fullfig[1] = image(theta, col = gray.colors(256),  main = "Adjacency Matrix")
    
    fullfig[2] = image(sigma, col = gray.colors(256), main = "Covariance Matrix")
    g = graph.adjacency(theta, mode="undirected", diag=FALSE)
    layout.grid = layout.fruchterman.reingold(g)
    
    fullfig[3] = plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA,main = "Graph Pattern")
    
    rm(fullfig,g,layout.grid)
    gc()
  }
  sim = list(data = x, sigma = sigma, sigmahat = sigmahat, omega = omega, theta = Matrix(theta,sparse = TRUE), sparsity= sum(theta)/(d*(d-1)))
  class(sim) = "sim" 
  return(sim)
}

