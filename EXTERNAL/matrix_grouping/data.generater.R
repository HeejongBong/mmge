source("util_funcs.R")
add.edges.on.subgraph <- function(theta,varying_cluster_size=NULL,prob=NULL)
{
  d = dim(theta)[1]
  if(is.null(varying_cluster_size)) varying_cluster_size = max(floor(d/10),3)
  #idx <- sample(1:d,varying_cluster_size)
  idx <- 1:varying_cluster_size
  if(is.null(prob))  prob = .5
  a <- theta[idx,idx]
  tmp = matrix(runif(varying_cluster_size^2,0,0.5),varying_cluster_size,varying_cluster_size)
  tmp = tmp + t(tmp)
  a[tmp < prob] = 1
  #theta[tmp >= tprob] = 0
  theta[idx,idx] <- a
  rm(tmp,a)
  gc()
  theta
}

add.random.edges <- function(theta, prob=NULL){
  d = dim(theta)[1]
  if(is.null(prob))  prob = 4/d
  tmp = matrix(runif(d^2,0,0.5),d,d)
  tmp = tmp + t(tmp)
  theta[tmp < prob] = 1
  rm(tmp)
  gc()
  theta
}


###########  main function generate a specified graph for one group, and add edges for  #########################
##########   other groups using various graph modification schemes                     ##########################
multi.graph.generator <- function(m.time,n.loc,temporal.graph=NULL,spatial.graph=NULL,graph.modification.scheme=NULL,sub.graph.size=NULL,prob=NULL){
  if (is.null(temporal.graph)) { temporal.graph = "band" }
  if (is.null(spatial.graph)) { spatial.graph="random" }
  if (spatial.graph == "random" || spatial.graph == "scale-free") {  graph.modification.scheme =  "densify.subgraph" }
  else if (spatial.graph == "band" || spatial.graph == "hub")  {  graph.modification.scheme = "add.random.edges"  }
  else { stop("graph type ", spatial.graph, " not supported!") }
  if (graph.modification.scheme == "densify.subgraph"){
    if (is.null(sub.graph.size)) { sub.graph.size = max(ceiling(n.loc/10),10) }
    if (is.null(prod)) { prob = .5 }
  } 
  if (graph.modification.scheme == "add.random.edges")  { if (is.null(prod))  { prob = 4/n.loc } } 
  K <- length(m.time)
  S.time <- list() # list of temporal cov matrix
  S.loc <- list() # list of spatial cov matrix
  omega.loc <- list()
  graph.loc <- list()
  common.graph.loc <- graph.generator(n.loc,graph=spatial.graph)
  for (k in 1:K){
    ## generate temporal cov matrix ##
    graph.time <- graph.generator(m.time[k],graph=temporal.graph)$theta
    S.time[[k]] <- omega.generator.simple(graph.time)$sigma
    ## generate spatial cov matrix ##
    if (k == 1) { 
      graph.loc[[1]] = common.graph.loc$theta
      tmp <- omega.generator.simple(graph.loc[[1]])
      S.loc[[1]] <- tmp$sigma
      omega.loc[[1]] <- tmp$omega
    }
    else {
      if (graph.modification.scheme == "densify.subgraph"){
        graph.loc[[k]] = add.edges.on.subgraph(common.graph.loc$theta,varying_cluster_size=sub.graph.size,prob=prob)
      }
      if (graph.modification.scheme == "add.random.edges"){
        graph.loc[[k]] = add.random.edges(common.graph.loc$theta,prob=prob)
      }
      tmp <- omega.generator.simple(graph.loc[[k]])
      S.loc[[k]] <- tmp$sigma
      omega.loc[[k]] <- tmp$omega
    }
  }
  list(omega.loc=omega.loc,graph.loc=graph.loc,S.time=S.time,S.loc=S.loc)
}


data.generator <- function(sample.size,S.time,S.loc){
  # if (length(sample.size) != length(m.time)) { stop("inconsistent input dimension! \n") }
  K <- length(sample.size)
  data = list()
  for (k in 1:K)
  {
    data[[k]] <- rmn_rep(sample.size[k],S.time[[k]],S.loc[[k]])
  }
  data
}