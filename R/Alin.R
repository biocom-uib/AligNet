#'Hungarian
#'@keywords internal
HungarianFinal <- function(mat,maxim = TRUE) {
  nodesnet1 <- rownames(mat)
  nodesnet2 <- colnames(mat)
  if (dim(mat)[1] > dim(mat)[2]) {
    alin <- solve_LSAP(t(mat),maximum = maxim)
    alin <- cbind(nodesnet1[alin],nodesnet2[seq_along(alin)])
  }
  else{
    alin <- solve_LSAP(mat,maximum = maxim)
    alin <- cbind(nodesnet1[seq_along(alin)],nodesnet2[alin])
  }
  return(alin)
}


#'Edge Correctness Score
#'
#' Given two networks, \code{net1}=(\eqn{V_1,E_1}) , \code{net2}=(\eqn{V_2,E_2})
#' , and an alignment \eqn{g:V_1 \rightarrow V_2}, then the edge correctnees is
#' defined as: \deqn{\frac{|\{(u,v)\in E_1 \: : \;
#' (g(u),g(v))\in E_2\}|}{|E_1|}}
#'@param alin  alignment, with the format of \code{alin.local} or
#'\code{alin.global}
#'@param net1  network
#'@param net2  network
#'@return EC score
EC.score <- function(alin,net1,net2) {
  alini <- function(x,y) {
    matrix(c(as.character(alin[x])
             ,as.character(alin[y])),nrow = 1)
  }
  E1 <- get.edgelist(net1)
  if (dim(E1)[1] == 0) {
    return(0)
  }
  eds <- t(mapply(alini,E1[,1],E1[,2]))
  nas <- unique(c(which(is.na(eds[,2])),which(is.na(eds[,1]))))
  if (length(nas) > 0) {
    eds <- eds[ - nas,]
  }
  if (is.null(dim(eds))) {
    if (length(eds) == 2) {
      eds <- cbind(eds[1],eds[2])
    }
  }
  if (dim(eds)[1] == 0) {
    return(0)
  }
  else{
    if (dim(eds)[2] == 1) {
      return(0)
    }
  }
  Gnet3 <- graph.edgelist(eds,directed = FALSE)
  ggg4 <- graph.intersection(Gnet3,net2,byname = TRUE)
  Gnet1 <- induced.subgraph(net1,vids = names(alin))
  Gnet2 <- induced.subgraph(net2,vids = alin)
  if (min(length(E(Gnet1)),length(E(Gnet2))) == 0) {
    return(0)
  }
  return(length(E(ggg4)) / min(length(E(Gnet1)),length(E(Gnet2))))

}

#'Functional Coherence Score
#'
#'Given an alignment \eqn{g:V_1 \rightarrow V_2}, and an ontologies \eqn{gos},
#'then the functional coherence score is defined as:
#'\deqn{\frac{1}{|V_1|} \sum_{v \in V_1} \frac{|gos(v) \cap gos(g(v))|}
#'{|gos(v)\cup gos(g(v))|}}
#'@param alin alignment, with the format of \code{alin.local} or
#'\code{alin.global}
#'@param gos a list of ontologies
#'@return FC score
FC.score=function(alin,gos=go){
    aa1 <- names(alin)
    aa2 <- alin
    fcs <- unlist(lapply(1:length(aa1), 
                         FUN=function(i) {length(intersect(gos[[aa1[i]]],
                           gos[[aa2[i]]]))/length(union(gos[[aa1[i]]], gos[[aa2[i]]]))}))
    
    return(mean(fcs, na.rm = TRUE))
  }

#'Local alignment
#'
#'Compute the local alignment between the networks \code{net1} and
#'\code{net2} with centers \code{p1} and \code{p2}
#'@param net1 network
#'@param net2 network
#'@param p1 center of net1
#'@param p2 center of net2
#'@param compute.ec compute ecscore (TRUE/FALSE)
#'@param mat a non dissimilarity matrix
#'@return alignment
align.local <-
  function(net1,net2,p1,p2,compute.ec = FALSE, mat = NULL) {
    mat1 <- compute.matrix.Degree(net1,net2)
    if (is.null(mat)) {
      mat <- mat1
    }
    else{
      mat <- mat[V(net1)$name,V(net2)$name] + mat1
    }
    dimnames(mat) <- list(V(net1)$name,V(net2)$name)
    neigh1 <- neighborhood(graph = net1, order = 1)
    neigh1 <-
      lapply(1:length(V(net1)),function(i)
        V(net1)$name[neigh1[[i]]])
    names(neigh1) <- V(net1)$name
    neigh2 <- neighborhood(graph = net2,order = 1)
    neigh2 <-
      lapply(1:length(V(net2)),function(i)
        V(net2)$name[neigh2[[i]]])
    names(neigh2) <- V(net2)$name
    assign <- p2
    names(assign) <- c(p1)
    completes <- c()
    incomplets <- c(p1)
    assignats <- c(p2)
    while (length(assign) < length(neigh1)) {
      q1 <- incomplets[1]
      q2 <- assign[q1]
      n1 <- setdiff(neigh1[[q1]],names(assign))
      n2 <- setdiff(neigh2[[q2]],assign)
      assign2 <- c()

      if (length(n2) > 1 && length(n1) > 1) {
        mat2 <- mat[n1,n2]
        dimnames(mat2) <- list(n1,n2)

      }
      else{
        mat2 <- c()
        if (length(n1) * length(n2) > 0) {
          if (length(n1) == 1 && length(n2) == 1) {
            assign2 <- n2
            names(assign2) <- n1
            assign <- c(assign,assign2)

          }
          else{
            if (length(n1) == 1) {
              assign2 <- names(which.min(mat[n1,n2]))
              names(assign2) <- n1
              assign <- c(assign,assign2)
            }
            if (length(n2) == 1) {
              assign2 <- n2
              names(assign2) <- names(which.min(mat[n1,n2]))
              assign <- c(assign,assign2)
            }
          }
        }
      }

      if (!is.null(dim(mat2))) {
        if (dim(mat2)[1] > dim(mat2)[2]) {
          hung <- HungarianFinal(as.matrix(t(mat2)),maxim = FALSE)
          assign2 <- hung[,1]
          names(assign2) <- hung[,2]
          inds <- sort(unlist(lapply(1:dim(hung)[1],
                                     function(i)
                                       t(mat2)[hung[i,1],hung[i,2]])),
                       decreasing = TRUE,index.return = TRUE)$ix
          assign2 <- assign2[inds]
          assign <- c(assign,assign2)

        }
        else{
          hung <- HungarianFinal(as.matrix(mat2),maxim = FALSE)
          assign2 <- hung[,2]
          names(assign2) <- hung[,1]
          inds <- sort(unlist(lapply(1:dim(hung)[1],
                                     function(i)
                                       (mat2)[hung[i,1],hung[i,2]])),
                       decreasing = TRUE,index.return = TRUE)$ix
          assign2 <- assign2[inds]
          assign <- c(assign,assign2)
        }
      }

      incomplets <- incomplets[ - 1]
      incomplets <- c(incomplets,names(assign2))


      if (length(incomplets) == 0) {
        break
      }
    }
    if (length(assign) == 0) {
      assign <- p2
      names(assign) <- p1
    }
    if (compute.ec) {
      return(list(align = assign,ec = EC.score(assign,net1,net2)))
    }
    return(list(align = assign))
  }

#'alin_aux
#'@keywords internal
alin_aux <- function(p1,mat,ll,clust1,clust2,mat2 = NULL) {
  pp2 <- intersect(names(which(mat[p1,] > ll)),names(clust2))
  protsc1 <- V(clust1[[p1]])$name
  if (length(pp2) == 0) {
    return(list())
  }
  return(lapply(pp2,
                function(p2)
                  align.local(clust1[[p1]],clust2[[p2]],p1,p2,mat = mat2)))

}

#'All local alignments
#'
#'Compute a list of local alignments between \code{clust1} and
#'\code{clust2}. The function compute all the local alignments of all pairs
#'of clusters whose centers have a similarity greather than \code{ll}
#'@param clust1 a list of clusters
#'@param clust2 a list of clusters
#'@param mat a similarity matrix
#'@param threshold a threshold
#'@param cores number of cores
#'@param dismat a disimilarity matrix to use in the local aligments
#'@return a list of alignments
align.local.all <-
  function(clust1,clust2,mat,threshold,cores = 2, dismat = NULL) {
    prots1 <- intersect(names(clust1),rownames(mat))
    aligns <- mclapply(prots1,
                       function(i)
                         alin_aux(i,mat,threshold,clust1,clust2,dismat),
                       mc.cores = cores)
    return(aligns)
  }


#'Size score
#'
#'Compute the size score of the list of local alignments. Given an alignment
#'\eqn{g:V_1 \rightarrow V_2}, the size score of \eqn{g} is defined as
#'\eqn{|V_1|}
#'@param localAligns a list of local alignments
#'@return sizeScore a table with 3 columns, the first and second ones
#'represents the local alignment and the third the score of alignment
size.score.all <- function(localAligns) {
  als <- localAligns
  align.size <- unlist(lapply(als,function(i)
    length(i)))
  align.name1 <- unlist(lapply(als,function(i)
    names(i)[1]))
  align.name2 <- unlist(lapply(als,function(i)
    i[1]))

  return(cbind(align.name1,align.name2,align.size))
}


#'Similarity Score
#'
#'Compute the similarity score of the list of local alignments. Given an
#'alignment \eqn{g:V_1\rightarrow V_2}, and a similarity matrix \eqn{sim}, the
#'similarity score of the alignment \eqn{g} is defined as:
#'\deqn{\frac{1}{|V_1|} \sum_{v \in V_1} sim[v,g(v)]}
#'@param localAligns a list of local alignments
#'@param sim a similarity matrix
#'@return simScore a table with 3 columns, the first and second ones
#'represents the local alignment and the third the score of alignment
sim.score.all <- function(localAligns,sim) {
  als <- localAligns
  align.sim <- unlist(lapply(als,function(i)
    sim.score(i,sim)))
  align.name1 <- unlist(lapply(als,function(i)
    names(i)[1]))
  align.name2 <- unlist(lapply(als,function(i)
    i[1]))

  return(cbind(align.name1,align.name2,align.sim))
}

#'simScore.aux
#'@keywords internal
sim.score <- function(align,sim) {
  sco <- sum(diag(sim[names(align),align]))
  n <- length(align)
  return(sco / n)
}

#'Edges score
#'
#'Compute the EC.score of the list of local alignments
#'@param localAligns, a list of local alignments
#'@param net1 an igraph object
#'@param net2 an igraph object
#'@return edgesScore  a table with 3 columns, the first and second ones
#'represents the local alignment and the third the score of alignment
EC.score.all <- function(localAligns,net1,net2) {
  als <- localAligns
  align.ec <- unlist(lapply(als,function(i)
    EC.score(i,net1,net2)))
  align.name <- unlist(lapply(als,function(i)
    names(i)[1]))
  align.name2 <- unlist(lapply(als,function(i)
    i[1]))

  return(cbind(align.name,align.name2,align.ec))
}
#'Select alignments
#'
#'Select a list of aligments that recover all the proteins, based on the
#'list of scores
#'@param localAligns a list of local alignments
#'@param scores a table of scores
#'@return selectAligns a list of local alignments
select.aligns <- function(localAligns, scores) {
  als <- localAligns
  als2 <- list()
  protsin <- 0
  prots <- c()
  n <- length(unique(scores[,1]))
  while (protsin < n) {
    i <- which.max(scores[,2])
    al1 <- als[[i]]
    prots <- append(prots,names(al1))
    prots <- unique(prots)
    protsin <- length(prots)
    for (p1 in names(al1)) {
      if (p1 %in% scores[,1]) {
        scores[which(scores[,1] == p1),2] <- - 1
      }
    }
    als2 <- append(als2,list(al1))
  }
  return(als2)
}

#'is.element2
#'@keywords internal
is.element2 <- function(i,j) {
  return(is.element(j,i))
}

#'Compute score
#'@keywords internal
compute.score <- function(als2,Sim) {
  blasts <- sim.score.all(als2,Sim)
  tamanys <- size.score.all(als2)
  score <- tamanys
  score[,3] <- as.numeric(score[,3]) / max(as.numeric(score[,3]))
  if (max(as.numeric(blasts[,3])) > 0) {
    score[,3] <-
      as.numeric(score[,3]) + as.numeric(blasts[,3]) /
      max(as.numeric(blasts[,3]))
  }
  colnames(score) <- c("V1","V2","V3")
  rownames(score) <- NULL
  sc3 <- floor(100 * as.numeric(score[,3]))
  score <- as.data.frame(score)
  score$V3 <- sc3
  return(score)
}


#'Global Alignment of two Protein Interaction Network
#'
#'Return a global alignment, from the list of local alignments and the table
#'of scores. The function first calculate a list of alignments with
#'\code{selectAligns}, then found a solution of the hypergraph mathching
#'problem. And finally extend this alignment to a global alignment
#'@param localAligns a list of local alignments
#'@param Sim a similarity matrix
#'@param AllSteps a boolean to determine whether return all intermediate
#'global alignments
#'@param getGlobal a boolean to set if the return is a global alignment or a better aligment but not strictly global
#'@return Global a global alignment and, if AllSteps is TRUE all intermediate alignments
align.global <- function(localAligns,Sim, AllSteps = TRUE, getGlobal = FALSE ) {
  global <- c()
  als <-
    unlist(unlist(localAligns,recursive = FALSE),recursive = FALSE)
  scores <- compute.score(als, Sim)
  Mat <-
    with(scores, sparseMatrix(
      i = as.numeric(V1), j = as.numeric(V2),
      x = V3, dimnames = list(levels(V1), levels(V2))
    ))
  globals <- list()
  Mat <- as.matrix(Mat)
  while (max(Mat) > 0) {
    hg <- HungarianFinal(Mat)
    als2 <- get.aligns(als, hg, Mat)
    score <- compute.score(als2, Sim)
    als2 <- select.aligns(als2, score[, c(1, 3)])
    hy <- hypergraph.solve(als2)
    global <- c(global, hy)
    globals <- append(globals, list(global))
    als <- aligns.update(als, global)
    if (length(als) > 0) {
      als <- remove.global(als, global)
      scores <- compute.score(als, Sim)
    }
    else {
      scores <- matrix(1:3, nrow = 1, ncol = 3)[ - 1,]
    }
    if (dim(scores)[1] > 0) {
      Mat <- with(scores, sparseMatrix(
        i = as.numeric(V1),
        j = as.numeric(V2),
        x = V3,
        dimnames = list(levels(V1), levels(V2))
      ))
      Mat <- as.matrix(Mat)
    }
    else {
      Mat <- matrix(0, nrow = 1, ncol = 1)
    }
  }
  if (getGlobal) {
  global2 <- align.end(localAligns, global)
  } else {
    global2 = global
  }
  if (AllSteps){
  return(list(globals, global2))
  }
  else{
    return(global2)
  }
}

#'Update aligns
#'@keywords internal
aligns.update <- function(als,global) {
  if (length(global) == 0) {
    return(als)
  }
  prots1 <- names(global)
  prots2 <- global
  als2 <- list()
  for (al in als) {
    if (!al[1] %in% prots2) {
      if (!names(al)[1] %in% prots1) {
        als2 <- append(als2,list(al))
      }
    }
  }
  return(als2)
}

#'get aligns
#'@keywords internal
get.aligns <- function(als,hg,Mat) {
  als2 <- list()
  prots1 <- hg[,1]
  prots2 <- hg[,2]
  for (al in als) {
    if (al[1] %in% prots2) {
      i <- which(al[1] == prots2)
      if (names(al)[1] == prots1[i]) {
        if (Mat[names(al[1]),al[1]] > 0) {
          als2 <- append(als2,list(al))
        }
      }
    }
  }
  return(als2)
}

#'remove global
#'@keywords internal
remove.global <- function(als2,global) {
  for (i in 1:length(als2)) {
    als2[[i]] <- als2[[i]][which(!als2[[i]] %in% global)]
    als2[[i]] <-
      als2[[i]][which(!names(als2[[i]]) %in% names(global))]
  }
  return(als2)
}

#'Solve hypergraph
#'@keywords internal
hypergraph.solve <- function(als) {
  E2 <- als
  E1 <- lapply(E2,names)
  scores <- unlist(lapply(E1, length))
  cprots1 <- count(unlist(E1))
  prots1 <- cprots1[,1]
  cprots2 <- count(unlist(E2))
  prots2 <- cprots2[,1]
  vars <- length(E1)
  constr1 <- length(prots1)
  constr2 <- length(prots2)
  constr <- constr1 + constr2
  lprec <- make.lp(constr,vars)
  constr11 <-
    lapply(prots1, function(i)
      as.numeric(unlist(lapply(E1,is.element2,i))))
  constr22 <-
    lapply(prots2, function(i)
      as.numeric(unlist(lapply(E2,is.element2,i))))

  ##Set constraints
  fcon1 <- function(i) {
    s <- sum(constr11[[i]])
    if (s > 0) {
      set.row(lprec,i,xt = rep(1,s),indices = which(constr11[[i]] == 1))
    }
  }

  aa <- lapply(1:constr1, fcon1)
  fcon2 <- function(i) {
    s <- sum(constr22[[i]])
    if (s > 0) {
      set.row(lprec,i + constr1,xt = rep(1,s),
              indices = which(constr22[[i]] == 1))
    }
  }

  bb <- lapply(1:constr2, fcon2)

  set.objfn(lprec,unlist(scores))
  set.constr.type(lprec,c(rep("<=",constr)))
  set.rhs(lprec,c(rep(1,constr)))
  set.bounds(
    lprec,lower = rep(0,vars),upper = rep(1,vars),columns = 1:vars
  )

  set.type(lprec,1:vars,"binary")
  break.value <- min(length(prots1),length(prots2))
  lp.control(
    lprec,sense = "max",verbose = "neutral",break.at.first = TRUE
  )
  solve(lprec)
  sols <- which(get.variables(lprec) > 0)

  getprots <- function(i,E2) {
    matrix(c(names(E2[[i]]),E2[[i]]),ncol = 2)
  }
  mmm <- getprots(sols[[1]],E2)
  nsol <- length(sols)
  if (nsol > 1) {
    for (i in 2:nsol) {
      mmm <- rbind(mmm,getprots(sols[[i]],E2))
    }
  }

  global <- mmm[,2]
  names(global) <- mmm[,1]
  return(global)
}


#'Update Matrix
#'@keywords internal
matrix.update <- function(mat,global) {
  mat2 <- mat
  rows.delete <- which(rownames(mat2) %in% names(global))
  cols.delete <- which(colnames(mat2) %in% global)
  mat2 <- mat2[ - rows.delete,]
  if (is.null(dim(mat2))) {
    if (length(mat2) > 0) {
      mat2 <- matrix(mat2,nrow = 1)
      rownames(mat2) <- setdiff(rownames(mat),names(global))
      colnames(mat2) <- colnames(mat)
      mat3 <- mat2[, - cols.delete]
      mat3 <- matrix(mat3,nrow = 1)
      colnames(mat3) <- setdiff(colnames(mat2),global)
      rownames(mat3) <- rownames(mat2)
      return(mat3)

    }
    else{
      return(matrix(0, nrow = 1,ncol = 1))
    }
  }
  mat3 <- mat2[, - cols.delete]
  if (is.null(dim(mat3))) {
    if (length(mat3) > 0) {
      mat3 <- matrix(mat3,ncol = 1)
      colnames(mat3) <- setdiff(colnames(mat2),global)
      rownames(mat3) <- rownames(mat2)
    }
    else{
      mat3 <- matrix(0, nrow = 1,ncol = 1)
    }
  }
  return(mat3)
}

#'End alignment
#'@keywords internal
align.end <- function(localAligns,global) {
  als <-
    unlist(unlist(localAligns,recursive = FALSE),recursive = FALSE)
  mmm <- cbind(names(global),global)
  E2 <- lapply(seq(1,length(als),2), function(i)
    als[i][[1]])
  E1 <- lapply(E2,names)

  prots1 <- unique(unlist(E1))

  rest2 <- unlist(E2)
  rest1 <- names(rest2)
  restm <- count(matrix(c(rest1,rest2),byrow = FALSE,ncol = 2))
  restm <- data.frame(restm,row.names = 1:dim(restm)[1])
  colnames(restm) <- c("V1","V2","V3")
  mat2 <- with(restm, sparseMatrix(
    i = as.numeric(V1),
    j = as.numeric(V2),
    x = V3,
    dimnames = list(levels(V1), levels(V2))
  ))
  mat2 = as.matrix(mat2)
  mat2[mat2 > 0] <- 1

  mat2[intersect(mmm[,1],rownames(mat2)),] <- -1

  mat2[,intersect(mmm[,2],colnames(mat2))] <- -1
  mat2 <- as.matrix(mat2) + 1
  hg <- HungarianFinal(as.matrix(mat2))


  for (i in 1:dim(hg)[1]) {
    if (mat2[hg[i,1],hg[i,2]] > 0) {
      mmm <- rbind(mmm,hg[i,])
    }
  }

  global2 <- mmm[,2]
  names(global2) <- mmm[,1]
  return(global2)
}

#'Plot the alignment
#'@param net1 an igraph object
#'@param net2 an igrpah object
#'@param global an alignment
#'@param k1 the width of the new edges of the alignment
#'@param k2 the width of the old edges
#'@param edge.curved Specifies whether to draw curved
#'edges, or not. This can be a logical or a numeric vector or scalar.
#'@param ... further arguments to be passed to igraph plotting
align.plot <-
  function(net1,net2,global,k1=1, k2=1, edge.curved = 0.5, ...) {
    coms1 <- fastgreedy.community(net1)
    coms2 <- fastgreedy.community(net2)
    newedges <- cbind(names(global),global)
    net3 <- graph.data.frame(newedges,directed = FALSE)
    net4 <- graph.union(net1,net2)
    net5 <- graph.union(net3,net4)
    num.eds <- ecount(net5)
    eds <- get.edgelist(net5)
    E(net5)$weight <- k1
    E(net5)$edge.curved <- 0
    eds1 <- get.edgelist(net4)
    ids <- get.edge.ids(net5,t(eds1),FALSE)

    net5 <- set.edge.attribute(net5,name = "weight",index = ids,k2)
    net5 <- set.edge.attribute(net5,name = "edge.curved",
                               index = ids,edge.curved)

    inds1 <- sort(coms1$membership,index.return = TRUE)

    lay1 <- 1:vcount(net1)
    lay2 <- rep(0,vcount(net2))
    lay12 <- lay1[inds1$ix]
    names(lay12) <- V(net1)$name
    lay <- cbind(10,rep(0,vcount(net5)))
    for (p1 in V(net1)$name) {
      i1 <- which(V(net5)$name == p1)
      p2 <- global[p1]
      i2 <- which(V(net5)$name == p2)
      lay[i1,] <- c(0,lay12[p1])
      lay[i2,] <- c(10,lay12[p1])
    }
    print("plot")
    print(net5)
    plot(
      net5,layout = layout.norm(lay),rescale = TRUE,
      edge.curved = E(net5)$edge.curved,edge.width = E(net5)$weight,...
    )
  }

#'Search alignments
#' Search alignments in als from p1 to p2
#' @param als a list of local alignments
#' @param p1 a protein
#' @param p2 a protein
#' @return list of alignments that includes p1 and p2
search.aligns <- function(als,p1,p2) {
  clustp1 <- c()
  clustp2 <- c()
  for (al in als) {
    if (p1 == names(al)[1]) {
      if (al[p1] == p2) {
        clustp1 <- names(al)
        clustp2 <- als
      }
    }
  }
  als2 <- list()
  for (al in als) {
    al.aux <- al[intersect(names(al),clustp1)]
    if (length(al.aux) > 0) {
      als2 <- append(als2,list(al.aux))
    }
  }
  return(als2)

}

#'Local alignment plot
#'plot a local alignment and all the local alignments that
#'intersect with him
#'@param localAligns a list of local alignments
#'@param global an alignment
#'@param p1 the center of cluster1
#'@param p2 the center of cluster2
#'@param net1 the first network
#'@param net2 the second network
#'@param ... further arguments to be passed to igraph plotting
align.local.plot <- function(localAligns,global,p1,p2,net1,net2,...) {
  als <- unlist(unlist(localAligns,recursive = FALSE),recursive = FALSE)
  als <- search.aligns(als,p1,p2)
  alp1p2 <- NULL
  for (al in als) {
    if (p1 == names(al)[1]) {
      if (al[p1] == p2) {
        alp1p2 <- al
      }
    }
  }
  if (is.null(alp1p2)) {
    print("Local alignment not found")
    return(NULL)
  }
  alini <- function(x,y) {
    matrix(c(as.character(alp1p2[x])
             ,as.character(alp1p2[y])),nrow = 1)
  }
  E1 <- get.edgelist(net1)
  if (dim(E1)[1] == 0) {
    return(0)
  }
  eds <- t(mapply(alini,E1[,1],E1[,2]))
  nas <- unique(c(which(is.na(eds[,2])),which(is.na(eds[,1]))))
  if (length(nas) > 0) {
    eds <- eds[ - nas,]
    eds1 <- E1[ - nas,]
  }
  if (is.null(dim(eds))) {
    if (length(eds) == 2) {
      eds <- cbind(eds[1],eds[2])
      eds1 <- cbind(eds1[1],eds1[2])
    }
  }
  eds <- rbind(eds,eds1)
  Gnet3 <- graph.edgelist(eds,directed = FALSE)

  cols <- rainbow(length(als))
  new.edges <- cbind(names(alp1p2),alp1p2)
  net3 <- graph.data.frame(new.edges,directed = FALSE)

  nets <- list(net3)
  for (al in als) {
    new.edges <- cbind(names(al),al)
    net3 <- graph.data.frame(new.edges,directed = FALSE)
    nets <- append(nets,list(net3))
  }
  net12 <- induced.subgraph(net1,vids = unique(unlist(lapply(als,names))))
  net22 <- induced.subgraph(net2,vids = unique(unlist(als)))
  net4 <- graph.union(net12,net22)
  net5 <- net4
  for (i in 1:length(nets)) {
    net5 <- graph.union(net5,nets[[i]])
  }
  num.eds <- ecount(net5)
  eds <- get.edgelist(net5)
  newedges <- cbind(names(global),global)
  net3 <- graph.data.frame(newedges,directed = FALSE)
  net3 <- induced.subgraph(net3, vids = intersect(V(net3)$name,V(net5)$name))
  E(net5)$color <- "black"
  E(net5)$lty <- 3
  E(net5)$width <- 0.5
  E(net5)$edge.curved <- 1
  ids <- c()
  eds1 <- get.edgelist(nets[[1]])
  ids2 <- get.edge.ids(net5,t(eds1),FALSE)
  ids2 <- setdiff(ids2,ids)
  net5 <- set.edge.attribute(net5,name = "color",index = ids2,cols[1])
  net5 <- set.edge.attribute(net5,name = "edge.curved",index = ids2,0)
  net5 <- set.edge.attribute(net5,name = "lty",index = ids2,2)

  ids <- c(ids,ids2)

  for (i in 2:length(nets)) {
    eds1 <- get.edgelist(nets[[i]])
    ids2 <- get.edge.ids(net5,t(eds1),FALSE)
    ids2 <- setdiff(ids2,ids)
    net5 <- set.edge.attribute(net5,name = "color",index = ids2,cols[i])
    net5 <- set.edge.attribute(net5,name = "edge.curved",index = ids2,0)
    ids <- c(ids,ids2)
  }
  eds1 <- get.edgelist(net3)
  ids <- get.edge.ids(net5,t(eds1),FALSE)

  net5 <- set.edge.attribute(net5,name = "lty",index = ids,1)
  net5 <- set.edge.attribute(net5,name = "edge.curved",index = ids,0)
  net5 <- set.edge.attribute(net5,name = "width",index = ids,1)

  eds1 <- get.edgelist(Gnet3)
  ids <- get.edge.ids(net5,t(eds1),FALSE)
  net5 <- set.edge.attribute(net5,name = "lty",index = ids,1)
  net5 <- set.edge.attribute(net5,name = "width",index = ids,1)
  net5 <- set.edge.attribute(net5,name = "color",index = ids,cols[1])

  eds1 <- get.edgelist(net22)
  ids <- get.edge.ids(net5,t(eds1),FALSE)
  net5 <- set.edge.attribute(net5,name = "edge.curved",index = ids, - 1)
  net5 <- induced.subgraph(net5,vids = union(alp1p2,names(alp1p2)))
  lay <- cbind(10,rep(0,vcount(net5)))
  cont1 <- 1
  cont2 <- 1
  for (i in 1:vcount(net5)) {
    p <- V(net5)$name[i]
    if (is.element(p,V(net12)$name)) {
      lay[i,] <- c(0,cont1)
      cont1 <- cont1 + 1
    }
    else{
      lay[i,] <- c(10,cont2)
      cont2 <- cont2 + 1

    }
  }
  e1 <- ecount(net5)
  v1 <- vcount(net5)
  print(paste(
    "plotting a graph with ",as.character(v1),
    " vertices and ",e1," edges",sep = ""
  ))
  net5 <- set.vertex.attribute(net5,"pos",
                               index = 1:vcount(net12),value = 1)
  net5 <- set.vertex.attribute(net5,"pos",
                               index = (vcount(net12) + 1):(vcount(net22) +
                                                              vcount(net12)),
                               value = - 1)

  plot(
    net5,layout = layout.norm(lay),rescale = TRUE,
    vertex.label.dist = V(net5)$pos,
    vertex.label.degree = pi,edge.curved = E(net5)$edge.curved,...
  )
}
