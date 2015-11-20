
#'Read Matrix
#'
#'Read a matrix in several formats
#'@param fileName directory of the matrix or data.frame (for mode col3)
#'@param mode table, col3, RData
#'@param sym True if the matrix must be symmetric.In modes table and RData,
#'if the matrix is not symmetric the value of M[i,j] will be (M[i,j]+M[j,i])/2.
#'In mode in 3col the value of M[i,j] will be the last that appears in the file
#'@param def the default value in mode 3col
#'@return matrix
read.matrix <- function(fileName,mode="table",sym="TRUE",def=0){
  switch(mode,
         table=return(read.matrix.table(fileName,sym)),
         col3=return(read.matrix.col3(fileName,def=def)),
         RData=return(read.matrix.RData(fileName,sym)),
         stop("Enter a valid mode"))
}

#'Intern
#'@keywords internal
read.matrix.table=function(fileName,sym="FALSE"){
  mat = as.matrix(read.table(fileName))
  if(sym){
    #Comprovar si es simetrica primer??????
         return((mat+t(mat))/2)
  }
  return(mat)
}

#'Intern
#'@keywords internal
read.matrix.RData=function(fileName,sym="FALSE"){
  load(fileName)
  mat=get(ls()[ls() != "fileName"])
  if(sym){
    #Comprovar si es simetrica primer??????
    return((mat+t(mat))/2)
  }
  return(mat)
}

#'Intern
#'@keywords internal
read.matrix.col3 = function(fileName,def=0){
  if(class(fileName)=="character"){
    MMM=read.table(fileName)
       }
       else{
         MMM=fileName
       }
  colnames(MMM)=c("V1","V2","V3")
  mat2=with(MMM, sparseMatrix(i = as.numeric(V1), j=as.numeric(V2), x=V3,dimnames=list(levels(V1), levels(V2))))
  return(mat2)
}


#'Read Network
#'
#'Read network from a file
#'@param fileName directory of network or a data.frame
#'@param mode tab : a file from DIP or BIOGRID of type tab or mitab;
#'edges: a list of edges with attributes (not necessary)
#'@param db database to select proteins of mode tab (uniprot,DIP,refseq,BIOGRID, gene)
#'@param cols which columns must be selected with mode edges, you can use "all" to
#'select all columns
#'@param sep separator between columns in edge mode
#'@param int.type boolean to decide if must be select the interaction type in tab mode.
read.network = function(fileName,mode="tab",db="uniprot",cols = "all",sep=" ",int.type=TRUE){
  switch(mode,
         tab=simplify(return(read.network.tab(fileName,db=db,int.type=int.type))),
         edges = simplify(return(read.network.edges(fileName,cols,sep=sep))),
         stop("No implemented yet"))
}

#' Network tab
#' @keywords internal
read.network.tab = function(fileName, db="uniprot", int.type=TRUE){
  if(class(fileName)=="character"){
  tab =read.table(fileName,sep="\t",row.names=NULL,stringsAsFactors=FALSE)
  }
  else{
    tab=fileName
  }
  tab2 = tab[,c(1,2)]
  get.int.db = function(i){
    switch(db,
           refseq = {x=splitrefseq(tab2[i,1])},
           uniprot = {x=splituniprot(tab2[i,1])},
           DIP = {x=splitDIP(tab2[i,1])},
           BIOGRID = {x=splitBIOGRID(tab2[i,1])},
           gene = {x=splitgene(tab2[i,1])},
           stop("Enter valid database"))
    switch(db,
           refseq = {y=splitrefseq(tab2[i,2])},
           uniprot = {y=splituniprot(tab2[i,2])},
           DIP = {y=splitDIP(tab2[i,2])},
           BIOGRID = {y=splitBIOGRID(tab2[i,2])},
           gene = {y=splitgene(tab2[i,2])},
           stop("Enter valid database"))

    return(list(x,y))
  }
  tab3=list()
  type=c()
  for(i in 1:dim(tab2)[1]){
    x=get.int.db(i)
    if(!is.na(x[[1]]) && !is.na(x[[2]])){
      tab3[[i]]=x[[1]]
      tab3[[i+dim(tab2)[1]]]=x[[2]]
      type[i]=tab[i,12]
    }
    else{
      tab3[[i]]=NA
      tab3[[i+dim(tab2)[1]]]=NA

    }
  }

  tab3=matrix(unlist(tab3)[!is.na(tab3)],ncol=2,byrow=FALSE)

  Net = graph.data.frame(tab3,directed=FALSE)
  if(int.type){
    type=type[!is.na(type)]
    Net=  set.edge.attribute(Net,name="Interaction type",value=type)
  }
  return(Net)
}

#' split refseq
#' @keywords internal
splitrefseq=function(x){
  x1 = strsplit(x,split="refseq:",fixed=TRUE)[[1]]
  if(length(x1)==2){
    return(strsplit(x1[2],split="|",fixed=TRUE)[[1]][1])
  }
  return(NA)
}

#' split uniprot
#' @keywords internal
splituniprot=function(x){
  x1 = strsplit(x,split="uniprotkb:",fixed=TRUE)[[1]]
  if(length(x1)==2){
    return(x1[2])
  }
  return(NA)
}

#' split DIP
#' @keywords internal
splitDIP=function(x){
  x1 = strsplit(x,split="DIP-",fixed=TRUE)[[1]]
  if(length(x1)==2){
      return(paste("DIP-",strsplit(x1[2],split="|",fixed=TRUE)[[1]][1],sep=""))

  }
  return(NA)
}

#' split BIOGRID
#' @keywords internal
splitBIOGRID = function(x){
  x1 = strsplit(x,split="BIOGRID:",fixed=TRUE)[[1]]
  if(length(x1)==2){
    return(x1[2])
  }
  return(NA)
}

#' split gene
#' @keywords internal
splitgene = function(x){
  x1 = strsplit(x,split="locuslink:",fixed=TRUE)[[1]]
  if(length(x1)==2){
   return(strsplit(x1[2],split="|",fixed=TRUE)[[1]][1])

  }
  return(NA)
}

#' Network edges
#' @keywords internal
read.network.edges = function(fileName,cols,sep){
  if(class(fileName)=="character"){
    tab = read.table(fileName,sep=sep,header = TRUE,stringsAsFactors=FALSE)

  }
  else{
    tab=fileName
  }
  if(cols[1] =="all"){
  return(graph.data.frame(tab,directed=FALSE))
  }
  return(graph.data.frame(tab[,cols],directed=FALSE))
}


#' Compute Matrix
#'
#' Compute several matrices
#' @param net1 an igraph object
#' @param net2 an igraph object for modes \code{BLAST} and \code{FC} if Net2 is not NULL,
#' computes the matrix between nodes in Net1 and Net2.
#' @param type what matrice do you want to compute (BLAST, DSD, Distance, FC, Degree)
#' @param mode for type \code{BLAST} : pident or bitscore.
#' For \code{DSD} and \code{Distance} distance, similarity or similarity by components.
#' For \code{FC} the list of ontologies
#' @param byComp for \code{DSD} and \code{Distance} if the similarity or distance must be normalized
#' by components
#' @param database for \code{BLAST} and \code{FC} the database which proteins belongs.
#' @param database2 for \code{BLAST} and \code{FC} the database which proteins belongs.
#' @param normalized if the matrix will be normalized
#' @return The matrix
compute.matrix = function(net1,net2=NULL, type="Distance", mode = "Similarity", database = NULL,database2=NULL, byComp=TRUE, normalized=TRUE){
  switch(type,
         BLAST =return(compute.matrix.Blast(net1,net2,mode,database,database2,normalized)),
         DSD =return(compute.matrix.DSD(net1,mode,byComp,normalized)),
         Distance =return(compute.matrix.Distance(net1,mode,byComp,normalized)),
         FC =return(compute.matrix.FC(net1,net2,mode)),
         Degree = return(compute.matrix.Degree(net1,net2)),
          stop("Enter a valid type")
  )
}

#'DSD
#'@keywords internal
compute.matrix.DSD=function(net,mode="Similarity",byComp=TRUE,normalized=TRUE){
  path <- paste(system.file(package="PINTAS"),"DSDmain.py", sep="/")
  n=length(V(net))
  prots=V(net)$name
  DSD = matrix(Inf,nrow=n,ncol=n)
  dimnames(DSD)=list(prots,prots)
  cc=decompose.graph(net)
  for(net in cc){
    tmp=tempfile()
    tmp2=tempfile()
    write.table(get.edgelist(net),quote=FALSE,file=tmp,row.names=FALSE,col.names=FALSE)
    command = paste("python",path,"-m 1 -o",tmp2,tmp)
    response= system(command, intern=T)
    table = as.matrix(read.table(paste(tmp2,"DSD1",sep=".")))
    diam=max(table)+1
    if(byComp){
      if(mode=="Similarity"){
        DSD[rownames(table),colnames(table)]=(max(table)+1-table)/(max(table)+1)
      }
      else{
        DSD[rownames(table),colnames(table)]= table/max(table)
      }
    }
    else{
      DSD[rownames(table),colnames(table)]=table
    }
  }
  mmm = max(DSD[DSD<Inf])
  if(!byComp){
    if(mode=="Similarity"){
      DSD = (mmm+1-DSD)/(mmm+1)
    }
    else{
      if(normalized){
        DSD = DSD/mmm
      }
    }
  }
  if(mode=="Similarity"){
    DSD[DSD==Inf]=0
  }
  DSD[DSD==-Inf]=0
  return(DSD)

  }

#'Distance
#'@keywords internal
compute.matrix.Distance=function(net,mode="Similarity",byComp=TRUE,normalized=TRUE){
  n=length(V(net))
  prots=V(net)$name
  if(!byComp){
    dist = shortest.paths(net)
    mmm = max(dist[dist<Inf])
    if(mode=="Similarity"){
      dist = (mmm+1-dist)/(mmm+1)
      dist[dist==-Inf]=0
      return(dist)
    }
    if(normalized){
      dist= dist/mmm
    }
    return(dist)
  }
  if(mode=="Similarity"){
    dist=matrix(0,nrow=n,ncol=n)
  }
  else{
    dist=matrix(Inf,nrow=n,ncol=n)
  }
  dimnames(dist)=list(prots,prots)
  cc=decompose.graph(net)
  for(nnn in cc){
  dist2 = shortest.paths(nnn)
    mmm2 = max(dist2[dist2<Inf])
    if(mode=="Similarity"){
      dist2 = (mmm2+1-dist2)/(mmm2+1)
      dist2[dist2==-Inf]=dist2[dist2==Inf]=0
      dist[V(nnn)$name,V(nnn)$name]=dist2
    }
    if(normalized){
      dist2= dist2/max(dist2[dist2<Inf])
    }
    dist[V(nnn)$name,V(nnn)$name]=dist2
  }
  return(dist)
}

#'FC
#'@keywords internal
compute.matrix.FC=function(net1,net2=NULL,gos){
  onenet=FALSE
  if(is.null(net2)){
    net2=net1
    onenet=TRUE
  }
  prots1=V(net1)$name
  prots2=V(net2)$name
  if(onenet) FSim=diag(1, nrow = length(prots1),ncol=length(prots2))
  else{
    FSim = matrix(0, nrow = length(prots1),ncol=length(prots2))
  }
  dimnames(FSim)=list(prots1,prots2)
  if(onenet){
    for(i in 1:(length(prots1)-1)){
      for(j in (i+1):length(prots1)){
        fc = length(intersect(gos[[prots1[i]]],gos[[prots1[j]]]))/length(
          union(gos[[prots1[i]]],gos[[prots1[j]]]))
        FSim[i,j]=FSim[j,i]=fc
      }
    }

  }
  else{
    for(i in 1:length(prots1)){
      for(j in 1:length(prots2)){
        fc = length(intersect(gos[[prots1[i]]],gos[[prots2[j]]]))/length(
          union(gos[[prots1[i]]],gos[[prots2[j]]]))
        FSim[i,j]=fc
      }
    }
  }
  return(FSim)
}

#'Degree
#'@keywords internal
compute.matrix.Degree=function(net1,net2=NULL){
  if(is.null(net2)){
    net2=net1
 }
 deg = degree(net1)
 deg2 = degree(net2)
 return(matrix(unlist(lapply(deg, function(i) lapply(deg2, function(j) abs(i-j)))),nrow=length(deg), byrow=TRUE))

}

#'IC
#'@keywords internal
computeIC = function(net1,net2){
  degree1 = degree(net1)
  degree2 = degree(net2)
  neigh1 = neighborhood(net1,order=1)
  neigh2 = neighborhood(net2,order=1)
  lens1 = unlist(lapply(neigh1,length))
  lens2 = unlist(lapply(neigh2,length))
  sums1 = unlist(lapply(neigh1,function(i) sum(1/lens1[i])))
  sums2 = unlist(lapply(neigh2,function(i) sum(1/lens2[i])))
  mat = matrix(0,nrow=vcount(net1),ncol=vcount(net2))
  maxneig = max(degree1,degree2)
  for(i in 1:vcount(net1)){
    for(j in 1:vcount(net2)){
      mat[i,j]=min(sums1[i],sums2[i])/maxneig
    }
  }
  dimnames(mat)=list(V(net1)$name,V(net2)$name)
  return(mat)
}

#'No funciona
#'@keywords internal
compute.matrix.Blast=function(Net1,Net2,mode,database,database2,normalized){
  tmp =tempfile()
  buildBlast(Net1,Net2,mode,database,database2,tmp)
  BlastM = read.matrix.table(tmp)
  if(normalized){
    return(BlastM/max(BlastM))
  }
  return(BlastM)
}

#'No funciona
#'@keywords internal
buildBlast = function(Net1,Net2,mode,database,database2,tmp){
  path <- paste(system.file(package="PINTAS"),"Blast.jar", sep="/")
  tmp1 = tempfile()
  write(V(Net1)$name,file=tmp1)
  tmp2 = tempfile()
  if(Net2==NULL){
    tmp2=tmp1
  }
  else{
    write(V(Net2)$name,file=tmp2)
  }
  if(database2==NULL){
    database=database2
  }
  command = paste("java",path,"-prot1",tmp1,"-prot2",tmp2,"-db1",
                  database,"-db2",database2,"-m",mode,"-outfmt mat -o",tmp)
  response= system(command, intern=T)

}

#'Cluster Network
#'
#'Compute the cluster matrix from the similarity matrix sigma, where all the nodes have more
#'similarity than lambda and the size of cluster is less than k.
#'
#'@param sigma a similarity matrix
#'@param lambda the similarity threshold
#'@param k the size threshold
#'@return The cluster matrix, where M[i,j] is 1 if node i belongs to the cluster of node j and
#'0 otherwise
cluster.network <- function(sigma, lambda=0, k=dim(sigma)[1]){
 if(is.data.frame(sigma)){
   sigma=read.matrix.col3(sigma,def=0)

 }
  n=dim(sigma)[1]
  clustmatrix=matrix(as.numeric(c(as.matrix(sigma))>=lambda),nrow=n,byrow=TRUE)
  sums = apply(clustmatrix,2,sum)
  ind=which(sums>k)
  if(length(ind)>0){
  lapply(ind,function(i)  clustmatrix[sort(as.matrix(sigma)[i,],index.return=TRUE)$ix[1:(n-k)],i]<<-0)
  }
  dimnames(clustmatrix)= dimnames(sigma)
  clustmatrix=as(clustmatrix,"sparseMatrix")
  return(clustmatrix)

}

# #'Compute i,j
# #'@keywords internal
# compute.i.j = function(mat,i,j){
#   if(is.data.table(mat)){
#     return(mat[list(i,j)]$V3)
#   }
#   return(mat[i,j])
# }

#'Extract clusters
#'
#'Compute the subnetworks of \code{Net} from a cluster matrix
#'@param Net an igraph object
#'@param ClustMat a cluster matrix (output of \code{cluster.matrix})
#'@return Clusters in igraph format
extract.clusters <- function(Net, ClustMat){
  prots=rownames(ClustMat)
  clusts = lapply(1:dim(ClustMat)[2], function(i) induced.subgraph(Net,prots[ClustMat[,i]==1]))
  return(clusts)

}

#'Display Clusters
#'
#'Given a cluster matrix, see \code{cluster.network} and a network
#'display the cluster matrix with the following colors for the position (i,j):
#'
#'- Yellow if the protein and protein i doesn't belongs to cluster of protein j and
#' proteins don't interact in the network
#'
#'-  Black if the protein and protein i belongs to cluster of protein j but
#' proteins don't interact in the network
#'
#'- Red if the protein and protein i doesn't belongs to cluster of protein j but
#' proteins interact in the network
#'
#'-  Green if the protein and protein i belongs to cluster of protein j and
#' proteins interact in the network
#'
#'@param clust a matrix which is the output1 of \code{cluster.network}
#'@param cols a list of 4 colors if you want to change the default colors
#'@param zoom an integer to define the size of plot or NA, to plot all clusters
#'@param type 0 = "yellow", 1 = "black", 2="red", 3 = "green"
#'@param col color to use with zoom
#'@param Net an igraph object
#'@param ... Additional plotting parameters
display.clusters <- function(clust, Net,zoom=NA, type = 1, cols=c("yellow","black","red","green"),col=cols[type+1],...){
  clust = as.matrix(clust)
  if(is.na(zoom)){
  mmm =   clust+2*as.matrix(get.adjacency(Net))
  dimnames(mmm) = list(1:dim(mmm)[1],1:dim(mmm)[2])
 # par(oma=c(1,1,1,1))
  image2D(mmm,
          col=cols,
          axes=F,
          xlab=NA,
          ylab=NA,
          colkey=FALSE,
          ...)
  axis(1,at=seq(0,1,by=1/(dim(clust)[2]-1)),
       labels=colnames(clust),las=2,
       pos=-1/(2*(dim(clust)[1]-1)))
  axis(2,at=seq(0,1,by=1/(dim(clust)[1]-1)),
       labels=rownames(clust),las=2,pos=-1/(2*(dim(clust)[2]-1)))
  legend(x=0.2,y=1.4,legend=0:3,fill=cols,horiz=TRUE,bty="n",xpd=TRUE,)

}
else{
  mat = clust+2*as.matrix(get.adjacency(Net))
  k = floor(length(V(Net))/zoom)
  conts = unlist(lapply(1:zoom,function(i) lapply(1:zoom, function(j) cont.type(mat,(i-1)*k+1,i*k,(j-1)*k+1,j*k,type)) ))
  mat2 = matrix(conts,nrow=zoom, byrow=TRUE)
  print("Min")
  print(min(mat2))
  print("Max")
  print(max(mat2))
  image2D(mat2, col = ramp.col(col=c("white",col),n=k*k),
          axes=F,
          xlab=NA,
          ylab=NA,
          ...)

}

}

#' Cont type
#' @keywords internal
cont.type = function(mat,init1, fin1,init2,fin2,type){
  return(length(which(mat[init1:fin1,init2:fin2]==type)))
}
