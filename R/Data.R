#'@name net1
#'@aliases net1
#'@title Protein Interaction Network
#'@description A small protein interaction network. Used as example
#'@docType data
#'@author Adria Alcala, 06/05/2015
NULL

#'@name net2
#'@aliases net2
#'@title Protein Interaction Network
#'@description A small protein interaction network. Used as example
#'@docType data
#'@author Adria Alcala, 06/05/2015
NULL

#'@name Sim1
#'@aliases Sim1
#'@title Example of similarity network
#'@description A similiarity matrix between proteins in net1
#'@format a \code{matrix} instance
#'@author Adria Alcala, 28/04/2015
NULL

#'@name Sim2
#'@aliases Sim2
#'@title Example of similarity network
#'@description A similarity matrix between proteins in net2
#'@format a \code{matrix} instance
#'@author Adria Alcala, 28/04/2015
NULL

#'@name Sim
#'@aliases Sim
#'@title Example of similarity network
#'@description A similarity matrix between proteins in net1 and net2
#'@format a \code{matrix} instance
#'@author Adria Alcala, 28/04/2015
NULL

#'@name dme
#'@aliases dme
#'@title Protein protein interaction network of dme
#'@description The biigest connected component of the protein protein interaction
#'network of drosophila melanogaster. Extracted from IsoBase
#'@format a \code{igraph} instance
#'@docType data
#'@author Adria Alcala, 07/12/2015
NULL

#'@name sce
#'@aliases sce
#'@title Protein protein interaction network of sce
#'@description The biigest connected component of the protein protein interaction
#'network of saccharomyces cerevisiae. Extracted from IsoBase
#'@format a \code{igraph} instance
#'@docType data
#'@author Adria Alcala, 07/12/2015
NULL

#'@name dme-dme
#'@aliases dme.dme
#'@title Similarity matrix between proteins of dme
#'@description A similarity matrix between the proteins of dme. The similarity was calculated
#'as an average between the blast bitscore normalized of the sequence of the proteins and the
#'similarity of the degrees. The second part was calculated using the function
#'\code{compute.matrix}
#'@format a \code{matrix} instance
#'@docType data
#'@author Adria Alcala, 07/12/2015
NULL

#'@name sce-sce
#'@aliases sce.sce
#'@title Similarity matrix between proteins of sce
#'@description A similarity matrix between the proteins of sce. The similarity was calculated
#'as an average between the blast bitscore normalized of the sequence of the proteins and the
#'similarity of the degrees. The second part was calculated using the function
#'\code{compute.matrix}
#'@format a \code{matrix} instance
#'@docType data
#'@author Adria Alcala, 07/12/2015
NULL

#'@name dme-sce
#'@aliases dme.sce
#'@title Similarity matrix between proteins of dme and sce
#'@description A similarity matrix between the proteins of dme. The similarity is
#' the blast bitscore normalized of the sequence of the proteins.
#'@format a \code{matrix} instance
#'@docType data
#'@author Adria Alcala, 07/12/2015
NULL

#'@name go
#'@aliases go
#'@title A list of go
#'@description A list of the gene ontology of the proteins in dme and sce protein
#'protein interaction network
#'@format a \code{list} instance
#'@source IsoBase
#'@author Adria Alcala, 07/12/2015
NULL

#'@name internal
#'@title A list of alignments for the vignette
#'@description A list of alignments used in the vignette dme-sce
#'@keywords internal
NULL

