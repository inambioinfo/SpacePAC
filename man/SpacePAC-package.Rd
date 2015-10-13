\name{SpacePAC-package}
\alias{SpacePAC-package}
\alias{SpacePAC}
\docType{package}
\title{
Identifying mutational clusters in 3D protein space using simulation.
}
\description{
The \emph{SpacePAC} package identifies non-random amino acid clusters in proteins in 3D space and is a sister package to \emph{iPAC} and \emph{GraphPAC}. \emph{SpacePAC} considers 1, 2 or 3 non-overlapping spheres with radii specified by the user and through simulation, attempts to identify spheres where there are more mutations than expected by random chance alone. These results are then outputted in the form of a list with p-values.
}
\details{
Please see \code{\link[iPAC]{get.Positions}} and \code{\link[iPAC]{get.AlignedPositions}} in the \pkg{iPAC} package for information about obtaining positional data. 
}
\author{
Gregory Ryslik, Yuwei Cheng, Hongyu Zhao

Maintainer: Gregory Ryslik <gregory.ryslik@yale.edu>

}
\references{
Gregory Ryslik and Hongyu Zhao (2012). iPAC: Identification of Protein
Amino acid Clustering. R package version 1.1.3. \url{http://www.bioconductor.org/}.

Gregory Ryslik and Hongyu Zhao (2013). GraphPAC: Identification of Mutational Clusters in Proteins via a Graph Theoretical Approach. R Package version 1.0.0 \url{http://www.bioconductor.org/}.

Bioconductor: Open software development for computational biology and bioinformatics R. Gentleman, V. J. Carey, D. M. Bates, B.Bolstad, M. Dettling, S.
  Dudoit, B. Ellis, L. Gautier, Y. Ge, and others 2004, Genome Biology, Vol. 5, R80

}

\seealso{
\code{\link[iPAC]{get.Positions}}
\code{\link[SpacePAC]{SpaceClust}}
}

\examples{
\dontrun{
CIF <- "http://www.pdb.org/pdb/files/3GFT.cif"
Fasta <- "http://www.uniprot.org/uniprot/P01116-2.fasta"
KRAS.Positions <- get.Positions(CIF, Fasta, "A")
data(KRAS.Mutations)

#Calculate the required clusters
SpaceClust(KRAS.Mutations, KRAS.Positions$Positions, radii.vector = c(1,2,3,4))

}
}
