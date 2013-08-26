\name{make.3D.Sphere}
\alias{make.3D.Sphere}
\title{
Plots a sphere centered at the specified amino acid with a specified radius.
}
\description{
Plots a sphere of radius r centered at a specific residue. Currently only 1 sphere can be plotted. The rgl pckage is required.
}
\usage{
make.3D.Sphere(position.matrix, center, radius, alpha = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{position.matrix}{
A dataframe consisting of six columns: 1) Residue Name, 2) Amino Acid number in the protein, 3) Side Chain, 4) X-coordinate, 5) Y-coordinate and 6) Z-coordinate.   Please see \code{\link[iPAC]{get.Positions}} and \code{\link[iPAC]{get.AlignedPositions}} in the \pkg{iPAC} package.
}

\item{center}{
  The residue number you want the sphere centered at. Use the number from the "Can.Count" column in the position matrix.
}
  \item{radius}{
The radius of the sphere you would like to draw.
}
  \item{alpha}{
The "darkness" of the sphere. 
}
}

\value{
The rgl package is called and a graph is shown.
}

\note{
This function is made for ease of use in preliminary analysis only. For more sophisticated graphing packages, consider using PyMol at:
http://www.pymol.org
}

\references{
Daniel Adler and Duncan Murdoch (2013). rgl: 3D visualization device system (OpenGL). R package version 0.93.935. http://CRAN.R-project.org/package=rgl


}
\seealso{
\code{\link[iPAC]{get.Positions}}
}


\examples{
\dontrun{
#Plots a sphere centered around amino acid 12 with radius 3.
library(rgl)

#loads the data
CIF <- "http://www.pdb.org/pdb/files/3GFT.cif"
Fasta <- "http://www.uniprot.org/uniprot/P01116-2.fasta"
KRAS.Positions <- get.Positions(CIF, Fasta, "A")

#generates the plot
make.3D.Sphere(KRAS.Positions$Positions, 12, 3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
