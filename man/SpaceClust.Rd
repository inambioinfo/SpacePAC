\name{SpaceClust}
\alias{SpaceClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SpaceClust
}
\description{
Finds mutational clusters via simulation. There are two options currently avaiable. The first is "SimMax" and the second is "Poisson".  The Poisson method is faster and finds the 1 sphere with the largest number of mutations at each radius. A bonferroni adjustment is then used to account for multiple radii. The SimMax method uses the ``simMaxspheres" parameter to find the 1, 2 or 3 (non-overlapping) spheres that together have the most number of mutations. A simulation approach is then used to find the most significant clusters. Please see the vignette for further details.
}
\usage{
SpaceClust(mutation.data, position.matrix, method = "SimMax", numsims = 1000,
		 simMaxSpheres = 3, radii.vector, multcomp = "bonferroni", alpha = 0.05)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

\item{mutation.data}{
A matrix of 0's (no mutation) and 1's (mutation) where each column represents an amino acid in the protein and each row represents an individual sample (test subject, cell line, etc). Thus if row i in column j had a 1, that would mean that the jth amino acid for person i had a nonsynonomous mutation. Please note that getting the mutation matrix is the responsibility of the user. Further, the column names of the matrix must be in the format V1, V2, ..., Vn where n is the total number of residues. One source of information is the COSMIC database \url{http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/}. However, extraction from this (or any other database) is not trivial and often requires pre-processing work. In the case of COSMIC, a local SQL server must be set up and query of the user's design must be run to pull the correct mutational data. This data must then be manipulated by the user into the matrix described. Please note, that the mutational data should come from a whole gene or a whole genome study and can not be selectively chosen as that will violate the uniformity assumption that the algorithm is based on.
}

  \item{position.matrix}{
A dataframe consisting of six columns: 1) Residue Name, 2) Amino Acid number in the protein, 3) Side Chain, 4) X-coordinate, 5) Y-coordinate and 6) Z-coordinate. Please see \code{\link[iPAC]{get.Positions}} and \code{\link[iPAC]{get.AlignedPositions}} in the \pkg{iPAC} package.

}
  \item{method}{
Either "SimMax" or "Poisson". Please see the vignette for further details on the difference.
}
  \item{numsims}{
The number of times to simulate the mutations on the protein. For each simulation, the mutations are uniformly distributed on the protein.
}
  \item{simMaxSpheres}{
The maximum number of spheres to consider. Currently, the implementation allows for simMaxspheres to be either 1, 2 or 3.
}
  \item{radii.vector}{
  A vector of radii. For each radius, we find the best sphere combination. If the positional data is obtained from the pdb, the radii are measured in Angstroms as the x,y,z coordinates in the PDB are in Angstroms. Thus, for instance, a radius of ``5" means that all residues with their carbon-alpha atom within 5 Angstroms of the center are included in the sphere.
}
  \item{multcomp}{
If the Poisson method is used, a multiple comparison adjustment is required to account for the multiple sphere. As the sphere iterates through the protein (centered at each amino acid), a p-value is calculated for each sphere. Options are: ``Bonferroni", "BH", or "none". The "BH" method stands for the Benjamini-Hochberg FDR correction. Please see \code{\link[stats]{p.adjust}} for a full description.
}
  \item{alpha}{
If the Poisson method is used, alpha is used as the cutoff value after the appropriate multiple comparison adjustment.
}
}
\details{
For the SimMax method, no multiple comparison is required for different radii sizes and sphere positions. See the vignette for more information. Furthermore, note that on average, residues are 3 Angstroms apart.
}
\value{
If the method is Poisson, the result is a list with the following components:
\item{result.poisson}{A data frame of the most significant clusters. The data frame has the following columns for each cluster (clusters shown as rows): 
Center:  The amino acid at which the sphere is centered.
Start: The smallest numbered residue in the sphere.
End: The largest numbered residue in the sphere.
Positions: The mutated positions in the sphere.
MutsCount: The total number of mutations in the sphere.
P.Value: The p-value for the cluster.
Within.Range: The residues within the sphere.
Line.Length: End-Start.}

\item{best.radii}{The radii at which the lowest p-value for the most significant cluster was found. Only the matrix for this p-value is shown.}

If the method is SimMax, the result is a list with the following components.
\item{p.value}{The smallest p.value when considering 1,2 or 3 spheres. Will match the p-value for the optimal sphere configuration.}
\item{optimal.num.spheres}{The number of spheres with the most statistically significant p-value.}
\item{optimal.radius}{The radius at which the most statistically significant p-value is identified.}
\item{optimal.sphere}{This presents the sphere results with the most statistically significant p-value. It will automatically display whether 1, 2 or 3 spheres is best. In the very unlikely event that more than one sphere result has the same z-score (for instance the z-score is the same whether you consider 2 or 3 spheres), the result that uses the minimum number of spheres will be displayed.}
\item{best.1.sphere}{This shows the orientation of the most statistically significant sphere. It will display the following items: 1) Center:  The amino acid at which the sphere is centered. 2) Start: The smallest numbered residue in the sphere. 3) End: The largest numbered residue in the sphere. 4) Positions: The mutated positions in the sphere. 5) MutsCount: The total number of mutations in the sphere. 6) Z-Score: The normalized z-score as defined in the vignette. 7) Within.Range: The residues within the sphere. 
8) Line.Length: End-Start.}
\item{best.2.sphere}{This shows the orientation of the most statistically significant 2 spheres. The entries are the same as the items for ``best.1.sphere" except for a "1" or "2" appended to each column name. A ``1" means that the information presented in the column belongs to the first sphere while a ``2" means that the information in the column belongs to the second sphere. The ``MutsCountTotal" column shows how many mutations are in both spheres and is just the sum of ``MutsCount1" and ``MutsCount2". Finally, the ``Intersection" column is the intersection of ``Within.Range1" and ``Within.Range2" and should be blank unless an error occurs.}
\item{best.3.sphere}{This shows the orientation of the most statistically significant 3 spheres. The entries are the same as in ``best2.sphere" except now there is a ``1", ``2" or ``3" appended to each column to signify whether the 1st, 2nd or 3rd sphere is being considered.}
\item{best.1.sphere.radius}{The radius that provides the most statistically significant result when only 1 sphere is considered.}
\item{best.2.sphere.radius}{The radius that provides the most statistically significant result when only 2 spheres are considered.}
\item{best.3.sphere.radius}{The radius that provides the most statistically signifificant result when only 3 spheres are considered.}
\item{bad.2.sphere.message}{If finding the optimal 2 spheres caused an error (possibly because no non-overlapping spheres or all the mutations are on one residue) a message is shown here with more details.}
\item{bad.3sphere.message}{If finding the optimal 3 spheres caused an error (possibly because no non-overlapping spheres or all the mutations are on one or two residues) a message is shown here with more details.}
\item{bad.2sphere.radii}{If finding the optimal 2 spheres caused an error, the radii at which errors occurred are displayed.}
\item{bad.3sphere.radii}{If finding the optimal 3 spheres caused an error, the radii at which errors occurred are displayed.}

}
\references{
Torsten Hothorn, Frank Bretz and Peter Westfall (2008). Simultaneous Inference in General Parametric Models. Biometrical Journal 50(3), 346--363.

}
\note{
See the 'multcomp' package on CRAN for a description of how the multiple comparison adjustment is made.

If you use the Poisson method, a Bonferroni correction is used to adjust for all the radii. As an example, supposing that the most significant cluster is found at radius 5, and the radii vector was (1,2,3,4,5), the p-values displayed in the result matrix would be the p-value_b*5 where p-value_b is the p-value if algorithm was run with radii vector= c(5).
}
\seealso{
\code{\link[iPAC]{get.Positions}}
\code{\link[SpacePAC]{SpaceClust}}
}


\examples{
CIF <- "http://www.pdb.org/pdb/files/3GFT.cif"
Fasta <- "http://www.uniprot.org/uniprot/P01116-2.fasta"
KRAS.Positions <- get.Positions(CIF, Fasta, "A")
data(KRAS.Mutations)


#Calculate the required clusters using SimMax
SpaceClust(KRAS.Mutations, KRAS.Positions$Positions, radii.vector = c(1,2,3,4))

#Calculate the required clusters using Poisson
SpaceClust(KRAS.Mutations, KRAS.Positions$Positions, radii.vector = c(1,2,3,4), method = "Poisson")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clusters }
\keyword{ Mutations }