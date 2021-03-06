% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/trimesh.R
\name{trimesh}
\alias{trimesh}
\title{Display triangles mesh (2D)}
\usage{
trimesh(T, p, p2, add = FALSE, axis = FALSE, boxed = FALSE, ...)
}
\arguments{
\item{T}{T is a \code{m}-by-3 matrix. A row of \code{T} contains
indices into \code{X} of the vertices of a triangle. \code{T} is
usually the output of \code{\link{delaunayn}}.}

\item{p}{A vector or a matrix.}

\item{p2}{if \code{p} is not a matrix \code{p} and \code{p2} are bind to a
matrix with \code{cbind}.}

\item{add}{Add to existing plot in current active device?}

\item{axis}{Draw axes?}

\item{boxed}{Plot box?}

\item{\dots}{Parameters to the rendering device. See the \link[rgl]{rgl}
package.}
}
\description{
\code{trimesh(T, p)} displays the triangles defined in the m-by-3
matrix \code{T} and points \code{p} as a mesh.  Each row of
\code{T} specifies a triangle by giving the 3 indices of its
points in \code{X}.
}
\examples{
#example trimesh
p = cbind(x=rnorm(30), y=rnorm(30))
tt = delaunayn(p)
trimesh(tt,p)
}
\author{
Raoul Grasman
}
\seealso{
\code{\link{tetramesh}}, \code{\link[rgl]{rgl}},
\code{\link{delaunayn}}, \code{\link{convhulln}},
\code{\link{surf.tri}}
}
\keyword{hplot}

