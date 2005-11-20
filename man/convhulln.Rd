\name{convhulln}
\alias{convhulln}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute smallest convex hull that encloses a set of points }
\description{
  Returns an index matrix to the points of simplices (\dQuote{triangles}) that
  form the smallest convex simplicial complex of a set of input points in
  N-dimensional space.
}
\usage{
convhulln(p, options = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ An \code{n}-by-\code{dim} matrix. The rows of \code{p} represent \code{n}
            points in \code{dim}-dimensional space.  }
  \item{options}{ Optional options, see details below. }
}
\details{
This function interfaces the qhull library, and intends to be a port from
Octave to R.

The input \code{n}-by-\code{dim} matrix contains \code{n} points of dimension
\code{dim}. If a second optional argument is given, it must be a string containing
extra options for the underlying qhull command. (See the Qhull
documentation for the available options - refer to \link{delaunayn}.)
}
\value{
 An \code{m}-by-\code{dim} index matrix of which each row defines
 a \code{dim}-dimensional \dQuote{triangle}. The indices refer to the
 rows in \code{p}.
}
\references{
\cite{Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T., \dQuote{The
Quickhull algorithm for convex hulls,} \emph{ACM Trans. on Mathematical
Software,} Dec 1996.}

\url{http://www.qhull.org}
}

\author{ Raoul Grasman }
\note{ This intents to be a port of the Octave's (\url{http://www.octave.org}) geometry
library. The sources originals were from Kai Habel.

Currently a call to convhulln unavoidably generates a diagnostic report
on an Rterm console. This is to be removed in future release. See \link{surf.tri}
for a quiet alternative, which computes the hull from a Delaunay
triangulation of the points.

See further notes in \link{delaunayn}.
}

\seealso{ \code{\link[tripack]{convex.hull}},
          \code{\link{delaunayn}},
          \code{\link{surf.tri}},
          \code{\link{distmesh2d}} }

\examples{
# example delaunayn
d = c(-1,1)
pc = as.matrix(rbind(expand.grid(d,d,d),0))
tc = delaunayn(pc)

# example tetramesh
\dontrun{
library(rgl)
rgl.viewpoint(60)
rgl.light(120,60)
tetramesh(tc,pc, alpha=0.9)    # render tetrahedron mesh
}

# example convhulln
# ==> see also surf.tri to avoid unwanted messages printed to the console by qhull
ps = matrix(rnorm(3000),ncol=3)                     # generate poinst on a sphere
ps = sqrt(3) * ps / drop(sqrt((ps^2) \%*\% rep(1,3)))
ts.surf = t( convhulln(ps,"QJ") ) # see the qhull documentations for the options
\dontrun{
rgl.triangles(ps[ts.surf,1],ps[ts.surf,2],ps[ts.surf,3],col="blue",alpha=.2)
for(i in 1:(8*360)) rgl.viewpoint(i/8)
}
}
\keyword{ math }
\keyword{ dplot }
\keyword{ graphs }
