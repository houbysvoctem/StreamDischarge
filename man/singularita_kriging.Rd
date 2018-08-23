\name{singularita_kriging}
\alias{singularita_kriging}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for use of interpolation methods in the stream cross-section
%%  ~~function to do ... ~~
}
\description{ Solution of possible matrix singularity in case of ordinary kriging.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
singularita_kriging(m, mesh, known, i, mynugget, mysill, myrange)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{number of points used to estiomation}
  \item{mesh}{points where the estimation should be done}
  \item{known}{measured velocities from which the estimation will be done}
  \item{i}{point to be estimated}
  \item{mynugget}{parameter of kriging Gauss semivariogram}
  \item{mysill}{parameter of kriging Gauss semivariogram}
  \item{myrange}{parameter of kriging Gauss semivariogram}
%%     ~~Describe \code{x} here~~

}
\details{Uses more points for estimation to created non-singular matrix.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Eliška Friedlová (Hasníková)
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
