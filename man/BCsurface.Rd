\name{BCsurface}
\alias{BCsurface}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for use of interpolation methods in the steam cross-section
%%  ~~function to do ... ~~
}
\description{It adds the information to the intersection of gauging verticals anf water surface. Value of velocity is assumed to be the same as in the top measured point in gauging vertical there.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BCsurface(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data must be in format: data.frame. In data there must be columns: Loc, Depth, MeasD, Vel. Means: Loc - "x" location of point from the bank of the stream [meters], Depth - depth of vertical of measuring point [meters], MeasD - depth of measuring point in vertical [meters], Vel - measured point velocity [m*s^(-1)]. Gauging verticals must be in order from the lowest value of "x" location to highest "x" location or from the highest value of "x" location to the lowest value of "x" location.
}
}
\details{
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
BCsurface(data_02)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
