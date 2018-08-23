\name{zjisteniP}
\alias{zjisteniP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for use of interpolation methods in the steam cross-section
%%  ~~function to do ... ~~
}
\description{It checks, how the boundarz conditions are set.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
zjisteniP(data, p1, p2, p3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data must be in format: data.frame. In data there must be columns: Loc, Depth, MeasD, Vel. Means: Loc - "x" location of point from the bank of the stream [meters], Depth - depth of vertical of measuring point [meters], MeasD - depth of measuring point in vertical [meters], Vel - measured point velocity [m*s^(-1)]. Gauging verticals must be in order from the lowest value of "x" location to highest "x" location or from the highest value of "x" location to the lowest value of "x" location.
}
  \item{p1}{boundary condition p1: if p1=1, velocity in the intersection of the river stream bed and gauging vertical is 0.}
  \item{p2}{boundary condition p2: if p2=1, velocity in the intersection of the river surface bed and gauging vertical is same as in the top measured point in gauging vertical.}
  \item{p3}{boundary condition p3: if p3=1, transfororm the "coordinate system" of measured depth data to the cartesian.}
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
