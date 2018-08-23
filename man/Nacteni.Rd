\name{Nacteni}
\alias{Nacteni}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Loading data
}
\description{Loading data from cross-section measurement of point velocities with FlowTracker HandHeld or hydrometric propeller.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Nacteni()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{}{
}
}
\details{
Function will open dialog for selection of files, that you want to use as input data (.txt file). It loads the data from the file. Input is data.frame. In data there must be columns: Loc, Depth, MeasD, Vel. Means: Loc - "x" location of point from the bank of the stream [meters], Depth - depth of vertical of measured point [meters], MeasD - depth of measured point in vertical [meters], Vel - measured point velocity [m*s^(-1)].  Measured verticals must be in order from the lowest value of "x" location to highest "x" location or from the highest value of "x" location to the lowest value of "x" location.

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
data = Nacteni()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
