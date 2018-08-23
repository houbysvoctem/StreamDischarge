\name{print_crossection}
\alias{print_crossection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphic view
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print_crossection(data, interpolace, oko)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data must be in format: data.frame. In data there must be columns: Loc, Depth, MeasD, Vel. Means: Loc - "x" location of point from the bank of the stream [meters], Depth - depth of vertical of measuring point [meters], MeasD - depth of measuring point in vertical [meters], Vel - measured point velocity [m*s^(-1)]. Gauging verticals must be in order from the lowest value of "x" location to highest "x" location or from the highest value of "x" location to the lowest value of "x" location.
}
\item{interpolace}{interpolace = 1 for the aritmetical mean, 2 for the inverse distance weighted, 3 for the ordinary kriging interpolation}
\item{oko}{size of the grid, in [m]}
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
  Deepayan Sarkar and Felix Andrews (2016). latticeExtra: Extra Graphical Utilities Based on
  Lattice. R package version 0.6-28. https://CRAN.R-project.org/package=latticeExtra


%% ~put references to the literature/web site here ~
}
\author{Eliška Friedlová (Hasníková)
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See Also \code{\link{levelplot}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
print_crossection(data_01,2,0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
