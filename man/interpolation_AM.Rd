\name{arith_mean}
\alias{arith_mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for use of interpolation methods in the steam cross-section
%%  ~~function to do ... ~~
}
\description{Interpolation of point velocities with the arithmetical mean method.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
arith_mean(data,sit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data must be in format: data.frame. In data there must be columns: Loc, Depth, MeasD, Vel. Means: Loc - "x" location of point from the bank of the stream [meters], Depth - depth of vertical of measuring point [meters], MeasD - depth of measuring point in vertical [meters], Vel - measured point velocity [m*s^(-1)]. Gauging verticals must be in order from the lowest value of "x" location to highest "x" location or from the highest value of "x" location to the lowest value of "x" location.
}
  \item{sit}{matrix of coordinates of points, were the estimation should be done
%%     ~~Describe \code{x} here~~
}
}
\details{Number of neighbouring points from that the interpolation is done = 2, the neighbourhood of estimated point is assumed elliptical, in the intersection of measured verticals and the water surface the velocita is assumed equal to the velocity of the top measured point velocity.
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
grid_size = 0.1
mesh = sit_sample(data_02, grid_size)
arith_mean(my_data, mesh)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
