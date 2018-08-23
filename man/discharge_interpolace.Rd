\name{discharge_interpolace}
\alias{discharge_interpolace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for use of interpolation methods in the stream cross-section
%%  ~~function to do ... ~~
}
\description{Discharge from the interpolated values of point velocities in the river cross-section.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
discharge_interpolace(data, interpolace)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data must be in format: data.frame. In data there must be columns: Loc, Depth, MeasD, Vel. Means: Loc - "x" location of point from the bank of the stream [meters], Depth - depth of vertical of measured point [meters], MeasD - depth of measured point in vertical [meters], Vel - measured point velocity [m*s^(-1)].  Measured verticals must be in order from the lowest value of "x" location to highest "x" location or from the highest value of "x" location to the lowest value of "x" location.
}
%%     ~~Describe \code{x} here~~
  \item{interpolace}{interpolace = 1 for the aritmetical mean, 2 for the inverse distance weighted, 3 for the ordinary kriging interpolation

%%     ~~Describe \code{x} here~~
}
}
\details{The parameters of the interpolation models were calibrated on 45 data sets of small Czech streams, with the DE/rand/i/bin optimizator. Interpolation is done in the ortoghonal regular mesh with calibrated size of grid. Volume of discharge is calculated as sum of blocks with the top in the point velocity and the bottom in the cross-section plane.
%%  ~~ If necessary, more details than the description above ~~
}
\value{Discharge in m^3.s^-1.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Eliska Friedlova (Hasnikova)
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See also \code{\link{arith_mean}}, \code{\link{ind_dis_w}},  \code{\link{ord_kriging}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Q = discharge_interpolace(data_01, 1)
Q
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
