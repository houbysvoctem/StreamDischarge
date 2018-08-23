\name{HorizontalV}
\alias{HorizontalV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calculation of mean velocities in verticals of cross-section by CSN/EU norm
}
\description{
The function calculates mean velocities in verticals of cross-section. Input data are originaly inteded for measurement of point velocity hydrometric measurement with FlowTracker HandHeld ADV, but if the data is provided in a right format, they can be from measurement with any other point velocity tool. Output is data.table with depths of verticals in [m] and mean velocities in verticals in [m*s^(-1)].
}
\usage{
HorizontalV(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data must be in format: data.frame. In data there must be columns: Loc, Depth, MeasD, Vel. Means: Loc - "x" location of point from the bank of the stream [meters], Depth - depth of vertical of measured point [meters], MeasD - depth of measured point in vertical [meters], Vel - measured point velocity [m*s^(-1)].  Measured verticals must be in order from the lowest value of "x" location to highest "x" location or from the highest value of "x" location to the lowest value of "x" location.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Number of points in vertical can be from 1 to 6. Used equations are:

1 point -- velocity of vertical = v1

2 points -- velocity of vertical = 0,5*(v(0,2) + v(0,8))

3 points -- velocity of vertical = 0.25*(v(0,2) + 2*v(0,4) + v(0,8))

4 points -- velocity of vertical = 0.2*(v(0,2) + v(0,4) + 2*v(0,6) + v(0,8))

5 points -- velocity of vertical = 0.1*(v(bottom) + 2*v(0,2) + 3*v(0,4) + 3*v(0,8) + v(top))

6 points -- velocity of vertical = 0.1*(v(bottom) + 2*v(0,2) + 2*v(0,4) + 2*v(0,6) + 2*v(0,8) + v(top))
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
technical norm ISO 748
}
\author{
Eliska Friedlova (Hasnikova)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
This function is used in functions \code{\link{SvisP}}, \code{\link{MeziSvisP}}.
}
\examples{
Velocities = HorizontalV(my_data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
