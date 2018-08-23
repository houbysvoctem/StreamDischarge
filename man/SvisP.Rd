\name{SvisP}
\alias{SvisP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calculation of discharge with method of mid-section segment method by CNS/EU norm
}
\description{
Input data are originaly inteded for measurement of point velocity hydrometric measurement with FlowTracker HandHeld ADV, but if the data is provided in a right format, they can be from measurement with any other point velocity tool. Output is discharge in [(m^3)*s^(-1)].
}
\usage{
SvisP(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data must be in format: data.frame. In data there must be columns: Loc, Depth, MeasD, Vel. Means: Loc - "x" location of point from the bank of the stream [meters], Depth - depth of vertical of measured point [meters], MeasD - depth of measured point in vertical [meters], Vel - measured point velocity [m*s^(-1)]. Measured verticals must be in order from the lowest value of "x" location to highest "x" location or from the highest value of "x" location to the lowest value of "x" location.
}
}
\details{
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
this function uses \code{\link{HorizontalV}}
}
\examples{
Qmids = SvisP(data_01)
Qmids
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
