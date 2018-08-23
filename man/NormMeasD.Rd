\name{NormMeasD}
\alias{NormMeasD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for use of interpolation methods in the stream cross-section
%%  ~~function to do ... ~~
}
\description{Transformation of measured depths during point velocity measurement of cross-section to the coordinates with the zero value in the river surface and minus values minimizing in the direction to the river bed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NormMeasD(aa, ba)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aa}{vector of depths of gauging verticals in the stream cross-section.
%%     ~~Describe \code{aa} here~~
}
  \item{ba}{vector of the depths of single points in gauging verticals in the stream cross-section
%%     ~~Describe \code{ba} here~~
}
}
\details{Length of the vectors must be equal.
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
Depth = data_01$Depth #verticals, but each vertical must be n times, smae as the number of single points
MeasD = data_01$MeasD #single points
Depth
MeasD
NormMeasD(Depth,MeasD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
