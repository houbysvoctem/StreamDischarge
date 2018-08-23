\name{Harlacher}
\alias{Harlacher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calculation of discharge with graphical Harlacher method
}
\description{
The function calculates discharge of river cross-section. Input data are originaly inteded for measurement of point velocity hydrometric measurement with FlowTracker HandHeld ADV, but if the data is provided in a right format, they can be from measurement with any other point velocity tool. Output is discharge in (m^3)*s^(-1).

}
\usage{
Harlacher(data, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data must be in format: data.frame. In data there must be columns: Loc, Depth, MeasD, Vel. Means: Loc - "x" location of point from the bank of the stream [meters], Depth - depth of vertical of measured point [meters], MeasD - depth of measured point in vertical [meters], Vel - measured point velocity [m*s^(-1)]. Measured verticals must be in order from the lowest value of "x" location to highest "x" location or from the highest value of "x" location to the lowest value of "x" location.
}
\item{m}{Coefficient depending on hydraulic friction for estimation of velocity in the area between the reiver bank an first/last measured vertical. Can be calculated according to technical norm ISO 748 with use of Chezy coeficient, or can be estimated as a value the usual from interval 5--7. According to ISO 748, value 4 us used for high roughness of the bottom and value 10 means the least rough bottom. If you do not want to use coefficient m, set it to 0.

}
}
\details{The function produces linear connection of end pointd of velocity vectors in gauging verticals. And linear connections between mean velocities in gauging vertivals times depth in gauging vertical.
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
}
\examples{
Qhar = Harlacher(data_01, 5)
Qhar
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")

