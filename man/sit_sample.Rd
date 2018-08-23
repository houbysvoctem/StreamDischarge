\name{sit_sample}
\alias{sit_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for use of interpolation methods in the stream cross-section
%%  ~~function to do ... ~~
}
\description{Creates orthogonal regular mesh inside cross-section.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sit_sample(data, grid_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data must be in format: data.frame. In data there must be columns: Loc, Depth, MeasD, Vel. Means: Loc - "x" location of point from the bank of the stream [meters], Depth - depth of vertical of measuring point [meters], MeasD - depth of measuring point in vertical [meters], Vel - measured point velocity [m*s^(-1)]. Gauging verticals must be in order from the lowest value of "x" location to highest "x" location or from the highest value of "x" location to the lowest value of "x" location.
}
  \item{grid_size}{ size of grid of the mesh, in [m]
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{SIT}{coordinates of points in mesh}
\item{POLYGONCS}{polygon of cross-section plane}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Pebesma, E.J., R.S. Bivand, 2005. Classes and methods for spatial data in R. R News 5 (2),
  https://cran.r-project.org/doc/Rnews/.

  Roger S. Bivand, Edzer Pebesma, Virgilio Gomez-Rubio, 2013. Applied spatial data analysis
  with R, Second edition. Springer, NY. http://www.asdar-book.org/

  Matt Dowle and Arun Srinivasan (2018). data.table: Extension of `data.frame`. R package
  version 1.11.0. https://CRAN.R-project.org/package=data.table
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
The function uses packages \code{\link{sp}} and \code{\link{data.table}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
grid_size = 0.1
sit_sample(data_02, grid_size)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
