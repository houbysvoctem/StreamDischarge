\name{StreamDischarge-package}
\alias{StreamDischarge-package}
\alias{StreamDischarge}
\docType{package}
\title{
\packageTitle{StreamDischarge}
}
\description{Package for evaluatuion of hydrometric measurement from the stream cross-sections.
\packageDescription{StreamDischarge}
}
\details{Package uncludes the functions for calculation of discharge and for visualistion of interpolated point velocitites in the cross-section. Package also includes some measurement data - data_01, data_02.

Functions for users:

\code{\link{Harlacher}} - calculation of discharge with the graphical Harlacher's method - usable for measurements with 3 and more points in the gauging vertical.

\code{\link{SvisP}} - calculation of discharge with the mid-section segments method (technical norm ISO 748)

\code{\link{MeziSvisP}} - calculation of discharge with the mean-section segments method (technical norm ISO 748)

\code{\link{discharge_interpolace}} - calculation of discharge with the interpolation, that was calibrated in Eliška´s dissertation thesis

\code{\link{print_crossection}} - from input measurement data interpolates the point velocites and shows then in the rectangular regular grid

}
\author{
\packageAuthor{StreamDischarge}

Maintainer: \packageMaintainer{StreamDischarge}
}
\references{
  ISO 748:2007 Hydrometry -- Measurement of liquid flow in open channels using current-meters or floats

  Pebesma, E.J., R.S. Bivand, 2005. Classes and methods for spatial data in R. R News 5 (2),
  https://cran.r-project.org/doc/Rnews/.

  Roger S. Bivand, Edzer Pebesma, Virgilio Gomez-Rubio, 2013. Applied spatial data analysis
  with R, Second edition. Springer, NY. http://www.asdar-book.org/

  Matt Dowle and Arun Srinivasan (2018). data.table: Extension of `data.frame`. R package
  version 1.11.0. https://CRAN.R-project.org/package=data.table

  Deepayan Sarkar and Felix Andrews (2016). latticeExtra: Extra Graphical Utilities Based on
  Lattice. R package version 0.6-28. https://CRAN.R-project.org/package=latticeExtra
}
%~~ Optionally other standard keywords, one per line, from file ~~
%~~ KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
#calculate discharge with the interpolation
interpolation = 1 #1 for the arithmetical mean
discharge_interpolace(data_01, interpolation)

#calculate discharge with Harlacher graphical method
m = 5
Harlacher(data_01, m)

#calculate discharge with mean-section segments ISO method
SvisP(data_01)

#calculate discharge with mid-section segments ISO method
m = 5
MeziSvisP(data_01, m)

#level plot of interpolated data set velocities
interpolation = 1
grid_size = 0.05
print_crossection(data_01,interpolation,grid_size)

%~~ simple examples of the most important functions ~~
}
