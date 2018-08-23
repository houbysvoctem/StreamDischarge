\name{testinv}
\alias{testinv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for use of interpolation methods in the stream cross-section
%%  ~~function to do ... ~~
}
\description{Function to testing singularity of matrix entered to the function.


%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
testinv(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{matrix to be tested
%%     ~~Describe \code{x} here~~
}
}
\details{Returns FALSE, if the matrix is not invertable, which means the matrix is singular.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{from webpage https://stackoverflow.com/questions/24961983/how-to-check-if-a-matrix-has-an-inverse-in-the-r-language
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
