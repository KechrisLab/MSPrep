\name{test}
\alias{test}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(test)}
\format{
  The format is:
List of 3
 $ sum_data1: num [1:9, 1:2654] 0 24885 23820 20730 19302 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:9] "1x_O1" "1x_O2" "1x_O3" "2x_O1" ...
  .. ..$ : chr [1:2654] "577.0322_0.5910416" "539.0207_0.58933336" "723.9585_0.599889" "525.3529_0.6085925" ...
 $ clinical :'data.frame':	9 obs. of  3 variables:
  ..$ SubjectID: Factor w/ 9 levels "1x_O1","1x_O2",..: 1 2 3 4 5 6 7 8 9
  ..$ Operator : int [1:9] 1 2 3 1 2 3 1 2 3
  ..$ Spike    : int [1:9] 1 1 1 2 2 2 4 4 4
 $ medians  : chr [1:23, 1:5] "251768" "101761" "79673" "468810" ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(test)
## maybe str(test) ; plot(test) ...
}
\keyword{datasets}
