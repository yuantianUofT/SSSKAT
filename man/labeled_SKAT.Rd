% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/labeled_SKAT.R
\name{labeled_SKAT}
\alias{labeled_SKAT}
\title{Standard SKAT function using the labeled data only}
\usage{
labeled_SKAT(dat, id.t)
}
\arguments{
\item{dat}{Basis matrix of data set to be analyzed.}

\item{id.t}{Row id of labeled data in the dat dataset.}
}
\value{
Parameter estimates, SS SKAT score and SS SKAT pvalue.
}
\description{
Standard SKAT function using the labeled data only
}
\examples{
analysis_data <- example_data$my_dat;
labeled_data_id <- example_data$id.t;
labeled_SKAT_result <- labeled_SKAT(dat = analysis_data, id.t = labeled_data_id);
}
