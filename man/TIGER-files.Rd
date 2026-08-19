\name{TIGER-files}
\alias{tiger_example_file}
\alias{tiger_template_file}
\title{Locate TIGER synthetic examples and input templates}
\description{
Return the path to a synthetic example file or hypothetical CSV template in an
installed TIGER package or repository checkout.
}
\usage{
tiger_example_file(filename, must_work = TRUE)

tiger_template_file(filename, must_work = TRUE)
}
\arguments{
  \item{filename}{One supplied filename, including its extension.}
  \item{must_work}{Whether to error when the requested file does not exist.}
}
\value{
A character path. Template filenames are
\code{individual_template.csv}, \code{rv_reference_template.csv},
\code{apoe_reference_template.csv}, and
\code{high_impact_reference_template.csv}.
}
\examples{
read.csv(tiger_template_file("individual_template.csv"))
read.csv(tiger_template_file("rv_reference_template.csv"))
}
