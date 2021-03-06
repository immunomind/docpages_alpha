% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis.R
\name{vis_box}
\alias{vis_box}
\title{Flexible box-plots for visualisation of distributions}
\usage{
vis_box(
  .data,
  .by = NA,
  .meta = NA,
  .melt = TRUE,
  .points = TRUE,
  .test = TRUE,
  .signif.label.size = 3.5,
  .defgroupby = "Sample",
  .grouping.var = "Group",
  .labs = c("X", "Y"),
  .title = "Boxplot (.title argument)",
  .subtitle = "Subtitle (.subtitle argument)",
  .legend = NA,
  .leg.title = "Legend (.leg.title argument)",
  .legend.pos = "right"
)
}
\arguments{
\item{.data}{Input matrix or data frame.}

\item{.by}{Pass NA if you want to plot samples without grouping.

You can pass a character vector with one or several column names from ".meta"
to group your data before plotting. In this case you should provide ".meta".

You can pass a character vector that exactly matches the number of samples in
your data, each value should correspond to a sample's property. It will be used
to group data based on the values provided. Note that in this case you should
pass NA to ".meta".}

\item{.meta}{A metadata object. An R dataframe with sample names and their properties,
such as age, serostatus or hla.}

\item{.melt}{If TRUE then apply \link{melt} to the ".data" before plotting.
In this case ".data" is supposed to be a data frame with the first character column reserved
for names of genes and other numeric columns reserved to counts or frequencies of genes.
Each numeric column should be associated with a specific repertoire sample.}

\item{.points}{A logical value defining whether points will be visualised or not.}

\item{.test}{A logical vector whether statistical tests should be applied. See "Details" for more information.}

\item{.signif.label.size}{An integer value defining the size of text for p-value.}

\item{.defgroupby}{A name for the column with sample names.}

\item{.grouping.var}{A name for the column to group by.}

\item{.labs}{Character vector of length two with names for x-axis and y-axis, respectively.}

\item{.title}{The text for the title of the plot.}

\item{.subtitle}{The The text for the plot's subtitle.}

\item{.legend}{If TRUE then displays a legend, otherwise removes legend from the plot.}

\item{.leg.title}{The The text for the plots's legend. Provide NULL to remove the legend's title completely.}

\item{.legend.pos}{Positions of the legend: either "top", "bottom", "left" or "right".}
}
\value{
A ggplot2 object.
}
\description{
Visualisation of distributions using ggplot2-based boxplots.
}
\examples{
vis_box(data.frame(Sample = sample(c("A", "B", "C"), 100, TRUE), Value = rnorm(100)), .melt = FALSE)
}
\seealso{
\link{vis.immunr_gene_usage}, \link{geneUsage}
}
\concept{vis}
