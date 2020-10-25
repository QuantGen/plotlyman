\name{plotlyman}
\alias{plotlyman}
\title{Create a Manhattan plot using Plotly}
\description{
    Create a Manhattan plot using Plotly.
}
\usage{
plotlyman(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", sigLine = 8)
}
\arguments{
    \item{x}{
        A data frame that contains the GWAS results.
    }
    \item{chr}{
        The name of the column in \code{x} that contains the chromosomes.
    }
    \item{bp}{
        The name of the column in \code{x} that contains the base-pair
        positions.
    }
    \item{p}{
        The name of the column in \code{x} that contains the p-values
    }
    \item{snp}{
        The name of the column in \code{x} that contains the variant names.
    }
    \item{sigLine}{
        Where to draw a genome-wide significant line in -log10 scale.
    }
}
