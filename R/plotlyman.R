plotlyman <- function(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", sigLine = 8) {

    # Check input
    if (!is.data.frame(x)) {
        stop("'x' has to be a data.frame.")
    }
    if (!chr %in% names(x)) {
        stop("'x' has to contain a '", chr, "' column")
    }
    if (!bp %in% names(x)) {
        stop("'x' has to contain a '", bp, "' column")
    }
    if (!p %in% names(x)) {
        stop("'x' has to contain a '", p, "' column")
    }
    if (!snp %in% names(x)) {
        stop("'x' has to contain a '", snp, "' column")
    }
    if (!is.numeric(x[[chr]])) {
        stop("'chr' column of 'x' has to be numeric.")
    }
    if (!is.numeric(x[[bp]])) {
        stop("'bp' column of 'x' has to be numeric.")
    }
    if (!is.numeric(x[[p]])) {
        stop("'p' column of 'x' has to be numeric.")
    }
    if (!is.character(x[[snp]])) {
        stop("'snp' column of 'x' has to be alphanumeric.")
    }

    # Reformat data
    assoc <- data.frame(
        CHR = x[[chr]],
        BP = x[[bp]],
        P = x[[p]],
        SNP = x[[snp]],
        stringsAsFactors = FALSE
    )

    # Sort data
    assoc <- assoc[order(assoc[["CHR"]], assoc[["BP"]]), ]

    # Convert p-values to -log10
    assoc[["P"]] <- -log10(assoc[["P"]])

    bpMax <- tapply(assoc[[bp]], assoc[["CHR"]], max)
    bpOffset <- cumsum(c(0, bpMax[seq(1, length(bpMax) - 1)]))
    assoc$POS <- bpOffset[assoc[["CHR"]]] + assoc[[bp]]

    # Determine ticks
    tickPositions <- tapply(assoc[["POS"]], assoc[["CHR"]], function(values) {
        min(values) + ((max(values) - min(values)) / 2)
    })
    tickNames <- unique(assoc[["CHR"]])

    plot <- suppressWarnings(plot_ly(
        data = assoc,
        x = ~POS,
        y = ~P,
        text = ~SNP,
        fillcolor = ~CHR,
        meta = ~BP, # for use in hovertemplate
        type = "scattergl",
        mode = "markers",
        marker = list(
            symbol = "circle-open"
        ),
        # There will be as many copies of this text as p-values...
        hovertemplate = paste(
            c(
                "ID: %{text}",
                "Chromosome: %{data.fillcolor}",
                "Base-pair position: %{meta:d}",
                "<extra></extra>" # hide extra box
            ),
            collapse = "<br>"
        )
    )) %>% layout(
        title = "Manhattan Plot",
        xaxis = list(
            title = "Chromosome",
            tickmode = "array",
            tickvals = tickPositions,
            ticktext = tickNames,
            tickangle = 0,
            showgrid = FALSE
        ),
        yaxis = list(
            title = "-log10(p-value)"
        ),
        showlegend = FALSE,
        shapes = list(
            list( # significance line
                type = "line",
                x0 = 0,
                x1 = ~max(POS),
                y0 = sigLine,
                y1 = sigLine,
                line = list(
                    color = "red"
                ),
                opacity = 0.5
            )
        ),
        # https://colorbrewer2.org/#type=diverging&scheme=PuOr&n=3
        colorway = list(
            "#f1a340",
            "#998ec3"
        )
    )

    return(plot)
}
