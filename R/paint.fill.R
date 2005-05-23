"paint.fill" <-
function (m, startloc = NULL, plot = TRUE) 
{
    k = list()
    if (plot) {
        k = lapply(1:length(dim(a)), function(i) if (i == 1 || 
            i == 2) 
            1:dim(a)[i]
        else round(dim(a)[i]/2))
        image(1:length(k[[1]]), 1:length(k[[2]]), do.call("[", 
            append(list(m), k)))
    }
    if (is.null(startloc) && plot) {
        startloc = lapply(locator(), round)
        k[[1]] = startloc[[1]]
        k[[2]] = startloc[[2]]
        startloc = k
    }
    if (is.null(startloc)) 
        stop("No starting element provided.")
    if (length(dim(m)) != length(startloc)) 
        stop(paste("Dimension of first argument `", deparse(substitute(m)), 
            "' incompatible with length of \nsecond argument '", 
            deparse(substitute(startloc)), "'", sep = ""))
    cl = do.call("cbind", startloc)
    testvalue = entry.value(m, cl)[1]
    neighb = as.matrix(do.call("expand.grid", lapply(dim(m), 
        function(b) -1:1)))
    checked = array(FALSE, dim = dim(m))
    boundary = array(dim = c(0, length(dim(m))))
    on.exit(return(invisible(list(boundary = boundary, painted = checked))))
    while (nrow(cl) > 0) {
        if (plot) 
            points(cl, pch = 20)

        # generate new elements
        nl = rep(1, nrow(neighb)) %x% cl + neighb %x% rep(1, 
            nrow(cl))
        nl = Unique(nl)
        outside = ((nl < 1) | t(t(nl) > dim(m))) %*% rep(1, ncol(nl))
        nl = nl[!outside, , drop = FALSE]                   # remove invalid elements
        nl = nl[!entry.value(checked, nl), , drop = FALSE]  # remove already checked elements
        if (nrow(nl) <= 0) 
            break

        # check image for boundary at new elements
        is.bound.elem = entry.value(m, nl) != testvalue
        if (any(is.bound.elem)) 
            boundary = rbind(boundary, nl[is.bound.elem, ])

        # make new elements the current elements
        entry.value(checked, nl) = TRUE                     # these elements have been checked
        cl = nl[!is.bound.elem, , drop = FALSE]
    }
}

