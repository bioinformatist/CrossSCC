#' Plot CrossSCC result by the interactive viewer
#'
#' @param x CrossSCC result.
#'
#' @return a interactive plot, as a data.tree object.
#' @export
#'
#' @examples
#' data('cl.b1')
#' handsome.zuo <- CrossSCC(cl.b1[1:500,], ncores = 10, mean.posterior.cutoff = 0.18)
#' plot_CrossSCC(handsome.zuo)
plot_CrossSCC <- function(x) {
  # Should not simplify here, for a matrix will be returned as sapply does
  contents <- x$Get('sampleNames', simplify = FALSE)
  element2label <- function(x) {
    left.side <- "<DIV align=center
        style='
        color: #ffffff;
        background-color: #000000;
        border: solid 2px black;
        width: 300px;
        height: 200px;
        overflow: scroll;
        scrollbar-face-color: #889B9F;
        scrollbar-shadow-color: #3D5054;
        scrollbar-highlight-color: #C3D6DA;
        scrollbar-3dlight-color: #3D5054;
        scrollbar-darkshadow-color: #85989C;
        scrollbar-track-color: #95A6AA;
        scrollbar-arrow-color: #FFD6DA;
        '>"

    middle <- ifelse(is(x, 'character'), paste0(paste0("<p>", x, "</p>"), collapse = ''), paste0(c('<p><b>Component 1:</b></p>', paste0("<p>", x[[1]], "</p>"), '<p><b>Component 2:</b></p>', paste0("<p>", x[[2]], "</p>")), collapse = ''))
    paste0(left.side, middle, '</DIV>', collapse = "")
  }
  titles <- vapply(contents, element2label, character(1))
  tmp.plot <- plot(x, output = 'visNetwork')
  tmp.plot$x$nodes$title <- factor(titles)
  tmp.plot %>% visNodes(image = 'https://reactome.org/icon/R-ICO-013594.png', shape = 'circularImage') %>% visHierarchicalLayout(sortMethod = 'directed') %>% visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
}
