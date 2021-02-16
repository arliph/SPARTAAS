mapclust_app <- function() {
  appDir <- system.file("shiny-examples", "map.app", package = "SPARTAAS")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `SPARTAAS`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
