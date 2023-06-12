cerardat_app <- function() {
  appDir <- system.file("shiny-examples", "cerardat.app", package = "SPARTAAS")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `SPARTAAS`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
