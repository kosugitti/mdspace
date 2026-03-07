#' mdspace Shiny アプリケーションを起動する
#'
#' ローカル環境でMDS分析用のGUIを起動します。
#'
#' @param ... \code{\link[shiny]{runApp}} に渡す追加引数
#' @return なし（Shinyアプリを起動）
#' @export
run_app <- function(...) {
  app_dir <- system.file("shiny", package = "mdspace")
  if (app_dir == "") {
    stop("Shiny app not found. Try re-installing `mdspace`.")
  }
  shiny::runApp(app_dir, ...)
}
