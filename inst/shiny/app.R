# mdspace Shiny Application
# このファイルは mdspace::run_app() から起動されます

# モジュール読み込み
module_dir <- file.path(dirname(sys.frame(1)$ofile %||% "."), "modules")
if (dir.exists(module_dir)) {
  module_files <- list.files(module_dir, pattern = "\\.R$", full.names = TRUE)
  for (f in module_files) source(f, local = TRUE)
}

ui <- shiny::navbarPage(
  title = "mdspace",
  # タブは modules/ の実装に合わせて追加していく
  shiny::tabPanel("Data", shiny::h3("Data import module (TODO)"))
)

server <- function(input, output, session) {
  # モジュールサーバーをここで呼び出す
}

shiny::shinyApp(ui = ui, server = server)
