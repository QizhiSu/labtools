#' An internal function to plot molecule
#'
#' \code{plot_molecule} is a simple wrapper to plot chemical structure based on
#' SMILES.
#'
#' @param molecule an object as returned by rcdk::load.molecules or rcdk::parse.smiles
#' @param name a character for the name of the molecule
#' @param sma a character with the smarts string as passed onto get.depictor
#' @param ... other arguments for get.depictor
#'
#' @return an image showing the chemical structure
#' @rawNamespace import(rcdk, except = "matches")

plot_molecule <- function(molecule, name = NULL, sma = NULL, ...){
  # Image aesthetics
  dep <- rcdk::get.depictor(
    width = 1000,
    height = 1000,
    zoom = 7,
    sma = sma,
    abbr = "off",
    ...
  )
  molecule_sdf <- rcdk::view.image.2d(molecule[[1]], depictor = dep)

  ## Remove extra margins around the molecule
  par(mar=c(0,0,0,0.5))
  plot(NA,
       xlim=c(0, 1),
       ylim=c(0, 1),
       # Remove the black bounding boxes around the molecule
       axes = F,
       type = "n")
  rasterImage(molecule_sdf, 0, 0, 1, 1)
  # Annotate the molecule
  # text(x = 0.5, y = 1.05,  deparse(substitute(molecule)))
}

#-------------------------------------------------------------------------------
#' A shinyApp to navigate through a chemical table
#'
#' \code{navigate_chem} is a shinyApp to navigate through a chemical table and
#' to see the structure of each one. The data must contain SMILES of each chemical.
#'
#'
#' @param data A dataframe contain at least SMILES of the chemical
#'
#' @export
#'
#' @import DT
#' @rawNamespace import(shiny, except = c("dataTableOutput", "renderDataTable"))
#' @rawNamespace import(shinyjs, except = c("show", "runExample"))
#' @rawNamespace import(rcdk, except = "matches")

navigate_chem <- function(data) {
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    warning("The 'rcdk' package is not installed. ",
            "Some functionality may be limited. ",
            "You can install it using: install.packages('rcdk')")
  }

  js_select_dt <- c(
    "var dt = table.table().node();",
    "var tblID = $(dt).closest('.datatables').attr('id');",
    "var inputName = tblID + '_rows_selected'",
    "var incrementName = tblID + '_rows_selected2_increment'",
    "table.on('key-focus', function(e, datatable, cell, originalEvent){",
    "  if (originalEvent.type === 'keydown'){",
    "    table.rows().deselect(); ",
    "    table.row(cell[0][0].row).select();",
    "    row = table.rows({selected: true})",
    "    Shiny.setInputValue(inputName, [parseInt(row[0]) + 1]);",
    "  }",
    "});"
  )

  ui <- fluidPage(
    shinyjs::useShinyjs(),  # Initialize shinyjs
    tags$head(
      tags$style(HTML("
        #structurePlot {
          position: absolute;
          left: 0;
          top: 0;
          width: 100%;
          height: 100%;
        }
      ")),
      tags$script(src = "https://code.jquery.com/ui/1.12.1/jquery-ui.js")
    ),
    textOutput("selectedRow"),
    fluidRow(
      column(width = 8,
             DTOutput("chemicalTable")
      ),
      column(width = 4,
             div(id = "structureContainer",
                 style = "display: flex; justify-content: flex-end;",
                 plotOutput("structurePlot", width = "100%", height = "400px")))
    )
  )

  server <- function(input, output) {

    output$chemicalTable <- renderDT({
      # Find the column index of the SMILES column
      smiles_col_index <- grep("SMILES", colnames(data), ignore.case = TRUE)

      # Hide the SMILES column if found
      if (length(smiles_col_index) > 0) {
        data[, smiles_col_index] <- ""
      }

      datatable(
        data,
        selection = "single",
        callback = JS(js_select_dt),
        extensions = c("KeyTable", "Select"),
        editable = TRUE,
        options = list(
          keys = TRUE,
          select = TRUE,
          columnDefs = list(list(visible = FALSE, targets = smiles_col_index))
        )
      )
    }, server = FALSE)

    output$structurePlot <- renderPlot({
      selectedRow <- input$chemicalTable_rows_selected

      if (!is.null(selectedRow)) {
        selectedRow <- as.integer(selectedRow)
        SMILES <- data$SMILES[selectedRow]

        tryCatch({
          mol <- rcdk::parse.smiles(SMILES)

          if (!is.null(mol)) {
            par(mar = c(0, 0, 0, 0))
            plot_molecule(mol)
          } else {
            plot.new()
            text(x = 0.5,
                 y = 0.5,
                 labels = "Error: Invalid SMILES string",
                 cex = 1.2)
          }
        }, error = function(e) {
          plot.new()
          text(x = 0.5,
               y = 0.5,
               labels = "Error: Failed to render the structure",
               cex = 1.2)
        })
      }
    })

    output$selectedRow <- renderText(input$chemicalTable_rows_selected)

    # Make the structurePlot resizable and draggable
    shinyjs::runjs('
      $("#structurePlot").resizable();
      $("#structurePlot").draggable();
    ')
  }

  shinyApp(ui, server)
}
