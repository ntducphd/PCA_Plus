
library(shiny)
library(FactoMineR)
library(factoextra)
library(DT)
library(ggplot2)
library(plotly)

ui <- fluidPage(
  titlePanel("PCA-Plus — Developed by Vietnamese Association at Kyushu University"),

  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", accept = ".csv"),
      uiOutput("varselect"),
      checkboxInput("scale", "Scale data", TRUE),
      numericInput("plotWidth", "Plot width (inches)", 8, min = 4),
      numericInput("plotHeight", "Plot height (inches)", 6, min = 4),
      selectInput("fontFamily", "Font family", choices = c("sans", "serif", "mono", "Times", "Arial", "Helvetica")),
      downloadButton("downloadPCA", "Download PCA Coordinates"),
      downloadButton("downloadPCAPlot", "Download PCA Biplot"),
      downloadButton("downloadHCPCPlot", "Download HCPC Dendrogram"),
      downloadButton("downloadScreePlot", "Download Scree Plot"),
      downloadButton("downloadVarPlot", "Download Variable Plot")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Eigenvalues", DTOutput("eigTable")),
        tabPanel("Scree Plot", plotOutput("screePlot", height = "600px")),
        tabPanel("Variables Plot", plotOutput("varPlot", height = "600px")),
        tabPanel("PCA Biplot", plotOutput("pcaPlot", height = "600px")),
        tabPanel("HCPC Dendrogram", plotOutput("hcpcPlot", height = "600px"))
      )
    )
  )
)

server <- function(input, output, session) {

  dataInput <- reactive({
    req(input$file)
    read.csv(input$file$datapath)
  })

  output$varselect <- renderUI({
    req(dataInput())
    numericVars <- names(Filter(is.numeric, dataInput()))
    checkboxGroupInput("selectedVars", "Select numeric traits:", 
                       choices = numericVars, selected = numericVars)
  })

  pcaResult <- reactive({
    req(input$selectedVars)
    df <- dataInput()[, input$selectedVars, drop = FALSE]
    PCA(df, scale.unit = input$scale, graph = FALSE)
  })

  hcpcResult <- reactive({
    req(pcaResult())
    HCPC(pcaResult(), graph = FALSE)
  })

  output$eigTable <- renderDT({
    eig <- as.data.frame(pcaResult()$eig)
    colnames(eig) <- c("Eigenvalue", "Proportion", "Cumulative")
    datatable(round(eig, 3))
  })

  output$screePlot <- renderPlot({
    fviz_eig(pcaResult(), addlabels = TRUE, barfill = "steelblue", barcolor = "black") +
      theme_minimal(base_size = 14, base_family = input$fontFamily)
  })

  output$varPlot <- renderPlot({
    fviz_pca_var(
      pcaResult(), col.var = "contrib",
      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
      repel = TRUE
    ) +
    theme_minimal(base_size = 14, base_family = input$fontFamily)
  })

  output$pcaPlot <- renderPlot({
    fviz_pca_biplot(
      pcaResult(), repel = TRUE,
      ggtheme = theme_classic(base_size = 14, base_family = input$fontFamily)
    )
  })

  output$hcpcPlot <- renderPlot({
    fviz_dend(
      hcpcResult(), 
      cex = 1.2, 
      palette = "jco", 
      rect = TRUE, 
      rect_fill = TRUE,
      ggtheme = theme_minimal(base_size = 14, base_family = input$fontFamily)
    )
  })

  output$downloadPCA <- downloadHandler(
    filename = function() {
      paste0("PCA_Coordinates_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(pcaResult()$ind$coord, file)
    }
  )

  output$downloadPCAPlot <- downloadHandler(
    filename = function() {
      paste0("PCA_Biplot_", Sys.Date(), ".png")
    },
    content = function(file) {
      g <- fviz_pca_biplot(
        pcaResult(), repel = TRUE,
        ggtheme = theme_classic(base_size = 16, base_family = input$fontFamily)
      )
      ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
    }
  )

  output$downloadHCPCPlot <- downloadHandler(
    filename = function() {
      paste0("HCPC_Dendrogram_", Sys.Date(), ".png")
    },
    content = function(file) {
      g <- fviz_dend(
        hcpcResult(), 
        cex = 1.2, 
        palette = "jco", 
        rect = TRUE, 
        rect_fill = TRUE,
        ggtheme = theme_minimal(base_size = 16, base_family = input$fontFamily)
      )
      ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
    }
  )

  output$downloadScreePlot <- downloadHandler(
    filename = function() {
      paste0("Scree_Plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      g <- fviz_eig(pcaResult(), addlabels = TRUE, barfill = "steelblue", barcolor = "black") +
        theme_minimal(base_size = 16, base_family = input$fontFamily)
      ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
    }
  )

  output$downloadVarPlot <- downloadHandler(
    filename = function() {
      paste0("PCA_Variable_Plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      g <- fviz_pca_var(
        pcaResult(), col.var = "contrib",
        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
        repel = TRUE
      ) +
      theme_minimal(base_size = 16, base_family = input$fontFamily)
      ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
    }
  )
}

shinyApp(ui, server)
