
library(shiny)
library(FactoMineR)
library(factoextra)
library(DT)
library(ggplot2)
library(plotly)
library(psych)
library(corrplot)

ui <- fluidPage(
  titlePanel("PCA & HCPC Analysis ver.1.0 — Developed by Vietnamese Association at Kyushu University"),

  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", accept = ".csv"),
      uiOutput("varselect"),
      checkboxInput("scale", "Scale data", TRUE),
      numericInput("plotWidth", "Plot width (inches)", 8),
      numericInput("plotHeight", "Plot height (inches)", 6),
      selectInput("fontFamily", "Font family", choices = c("sans", "serif", "mono", "Times", "Arial", "Helvetica")),
      numericInput("fontSize", "Font size (pt)", 14, min = 8),
      numericInput("topN", "Top variables/individuals to show (contrib)", 10, min = 1),
      
      radioButtons("exportFormat", "Export Format", choices = c("png", "tiff"), inline = TRUE),
      downloadButton("downloadPCAplotImg", "Download PCA Biplot Image"),
      downloadButton("download_screePlot", "Download Scree Plot"),
      downloadButton("download_varPlot", "Download Variable Plot"),
      downloadButton("download_indPlot", "Download Individuals Plot"),
      downloadButton("download_pcaPlot", "Download PCA Biplot"),
      downloadButton("download_contribVarPlot", "Download Contribution Variables"),
      downloadButton("download_contribIndPlot", "Download Contribution Individuals"),
      downloadButton("download_cos2Plot", "Download Cos2 Heatmap"),
      downloadButton("download_hcpcPlot", "Download HCPC Dendrogram")
,
      downloadButton("downloadAllTables", "Download All PCA Tables (.csv)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Eigenvalues", DTOutput("eigTable")),
        tabPanel("Scree Plot", plotOutput("screePlot", height = "600px")),
        tabPanel("Individuals", plotOutput("indPlot", height = "600px")),
        tabPanel("Variables", plotOutput("varPlot", height = "600px")),
        tabPanel("Biplot", plotOutput("pcaPlot", height = "600px")),
        tabPanel("Contribution (Variables)", plotOutput("contribVarPlot", height = "600px")),
        tabPanel("Contribution (Individuals)", plotOutput("contribIndPlot", height = "600px")),
        tabPanel("Cos2 Heatmap", plotOutput("cos2Plot", height = "600px")),
        tabPanel("HCPC Dendrogram", plotOutput("hcpcPlot", height = "600px")),
        tabPanel("KMO & Bartlett Test", verbatimTextOutput("testOutput"))
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
    checkboxGroupInput("selectedVars", "Select numeric traits:", choices = numericVars, selected = numericVars)
  })

  pcaResult <- reactive({
    req(input$selectedVars)
    df <- dataInput()[, input$selectedVars]
    PCA(df, scale.unit = input$scale, graph = FALSE)
  })

  hcpcResult <- reactive({
    req(pcaResult())
    HCPC(pcaResult(), graph = FALSE)
  })

  base_theme <- reactive({
    theme_minimal(base_size = input$fontSize, base_family = input$fontFamily)
  })

  output$eigTable <- renderDT({
    eig <- as.data.frame(pcaResult()$eig)
    colnames(eig) <- c("Eigenvalue", "Proportion", "Cumulative")
    datatable(round(eig, 3))
  })

  output$screePlot <- renderPlot({
    fviz_eig(pcaResult(), addlabels = TRUE, barfill = "#0073C2FF", barcolor = "black") +
      base_theme()
  })

  output$indPlot <- renderPlot({
    fviz_pca_ind(pcaResult(), col.ind = "cos2",
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
      base_theme()
  })

  output$varPlot <- renderPlot({
    fviz_pca_var(pcaResult(), 
                 col.var = "contrib", 
                 gradient.cols = c("deepskyblue3", "orange", "red3"),
                 repel = TRUE, labelsize = 5, arrowsize = 0.8) +
      theme_minimal(base_size = input$fontSize, base_family = input$fontFamily) +
      labs(title = "PCA Variable Contributions")
  })

  output$pcaPlot <- renderPlot({
    fviz_pca_biplot(pcaResult(), repel = TRUE, labelsize = 5) +
      theme_classic(base_size = input$fontSize, base_family = input$fontFamily)
  })

  output$contribVarPlot <- renderPlot({
    fviz_contrib(pcaResult(), choice = "var", axes = 1, top = input$topN) +
      base_theme()
  })

  output$contribIndPlot <- renderPlot({
    fviz_contrib(pcaResult(), choice = "ind", axes = 1, top = input$topN) +
      base_theme()
  })

  output$cos2Plot <- renderPlot({
    cos2 <- pcaResult()$var$cos2
    corrplot(as.matrix(cos2), is.corr = FALSE, method = "color", tl.col = "black", mar = c(1,1,1,1))
  })

  output$hcpcPlot <- renderPlot({
    fviz_dend(hcpcResult(), rect = TRUE, show_labels = TRUE, cex = 1.2)
  })

  output$testOutput <- renderPrint({
    df <- dataInput()[, input$selectedVars]
    kmo <- KMO(cor(df))
    bartlett <- cortest.bartlett(cor(df), n = nrow(df))
    list(KMO = kmo, Bartlett = bartlett)
  })

  
  output$downloadPCAplotImg <- downloadHandler(
    filename = function() {
      paste0("PCA_Biplot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      g <- fviz_pca_biplot(pcaResult(), repel = TRUE, labelsize = 5) +
        theme_classic(base_size = input$fontSize, base_family = input$fontFamily)
      if (input$exportFormat == "png") {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
      } else {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
      }
    }
  )


  
  output$download_screePlot <- downloadHandler(
    filename = function() {
      paste0("Scree_Plot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      g <- isolate(output$screePlot)$func()
      if (input$exportFormat == "png") {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
      } else {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
      }
    }
  )

  output$download_varPlot <- downloadHandler(
    filename = function() {
      paste0("Variable_Plot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      g <- isolate(output$varPlot)$func()
      if (input$exportFormat == "png") {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
      } else {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
      }
    }
  )

  output$download_indPlot <- downloadHandler(
    filename = function() {
      paste0("Individuals_Plot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      g <- isolate(output$indPlot)$func()
      if (input$exportFormat == "png") {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
      } else {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
      }
    }
  )

  output$download_pcaPlot <- downloadHandler(
    filename = function() {
      paste0("PCA_Biplot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      g <- isolate(output$pcaPlot)$func()
      if (input$exportFormat == "png") {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
      } else {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
      }
    }
  )

  output$download_contribVarPlot <- downloadHandler(
    filename = function() {
      paste0("Contribution_Variables_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      g <- isolate(output$contribVarPlot)$func()
      if (input$exportFormat == "png") {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
      } else {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
      }
    }
  )

  output$download_contribIndPlot <- downloadHandler(
    filename = function() {
      paste0("Contribution_Individuals_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      g <- isolate(output$contribIndPlot)$func()
      if (input$exportFormat == "png") {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
      } else {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
      }
    }
  )

  output$download_cos2Plot <- downloadHandler(
    filename = function() {
      paste0("Cos2_Heatmap_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      g <- isolate(output$cos2Plot)$func()
      if (input$exportFormat == "png") {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
      } else {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
      }
    }
  )

  output$download_hcpcPlot <- downloadHandler(
    filename = function() {
      paste0("HCPC_Dendrogram_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      g <- isolate(output$hcpcPlot)$func()
      if (input$exportFormat == "png") {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
      } else {
        ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
      }
    }
  )

  
  output$download_screePlot <- downloadHandler(
    filename = function() {
      paste0("screePlot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      if ("corrplot" %in% class(fviz_eig)) {
        tiff(file, width = input$plotWidth * 100, height = input$plotHeight * 100, res = 100)
        fviz_eig(pcaResult(), addlabels = TRUE, barfill = "#0073C2FF", barcolor = "black") + theme_minimal(base_size = input$fontSize, base_family = input$fontFamily)
        dev.off()
      } else {
        g <- fviz_eig(pcaResult(), addlabels = TRUE, barfill = "#0073C2FF", barcolor = "black") + theme_minimal(base_size = input$fontSize, base_family = input$fontFamily)
        if (input$exportFormat == "png") {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
        } else {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
        }
      }
    }
  )

  output$download_varPlot <- downloadHandler(
    filename = function() {
      paste0("varPlot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      if ("corrplot" %in% class(fviz_pca_var)) {
        tiff(file, width = input$plotWidth * 100, height = input$plotHeight * 100, res = 100)
        fviz_pca_var(pcaResult(), col.var = "contrib", gradient.cols = c("deepskyblue3", "orange", "red3"), repel = TRUE, labelsize = 5, arrowsize = 0.8) + theme_minimal(base_size = input$fontSize, base_family = input$fontFamily) + labs(title = "PCA Variable Contributions")
        dev.off()
      } else {
        g <- fviz_pca_var(pcaResult(), col.var = "contrib", gradient.cols = c("deepskyblue3", "orange", "red3"), repel = TRUE, labelsize = 5, arrowsize = 0.8) + theme_minimal(base_size = input$fontSize, base_family = input$fontFamily) + labs(title = "PCA Variable Contributions")
        if (input$exportFormat == "png") {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
        } else {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
        }
      }
    }
  )

  output$download_indPlot <- downloadHandler(
    filename = function() {
      paste0("indPlot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      if ("corrplot" %in% class(fviz_pca_ind)) {
        tiff(file, width = input$plotWidth * 100, height = input$plotHeight * 100, res = 100)
        fviz_pca_ind(pcaResult(), col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) + theme_minimal(base_size = input$fontSize, base_family = input$fontFamily)
        dev.off()
      } else {
        g <- fviz_pca_ind(pcaResult(), col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) + theme_minimal(base_size = input$fontSize, base_family = input$fontFamily)
        if (input$exportFormat == "png") {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
        } else {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
        }
      }
    }
  )

  output$download_pcaPlot <- downloadHandler(
    filename = function() {
      paste0("pcaPlot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      if ("corrplot" %in% class(fviz_pca_biplot)) {
        tiff(file, width = input$plotWidth * 100, height = input$plotHeight * 100, res = 100)
        fviz_pca_biplot(pcaResult(), repel = TRUE, labelsize = 5) + theme_classic(base_size = input$fontSize, base_family = input$fontFamily)
        dev.off()
      } else {
        g <- fviz_pca_biplot(pcaResult(), repel = TRUE, labelsize = 5) + theme_classic(base_size = input$fontSize, base_family = input$fontFamily)
        if (input$exportFormat == "png") {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
        } else {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
        }
      }
    }
  )

  output$download_contribVarPlot <- downloadHandler(
    filename = function() {
      paste0("contribVarPlot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      if ("corrplot" %in% class(fviz_contrib)) {
        tiff(file, width = input$plotWidth * 100, height = input$plotHeight * 100, res = 100)
        fviz_contrib(pcaResult(), choice = "var", axes = 1, top = input$topN) + theme_minimal(base_size = input$fontSize, base_family = input$fontFamily)
        dev.off()
      } else {
        g <- fviz_contrib(pcaResult(), choice = "var", axes = 1, top = input$topN) + theme_minimal(base_size = input$fontSize, base_family = input$fontFamily)
        if (input$exportFormat == "png") {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
        } else {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
        }
      }
    }
  )

  output$download_contribIndPlot <- downloadHandler(
    filename = function() {
      paste0("contribIndPlot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      if ("corrplot" %in% class(fviz_contrib)) {
        tiff(file, width = input$plotWidth * 100, height = input$plotHeight * 100, res = 100)
        fviz_contrib(pcaResult(), choice = "ind", axes = 1, top = input$topN) + theme_minimal(base_size = input$fontSize, base_family = input$fontFamily)
        dev.off()
      } else {
        g <- fviz_contrib(pcaResult(), choice = "ind", axes = 1, top = input$topN) + theme_minimal(base_size = input$fontSize, base_family = input$fontFamily)
        if (input$exportFormat == "png") {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
        } else {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
        }
      }
    }
  )

  output$download_cos2Plot <- downloadHandler(
    filename = function() {
      paste0("cos2Plot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      if ("corrplot" %in% class(corrplot)) {
        tiff(file, width = input$plotWidth * 100, height = input$plotHeight * 100, res = 100)
        corrplot(as.matrix(pcaResult()$var$cos2), is.corr = FALSE, method = "color", tl.col = "black", mar = c(1,1,1,1)); return(invisible())
        dev.off()
      } else {
        g <- corrplot(as.matrix(pcaResult()$var$cos2), is.corr = FALSE, method = "color", tl.col = "black", mar = c(1,1,1,1)); return(invisible())
        if (input$exportFormat == "png") {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
        } else {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
        }
      }
    }
  )

  output$download_hcpcPlot <- downloadHandler(
    filename = function() {
      paste0("hcpcPlot_", Sys.Date(), ".", input$exportFormat)
    },
    content = function(file) {
      if ("corrplot" %in% class(fviz_dend)) {
        tiff(file, width = input$plotWidth * 100, height = input$plotHeight * 100, res = 100)
        fviz_dend(hcpcResult(), rect = TRUE, show_labels = TRUE, cex = 1.2)
        dev.off()
      } else {
        g <- fviz_dend(hcpcResult(), rect = TRUE, show_labels = TRUE, cex = 1.2)
        if (input$exportFormat == "png") {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600)
        } else {
          ggsave(file, plot = g, width = input$plotWidth, height = input$plotHeight, dpi = 600, device = "tiff", compression = "lzw")
        }
      }
    }
  )

  output$downloadAllTables <- downloadHandler(
    filename = function() {
      paste0("PCA_Tables_", Sys.Date(), ".zip")
    },
    content = function(file) {
      tmpdir <- tempdir()
      eig <- as.data.frame(pcaResult()$eig)
      coords <- pcaResult()$ind$coord
      vars <- pcaResult()$var$coord
      write.csv(eig, file.path(tmpdir, "Eigenvalues.csv"))
      write.csv(coords, file.path(tmpdir, "Individuals_Coordinates.csv"))
      write.csv(vars, file.path(tmpdir, "Variables_Coordinates.csv"))
      zip::zipr(zipfile = file, files = list.files(tmpdir, full.names = TRUE))
    }
  )
}

shinyApp(ui, server)
