
library(shiny)
library(FactoMineR)
library(factoextra)
library(DT)
library(ggplot2)
library(plotly)
library(psych)
library(corrplot)
library(car)
library(emmeans)
library(multcomp)
library(readxl)

ui <- navbarPage(
  title = "PCA & ANOVA Analysis ver.2.0 — Developed by Vietnamese Association at Kyushu University",
  
  # PCA Tab
  tabPanel("PCA Analysis",
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
        downloadButton("download_hcpcPlot", "Download HCPC Dendrogram"),
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
  ),
  
  # ANOVA Tab
  tabPanel("ANOVA Analysis",
    sidebarLayout(
      sidebarPanel(
        h4("Step 1: Upload Data"),
        fileInput("anovaFile", "Upload Data File", 
                  accept = c(".csv", ".xlsx", ".xls")),
        actionButton("useSampleData", "Use Sample Data"),
        hr(),
        
        h4("Step 2: Select ANOVA Model"),
        radioButtons("anovaType", "ANOVA Type:",
                    choices = c("One-way ANOVA" = "oneway",
                               "Two-way ANOVA" = "twoway"),
                    selected = "oneway"),
        
        uiOutput("anovaFactorSelect"),
        uiOutput("anovaResponseSelect"),
        
        actionButton("runANOVA", "Run ANOVA", class = "btn-primary"),
        hr(),
        
        h4("Plot Settings"),
        numericInput("anovaPlotWidth", "Plot width (inches)", 8),
        numericInput("anovaPlotHeight", "Plot height (inches)", 6),
        selectInput("anovaFontFamily", "Font family", 
                   choices = c("sans", "serif", "mono", "Times", "Arial", "Helvetica"),
                   selected = "sans"),
        numericInput("anovaFontSize", "Font size (pt)", 12, min = 8),
        
        radioButtons("anovaExportFormat", "Export Format", 
                    choices = c("png", "tiff"), inline = TRUE),
        downloadButton("downloadAnovaBoxplot", "Download Boxplot")
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel("Data Preview", 
                  p("Preview of uploaded data:"),
                  DTOutput("anovaDataPreview")),
          tabPanel("ANOVA Table", 
                  h4("Analysis of Variance Results"),
                  p("ANOVA table shows the significance of factors."),
                  verbatimTextOutput("anovaTableOutput"),
                  DTOutput("anovaTableFormatted")),
          tabPanel("Post-hoc Tests", 
                  h4("Post-hoc Comparisons"),
                  p("If ANOVA is significant (p < 0.05), post-hoc tests compare groups."),
                  uiOutput("posthocMethodSelect"),
                  verbatimTextOutput("posthocOutput"),
                  DTOutput("posthocTableFormatted")),
          tabPanel("Boxplot", 
                  h4("Group Comparison Boxplot"),
                  p("Boxplot in Nature journal style (white background, black borders)."),
                  plotOutput("anovaBoxplot", height = "600px"))
        )
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
  
  # ===== ANOVA Module =====
  
  # Reactive: Load ANOVA data
  anovaDataInput <- reactive({
    if (!is.null(input$anovaFile)) {
      ext <- tools::file_ext(input$anovaFile$name)
      if (ext == "csv") {
        return(read.csv(input$anovaFile$datapath, stringsAsFactors = TRUE))
      } else if (ext %in% c("xlsx", "xls")) {
        return(as.data.frame(read_excel(input$anovaFile$datapath)))
      }
    } else if (input$useSampleData > 0) {
      # Load sample data
      df <- read.csv("Simulated_Corn_Trait_Data.csv", stringsAsFactors = TRUE)
      # Add a categorical factor for demonstration
      df$Treatment <- factor(rep(c("T1", "T2", "T3", "T4"), length.out = nrow(df)))
      df$Location <- factor(rep(c("L1", "L2"), each = 2, length.out = nrow(df)))
      return(df)
    }
    return(NULL)
  })
  
  # UI: Select factors and response
  output$anovaFactorSelect <- renderUI({
    req(anovaDataInput())
    df <- anovaDataInput()
    factorVars <- names(df)[sapply(df, function(x) is.factor(x) || is.character(x))]
    
    if (input$anovaType == "oneway") {
      selectInput("factor1", "Select Factor:", choices = factorVars)
    } else {
      tagList(
        selectInput("factor1", "Select Factor 1:", choices = factorVars),
        selectInput("factor2", "Select Factor 2:", choices = factorVars)
      )
    }
  })
  
  output$anovaResponseSelect <- renderUI({
    req(anovaDataInput())
    df <- anovaDataInput()
    numericVars <- names(df)[sapply(df, is.numeric)]
    selectInput("response", "Select Response Variable:", choices = numericVars)
  })
  
  # Data preview
  output$anovaDataPreview <- renderDT({
    req(anovaDataInput())
    datatable(head(anovaDataInput(), 20), options = list(scrollX = TRUE))
  })
  
  # ANOVA result
  anovaResult <- eventReactive(input$runANOVA, {
    req(anovaDataInput(), input$factor1, input$response)
    
    df <- anovaDataInput()
    # Convert to factor if needed
    df[[input$factor1]] <- as.factor(df[[input$factor1]])
    
    if (input$anovaType == "oneway") {
      formula_str <- paste(input$response, "~", input$factor1)
      model <- aov(as.formula(formula_str), data = df)
    } else {
      req(input$factor2)
      df[[input$factor2]] <- as.factor(df[[input$factor2]])
      formula_str <- paste(input$response, "~", input$factor1, "*", input$factor2)
      model <- aov(as.formula(formula_str), data = df)
    }
    
    list(
      model = model,
      summary = summary(model),
      df = df,
      formula = formula_str
    )
  })
  
  # ANOVA table output
  output$anovaTableOutput <- renderPrint({
    req(anovaResult())
    print(anovaResult()$summary)
  })
  
  output$anovaTableFormatted <- renderDT({
    req(anovaResult())
    anova_table <- as.data.frame(anova(anovaResult()$model))
    anova_table <- round(anova_table, 4)
    datatable(anova_table, options = list(scrollX = TRUE))
  })
  
  # Post-hoc test method selection
  output$posthocMethodSelect <- renderUI({
    req(anovaResult())
    selectInput("posthocMethod", "Select Post-hoc Test:",
               choices = c("Tukey HSD" = "tukey",
                          "Bonferroni" = "bonferroni"),
               selected = "tukey")
  })
  
  # Post-hoc test results
  posthocResult <- reactive({
    req(anovaResult(), input$posthocMethod)
    
    model <- anovaResult()$model
    
    if (input$posthocMethod == "tukey") {
      # Tukey HSD
      emm <- emmeans(model, as.formula(paste("~", input$factor1)))
      result <- pairs(emm, adjust = "tukey")
    } else if (input$posthocMethod == "bonferroni") {
      # Bonferroni
      emm <- emmeans(model, as.formula(paste("~", input$factor1)))
      result <- pairs(emm, adjust = "bonferroni")
    }
    
    return(result)
  })
  
  output$posthocOutput <- renderPrint({
    req(posthocResult())
    print(posthocResult())
  })
  
  output$posthocTableFormatted <- renderDT({
    req(posthocResult())
    posthoc_df <- as.data.frame(posthocResult())
    posthoc_df <- round(posthoc_df, 4)
    datatable(posthoc_df, options = list(scrollX = TRUE))
  })
  
  # Boxplot - Nature style
  output$anovaBoxplot <- renderPlot({
    req(anovaResult())
    
    df <- anovaResult()$df
    
    if (input$anovaType == "oneway") {
      p <- ggplot(df, aes_string(x = input$factor1, y = input$response, fill = input$factor1)) +
        geom_boxplot(outlier.shape = 21, outlier.size = 2, color = "black", alpha = 0.7) +
        geom_jitter(width = 0.2, size = 2, alpha = 0.5, color = "black") +
        labs(x = input$factor1, y = input$response, title = paste("Boxplot:", input$response, "by", input$factor1)) +
        theme_classic(base_size = input$anovaFontSize, base_family = input$anovaFontFamily) +
        theme(
          panel.background = element_rect(fill = "white", color = "black", size = 1),
          panel.grid = element_blank(),
          axis.line = element_line(color = "black", size = 0.5),
          axis.ticks = element_line(color = "black"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
    } else {
      req(input$factor2)
      p <- ggplot(df, aes_string(x = input$factor1, y = input$response, fill = input$factor2)) +
        geom_boxplot(outlier.shape = 21, outlier.size = 2, color = "black", alpha = 0.7, position = position_dodge(0.8)) +
        geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
                  size = 1.5, alpha = 0.5, color = "black") +
        labs(x = input$factor1, y = input$response, fill = input$factor2,
             title = paste("Boxplot:", input$response, "by", input$factor1, "and", input$factor2)) +
        theme_classic(base_size = input$anovaFontSize, base_family = input$anovaFontFamily) +
        theme(
          panel.background = element_rect(fill = "white", color = "black", size = 1),
          panel.grid = element_blank(),
          axis.line = element_line(color = "black", size = 0.5),
          axis.ticks = element_line(color = "black"),
          legend.background = element_rect(fill = "white", color = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
    }
    
    print(p)
  })
  
  # Download boxplot
  output$downloadAnovaBoxplot <- downloadHandler(
    filename = function() {
      paste0("ANOVA_Boxplot_", Sys.Date(), ".", input$anovaExportFormat)
    },
    content = function(file) {
      df <- anovaResult()$df
      
      if (input$anovaType == "oneway") {
        p <- ggplot(df, aes_string(x = input$factor1, y = input$response, fill = input$factor1)) +
          geom_boxplot(outlier.shape = 21, outlier.size = 2, color = "black", alpha = 0.7) +
          geom_jitter(width = 0.2, size = 2, alpha = 0.5, color = "black") +
          labs(x = input$factor1, y = input$response, title = paste("Boxplot:", input$response, "by", input$factor1)) +
          theme_classic(base_size = input$anovaFontSize, base_family = input$anovaFontFamily) +
          theme(
            panel.background = element_rect(fill = "white", color = "black", size = 1),
            panel.grid = element_blank(),
            axis.line = element_line(color = "black", size = 0.5),
            axis.ticks = element_line(color = "black"),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold")
          )
      } else {
        p <- ggplot(df, aes_string(x = input$factor1, y = input$response, fill = input$factor2)) +
          geom_boxplot(outlier.shape = 21, outlier.size = 2, color = "black", alpha = 0.7, position = position_dodge(0.8)) +
          geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
                    size = 1.5, alpha = 0.5, color = "black") +
          labs(x = input$factor1, y = input$response, fill = input$factor2,
               title = paste("Boxplot:", input$response, "by", input$factor1, "and", input$factor2)) +
          theme_classic(base_size = input$anovaFontSize, base_family = input$anovaFontFamily) +
          theme(
            panel.background = element_rect(fill = "white", color = "black", size = 1),
            panel.grid = element_blank(),
            axis.line = element_line(color = "black", size = 0.5),
            axis.ticks = element_line(color = "black"),
            legend.background = element_rect(fill = "white", color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold")
          )
      }
      
      if (input$anovaExportFormat == "png") {
        ggsave(file, plot = p, width = input$anovaPlotWidth, height = input$anovaPlotHeight, dpi = 600)
      } else {
        ggsave(file, plot = p, width = input$anovaPlotWidth, height = input$anovaPlotHeight, 
               dpi = 600, device = "tiff", compression = "lzw")
      }
    }
  )
}

shinyApp(ui, server)
