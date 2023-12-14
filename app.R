# Shiny app: multi-omics data analysis
# accepts either 1 or 2 input groups
# if 2 groups, performs differential expression analysis (limma::eBayes)
# also performs enrichment analyses, including drug mechanism enrichment
# analysis for each input data frame

library(shiny)
library(DT)
library(panSEA)
library(DMEA)
library(data.table)
library(dplyr)
library(msigdbr)
library(visNetwork)

# set limit for upload file size
MB.limit <- 180

# Define UI for application
ui <- fluidPage(
  # Application title
  titlePanel("panSEA: Multi-omic Set Enrichment Analysis"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      # select 1 or 2 groups
      radioButtons("Ngroups",
        "How many groups/conditions would you like to analyze?",
        choices = list(1, 2), selected = 1
      ),

      # keep track of last selection to add inputs if Ntypes is changed
      # source: https://stackoverflow.com/questions/40631788/shiny-observe-triggered-by-dynamicaly-generated-inputs?noredirect=1
      tags$script("$(document).on('change', '.dynamicSI select', function () {
                              Shiny.onInputChange('lastSelectId',this.id);
                              Shiny.onInputChange('lastSelect', Math.random());
                             });"),
      
      # select # of input data sets
      numericInput("Ntypes",
        "How many omics data sets would you like to process?",
        value = 1, min = 1, max = 6
      ),
      
      # allow input of multiple data sets
      uiOutput("omics"),

      # offer advanced settings
      checkboxInput(inputId = "advanced", label = "Optional: advanced settings", value = FALSE),
      conditionalPanel(
        condition = "input.advanced",
        textInput("g1Samples", "If 2 groups: columns containing samples belonging to group 1 (e.g., 2:4)"),
        textInput("g2Samples", "If 2 groups: columns containing samples belonging to group 2 (e.g., 5:6)"),
        selectInput(
          "GSEA.rank", "If 2 groups: rank for GSEA",
          c("Log2FC", "t"), "Log2FC"
        ),
        selectInput(
          "DMEA.rank", "Rank for DMEA",
          c("Pearson.est", "Spearman.est"), "Pearson.est"
        ),
        numericInput("n.min.per.set", "Minimum features or drugs per set", 6),
        numericInput("n.dot.sets", "Number of sets for dot plots", 10),
        numericInput("n.net.sets", "Number of sets for network graphs", 10),
        sliderInput("p.cutoff", "P-value cutoff",
          min = 0, max = 1, value = 0.05
        ),
        sliderInput("FDR.cutoff", "False discovery rate cutoff",
          min = 0, max = 1, value = 0.25
        ),
        selectInput("drugSensType",
          "Drug sensitivity data set",
          choices = c(
            "PRISM (adherent cancer cell lines)" = "prism",
            "Other" = "other"
          )
        ),
        conditionalPanel(
          condition = "input.drugSensType == 'other'",
          fileInput("drug.sens",
            paste0("Upload a CSV file with sample names in column 1 and drug names as the rest of the column names (", MB.limit, " MB limit)"),
            accept = ".csv"
          ),
          fileInput("drug.info",
            paste0("Upload a CSV file with a column 'Drug' containing drug names and 'moa' containing set annotations (", MB.limit, " MB limit)"),
            accept = ".csv"
          )
        )
      ),

      # include "Run" button
      actionButton(inputId = "run", label = "Run"),

      # including loading message
      conditionalPanel(
        condition = "input.run && !output.msg",
        textOutput(outputId = "fyi")
      ),

      # indicate when run is completed
      textOutput("msg")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      # display plots in tabs
      tabsetPanel(
        tabPanel(
          "Gene Set Enrichment Analysis: Dot Plot", plotOutput("GSEA.dot"),
          conditionalPanel(
            condition = "output.msg=='Run completed'",
            downloadButton(outputId = "GSEA.dot.dwnld", 
                           label = "Download dot plot of GSEA results")
          )
        ),
        tabPanel(
          "Drug Mechanism Enrichment Analysis: Dot Plot", 
          plotOutput("DMEA.dot"),
          conditionalPanel(
            condition = "output.msg=='Run completed'",
            downloadButton(outputId = "DMEA.dot.dwnld", 
                           label = "Download dot plot of DMEA results")
          )
        ),
        tabPanel(
          "Gene Set Enrichment Analysis: Interactive Network Graph", 
          visNetworkOutput("GSEA.int.net"),
          conditionalPanel(
            condition = "output.msg=='Run completed'",
            downloadButton(outputId = "GSEA.int.net.dwnld", 
                           label = "Download network graph of GSEA results")
          )
        ),
        tabPanel(
          "Drug Mechanism Enrichment Analysis: Interactive Network Graph", 
          visNetworkOutput("DMEA.int.net"),
          conditionalPanel(
            condition = "output.msg=='Run completed'",
            downloadButton(outputId = "DMEA.int.net.dwnld", 
                           label = "Download network graph of GSEA results")
          )
        ),
        tabPanel(
          "Gene Set Enrichment Analysis: Static Network Graph", 
          plotOutput("GSEA.net"),
          conditionalPanel(
            condition = "output.msg=='Run completed'",
            downloadButton(outputId = "GSEA.net.dwnld", 
                           label = "Download network graph of GSEA results")
          )
        ),
        tabPanel(
          "Drug Mechanism Enrichment Analysis: Static Network Graph", 
          plotOutput("DMEA.net"),
          conditionalPanel(
            condition = "output.msg=='Run completed'",
            downloadButton(outputId = "DMEA.net.dwnld", 
                           label = "Download network graph of GSEA results")
          )
        ),
        tabPanel(
          "Gene Set Enrichment Analysis: Averaged Results", 
          DT::dataTableOutput("GSEAResults"),
          conditionalPanel(
            condition = "output.msg=='Run completed'",
            downloadButton(outputId = "GSEAResults.dwnld", 
                           label = "Download averaged GSEA results")
          )),
        tabPanel(
          "Drug Mechanism Enrichment Analysis: Averaged Results", 
          DT::dataTableOutput("DMEAResults"),
          conditionalPanel(
            condition = "output.msg=='Run completed'",
            downloadButton(outputId = "DMEAResults.dwnld", 
                           label = "Download averaged DMEA results"))
        )
      ),

      # output result files (.zip)
      conditionalPanel(
        condition = "output.msg=='Run completed'",
        downloadButton(outputId = "results", label = "Download all results")
      ),
      uiOutput("info"),
      textOutput("private"),
      textOutput("refresh")
    )
  )
)

# Define backend
server <- function(input, output) {
  url <- a("https://github.com/belindabgarana/panSEA", href = "https://github.com/belindabgarana/panSEA")
  output$info <- renderUI({
    tagList("For more information or to contact us, please visit: ", url)
  })
  output$private <- renderText({
    "No user data is stored on our secure server, so your data will remain private."
  })
  output$refresh <- renderText({
    "Please refresh your web browser after each analysis and format your inputs to match the examples on the 'How to Use' page at the URL above to avoid errors. You will also need to refresh this webpage after 5 minutes of inactivity."
  })
  output$fyi <- renderText({
    "Running analysis... Please allow 1 to 5 minutes of run time"
  })
  options(shiny.maxRequestSize = MB.limit * 1024^2)
  
  # accept data inputs
  output$omics <- renderUI({
    # create list based on Ntypes
    buttons <- as.list(1:input$Ntypes)

    # create wellPanel for each input
    div( class = "dynamicSI",
         lapply(buttons, function(i) {
           wellPanel(
             helpText(paste0("Input omics data set #",i)),
             selectInput(paste0("type",i),
                         "Omics data type",
                         choices = c(
                           "Example (1 group): GSE66539 SKCM sensitive vs. resistant to vemurafenib",
                           "Example (1 group): GSE31625 NSCLC sensitive vs. resistant to erlotinib",
                           "Transcriptomics",
                           "Proteomics",
                           "Phospho-proteomics",
                           "Acetyl-proteomics",
                           "Metabolomics",
                           "Lipidomics"
                         )
             ),
             fileInput(paste0("data",i),
                       paste0("Upload a CSV file (", MB.limit, " MB limit)"),
                       accept = ".csv"
             ),
             textInput(paste0("featureName",i),
                       "Column name containing feature names", "Gene"),
              selectInput(paste0("gmtType",i),
                         "Feature set annotations",
                         choices = c(
                           "MSigDB (gene sets)",
                           "Kinase-substrate (substrate_site sets)",
                           "Other"
                         )
              ),
             conditionalPanel(
               condition = paste0("input.gmtType",i," == 'MSigDB (gene sets)'"),
               textInput(paste0("species",i), "Species", "Homo sapiens"),
               textInput(paste0("cat",i), "Category (e.g., H, C1, C2, ..., C8)", "C2"),
               textInput(paste0("subcat",i), "Optional: sub-category (e.g., CGP, CP:KEGG, or CP:REACTOME for category C2)", "CP:KEGG"),
               helpText("More info about MSigDB gene sets is available here: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp")
             ),
             conditionalPanel(
               condition = paste0("input.gmtType",i," == 'Kinase-substrate (substrate_site sets)'"),
               textInput(paste0("species",i), "Species", "human")
              ),
             conditionalPanel(
               condition = paste0("input.gmtType",i," == 'Other'"),
               fileInput(paste0("gmtInfo",i),
                         paste0("Upload a CSV file (", MB.limit, " MB limit)"),
                         accept = ".csv"
               ),
               textInput(paste0("feature",i), "Column name containing feature names (e.g., gene_symbol)"),
               textInput(paste0("set",i), "Column name containing set names (e.g., gs_name)"),
             ),
             selectInput(paste0("exprType",i),
                         "Expression data set to pair with drug sensitivity data",
                         choices = c(
                           "Adherent CCLE v19Q4" = "ccle",
                           "Other" = "other"
                         )
             ),
             conditionalPanel(
               condition = paste0("input.exprType",i," == 'other'"),
               fileInput(paste0("exprInfo",i),
                         paste0("Upload a CSV file with sample names in column 1 and feature names as the rest of the column names (", MB.limit, " MB limit)"),
                         accept = ".csv"
               )
             )
           )
         }
         )
        )
  })
  
  observeEvent(input$run, {
    # get input parameters
    cat(file = stderr(), "About to get n.groups", "\n")
    n.groups <- input$Ngroups
    if (n.groups == 1) {
      group.names <- "Diseased"
    } else {
      group.names <- c("Diseased", "Healthy")
    }
    n.types <- input$Ntypes
    GSEA.rank.var <- rep(input$GSEA.rank, n.types)
    DMEA.rank.var <- rep(input$DMEA.rank, n.types)
    p <- input$p.cutoff
    FDR <- input$FDR.cutoff
    min.per.set <- input$n.min.per.set
    drugSensType <- input$drugSensType
    cat(file = stderr(), "About to get drug.sensitivity", "\n")
    if (is.character(drugSensType)) {
      if (drugSensType == "prism") {
        drug.sensitivity <- "PRISM"
      }
    } else {
      drug.sensitivity <- input$drug.sens
    }
    n.network.sets <- input$n.net.sets
    n.dot.plot.sets <- input$n.dot.sets
    cat(file = stderr(), "About to get groupSamples", "\n")
    if (input$g1Samples != "" & input$g2Samples != "" &
        n.groups == 2) {
      groupSamples <- list(input$g1Samples, input$g2Samples)
    } else if (input$g1Samples != "" & input$g2Samples == "" &
               n.groups == 1) {
      groupSamples <- list(input$g1Samples)
    } else if (input$g2Samples != "" & input$g1Samples == "" &
               n.groups == 1) {
      groupSamples <- list(input$g2Samples)
    } else {
      groupSamples <- list(2)
    }
    cat(file = stderr(), "About to get gmt.drugs", "\n")
    if (input$drugSensType == "prism") {
      gmt.drugs <- "PRISM"
    } else {
      gmt.drugs <- DMEA::as_gmt(drug.info, min.per.set = min.per.set)
    }
    
    # get input data
    cat(file = stderr(), "About to get input", "\n")
    gmt.features <- list()
    data.list <- list()
    expression <- list()
    types <- c()
    feature.names <- c()
    for (i in 1:n.types){
      # reformat gmt.features
      cat(file = stderr(), "About to get info for MSigDB", "\n")
      if (input[[paste0("gmtType",i)]] == "MSigDB (gene sets)") {
        gmt.features[[i]] <- paste0(
          "msigdb_", input[[paste0("species",i)]], "_",
          input[[paste0("cat",i)]], "_", input[[paste0("subcat",i)]]
        )
        cat(file = stderr(), "About to make ksdb gmt", "\n")
      } else if (input[[paste0("gmtType",i)]] == "Kinase-substrate (substrate_site sets)") {
        ksdb <- read.csv(paste0("https://raw.githubusercontent.com/",
                                "BelindaBGarana/panSEA/shiny-app/data/",
                                "ksdb_20231101.csv"))
        ksdb <- ksdb[ksdb$KIN_ORGANISM == input[[paste0("species",i)]], ]
        ksdb$SUB_SITE <- paste0(ksdb$SUBSTRATE, ksdb$SUB_MOD_RSD, collapse = "_")
        gmt.features[[i]] <- DMEA::as_gmt(ksdb, "SUB_SITE", "KINASE", min.per.set = 6, 
                            descriptions = "KIN_ACC_ID")
        cat(file = stderr(), "About to make custom gmt", "\n")
      } else if (input[[paste0("gmtType",i)]] == "Other") {
        custom.gmt.features <- read.csv(input[[paste0("gmtInfo",i)]]$datapath) 
        gmt.features[[i]] <- 
          DMEA::as_gmt(custom.gmt.features, 
                       input[[paste0("feature",i)]], input[[paste0("set",i)]], 
                       min.per.set)
      }
      
      # read input files
      cat(file = stderr(), "About to check for GSE66539 example", "\n")
      if (input[[paste0("type",i)]] == "Example (1 group): GSE66539 SKCM sensitive vs. resistant to vemurafenib") {
        data.list[[i]] <- read.csv(
          paste0("https://raw.githubusercontent.com/BelindaBGarana/DMEA/",
                 "shiny-app/Examples/Gene_signature/",
                 "GSE66539_SKCM_sensitive_vs_resistant_to_vemurafenib",
                 "/Filtered_gene_signature_no_duplicates.csv"))
        cat(file = stderr(), "About to check for GSE31625 example", "\n")
      } else if (input[[paste0("type",i)]] == "Example (1 group): GSE31625 NSCLC sensitive vs. resistant to erlotinib") {
        data.list[[i]] <- read.csv(
          paste0("https://raw.githubusercontent.com/BelindaBGarana/DMEA/",
                 "shiny-app/Examples/Gene_signature/",
                 "GSE31625_NSCLC_sensitive_vs_resistant_to_erlotinib",
                 "/Filtered_gene_signature_no_duplicates.csv"))
        data.list[[i]][ , 3] <- NULL
        cat(file = stderr(), "About to get check for custom input", "\n")
      } else {
        data.list[[i]] <- read.csv(input[[paste0("data",i)]]$datapath)
      }
      
      # adjust groupSamples if necessary
      cat(file = stderr(), "About to adjust groupSamples if necessary", "\n")
      if (groupSamples == 2 & n.groups == 2) {
        groupSamples <- list(2:(0.5*(ncol(data.list[[1]])+1)), 
                             (0.5*(ncol(data.list[[1]])+1)+1):ncol(data.list[[1]]))
      }
      cat(file = stderr(), "About to get expression data", "\n")
      if (input[[paste0("exprType",i)]] == "ccle") {
        expression[[i]] <- "adherent CCLE"
      } else {
        expression[[i]] <- read.csv(input[[paste0("exprInfo",i)]]$datapath)
      }
  
      # compile other input parameters
      types <- c(input[[paste0("type",i)]], types)
      feature.names <- c(input[[paste0("featureName",i)]], feature.names)
    }

    # run panSEA
    cat(file = stderr(), "About to run panSEA", "\n")
    results <- panSEA::panSEA(
      data.list, types, feature.names,
      GSEA.rank.var, DMEA.rank.var, group.names,
      groupSamples, gmt.features, gmt.drugs, p,
      FDR, min.per.set, drug.sensitivity, expression,
      n.network.sets, n.dot.plot.sets
    )
    deg <- data.table::rbindlist(results$DEGs, use.names = TRUE, idcol = "type")

    # output results
    cat(file = stderr(), "About to output results", "\n")
    if (length(group.names) == 2) {
      output$results <- downloadHandler(
        filename = function() {
          paste0("panSEA_results_", Sys.Date(), ".zip")
        },
        content = function(file) {
          all.files <- c(
            "Differential_expression_results.csv",
            "GSEA_results.csv",
            "Mean_GSEA_results.csv",
            "GSEA_dot_plot.pdf",
            "GSEA_correlation_matrix.pdf",
            "GSEA_static_network_graph.pdf",
            "GSEA_interactive_network.graph.html",
            "DMEA_results.csv",
            "Mean_DMEA_results.csv",
            "DMEA_dot_plot.pdf",
            "DMEA_correlation_matrix.pdf",
            "DMEA_static_network_graph.pdf",
            "DMEA_interactive_network.graph.html",
            "All_panSEA_results.rds"
          )
          write.csv(deg, all.files[1], row.names = FALSE)
          write.csv(results$mGSEA.results[[1]]$compiled.results$results,
                    all.files[2],
                    row.names = FALSE
          )
          write.csv(results$mGSEA.results[[1]]$compiled.results$mean.results,
                    all.files[3],
                    row.names = FALSE
          )
          ggsave(all.files[4], 
                 results$mGSEA.results[[1]]$compiled.results$dot.plot,
                 device = "pdf"
          )
          ggsave(all.files[5], 
                 results$mGSEA.results[[1]]$compiled.results$corr.matrix,
                 device = "pdf"
          )
          results$mGSEA.network[[1]]$static
          dev.print(pdf, all.files[6])
          visNetwork::visSave(results$mGSEA.network[[1]]$interactive, 
                              all.files[7])
          write.csv(results$mDMEA.results[[1]]$compiled.results$results,
                    all.files[8],
                    row.names = FALSE
          )
          write.csv(results$mDMEA.results[[1]]$compiled.results$mean.results,
                    all.files[9],
                    row.names = FALSE
          )
          ggsave(all.files[10], 
                 results$mDMEA.results[[1]]$compiled.results$dot.plot,
                 device = "pdf"
          )
          ggsave(all.files[11], 
                 results$mDMEA.results[[1]]$compiled.results$corr.matrix,
                 device = "pdf"
          )
          results$mDMEA.network$static[[1]]
          dev.print(pdf, all.files[12])
          visNetwork::visSave(results$mDMEA.network[[1]]$interactive, 
                              all.files[13])
          saveRDS(results, all.files[14])
          zip(zipfile = file, files = all.files)
        }
      )

      output$GSEA.dot.dwnld <- downloadHandler(
        filename = function() {
          paste0("GSEA_dot_plot_", Sys.Date(), ".pdf")
        },
        content = function(file) {
          ggsave(file, results$mGSEA.results$compiled.results$dot.plot, 
                 device = "pdf")
        }
      )
      output$GSEA.dot <-
        renderPlot({
          results$mGSEA.results$compiled.results$dot.plot
        })

      output$DMEA.dot.dwnld <- downloadHandler(
        filename = function() {
          paste0("DMEA_dot_plot_", Sys.Date(), ".pdf")
        },
        content = function(file) {
          ggsave(file, results$mDMEA.results$compiled.results$dot.plot, 
                 device = "pdf")
        }
      )
      output$DMEA.dot <-
        renderPlot({
          results$mDMEA.results$compiled.results$dot.plot
        })
      
      output$GSEA.int.net.dwnld <- downloadHandler(
        filename = function() {
          paste0("GSEA_interactive_network_graph_", Sys.Date(), ".html")
        },
        content = function(file) {
          visSave(results$mGSEA.network$interactive, file)
        }
      )
      output$GSEA.int.net <-
        renderVisNetwork({
          results$mGSEA.network$interactive
        })
      
      output$DMEA.int.net.dwnld <- downloadHandler(
        filename = function() {
          paste0("DMEA_interactive_network_graph_", Sys.Date(), ".html")
        },
        content = function(file) {
          visSave(results$mDMEA.network$interactive, file)
        }
      )
      output$DMEA.int.net <-
        renderVisNetwork({
          results$mDMEA.network$interactive
        })
      

      output$GSEA.net.dwnld <- downloadHandler(
        filename = function() {
          paste0("GSEA_static_network_graph_", Sys.Date(), ".pdf")
        },
        content = function(file) {
          results$mGSEA.network$static
          dev.print(pdf, file)
        }
      )
      output$GSEA.net <-
        renderPlot({
          results$mGSEA.network$static
        })

      output$DMEA.net.dwnld <- downloadHandler(
        filename = function() {
          paste0("DMEA_static_network_graph_", Sys.Date(), ".pdf")
        },
        content = function(file) {
          results$mDMEA.network
          dev.print(pdf, file)
        }
      )
      output$DMEA.net <-
        renderPlot({
          results$mDMEA.network
        })

      output$GSEAResults.dwnld <- downloadHandler(
        filename = function() {
          paste0("GSEA_mean_results_", Sys.Date(), ".csv")
        },
        content = function(file) {
          write.csv(results$mGSEA.results$compiled.results$mean.results, file)
        }
      )
      output$GSEAResults <- 
        results$mGSEA.results$compiled.results$mean.results %>% 
        dplyr::mutate_if(is.numeric, ~round(., 3))
      
      output$DMEAResults.dwnld <- downloadHandler(
        filename = function() {
          paste0("DMEA_mean_results_", Sys.Date(), ".csv")
        },
        content = function(file) {
          write.csv(results$mDMEA.results$compiled.results$mean.results, file)
        }
      )
      output$DMEAResults <- 
        results$mDMEA.results$compiled.results$mean.results %>% 
        dplyr::mutate_if(is.numeric, ~round(., 3))
      
      output$msg <- renderText({
        "Run completed"
      })
    } else if (n.groups == 1) {
      output$results <- downloadHandler(
        filename = function() {
          paste0("panSEA_results_", Sys.Date(), ".zip")
        },
        content = function(file) {
          all.files <- c(
            "GSEA_results.csv",
            "Mean_GSEA_results.csv",
            "GSEA_dot_plot.pdf",
            "GSEA_correlation_matrix.pdf",
            "GSEA_static_network_graph.pdf",
            "GSEA_interactive_network.graph.html",
            "DMEA_results.csv",
            "Mean_DMEA_results.csv",
            "DMEA_dot_plot.pdf",
            "DMEA_correlation_matrix.pdf",
            "DMEA_static_network_graph.pdf",
            "DMEA_interactive_network.graph.html",
            "All_panSEA_results.rds"
          )
          write.csv(results$mGSEA.results[[1]]$compiled.results$results,
            all.files[1],
            row.names = FALSE
          )
          write.csv(results$mGSEA.results[[1]]$compiled.results$mean.results,
                    all.files[2],
                    row.names = FALSE
          )
          ggsave(all.files[3], 
                 results$mGSEA.results[[1]]$compiled.results$dot.plot,
            device = "pdf"
          )
          ggsave(all.files[4], 
                 results$mGSEA.results[[1]]$compiled.results$corr.matrix,
            device = "pdf"
          )
          results$mGSEA.network[[1]]$static
          dev.print(pdf, all.files[5])
          visNetwork::visSave(results$mGSEA.network[[1]]$interactive, 
                              all.files[6])
          write.csv(results$mDMEA.results[[1]]$compiled.results$results,
            all.files[7],
            row.names = FALSE
          )
          write.csv(results$mDMEA.results[[1]]$compiled.results$mean.results,
                    all.files[8],
                    row.names = FALSE
          )
          ggsave(all.files[9], 
                 results$mDMEA.results[[1]]$compiled.results$dot.plot,
            device = "pdf"
          )
          ggsave(all.files[10], 
                 results$mDMEA.results[[1]]$compiled.results$corr.matrix,
            device = "pdf"
          )
          results$mDMEA.network$static[[1]]
          dev.print(pdf, all.files[11])
          visNetwork::visSave(results$mDMEA.network[[1]]$interactive, 
                              all.files[12])
          saveRDS(results, all.files[13])
          zip(zipfile = file, files = all.files)
        }
      )

      output$GSEA.dot.dwnld <- downloadHandler(
        filename = function() {
          paste0("GSEA_dot_plot_", Sys.Date(), ".pdf")
        },
        content = function(file) {
          ggsave(file, results$mGSEA.results[[1]]$compiled.results$dot.plot, 
                 device = "pdf")
        }
      )
      output$GSEA.dot <-
        renderPlot({
          results$mGSEA.results[[1]]$compiled.results$dot.plot
        })

      output$DMEA.dot.dwnld <- downloadHandler(
        filename = function() {
          paste0("DMEA_dot_plot_", Sys.Date(), ".pdf")
        },
        content = function(file) {
          ggsave(file, results$mDMEA.results[[1]]$compiled.results$dot.plot, 
                 device = "pdf")
        }
      )
      output$DMEA.dot <-
        renderPlot({
          results$mDMEA.results[[1]]$compiled.results$dot.plot
        })
      
      output$GSEA.int.net.dwnld <- downloadHandler(
        filename = function() {
          paste0("GSEA_interactive_network_graph_", Sys.Date(), ".html")
        },
        content = function(file) {
          visSave(results$mGSEA.network[[1]]$interactive, file)
        }
      )
      output$GSEA.int.net <-
        renderVisNetwork({
          results$mGSEA.network[[1]]$interactive
        })

      output$DMEA.int.net.dwnld <- downloadHandler(
        filename = function() {
          paste0("DMEA_interactive_network_graph_", Sys.Date(), ".html")
        },
        content = function(file) {
          visSave(results$mDMEA.network[[1]]$interactive, file)
        }
      )
      output$DMEA.int.net <-
        renderVisNetwork({
          results$mDMEA.network[[1]]$interactive
        })
      
      
      output$GSEA.net.dwnld <- downloadHandler(
        filename = function() {
          paste0("GSEA_static_network_graph_", Sys.Date(), ".pdf")
        },
        content = function(file) {
          results$mGSEA.network[[1]]$static
          dev.print(pdf, file)
        }
      )
      output$GSEA.net <-
        renderPlot({
          results$mGSEA.network[[1]]$static
        })

      output$GSEAResults.dwnld <- downloadHandler(
        filename = function() {
          paste0("GSEA_mean_results_", Sys.Date(), ".csv")
        },
        content = function(file) {
          write.csv(results$mGSEA.results[[1]]$compiled.results$mean.results, 
                    file)
        }
      )
      output$GSEAResults <- 
        results$mGSEA.results[[1]]$compiled.results$mean.results %>% 
        dplyr::mutate_if(is.numeric, ~round(., 3))
      
      output$DMEAResults.dwnld <- downloadHandler(
        filename = function() {
          paste0("DMEA_mean_results_", Sys.Date(), ".csv")
        },
        content = function(file) {
          write.csv(results$mDMEA.results[[1]]$compiled.results$mean.results, 
                    file)
        }
      )
      output$DMEAResults <- 
        results$mDMEA.results[[1]]$compiled.results$mean.results %>% 
        dplyr::mutate_if(is.numeric, ~round(., 3))
      
      output$msg <- renderText({
        "Run completed"
      })
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
