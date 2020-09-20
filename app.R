#' CSV formatted file
#' 1st column should be samples and named "Group"
#' Please ensure the control group is at the top of table
#' Please ensure the sample names contain alphabetical characters,
#' and its alphabetical part in a group are consistent
#' @param pcr_rawfile is a file name of tidy PCR table
#' @param housekeeping_gene is a charactor name
#' One PCR run one input file
#' WeChat Official Account: HalfAppleTuring

library(shiny)
library(tidyverse)
library(DT)
library(shinythemes)

# Define UI ---------------------------------------------------------------
ui <- fluidPage(
  theme = "litera_bootswatch.scss",
  navbarPage(
    "Turing's Half Apple",
    tabPanel(
      "PCR Calc",

      # Application title
      titlePanel(
        h2("PCR Results: From CT to Fold Change", align = "center")
      ),
      sidebarLayout(
        sidebarPanel(
          fileInput(
            inputId = "pcr_rawfile",
            label = "Choose A CSV File",
            accept = c("text/csv", "text/comma-separated-values,text/plain")
          ),
          textInput(
            inputId = "housekeeping_gene",
            label = "Housekeeping Gene, CaSe SeNsItIvE",
            value = "",
            placeholder = "A Gene Name"
          ),
          actionButton("submitButton", label = "Submit"),
          downloadButton("downloadData", "Download")
        ), # sidebarPanel()
        mainPanel(
          h3("Showing Results"),
          dataTableOutput(
            outputId = "fc"
          ) # Result table
        ) # mainPanel()
      ) # sidebarLayout()
    ), # tabPanel(), PCR Calc
    tabPanel(
      "About",
      div(includeMarkdown("README.md"),
        align = "justify"
      )
    ) # tabPanel(), About
  ) # navbarPage()
) # fluidPage()

# Define server -----------------------------------------------------------
server <- function(input, output, session) {

  # Submit input
  submitButton <- eventReactive(input$submitButton, {
    input$pcr_rawfile
    input$housekeeping_gene
  })

  # Run calc
  run <- reactive({
    submitButton()

    # Environment -------------------------------------------------------------

    pcr_rawfile <- input$pcr_rawfile
    hk_gene <- input$housekeeping_gene
    pcr_raw <- pcr_rawfile$datapath %>%
      read_csv(col_names = TRUE) %>%
      group_by(Group)
    hk_gene_col <- pcr_raw %>%
      colnames() %>%
      str_which(hk_gene)
    colnames(pcr_raw)[hk_gene_col] <- "hk_gene"

    # Calc --------------------------------------------------------------------

    # Step 1
    ct_mean <- pcr_raw %>%
      summarise_at(
        .vars = names(.)[-1],
        tibble::lst(mean)
      )

    # Step 2
    d_ct <-
      {
        ct_mean %>% select(-c(Group, hk_gene_mean)) - ct_mean$hk_gene_mean
      } %>%
      mutate(
        Group = ct_mean$Group %>%
          str_match("[A-Za-z]*")
      ) %>%
      group_by(Group) %>%
      group_split()
    names(d_ct) <- ct_mean$Group %>%
      str_match("[A-Za-z]*") %>%
      as.character() %>%
      unique()

    # Step 3
    control_mean <- d_ct[[1]] %>%
      summarise_at(
        .vars = names(.)[1:(ncol(d_ct$CON) - 1)],
        tibble::lst(mean)
      )

    for (n_treat_group in 1:length(d_ct)) {
      # Multiple bio rep
      for (n_bio_rep in 2:nrow(d_ct[[n_treat_group]])) {
        if (exists("control_mean_rep")) {
          control_mean_rep <- control_mean_rep %>%
            add_row(control_mean)
        }
        if (!exists("control_mean_rep")) {
          control_mean_rep <- control_mean %>%
            add_row(control_mean)
        }
      }

      dd_ct_temp <- d_ct[[n_treat_group]][-{
        ncol(ct_mean) - 1
      }] - control_mean_rep
      row.names(dd_ct_temp) <- ct_mean$Group[which(ct_mean$Group %>%
        str_detect(
          names(d_ct)[n_treat_group]
        ))]
      dd_ct_temp <- list(dd_ct_temp)

      if (exists("dd_ct")) {
        dd_ct <- append(dd_ct, dd_ct_temp)
      }
      if (!exists("dd_ct")) {
        dd_ct <- dd_ct_temp
      }
      rm(control_mean_rep, dd_ct_temp)
    }

    names(dd_ct) <- ct_mean$Group %>%
      str_match("[A-Za-z]*") %>%
      as.character() %>%
      unique()

    # Step 4
    for (n_treat_group in 1:length(dd_ct)) {
      if (exists("fc_temp")) {
        fc <- rbind(
          fc_temp,
          2^-dd_ct[[n_treat_group]]
        )
        fc_temp <- fc
      }
      if (!exists("fc_temp")) {
        fc_temp <- 2^-dd_ct[[n_treat_group]]
      }
    }
    rm(fc_temp)

    colnames(fc) <- colnames(pcr_raw) %>%
      setdiff(c("Group", "hk_gene"))

    datatable(fc, rownames = TRUE)
    return(fc)
  })

  # Status/Output Text Box
  output$contents <- renderPrint({
    if (input$submitButton > 0) {
      isolate("Calculation complete.")
    } else {
      return("Server is ready for calculation.")
    }
  })

  # Prediction results table
  output$fc <- renderDataTable({
    if (input$submitButton > 0) {
      isolate(run())
    }
  })

  # Export data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("PCR Results_Fold Change_", input$pcr_rawfile)
    },
    content = function(file) {
      write.csv(
        run(),
        file,
        row.names = TRUE
      )
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
