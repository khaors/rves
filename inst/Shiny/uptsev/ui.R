#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(rves)
#
dbHeader <- dashboardHeader(title = "UPTSE-V",
                            tags$li(a(href = 'http://www.uptc.edu.co',
                                      icon("power-off"),
                                      title = "Back to Apps Home"),
                                    class = "dropdown"),
                            tags$li(a(href = 'http://www.uptc.edu.co',
                                      img(src = 'uptc_jpg.jpg',
                                          title = "Company Home", height = "30px"),
                                      style = "padding-top:10px; padding-bottom:10px;"),
                                    class = "dropdown"))



#
dbSidebar <- dashboardSidebar(
  # DEFINE SIDEBAR ITEMS
  #uiOutput("selectModelName"),
  #uiOutput("resettableFileInput"),
  #uiOutput("selectSlotName"),
  sidebarMenu(
    menuItem("Home", tabName = "home", icon = icon("home")),
    menuItem("Data", tabName = "data", icon = icon("table")),
    menuItem("Graphical Inversion", tabName = "manual", icon = icon("hand-spock-o")),
    menuItem("Automatic Inversion", tabName = "automatic", icon = icon("fighter-jet")),
    #menuItem("Reports", tabName = "reports", icon = icon("cogs")),
    menuItem("Model Diagnostic", tabName = "diagnostic", icon = icon("wrench")),
    menuItem("Source Code", icon = icon("code"), href = "https://github.com/khaors/rves")
  )
)
#########################################################################################
#                                    homeTab
#########################################################################################
homeTab <- tabItem(
  tabName = "home",
  h1("UPTSE-V: Vertical Electrical Sounding Analyst"),
  fluidRow(
    tags$em(
      "Copyright (c) 2017 Oscar Garcia-Cabrejo <khaors@gmail.com>",
      br(),
      br(),
        "Permission to use, copy, modify, and distribute this software for any
      purpose with or without fee is hereby granted, provided that the above
      copyright notice and this permission notice appear in all copies.",
      br(),
      br(),
      'THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES',
      "WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF",
      "MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR",
      "ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES",
      "WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN",
      "ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF",
      "OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE."
    ),
    br(),br(),
    box(
      h2("Information"),
      "This interface allows users to load, plot and analyze VES acquired using Schlumberger arrays.",
      "This tool is primarily meant to support the teaching of Hydrogeology at School of Geological ",
      "Engineering, UPTC - Sogamoso, Boyaca, Colombia. Although the main target is only for educational",
      "purposes, this tool can be used by other users to analyze their VES but the accuracy ",
      "and the geological interpretation derived from these results are not warranted",
      br(),br(),
      "The dashboard uses the following R libraries below and is being developed in RStudio. ",
      br(),
      tags$a(href="http://shiny.rstudio.com/", "shiny"),
      " and ",
      tags$a(href="https://rstudio.github.io/shinydashboard", "shinydashboard"),
      ". ",
      br(),br(),
      "The source code is available on ",
      tags$a(href="https://github.com/khaors/rves", "GitHub. "),
      "Contact the developer, Oscar Garcia-Cabrejo, via e-mail at khaors@gmail.com for
      questions or feedback."
    ),
    box(
      h2("Instructions"),
      "1. Load the results of the VES using the Data Tab on the right.",
      br(),br(),
      "2. Use the Graphical Inversion Tab to define a model using your 'Geophysical Skill'.",
      br(),br(),
      "3. Use the Automatic Inversion Tab to estimate the real resistivities and thicknesses using ",
      "a Nonlinear Least-Squares approach."
    )

  )
)
#########################################################################################
#                                    dataTab
#########################################################################################
dataTab <- tabItem(
  tabName = "data",
  h2("Import VES Data"),
  br(),
  br(),
  fluidRow(
    tags$em("This tab is to import data. VES data must be in a text file and you can use
            all the options below to make sure that you can import your data.")
  ),
  br(),
  br(),
  checkboxInput('header', ' Header?', TRUE),
  checkboxInput('rownames', ' Row names?', FALSE),
  selectInput('sep', 'Separator:',
              c("Comma","Semicolon","Tab","Space"), 'Comma'),
  selectInput('quote', 'Quote:',
              c("None","Double Quote","Single Quote"),
              'Double Quote'),
  selectInput('dec', 'Decimal mark', c("Period", "Comma"),
              'Period'),
  numericInput('nrow.preview','Number of rows in the preview:',20),
  numericInput('ncol.preview', 'Number of columns in the preview:',
               10),
  fileInput('file1', 'Choose CSV/TXT File'),
  helpText("Note: Even if the preview only shows a restricted
                          number of observations, the pumping_test object will be created
                          based on the full dataset."),
  tableOutput("view"),
  br()
  #downloadLink('downloadDump', 'Download source'),
  #br(),
  #downloadLink('downloadSave', 'Download binary')
)
#########################################################################################
#                                    Manual Inversion Tab
#########################################################################################
manualTab <- tabItem(
  tabName = "manual",
  h3("VES Inversion: Graphical Method"),
  br(),
  fluidRow(
    tags$em("On this tab you can specify an Earth model (real resisivities and thicknesses),
            calculate the apparent resisivities of this model, and compare them with the
            apparent resistivities measured on the field. The goal is to perform an inversion
            by hand/eye to help the students to get familiar with the interpretation of VES.
            Some fitness measures of the specified model are presented below the plot.")
  ),

  br(),
  br(),
  sidebarLayout(
    sidebarPanel(
      textInput(inputId = "manual_nlayers", label = "Number Layers", value = 1),
      textInput(inputId = "manual_res", label = "Real Resistivities", value = 10),
      textInput(inputId = "manual_thick", label = "Thicknesses", value = 10),
      uiOutput("manual_run")
    ),
    mainPanel(
      plotOutput(outputId="manual_plot"),
      br(),
      uiOutput("manual_results")

    )
  )
)
#########################################################################################
#                                    Automatic Inversion Tab
#########################################################################################
automaticTab <- tabItem(
  tabName = "automatic",
  h3("VES Inversion: Automatic Method"),
  br(),
  br(),
  fluidRow(
    tags$em("On this tab you can specify an initial Earth model (real resisivities and thicknesses),
             and using different optimization techniques an Earth model with minimum error is found.
            Some fitness measures of the estimated model are presented below the plot.")
  ),
  br(),
  br(),
  sidebarLayout(
    sidebarPanel(
      actionButton(inputId = "automatic_import", label = "Import Graphical Model"),
      textInput(inputId = "automatic_nlayers", label = "Number Layers", value = "1"),
      textInput(inputId = "automatic_res", label = "Initial Solution: Real Resistivities", value = "10"),
      textInput(inputId = "automatic_thick", label = "Initial Solution: Thicknesses", value = "10"),
      selectInput(inputId = "automatic_method", label = "Optimization Method: ",
                  choices = c("None", "Nonlinear Least Squares",
                              "L-BFGS-B",
                              "Simulated Annealing",
                              "Genetic Algorithms",
                              "Particle Swarm Optimization",
                              "Differential Evolution"),
                  selected = "None"),
      checkboxInput(inputId = "automatic_options1", label = "Specify Optimization Options?", value = FALSE),
      uiOutput("automatic_options2"),
      actionButton(inputId = "automatic_plot", label = "Plot"),
      actionButton(inputId="auto_run",label= "run, SEV run", icon = icon("bullseye"))
    ),
    mainPanel(
      plotOutput(outputId = "automatic_plot"),
      br(),
      h4(HTML("<b>Message Window</b>")),
      br(),
      textOutput(outputId = "automatic_msg"),
      br(),
      uiOutput(outputId = 'automatic_results'),
      br(),
      uiOutput(outputId = "automatic_table")
    )
  )
)
#########################################################################################
#                                    diagnosticTab
#########################################################################################
diagnosticTab <- tabItem(
  tabName = "diagnostic",
  h3("VES Inversion: Estimated Model Diagnostic"),
  br(),
  fluidRow(
    tags$em("On this tab you can check if the assumptions behind the estimation procedure
            hold or are violated. If these assumption hold then there is a greater confidence
            in the model results, otherwise the estimated parameters (real resistivities and
            thicknesses) must be used with care.")
  ),
  br(),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      selectInput(inputId = "diagnostic.type", label = "Select the Diagnostic Type",
                  choices = c("None", "Model Diagnostic"), selected = "None")
    ),
    mainPanel(
      # Add Plot
      plotOutput(outputId = "model_diagnostic") #, height = 500, width = 500*1.5)
    )
  )
)

#########################################################################################
#                                    reportsTab
#########################################################################################
reportsTab <- tabItem(
  tabName = "reports",
  h2("Reports"),
  "Notes:",
  br(),
  "1. This tab allows the automatic generation and download of some preconfigured graphs and data ",
  "tables given the selected RDF and slot on the sidebar menu.",
  br(),
  "2. Automatic generation of the reports for all the slots in a given RDF will be programmed at a ",
  "later time. ",
  br(),
  "3. This tab will also have options to generate a preconfigured set of graphs and analysis ",
  "typicaly produced for the USBR UC & LC operations and planning models. The idea is to ",
  "make the report generator smart enough to know from the RDF metadata whether to allow ",
  "the generation of 24MS, MTOM, or CRSS reports.",
  br(),br(),
  fluidRow(
    box(
      "Using the selected RDF and slot, ",
      "generate and download a zip file of the raw time-series, percentile envelope, ",
      "and percent exceedance graphs shown on the Graphs tab along with a csv file of the raw",
      "data.",
      br(),br(),
      "The selected graph options on the Graphs tab is used to generate the saved graphs. ",
      br(),br(),
      uiOutput("reportDownloadButton")
    )
  ),
  br(),
  div(
    style = 'overflow-x: scroll',
    DT::dataTableOutput("tableProbabilityData")
  )
)
#
userInterface <- dashboardPage(
  skin = "blue",
  # DASHBOARD HEADER
  dbHeader,
  # DASHBOARD SIDEBAR
  dbSidebar,
  # DASHBOARD BODY
  dashboardBody(
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "rdfTool.css")),
    tabItems(
      homeTab ,
      dataTab,
      manualTab,
      automaticTab, #,reportsTab
      diagnosticTab
    )
  )
)
