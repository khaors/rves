#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(rves)
library(ggplot2)
library(gridExtra)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  server.env <- environment() # used to allocate in functions
  current.table <- NULL
  current.ves.manual <- NULL
  current.ves.auto <- NULL
  first <- FALSE
  # Panel 'Import data'
  dInput <- reactive({
    in.file <- input$file1
    #
    if (is.null(in.file))
      return(NULL)
    #
    the.sep <- switch(input$sep, "Comma"=",", "Semicolon"=";", "Tab"="\t",
                      "Space"="")
    the.quote <- switch(input$quote, "None"="","Double Quote"='"',
                        "Single Quote"="'")
    the.dec <- switch(input$dec, "Period"=".", "Comma"=",")
    if (input$rownames) {
      the.table <- read.table(in.file$datapath, header=input$header,
                              sep=the.sep, quote=the.quote, row.names=1,
                              dec=the.dec)
    } else {
      the.table <- read.table(in.file$datapath, header=input$header,
                              sep=the.sep, quote=the.quote, dec=the.dec)
      server.env$first <- TRUE
    }
    # return the table
    server.env$current.table <- the.table
    the.table
  })
  #
  # data preview table
  output$view <- renderTable({
    d.input <- dInput()
    if (is.null(d.input))
      return(NULL)
    if (ncol(d.input)>input$ncol.preview)
      d.input <- d.input[,1:input$ncol.preview]
    ves <- ves("VES1", ab2 =d.input[,1],apprho= d.input[,2])
    server.env$current.ves <- ves
    #print(ves)
    head(d.input, n=input$nrow.preview)
  })
  ########################################################################################
  #                           Manual Inversion Tab
  ########################################################################################
  output$manual_run <- renderUI({
    tmp <- actionButton(inputId = "manual_run1", label = "Press to Plot")
  })
  #
  manual.model.results <- function(){
    current.ves.manual <- server.env$current.ves
    if(is.null(current.ves.manual))
      return(NULL)
    nlayers <- isolate(as.numeric(input$manual_nlayers))
    rho <- isolate(as.numeric(unlist(strsplit(input$manual_res,","))))
    thick <- isolate(as.numeric(unlist(strsplit(input$manual_thick,","))))
    if(nlayers == 1)
      return(NULL)
    if(length(rho) != nlayers)
      return(NULL)
    if(length(thick) != nlayers)
      return(NULL)
    current.ves.manual$rhopar <- rho
    current.ves.manual$thickpar <- thick
    current.ves.manual$interpreted <- TRUE
    server.env$current.ves.manual <- current.ves.manual
    return(current.ves.manual)
  }
  #
  observeEvent(input$manual_run1,{
    output$manual_plot <- renderPlot({
      nlayers <- isolate(as.numeric(input$manual_nlayers))
      rho <- isolate(as.numeric(unlist(strsplit(input$manual_res,","))))
      thick <- isolate(as.numeric(unlist(strsplit(input$manual_thick,","))))
      p <- NULL
      #
      if(server.env$first){
        p <- plot(server.env$current.ves)
        server.env$first <- FALSE
      }
      else{
        tmp <- manual.model.results()
        p <- plot(tmp)
      }
      #
      output$manual_results <- renderUI({
        current.ves.manual <- server.env$current.ves.manual
        if(is.null(current.ves.manual))
          return(NULL)
        #print(names(current.ves.manual))
        rho <- current.ves.manual$rhopar
        thick <- current.ves.manual$thickpar
        spacing <- current.ves.manual$ab2
        meas.app.rho <- current.ves.manual$appres
        cal.app.rho <- apparent_resistivities(rho, thick, filt = rves::filt$V1,
                                              spacing = spacing)
        #print(cal.app.rho$appres)
        #print(meas.app.rho)
        rel.err <- 100*mean(abs(cal.app.rho$appres-meas.app.rho)/meas.app.rho)
        mse <- mean((cal.app.rho$appres-meas.app.rho)^2)
        str1 <- "<h3>Results Parameter Estimation</h3><br>"
        str2 <- paste("<b>Relative Error(%)= </b>", format(rel.err, digits = 3), "<br>", sep = " ")
        str3 <- paste("<b>Mean Squared Error= </b>", format(mse, digits = 3), sep = " ")
        HTML(paste(str1, str2, str3))
      })
      #
      return(p)
    })
  })
  #
  #
  ########################################################################################
  #                           Automatic Inversion Tab
  ########################################################################################
  # output$video <- renderUI({
  #   link <- "I5Z9WtTBZ_w"
  #     HTML(paste0('<iframe width="600" height="400" src="https://www.youtube.com/embed/',
  #                 link ,'" frameborder="0" allowfullscreen></iframe>'))
  # })
  #
  observeEvent(input$automatic_import, {
    nlayers <- as.numeric(input$manual_nlayers)
    rho <- as.numeric(unlist(strsplit(input$manual_res,",")))
    thick <- as.numeric(unlist(strsplit(input$manual_thick,",")))
    print(nlayers)
    print(rho)
    print(thick)
    updateTextInput(session, inputId = "automatic_nlayers", value = as.character(nlayers))
    updateTextInput(session, inputId = "automatic_res", value = as.character(rho))
    updateTextInput(session, inputId = "automatic_thick", value = as.character(thick))
  })
  #
  output$automatic_options2 <- renderUI({
    current.ves.auto <- server.env$current.ves
    tmp <- NULL
    if(is.null(current.ves.auto)){
      #output$automatic_msg <- "A VES has not been defined"
      return(NULL)
    }
    if(input$automatic_options1){
      if(input$automatic_method == "Nonlinear Least Squares"){
        tmp <- wellPanel(
          h3("NLS Options"),
          br(),
          textInput(inputId = "nls_niter", label = "Number Iterations= ", value = 100),
          textInput(inputId = "nls_nreport", label = "Number Iter Report= ", value = 10)
        )
      }
      else if(input$automatic_method == "Simmulated Annealing"){
        tmp <- wellPanel(
          h3("Simmulated Annealing Options"),
          br(),
          textInput(inputId = "sa_niter", label = "Max. Number Iterations", value = 100),
          checkboxInput(inputId = "sa_smooth", label = "Smooth Function?", value = TRUE),
          checkboxInput(inputId = "sa_simplefn", label = "Simple Function?", value = FALSE)
        )
      }
      else if(input$automatic_method == "Genetic Algorithms"){
        tmp <- wellPanel(
          h3("Genetic Algorithms Options"),
          br(),
          textInput(inputId = "ga_niter", label = "Number Iterations", value = 100)
        )
      }
    }
    return(tmp)
  })
  #
  calibrate.results <- function(){
    #
    current.ves.auto <- server.env$current.ves
    if(is.null(current.ves.auto)){
      output$automatic_msg <- renderText("VES has not been defined")
      return(NULL)
    }
    #
    nlayers <- isolate(as.numeric(input$automatic_nlayers))
    if(nlayers == 1){
      output$automatic_msg <- renderText("Number of layers in Earth model is 1")
      return(NULL)
    }
    #
    rho <- isolate(as.numeric(unlist(strsplit(input$automatic_res,","))))
    thick <- isolate(as.numeric(unlist(strsplit(input$automatic_thick,","))))
    #
    if(length(rho) != nlayers | length(thick) != nlayers){
      output$automatic_msg <- renderText("Incorrect dimensions of rho or thick")
      return(NULL)
    }
    #
    check_options <- isolate(input$automatic_options1)
    nls_niter <- 100
    nls_nreport <- 10
    if(check_options){
      if(input$automatic_method == "Nonlinear Least Squares"){
        nls_niter <- isolate(as.numeric(input$nls_niter))
        nls_nreport <- isolate(as.numeric(nls_nreport))
      }
    }
    # Define Initial Solution
    par0 <- c(rho, thick)
    #print(par0)
    #print(class(par0))
    # Estimate model parameters
    output$automatic_msg <- renderText({"Working on estimation..."})
    current.res <- calibrate_nls(current.ves.auto, par0 = par0,
                                 iterations = nls_niter,
                                 ireport = nls_nreport)
    output$automatic_msg <- renderText({"Working on estimation...Finished"})
    #print(names(current.res))
    return(current.res)
  }
  #
  observeEvent(input$automatic_plot,{
    output$automatic_plot <- renderPlot({
      current.ves <- server.env$current.ves
      current.ves$interpreted <- FALSE
      plot(current.ves)
    })
  })
  #
  observeEvent(input$auto_run, {
    output$automatic_plot <- renderPlot({
      current.ves <- server.env$current.ves
      nlayers <- isolate(as.numeric(input$automatic_nlayers))
      rho <- isolate(as.numeric(unlist(strsplit(input$automatic_res,","))))
      thick <- isolate(as.numeric(unlist(strsplit(input$automatic_thick,","))))
      p <- NULL
      if(length(rho) != nlayers)
        return(NULL)
      if(length(thick) != nlayers)
        return(NULL)
      if(is.null(current.ves))
        return(NULL)
      if(nlayers == 1){
        output$automatic_msg <- renderText("Initial Solution with only 1 layer")
        return(NULL)
      }
      #
      current.res <- calibrate.results()
      #print(names(current.res))
      #print(current.res$rho)
      #print(current.res$thickness)
      #
      current.ves$rhopar <- current.res$rho
      current.ves$thickpar <- current.res$thickness
      current.ves$interpreted <- TRUE
      server.env$current.ves <- current.ves
      #
      output$automatic_results <- renderUI({
        current.ves <- server.env$current.ves
        if(is.null(current.ves))
          return(NULL)
        #print(names(current.ves.manual))
        rho <- current.ves$rhopar
        thick <- current.ves$thickpar
        spacing <- current.ves$ab2
        meas.app.rho <- current.ves$appres
        cal.app.rho <- apparent_resistivities(rho, thick, filt = rves::filt$V1,
                                              spacing = spacing)
        #print(cal.app.rho$appres)
        #print(meas.app.rho)
        rel.err <- 100*mean(abs(cal.app.rho$appres-meas.app.rho)/meas.app.rho)
        mse <- mean((cal.app.rho$appres-meas.app.rho)^2)
        str1 <- "<h3>Results Parameter Estimation</h3><br>"
        str2 <- paste("<b>Relative Error(%)= </b>", format(rel.err, digits = 3), "<br>", sep = " ")
        str3 <- paste("<b>Mean Squared Error= </b>", format(mse, digits = 3), "<br><br>", sep = " ")
        HTML(paste(str1, str2, str3)) #, str4, str5))
      })
      plot(current.ves)
    }) #renderPlot
    #
    output$automatic_table <- renderTable({
      current.ves <- server.env$current.ves
      rho <- current.ves$rhopar
      thick <- current.ves$thickpar
      res.df <- data.frame('Real_Resistivity(Ohm_m)' = rho, 'Thickness(m)' = thick)
      nlay <- length(rho)
      layers <- vector('character', nlay)
      for(i in 1:nlay){
        layers[i] <- paste0("Layer", as.character(i))
      }
      row.names(res.df) <- layers
      res.df
    }, rownames = TRUE, digits = 3)# renderTable

  }) #observeEvent
  ########################################################################################
  #                           Model Diagnostic Tab
  ########################################################################################
  output$model_diagnostic <- renderPlot({
    current.ves <- server.env$current.ves
    if(is.null(current.ves))
      return(NULL)
    if(!current.ves$interpreted)
      return(NULL)
    print(input$diagnostic.type)
    if(input$diagnostic.type == "None")
      return(NULL)
    else if(input$diagnostic.type == "Model Diagnostic"){
      rho <- current.ves$rhopar
      thick <- current.ves$thickpar
      spacing <- current.ves$ab2
      meas.app.rho <- current.ves$appres
      cal.app.rho <- apparent_resistivities(rho, thick, filt = rves::filt$V1,
                                            spacing = spacing)
      print(names(cal.app.rho))
      residuals.rho <- meas.app.rho-cal.app.rho$appres
      abs.residuals.rho <- sqrt(abs(residuals.rho))
      # Plot measured rho vs calculated rho
      df1 <- data.frame(measured = meas.app.rho, calculated = cal.app.rho$appres,
                        residuals = residuals.rho,
                        abs.residuals = abs.residuals.rho)
      p1 <- ggplot() + geom_point(aes(x=measured, y=calculated), data = df1) +
        ggtitle("a) Measured vs Calculated")
      #
      p2 <- ggplot(data = df1, aes(x = calculated, y = residuals)) + geom_point() +
        geom_smooth() +
        ggtitle("b) Residuals")
      #
      p3 <- ggplot(data = df1, aes(x = calculated, y = abs.residuals)) + geom_point() +
        geom_smooth() +
        ggtitle("c) Absolute Residuals")
      #
      p4 <- ggplot(data = df1, aes(sample = residuals)) + geom_qq() +
        ggtitle("d) QQ plot")
      ptot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
      print(ptot)
    }
  })

})
