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
library(DT)
library(pracma)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  server.env <- environment() # used to allocate in functions
  current.table <- NULL
  current.ves.manual <- NULL
  current.ves.auto <- NULL
  original.ves <- NULL
  first <- FALSE
  filtered <- FALSE
  rmd.report.file <- NULL
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
    fname <- input$file1
    pos.valid <- !duplicated(d.input[,1])
    d.input <- d.input[pos.valid,]
    ves <- ves(fname, ab2 =d.input[,1],apprho= d.input[,2])
    server.env$current.ves <- ves
    server.env$original.ves <- ves
    #print(ves)
    head(d.input, n=input$nrow.preview)
  })
  ########################################################################################
  #                               Filter VES Tab
  ########################################################################################
  filter.ves <- function(){
    current.ves <- isolate(server.env$current.ves)
    if(is.null(current.ves))
      return(NULL)
    if(server.env$filtered){
      return(NULL)
    }
    filterMethod <- isolate(input$filterMethod)
    validate(
      need(filterMethod != "None", "Please select a filter method")
    )
    res <- NULL
    if(filterMethod == "None")
      return(NULL)
    else if(filterMethod == "smooth.spline"){
      res <- smoothing_ves(current.ves, method = filterMethod)
    }
    else if(filterMethod == "kernel.regression"){
      bw <- isolate(as.numeric(input$kernel_bw))
      res <- smoothing_ves(current.ves, method = filterMethod, bw = bw)
    }
    else if(filterMethod == "wavelet"){
      res <- smoothing_ves(current.ves, method = filterMethod)
    }
    #
    current.ves$ab2 <- res$ab2
    current.ves$appres <- res$apprho
    server.env$current.ves <- current.ves
    server.env$filtered <- TRUE
    return(current.ves)
  }
  #
  observeEvent(input$filterRun, {
    output$filterResultsPlot <- renderPlot({
      current.ves <- isolate(server.env$current.ves)
      validate(
        need(!is.null(current.ves), "The VES is not defined.  No filtering applied.")
      )
      if(is.null(current.ves))
        return(NULL)
      validate(
        need(!server.env$filtered, "The VES is already filtered")
      )
      if(server.env$filtered){
        return(NULL)
      }
      ab2.original <- isolate(current.ves$ab2)
      appres.original <- isolate(current.ves$appres)
      current.ves <- filter.ves()
      p1 <- plot(current.ves, type = "ves")
      original.df <- data.frame(ab2 = ab2.original, appres = appres.original)
      current.title <- paste0(current.ves$id, " Filtered (", input$filterMethod,")")
      p1 <- p1 + geom_point(aes(x = ab2, y = appres), data = original.df, show.legend = T) +
        ggtitle(current.title) #+
        #scale_colour_manual(name="Line.Color",
        #                    values=c(red="red", black="black"))
      print(p1)
    })
    #
    output$filterResultsTable <- renderDataTable({
      current.ves <- server.env$current.ves
      validate(
        need(!is.null(current.ves), "The VES is not defined")
      )
      if(is.null(current.ves))
        return(NULL)
      validate(
        need(server.env$filtered, "The VES has not been filtered")
      )
      #
      if(!server.env$filtered)
        return(NULL)
      ab2 <- current.ves$ab2
      filt.rho <- current.ves$appres
      real.rho <- server.env$original.ves$appres
      #print(filt.rho)
      res <- data.frame(spacing = ab2, filtered.rho = filt.rho,
                        real.rho = real.rho)
      return(res)
    },
    extensions = c('Buttons'),
    options = list(
      pageLength = length(server.env$current.ves$ab2),
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      text = 'Download',
      scrollY = 200,
      scroller = TRUE))
  })
  #
  observeEvent(input$filterRestore, ({
    server.env$current.ves <- server.env$original.ves
    server.env$filtered <- FALSE
    print("VES Restored.")
  }))
  ########################################################################################
  #                             Transformation Tab
  ########################################################################################
  observeEvent(input$transformationRun, {
    #
    output$transform_results <- renderDataTable({
      current.ves <- server.env$current.ves
      if(input$transform_results_plot){
        validate(
          need(!is.null(current.ves), "The VES is not defined")
        )
        #
        if(is.null(current.ves))
          return(NULL)
        #
        current.transformation <- isolate(input$transformation.type)
        validate(
          need(current.transformation != "None", "Please choose a Transformation")
        )
        current.transf <- NULL
        if(current.transformation == "Direct"){
          current.transf <- "direct"
        }
        else if(current.transformation == "Scaling"){
          current.transf <- "scaling"
        }
        else if(current.transformation == "Zohdy"){
          current.transf <- "zohdy"
        }
        else if(current.transformation == "Smoothed.Zohdy"){
          current.transf <- "smoothed_zohdy"
        }
        #
        base_transform <- "transform_"
        def_transform <- paste0(base_transform, current.transf)
        args <- list(ves = current.ves)
        res_transform <- do.call(def_transform, args)
        res <- data.frame(depth = res_transform$depth,
                          real.resistivity = res_transform$real.res)
        return(res)
      }
    }, extensions = c('Buttons'),
    options = list(
      pageLength = length(server.env$current.ves$ab2),
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      text = 'Download',
      scrollY = 200,
      scroller = TRUE))
    #
    output$transformSampleTable <- renderDataTable({
      current.ves <- server.env$current.ves
      if(input$transform_results_plot){
        validate(
          need(!is.null(current.ves), "The VES is not defined")
        )
        #
        if(is.null(current.ves))
          return(NULL)
        #
        if(!input$transform_sample_plot)
          return(NULL)
        #
        current.transformation <- isolate(input$transformation.type)
        validate(
          need(current.transformation != "None", "Please choose a Transformation")
        )
        current.transf <- NULL
        if(current.transformation == "Direct"){
          current.transf <- "direct"
        }
        else if(current.transformation == "Scaling"){
          current.transf <- "scaling"
        }
        else if(current.transformation == "Zohdy"){
          current.transf <- "zohdy"
        }
        else if(current.transformation == "Smoothed.Zohdy"){
          current.transf <- "smoothed_zohdy"
        }
        #
        base_transform <- "transform_"
        def_transform <- paste0(base_transform, current.transf)
        args <- list(ves = current.ves)
        res_transform <- do.call(def_transform, args)
        depth.mn <- min(res_transform$depth)
        depth.mx <- max(res_transform$depth)
        depth.seq <- seq(depth.mn, depth.mx, by = 1)
        #print(depth.seq)
        depth.int <- interp1(x = res_transform$depth,
                             y = res_transform$real.res,
                             xi = depth.seq,
                             method = "linear")
        #print(depth.int)
        #print(names(depth.int))
        res <- data.frame(depth = depth.seq,
                          real.resistivity = depth.int)
        return(res)
    }
  }, extensions = c('Buttons'),
     options = list(
     pageLength = length(server.env$current.ves$ab2),
     dom = 'Bfrtip',
     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
     text = 'Download',
     scrollY = 200,
     scroller = TRUE))
    #
    output$transformationPlot <- renderPlot({
      current.ves <- server.env$current.ves
      validate(
        need(!is.null(current.ves), "The VES is not defined")
      )
      if(is.null(current.ves))
        return(NULL)
      current.transformation <- isolate(input$transformation.type)
      validate(
        need(current.transformation != "None", "Please choose a Transformation")
      )
      current.transf <- NULL
      if(current.transformation == "Direct"){
        current.transf <- "direct"
      }
      else if(current.transformation == "Scaling"){
        current.transf <- "scaling"
      }
      else if(current.transformation == "Zohdy"){
        current.transf <- "zohdy"
      }
      else if(current.transformation == "Smoothed.Zohdy"){
        current.transf <- "zohdy.smoothed"
      }
      p1 <- plot(current.ves, type = "transformation", trans.type = current.transf)
      print(p1)
      return(p1)
    })
  })
  ########################################################################################
  #                           Graphical Inversion Tab
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
    validate(
      need(nlayers != 1, "Earth model with a single layer"),
      need(length(rho) == nlayers, "Resisitivities do not match number of layers"),
      need(length(thick) == nlayers, "Thicknesses do not match number of layers")
    )
    if(nlayers == 1 || length(rho) != nlayers || length(thick) != nlayers){
      current.ves.manual$interpreted <- FALSE
      server.env$current.ves.manual <- current.ves.manual
      return(current.ves.manual)
    }
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
      #validate(
      #  need(nlayers !=1, "Earth model with a single layer"),
      #  need(length(rho) == nlayers, "Resistivities do not match"),
      #  need(length(thick) == nlayers, "Thicknesses do not match")
      #)
      if(nlayers != length(rho)){
        return(NULL)
      }
      if(nlayers != length(thick)){
        return(NULL)
      }
      if(length(rho) != length(thick)){
        return(NULL)
      }
      p <- NULL
      #
      if(server.env$first){
        p <- plot(server.env$current.ves, type = "ves")
        server.env$first <- FALSE
      }
      else{
        tmp <- manual.model.results()
        p <- plot(tmp, type = "ves")
      }
      #
      output$manual_results <- renderUI({
        current.ves.manual <- server.env$current.ves.manual
        validate(
          need(!is.null(current.ves.manual), "The VES object is not defined.")
        )
        if(is.null(current.ves.manual))
          return(NULL)
        validate(
          need(input$manual_nlayers > 1, "Earth Model with a single layer")
        )
        if(isolate(input$manual_nlayers) == 1)
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
        rerr.num <- ((log10(cal.app.rho$appres)-log10(meas.app.rho))/log10(meas.app.rho))^2
        rel.err <- 100*sqrt(mean(rerr.num))
        mse <- mean((log10(cal.app.rho$appres)-log10(meas.app.rho))^2)
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
    updateTextInput(session, inputId = "automatic_nlayers", value = as.character(nlayers))
    updateTextInput(session, inputId = "automatic_res", value = as.character(rho))
    updateTextInput(session, inputId = "automatic_thick", value = as.character(thick))
  })
  #
  output$automatic_options2 <- renderUI({
    current.ves.auto <- server.env$current.ves
    tmp <- NULL
    if(is.null(current.ves.auto)){
      return(NULL)
    }
    if(input$automatic_options1){
      if(input$automatic_method == "Nonlinear Least Squares"){
        tmp <- wellPanel(
          h3("NLS Options"),
          br(),
          textInput(inputId = "nls_niter", label = "Number Iterations= ", value = 30),
          textInput(inputId = "nls_nreport", label = "Number Iter Report= ", value = 5)
        )
      }
      else if(input$automatic_method == "L-BFGS-B"){
        tmp <- wellPanel(
          h3("L-BFGS-B Options"),
          br(),
          textInput(inputId = "lbfgs_low_rho", label = "Min value rho= ",
                    value = "1.0"),
          textInput(inputId = "lbfgs_low_thick", label = "Min value thickness= ",
                    value = "1.0"),
          textInput(inputId = "lbfgs_up_rho", label = "Max Value rho= ", value = "1000"),
          textInput(inputId = "lbfgs_up_thick", label = "Max Value thickness= ",
                    value = "1000")
        )
      }
      else if(input$automatic_method == "Simulated Annealing"){
        tmp <- wellPanel(
          h3("Simulated Annealing Options"),
          br(),
          textInput(inputId = "sa_low_rho", label = "Min value rho= ",
                    value = "1.0"),
          textInput(inputId = "sa_low_thick", label = "Min value thicknness= ",
                    value = "1.0"),
          textInput(inputId = "sa_up_rho", label = "Max Value rho= ", value = "1000"),
          textInput(inputId = "sa_up_thick", label = "Max Value thickness= ",
                    value = "1000")
        )
      }
      else if(input$automatic_method == "Genetic Algorithms"){
        tmp <- wellPanel(
          h3("Genetic Algorithms Options"),
          br(),
          textInput(inputId = "ga_low_rho", label = "Min value rho= ",
                    value = "1.0"),
          textInput(inputId = "ga_low_thick", label = "Min value thickness= ",
                    value = "1.0"),
          textInput(inputId = "ga_up_rho", label = "Max Value rho= ", value = "1000"),
          textInput(inputId = "ga_up_thick", label = "Max Value thickness= ",
                    value = "1000")
        )
      }
      else if(input$automatic_method == "Particle Swarm Optimization"){
        tmp <- wellPanel(
          h3("Particle Swarm Optimization"),
          br(),
          textInput(inputId = "pso_low_rho", label = "Min value rho= ",
                    value = "1.0"),
          textInput(inputId = "pso_low_thick", label = "Min value thickness= ",
                    value = "1.0"),
          textInput(inputId = "pso_up_rho", label = "Max Value rho= ", value = "1000"),
          textInput(inputId = "pso_up_thick", label = "Max Value thickness= ",
                    value = "1000")
        )
      }
      else if(input$automatic_method == "Differential Evolution"){
        tmp <- wellPanel(
          h3("Differential Evolution"),
          br(),
          textInput(inputId = "de_low_rho", label = "Min value rho = ",
                    value = "1.0"),
          textInput(inputId = "de_low_thick", label = "Min value thickness = ",
                    value = "1.0"),
          textInput(inputId = "de_up_rho", label = "Max Value rho= ", value = "1000"),
          textInput(inputId = "de_up_thick", label = "Max Value thickness= ",
                    value = "1000")
        )
      }
    }
    return(tmp)
  })
  #
  calibrate.results <- function(){
    #
    current.ves.auto <- server.env$current.ves
    validate(
      need(!is.null(current.ves.auto), "VES has not been defined")
    )
    if(is.null(current.ves.auto)){
      return(NULL)
    }
    #
    nlayers <- isolate(as.numeric(input$automatic_nlayers))
    validate(
      need(nlayers > 1, "Earth Model with a single layer")
    )
    if(nlayers == 1){
      return(NULL)
    }
    #
    rho <- isolate(as.numeric(unlist(strsplit(input$automatic_res,","))))
    validate(
      need(length(rho) == nlayers & length(rho) > 1, "The number of resistivities
           does not match with the number of layers")
    )
    thick <- isolate(as.numeric(unlist(strsplit(input$automatic_thick,","))))
    validate(
      need(length(thick) == nlayers & length(thick) > 1, "The number of thicknesses
           does not match with the number of layers")    )

    automatic_method <- isolate(input$automatic_method)
    validate(
      need(automatic_method != "None", "Select a valid optimization method")
    )
    check_options <- isolate(input$automatic_options1)
    #
    if(length(rho) != nlayers | length(thick) != nlayers){
      output$automatic_msg <- renderText("Incorrect dimensions of rho or thick")
      return(NULL)
    }
    #
    nls_niter <- 100
    nls_nreport <- 10
    if(check_options){
      if(automatic_method == "Nonlinear Least Squares"){
        nls_niter <- isolate(as.numeric(input$nls_niter))
        nls_nreport <- isolate(as.numeric(nls_nreport))
      }
    }
    # Define Initial Solution
    par0 <- c(rho, thick)
    #print(par0)
    #print(class(par0))
    # Estimate model parameters
    if(automatic_method == "Nonlinear Least Squares"){
      current.res <- calibrate_nls(current.ves.auto, par0 = par0,
                                   iterations = nls_niter,
                                   ireport = nls_nreport,
                                   trace = FALSE)
    }
    else if(automatic_method == "L-BFGS-B"){
      lower_rho <- isolate(as.numeric(input$lbfgs_low_rho))
      lower_thick <- isolate(as.numeric(input$lbfgs_low_thick))
      upper_rho <- isolate(as.numeric(input$lbfgs_up_rho))
      upper_thick <- isolate(as.numeric(input$lbfgs_up_thick))
      lower <- c(rep(lower_rho, nlayers), rep(lower_thick, nlayers))
      upper <- c(rep(upper_rho, nlayers), rep(upper_thick, nlayers))
      #print(lower)
      #print(upper)
      print("L-BFGS-B: Working...")
      current.res <- calibrate(current.ves.auto, opt.method = "L-BFGS-B",
                               obj.fn = "log_rss",
                               par0 = par0,
                               lower = lower, upper = upper)
      print("L-BFGS-B: Finished")
    }
    else if(automatic_method == "Simulated Annealing"){
      lower_rho <- isolate(as.numeric(input$sa_low_rho))
      lower_thick <- isolate(as.numeric(input$sa_low_thick))
      upper_rho <- isolate(as.numeric(input$sa_up_rho))
      upper_thick <- isolate(as.numeric(input$sa_up_thick))
      lower <- c(rep(lower_rho, nlayers), rep(lower_thick, nlayers))
      upper <- c(rep(upper_rho, nlayers), rep(upper_thick, nlayers))
      #print(lower)
      #print(upper)
      print("Simulated Annealing: Working...")
      current.res <- calibrate(current.ves.auto, opt.method = "SA",
                               obj.fn = "log_rss",
                               par0 = par0,
                               lower = lower, upper = upper)
      print("Simulated Annealing: Finished")
    }
    else if(automatic_method == "Genetic Algorithms"){
      lower_rho <- isolate(as.numeric(input$ga_low_rho))
      lower_thick <- isolate(as.numeric(input$ga_low_thick))
      upper_rho <- isolate(as.numeric(input$ga_up_rho))
      upper_thick <- isolate(as.numeric(input$ga_up_thick))
      lower <- c(rep(lower_rho, nlayers), rep(lower_thick, nlayers))
      upper <- c(rep(upper_rho, nlayers), rep(upper_thick, nlayers))
      print("Genetic Algorithms: Working...")
      current.res <- calibrate(current.ves.auto, opt.method = "GA",
                               obj.fn = "log_rss",
                               par0 = par0,
                               lower = lower, upper = upper)
      print("Genetic Algorithms: Finished")
    }
    else if(automatic_method == "Particle Swarm Optimization"){
      lower_rho <- isolate(as.numeric(input$pso_low_rho))
      lower_thick <- isolate(as.numeric(input$pso_low_thick))
      upper_rho <- isolate(as.numeric(input$pso_up_rho))
      upper_thick <- isolate(as.numeric(input$pso_up_thick))
      lower <- c(rep(lower_rho, nlayers), rep(lower_thick, nlayers))
      upper <- c(rep(upper_rho, nlayers), rep(upper_thick, nlayers))
      print("PSO: Working...")
      current.res <- calibrate(current.ves.auto, opt.method = "PSO",
                               obj.fn = "log_rss",
                               par0 = par0,
                               lower = lower, upper = upper)
      print("PSO: Finished")
    }
    else if(automatic_method == "Differential Evolution"){
      lower_rho <- isolate(as.numeric(input$de_low_rho))
      lower_thick <- isolate(as.numeric(input$de_low_thick))
      upper_rho <- isolate(as.numeric(input$de_up_rho))
      upper_thick <- isolate(as.numeric(input$de_up_thick))
      lower <- c(rep(lower_rho, nlayers), rep(lower_thick, nlayers))
      upper <- c(rep(upper_rho, nlayers), rep(upper_thick, nlayers))
      print("Differential Evolution: Working...")
      current.res <- calibrate(current.ves.auto, opt.method = "DE",
                               obj.fn = "log_rss",
                               par0 = par0,
                               lower = lower, upper = upper)
      print("Differential Evolution: Finished")
    }
    updateSelectInput(session, inputId = "diagnostic.type", selected = "None")
    #print(names(current.res))
    return(current.res)
  }
  #
  observeEvent(input$automatic_plot,{
    output$automatic_plot <- renderPlot({
      current.ves <- server.env$current.ves
      validate(
        need(!is.null(current.ves), "The VES object is not defined")
      )
      current.ves$interpreted <- FALSE
      plot(current.ves, type = "ves")
    })
  })
  #
  observeEvent(input$auto_run, {
    output$automatic_plot <- renderPlot({
      current.ves <- server.env$current.ves
      validate(
        need(!is.null(current.ves), "The VES object is not defined")
      )
      nlayers <- isolate(as.numeric(input$automatic_nlayers))
      validate(
        need(nlayers > 1, "Earth model with a single layer")
      )
      rho <- isolate(as.numeric(unlist(strsplit(input$automatic_res,","))))
      thick <- isolate(as.numeric(unlist(strsplit(input$automatic_thick,","))))
      p <- NULL
      validate(
        need(length(rho) == nlayers & length(rho) > 1, "Number of resistivities
             does not match with number of layers")
      )
      if(length(rho) != nlayers)
        return(NULL)
      validate(
        need(length(thick) == nlayers & length(thick) > 1, "Number of thicknesses
             does not match with number of layers")
      )
      if(length(thick) != nlayers)
        return(NULL)
      if(is.null(current.ves))
        return(NULL)
      validate(
        need(nlayers > 1, "VES model with a single layer")
      )
      if(nlayers == 1){
        return(NULL)
      }
      #
      current.res <- calibrate.results()
      #print(current.res)
      #
      total.depth.model <- sum(current.res$thickness)
      total.depth.ves <- max(current.ves$ab2)
      if(total.depth.model < total.depth.ves){
        depth.corr <- total.depth.ves - total.depth.model
        current.res$thickness[nlayers] <- current.res$thickness[nlayers] + depth.corr
      }
      current.ves$rhopar <- current.res$rho
      current.ves$thickpar <- current.res$thickness
      current.ves$interpreted <- TRUE
      server.env$current.ves <- current.ves
      #
      output$automatic_results <- renderUI({
        current.ves <- server.env$current.ves
        validate(
          need(!is.null(current.ves), "The VES object is not defined")
        )
        if(is.null(current.ves))
          return(NULL)
        #print(names(current.ves.manual))
        rho <- isolate(current.ves$rhopar)
        thick <- isolate(current.ves$thickpar)
        spacing <- current.ves$ab2
        meas.app.rho <- isolate(current.ves$appres)
        automatic_method <- isolate(input$automatic_method)
        cal.app.rho <- apparent_resistivities(rho, thick, filt = rves::filt$V1,
                                              spacing = spacing)
        #print(cal.app.rho$appres)
        #print(meas.app.rho)
        rerr.num <- ((log10(cal.app.rho$appres)-log10(meas.app.rho))/log10(meas.app.rho))^2
        rel.err <- 100*sqrt(mean(rerr.num))
        mse <- mean((cal.app.rho$appres-meas.app.rho)^2)
        str1 <- "<h3>Results Parameter Estimation</h3><br>"
        str2 <- paste("<b>Relative Error(%)= </b>", format(rel.err, digits = 3), "<br>", sep = " ")
        str3 <- paste("<b>Mean Squared Error= </b>", format(mse, digits = 3), "<br>", sep = " ")
        str4 <- paste("<b>Optimization Method= </b>", automatic_method, "<br><br>", sep = " ")
        HTML(paste(str1, str2, str3, str4))
      })
      plot(current.ves, type = "ves")
    }) #renderPlot
    #
    output$automatic_table <- renderDataTable({
      current.ves <- server.env$current.ves
      validate(
        need(!is.null(current.ves), "The VES object is not defined")
      )
      #
      validate(
        need(current.ves$interpreted, "The VES object is not interpreted")
      )
      rho <- current.ves$rhopar
      thick <- current.ves$thickpar
      depth <- cumsum(current.ves$thickpar)
      res.df <- data.frame('Real_Resistivity(Ohm_m)' = rho, 'Thickness(m)' = thick,
                           'Depth(m)' = depth)
      nlay <- length(rho)
      layers <- vector('character', nlay)
      for(i in 1:nlay){
        layers[i] <- paste0("Layer", as.character(i))
      }
      row.names(res.df) <- layers
      res.df
    },
    extensions = c('Buttons'),
    options = list(
      pageLength = length(server.env$current.ves$ab2),
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      text = 'Download',
      scrollY = 200,
      scroller = TRUE))# renderDataTable

  }) #observeEvent
  ########################################################################################
  #                           Sequential Estimation Tab
  ########################################################################################
  calibrate.seq.results <- function(){
    #
    current.ves<- server.env$current.ves
    validate(
      need(!is.null(current.ves), "VES has not been defined")
    )
    if(is.null(current.ves)){
      return(NULL)
    }
    #
    niterations <- isolate(as.numeric(input$seqIterations))
    validate(
      need(niterations > 9, "Stepwise Estimation: Number of iterations too low")
    )
    if(niterations < 9){
      return(NULL)
    }
    #
    nreport <- isolate(as.numeric(input$seqReport))
    validate(
      need(nreport > 0 & nreport < niterations, "Stepwise Estimation: Number of report iterations wrong")
    )
    if(nreport == 1 & nreport == niterations){
      return(NULL)
    }
    #
    max.layers <- isolate(as.numeric(input$seqMaxlayers))
    validate(
      need(max.layers > 1, "Stepwise Estimation: Number of layers must be greater than 1")
    )
    if(max.layers == 1){
      return(NULL)
    }
    # Stepwise estimation using NLS
    #print(current.ves)
    current.optMethod <- isolate(input$seqOptMethod)
    current.res <- NULL
    if(current.optMethod == "None"){
      return(NULL)
    }
    else if(current.optMethod == "NLS"){
      current.res <- calibrate_step_nls(current.ves,
                                        iterations = niterations,
                                        ireport = nreport,
                                        max.layers = max.layers,
                                        trace = FALSE)
      current.res <- current.res$best.res.rss
    }
    else if(current.optMethod != "NLS"){
      lower.lim <- isolate(as.numeric(unlist(strsplit(input$seqLowerLim,",")))) #isolate(as.numeric(input$seqLowerLim))
      upper.lim <- isolate(as.numeric(unlist(strsplit(input$seqUpperLim,","))))#isolate(as.numeric(input$seqUpperLim))
      current.res1 <- calibrate_step(current.ves, opt.method = current.optMethod,
                                    max.layers = max.layers,
                                    lower = lower.lim,
                                    upper = upperlim)
      #print(names(current.res1))
      #print(current.res1$all.measures)
      current.res <- current.res1$best.res.rss
    }
    #
    current.ves$rhopar <- current.res$rho
    current.ves$thickpar <- current.res$thickness
    current.ves$interpreted <- TRUE
    #print(current.ves)
    #
    server.env$current.ves <- current.ves
    return(current.res)
  }
  #
  observeEvent(input$seq_plot,{
    output$seq_plot <- renderPlot({
      current.ves <- server.env$current.ves
      validate(
        need(!is.null(current.ves), "The VES object is not defined")
      )
      current.ves$interpreted <- FALSE
      plot(current.ves, type = "ves")
    })
  })
  #
  observeEvent(input$seqRun, {
    output$seq_plot <- renderPlot({
      current.ves <- server.env$current.ves
      validate(
        need(!is.null(current.ves), "The VES object is not defined")
      )
      #
      current.res <- calibrate.seq.results()
      #print(current.res)
      #current.res <- current.res1$best.res.rss
      nlayers <- length(current.res$rho)
      #print(nlayers)
      #
      total.depth.model <- sum(current.res$thickness)
      total.depth.ves <- max(current.ves$ab2)
      if(total.depth.model < total.depth.ves){
        depth.corr <- total.depth.ves - total.depth.model
        current.res$thickness[nlayers] <- current.res$thickness[nlayers] + depth.corr
      }
      current.ves$rhopar <- current.res$rho
      current.ves$thickpar <- current.res$thickness
      current.ves$interpreted <- TRUE
      server.env$current.ves <- current.ves
      #
      plot(current.ves, type = "ves")
    })
    #
    output$seq_results <- renderUI({
      current.ves <- server.env$current.ves
      validate(
        need(!is.null(current.ves), "The VES object is not defined")
      )
      if(is.null(current.ves))
        return(NULL)
      #print(names(current.ves.manual))
      rho <- isolate(current.ves$rhopar)
      thick <- isolate(current.ves$thickpar)
      spacing <- current.ves$ab2
      meas.app.rho <- isolate(current.ves$appres)
      automatic_method <- isolate(input$automatic_method)
      cal.app.rho <- apparent_resistivities(rho, thick, filt = rves::filt$V1,
                                            spacing = spacing)
      #print(cal.app.rho$appres)
      #print(meas.app.rho)
      rerr.num <- ((log10(cal.app.rho$appres)-log10(meas.app.rho))/log10(meas.app.rho))^2
      rel.err <- 100*sqrt(mean(rerr.num))
      mse <- mean((cal.app.rho$appres-meas.app.rho)^2)
      str1 <- "<h3>Results Parameter Estimation</h3><br>"
      str2 <- paste("<b>Relative Error(%)= </b>", format(rel.err, digits = 3), "<br>", sep = " ")
      str3 <- paste("<b>Mean Squared Error= </b>", format(mse, digits = 3), "<br>", sep = " ")
      opt.m <- paste(input$seqOptMethod, " Stepwise")
      str4 <- paste("<b>Optimization Method= </b>", opt.m, "<br><br>", sep = " ")
      HTML(paste(str1, str2, str3, str4))
    })
    #
    output$seq_table <- renderDataTable({
      current.ves <- server.env$current.ves
      validate(
        need(!is.null(current.ves), "The VES object is not defined")
      )
      #
      validate(
        need(current.ves$interpreted, "The VES object is not interpreted")
      )
      rho <- current.ves$rhopar
      nlay <- length(rho)
      thick <- current.ves$thickpar
      max.depth <- max(current.ves$ab2)/2.3
      thick[nlay] <- max.depth-thick[nlay-1]
      depth <- cumsum(thick)
      res.df <- data.frame('Real_Resistivity(Ohm_m)' = rho,
                           'Thickness(m)' = thick,
                           'Depth(m)' = depth)
      layers <- vector('character', nlay)
      for(i in 1:nlay){
        layers[i] <- paste0("Layer", as.character(i))
      }
      row.names(res.df) <- layers
      res.df
    },
    extensions = c('Buttons'),
    options = list(
      pageLength = length(server.env$current.ves$ab2),
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      text = 'Download',
      scrollY = 200,
      scroller = TRUE))# renderDataTable
  })


  ########################################################################################
  #                           Model Diagnostic Tab
  ########################################################################################
  output$model_diagnostic <- renderPlot({
    current.ves <- server.env$current.ves
    validate(
      need(!is.null(current.ves), "The VES object is not defined")
    )
    if(is.null(current.ves))
      return(NULL)
    validate(
      need(current.ves$interpreted, "An interpreted VES object is required")
    )
    if(!current.ves$interpreted)
      return(NULL)
    #print(input$diagnostic.type)
    validate(
      need(input$diagnostic.type != "None", "A valid diagnostic plot type is required")
    )
    if(input$diagnostic.type == "None")
      return(NULL)
    else if(input$diagnostic.type == "Model Diagnostic"){
      rho <- current.ves$rhopar
      thick <- current.ves$thickpar
      spacing <- current.ves$ab2
      meas.app.rho <- current.ves$appres
      cal.app.rho <- apparent_resistivities(rho, thick, filt = rves::filt$V1,
                                            spacing = spacing)
      #print(names(cal.app.rho))
      residuals.rho <- meas.app.rho-cal.app.rho$appres
      abs.residuals.rho <- sqrt(abs(residuals.rho))
      # Plot measured rho vs calculated rho
      df1 <- data.frame(measured = meas.app.rho, calculated = cal.app.rho$appres,
                        residuals = residuals.rho,
                        abs.residuals = abs.residuals.rho)
      p1 <- ggplot() + geom_point(aes(x=measured, y=calculated), data = df1,
                                  color = "red") +
        coord_equal() +
        ggtitle("a) Measured vs Calculated")
      #
      p2 <- ggplot(data = df1, aes(x = calculated, y = residuals)) +
        geom_point(color = "red") +
        geom_smooth() +
        ggtitle("b) Residuals")
      #
      p3 <- ggplot(data = df1, aes(x = calculated, y = abs.residuals)) +
        geom_point(color = "red") +
        geom_smooth() +
        ggtitle("c) Absolute Residuals")
      #
      p4 <- ggplot(data = df1, aes(sample = residuals)) + geom_qq(color = "red") +
        ggtitle("d) QQ plot")
      ptot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
      print(ptot)
    }
  })
  ########################################################################################
  #                           Report Generator Tab
  ########################################################################################
  #
  output$report.html.eng <- downloadHandler(
    filename = "report_eng.html",
    content = function(file) {
      tempReport <- file.path(tempdir(), "report_html_eng.Rmd")
      file.copy("report_html_eng.Rmd", tempReport, overwrite = TRUE)
      # Set up parameters to pass to Rmd document
      if(!is.null(server.env$current.ves)){
        params <- list(current.ves = server.env$current.ves)
      }
      else{
        #params <- list(current.ves = NULL)
        return(NULL)
      }
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  #
  output$report.html.spa <- downloadHandler(
    filename = "report_spa.html",
    content = function(file) {
      tempReport <- file.path(tempdir(), "report_html_spa.Rmd")
      file.copy("report_html_spa.Rmd", tempReport, overwrite = TRUE)
      # Set up parameters to pass to Rmd document
      if(!is.null(server.env$current.ves)){
        params <- list(current.ves = server.env$current.ves)
      }
      else{
        #params <- list(current.ves = NULL)
        return(NULL)
      }
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  #
  output$report.word.eng <- downloadHandler(
    filename = "report_eng.doc",
    content = function(file) {
      tempReport <- file.path(tempdir(), "report_word_eng.Rmd")
      file.copy("report_word_eng.Rmd", tempReport, overwrite = TRUE)
      # Set up parameters to pass to Rmd document
      if(!is.null(server.env$current.ves)){
        params <- list(current.ves = server.env$current.ves)
      }
      else{
        #params <- list(current.ves = NULL)
        return(NULL)
      }
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  #
  output$report.word.spa <- downloadHandler(
    filename = "report_spa.doc",
    content = function(file) {
      tempReport <- file.path(tempdir(), "report_word_spa.Rmd")
      file.copy("report_word_spa.Rmd", tempReport, overwrite = TRUE)
      # Set up parameters to pass to Rmd document
      if(!is.null(server.env$current.ves)){
        params <- list(current.ves = server.env$current.ves)
      }
      else{
        #params <- list(current.ves = NULL)
        return(NULL)
      }
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  #
  define_report_files1 <- reactive({
    current.format <- input$report.format
    current.language <- input$report.lang
    #print(current.language)
    current.language1 <- NULL
    if(current.format == "None"){
      return(NULL)
    }
    if(current.language == "None"){
      return(NULL)
    }
    else if(current.language == "English"){
      current.language1 <- "eng"
    }
    else{
      current.language1 <- "spa"
    }
    output.report.file <- paste0("report_", current.format, "_", current.language1, ".",
                                 current.format)
    print(output.report.file)
    rmd.report.file <- paste0("report_", current.format, "_", current.language1, ".Rmd")
    server.env$rmd.report.file <- rmd.report.file
    return(output.report.file)
  })
  #
  define_report_files <- function(){
    res <- define_report_files1()
    return(res)
  }
})
