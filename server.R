server <- function(input, output, session) {
  
  
  clouldscol1<-rainbow(n=30,start=0.5001, end=1)
  clouldscol2<-rainbow(n=15,start=0.1, end=0.4)
  clouldscol<-c("darkgrey",clouldscol1,clouldscol2)
  
  
  position <- reactiveValues(label=1)
  
  Acceptable <<- c()
  Unacceptable <<- c()
  
  ###page1
  output$distPlot0 <- renderPlot({
  image.plot(1:clouds$timeN,clouds$zgrid, matrix(t(clouds$obs),ncol=50), ylim=c(0,max(clouds$zgrid )),col = clouldscol,
             xlab="Time",ylab="Level")
  title(main = "Observed field")
  })
  
  output$distPlot00 <- renderPlot({
    image.plot(1:clouds$timeN,clouds$zgrid, matrix(t(clouds$obs),ncol=50), ylim=c(0,max(clouds$zgrid )),col = clouldscol,
               xlab="Time",ylab="Level")
    title(main = "Observed field")
  })
  
  output$distPlot000 <- renderPlot({
    image.plot(1:clouds$timeN,clouds$zgrid, matrix(t(clouds$obs),ncol=50), ylim=c(0,max(clouds$zgrid )),col = clouldscol,
               xlab="Time",ylab="Level")
    title(main = "Observed field")
  })
  
  output$plots1 <- renderUI({
    # radioButtons(inputId = paste0("mVariable",i), label = paste0("mVariable",i), choices = c("A","B","C"))
    plot_output_list <- lapply(seq(1,90,4), function(i) {
      plotname1 <- paste("ui_plot_", i, sep="")
      plotname2 <- paste("ui_plot_", i+1, sep="")
      plotname3 <- paste("ui_plot_", i+2, sep="")
      plotname4 <- paste("ui_plot_", i+3, sep="")
      fluidRow(splitLayout(cellWidths = c("25%", "25%", "25%", "25%"),
                           plotOutput(plotname1), plotOutput(plotname2),plotOutput(plotname3),
                           plotOutput(plotname4)))
    })
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  
  for (i in 1:90) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_data = clouds$data[,i]
      plotname = sprintf('ui_plot_%i',i)
      titlename = paste0("Ensemble member ",i)
      output[[plotname]] = renderPlot({
        image.plot(1:clouds$timeN,clouds$zgrid, matrix(t(my_data),ncol=50), ylim=c(0,max(clouds$zgrid )),col = clouldscol,
                   xlab="Time",ylab="Level")
        title(main = titlename)})
    })
  }
  
  ###page2
  output$label_number <- renderText({
    loc<-position[['label']]
    if(loc < 90){ paste0("The data is ensemble member ",position[['label']],".")}
    else{ paste0("The data is ensemble member 90.
You have made all the selections. Please go to page 3, all of your selections will be presented.")}
    })
  
  output$label_number2 <- renderText({"Please see the ensemble plots and select your acceptable runs.
The acceptable runs cloud be the perfect runs which are consistent with observation. 
When there are no perfect runs, the data does not quite match with observation but still acceptable 
could be also seen as acceptable runs. 
In the extreme, without any good runs, the acceptable runs could be set at the best in all the 
training data."
  })
  
  output$distPlot_choose <- renderPlot({
    image.plot(1:clouds$timeN,clouds$zgrid, matrix(t(clouds$data[,position[['label']]]),ncol=50), ylim=c(0,max(clouds$zgrid )),col = clouldscol,
               xlab="Time",ylab="Level")
    title(main = "Ensemble member")
  })
  
  observeEvent(input$goButton1,{
    Acceptable <<- c(Acceptable,position[['label']])
    Acceptable <<- sort(Acceptable)
    Acceptable <<- unique(Acceptable)
    if (position[['label']] %in% Unacceptable) {
      Unacceptable <<- Unacceptable[-which(Unacceptable==position[['label']])]
    }
    position[['label']] <- position[['label']] + 1
    if( position[['label']]>90){ position[['label']]=90}
    })
  
  observeEvent(input$goButton2,{
    Unacceptable <<- c(Unacceptable,position[['label']])
    Unacceptable <<- sort(Unacceptable)
    Unacceptable <<- unique(Unacceptable)
    
    if (position[['label']] %in% Acceptable) {
      Acceptable <<- Acceptable[-which(Acceptable==position[['label']])]
    }
    position[['label']] <- position[['label']] + 1
    if( position[['label']]>90){ position[['label']]=90}
  })
  
  observeEvent(input$goButton3,{
    position[['label']] <- position[['label']] - 1
    if( position[['label']]>90){ position[['label']]=90}
  })
  
  observeEvent(input$goButton4,{
    position[['label']] <- position[['label']] + 1
    if( position[['label']]>90){ position[['label']]=90}
  })
  
  observeEvent(input$goButton5,{
    position[['label']] <- as.numeric(input$input_label)
    if( position[['label']]>90){ position[['label']]=90}
  })
  
  ####page3
  observeEvent(input$goButton7,{
    csvfile<-rep(0,90)
    csvfile[Acceptable]=1
    csvfile[Unacceptable]=2
    write.csv(csvfile,"Accepttable.csv")
  })
 
  
  output$plots2 <- renderUI({
    input$goButton6
    for (i in Acceptable) {
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderPlot() will be the same across all instances, because
      # of when the expression is evaluated.
      local({
        my_data = clouds$data[,i]
        plotname = sprintf('ui2_plot_%i',i)
        titlename = paste0("Data",i)
        output[[plotname]] = renderPlot({
          image.plot(1:clouds$timeN,clouds$zgrid, matrix(t(my_data),ncol=50), ylim=c(0,max(clouds$zgrid )),col = clouldscol,
                     xlab="Time",ylab="Level")
          title(main = titlename)})
      })
    }
    # radioButtons(inputId = paste0("mVariable",i), label = paste0("mVariable",i), choices = c("A","B","C"))
    plot_output_list <- lapply(seq(1,length(Acceptable),4), function(i) {
        plotname1 <- paste("ui2_plot_", Acceptable[i], sep="")
        plotname2 <- paste("ui2_plot_", Acceptable[i+1], sep="")
        plotname3 <- paste("ui2_plot_", Acceptable[i+2], sep="")
        plotname4 <- paste("ui2_plot_", Acceptable[i+3], sep="")
        fluidRow(splitLayout(cellWidths = c("25%", "25%", "25%", "25%"),
                             plotOutput(plotname1), plotOutput(plotname2),plotOutput(plotname3),
                             plotOutput(plotname4)))
    })
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  
  output$plots3 <- renderUI({
    input$goButton8
    for (i in Unacceptable) {
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderPlot() will be the same across all instances, because
      # of when the expression is evaluated.
      local({
        my_data = clouds$data[,i]
        plotname = sprintf('ui3_plot_%i',i)
        titlename = paste0("Data",i)
        output[[plotname]] = renderPlot({
          image.plot(1:clouds$timeN,clouds$zgrid, matrix(t(my_data),ncol=50), ylim=c(0,max(clouds$zgrid )),col = clouldscol,
                     xlab="Time",ylab="Level")
          title(main = titlename)})
      })
    }
    # radioButtons(inputId = paste0("mVariable",i), label = paste0("mVariable",i), choices = c("A","B","C"))
    plot_output_list <- lapply(seq(1,length(Unacceptable),4), function(i) {
      plotname1 <- paste("ui3_plot_", Unacceptable[i], sep="")
      plotname2 <- paste("ui3_plot_", Unacceptable[i+1], sep="")
      plotname3 <- paste("ui3_plot_", Unacceptable[i+2], sep="")
      plotname4 <- paste("ui3_plot_", Unacceptable[i+3], sep="")
      
      fluidRow(splitLayout(cellWidths = c("25%", "25%", "25%", "25%"),
                           plotOutput(plotname1), plotOutput(plotname2),plotOutput(plotname3),
                           plotOutput(plotname4)))
    })
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
}



