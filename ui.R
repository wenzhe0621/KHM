ui <-fluidPage(
  #shinythemes::themeSelector(),
  navbarPage(
    title = "LMDZ_SANDU/REF Select acceptable data",
    tabsetPanel(
    tabPanel("Page1",h4("This app is created for the LMDZ_SANDU/REF model calibration. Please see the ensemble plots and go to page 2 to select your acceptable runs."),
             h5("Page1 shows the observed field and 90 ensemble member plots."),  h5("Page2 is the Selection page where you need to choose your acceptable runs."), 
             h5("Page3 is used to do a final check and save your selection."),
             fluidRow(column(6,plotOutput('distPlot0',width = 500))),
             uiOutput("plots1")),
    tabPanel("Page2",sidebarPanel(  
      fluidRow(actionButton("goButton1", "Acceptable"),
      actionButton("goButton2", "Unacceptable")),
      br(),
      fluidRow(actionButton("goButton3", "Back"),actionButton("goButton4", "Next")),
      br(),
      fluidRow(textInput('input_label',"Select the ensemble member"),actionButton("goButton5", "Jump")),
    width = 3),mainPanel(
      column(6,plotOutput('distPlot00')),
      column(6,plotOutput('distPlot_choose')),width = 9,
      verbatimTextOutput("label_number"),
      verbatimTextOutput("label_number2")
    )),
    tabPanel("Page3",
             h3("You choose the following as acceptable"),
             fluidRow(actionButton("goButton6", "Refresh acceptable data"),
                      actionButton("goButton8", "Refresh unacceptable data")),
             br(),
             fluidRow(actionButton("goButton7", "Save your selection")),
             br(),
             fluidRow(column(width = 5,
                    box(title = "",
                        width = 12,
                        plotOutput(outputId="distPlot000")
                    ))),             br(),
             uiOutput("plots2"),
             h3("You choose the following as unacceptable"),
             uiOutput("plots3")
             )
    )
  )
)
