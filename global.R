##   This is a global file to run the app.
##   This app is created for the LMDZ model calibration,
##   we are aiming to consider the modellers information in the calibration process, 
##   please see the ensemble plots and select your acceptable runs.

#  Acceptable runs
#  The acceptable runs cloud be the training data which are consistent with observation. 
#  When there is no perfect runs, the data does not quite match with observation but still acceptable 
#  could be also seen as acceptable runs. In the extreme, without any good runs, 
#  the acceptable runs could be set at the best in all the training data.


#    Components for this file
##   server.r, ui.r  :     Shiny app files 
##   clouds.RData    :     SANDU/REF case 90 ensemble runs, observations
##   Selection file  :     Responses are saved to this file 
##                         A exeter selection csv file is include in this file (as an example). 

##   STEP 1 download the required packages
twd <- getwd()

if (!require("shiny")) install.packages("shiny")
if (!require("fields")) install.packages("fields")
#if (!require("shinyforms")) install.packages("shinyforms")
if (!require("shinydashboard")) install.packages("shinydashboard")

# if you have download the packages before, just run the following code.
packages <- c('shiny', 'fields', 'shinyforms', 'shinydashboard')
sapply(packages, require, character.only = TRUE, quietly = TRUE)

##   STEP 2 set the working directory,load the R file
# getwd()
#setwd("you filepath")
load("clouds_wave1_v2.RData")

##   STEP 3 run the Shiny app
# Attention:  The app shows two features.
#             The first one shows the plots for all ensembles and also the observation (for compare), 
#             you can check all the ensemble plots by draging the bar in page 1.
source('server.R')
source('ui.R')
shinyApp(ui,server)

##   STEP 4 Responses
# Responses are saved to the Selection file. 
# Noted: Please send us your response file by email.


