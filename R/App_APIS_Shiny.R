#' Shiny App for interactive session of APIS
#'
#' Launch the shiny interface to use APIS interactively
#'
#' @import dplyr
#' @import ggplot2
#' @import cowplot
#' @import shinyBS
#' @import shinythemes
#' @import htmltools
#' @rawNamespace import(shiny, except=c(dataTableOutput, renderDataTable))
#' @importFrom plotly plotlyOutput renderPlotly ggplotly
#' @importFrom parallel detectCores
#' @importFrom utils read.table
#' @importFrom DT datatable dataTableOutput renderDataTable
#' @importFrom gridExtra grid.table
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot.new
#' @importFrom utils head write.table read.table
#' @importFrom rlang .data
#'
#' @return void : most results are automatically saved
#'
#' @export

launch_APIShiny = function(){
  addResourcePath("www",system.file("www",package = "APIS"))
  ui <- fluidPage(
    theme = shinytheme("spacelab"),
    # Application title
    fluidRow(column(12,titlePanel(title = div("Parentage assignment with APIS",style="font-family:Arial;font-weight:bold;margin-bottom:-1em",
                                              img(src="/www/sysaaf.png",height=70,width=70),
                                              img(src="/www/inrae.png",height=50,width=120),
                                              img(src="/www/Ifremer.png",height=70,width=120),
                                              img(src="/www/feamp2.png",height=70,width=110))))),
    fluidRow(column(12,titlePanel(title = div("By J.Roche",style="color:#000000;font-size:10px;height:20px;")))),
    fluidRow(
      column(6,
             h2("Choose your step",style="font-size:20px;font-weight:bold"),
             navlistPanel(widths = c(3,9),id = 'nav',selected = 0,
                          tabPanel("Formatting",value=0,
                                   fileInput(inputId = "to_format",
                                             label = div("File of genotype to format for APIS (.ped from AXAS or .vcf)",
                                                         bsButton(inputId = "qx1",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                             accept = c(".txt",".ped",".vcf",".csv")),
                                   conditionalPanel(condition = "output.colFormat",
                                                    selectInput(inputId = 'header_file_format',
                                                                label = div("Is there a header in the file to format ?",
                                                                            bsButton(inputId = "q1110",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                choices = c('Yes','No'),selected = 'No'),
                                                    uiOutput(outputId = 'uiNoColSN'),
                                                    uiOutput(outputId = 'uiNoColGeno'),
                                                    selectInput(inputId = "what",
                                                                label = div("Which individuals are in the file ?",
                                                                            bsButton(inputId = "qx2",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                choices = c("Only parents","Only offspring","Both"),selected = "Only parents",multiple = FALSE)),
                                   conditionalPanel(condition = "input.what=='Both' | input.what=='Only offspring'",
                                                    selectInput(inputId = "ploidy_format",
                                                                label = div("Choose the ploidy level of the offspring",
                                                                            bsButton(inputId = "q52",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                choices = c(2,3),selected = 2)),
                                   conditionalPanel(condition = "input.what=='Both'",
                                                    fileInput(inputId = "list_par",
                                                              label = div("File with the name of the Parents",
                                                                          bsButton(inputId = "qx3",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                              accept = c(".txt",".csv")),
                                                    selectInput(inputId = 'header_list_par',
                                                                label = div("Is there a header in the file with the parent names ?",
                                                                            bsButton(inputId = "q111",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                choices = c('Yes','No'),selected = 'No'),
                                                    uiOutput(outputId = 'col_head0')),
                                   fileInput(inputId = "snp_map",
                                             label = div("File with the marker names (ex : .map from AXAS) (optional)",
                                                         bsButton(inputId = "qx4",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                             accept = c(".txt",".map",".csv")),
                                   conditionalPanel(condition = "output.colMap",
                                                    selectInput(inputId = 'header_snp_map',
                                                                label = div("Is there a header in the file with the marker names ?",
                                                                            bsButton(inputId = "q1111",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                choices = c('Yes','No'),selected = 'No'),
                                                    uiOutput(outputId = 'uiNoColMap')),
                                   selectInput(inputId = "markerType",
                                               label = div("Choose the type of marker",
                                                           bsButton(inputId = "q53",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                               choices = c("SNP","microsat"),selected = "SNP"),
                                   uiOutput(outputId = "go_format")
                          ),
                          tabPanel("APIS",value = 1,
                                   conditionalPanel(condition = "output.APIS_launched == 0",
                                                    checkboxInput(inputId = "both_parents_in",value = TRUE,
                                                                  label = div("Male and female parents are in the same dataset",
                                                                              bsButton(inputId = "q131",label = "",icon = icon("question"), style = "info", size = "extra-small"))),
                                                    selectInput(inputId = "ploidy",
                                                                label = div("Choose the ploidy level of the offspring",
                                                                            bsButton(inputId = "q5",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                choices = c(2,3),selected = 2),
                                                    fileInput(inputId = "data_off",
                                                              label = div("Dataset of OFFSPRING",
                                                                          bsButton(inputId = "q1",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                              accept = c(".Rdata",".txt")),
                                                    conditionalPanel(condition = "input.both_parents_in==true",
                                                                     fileInput(inputId = "data_par",
                                                                               label = div("Dataset of PARENTS",
                                                                                           bsButton(inputId = "q2",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                               accept = c(".Rdata",".txt")),
                                                                     fileInput(inputId = "sexe_par",
                                                                               label = div("File (.txt) with the sex of each parents",
                                                                                           bsButton(inputId = "q3",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                               accept = c(".txt",".csv")),
                                                                     conditionalPanel(condition = "output.dta_sex_load",
                                                                                      selectInput(inputId = 'header_sexe_par',
                                                                                                  label = div("Is there a header in the file with the marker names ?",
                                                                                                              bsButton(inputId = "q11111",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                                                  choices = c('Yes','No'),selected = 'Yes'),
                                                                                      uiOutput(outputId = "uiChangeSN"),
                                                                                      # conditionalPanel(condition = "output.dta_sex_load",
                                                                                      uiOutput(outputId = "uiChangeSe"))),
                                                    conditionalPanel(condition = "input.both_parents_in!=true",
                                                                     fileInput(inputId = "data_par1",
                                                                               label = div("Dataset of male (sire)",
                                                                                           bsButton(inputId = "q2.1",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                               accept = c(".Rdata",".txt")),
                                                                     fileInput(inputId = "data_par2",
                                                                               label = div("Dataset of female (dam)",
                                                                                           bsButton(inputId = "q2.2",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                               accept = c(".Rdata",".txt"))),
                                                    fileInput(inputId = "snp_in",
                                                              label = div("File (.txt) with the markers to use (optional)",
                                                                          bsButton(inputId = "q4",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                              accept = ".txt"),
                                                    uiOutput(outputId = "nb_snp_poss"),
                                                    selectInput(inputId = "markerType2",
                                                                label = div("Choose the type of marker",
                                                                            bsButton(inputId = "q54",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                choices = c("SNP","microsat"),selected = "SNP"),
                                                    sliderInput(inputId = "core",
                                                                label = div("Choose number of cores for parallelization",
                                                                            bsButton(inputId = "q7",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                min = 1,max = detectCores() - 1,value = detectCores() - 2,step = 1),
                                                    uiOutput(outputId = "go_apis")),
                                   conditionalPanel(condition = "output.APIS_launched == 1",
                                                    selectInput(inputId = "method",
                                                                label = div("Choose the method for assignment",
                                                                            bsButton(inputId = "q6",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                choices = c("exclusion","likelihood"),selected = "likelihood"),
                                                    conditionalPanel(condition = "input.method=='exclusion'",
                                                                     radioButtons(inputId = "exclu_thres",
                                                                                  label = div("Choose threshold method for exclusion",
                                                                                              bsButton(inputId = "q61",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                                  choices = c("Mismatch number","Error rate"),
                                                                                  selected = "Error rate")
                                                    ),
                                                    conditionalPanel(condition = "input.method=='likelihood' || (input.method=='exclusion' && input.exclu_thres=='Error rate')",
                                                                     sliderInput(inputId = "acceptError",label = "Error rate allowed",min = 0,max = 0.25,value = 0.05,step = 0.005)
                                                    ),
                                                    conditionalPanel(condition = "input.method=='exclusion' && input.exclu_thres=='Mismatch number'",
                                                                     uiOutput(outputId = "acceptMismatch0")
                                                                     # sliderInput(inputId = "acceptMismatch",label = "Number of mismatch allowed",min = 0,max = 50,value = 5,step = 1)
                                                    ),
                                                    # sliderInput(inputId = "acceptError",label = "Error rate allowed",min = 0,max = 0.25,value = 0.05,step = 0.01),
                                                    textInput(inputId = "save_name",
                                                              label = div("Name of the file to save (click the '?' for help)",
                                                                          bsButton(inputId = "q8",label = "",icon = icon("question"), style = "info", size = "extra-small"))),
                                                    actionButton(inputId = "SaveAPIS",label = "Save pedigree & log")
                                   )
                          ),
                          tabPanel("Verification",value = 2,
                                   fileInput(inputId = "data_res",
                                             label = div("Dataset created after the assignment with APIS in the first part (.Rdata)",
                                                         bsButton(inputId = "q9",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                             accept = ".Rdata"),
                                   fileInput(inputId = "tab_accou",
                                             label = div("The file with the mating plan",
                                                         bsButton(inputId = "q10",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                             accept = c(".txt",".csv")),
                                   conditionalPanel(condition = "output.dta_fac_load",
                                                    selectInput(inputId = 'header_tab_accou',
                                                                label = div("Is there a header in the file ?",
                                                                            bsButton(inputId = "q11100",label = "",icon = icon("question"), style = "info", size = "extra-small")),
                                                                choices = c('Yes','No'),selected = 'No'),
                                   ),
                                   conditionalPanel(condition = "output.dta_fac_load",
                                                    uiOutput(outputId = "uiChangeSN2")),
                                   conditionalPanel(condition = "output.dta_fac_load",
                                                    uiOutput(outputId = "uiChangeFa")),
                                   textInput(inputId = "save_name2",
                                             label = div("Name of the file to save (optional)",
                                                         bsButton(inputId = "q18",label = "",icon = icon("question"), style = "info", size = "extra-small"))),
                                   actionButton(inputId = "launch_verif",label = "Launch the verification"),
                                   conditionalPanel(condition = "input.launch_verif>0",
                                                    actionButton(inputId = "SavePlot",label = "Save output"))
                                   
                          )
             ),
      ),
      column(5,
             h5(div("Warnings : Make sure your app is launched on the desired directory. If not, quit the app and use setwd() !",
                    style="color:red;font-weight:bold;margin-bottom:3em")),
             conditionalPanel(condition = "input.nav==0",
                              h5(div("Warnings : this step is only required if you did not use the genotyping shiny application AND works only for diploid genotype. Else the format is already good OR you should use the genotyping application for your triploids.",
                                     style="color:red;font-weight:bold;margin-bottom:3em")),
                              a("Input .ped format example",href="/www/Ped_Example.txt",target="_blank"), # for hyperlink
                              a("Input .vcf format example",href="/www/Vcf_Example.txt",target="_blank"), # for hyperlink
                              h6("You can find 'log_Formating.txt' in './log' directory for more info."),
                              # a("log Formating",href="/www/log/log_Formating.txt",target="_blank"), # for hyperlink
                              dataTableOutput(outputId = "head"),
                              dataTableOutput(outputId = "headLpar"),
                              dataTableOutput(outputId = "headMap"),
                              textOutput(outputId = "end_format")),
             conditionalPanel(condition = "input.nav==1",
                              a("Input .Rdata or .txt format example",href="/www/Data_Example.txt",target="_blank"), # for hyperlink
                              h6("You can find 'log_APIS.txt' in './log' directory for more info."),
                              # a("log APIS",href="/www/log/log_APIS.txt",target="_blank"), # for hyperlink
                              dataTableOutput(outputId = "head_sex"),
                              textOutput(outputId = "WarningMarker"),
                              tags$head(tags$style("#WarningMarker{color: red;font-size: 20px;}")),
                              conditionalPanel(condition = "output.APIS_launched == 1",
                                               h4("Plot result from APIS : Select accepted error or mismatch threshold"),
                                               plotOutput(outputId = "graph_res"))),
             conditionalPanel(condition = "input.nav==2",
                              dataTableOutput(outputId = "head_fac"),
                              navbarPage("Choose section",
                                         tabPanel(title = "Summary",
                                                  textOutput(outputId = "txt1"),
                                                  textOutput(outputId = "txt2"),
                                                  textOutput(outputId = "txt3"),
                                                  textOutput(outputId = "txt4"),
                                                  dataTableOutput(outputId = "tab"),
                                                  plotOutput(outputId = "plot")),
                                         tabPanel(title = "Parents fertility",
                                                  textOutput(outputId = "sent1.1"),
                                                  textOutput(outputId = "sent1.2"),
                                                  textOutput(outputId = "sent2.1"),
                                                  textOutput(outputId = "sent2.2"),
                                                  plotlyOutput(outputId = "barplot")),
                                         tabPanel(title = "Mating plan",
                                                  plotlyOutput(outputId = "heatmap"))))
      ),
      column(1)
    ),
    bsTooltip(id = "qx1",title = "File with in row, individuals and in column, genotype. Each allele have to be separated by a space if .txt<br/>Else, file created by AXAS with the extension .ped<br/>.vcf format also accepted.<br/>See examples for more details"),
    bsTooltip(id = "qx2",title = "Select whether there is only parents/offspring or if there is both in the .ped file<br/>If it is both, a file with parents names will be asked to separate parents from offspring"),
    bsTooltip(id = "qx3",title = "The .txt file with the parents names must be an unique column as follow (no header) :<br/>Indiv1<br/>Indiv2<br/>Indiv3<br/>..."),
    bsTooltip(id = "qx4",title = "File with MarkerName in the same order as the genotype in the input file<br/>Can be .txt file with all marker names as a single column (no header) or the file created by AXAS (for example) with the markers names with the extension .map<br/>It is used to have the names of the markers (no provided in .ped file)"),
    bsTooltip(id = "q1",title = "File created by the genotyping application (or during formating phase) with extension _genoAPIS.Rdata<br/>Can also be a .txt file with marker as columns and individuals as rows with genotype format as A/A/A or A/A.<br/>This is the offspring dataset"),
    bsTooltip(id = "q2",title = "File created by the genotyping application (or during formating phase) with extension _genoAPIS.Rdata<br/>Can also be a .txt file with marker as columns and individuals as rows with genotype format as A/A.<br/>This is the parents dataset"),
    bsTooltip(id = "q3",title = "A file containing at least two rows : SampleName (or CodeBarre) & Sex like so :<br/>SampleName  Sex<br/>Name1  Sex1<br/>Name2  Sex2<br/>Name3  Sex3<br/>...  ...<br/>The SampleName variable is the name of the sample as in the genotyping application.<br/>The Sex variable should take values from 1 to 3 : 1, male ; 2, female ; 3, neo-male ; 4, neo-female<br/>It is used to separate male from female"),
    bsTooltip(id = "q4",title = "A file containing one row with NO header like so :<br/>Marker1<br/>Marker2<br/>Marker3<br/>...  ...<br/>Only those markers will be used ofr assignment at the condition that they are in the dataset."),
    bsTooltip(id = "q5",title = "Select 2 if offspring are diploids and 3 if they are triploids"),
    bsTooltip(id = "q131",title = "Tick if male (sire) and female (dam) parents are in the same dataset"),
    bsTooltip(id = "q2.1",title = "File created by the genotyping application (or during formating phase) with extension _genoAPIS.Rdata<br/>File with only the male (sire) parents.<br/>Can also be a .txt file with marker as columns and individuals as rows with genotype format as A/A."),
    bsTooltip(id = "q2.2",title = "File created by the genotyping application (or during formating phase) with extension _genoAPIS.Rdata<br/>File with only the female (dam) parents.<br/>Can also be a .txt file with marker as columns and individuals as rows with genotype format as A/A."),
    bsTooltip(id = "q111",title = "If the loaded file has a header, select Yes. If it is a single column file with no header, select No."),
    bsTooltip(id = "q1110",title = "If the loaded file has a header, select Yes. If it is a single column file with no header, select No."),
    bsTooltip(id = "q11100",title = "If the loaded file has a header, select Yes. If it is a single column file with no header, select No."),
    bsTooltip(id = "q1111",title = "If the loaded file has a header, select Yes. If it is a single column file with no header, select No."),
    bsTooltip(id = "q11111",title = "If the loaded file has a header, select Yes. If it is a single column file with no header, select No."),
    bsTooltip(id = "q52",title = "Select 2 if individuals are diploids and 3 if they are triploids"),
    bsTooltip(id = "q53",title = "Select SNP if markers in your dataset are SNP and microsat if it is microsatellite markers"),
    bsTooltip(id = "q54",title = "Select SNP if markers in your dataset are SNP and microsat if it is microsatellite markers"),
    bsTooltip(id = "q6",title = "APIS can use two differents method for assignment : likelihood and exclusion.<br/>It selects the best pair of parents based on differents criteria<br/>See the documentation of APIS for more details"),
    bsTooltip(id = "q61",title = "APIS can use two threshold for exclusion method.<br/>Error rate: automatically estimate the threshold given the accepted error rate.<br/>Mismatch number: select the number of maximum mismatch."),
    bsTooltip(id = "q7",title = "APIS uses paralellization to reduce running time. Select the number of core<br/>Max is number of core of the computer -1 so you can t block your computer<br/>Max-1 recommended (default)"),
    bsTooltip(id = "q8",title = "Type the name you want for saving ;<br/>Will automatically be added : _MethodSelected.Rdata for the file that contains the result of APIS and the list of markers selected ; and _ped.txt for the pedigree result",trigger = "click"),
    bsTooltip(id = "q18",title = "Type the name you want for saving plots<br/>Will automatically be added an extension to identify each plots<br/>The same name will be use if you want an output for INFAQUA"),
    bsTooltip(id = "q9",title = "Select the file finishing by _MethodSelected.Rdata created in the APIS part"),
    bsTooltip(id = "q10",title = "A file containing at least two rows : SampleName (or CodeBarre) & Facto like so :<br/>SampleName  Facto<br/>Name1  Facto1<br/>Name2  Facto2<br/>Name3  Facto3<br/>...  ...<br/>The SampleName (or CodeBarre) variable is the name of the sample as in the genotyping application.<br/>The facto variable means factorial. It is used to verify whether parents of an idividual have met or not.")
  )
  
  # Define server logic
  server <- function(input, output,session) {
    ##### Event from APIS #####
    dataset = reactiveValues(off=data.frame(),par=data.frame(),sexe=data.frame(),snp_off=data.frame(),snp_par=data.frame(),snp_to_keep=c(),snp_kept=c(),
                             nbMarker=c(),AlertNbMarker=NULL,
                             APIS_launched=0,
                             apis_likelihood=data.frame(),apis_exclusion=data.frame(),
                             dta_sex_load=FALSE,tmp_sexe=data.frame(),choices=NULL,displayed=data.frame(),
                             path_log='')
    
    #---Load dataset of offspring and stock the data
    observeEvent(input$data_off,{
      if (!is.null(input$data_off$datapath)){
        if (grepl(pattern = ".txt",x = input$data_off$name)){
          dataset$off = read.table(file=input$data_off$datapath)
          allele_freq = as.data.frame(get_allele_frequencies(dataset$off,ploidy_level = as.numeric(input$ploidy)))
          min_non_0 = function(x){
            x[x==min(x[x!=0])][1]
          }
          dataset$snp_off = data.frame(MarkerName=rownames(allele_freq),toKeep=TRUE,MAF=apply(X = allele_freq %>% select(.data$Freq_A,.data$Freq_T,.data$Freq_C,.data$Freq_G),MARGIN = 1,FUN = min_non_0,simplify = TRUE))
        } else if (grepl(pattern = ".Rdata",x = input$data_off$name)){
          tmp=load(file = input$data_off$datapath)
          eval(parse(text = paste0("dataset$off = ",tmp[1])))
          eval(parse(text = paste0("dataset$snp_off = ",tmp[2])))
        } else {
          stop("File extension not supported !")
        }
      }
    })
    #---Load dataset of parents and stock the data
    observeEvent(input$data_par,{
      if (!is.null(input$data_par$datapath)){
        if (grepl(pattern = ".txt",x = input$data_par$name)){
          dataset$par = read.table(file=input$data_par$datapath)
          
          allele_freq = as.data.frame(get_allele_frequencies(dataset$par,ploidy_level = 2))
          
          if (!is.null(allele_freq$Freq_NA)){
            dataset$snp_par = data.frame(MarkerName=rownames(allele_freq),toKeep=TRUE,CR=1-allele_freq$Freq_NA)
          } else {
            dataset$snp_par = data.frame(MarkerName=rownames(allele_freq),toKeep=TRUE,CR=1)
          }
          
        } else if (grepl(pattern = ".Rdata",x = input$data_par$name)){
          tmp=load(file = input$data_par$datapath)
          eval(parse(text = paste0("dataset$par = ",tmp[1])))
          eval(parse(text = paste0("dataset$snp_par = ",tmp[2])))
        } else {
          stop("File extension not supported !")
        }
      }
    })
    
    observeEvent(c(input$data_par1,input$data_par2),{
      if (!is.null(input$data_par1$datapath) & !is.null(input$data_par2$datapath)){
        if (grepl(pattern = ".txt",x = input$data_par1$name)){
          tmp1.1 = read.table(file=input$data_par1$datapath)
          
          allele_freq1 = as.data.frame(get_allele_frequencies(tmp1.1,ploidy_level = 2))
          
          if (!is.null(allele_freq1$Freq_NA)){
            tmp1.2 = data.frame(MarkerName=rownames(allele_freq1),toKeep=TRUE,CR=1-allele_freq1$Freq_NA)
          } else {
            tmp1.2 = data.frame(MarkerName=rownames(allele_freq1),toKeep=TRUE,CR=1)
          }
          
        } else if (grepl(pattern = ".Rdata",x = input$data_par1$name)){
          tmp=load(file = input$data_par1$datapath)
          eval(parse(text = paste0("tmp1.1 = ",tmp[1])))
          eval(parse(text = paste0("tmp1.2 = ",tmp[2])))
        } else {
          stop("File extension not supported !")
        }
        
        if (grepl(pattern = ".txt",x = input$data_par2$name)){
          tmp2.1 = read.table(file=input$data_par2$datapath)
          
          allele_freq2 = as.data.frame(get_allele_frequencies(tmp2.1,ploidy_level = 2))
          
          if (!is.null(allele_freq2$Freq_NA)){
            tmp2.2 = data.frame(MarkerName=rownames(allele_freq2),toKeep=TRUE,CR=1-allele_freq2$Freq_NA)
          } else {
            tmp2.2 = data.frame(MarkerName=rownames(allele_freq2),toKeep=TRUE,CR=1)
          }
          
        } else if (grepl(pattern = ".Rdata",x = input$data_par2$name)){
          tmp=load(file = input$data_par2$datapath)
          eval(parse(text = paste0("tmp2.1 = ",tmp[1])))
          eval(parse(text = paste0("tmp2.2 = ",tmp[2])))
        } else {
          stop("File extension not supported !")
        }
        
        marker_shared = tmp1.2$MarkerName[which(tmp1.2$MarkerName %in% tmp2.2$MarkerName)]
        
        dataset$snp_par = data.frame(MarkerName=marker_shared,
                                     CR=(tmp1.2$CR[tmp1.2$MarkerName %in% marker_shared]*nrow(tmp1.1)+tmp2.2$CR[tmp2.2$MarkerName %in% marker_shared]*nrow(tmp2.1))/(nrow(tmp1.1)+nrow(tmp2.1)),
                                     toKeep=TRUE)
        
        dataset$par = rbind(tmp1.1 %>% select(all_of(marker_shared)),tmp2.1 %>% select(all_of(marker_shared)))
        
        dataset$sexe = rbind(data.frame(SampleName=rownames(tmp1.1),Sexe=1),data.frame(SampleName=rownames(tmp2.1),Sexe=2))
      }
    })
    
    #---Load dataset of parents with their sexe and stock the data
    # Verify that 'SampleName' and 'Sexe' is a variable : if not, detection to change the colname of the samples name to match with SampleName
    observeEvent(c(input$sexe_par,input$header_sexe_par),{
      if (!is.null(input$sexe_par$datapath) & ! is.null(input$header_sexe_par)){
        dataset$dta_sex_load=TRUE
        header=ifelse(input$header_sexe_par=='Yes',TRUE,FALSE)
        dataset$sexe=read.table(file = input$sexe_par$datapath,header=header)
        if (ncol(dataset$sexe)==1){
          delim=find_delim(readLines(con=input$sexe_par$datapath,n = 1))
          dataset$sexe=read.table(file = input$sexe_par$datapath,header=header,sep=delim)
        }
        dataset$choices = colnames(dataset$sexe)
        dataset$tmp_sexe=dataset$sexe
      }
    })
    
    #--- Change SampleName and Sex column -----
    output$dta_sex_load <- reactive({
      dataset$dta_sex_load
    })
    outputOptions(output, 'dta_sex_load', suspendWhenHidden=FALSE)
    
    #---SelectInput with the names of the different columns : select the sample name
    output$uiChangeSN = renderUI({
      if (dataset$dta_sex_load){
        selectInput(inputId = "newSN",
                    label = div("Choose the variable corresponding to SampleName",
                                tipify(el = bsButton(inputId = "1",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "Please select the corresponding SampleName column so that APIS could run. See the table nearby to help you.")),
                    choices = dataset$choices,selected = dataset$choices[1],multiple = FALSE)
      }
    })
    
    #---SelectInput with the names of the different columns : select the sample name
    output$uiChangeSe = renderUI({
      if (dataset$dta_sex_load){
        selectInput(inputId = "newSe",
                    label = div("Choose the variable corresponding to Sex",
                                tipify(el = bsButton(inputId = "2",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "Please select the corresponding Sex column so that APIS could run. See the table nearby to help you.")),
                    choices = dataset$choices,selected = dataset$choices[2],multiple = FALSE)
      }
    })
    #---Change the name in regard to what is selected as samplename / sexe
    observeEvent(c(input$newSN,input$newSe),{
      if (!is.null(input$newSN) & !is.null(input$newSe)){
        dataset$sexe = dataset$tmp_sexe
        dataset$displayed = head(dataset$sexe)
        colnames(dataset$sexe)[which(colnames(dataset$sexe)==input$newSN)]='SampleName'
        colnames(dataset$sexe)[which(colnames(dataset$sexe)==input$newSe)]='Sexe'
        
        colnames(dataset$displayed)[which(colnames(dataset$displayed)==input$newSN)]=paste0(input$newSN,' (SampleName)')
        colnames(dataset$displayed)[which(colnames(dataset$displayed)==input$newSe)]=paste0(input$newSe,' (Sex)')
      }
    })
    #---Head of the dataset to help select corresponding columns
    output$head_sex = renderDataTable({
      datatable(head(dataset$displayed),rownames = FALSE,options = list(dom = 't'),
                caption = htmltools::tags$caption( style = 'caption-side: top; text-align: center; color:black;  font-size:150% ;','Head data sex'))
    })
    
    #---Find the SNPs in common between parents and offspring => keeps only good snp !!
    observeEvent(c(dataset$snp_off,dataset$snp_par),{
      if (length(dataset$snp_off)>0 & length(dataset$snp_par)>0){
        m_par_true = dataset$snp_par$MarkerName[dataset$snp_par$toKeep]
        m_off_true = dataset$snp_off$MarkerName[dataset$snp_off$toKeep]
        dataset$snp_to_keep = m_off_true[which(m_off_true %in% m_par_true)]
      }
    })
    #---textInput to select the number of SNP for assignment
    output$nb_snp_poss = renderUI({
      if (length(dataset$snp_to_keep)>0 & is.null(input$snp_in$datapath)){
        textInput(inputId = "nb_snp",value = paste0(min(100,length(dataset$snp_to_keep))),
                  label = div("Choose the number of markers to be used for the assignment",
                              tipify(el = bsButton(inputId = "3",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                     title = "Number of markers to select for the assignment after ordering them by call rate of parents and MAF of offspring (only for SNP markers).")))
      }
    })
    #---Stock the number of markers asked by user and display warning message if necessary
    observeEvent(input$nb_snp,{
      suppressWarnings(a <- as.numeric(input$nb_snp))
      if (is.na(a)){
        dataset$AlertNbMarker = "Warning : The input must be a numeric !"
      } else {
        if (a>length(dataset$snp_to_keep)){ # a>nb marker shared between par and off
          if (length(dataset$snp_to_keep)>1500){
            dataset$AlertNbMarker = paste0("Warning : There is only ",length(dataset$snp_to_keep)," markers shared between parents and offspring ! Number max of markers will be used. Additional warning : The computation time might be high ! 100 to 500 markers is usually enough for assignment.")
          } else {
            dataset$AlertNbMarker = paste0("Warning : There is only ",length(dataset$snp_to_keep)," markers shared between parents and offspring ! Number max of markers will be used.")
          }
          dataset$nbMarker = length(dataset$snp_to_keep)
        } else if (a>1500){
          dataset$AlertNbMarker = "Warning : The computation time might be high ! 100 to 500 markers is usually enough for assignment."
          dataset$nbMarker = a
        } else if (a<1){
          dataset$AlertNbMarker = "Warning : The number of marker cannot be null or negative. Automatically set to 10 and less (if 10 not possible)."
          dataset$nbMarker = ifelse(length(dataset$snp_to_keep)>9,10,1)
        } else {
          dataset$AlertNbMarker = NULL
          dataset$nbMarker = a
        }
      }
    })
    #--- Warning message du to marker number
    output$WarningMarker = renderText({
      dataset$AlertNbMarker
    })
    
    # #---RenderUI for exclusion threshold : to control max value in function of dataset length
    # output$exclusionErrorUI = renderUI({
    #   if (length(dataset$snp_kept)>0){
    #     sliderInput(inputId = "exclusionError",
    #                 label = div("Choose the threshold for mismatch",
    #                             tipify(el = bsButton(inputId = "4",label = "",icon = icon("question"), style = "info", size = "extra-small"),
    #                                    title = "Choose the number of mismatch allowed for the assignment for the best pair of parents.")),
    #                 min = 0,max = max(10,floor(0.2*length(dataset$snp_kept))),value = 0,step = 1)
    #   }
    # })
    graphApis = reactiveValues(p1=ggplot(),p2=ggplot(),p3=ggplot(),tot=ggplot())
    ##### Graph APIS #####
    output$graph_res = renderPlot({
      if (length(dataset$apis_likelihood)>0){
        p1=graphApis$p1
        if (input$method=='likelihood'){
          p3=graphApis$p3
          
          THRESHOLD=estimate_mendel_threshold(dataset$apis_likelihood,as.numeric(input$acceptError))
          p2=graphApis$p2 + geom_vline(xintercept=THRESHOLD)
        } else if (input$method=='exclusion'){
          p2=graphApis$p2
          if (input$exclu_thres=='Error rate'){
            THRESHOLD=estimate_exclusion_threshold(dataset$apis_exclusion,as.numeric(input$acceptError))
          } else {
            THRESHOLD=as.numeric(input$acceptMismatch)
          }
          p3=graphApis$p3 + geom_vline(xintercept=THRESHOLD+0.5)
        }
        graphApis$tot = plot_grid(plot_grid(p1,p2),p3,nrow = 2)
        graphApis$tot
      }
    })
    
    #---To specify that APIS has been launched and that the rest of the ui can be displayed
    output$APIS_launched <- reactive({
      dataset$APIS_launched
    })
    outputOptions(output, 'APIS_launched', suspendWhenHidden=FALSE)
    
    #---To launch APIS by beeing sure that there is everything required
    output$go_apis = renderUI({
      if (!is.null(input$data_off) & (!is.null(input$data_par) | (!is.null(input$data_par1) & !is.null(input$data_par2))) & (length(dataset$nbMarker)==1 | !is.null(input$snp_in$datapath)) & !is.null(input$save_name)){
        actionButton(inputId = "launch1",label = "Launch APIS assignment")
      }
    })
    
    output$acceptMismatch0 = renderUI({
      if (input$exclu_thres=='Mismatch number' & input$method=='exclusion' & !is.null(dataset$apis_exclusion)){
        sliderInput(inputId = "acceptMismatch",label = "Number of mismatch allowed",min = 0,max = max(dataset$apis_exclusion$mismatch_2,na.rm = TRUE)+5,value = 5,step = 1)
      }
    })
    
    ##### Launch APIS #####
    observeEvent(input$launch1,{
      if (length(dataset$off)>0 & length(dataset$par)>0 & (length(dataset$nbMarker)==1 | !is.null(input$snp_in$datapath)) & !is.null(input$save_name)){
        if (!dir.exists("./log")){
          dir.create("./log")
        }
        date_time = Sys.time()
        date_time=gsub(pattern = "-",replacement = "",x = date_time)
        date_time=gsub(pattern = ":",replacement = "",x = date_time)
        # date_time=gsub(pattern = " CEST",replacement = "_",x = date_time)
        date_time=gsub(pattern = " ",replacement = "_",x = date_time)
        date_time=substr(x = date_time,start = 3,stop = 16)
        path_log = paste0("./log/",date_time,"_log_APIS.txt")
        dataset$path_log=path_log
        write(x = "-----Launching of APIS-----",file = path_log)
        print("-----Launching of APIS-----")
        write(x = paste0("Offspring file : ",input$data_off$name),file = path_log,append = TRUE)
        if (input$both_parents_in){
          if (is.null(input$sexe_par$datapath)){
            stop("Must provide dataset with sexe of each parents !")
          }
          write(x = paste0("Parents file : ",input$data_par$name),file = path_log,append = TRUE)
          write(x = paste0("Parental sex file : ",input$sexe_par$name),file = path_log,append = TRUE)
        } else {
          write(x = paste0("Parents file sire : ",input$data_par1$name),file = path_log,append = TRUE)
          write(x = paste0("Parents file dam : ",input$data_par2$name),file = path_log,append = TRUE)
        }
        
        if (!is.null(input$snp_in$datapath)){
          snp_kept = read.table(file = input$snp_in$datapath,header = FALSE)
          snp_kept = snp_kept[,1]
          write(x = paste0("Marker file : ",input$snp_in$name),file = path_log,append = TRUE)
          snp_kept=snp_kept[which(snp_kept %in% colnames(dataset$off) & snp_kept %in% colnames(dataset$par))]
        } else {
          if (input$markerType2=='SNP'){
            snp_kept = dataset$snp_off %>%
              select(.data$MarkerName,.data$MAF) %>%
              left_join(dataset$snp_par %>% select(.data$MarkerName,.data$CR),by=c("MarkerName"="MarkerName")) %>%
              arrange(desc(.data$MAF),desc(.data$CR)) %>%
              filter(!is.na(.data$MAF),!is.na(.data$CR)) %>%
              select(.data$MarkerName)
            
          } else { # input$markerType2=='microsat'
            snp_kept = dataset$snp_off %>%
              select(.data$MarkerName) %>%
              left_join(dataset$snp_par %>% select(.data$MarkerName,.data$CR),by=c("MarkerName"="MarkerName")) %>%
              arrange(desc(.data$CR)) %>%
              filter(!is.na(.data$CR)) %>%
              select(.data$MarkerName)
          }
          snp_kept = snp_kept$MarkerName[1:min(as.numeric(dataset$nbMarker),length(snp_kept$MarkerName))]
        }
        if (length(snp_kept) < 1){
          write(x = "Error : 0 marker kept !",file = path_log,append = TRUE)
          write(x = "If a list is provided: dataset doesnt have those markers ; if not, dataset doesnt share marker !",file = path_log,append = TRUE)
          stop("If a list is provided: dataset doesnt have those markers ; if not, dataset doesnt share marker !")
        }
        write(x = paste0("Marker type : ",input$markerType2),file = path_log,append = TRUE)
        write(x = paste0("Number of snp kept : ",length(snp_kept)),file = path_log,append = TRUE)
        
        # New 2.0.5
        dataset$sexe = dataset$sexe %>% select(.data$SampleName,.data$Sexe) %>% distinct()
        
        offspring = dataset$off %>%
          as.data.frame() %>%
          select(all_of(snp_kept)) %>%
          as.matrix()
        sire0 = dataset$par %>%
          as.data.frame() %>%
          select(all_of(snp_kept))
        sire = sire0[which(toupper(rownames(sire0)) %in% toupper(dataset$sexe$SampleName[dataset$sexe$Sexe==1 | dataset$sexe$Sexe==3])),] %>%
          as.matrix()
        dam0 = dataset$par %>%
          as.data.frame() %>%
          select(all_of(snp_kept))
        dam = dam0[which(toupper(rownames(dam0)) %in% toupper(dataset$sexe$SampleName[dataset$sexe$Sexe==2 | dataset$sexe$Sexe==4])),] %>%
          as.matrix()
        
        par_nam_tot=toupper(rownames(dataset$par))
        par_nam_rest=toupper(c(rownames(sire),rownames(dam)))
        par_nam_list=toupper(dataset$sexe$SampleName)
        
        no_sex = which(! par_nam_tot %in% par_nam_rest)
        sex_but_no = which(! par_nam_list %in% par_nam_rest)
        if (length(no_sex)>0){
          warning("Some parent(s) does not have a sex assigned !")
          print(par_nam_tot[no_sex])
          write(x = paste0("WARNING -- Sex missing for some parents : ",par_nam_tot[no_sex]),file = path_log,append = TRUE)
        }
        if (length(sex_but_no)>0){
          warning("Some parent(s) are missing in the dataset !")
          print(par_nam_list[sex_but_no])
          write(x = paste0("WARNING -- Parents missing : ",par_nam_list[sex_but_no]),file = path_log,append = TRUE)
        }
        write(x = paste0("Offspring ploidy : ",input$ploidy),file = path_log,append = TRUE)
        write(x = paste0("Number of cores (parallelization) : ",input$core),file = path_log,append = TRUE)
        write(x = paste0("Launch time : ",Sys.time()),file = path_log,append = TRUE)
        if (as.numeric(input$ploidy)==3){
          if (FALSE){ # input$try_recom
            v_recom = seq(0.1,0.9,0.1)
            v_eff = NULL
            m_re = 0
            for (v in v_recom){
              to_save=APIS_3n(offspring_genotype = offspring[,1:min(100,ncol(offspring))],
                              sire_genotype = sire[,1:min(100,ncol(offspring))],
                              dam_genotype = dam[,1:min(100,ncol(offspring))],
                              number_cores = 10,
                              verbose=TRUE,
                              simulation_if_small = FALSE,
                              t_recom=v)
              if (m_re<mean(to_save$log_file$delta_1_2,na.rm = T)){
                v_eff = v
                m_re = mean(to_save$log_file$delta_1_2,na.rm = T)
              }
            }
            res_APIS=APIS_3n(offspring_genotype = offspring,
                             sire_genotype = sire,
                             dam_genotype = dam,
                             number_cores = as.numeric(input$core),
                             verbose=TRUE,
                             simulation_if_small = FALSE,
                             t_recom = v_eff,
                             method = "mendel")
          } else {
            res_APIS=APIS_3n(offspring_genotype = offspring,
                             sire_genotype = sire,
                             dam_genotype = dam,
                             number_cores = as.numeric(input$core),
                             verbose=TRUE,
                             simulation_if_small = FALSE,
                             t_recom = 0.5,
                             method = "")
          }
        } else if(as.numeric(input$ploidy)==2){
          res_APIS=APIS_2n(offspring_genotype = offspring,
                           sire_genotype = sire,
                           dam_genotype = dam,
                           number_cores = as.numeric(input$core),
                           verbose=TRUE,
                           simulation_if_small = FALSE,
                           method = "")
        } else {
          stop("Incorrect number of ploidy !")
        }
        write(x = paste0("End time : ",Sys.time()),file = path_log,append = TRUE)
        dataset$snp_kept = snp_kept
        dataset$apis_likelihood = res_APIS$log_file_likelihood
        dataset$apis_exclusion = res_APIS$log_file_exclusion
        dataset$APIS_launched = 1
        dataset$displayed = data.frame() # reinitialisation so that it does not overcharge the user experience
        #---Generate graph
        tmp = data.frame(Value1 = c(dataset$apis_likelihood$probability_1,dataset$apis_likelihood$probability_2),
                         Value2 = c(dataset$apis_likelihood$delta_1_2,dataset$apis_likelihood$delta_2_3),
                         Which = as.factor(rep(c(1,2),each=length(dataset$apis_likelihood$probability_1))))
        graphApis$p1=ggplot(data=tmp,aes(x = .data$Value1,color = .data$Which,fill = .data$Which))+
          geom_histogram(alpha = 0.25,position = "identity")+
          scale_fill_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
          scale_color_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
          theme_bw()+
          labs(x = "Mendelian probability",y="Count",fill="Assignment",col="Assignment")
        graphApis$p2=ggplot(data=tmp,aes(x = .data$Value2,color = .data$Which,fill = .data$Which))+
          geom_histogram(alpha = 0.25,position = "identity")+
          scale_fill_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
          scale_color_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
          theme_bw()+
          labs(x = "Delta of mendelian probability",y="Count",fill="Assignment",col="Assignment")
        
        tmp = data.frame(Value = c(dataset$apis_exclusion$mismatch_1,dataset$apis_exclusion$mismatch_2),
                         Which = as.factor(rep(c(1,2),each=length(dataset$apis_exclusion$mismatch_1))))
        
        graphApis$p3=ggplot(data=tmp,aes(x = .data$Value,color = .data$Which,fill = .data$Which))+
          geom_histogram(alpha = 0.25,position = "identity")+
          scale_fill_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
          scale_color_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
          theme_bw()+
          labs(x = "Number of mismatch",y="Count",fill="Assignment",col="Assignment")
        print("-----End APIS-----")
      }
    })
    
    observeEvent(input$SaveAPIS,{
      if (length(dataset$apis_likelihood)>0){
        write(x = "----Saving APIS assignment----",file = dataset$path_log,append = TRUE)
        write(x = paste0("Method : ",input$method),file = dataset$path_log,append = TRUE)
        mismatch_error=FALSE
        if (input$method=='exclusion'){
          if (input$exclu_thres=='Error rate'){
            THRESHOLD=as.numeric(input$acceptError)
            estiThreshold=estimate_exclusion_threshold(dataset$apis_exclusion,as.numeric(input$acceptError))
            write(x = paste0("Error accepted : ",as.numeric(input$acceptError*100),"%"),file = dataset$path_log,append = TRUE)
            write(x = paste0("Mismatch threshold : ",estiThreshold),file = dataset$path_log,append = TRUE)
          } else {
            mismatch_error=TRUE
            THRESHOLD=as.numeric(input$acceptMismatch)
            estiThreshold=THRESHOLD
            write(x = paste0("Mismatch threshold : ",estiThreshold),file = dataset$path_log,append = TRUE)
          }
          ped = dataset$apis_exclusion %>% select(.data$offspring,.data$sire_1,.data$dam_1) %>% rename(sire=.data$sire_1,dam=.data$dam_1)
          ind_na=unique(which(dataset$apis_exclusion$mismatch_1>estiThreshold | dataset$apis_exclusion$mismatch_1==dataset$apis_exclusion$mismatch_2))
          ped[ind_na,2:3]=c(NA,NA)
          log_APIS = dataset$apis_exclusion
        } else if (input$method=='likelihood'){
          THRESHOLD=as.numeric(input$acceptError)
          estiThreshold=estimate_mendel_threshold(dataset$apis_likelihood,as.numeric(input$acceptError))
          ped = dataset$apis_likelihood %>% select(.data$offspring,.data$sire_1,.data$dam_1) %>% rename(sire=.data$sire_1,dam=.data$dam_1)
          ind_na=unique(which(dataset$apis_likelihood$delta_1_2<estiThreshold))
          ped[ind_na,2:3]=c(NA,NA)
          write(x = paste0("Error accepted : ",as.numeric(input$acceptError*100),"%"),file = dataset$path_log,append = TRUE)
          write(x = paste0("Threshold of delta probability estimated : ",estiThreshold),file = dataset$path_log,append = TRUE)
          log_APIS = dataset$apis_likelihood
        }
        write(x = paste0("Saving name (in './Results_APIS' folder) : ",input$save_name),file = dataset$path_log,append = TRUE)
        if (!dir.exists("./Results_APIS")){
          dir.create("./Results_APIS")
        }
        snp_kept = dataset$snp_kept
        df_par = dataset$sexe
        save(log_APIS,ped,snp_kept,df_par,THRESHOLD,estiThreshold,mismatch_error,file = paste0("./Results_APIS/",input$save_name,"_",input$method,".Rdata"))
        ggsave(graphApis$tot,filename = paste0("./Results_APIS/",input$save_name,"_",input$method,".png"),width = 7,height = 7)
        write.table(x = ped,file = paste0("./Results_APIS/",input$save_name,"_ped.txt"),quote = FALSE,row.names = FALSE)
        write.table(x = log_APIS,file = paste0("./Results_APIS/",input$save_name,"_logfile.csv"),sep=";",quote = FALSE,row.names = FALSE)
        print("-----APIS files saved ! -----")
      }
    })
    
    ##### Event from Verification #####
    verif = reactiveValues(out1=NULL,out2=NULL,out3=NULL,out4=NULL,tab=data.frame(),
                           threshold=NULL,
                           dta_accou_load=FALSE,
                           sentence_dam1=NULL,sentence_sire1=NULL,sentence_dam2=NULL,sentence_sire2=NULL,
                           data=data.frame(),accou=data.frame(),ped=data.frame(),accou_tmp=data.frame(),
                           dta_fac_load=FALSE,choices=NULL,displayed=data.frame())
    
    #---Summary of the assignment
    output$txt1 = renderText({
      verif$out1
    })
    output$txt2 = renderText({
      verif$out2
    })
    output$txt3 = renderText({
      verif$out3
    })
    output$txt4 = renderText({
      verif$out4
    })
    #---Summary number of offspring by sexe
    output$sent1.1 = renderText({
      verif$sentence_dam1
    })
    output$sent2.1 = renderText({
      verif$sentence_sire1
    })
    output$sent1.2 = renderText({
      verif$sentence_dam2
    })
    output$sent2.2 = renderText({
      verif$sentence_sire2
    })
    #---Table with assignment where parents should not have met
    output$tab = renderDataTable({
      if (is.null(input$tab_accou$datapath)){
        datatable(verif$tab,rownames = FALSE,#options = list(dom = 't'),
                  caption = htmltools::tags$caption( style = 'caption-side: top; text-align: center; color:black;  font-size:150% ;','Un-assigned offspring'))
      } else {
        datatable(verif$tab,rownames = FALSE,#options = list(dom = 't'),
                  caption = htmltools::tags$caption( style = 'caption-side: top; text-align: center; color:black;  font-size:150% ;','Out-of-plan assignment'))
      }
    })
    
    #---Load dataset of parents with their sexe and stock the data
    # Verify that 'SampleName' is a variable : if not, detection to change the colname of the samples name to match with SampleName
    observeEvent(c(input$tab_accou,input$header_tab_accou),{
      if (!is.null(input$tab_accou$datapath) & ! is.null(input$header_tab_accou)){
        verif$dta_fac_load=TRUE
        header=ifelse(input$header_tab_accou=='Yes',TRUE,FALSE)
        verif$accou=read.table(file = input$tab_accou$datapath,header=header)
        if (length(verif$accou)[1]==1){
          delim = find_delim(readLines(con=input$tab_accou$datapath,n=1))
          verif$accou=read.table(file = input$tab_accou$datapath,header=header,sep=delim)
        }
        verif$choices = colnames(verif$accou)
        verif$tmp_accou=verif$accou
      }
    })
    
    #--- Change SampleName and Facto when data loaded VERIF -----
    output$dta_fac_load <- reactive({
      verif$dta_fac_load
    })
    outputOptions(output, 'dta_fac_load', suspendWhenHidden=FALSE)
    
    #--- SelectInput with the names of the different columns : select the sample name
    output$uiChangeSN2 = renderUI({
      if (verif$dta_fac_load){
        selectInput(inputId = "newSN2",
                    label = div("Choose the variable corresponding to SampleName",
                                tipify(el = bsButton(inputId = "6",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "Please select the corresponding SampleName column so that APIS could run. See the table nearby to help you.")),
                    choices = verif$choices,selected = verif$choices[1],multiple = FALSE)
      }
    })
    
    #--- SelectInput with the names of the different columns : select the Facto
    output$uiChangeFa = renderUI({
      if (verif$dta_fac_load){
        selectInput(inputId = "newFa",
                    label = div("Choose the variable corresponding to Facto",
                                tipify(el = bsButton(inputId = "7",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "Please select the corresponding Facto column so that APIS could run. See the table nearby to help you.")),
                    choices = verif$choices,selected = verif$choices[2],multiple = FALSE)
      }
    })
    #--- Change the name in regard to what is selected as samplename and facto
    observeEvent(c(input$newSN2,input$newFa,input$tab_accou),{
      if (!is.null(input$newSN2) & !is.null(input$newFa)){
        verif$accou = verif$tmp_accou
        verif$displayed = verif$accou
        colnames(verif$accou)[which(colnames(verif$accou)==input$newSN2)]='SampleName'
        colnames(verif$accou)[which(colnames(verif$accou)==input$newFa)]='Facto'
        
        colnames(verif$displayed)[which(colnames(verif$displayed)==input$newSN2)]=paste0(input$newSN2,' (SampleName)')
        colnames(verif$displayed)[which(colnames(verif$displayed)==input$newFa)]=paste0(input$newFa,' (Facto)')
      }
    })
    
    #---Head of the dataset to help select corresponding columns
    output$head_fac = renderDataTable({
      datatable(head(verif$displayed),rownames = FALSE,options = list(dom = 't'),
                caption = htmltools::tags$caption( style = 'caption-side: top; text-align: center; color:black;  font-size:150% ;','Head data factorial'))
    })
    
    ##### Graph Verification #####
    #--- Plot of proba/mismatch between couple 1 and 2
    output$plot = renderPlot({
      if (length(verif$data)>0){
        if (regexpr("likelihood",input$data_res$name)==-1){ # si il ny a pas likelihood dans le nom => donc exclusion
          tmp = data.frame(Value = c(verif$data$mismatch_1,verif$data$mismatch_2),
                           Which = as.factor(rep(c(1,2),each=length(verif$data$mismatch_1))))
          ggplot(data=tmp,aes(x = .data$Value,color = .data$Which,fill = .data$Which))+
            geom_histogram(alpha = 0.25,position = "identity")+
            scale_fill_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
            scale_color_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
            geom_vline(xintercept = verif$threshold)+
            theme_bw()+
            labs(x = "Number of mismatch",y="Count",fill="Assignment",col="Assignment")
        } else {
          tmp = data.frame(Value1 = c(verif$data$probability_1,verif$data$probability_2),
                           Value2 = c(verif$data$delta_1_2,verif$data$delta_2_3),
                           Which = as.factor(rep(c(1,2),each=length(verif$data$probability_1))))
          p1=ggplot(data=tmp,aes(x = .data$Value1,color = .data$Which,fill = .data$Which))+
            geom_histogram(alpha = 0.25,position = "identity")+
            scale_fill_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
            scale_color_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
            theme_bw()+
            labs(x = "Mendelian probability",y="Count",fill="Assignment",col="Assignment")
          p2=ggplot(data=tmp,aes(x = .data$Value2,color = .data$Which,fill = .data$Which))+
            geom_histogram(alpha = 0.25,position = "identity")+
            scale_fill_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
            scale_color_manual(values = c("1"="#56B4E9","2"="#D55E00"))+
            theme_bw()+
            geom_vline(xintercept = verif$threshold)+
            labs(x = "Delta of mendelian probability",y="Count",fill="Assignment",col="Assignment")
          plot_grid(p1,p2)
        }
      }
    })
    #--- Barplot with the number of offspring for a given parents
    output$barplot = renderPlotly({
      if(length(to_plot$ggbar)>0){
        ggplotly(to_plot$ggbar)
      }
    })
    #--- Heatmap with color for the number of offspring between sireX and damX
    output$heatmap = renderPlotly({
      if (length(to_plot$ggheat)>0){
        ggplotly(p=to_plot$ggheat,tooltip=c('x','y','fill','text'))
      }
    })
    observeEvent(input$SavePlot,{
      if (length(to_plot$ggbar$data)>0 & length(to_plot$ggheat$data)>0){
        if (!dir.exists("./Results_verif")){
          dir.create("./Results_verif")
        }
        ggsave(plot = to_plot$ggbar,filename = paste0("./Results_verif/",input$save_name2,"_barplot.png"),width = 18,height = 9)
        ggsave(plot = to_plot$ggheat+coord_fixed(),filename = paste0("./Results_verif/",input$save_name2,"_heatmap.png"),width = 12,height = 12)
        
        txt=ggplot(data=data.frame(x=0,y=0),aes(x=.data$x,y=.data$y))+
          geom_text(x=0,y=0.9,label="Summary of the APIS assignment",size=10)+
          geom_text(x=0,y=0,label=paste(verif$out1,verif$out2,verif$out3,verif$out4,paste0("File name : ",input$data_res$name),sep = "\n"),size=3)+
          geom_text(x=0,y=-0.9,label=paste(verif$sentence_dam1,verif$sentence_dam2,verif$sentence_sire1,verif$sentence_sire2,sep = "\n"),size=3)+
          xlim(-1,1)+ylim(-1,1)+
          theme(axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                panel.background = element_blank(),
                axis.text = element_blank())
        plot1 = plot_grid(txt,to_plot$ggbar,nrow=2)
        
        if (is.null(input$tab_accou$datapath)){
          txt2=ggplot(data=data.frame(x=0,y=0),aes(x=.data$x,y=.data$y))+
            geom_text(x=0,y=0.9,label="Un-assigned offspring",size=5)+
            xlim(-1,1)+ylim(-1,1)+
            theme(axis.line = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  panel.background = element_blank(),
                  axis.text = element_blank())
        } else {
          txt2=ggplot(data=data.frame(x=0,y=0),aes(x=.data$x,y=.data$y))+
            geom_text(x=0,y=0.9,label="Out-of-plan assignment",size=5)+
            xlim(-1,1)+ylim(-1,1)+
            theme(axis.line = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  panel.background = element_blank(),
                  axis.text = element_blank())
        }
        
        
        nMAX = 15
        pdf(paste0("./Results_verif/",input$save_name2,"_summary.pdf"))
        print(plot1)
        print(to_plot$ggheat)
        
        if (nrow(verif$tab)==0 & !is.null(input$tab_accou$datapath)){
          print("No out-of-plan assignment")
        } else if (nrow(verif$tab)==0 & is.null(input$tab_accou$datapath)){
          print("No un-assigned offspring")
        } else if (nrow(verif$tab)<=nMAX){
          print(txt2)
          grid.table(verif$tab,rows=NULL)
        } else {
          print(txt2)
          grid.table(verif$tab[1:nMAX,],rows=NULL)
          n = (nrow(verif$tab)-nMAX)%/%nMAX
          leftover = (nrow(verif$tab)-nMAX)%%nMAX
          for (k in 1:n){
            plot.new()
            grid.table(verif$tab[(nMAX*k+1):(nMAX*(k+1)),],rows=NULL)
          }
          if (leftover>0){
            plot.new()
            grid.table(verif$tab[(nMAX*(n+1)+1):nrow(verif$tab),],rows=NULL)
          }
        }
        dev.off()
        print("----- Plot saved ! -----")
      }
    })
    
    to_plot=reactiveValues(heatmap=data.frame(),ggbar=ggplot(),ggheat=ggplot())
    ##### Launch Verification #####
    observeEvent(input$launch_verif,{
      if (!is.null(input$data_res)){
        load(file = input$data_res$datapath) # log_APIS , ped , snp_kept , df_par , THRESHOLD , estiThreshold , mismatch_error
        verif$data = log_APIS
        verif$ped = ped
        verif$threshold = estiThreshold
        
        # Barplot
        nb_sire = verif$ped %>% group_by(.data$sire) %>% count()
        nb_dam = verif$ped %>% group_by(.data$dam) %>% count()
        tmp=data.frame(par=c(nb_sire$sire,nb_dam$dam),n=c(nb_sire$n,nb_dam$n))
        tmp$par=toupper(tmp$par)
        
        # New 2.0.5
        df_par = df_par %>% select(.data$SampleName,.data$Sexe) %>% distinct()
        
        df_par$SampleName=toupper(df_par$SampleName)
        barplot = left_join(df_par,tmp,by=c("SampleName"="par"))
        barplot$n[is.na(barplot$n)]=0
        count0 = function(n){
          length(which(n==0))
        }
        countDif0 = function(n){
          length(which(n!=0))
        }
        suma_dam = barplot %>% filter(.data$Sexe %in% c(2,4)) %>% select(.data$n) %>%
          summarise(min=min(.data$n,na.rm=T),max=max(.data$n,na.rm=T),mean=mean(.data$n,na.rm=T),
                    median=median(.data$n,na.rm=T),nb0=count0(.data$n),nbdif0=countDif0(.data$n)) %>% round(digits = 2)
        suma_sire = barplot %>% filter(.data$Sexe %in% c(1,3)) %>% select(.data$n) %>%
          summarise(min=min(.data$n,na.rm=T),max=max(.data$n,na.rm=T),mean=mean(.data$n,na.rm=T),
                    median=median(.data$n,na.rm=T),nb0=count0(.data$n),nbdif0=countDif0(.data$n)) %>% round(digits = 2)
        verif$sentence_dam1 = paste0("Dam Number of Offspring -- Min : ",suma_dam$min," ; Max : ",suma_dam$max," ; Mean : ",suma_dam$mean," ; Median : ",suma_dam$median)
        verif$sentence_dam2 = paste0("Nb no off : ",suma_dam$nb0," -- Nb with off : ",suma_dam$nbdif0)
        
        verif$sentence_sire1 = paste0("Sire Number of Offspring -- Min : ",suma_sire$min," ; Max : ",suma_sire$max," ; Mean : ",suma_sire$mean," ; Median : ",suma_sire$median)
        verif$sentence_sire2 = paste0("Nb no off : ",suma_sire$nb0," -- Nb with off : ",suma_sire$nbdif0)
        
        # Def ggbar
        order_ind = barplot %>%
          arrange(.data$Sexe,desc(.data$n))%>%
          select(SampleName)
        barplot$SampleName=factor(barplot$SampleName,levels=as.character(order_ind$SampleName))
        barplot$Sexe[barplot$Sexe==1 | barplot$Sexe==3]="Sire"
        barplot$Sexe[barplot$Sexe==2 | barplot$Sexe==4]="Dam"
        barplot$Sexe=as.factor(barplot$Sexe)
        to_plot$ggbar=ggplot(data=barplot,aes(x=.data$SampleName,y=.data$n,fill=.data$Sexe))+
          geom_bar(stat='identity')+
          scale_fill_manual(values = c("Dam" = "#009E73", "Sire" = "#56B4E9"))+
          theme_bw()+
          geom_hline(yintercept = 0)+
          labs(x="Individuals",y="Number of offspring")+
          theme(panel.grid.major.x = element_blank(),
                axis.ticks.x = element_blank(),
                # axis.text.x = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, size=5,face = "bold"))
        
        # Heatmap
        tmp=verif$ped %>% group_by(sire,dam) %>% count() %>% arrange(desc(.data$n))
        ind_sire=df_par$SampleName[df_par$Sexe==1 | df_par$Sexe==3]
        ind_dam=df_par$SampleName[df_par$Sexe==2 | df_par$Sexe==4]
        tmp2=expand.grid(Sire=ind_sire,Dam=ind_dam)
        heatmap = left_join(tmp2,tmp,by=c("Sire"="sire","Dam"="dam"))
        
        # Def ggheat
        if (!is.null(input$tab_accou)){
          verif$accou$SampleName=toupper(verif$accou$SampleName)
          tmpSire = df_par %>%
            filter(.data$Sexe==1 | .data$Sexe==3) %>%
            select(.data$SampleName) %>%
            left_join(verif$accou,by=c("SampleName"="SampleName"),multiple="all") %>%
            arrange(.data$Facto) %>%
            select(.data$SampleName)
          tmpDam = df_par %>%
            filter(.data$Sexe==2 | .data$Sexe==4) %>%
            select(.data$SampleName) %>%
            left_join(verif$accou,by=c("SampleName"="SampleName"),multiple="all") %>%
            arrange(.data$Facto) %>%
            select(.data$SampleName)
          
          f_sire=data.frame(SN=unique(tmpSire$SampleName),Fa_sire=NA)
          for (k in 1:length(f_sire$SN)){
            sire_k=f_sire$SN[k]
            f_k=paste0(sort(verif$accou$Facto[verif$accou$SampleName==sire_k]),collapse = "/")
            f_sire$Fa_sire[k]=f_k
          }
          f_dam=data.frame(SN=unique(tmpDam$SampleName),Fa_dam=NA)
          for (k in 1:length(f_dam$SN)){
            dam_k=f_dam$SN[k]
            f_k=paste0(sort(verif$accou$Facto[verif$accou$SampleName==dam_k]),collapse = "/")
            f_dam$Fa_dam[k]=f_k
          }
          
          heatmap = heatmap %>%
            left_join(f_sire,by=c("Sire"="SN")) %>%
            left_join(f_dam,by=c("Dam"="SN"))
          
          
          heatmap$Sire=factor(heatmap$Sire,levels = unique(as.character(tmpSire$SampleName)))
          heatmap$Dam=factor(heatmap$Dam,levels = unique(as.character(tmpDam$SampleName)))
          
          if (length(unique(tmpDam$SampleName))==length(tmpDam$SampleName) & length(unique(tmpSire$SampleName))==length(tmpDam$SampleName)){
            heatmap$Same = heatmap$Fa_sire==heatmap$Fa_dam
          } else {
            heatmap$Same = NA
          }
          to_plot$ggheat=ggplot(data = heatmap,aes(x = .data$Sire,y = .data$Dam,fill = .data$n,text = paste0("Fac Sire : ",.data$Fa_sire,"<br>Fac Dam : ",.data$Fa_dam)))+
            geom_tile(col=heatmap$Same,linewidth=0.05,width=0.8,height=0.8)+
            scale_fill_gradient2(low = "#66CCFF", high = "#D55E00", mid="#F0E442",
                                 na.value="white",name="Nombre de \ndescendants",midpoint =5,space = "Lab")+
            theme_bw()+
            theme(panel.grid = element_blank())+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                             size = 5, hjust = 1,face = "bold"))+
            theme(axis.text.y = element_text(angle = 0, vjust = 0.5,
                                             size = 5, hjust = 1,face="bold"))
          
          # Indiv Hors Plan
          if (length(verif$accou$SampleName)==length(unique(verif$accou$SampleName))){
            tab_assign2 = left_join(verif$ped,verif$accou,by=c('sire'='SampleName')) %>% rename(FactoSire=.data$Facto) %>% select(.data$offspring,.data$sire,.data$dam,.data$FactoSire)
            # names(tab_assign2)[6]='FactoSire'
            # tab_assign2 = tab_assign2[,-c(4,5)]
            tab_assign3 = left_join(tab_assign2,verif$accou,by=c('dam'='SampleName')) %>% rename(FactoDam=.data$Facto) %>% select(.data$offspring,.data$sire,.data$dam,.data$FactoSire,.data$FactoDam)
            # names(tab_assign3)[7]='FactoDam'
            # tab_assign3 = tab_assign3[,-c(5,6)]
            indiv_pb=tab_assign3$offspring[which(tab_assign3$FactoSire!=tab_assign3$FactoDam)]
            
            n_pb=length(indiv_pb)
            n_na=length(which(is.na(tab_assign3$sire)))
            n_tot=length(verif$ped$offspring)
          } else { # si il y a des individus avec plusieurs facto
            tab_assign2 = left_join(verif$ped,verif$accou,by=c('sire'='SampleName'),multiple = "all") %>% rename(FactoSire=.data$Facto) %>% select(.data$offspring,.data$sire,.data$dam,.data$FactoSire)
            tab_assign3 = left_join(tab_assign2,verif$accou,by=c('dam'='SampleName'),multiple = "all") %>% rename(FactoDam=.data$Facto) %>% select(.data$offspring,.data$sire,.data$dam,.data$FactoSire,.data$FactoDam)
            
            n_pb=0
            n_na=0
            indiv_pb=c()
            n_tot=length(verif$ped$offspring)
            for (off in unique(verif$ped$offspring)){
              tmp=tab_assign3[tab_assign3$offspring==off,]
              n.row=nrow(tmp)
              if (n.row==1){ # parents have only 1 facto
                if (!is.na(tmp$sire)){
                  if (tmp$FactoSire!=tmp$FactoSire){
                    n_pb=n_pb+1
                    indiv_pb=c(indiv_pb,off)
                  }
                } else {
                  n_na=n_na+1
                }
              } else { # un des parents au moins fait parti de plusieurs factorielles
                # parents ne peuvent pas etre NA si plusieurs lignes (par construction)
                same_facto = FALSE
                for (k in 1:n.row){
                  if (tmp$FactoDam[k]==tmp$FactoSire[k]){
                    same_facto = TRUE
                  }
                }
                if (!same_facto){
                  n_pb=n_pb+1
                  indiv_pb=c(indiv_pb,off)
                }
              }
            }
          }
          verif$out1=paste0("There is/are ",n_na," no assigned offspring(s) (",round((n_na/n_tot)*100,2),"%)\nand among assigned offspring(s) ",n_pb," have parents that are not in the same factorial (",round((n_pb/n_tot)*100,2),"%).")
          verif$out2=paste0("Real assignment rate : ",n_tot-n_pb-n_na,"/",n_tot,"=",round((n_tot-n_pb-n_na)*100/n_tot,2),"%.")
          verif$out3=paste0("The assignment was done using ",length(snp_kept)," markers.")
          if (!exists("mismatch_error")){
            verif$out4=paste0("The maximum theoretical error rate for this assignment is ",round(THRESHOLD,2)*100,"%.")
          } else {
            if (!mismatch_error){
              verif$out4=paste0("The maximum theoretical error rate for this assignment is ",round(THRESHOLD,2)*100,"%.")
            }
          }
          verif$tab=tab_assign3[tab_assign3$offspring %in% indiv_pb,]
        } else {
          to_plot$ggheat=ggplot(data = heatmap,aes(x = .data$Sire,y = .data$Dam,fill = .data$n))+
            geom_tile(linewidth=0.05,width=0.8,height=0.8)+
            scale_fill_gradient2(low = "#66CCFF", high = "#D55E00", mid="#F0E442",
                                 na.value="white",name="Nombre de \ndescendants",midpoint =5,space = "Lab")+
            theme_bw()+
            theme(panel.grid = element_blank())+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                             size = 5, hjust = 1,face = "bold"))+
            theme(axis.text.y = element_text(angle = 0, vjust = 0.5,
                                             size = 5, hjust = 1,face="bold"))
          
          n_na=length(which(is.na(verif$ped$sire)))
          n_tot=length(verif$ped$sire)
          
          verif$out1=paste0("There is/are ",n_na," no assigned offspring(s) (",round((n_na/n_tot)*100,2),"%).")
          verif$out2=paste0("Real assignment rate : ",n_tot-n_na,"/",n_tot,"=",round((n_tot-n_na)*100/n_tot,2),"%.")
          verif$out3=paste0("The assignment was done using ",length(snp_kept)," markers.")
          if (!exists("mismatch_error")){
            verif$out4=paste0("The maximum theoretical error rate for this assignment is ",round(THRESHOLD,2)*100,"%.")
          } else {
            if (!mismatch_error){
              verif$out4=paste0("The maximum theoretical error rate for this assignment is ",round(THRESHOLD,2)*100,"%.")
            }
          }
          verif$tab=verif$ped[which(is.na(verif$ped$sire)),]
        }
        verif$displayed = data.frame() # reinitialisation so that it does not overcharge the user experience
      }
    })
    
    ##### Event from Formating #####
    formating=reactiveValues(data=NULL,head=NULL,end=NULL,
                             dataMap=NULL,
                             colFormat=FALSE,colMap=FALSE,
                             Lpar=NULL,import_vcf = FALSE)
    
    #---Display launch button when ready
    output$go_format = renderUI({
      if (length(formating$head)!=0 | formating$import_vcf){
        actionButton(inputId = "format",label = "Launch formatting")
      }
    })
    
    #---Load list with parents
    observeEvent(c(input$list_par,input$header_list_par),{
      if (!is.null(input$list_par$datapath) & !is.null(input$header_list_par)){
        if (input$header_list_par=='Yes'){
          if (grepl(pattern = ".csv",x = input$list_par$datapath)){
            tmp1 = read.table(file = input$list_par$datapath,header = TRUE,comment.char = "#",sep=";")
            tmp2 = read.table(file = input$list_par$datapath,header = TRUE,comment.char = "#",sep=",")
            if (ncol(tmp1)>ncol(tmp2)){
              formating$Lpar=tmp1
            } else {
              formating$Lpar=tmp2
            }
            rm(tmp1,tmp2)
            gc()
          } else {
            formating$Lpar = read.table(file = input$list_par$datapath,header = TRUE,comment.char = "#")
          }
        } else {
          if (grepl(pattern = ".csv",x = input$list_par$datapath)){
            tmp1 = read.table(file = input$list_par$datapath,header = FALSE,comment.char = "#",sep=";")
            tmp2 = read.table(file = input$list_par$datapath,header = FALSE,comment.char = "#",sep=",")
            if (ncol(tmp1)>ncol(tmp2)){
              formating$Lpar=tmp1
            } else {
              formating$Lpar=tmp2
            }
            rm(tmp1,tmp2)
            gc()
          } else {
            formating$Lpar = read.table(file = input$list_par$datapath,header = FALSE,comment.char = "#")
          }
        }
      }
    })
    
    #---Check if there is a header and ask for the good column
    output$col_head0 = renderUI({
      if (!is.null(formating$Lpar)){
        sliderInput(inputId = "col_head",
                    label = div("Choose the number of the column corresponding to sample names",
                                tipify(el = bsButton(inputId = "5",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "Select the number corresponding to the column with the sample names. Look at the table displayed to be sure.")),
                    min = 1,max = min(20,ncol(formating$Lpar)),value = 1,step = 1)
      }
    })
    #---Display head of the list with parents
    # output$headLpar = renderTable({
    #   head(formating$Lpar)
    # })
    output$headLpar = renderDataTable({
      datatable(head(formating$Lpar),rownames = FALSE,options = list(dom = 't'),
                caption = htmltools::tags$caption( style = 'caption-side: top; text-align: center; color:black;  font-size:150% ;','Head list parents'))
    })
    
    #---Load genotype file to format and display the head of df to help select the columns
    observeEvent(c(input$to_format,input$header_file_format),{
      if (! is.null(input$to_format) & ! is.null(input$header_file_format)){
        header = ifelse(input$header_file_format=='Yes',TRUE,FALSE)
        if (grepl(pattern = ".vcf",x = input$to_format$datapath)){
          formating$import_vcf = TRUE
        } else {
          formating$import_vcf = FALSE
          if (grepl(pattern = ".csv",x = input$to_format$datapath)) {
            dta = read.table(file = input$to_format$datapath,header = header,comment.char = "#",sep = ";")
            if (ncol(dta)==1){
              dta = read.table(file = input$to_format$datapath,header = header,comment.char = "#",sep=",")
            }
          } else {
            dta = read.table(file = input$to_format$datapath,header = header,comment.char = "#")
          }
          colnames(dta) = 1:ncol(dta) # change the colnames for easier use when choosing the number of the column
          formating$data = dta
          formating$colFormat=TRUE
          if (ncol(formating$data)>20){
            formating$head = head(formating$data,n = c(5,20))
          } else {
            formating$head = head(formating$data)
          }
        }
      }
    })
    # output$head = renderTable({
    #   formating$head
    # })
    output$head = renderDataTable({
      datatable(formating$head,rownames = FALSE,options = list(dom = 't'),
                caption = htmltools::tags$caption( style = 'caption-side: top; text-align: center; color:black;  font-size:150% ;','Head dataset to format'))
    })
    output$uiNoColSN = renderUI({
      if (length(formating$data)>0){
        sliderInput(inputId = "col_SN",
                    label = div("Choose the number of the column corresponding to sample names",
                                tipify(el = bsButton(inputId = "10",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "Select the number corresponding to the column with the sample names. Look at the table displayed to be sure.")),
                    min = 1,max = min(20,ncol(formating$data)),value = 1,step = 1)
      }
    })
    output$uiNoColGeno = renderUI({
      if (length(formating$data)>0){
        sliderInput(inputId = "col_geno",
                    label = div("Choose the number of the first column with genotype",
                                tipify(el = bsButton(inputId = "11",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "Select the number corresponding to the first column genotype. Look at the table displayed to be sure.")),
                    min = 1,max = min(20,ncol(formating$data)),value = 2,step = 1)
      }
    })
    output$colFormat <- reactive({
      formating$colFormat
    })
    outputOptions(output, 'colFormat', suspendWhenHidden=FALSE)
    
    #---Load map or txt file with marker names and display the head of df to help select the column
    observeEvent(c(input$snp_map,input$header_snp_map),{
      if (!is.null(input$snp_map) & !is.null(input$header_snp_map)){
        if (input$header_snp_map=='Yes'){
          if (grepl(pattern = ".csv",x = input$snp_map$datapath)){
            tmp1 = read.table(file = input$snp_map$datapath,header = TRUE,comment.char = "#",sep=";")
            tmp2 = read.table(file = input$snp_map$datapath,header = TRUE,comment.char = "#",sep=",")
            if (ncol(tmp1)>ncol(tmp2)){
              formating$dataMap=tmp1
            } else {
              formating$dataMap=tmp2
            }
            rm(tmp1,tmp2)
            gc()
          } else {
            formating$dataMap = read.table(file = input$snp_map$datapath,header = TRUE,comment.char = "#")
          }
        } else {
          if (grepl(pattern = ".csv",x = input$snp_map$datapath)){
            tmp1 = read.table(file = input$snp_map$datapath,header = FALSE,comment.char = "#",sep=";")
            tmp2 = read.table(file = input$snp_map$datapath,header = FALSE,comment.char = "#",sep=",")
            if (ncol(tmp1)>ncol(tmp2)){
              formating$dataMap=tmp1
            } else {
              formating$dataMap=tmp2
            }
            rm(tmp1,tmp2)
            gc()
          } else {
            formating$dataMap = read.table(file = input$snp_map$datapath,header = FALSE,comment.char = "#")
          }
        }
        formating$colMap=TRUE
      }
    })
    # output$headMap = renderTable({
    #   head(formating$dataMap)
    # })
    output$headMap = renderDataTable({
      datatable(head(formating$dataMap),rownames = FALSE,options = list(dom = 't'),
                caption = htmltools::tags$caption( style = 'caption-side: top; text-align: center; color:black;  font-size:150% ;','Head list markers'))
    })
    output$uiNoColMap = renderUI({
      if (length(formating$dataMap)>0){
        sliderInput(inputId = "col_marker",
                    label = div("Choose the number of the column correspoonding to marker names",
                                tipify(el = bsButton(inputId = "12",label = "",icon = icon("question"), style = "info", size = "extra-small"),
                                       title = "Select the number corresponding to the first column genotype. Look at the table displayed to be sure.")),
                    min = 1,max = min(20,ncol(formating$dataMap)),value = 1,step = 1)
      }
    })
    output$colMap <- reactive({
      formating$colMap
    })
    outputOptions(output, 'colMap', suspendWhenHidden=FALSE)
    ###### Launch Formating #####
    observeEvent(input$format,{
      if (!is.null(input$to_format$name)){
        if (!dir.exists("./log")){
          dir.create("./log")
        }
        
        # # Set the saving name
        # indice1 = regexpr(pattern = ".ped",text = input$to_format$name,fixed = TRUE)
        # indice2 = regexpr(pattern = ".txt",text = input$to_format$name,fixed = TRUE)
        indice_vcf = regexpr(pattern = ".vcf",text = input$to_format$name,fixed = TRUE)
        # if (indice1!=-1 | indice2!=-1 | indice_vcf!=-1){
        #   saving_name = substr(x = input$to_format$name,start = 1,stop = nchar(input$to_format$name)-4)
        # } else {
        indice=gregexpr(pattern = ".",text = input$to_format$name,fixed = TRUE)[[1]]
        if (indice[1] !=-1){
          saving_name = substr(x = input$to_format$name,start = 1,stop = indice[length(indice)]-1)
        } else {
          saving_name = input$to_format$name
        }
        # }
        print("-----Formating dataset-----")
        date_time = Sys.time()
        date_time=gsub(pattern = "-",replacement = "",x = date_time)
        date_time=gsub(pattern = ":",replacement = "",x = date_time)
        # date_time=gsub(pattern = " CEST",replacement = "_",x = date_time)
        date_time=gsub(pattern = " ",replacement = "_",x = date_time)
        date_time=substr(x = date_time,start = 3,stop = 16)
        path_log = paste0("./log/",date_time,"_log_Formating.txt")
        write(x = "-----Formating dataset-----",file = path_log)
        if (!dir.exists("./data")){ # result folder
          dir.create("./data")
        }
        write(x = paste0("Dataset : ",input$to_format$name),file = path_log,append = TRUE)
        if (! formating$import_vcf){
          # dta = read.table(file = input$to_format$datapath,header = FALSE,comment.char = "#")
          dta = formating$data
          if(!is.null(input$snp_map$name)){
            write(x = paste0("Marker file : ",input$snp_map$name),file = path_log,append = TRUE)
            marker_name=formating$dataMap[,as.numeric(input$col_marker)]
          } else {
            marker_name=paste0("Marker",seq(1,(ncol(dta)-1)/as.numeric(input$ploidy_format),1))
            write(x = "No marker file provided : markers names are Marker1, Marker2, ...",file = path_log,append = TRUE)
          }
        }
        
        if (formating$import_vcf){
          dta = import_from_vcf(input$to_format$datapath)
          rownames(dta) = toupper(rownames(dta))
          SampleName = rownames(dta)
        } else {
          SampleName = dta[,as.numeric(input$col_SN)]
          indice1=regexpr(pattern="_[A-Z][0-9][0-9].CEL$",SampleName)
          indice2=regexpr(pattern="_[A-Z][0-9].CEL$",SampleName)
          SampleName[indice1!=-1]=substr(x = SampleName[indice1!=-1],start = 1,stop = (indice1[indice1!=-1]-1))
          SampleName[indice2!=-1]=substr(x = SampleName[indice2!=-1],start = 1,stop = (indice2[indice2!=-1]-1))
          SampleName=toupper(SampleName) # au cas ou 'en' en minuscule
          dta=dta[,-c(1:(as.numeric(input$col_geno)-1))]
          dta[dta==0]=NA
          if (ncol(dta)%%as.numeric(input$ploidy_format)!=0){
            stop(paste0("Invalide number of columns ! There should be ",as.numeric(input$ploidy_format)," columns by markers"))
          }
        }
        if (input$what=='Both'){
          write(x = "Parents and offspring are in the same dataset",file = path_log,append = TRUE)
          if (!is.null(input$list_par)){
            write(x = paste0("File with parents names : ",input$list_par$name),file = path_log,append = TRUE)
            par_nam = formating$Lpar[,as.numeric(input$col_head)]
            indice1=regexpr(pattern="_[A-Z][0-9][0-9].CEL$",par_nam)
            indice2=regexpr(pattern="_[A-Z][0-9].CEL$",par_nam)
            par_nam[indice1!=-1]=substr(x = par_nam[indice1!=-1],start = 1,stop = (indice1[indice1!=-1]-1))
            par_nam[indice2!=-1]=substr(x = par_nam[indice2!=-1],start = 1,stop = (indice2[indice2!=-1]-1))
            par_nam=toupper(par_nam) # au cas ou 'en' en minuscule
            
            data_par = dta[SampleName %in% par_nam,]
            data_off = dta[! SampleName %in% par_nam,]
            
            # For parents
            if (! grepl(pattern = ".vcf",x = input$to_format$name,fixed = TRUE)){
              res=Run_formating(data = data_par,
                                SampleName = SampleName[SampleName %in% par_nam],
                                marker_name = marker_name,
                                ploidy = as.numeric(input$ploidy_format),
                                marker_type = input$markerType)
              data_par=res[[1]]
              df_SNP=res[[2]]
            } else {
              allele_freq = as.data.frame(get_allele_frequencies(data_par,ploidy_level = 2))
              
              if (!is.null(allele_freq$Freq_NA)){
                df_SNP = data.frame(MarkerName=rownames(allele_freq),toKeep=TRUE,CR=1-allele_freq$Freq_NA)
              } else {
                df_SNP = data.frame(MarkerName=rownames(allele_freq),toKeep=TRUE,CR=1)
              }
            }
            save(data_par,df_SNP,file=paste0("./data/",saving_name,"_Parents_genoAPIS.Rdata"))
            
            # For offspring
            if (! grepl(pattern = ".vcf",x = input$to_format$name,fixed = TRUE)){
              res=Run_formating(data = data_off,
                                SampleName = SampleName[! SampleName %in% par_nam],
                                marker_name = marker_name,
                                ploidy = as.numeric(input$ploidy_format),
                                marker_type = input$markerType)
              data_off=res[[1]]
              df_SNP=res[[2]]
            } else {
              allele_freq = as.data.frame(get_allele_frequencies(data_off,ploidy_level = as.numeric(input$ploidy)))
              min_non_0 = function(x){
                x[x==min(x[x!=0])][1]
              }
              df_SNP = data.frame(MarkerName=rownames(allele_freq),toKeep=TRUE,MAF=apply(X = allele_freq %>% select(.data$Freq_0,.data$Freq_1),MARGIN = 1,FUN = min_non_0,simplify = TRUE))
            }
            save(data_off,df_SNP,file=paste0("./data/",saving_name,"_Offspring_genoAPIS.Rdata"))
          } else {
            stop("There must be a list with parents names !")
          }
        } else {
          # If only parents or offspring
          if (! grepl(pattern = ".vcf",x = input$to_format$name,fixed = TRUE)){
            res=Run_formating(data = dta,
                              SampleName = SampleName,
                              marker_name = marker_name,
                              ploidy = as.numeric(input$ploidy_format),
                              marker_type = input$markerType)
            new_data=res[[1]]
            df_SNP=res[[2]]
          } else {
            allele_freq = as.data.frame(get_allele_frequencies(dta,ploidy_level = as.numeric(input$ploidy)))
            min_non_0 = function(x){
              x[x==min(x[x!=0])][1]
            }
            if (is.null(allele_freq$Freq_NA)){
              df_SNP = data.frame(MarkerName=rownames(allele_freq),toKeep=TRUE,
                                  MAF=apply(X = allele_freq %>% select(.data$Freq_0,.data$Freq_1),MARGIN = 1,FUN = min_non_0,simplify = TRUE),
                                  CR=1)
            } else {
              df_SNP = data.frame(MarkerName=rownames(allele_freq),toKeep=TRUE,
                                  MAF=apply(X = allele_freq %>% select(.data$Freq_0,.data$Freq_1),MARGIN = 1,FUN = min_non_0,simplify = TRUE),
                                  CR=1-allele_freq$Freq_NA)
            }
            new_data = dta
          }
          save(new_data,df_SNP,file=paste0("./data/",saving_name,"_genoAPIS.Rdata"))
        }
        write(x = "-----End formating-----",file = path_log,append = TRUE)
        write(x = paste0("Saving name : ",saving_name, " (same name as input file but with extension _genoAPIS.Rdata => in ./data folder)"),file = path_log,append = TRUE)
        formating$end="Format OK"
        print("-----End formating-----")
      }
    })
    output$end_format=renderText({
      formating$end
    })
  }
  options(shiny.maxRequestSize=10000*1024^2) # 10000 pour 10Go => augmenter si besoin
  shinyApp(ui = ui, server = server)
}
