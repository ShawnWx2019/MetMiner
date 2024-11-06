#' Compound annotation
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom bsicons bs_icon
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyFiles shinyDirButton shinyFilesButton
#' @importFrom DT dataTableOutput
#' @noRd


annotation_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    useShinyjs(),
    title = 'Annotation',
    icon = icon("stethoscope"),
    sidebarLayout(
      div(id = ns('Sidebar'),
          #### sidebarPanel ======
          sidebarPanel(width = 2,
                       tags$h3("Main parameters",style = 'color: #008080'),
                       hr_main(),
                       textInput_div(
                         inputId = ns('anno_ms1.match.ppm'),
                         label = "ms1.match.ppm ",
                         value = 25,
                         placeholder = "Only accept number, Precursor match ppm tolerance.",
                         title = "Precursor match ppm tolerance."
                       ),
                       textInput_div(
                         inputId = ns('anno_ms2.match.ppm'),
                         label = "ms2.match.ppm ",
                         value = 30,
                         placeholder = "Only accept number, Fragment ion match ppm tolerance.",
                         title = "Fragment ion match ppm tolerance."
                       ),
                       textInput_div(
                         inputId = ns('anno_rt.match.tol'),
                         label = "rt.match.tol",
                         value = 30,
                         placeholder = "Only accept number, RT match tolerance.",
                         title = "RT match tolerance."
                       ),
                       textInput_div(
                         inputId = ns('anno_candidate.num'),
                         label = "candidate.num",
                         value = 3,
                         placeholder = "Only accept number, The number of candidate.",
                         title = "The number of candidate."
                       ),
                       selectInput_div(
                         inputId = ns('anno_column'),
                         label = "column",choices = c("rp","hilic"),
                         selected = "rp",multiple = FALSE,
                         title = "rp: reverse phase \nhilic: HILIC column"
                       ),
                       textInput_div(
                         inputId = ns('anno_threads'),
                         label = "threads",
                         value = 3,
                         placeholder = "Only accept number, The number of threads",
                         title = "Number of threads"
                       ),
                       tags$h3("Public MS1 and MS/MS Database",style = 'color: #008080'),
                       hr_main(),
                       selectInput_div(
                         inputId = ns('norm_db'),
                         label = "Select database",
                         choices = c(
                           "MoNA","Massbank","ReSpect","PlaSMA","MetaboBASE","KEGG","KNApSAcK","Ath_Cyc","Orbitrap","Zma_Cyc"
                         ),selected = c(
                           "MoNA","Massbank","ReSpect","PlaSMA","KEGG","KNApSAcK","Ath_Cyc","MetaboBASE","Orbitrap","Zma_Cyc"
                         ),multiple = T,
                         title = "Select database"
                       ),
                       tags$h3(ns("Customized database"),style = 'color: #008080'),
                       hr_main(),
                       shinyDirButton(id = ns("norm_customized_db"), label = "Choose folder",
                                      title = "Customized database path:",
                                      buttonType = "default", class = NULL,
                                      icon = bs_icon("folder"), multiple = FALSE),
                       tags$span(textOutput(outputId = ns("ms_db_folder_selected")), class = "text-wrap"),
          )# sidebarPanel
      ),# div
      #### MainPanel ========
      mainPanel(
        fluidPage(
          actionButton(ns('toggleSidebar'),"Toggle sidebar"),
          actionButton(inputId = ns('anno_start'),label = "Start annotation",icon = icon("play")),
          hr_bar(),
          tabsetPanel(
            tabPanel(
              title = 'Positive',height = '500px',width = "100%",
              icon = icon('plus'),
              div(id = ns("parameters2"),
                  style = "position: fixed;
                           top: 120px; right: -300px; width: 300px; height: 80%;
                           background-color: #f7f7f7; border: 1.5px solid #008080; overflow-y: scroll;
                           ; border-top-left-radius: 10px; border-bottom-left-radius: 10px; padding: 20px;",
                  h3("Rest parameters for annotation",style = "color:#darkgreen"),
                  hr_main(),
                  textInput_div(
                    inputId = ns('anno_mz.ppm.thr'),label = "mz.ppm.thr",value = 400,title = "Accurate mass tolerance for m/z error calculation.",placeholder = "numeric"
                  ),
                  hr_head(),
                  textInput_div(
                    inputId = ns('anno_ms2.match.tol'),label = 'ms2.match.tol',value = 0.5,placeholder = "numeric",title = "MS2 match (MS2 similarity) tolerance."
                  ),
                  hr_head(),
                  textInput_div(
                    inputId = ns('anno_fraction.weight'),label = 'fraction.weight',value = 0.3,title = "The weight for matched fragments.",placeholder = "numeric"
                  ),
                  hr_head(),
                  textInput_div(
                    inputId = ns('anno_dp.forward.weight'),label = 'dp.forward.weight',value = 0.6,title = "Forward dot product weight.",placeholder = "numeric"
                  ),
                  hr_head(),
                  textInput_div(
                    inputId = ns('anno_dp.reverse.weight'),label = 'dp.reverse.weight',value = 0.1,title = "Reverse dot product weight.",placeholder = "numeric"
                  ),
                  hr_head(),
                  textInput_div(
                    inputId = ns('anno_remove_fragment_intensity_cutoff'),label = 'remove_fragment_intensity_cutoff',value = 0,title = "remove_fragment_intensity_cutoff",placeholder = "numeric"
                  ),
                  hr_head(),
                  textInput_div(
                    inputId = ns('anno_ce'),label = 'ce',value = "all",title = "Collision energy. Please confirm the CE values in your database. Default is all",placeholder = "CE model"
                  ),
                  hr_head(),
                  textInput_div(
                    inputId = ns('anno_ms1.match.weight'),label = 'ms1.match.weight',value = 0.25,title = "The weight of MS1 match for total score calculation.",placeholder = "numeric"
                  ),
                  hr_head(),
                  textInput_div(
                    inputId = ns('anno_rt.match.weight'),label = 'rt.match.weight',value = 0.25,title = "The weight of RT match for total score calculation.",placeholder = "numeric"
                  ),
                  hr_head(),
                  textInput_div(
                    inputId = ns('anno_ms2.match.weight'),label = 'ms2.match.weight',value = 0.5,title = "The weight of MS2 match for total score calculation.",placeholder = "numeric"
                  ),
                  hr_head(),
                  textInput_div(
                    inputId = ns('anno_total.score.tol'),label = 'total.score.tol',value = 0.5,title = "Total score tolerance. The total score are referring to MS-DIAL.",placeholder = "numeric"
                  )
              ),
              tags$h3("Parameter settings",style = 'color: #008080'),
              hr_main(),
              div(actionButton(ns("para_select2"), "Show | Hide",icon = icon('computer-mouse')),
                  style = "margin-bottom: 15px;"),
              tags$h3("Summary",style = 'color: #008080'),
              hr_main(),
              htmlOutput(ns("anno_check1_pos")),
              tags$h3("Annotation table",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("Pos_anno")),
              tags$h3("Status",style = 'color: #008080'),
              hr_main(),
              verbatimTextOutput(ns("object_pos.anno")),
              tags$script("
                        function toggleParameters2() {
                          var div = document.getElementById('parameters2');
                          var right = div.style.right;
                          if (right === '-300px') {
                            div.style.right = '0px';
                            div.style.transition = 'right 0.8s ease-in-out';
                            return;
                          } else if (right === '0px') {
                            div.style.right = '-300px';
                            div.style.transition = 'right 0.8s ease-in-out';
                            return;
                          }
                        }
                      ")
            ),
            tabPanel(
              title = 'Negative',height = '500px',width = "100%",
              icon = icon('minus'),
              tags$h3(ns("Annotation table"),style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("Neg_anno")),
              tags$h3("Status",style = 'color: #008080'),
              hr_main(),
              verbatimTextOutput(ns("object_neg.anno"))
            )
          )
        )
      )
    )
  )


}


#' Compound annotation
#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs toggle runjs
#' @importFrom dplyr select left_join
#' @importFrom stringr str_remove
#' @importFrom tibble rownames_to_column
#' @importFrom massdataset activate_mass_dataset extract_annotation_table
#' @importFrom metid annotate_metabolites_mass_dataset
#' @import MDAtoolkits
#' @param id module of server
#' @param volumes shinyFiles volumes
#' @param prj_init use project init variables.
#' @param data_clean_rv reactivevalues p2 dataclean
#' @noRd


annotation_server <- function(id,volumes,prj_init,data_clean_rv) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    observeEvent(input$toggleSidebar, {
      shinyjs::toggle(id = "Sidebar")
    })

    #> parameters - for compound annotation parameters.
    observeEvent(input$para_select2,{
      shinyjs::runjs(sprintf("toggleParameters('%s')", ns("parameters2")))
    })

    p2_anno <- reactiveValues(data = NULL)



    #> MS database dir
    observe({
      shinyDirChoose(input = input,
                     id = "norm_customized_db",
                     roots =  volumes,
                     session = session)
      if(!is.null(input$norm_customized_db)){
        # browser()
        ms_db_folder_selected<-parseDirPath(roots = volumes, input$norm_customized_db)
        output$ms_db_folder_selected <- renderText(ms_db_folder_selected)
      }})

    observeEvent(
      input$anno_start,
      {
        if(!is.null(prj_init$object_negative.init) & !is.null(prj_init$object_positive.init) & prj_init$steps == "Annotation"){
          p2_anno$object_neg.norm= prj_init$object_negative.init %>%
            activate_mass_dataset('sample_info') %>%
            dplyr::select('sample_id') %>%
            dplyr::left_join(prj_init$sample_info)
          p2_anno$object_pos.norm = prj_init$object_positive.init %>%
            activate_mass_dataset('sample_info') %>%
            dplyr::select('sample_id') %>%
            dplyr::left_join(prj_init$sample_info)
        } else {
          if(is.null(data_clean_rv$object_neg.norm)){return()}
          if(is.null(data_clean_rv$object_pos.norm)){return()}
          p2_anno$object_neg.norm = data_clean_rv$object_neg.norm
          p2_anno$object_pos.norm = data_clean_rv$object_pos.norm
        }


        ## buildin database
        p2_anno$buildin_db <-
          list(
            MoNA = MDAtoolkits::mona_database0.0.4,
            Massbank = MDAtoolkits::massbank_database0.0.4,
            ReSpect = MDAtoolkits::respect_database0.0.1,
            PlaSMA = MDAtoolkits::plasma_database0.0.1,
            KEGG = MDAtoolkits::kegg_plant_database0.0.1,
            KNApSAcK = MDAtoolkits::knapsack_agri_database0.0.1,
            Ath_Cyc = MDAtoolkits::ath_plantcyc.database0.0.1,
            Zma_Cyc = MDAtoolkits::zma_plantcyc_database0.0.1,
            MetaboBASE = MDAtoolkits::metabobase_database0.0.1,
            Orbitrap = MDAtoolkits::orbitrap_database0.0.3
          )
        ##
        p2_anno$buildin_name = input$norm_db %>% as.character()

        if(length(p2_anno$buildin_name) == 0) {
          p2_anno$buildin_db = NULL
        } else {
          temp_anno_idx = match(p2_anno$buildin_name,names(p2_anno$buildin_db))
          p2_anno$buildin_db = p2_anno$buildin_db[temp_anno_idx]
        }

        ## Customized ms database
        p2_anno$norm_customized_db <- parseDirPath(volumes, input$norm_customized_db)
        p2_anno$cuz_db_path <- p2_anno$norm_customized_db %>% as.character()
        temp_file_name = dir(p2_anno$cuz_db_path,"*.rda")

        if(length(temp_file_name) == 0) {
          p2_anno$db = p2_anno$buildin_db
        } else {
          p2_anno$cuz_db = list()
          for (i in 1:length(temp_file_name)) {
            xx = load(file = paste0(p2_anno$cuz_db_path,"/",temp_file_name[[i]]))
            p2_anno$cuz_db[[i]] = get(xx)
          }
          p2_anno$cuz_name = str_remove(string = temp_file_name,pattern = "\\.rda")

          names(p2_anno$cuz_db) = p2_anno$cuz_name

          if(is.null(p2_anno$buildin_db)){
            p2_anno$db = p2_anno$cuz_db
          } else {
            p2_anno$db <- c(p2_anno$buildin_db,p2_anno$cuz_db)
          }
        }
        dir.create(path = paste0(prj_init$wd,"/Result/Database/"),showWarnings = F,recursive = T)
        temp_db <- p2_anno$db
        data_clean_rv$db <- p2_anno$db
        save(temp_db,file =  paste0(prj_init$wd,"/Result/Database/auto_saved.dblist"))

        ## parameters

        p2_anno$ms1.match.ppm = input$anno_ms1.match.ppm %>% as.numeric()
        p2_anno$ms2.match.ppm = input$anno_ms2.match.ppm %>% as.numeric()
        p2_anno$rt.match.tol = input$anno_rt.match.tol %>% as.numeric()
        p2_anno$candidate.num = input$anno_candidate.num %>% as.numeric()
        p2_anno$column = input$anno_column %>% as.character()
        p2_anno$threads = input$anno_threads %>% as.numeric()

        p2_anno$mz.ppm.thr = input$anno_mz.ppm.thr %>% as.numeric()
        p2_anno$ms2.match.tol = input$anno_ms2.match.tol %>% as.numeric()
        p2_anno$fraction.weight = input$anno_fraction.weight %>% as.numeric()
        p2_anno$dp.forward.weight = input$anno_dp.forward.weight %>% as.numeric()
        p2_anno$dp.reverse.weight = input$anno_dp.reverse.weight %>% as.numeric()
        p2_anno$remove_fragment_intensity_cutoff = input$anno_remove_fragment_intensity_cutoff %>% as.numeric()
        p2_anno$ce = input$anno_ce %>% as.character()
        p2_anno$ms1.match.weight = input$anno_ms1.match.weight %>% as.numeric()
        p2_anno$rt.match.weight = input$anno_rt.match.weight %>% as.numeric()
        p2_anno$ms2.match.weight = input$anno_ms2.match.weight %>% as.numeric()
        p2_anno$total.score.tol = input$anno_total.score.tol %>% as.numeric()

        if(length(p2_anno$db) == 0) {
          output$anno_check1_pos = renderUI({
            isolate(HTML(paste0(
              '<font color = red>"No database were selected!"<b></font></b>',
            )))
          })
        } else {
          tags = names(p2_anno$db)
          ##> compound annotation
          pro_steps_anno = c(paste0("Database ",tags," in progress..."),"Finish!")

          anno_steps = length(pro_steps_anno)
          withProgress(message = 'Compound annoation', value = 0,
                       expr = {
                         for (i in 1:(anno_steps)) {
                           incProgress(1/anno_steps,detail = pro_steps_anno[i])
                           if(i == 1) {
                             p2_anno$object_neg.anno = annotate_metabolites_mass_dataset(
                               object = p2_anno$object_neg.norm,polarity = "negative",
                               database = p2_anno$db[[i]] ,
                               ms1.match.ppm = p2_anno$ms1.match.ppm,
                               ms2.match.ppm = p2_anno$ms2.match.ppm,
                               rt.match.tol = p2_anno$rt.match.tol,
                               candidate.num = p2_anno$candidate.num,
                               column = p2_anno$column,
                               threads = p2_anno$threads,
                               mz.ppm.thr = p2_anno$mz.ppm.thr,
                               ms2.match.tol = p2_anno$ms2.match.tol,
                               fraction.weight = p2_anno$fraction.weight,
                               dp.forward.weight = p2_anno$dp.forward.weight,
                               dp.reverse.weight = p2_anno$dp.reverse.weight,
                               remove_fragment_intensity_cutoff = p2_anno$remove_fragment_intensity_cutoff,
                               ce = p2_anno$ce,
                               ms1.match.weight = p2_anno$ms1.match.weight,
                               rt.match.weight = p2_anno$rt.match.weight,
                               ms2.match.weight = p2_anno$ms2.match.weight,
                               total.score.tol = p2_anno$total.score.tol
                             )
                             p2_anno$object_pos.anno = annotate_metabolites_mass_dataset(
                               object = p2_anno$object_pos.norm,polarity = "positive",
                               database = p2_anno$db[[i]] ,
                               ms1.match.ppm = p2_anno$ms1.match.ppm,
                               ms2.match.ppm = p2_anno$ms2.match.ppm,
                               rt.match.tol = p2_anno$rt.match.tol,
                               candidate.num = p2_anno$candidate.num,
                               column = p2_anno$column,
                               threads = p2_anno$threads,
                               mz.ppm.thr = p2_anno$mz.ppm.thr,
                               ms2.match.tol = p2_anno$ms2.match.tol,
                               fraction.weight = p2_anno$fraction.weight,
                               dp.forward.weight = p2_anno$dp.forward.weight,
                               dp.reverse.weight = p2_anno$dp.reverse.weight,
                               remove_fragment_intensity_cutoff = p2_anno$remove_fragment_intensity_cutoff,
                               ce = p2_anno$ce,
                               ms1.match.weight = p2_anno$ms1.match.weight,
                               rt.match.weight = p2_anno$rt.match.weight,
                               ms2.match.weight = p2_anno$ms2.match.weight,
                               total.score.tol = p2_anno$total.score.tol
                             )
                           } else if(i > 1 & i < anno_steps) {
                             p2_anno$object_neg.anno = annotate_metabolites_mass_dataset(
                               object = p2_anno$object_neg.anno,polarity = "negative",
                               database = p2_anno$db[[i]] ,
                               ms1.match.ppm = p2_anno$ms1.match.ppm,
                               ms2.match.ppm = p2_anno$ms2.match.ppm,
                               rt.match.tol = p2_anno$rt.match.tol,
                               candidate.num = p2_anno$candidate.num,
                               column = p2_anno$column,
                               threads = p2_anno$threads,
                               mz.ppm.thr = p2_anno$mz.ppm.thr,
                               ms2.match.tol = p2_anno$ms2.match.tol,
                               fraction.weight = p2_anno$fraction.weight,
                               dp.forward.weight = p2_anno$dp.forward.weight,
                               dp.reverse.weight = p2_anno$dp.reverse.weight,
                               remove_fragment_intensity_cutoff = p2_anno$remove_fragment_intensity_cutoff,
                               ce = p2_anno$ce,
                               ms1.match.weight = p2_anno$ms1.match.weight,
                               rt.match.weight = p2_anno$rt.match.weight,
                               ms2.match.weight = p2_anno$ms2.match.weight,
                               total.score.tol = p2_anno$total.score.tol
                             )
                             p2_anno$object_pos.anno = annotate_metabolites_mass_dataset(
                               object = p2_anno$object_pos.anno,polarity = "positive",
                               database = p2_anno$db[[i]] ,
                               ms1.match.ppm = p2_anno$ms1.match.ppm,
                               ms2.match.ppm = p2_anno$ms2.match.ppm,
                               rt.match.tol = p2_anno$rt.match.tol,
                               candidate.num = p2_anno$candidate.num,
                               column = p2_anno$column,
                               threads = p2_anno$threads,
                               mz.ppm.thr = p2_anno$mz.ppm.thr,
                               ms2.match.tol = p2_anno$ms2.match.tol,
                               fraction.weight = p2_anno$fraction.weight,
                               dp.forward.weight = p2_anno$dp.forward.weight,
                               dp.reverse.weight = p2_anno$dp.reverse.weight,
                               remove_fragment_intensity_cutoff = p2_anno$remove_fragment_intensity_cutoff,
                               ce = p2_anno$ce,
                               ms1.match.weight = p2_anno$ms1.match.weight,
                               rt.match.weight = p2_anno$rt.match.weight,
                               ms2.match.weight = p2_anno$ms2.match.weight,
                               total.score.tol = p2_anno$total.score.tol
                             )
                           } else if(i == anno_steps) {
                             Sys.sleep(2)
                           }
                         }
                       }

          )
          data_clean_rv$object_pos.anno = p2_anno$object_pos.anno
          data_clean_rv$object_neg.anno = p2_anno$object_neg.anno
          save_massobj(
            polarity = 'positive',
            file_path = paste0(prj_init$wd,"/Result/POS/Objects/"),
            stage = 'anno',
            obj = p2_anno$object_pos.anno)

          save_massobj(
            polarity = 'negative',
            file_path = paste0(prj_init$wd,"/Result/NEG/Objects/"),
            stage = 'anno',
            obj = p2_anno$object_neg.anno)


          #> information of mass datasets
          output$object_pos.anno = renderPrint({
            print(p2_anno$object_pos.anno)
          })
          output$object_neg.anno = renderPrint({
            print(p2_anno$object_neg.anno)
          })

          #> data table
          #>
          output$Pos_anno = renderDataTable_formated(
            actions = input$anno_start,
            condition1 = p2_anno$object_pos.anno,filename.a = "3.6.6.annotation_pos",
            tbl = p2_anno$object_pos.anno %>% extract_annotation_table()
          )

          output$Neg_anno = renderDataTable_formated(
            actions = input$anno_start,
            condition1 = p2_anno$object_neg.anno,filename.a = "3.6.6.annotation_neg",
            tbl = p2_anno$object_neg.anno %>% extract_annotation_table()
          )


          #> Summary
          temp_db_name = paste(tags,collapse = " | ")
          output$anno_check1_pos = renderUI({
            isolate(HTML(paste0(
              '<font color = blue> <b>Selected database: </b> </font> <font color=red>',temp_db_name,'</font> ')))
          })
        }
      }
    )

  })
}

