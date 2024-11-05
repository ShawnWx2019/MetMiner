#' class
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyjqui jqui_resizable
#' @importFrom shinyWidgets materialSwitch
#' @noRd


data_class_ui <- function(id) {
  ns <- NS(id)
  ### Part1.4.2 Remove noisy features -----------------------------------------

  tabPanel(
    useShinyjs(),
    title = 'Classification',
    icon = icon("list"),
    sidebarLayout(
      div(id = ns('Sidebar'),
          #### sidebarPanel ======
          sidebarPanel(width = 2,
                       tags$h3("Parameters",style = 'color: #008080'),
                       hr_main(),
                       radioButtons(
                         inputId = ns('class_method'),
                         label = 'Method',
                         choices = c("from new upload annotation file","from mass_dataset"),
                         selected = "from mass_dataset"
                       ),
                       hr_head(),
                       fileInput(
                         inputId = ns('class_annotation_file'),
                         label = "upload annotation file (option)",
                         accept = ".csv"
                       ),
                       hr_head(),
                       radioButtons(
                         inputId = ns('class_rest'),
                         label = 'Search unclassified metabolites',choices = c("TRUE","FALSE"),
                         selected = 'TRUE'
                       ),
                       tags$h3("Pie plot",style = 'color: #008080'),
                       hr_main(),
                       radioButtons(
                         inputId = ns('pie_tag'),
                         label = 'choise',choices = c("superclass","class","subclass"),
                         selected = 'class'
                       ),
                       hr_head(),
                       textInput(
                         inputId = ns('pie_cut'),
                         label = "other group size",
                         value = 10
                       ),
                       hr_head(),
                       textInput(
                         inputId = ns('pie_ntop'),
                         label = "Show top x main groups",
                         value = 15
                       )
          )# sidebarPanel
      ),# div
      #### MainPanel ========
      mainPanel(
        fluidPage(
          actionButton(ns('toggleSidebar'),"Toggle sidebar"),
          actionButton(inputId = ns('class_start'),label = "Start",icon = icon("play")),
          actionButton(inputId = ns('class_show'),"Show pie plot",icon = icon("angles-up")),
          hr_head(),
          materialSwitch(inputId = ns("class_plt_format"),label = "Interactive plot", status = "primary"),
          tabsetPanel(
            tabPanel(
              title = 'Table',height = '500px',width = "100%",
              icon = icon('chart'),
              tags$h3("Classified",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("class_table")),
              tags$h3("Unclassified",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("unclass_table")),
              hr_head(),
              actionButton(inputId = ns('class_rest_try'),label = "Search for the rest",icon = icon("play")),
              tags$h3("New classified",style = 'color: #008080'),
              DT::dataTableOutput(ns("class_table_new")),
              hr_main(),
            ),

            tabPanel(
              title = 'Plot',height = '500px',width = "100%",
              icon = icon('image'),
              tags$h3("pie plot",style = 'color: #008080'),
              hr_main(),
              jqui_resizable(
                uiOutput(ns("class_pie_plot"),fill = T)
              ),
            )
          )
        )
      )
    )
  )


}


#' classification
#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs toggle runjs useShinyjs
#' @importFrom dplyr select all_of left_join pull filter rename
#' @importFrom tibble rownames_to_column
#' @importFrom massdataset activate_mass_dataset extract_expression_data extract_sample_info extract_annotation_table
#' @importFrom plotly renderPlotly plotlyOutput
#' @import MDAtoolkits
#' @param id module of server
#' @param volumes shinyFiles volumes
#' @param prj_init use project init variables.
#' @param data_clean_rv reactivevalues p2 dataclean
#' @noRd


data_class_server <- function(id,volumes,prj_init,data_clean_rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$toggleSidebar, {
      shinyjs::toggle(id = "Sidebar")
    })

    ### 3.6.9 classification ---------------------------------------------------------
    p3_class <- reactiveValues(data = NULL)
    #> variable table

    temp_class_anno <- reactive({
      file1 <- input$class_annotation_file
      if(is.null(file1)){return()}
      read.csv(file = file1$datapath,
               sep=",",header = T,stringsAsFactors = F)
    })

    #> File check event
    observeEvent(
      input$class_start,
      {
        if (!is.null(temp_class_anno())) {
          temp_anno_file = temp_class_anno()
        } else if(!is.null(prj_init$object_positive.init) & prj_init$steps == "DAM and rest"){
          p3_class$object_merge= prj_init$object_positive.init;
        } else if (!is.null(data_clean_rv$object_merge)) {
          p3_class$object_merge = data_clean_rv$object_merge
        } else {return()}

        ## para
        temp_method =input$class_method %>% as.character()
        p3_class$rest_method = input$class_rest %>% as.logical()

        if(temp_method == "from new upload annotation file") {
          temp_upload_file = temp_anno_file %>% select(1:3) %>% setNames(c("variable_id","Compound.name","Lab.ID"))
        } else {
          temp_upload_file = p3_class$object_merge %>% extract_annotation_table() %>%
            dplyr::filter(filter_tag_addcut == 'retain') %>%
            select(variable_id,Compound.name,Lab.ID)
        }

        ## add classification info
        p3_class$class_inner <- temp_upload_file %>% left_join(class.database,by = "Lab.ID")
        p3_class$class_inner_clean <- p3_class$class_inner %>% filter(!is.na(superclass))
        p3_class$temp_rest <- p3_class$class_inner %>% filter(is.na(superclass)) %>%
          select(variable_id,Compound.name,Lab.ID)

        output$class_table = renderDataTable_formated(
          actions = input$class_start,
          condition1 = p3_class$class_inner_clean,
          filename.a = 'classyfire_1st.csv',
          tbl = p3_class$class_inner_clean
        )

        output$unclass_table = renderDataTable_formated(
          actions = input$class_start,
          condition1 = p3_class$temp_rest,
          filename.a = 'unclassyfire_1st.csv',
          tbl = p3_class$temp_rest
        )

      }
    )


    observeEvent(
      input$class_rest_try,
      {
        if(isTRUE(p3_class$rest_method) & !is.null(p3_class$temp_rest)){
          query = p3_class$temp_rest %>% pull(Compound.name) %>% unique()
          tags = c("Get InChIKey via PUG REST","Get classification via CBF","Get KEGG.ID via KEGG REST")
          ##> compound annotation
          pro_steps_class = c(paste0(tags," in progress..."),"Finish!")

          class_steps = length(pro_steps_class)
          withProgress(message = 'Compound classification', value = 0,
                       expr = {
                         for (i in 1:(class_steps)) {
                           incProgress(1/class_steps,detail = pro_steps_class[i])
                           if(i == 1 ) {
                             temp_inchikey = MDAtoolkits::mda_get_cid_fast(query = query,probe = 'name',core_num = 1)
                           } else if (i == 2) {
                             query2 = temp_inchikey %>% pull(InChIKey) %>% unique()
                             temp_class = MDAtoolkits::cfb_crawler(query = query2,delay_max = .5,ssl.verifypeer = F)
                           } else if(i == 3) {
                             temp_keggid = MDAtoolkits::mda_name2kegg(query = query,type = 'multiple',core_num = 1)
                           } else if(i == 4){
                             p3_class$class_new = p3_class$temp_rest %>%
                               left_join(temp_inchikey %>% dplyr::rename('Compound.name' = 'query')) %>%
                               left_join(temp_class,by = "InChIKey") %>%
                               left_join(temp_keggid,by = "Compound.name")
                           }
                         }})

          output$class_table_new = renderDataTable_formated(
            actions = input$class_start,
            condition1 = p3_class$class_new,
            filename.a = 'new_classyfire.csv',
            tbl = p3_class$class_new
          )

          temp_new = p3_class$class_new %>% select(colnames( p3_class$class_inner_clean))
          p3_class$class_inner_clean = rbind(p3_class$class_inner_clean,temp_new)
          ##> update class and KEGG result
          data_clean_rv$object_merge =
            data_clean_rv$object_merge %>%
            activate_mass_dataset('annotation_table') %>%
            select(-all_of(c("superclass","class","subclass","parent_levels","description","KEGG.ID"))) %>%
            dplyr::left_join(p3_class$class_inner_clean,by = 'variable_id')
          temp_filepath = paste0(prj_init$wd,"/Result/Data_Export/")
          object_merge = data_clean_rv$object_merge
          save(object_merge,file = paste0(temp_filepath,"object_merge.rda"))
          output$class_table = renderDataTable_formated(
            actions = input$class_start,
            condition1 = p3_class$class_inner_clean,
            filename.a = 'classyfire_1st.csv',
            tbl = p3_class$class_inner_clean
          )

        }
      }
    )


    observeEvent(
      input$class_show,
      {
        if(is.null(p3_class$class_inner_clean)) {return()}

        temp_pie_tag = input$pie_tag %>% as.character()
        temp_pie_cut = input$pie_cut %>% as.numeric()
        temp_pie_ntop = input$pie_ntop %>% as.numeric()
        p3_class$plot_pie = pie_plot(x = p3_class$class_inner_clean,tag = temp_pie_tag,cut = temp_pie_cut,ntop = temp_pie_ntop)
        p3_class$plotly_pie = pie_plot_plotly(x = p3_class$class_inner_clean,tag = temp_pie_tag,cut = temp_pie_cut,ntop = temp_pie_ntop)
        output$class_pie_plot <- renderUI({
          plot_type <- input$class_plt_format
          if (plot_type) {
            plotlyOutput(outputId = ns("pie_plotly"))
          } else {
            plotOutput(outputId = ns("pie_plot"))
          }
        })

        output$pie_plot <- renderPlot({
          if(is.null(p3_class$plot_pie)){return()}
          p3_class$plot_pie
        })

        output$pie_plotly <- renderPlotly({
          if(is.null(p3_class$plotly_pie)){return()}
          p3_class$plotly_pie
        })
      })

  })
}

