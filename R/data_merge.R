#' Merge ESI+ and ESI-
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyjqui jqui_resizable
#' @importFrom DT dataTableOutput
#' @importFrom shinyWidgets materialSwitch
#' @noRd


data_merge_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    useShinyjs(),
    title = 'Data integration',
    icon = icon("object-group"),
    sidebarLayout(
      div(id = ns('Sidebar'),
          #### sidebarPanel ======
          sidebarPanel(width = 2,
                       tags$h3("Merge method",style = 'color: #008080'),
                       hr_main(),
                       radioButtons(
                         inputId = ns("data_merge_s_direction"),
                         label = "sample direction",
                         choices = c("left","right","inner","full"),
                         selected = "inner"
                       ),
                       radioButtons(
                         inputId = ns("data_merge_v_direction"),
                         label = "variable direction",
                         choices = c("left","right","inner","full"),
                         selected = "full"
                       ),
          )# sidebarPanel
      ),# div
      #### MainPanel ========
      mainPanel(
        fluidPage(
          actionButton(ns('toggleSidebar'),"Toggle sidebar"),
          actionButton(inputId = ns('data_merge_show'),label = "Start intergration",icon = icon("play")),
          actionButton(inputId = ns('data_merge_vis'),label = "Start visualization",icon = icon("feather")),
          actionButton(inputId = ns('data_merge_export'),label = "Export data",icon = icon("file-export")),
          hr_bar(),
          materialSwitch(inputId = ns("merge_plt_format"),label = "Interactive plot", status = "primary"),
          tabsetPanel(
            tabPanel(
              title = 'Tables',height = '500px',width = "100%",
              icon = icon('table'),
              tags$h3("Accumulation matrix",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("data_merge_expmat")),
              tags$h3("Sample information",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("data_merge_sample_info")),
              tags$h3("variable information",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("data_merge_variable_info")),
              tags$h3("Compound annotation",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("data_merge_annotation"))
            ),
            tabPanel(
              title = 'Figures',height = '500px',width = "100%",
              icon = icon('image'),
              tags$h3("PCA",style = 'color: #008080'),
              hr_main(),
              selectInput_div(inputId = ns("merge_colorby"),label = "colored by",choices = c("group",'batch'),selected = "group",multiple = F,title = "PCA color by"),
              jqui_resizable(
                uiOutput(ns("merge_pca"),fill = T)
              ),
              hr_head(),
              textInput(inputId = ns("width4.9.1"),
                        label = "width",
                        value = 10),
              textInput(inputId = ns("height4.9.1"),
                        label = "height",
                        value = 10),
              actionButton(ns("adjust4.9.1"),"Set fig size"),
              downloadButton(ns("downfig4.9.1"),"Download"),
              tags$h3("Correlation",style = 'color: #008080'),
              hr_main(),
              fluidRow(
                column(width = 4,
                       sliderInput(
                         inputId = ns("merge_break_min"),step = 0.1,
                         label = "break min",min = -1,max = 1,value = -1
                       ),
                       sliderInput(
                         inputId = ns("merge_break_mid"),step = 0.1,
                         label = "break min",min = -1,max = 1,value = 0
                       ),
                       sliderInput(
                         inputId = ns("merge_break_max"),step = 0.1,
                         label = "break min",min = -1,max = 1,value = 1
                       )

                ),
                column(width = 4,

                       colourpicker::colourInput(
                         inputId = ns("merge_min"),
                         label = "corr min",
                         value = "blue",
                       ),
                       colourpicker::colourInput(
                         inputId = ns("merge_mid"),
                         label = "corr mid",
                         value = "white",
                       ),
                       colourpicker::colourInput(
                         inputId = ns("merge_max"),
                         label = "corr max",
                         value = "yellow",
                       )
                ),
                column(width = 4,
                       selectInput(
                         inputId = ns("merge_showname"),
                         label = "show names",
                         choices = c("TRUE","FALSE"),
                         selected = "TRUE",multiple = F
                       ),
                       selectInput(
                         inputId = ns("merge_cluster"),
                         label = "show cluster",
                         choices = c("TRUE","FALSE"),
                         selected = "TRUE",multiple = F
                       ),
                       hr_head(),
                       selectInput_div(inputId = ns("merge_anno_key"),label = "colored by",choices = c("group",'batch'),selected = "group",multiple = F,title = "b_anno by"),
                       textAreaInput(
                         inputId = ns("merge_anno_color"),
                         label = "bottom annotation color",resize = "both",placeholder = NULL,
                         value = NULL
                       )
                )
              ),
              hr_head(),
              jqui_resizable(
                plotOutput(ns("merge_corr"))
              ),
            )
          )
        )
      )
    )
  )
}


#' Merge ESI+ and ESI-
#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs toggle runjs
#' @importFrom dplyr select left_join filter inner_join rename pull
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom massdataset activate_mass_dataset extract_annotation_table extract_sample_info extract_variable_info extract_expression_data export_mass_dataset4metdna export_ms2_data
#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom stringr str_split
#' @import MDAtoolkits
#' @param id module of server
#' @param volumes shinyFiles volumes
#' @param prj_init use project init variables.
#' @param data_clean_rv reactivevalues p2 dataclean
#' @noRd


data_merge_server <- function(id,volumes,prj_init,data_clean_rv) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    observeEvent(input$toggleSidebar, {
      shinyjs::toggle(id = "Sidebar")
    })
    p2_data_merge <- reactiveValues(data = NULL)

    observeEvent(
      input$data_merge_show,
      {
        if(!is.null(prj_init$object_negative.init) & !is.null(prj_init$object_positive.init) & prj_init$steps == "Data integrate"){
          p2_data_merge$object_neg.af= prj_init$object_negative.init %>%
            activate_mass_dataset('sample_info') %>%
            dplyr::select('sample_id') %>%
            dplyr::left_join(prj_init$sample_info)
          p2_data_merge$object_pos.af = prj_init$object_positive.init %>%
            activate_mass_dataset('sample_info') %>%
            dplyr::select('sample_id') %>%
            dplyr::left_join(prj_init$sample_info)
        } else {
          if(is.null(data_clean_rv$object_neg.af)){return()}
          if(is.null(data_clean_rv$object_pos.af)){return()}
          p2_data_merge$object_neg.af = data_clean_rv$object_neg.af
          p2_data_merge$object_pos.af = data_clean_rv$object_pos.af
        }
        p2_data_merge$object_pos.mrm = p2_data_merge$object_pos.af
        p2_data_merge$object_neg.mrm = p2_data_merge$object_neg.af
        p2_data_merge$object_neg.af =
          p2_data_merge$object_neg.af %>%
          activate_mass_dataset("sample_info") %>%
          filter(class != "QC")

        p2_data_merge$object_pos.af =
          p2_data_merge$object_pos.af %>%
          activate_mass_dataset("sample_info") %>%
          filter(class != "QC")
        ##> merge massdataset
        s_direction = input$data_merge_s_direction %>% as.character()
        v_direction = input$data_merge_v_direction %>% as.character()

        p2_data_merge$object_merge <- merge_mass_dataset_fix(
          x = p2_data_merge$object_neg.af,
          y = p2_data_merge$object_pos.af,
          sample_direction = s_direction,
          variable_direction = v_direction
        )
        ##> add inbuild classyFire database
        p2_data_merge$object_merge =
          p2_data_merge$object_merge %>% activate_mass_dataset('annotation_table') %>%
          left_join(MDAtoolkits::class.database)

        ##> correct sample info and variable_info
        ##> temp table
        neg_sample_info <-
          p2_data_merge$object_neg.af %>%
          extract_sample_info()
        pos_sample_info <-
          p2_data_merge$object_pos.af %>%
          extract_sample_info()

        neg_vari_info <-
          p2_data_merge$object_neg.af %>%
          extract_variable_info()
        pos_vari_info <-
          p2_data_merge$object_pos.af %>%
          extract_variable_info()

        ##> Output tables
        ##> sample info
        p2_data_merge$sample_info_full <-
          rbind(neg_sample_info,pos_sample_info) %>%
          unique() %>%
          inner_join(p2_data_merge$object_merge %>%
                       extract_sample_info() %>%
                       dplyr::select(sample_id))

        output$data_merge_sample_info = renderDataTable_formated(
          actions = input$data_merge_show,
          condition1 = p2_data_merge$sample_info_full,
          filename.a = "3.6.8.Merge_sample_info",
          tbl = p2_data_merge$sample_info_full
        )

        ##> variable info
        p2_data_merge$vari_info_full <-
          rbind(neg_vari_info,pos_vari_info)%>%
          inner_join(p2_data_merge$object_merge %>%
                       extract_variable_info() %>%
                       dplyr::select(variable_id))

        output$data_merge_variable_info = renderDataTable_formated(
          actions = input$data_merge_show,
          condition1 = p2_data_merge$vari_info_full,
          filename.a = "3.6.8.Merge_variable_info",
          tbl = p2_data_merge$vari_info_full
        )

        ##> annotation info
        p2_data_merge$anno_table <-
          p2_data_merge$object_merge %>%
          extract_annotation_table() %>% dplyr::filter(filter_tag_addcut == 'retain')

        output$data_merge_annotation = renderDataTable_formated(
          actions = input$data_merge_show,
          condition1 = p2_data_merge$anno_table,
          filename.a = "3.6.8.Merge_anno_table",
          tbl = p2_data_merge$anno_table
        )

        ##> expmat info
        p2_data_merge$expmat<-
          p2_data_merge$object_merge %>%
          extract_expression_data() %>%
          rownames_to_column("variable_id")

        output$data_merge_expmat= renderDataTable_formated(
          actions = input$data_merge_show,
          condition1 = p2_data_merge$expmat,
          filename.a = "3.6.8.Merge_expmat",
          tbl = p2_data_merge$expmat
        )

      }

    )

    observe({
      updateSelectInput(session, "merge_colorby",choices = colnames(p2_data_merge$sample_info_full),selected = "group")
      updateSelectInput(session, "merge_anno_key",choices = colnames(p2_data_merge$sample_info_full),selected = "group")
    })


    observeEvent(
      input$data_merge_vis,
      {
        temp_tag = input$merge_colorby %>% as.character()
        output$merge_pca <- renderUI({
          plot_type <- input$merge_plt_format
          if (plot_type) {
            plotlyOutput(outputId = ns("plotly_pca"))
          } else {
            plotOutput(outputId = ns("plot_pca"))
          }
        })

        output$plot_pca <- renderPlot({
          if(is.null(p2_data_merge$object_merge)){return()}
          pca_plot(
            object = p2_data_merge$object_merge,
            colby = temp_tag,
            center = T,
            scale = T,
            removeVar = .1,interactive = F
          )
        })

        output$plotly_pca <- renderPlotly({

          if(is.null(p2_data_merge$object_merge)){return()}
          pca_plot(
            object = p2_data_merge$object_merge,
            colby = temp_tag,
            center = T,
            scale = T,
            removeVar = .1,interactive = T
          )

        })
        ##> correlation plot breaks and colors
        temp_break_min = input$merge_break_min %>% as.numeric()
        temp_break_mid = input$merge_break_mid %>% as.numeric()
        temp_break_max = input$merge_break_max %>% as.numeric()
        temp_break = c(temp_break_min,temp_break_mid,temp_break_max)
        temp_min = input$merge_min %>% as.character()
        temp_mid = input$merge_mid %>% as.character()
        temp_max = input$merge_max %>% as.character()
        temp_col = c(temp_min,temp_mid,temp_max)
        col_fun = circlize::colorRamp2(breaks = temp_break,colors = temp_col)
        ##> show name and do cluster
        temp_show_name = input$merge_showname %>% as.logical()
        temp_cluster = input$merge_cluster %>% as.logical()
        ##> add b_annotation
        temp_lab_key = input$merge_anno_key %>% as.character()
        temp_lab_color = input$merge_anno_color %>% as.character() %>%
          stringr::str_split(pattern  = "\\\n|,| ",n = Inf,simplify = T) %>% unlist()
        print(temp_lab_color)
        ##> cor mat
        temp_cor = cor(p2_data_merge$expmat %>% column_to_rownames("variable_id"))

        ##> generate heatmap
        if(is.null(temp_lab_key)){
          ht_cor = ComplexHeatmap::Heatmap(
            matrix = temp_cor,col = col_fun,name = "r",border = T,show_row_names = temp_show_name,show_column_names = temp_show_name,
            cluster_rows = temp_cluster,cluster_columns = temp_cluster
          )
        } else {
          ##> extract tags
          temp_lab_tags = p2_data_merge$sample_info_full %>% dplyr::rename("tag" = temp_lab_key) %>% pull(tag) %>% unique()
          ##> match col and tag
          if(is.null(temp_lab_color)) {
            b_anno = HeatmapAnnotation(
              sample_tag = p2_data_merge$sample_info_full %>% select(sample_id,temp_lab_key) %>% column_to_rownames("sample_id") %>% as.matrix()
            )
          } else if(length(temp_lab_color) == length(temp_lab_tags)) {
            b_anno = HeatmapAnnotation(
              sample_tag = p2_data_merge$sample_info_full %>% select(sample_id,temp_lab_key) %>% column_to_rownames("sample_id") %>% as.matrix(),
              col = list(
                sample_tag = temp_lab_color %>% setNames(temp_lab_tags)
              )
            )
          } else {
            b_anno = HeatmapAnnotation(
              sample_tag = p2_data_merge$sample_info_full %>% select(sample_id,temp_lab_key) %>% column_to_rownames("sample_id") %>% as.matrix()
            )
          }
          ##> plot
          ht_cor = ComplexHeatmap::Heatmap(
            matrix = temp_cor,col = col_fun,name = "r",border = T,show_row_names = temp_show_name,show_column_names = temp_show_name,
            cluster_rows = temp_cluster,cluster_columns = temp_cluster,bottom_annotation = b_anno
          )
        }
        #> draw ht
        p2_data_merge$ht_cor = draw(ht_cor)

        output$merge_corr <- renderPlot({
          if(is.null(p2_data_merge$ht_cor)){return()}
          p2_data_merge$ht_cor
        })

        #      InteractiveComplexHeatmapWidget(input, output, session, p2_data_merge$ht_cor,output_id = "merge_corr_ht")
      }

    )

    observeEvent(input$data_merge_export,{
      if(is.null(p2_data_merge$object_merge)){return()}
      if(is.null(p2_data_merge$sample_info_full)){return()}
      if(is.null(p2_data_merge$vari_info_full)){return()}
      if(is.null(p2_data_merge$anno_table)){return()}
      if(is.null(p2_data_merge$expmat)){return()}

      dir.create(path = paste0(prj_init$wd,"/Result/Data_Export/"),showWarnings = F,recursive = T)
      temp_filepath = paste0(prj_init$wd,"/Result/Data_Export/")
      object_merge <- p2_data_merge$object_merge
      data_clean_rv$object_merge = p2_data_merge$object_merge
      object_pos = p2_data_merge$object_pos.af
      object_neg = p2_data_merge$object_neg.af
      print('check point mrm')
      temp_pos_mrm = p2_data_merge$object_pos.mrm %>%  MDAtoolkits::oneStepMRMselection() %>%
        extract_variable_info() %>%
        dplyr::select(variable_id,rt,precursor,product) %>%
        dplyr::mutate(
          Compound = variable_id,
          `Retention Time (min)` = rt,
          `RT Window (min)` = 1,
          Polarity = "Positive",
          `Precursor (m/z)` = precursor,
          `Product (m/z)` = product
        ) %>%
        dplyr::filter(!is.na(precursor)) %>%
        dplyr::mutate(
          `Collision Energy (V)` = dplyr::case_when(
            `Precursor (m/z)` < 200 ~ 15,
            `Precursor (m/z)` >= 200 & `Precursor (m/z)` < 400 ~ 20,
            `Precursor (m/z)` >= 400 & `Precursor (m/z)` < 600 ~ 25,
            TRUE ~ 30
          ),
          `Min Dwell Time (ms)` = 1
        ) %>%
        dplyr::select(Compound,`Retention Time (min)`,`RT Window (min)`,Polarity,`Precursor (m/z)`,`Product (m/z)`,
                      `Collision Energy (V)`,`Min Dwell Time (ms)`)
      temp_neg_mrm =  p2_data_merge$object_neg.mrm %>% MDAtoolkits::oneStepMRMselection() %>%
        extract_variable_info() %>%
        dplyr::select(variable_id,rt,precursor,product) %>%
        dplyr::mutate(
          Compound = variable_id,
          `Retention Time (min)` = rt,
          `RT Window (min)` = 1,
          Polarity = "Negative",
          `Precursor (m/z)` = precursor,
          `Product (m/z)` = product
        ) %>%
        dplyr::filter(!is.na(precursor)) %>%
        dplyr::mutate(
          `Collision Energy (V)` = dplyr::case_when(
            `Precursor (m/z)` < 200 ~ 15,
            `Precursor (m/z)` >= 200 & `Precursor (m/z)` < 400 ~ 20,
            `Precursor (m/z)` >= 400 & `Precursor (m/z)` < 600 ~ 25,
            TRUE ~ 30
          ),
          `Min Dwell Time (ms)` = 1
        ) %>%
        dplyr::select(Compound,`Retention Time (min)`,`RT Window (min)`,Polarity,`Precursor (m/z)`,`Product (m/z)`,
                      `Collision Energy (V)`,`Min Dwell Time (ms)`)


      mrm_selection = rbind(temp_pos_mrm,temp_neg_mrm)

      temp_datalist = list(
        sample_info = p2_data_merge$sample_info_full,
        variable_info = p2_data_merge$vari_info_full,
        annotation_table = p2_data_merge$anno_table,
        expmat = p2_data_merge$expmat,
        mrm_table = mrm_selection
      )
      tags = c("mass_dataset.rda","table.xlsx","ms2_file.mgf","MRM_selection.xlsx")
      merge_steps_tag = c(paste0("Save ",tags," in progress..."),"Finish!")

      merge_steps = length(merge_steps_tag)
      withProgress(message = 'Export Intergrated data', value = 0,
                   expr = {
                     for (i in 1:(merge_steps)) {
                       incProgress(1/merge_steps,detail = merge_steps_tag[i])
                       if(i == 1){
                         save(object_merge,file = paste0(temp_filepath,"object_merge.rda"))
                       } else if(i == 2) {
                         writexl::write_xlsx(x = temp_datalist,path = paste0(temp_filepath,"table_merge.xlsx"))
                       } else if(i == 3) {
                         dir.create(paste0(temp_filepath,"/MS2_match/NEG/"),showWarnings = F,recursive = T)
                         dir.create(paste0(temp_filepath,"/MS2_match/POS/"),showWarnings = F,recursive = T)
                         dir.create(paste0(temp_filepath,"/MetDNA/POS/"),showWarnings = F,recursive = T)
                         dir.create(paste0(temp_filepath,"/MetDNA/NEG/"),showWarnings = F,recursive = T)
                         export_ms2_data(object = object_pos,file_type = "mgf",path = paste0(temp_filepath,"/MS2_match/POS/"))
                         export_ms2_data(object = object_neg,file_type = "mgf",path = paste0(temp_filepath,"/MS2_match/NEG/"))
                         export_mass_dataset4metdna(object = object_pos,path = paste0(temp_filepath,"/MetDNA/POS/"))
                         export_mass_dataset4metdna(object = object_neg,path = paste0(temp_filepath,"/MetDNA/NEG/"))
                       } else {
                         Sys.sleep(1)
                       }
                     }
                   })

    })

  })
}

