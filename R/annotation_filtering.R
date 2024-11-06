#' Compound annotation filtering
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom bsicons bs_icon
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyjqui jqui_resizable
#' @importFrom shinyWidgets materialSwitch
#' @importFrom DT dataTableOutput
#' @noRd


annotation_filter_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    useShinyjs(),
    title = 'Annotation filtering',
    icon = icon("wind"),
    sidebarLayout(
      div(id = ns('Sidebar'),
          #### sidebarPanel ======
          sidebarPanel(width = 2,
                       tags$h3("Addcut for only ms1 matched features",style = 'color: #008080'),
                       hr_main(),
                       HTML(
                         "<span style='color: #7a8788; font-size: 12px; font-style: italic;'>For features ONLY annotated with MS1 database (level 3), the most common adduct are [M+H]+ or [M-H]-. We do not recommend presuming other adducts as highly accurate annotations in the absence of MS/MS matching. Typically, adduct filtering is applied to level 3 annotations here. However, other adducts can be retained based on your preferences. </span>"
                       ),
                       selectInput_div(
                         inputId = ns('af_column'),
                         label = "column type",
                         choices = c("rp","hilic"),
                         selected = "rp",
                         multiple = F,
                         title = "column type"
                       ),
                       selectInput_div(
                         inputId = ns('af_Adduct_pos'),
                         label = "Addcut positive model",
                         choices = "(M+H)+",
                         selected = "(M+H)+",
                         multiple = T,
                         title = "Addcut based on column type"
                       ),
                       selectInput_div(
                         inputId = ns('af_Adduct_neg'),
                         label = "Addcut negative model",
                         choices = "(M-H)-",
                         selected = "(M-H)-",
                         multiple = T,
                         title = "Addcut based on column type"
                       ),
                       tags$h3("Parameters for filtering",style = 'color: #008080'),
                       hr_main(),
                       radioButtons(
                         inputId = ns("af_multi_anno"),
                         label = "Multiple annotations",
                         choices = c("keep all","keep top total score","keep the first one"),
                         selected = "keep the first one"
                       ),
                       radioButtons(
                         inputId = ns("af_redundancy"),
                         label = "Remove redundancy",
                         choices = c("keep all","keep the first one"),
                         selected = "keep all"
                       ),
                       # radioButtons(
                       #   inputId = ns("af_levels"),
                       #   label = "Annotation level",
                       #   choices = c("keep all","keep level 1 and 2 only"),
                       #   selected = "keep all"
                       # ),
                       tags$h3("Parameters for MS/MS match plot",style = 'color: #008080'),
                       hr_main(),
                       radioButtons(
                         inputId = ns("show_mz"),
                         label = "show mz value on plot",
                         choices = c("TRUE","FALSE"),
                         selected = "TRUE"
                       ),
                       radioButtons(
                         inputId = ns("show_detail"),
                         label = "show match details",
                         choices = c("TRUE","FALSE"),
                         selected = "TRUE"
                       ),
                       tags$h3("Export features",style = 'color: #008080'),
                       hr_main(),
                       radioButtons(
                         inputId = ns("feature_remove"),
                         label = "Method",
                         choices = c("Only annotated features","Only features with MS2 spectra","Both","Keep unknown features"),
                         selected = "Both"
                       ),
          )# sidebarPanel
      ),# div
      #### MainPanel ========
      mainPanel(
        fluidPage(
          actionButton(ns('toggleSidebar'),"Toggle sidebar"),
          actionButton(inputId = ns('af_start'),label = "Start annotation filtering",icon = icon("play")),
          hr_bar(),
          materialSwitch(inputId = ns("af_plt_format"),label = "Interactive plot", status = "primary"),
          tabsetPanel(
            tabPanel(
              title = 'Positive',height = '500px',width = "100%",
              icon = icon('plus'),
              tags$h3("Filtered compound annotation",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("af_pos")),
              tags$h3("Status",style = 'color: #008080'),
              hr_main(),
              verbatimTextOutput(ns("object_pos.af")),
              tags$h3("MS/MS match plot",style = 'color: #008080'),
              hr_main(),
              actionButton(inputId = ns("af_pos_show_plot"),label = "show plot",icon = icon("play")),
              hr_head(),
              fluidRow(
                column(width = 6,
                       dataTableOutput(outputId = ns('af_pos_ms2_tbl'))
                ),
                column(width = 6,
                       jqui_resizable(
                         uiOutput(ns("pos_match_mz"))
                       )
                ),
                textInput(inputId = ns("fig1_width"),
                          label = "width",
                          value = 10),
                textInput(inputId = ns("fig1_height"),
                          label = "height",
                          value = 10),
                selectInput(
                  inputId = ns("fig1_format"),label = "format",
                  choices = c("jpg","pdf","png","tiff"),
                  selected = "pdf",selectize = F
                ),
                downloadButton(ns("fig1_download"),"Download"),
              )
            ),
            tabPanel(
              title = 'Negative',height = '500px',width = "100%",
              icon = icon('minus'),
              tags$h3("Filtered compound annotation",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("af_neg")),
              tags$h3("Status",style = 'color: #008080'),
              hr_main(),
              verbatimTextOutput(ns("object_neg.af")),
              tags$h3("MS/MS match plot",style = 'color: #008080'),
              hr_main(),
              actionButton(inputId = ns("af_neg_show_plot"),label = "show plot",icon = icon("play")),
              hr_head(),
              fluidRow(
                column(width = 6,
                       dataTableOutput(outputId = ns('af_neg_ms2_tbl'))
                ),
                column(width = 6,
                       jqui_resizable(
                         uiOutput(ns("neg_match_mz"))
                       ),
                       textInput(inputId = ns("fig2_width"),
                                 label = "width",
                                 value = 10),
                       textInput(inputId = ns("fig2_height"),
                                 label = "height",
                                 value = 10),
                       selectInput(
                         inputId = ns("fig2_format"),label = "format",
                         choices = c("jpg","pdf","png","tiff"),
                         selected = "pdf",selectize = F
                       ),
                       downloadButton(ns("fig2_download"),"Download")
                )
              )
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
#' @importFrom shinyjs toggle runjs useShinyjs
#' @importFrom shinyFiles shinyDirChoose parseDirPath
#' @importFrom dplyr select all_of left_join pull case_when group_by filter mutate slice_head
#' @importFrom stringr str_remove str_detect str_to_lower str_split
#' @importFrom purrr map_chr
#' @importFrom tibble rownames_to_column
#' @importFrom massdataset activate_mass_dataset extract_expression_data extract_sample_info extract_annotation_table
#' @importFrom metid annotate_metabolites_mass_dataset
#' @import MDAtoolkits
#' @param id module of server
#' @param volumes shinyFiles volumes
#' @param prj_init use project init variables.
#' @param data_clean_rv reactivevalues p2 dataclean
#' @noRd


annotation_filter_server <- function(id,volumes,prj_init,data_clean_rv,data_download) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$toggleSidebar, {
      shinyjs::toggle(id = "Sidebar")
    })

    ns <- session$ns

    #> parameters
    ##> download parameters ================
    download_para = reactive({
      list(
        ##> fig1
        fig1_width = as.numeric(input$fig1_width),
        fig1_height = as.numeric(input$fig1_height),
        fig1_format = as.character(input$fig1_format),
        ##> fig2
        fig2_width = as.numeric(input$fig2_width),
        fig2_height = as.numeric(input$fig2_height),
        fig2_format = as.character(input$fig2_format)
      )
    })

    ### 3.6.7 Annotation filtering-----------------------------------------------------

    p2_af_filter <- reactiveValues(data = NULL)

    p2_af_column <- reactive({
      input$af_column %>% as.character()
    })

    observe({
      if (p2_af_column() == "rp") {
        add_pos <- c(
          "(M+H)+","(M+H-H2O)+","(M+H-2H2O)+","(M+NH4)+","(M+Na)+","(M-H+2Na)+","(M-2H+3Na)+","(M+K)+",
          "(M-H+2K)+","(M-2H+3K)+","(M+CH3CN+H)+","(M+CH3CN+Na)+","(2M+H)+","(2M+NH4)+","(2M+Na)+","(M+HCOO+2H)"
        )
        add_neg <- c(
          "(M-H)-","(M-H2O-H)-","(M+Na-2H)-","(M+K-2H)-","(M+NH4-2H)-","(2M-H)-","(M+F)-"
        )
      } else {
        add_pos <- c(
          "(M+H)+","(M+H-H2O)+","(M+H-2H2O)+","(M+NH4)+","(M+Na)+","(M-H+2Na)+","(M-2H+3Na)+",
          "(M+K)+","(M-H+2K)+","(M-2H+3K)+","(M+CH3CN+H)+","(M+CH3CN+Na)+","(2M+H)+","(2M+NH4)+",
          "(2M+Na)+","(2M+K)+","(M+CH3COO+2H)"
        )
        add_neg <- c("(M-H)-","(M-H2O-H)-","(M+Na-2H)-","(M+K-2H)-","(M+NH4-2H)-","(2M-H)-","(M+CH3COO)-")
      }

      updateSelectInput(session, "af_Adduct_pos", choices = add_pos,selected = add_pos[1])
      updateSelectInput(session, "af_Adduct_neg",choices = add_neg,selected = add_neg[1])

    })


    observeEvent(
      input$af_start,
      {
        if(!is.null(prj_init$object_negative.init) & !is.null(prj_init$object_positive.init) & !is.null(prj_init$dblist) & prj_init$steps == "Annotation filtering"){
          p2_af_filter$object_neg.anno= prj_init$object_negative.init %>%
            activate_mass_dataset('sample_info') %>%
            dplyr::select('sample_id') %>%
            dplyr::left_join(prj_init$sample_info)
          p2_af_filter$object_pos.anno = prj_init$object_positive.init %>%
            activate_mass_dataset('sample_info') %>%
            dplyr::select('sample_id') %>%
            dplyr::left_join(prj_init$sample_info)
          p2_af_filter$dblist = prj_init$dblist
        } else {
          if(is.null(data_clean_rv$object_neg.anno)){return()}
          if(is.null(data_clean_rv$object_pos.anno)){return()}
          if(is.null(data_clean_rv$db)){return()}
          p2_af_filter$object_neg.anno = data_clean_rv$object_neg.anno
          p2_af_filter$object_pos.anno = data_clean_rv$object_pos.anno
          p2_af_filter$dblist = data_clean_rv$db
        }
        ##> database

        p2_af_filter$db.name <- purrr::map_chr(.x = 1:length(p2_af_filter$dblist),.f = function(.x) {
          temp_db = p2_af_filter$dblist[[.x]]
          database.name <- paste(temp_db@database.info$Source, temp_db@database.info$Version, sep = "_")
        })

        af_Adduct_pos = input$af_Adduct_pos %>% as.character()
        af_Adduct_neg = input$af_Adduct_neg %>% as.character()

        ##> parameters
        p2_af_filter$af_multi_anno = input$af_multi_anno %>% as.character()
        p2_af_filter$af_redundancy = input$af_redundancy
        # p2_af_filter$af_levels = input$af_levels
        p2_af_filter$feature_remove = input$feature_remove


        ##> reform adduct
        if(length(af_Adduct_neg > 1)) {
          af_Adduct_neg = paste0(af_Adduct_neg,collapse = "|")
        }


        if(length(af_Adduct_pos > 1)) {
          af_Adduct_pos = paste0(af_Adduct_pos,collapse = "|")
        }

        ##> addcut filtring
        ##> pos
        p2_af_filter$object_pos.anno <-
          p2_af_filter$object_pos.anno %>%
          activate_mass_dataset("annotation_table") %>%
          mutate(filter_tag_addcut =
                   dplyr::case_when(
                     Level == 3 & str_detect(Adduct,re_form_reg(af_Adduct_pos)) ~ "retain",
                     Level == 1 | Level == 2 ~ "retain",
                     TRUE ~ "remove"
                   )
          )
        ##> neg
        p2_af_filter$object_neg.anno  <-
          p2_af_filter$object_neg.anno %>%
          activate_mass_dataset("annotation_table") %>%
          mutate(filter_tag_addcut =
                   dplyr::case_when(
                     Level == 3 & str_detect(Adduct,re_form_reg(af_Adduct_neg)) ~ "retain",
                     Level == 1 | Level == 2 ~ "retain",
                     TRUE ~ "remove"
                   )
          )
        ##> add tags
        p2_af_filter$object_pos_temp.af <- p2_af_filter$object_pos.anno
        p2_af_filter$object_neg_temp.af <- p2_af_filter$object_neg.anno

        ##> Annotation cleaning for multi-matched annotation
        if(p2_af_filter$af_multi_anno == "keep top total score") {
          p2_af_filter$object_neg_temp.af =
            p2_af_filter$object_neg_temp.af %>%
            activate_mass_dataset("annotation_table") %>%
            dplyr::group_by(variable_id) %>%
            dplyr::filter(Level == min(Level)) %>%
            dplyr::filter(Total.score == max(Total.score))
          p2_af_filter$object_pos_temp.af =
            p2_af_filter$object_pos_temp.af %>%
            activate_mass_dataset("annotation_table") %>%
            dplyr::group_by(variable_id) %>%
            dplyr::filter(Level == min(Level)) %>%
            dplyr::filter(Total.score == max(Total.score))
        } else if(p2_af_filter$af_multi_anno == "keep the first one") {
          p2_af_filter$object_neg_temp.af =
            p2_af_filter$object_neg_temp.af %>%
            activate_mass_dataset("annotation_table") %>%
            dplyr::group_by(variable_id) %>%
            dplyr::filter(Level == min(Level)) %>%
            dplyr::filter(Total.score == max(Total.score)) %>%
            dplyr::mutate(order = 1:length(variable_id)) %>%
            dplyr::filter(order == 1) %>% dplyr::select(-order)
          p2_af_filter$object_pos_temp.af =
            p2_af_filter$object_pos_temp.af %>%
            activate_mass_dataset("annotation_table") %>%
            dplyr::group_by(variable_id) %>%
            dplyr::filter(Level == min(Level)) %>%
            dplyr::filter(Total.score == max(Total.score)) %>%
            dplyr::mutate(order = 1:length(variable_id)) %>%
            dplyr::filter(order == 1) %>% dplyr::select(-order)
        }
        ##> Remove redundancy
        if(p2_af_filter$af_redundancy == "keep the first one") {
          p2_af_filter$object_neg_temp.af =
            p2_af_filter$object_neg_temp.af %>%
            activate_mass_dataset("annotation_table") %>%
            mutate(Compound.name.fix = str_split(Compound.name,";",Inf,T)[,1] %>% stringr::str_to_lower()) %>%
            group_by(Compound.name.fix) %>%
            slice_head(n = 1)
          p2_af_filter$object_pos_temp.af =
            p2_af_filter$object_pos_temp.af %>%
            activate_mass_dataset("annotation_table") %>%
            mutate(Compound.name.fix = str_split(Compound.name,";",Inf,T)[,1] %>% stringr::str_to_lower()) %>%
            group_by(Compound.name.fix) %>%
            slice_head(n = 1)
        }

        ##> plot
        p2_af_filter$pos_clean_anno =
          p2_af_filter$object_pos_temp.af %>%
          extract_annotation_table() %>%
          filter(filter_tag_addcut == "retain")

        p2_af_filter$neg_clean_anno =
          p2_af_filter$object_neg_temp.af %>%
          extract_annotation_table() %>%
          filter(filter_tag_addcut == "retain")

        p2_af_filter$pos_var_id =
          p2_af_filter$pos_clean_anno %>%
          filter(Level == 2) %>% pull(variable_id) %>% unique()

        p2_af_filter$neg_var_id =
          p2_af_filter$neg_clean_anno %>%
          filter(Level == 2) %>% pull(variable_id) %>% unique()

        output$af_pos = renderDataTable_formated(
          actions = input$af_start,
          condition1 = p2_af_filter$object_pos_temp.af,
          condition2 = p2_af_filter$object_pos.anno,
          filename.a = "3.6.7.AnnoFiltering_pos",
          tbl = p2_af_filter$pos_clean_anno
        )

        p2_af_filter$temp_af_pos_tbl = p2_af_filter$pos_clean_anno %>% dplyr::filter(Level < 3) %>%
          select(variable_id,Compound.name,Lab.ID,Database)

        output$af_pos_ms2_tbl = renderDataTable_formated(
          actions = input$af_start,
          filename.a = "3.6.7.AnnoFiltering_pos_ms2",
          tbl = p2_af_filter$temp_af_pos_tbl
        )

        p2_af_filter$temp_af_neg_tbl = p2_af_filter$neg_clean_anno %>% dplyr::filter(Level < 3) %>%
          select(variable_id,Compound.name,Lab.ID,Database)

        output$af_neg_ms2_tbl = renderDataTable_formated(
          actions = input$af_start,
          filename.a = "3.6.7.AnnoFiltering_neg_ms2",
          tbl = p2_af_filter$temp_af_neg_tbl
        )

        output$af_neg = renderDataTable_formated(
          actions = input$af_start,
          condition1 = p2_af_filter$object_neg_temp.af,
          condition2 = p2_af_filter$object_neg.anno,
          filename.a = "3.6.7.AnnoFiltering_neg",
          tbl = p2_af_filter$neg_clean_anno
        )

        ##> save object

        temp_anno.neg = p2_af_filter$object_neg_temp.af %>%
          extract_annotation_table()

        temp_anno.pos = p2_af_filter$object_pos_temp.af %>%
          extract_annotation_table()

        if(p2_af_filter$feature_remove == "Both") {
          p2_af_filter$object_neg.af = filter_annotations_massdataset(object = p2_af_filter$object_neg.anno,annotate_tbl = temp_anno.neg,method = 'both')
          p2_af_filter$object_pos.af = filter_annotations_massdataset(object = p2_af_filter$object_pos.anno,annotate_tbl = temp_anno.pos,method = 'both')
        } else if(p2_af_filter$feature_remove == "Only features with MS2 spectra") {
          p2_af_filter$object_neg.af = filter_annotations_massdataset(object = p2_af_filter$object_neg.anno,method = 'only ms2')
          p2_af_filter$object_pos.af = filter_annotations_massdataset(object = p2_af_filter$object_pos.anno,method = 'only ms2')
        } else if(p2_af_filter$feature_remove == "Only annotated features") {
          p2_af_filter$object_neg.af = filter_annotations_massdataset(object = p2_af_filter$object_neg.anno,annotate_tbl = temp_anno.neg,method = 'only annotation')
          p2_af_filter$object_pos.af = filter_annotations_massdataset(object = p2_af_filter$object_pos.anno,annotate_tbl = temp_anno.pos,method = 'only annotation')
        } else if(p2_af_filter$feature_remove == "Keep unknown features"){
          p2_af_filter$object_neg.af = p2_af_filter$object_neg_temp.af
          p2_af_filter$object_pos.af = p2_af_filter$object_pos_temp.af
        }

        data_clean_rv$object_neg.af = p2_af_filter$object_neg.af
        data_clean_rv$object_pos.af = p2_af_filter$object_pos.af

        save_massobj(
          polarity = 'positive',
          file_path = paste0(prj_init$wd,"/Result/POS/Objects/"),
          stage = 'af',
          obj = data_clean_rv$object_pos.af)

        save_massobj(
          polarity = 'negative',
          file_path = paste0(prj_init$wd,"/Result/NEG/Objects/"),
          stage = 'af',
          obj = data_clean_rv$object_neg.af)

        ##> status
        output$object_pos.af = renderPrint({
          print(data_clean_rv$object_pos.af)
        })
        output$object_neg.af = renderPrint({
          print(data_clean_rv$object_neg.af)
        })
      }

    )

    observeEvent(input$af_pos_show_plot, {
      tryCatch({
        if (is.null(p2_af_filter$object_pos_temp.af)) {
          return()
        }
        if (is.null(p2_af_filter$pos_clean_anno)) {
          return()
        }

        # plot pos
        show_mz = input$show_mz
        if (show_mz == "TRUE") { show_mz = TRUE } else { show_mz = FALSE }
        show_detail = input$show_detail
        if (show_detail == "TRUE") { show_detail = TRUE } else { show_detail = FALSE }

        # get index
        af_pos_row_idx = input$af_pos_ms2_tbl_rows_selected
        # extract info
        af_pos_row = p2_af_filter$temp_af_pos_tbl[af_pos_row_idx, ]
        p2_af_filter$pos_vari_id = af_pos_row[[1]]
        p2_af_filter$pos_db_name = af_pos_row[[4]]

        temp_idx.pos = match(p2_af_filter$pos_db_name, p2_af_filter$db.name)
        temp_db.pos = p2_af_filter$dblist[[temp_idx.pos]]

        # plot
        p2_af_filter$temp_ms2_match.pos = ms2_plot_mass_dataset_mz(
          object = p2_af_filter$object_pos_temp.af,
          polarity = "positive",
          variable_id = p2_af_filter$pos_vari_id,
          database = temp_db.pos,
          show_mz = show_mz,
          show_detail = show_detail
        )

        # vis
        output$pos_match_mz <- renderUI({
          if (is.null(p2_af_filter$temp_ms2_match.pos)) { return() }
          plot_type <- input$af_plt_format

          if (plot_type) {
            plotlyOutput(outputId = ns("plotly_pos_match_mz"))
          } else {
            plotOutput(outputId = ns("plot_pos_match_mz"))
          }
        })

        output$plot_pos_match_mz <- renderPlot({
          if (is.null(p2_af_filter$temp_ms2_match.pos)) { return() }
          p2_af_filter$temp_ms2_match.pos[[1]]
        })

        output$plotly_pos_match_mz <- renderPlotly({
          if (is.null(p2_af_filter$temp_ms2_match.pos)) { return() }
          plotly::ggplotly(p2_af_filter$temp_ms2_match.pos[[1]])
        })
      }, error = function(e) {
        message("Error occurred: ", e$message)
      })
    })


    observeEvent(input$af_neg_show_plot, {
      tryCatch({
        if (is.null(p2_af_filter$object_neg_temp.af)) {
          return()
        }
        if (is.null(p2_af_filter$neg_clean_anno)) {
          return()
        }

        # plot neg
        show_mz = input$show_mz
        if (show_mz == "TRUE") { show_mz = TRUE } else { show_mz = FALSE }
        show_detail = input$show_detail
        if (show_detail == "TRUE") { show_detail = TRUE } else { show_detail = FALSE }

        af_neg_row_idx = input$af_neg_ms2_tbl_rows_selected
        af_neg_row = p2_af_filter$temp_af_neg_tbl[af_neg_row_idx, ]
        p2_af_filter$neg_vari_id = af_neg_row[[1]]
        p2_af_filter$neg_db_name = af_neg_row[[4]]

        temp_idx.neg = match(p2_af_filter$neg_db_name, p2_af_filter$db.name)
        temp_db.neg = p2_af_filter$dblist[[temp_idx.neg]]

        p2_af_filter$temp_ms2_match.neg = ms2_plot_mass_dataset_mz(
          object = p2_af_filter$object_neg_temp.af,
          polarity = "negative",
          variable_id = p2_af_filter$neg_vari_id,
          database = temp_db.neg,
          show_mz = show_mz,
          show_detail = show_detail
        )

        # mv plot original neg
        output$neg_match_mz <- renderUI({
          plot_type <- input$af_plt_format

          if (plot_type) {
            plotlyOutput(outputId = ns("plotly_neg_match_mz"))
          } else {
            plotOutput(outputId = ns("plot_neg_match_mz"))
          }
        })

        output$plot_neg_match_mz <- renderPlot({
          if (is.null(p2_af_filter$temp_ms2_match.neg)) { return() }
          p2_af_filter$temp_ms2_match.neg[[1]]
        })

        output$plotly_neg_match_mz <- renderPlotly({
          if (is.null(p2_af_filter$temp_ms2_match.neg)) { return() }
          plotly::ggplotly(p2_af_filter$temp_ms2_match.neg[[1]])
        })

      }, error = function(e) {
        message("Error occurred: ", e$message)
      })
    })

    # download ----------------------------------------------------------------
    ###> fig1 =====
    output$fig1_download = downloadHandler(
      filename = function() {
        paste0(p2_af_filter$pos_vari_id,"_ms_ms_mirror_plot.", download_para()$fig1_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # draw condition
        p = p2_af_filter$temp_ms2_match.pos[[1]]
        # save plot
        ggsave(
          filename = file,
          plot = p,
          width = para_d$fig1_width,
          height = para_d$fig1_height,
          device = para_d$fig1_format
        )
      }
    )
    ###> fig2 ====
    output$fig2_download = downloadHandler(
      filename = function() {
        paste0(p2_af_filter$neg_vari_id,"_ms_ms_mirror_plot.", download_para()$fig2_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # draw condition
        p = p2_af_filter$temp_ms2_match.neg

        # save plot

        ggsave(
          filename = file,
          plot = p,
          width = para_d$fig2_width,
          height = para_d$fig2_height,
          device = para_d$fig2_format
        )
      }
    )

  })
}

