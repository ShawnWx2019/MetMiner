#' data normalization
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyjqui jqui_resizable
#' @importFrom shinyWidgets materialSwitch
#' @importFrom DT dataTableOutput
#' @noRd


data_normalize_ui <- function(id) {
  ns <- NS(id)
  ### Part1.4.2 Remove noisy features -----------------------------------------

  tabPanel(
    useShinyjs(),
    title = 'Normalization',
    icon = icon("align-justify"),
    sidebarLayout(
      div(id = ns('Sidebar'),
          #### sidebarPanel ======
          sidebarPanel(width = 2,
                       tags$h3("parameters",style = 'color: #008080'),
                       hr_main(),
                       selectInput(
                         inputId = ns('norm_method'),
                         label = "method ",multiple = F,
                         choices = c("svr", "total", "median", "mean", "pqn", "loess",
                                     "ppca"),
                         selected = 'svr'
                       ),
                       radioButtons(
                         inputId = ns('norm_keep_scale'),
                         label = "keep_scale",choices = c("TRUE","FALSE"),selected = "TRUE"
                       ),
                       radioButtons(
                         inputId = ns('norm_optimization'),
                         label = "optimization",choices = c("TRUE","FALSE"),selected = "TRUE"
                       ),
                       radioButtons(
                         inputId = ns('norm_pqn_reference'),
                         label = "pqn_reference",choices = c("median","mean"),selected = "median"
                       ),
                       textInput(
                         inputId = ns('norm_begin'),
                         label = "begin",
                         value = 0.5
                       ),
                       textInput(
                         inputId = ns('norm_end'),
                         label = "end",
                         value = 1
                       ),
                       textInput(
                         inputId = ns('norm_step'),
                         label = "step",
                         value = 0.2
                       ),
                       textInput(
                         inputId = ns('norm_multiple'),
                         label = "multiple",
                         value = 1
                       ),
                       textInput(
                         inputId = ns('norm_threads'),
                         label = "threads",
                         value = 1
                       ),
                       tags$h3("Visualize parameter",style = 'color: #008080'),
                       hr_main(),
                       selectInput(
                         inputId = ns('pca_col_by'),
                         label = "pca colored by",
                         choices = "class",selected = "class",multiple = F
                       )
          )# sidebarPanel
      ),# div
      #### MainPanel ========
      mainPanel(
        fluidPage(
          actionButton(ns('toggleSidebar'),"Toggle sidebar"),
          actionButton(inputId = ns('norm_start'),label = "Start normalization",icon = icon("play")),
          actionButton(inputId = ns('norm_plot_start'),label = "Visualize",icon = icon("images")),
          hr_bar(),
          materialSwitch(inputId = ns("norm_plt_format"),label = "Interactive plot", status = "primary"),
          tabsetPanel(
            tabPanel(
              title = 'Positive',height = '500px',width = "100%",
              icon = icon('plus'),
              tags$h3("Normalized accumulation table",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("pos_norm_tbl")),
              tags$h3("PCA normalized",style = 'color: #008080'),
              hr_main(),
              jqui_resizable(
                uiOutput(ns("pca.pos_plt"),fill = T)
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
              tags$h3("PCA raw",style = 'color: #008080'),
              hr_main(),
              jqui_resizable(
                uiOutput(ns("pca.pos_plt_raw"),fill = T)
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
              downloadButton(ns("fig2_download"),"Download"),
              tags$h3("RSD",style = 'color: #008080'),
              hr_main(),
              jqui_resizable(
                uiOutput(ns("rsd.pos_plt"))
              ),
              textInput(inputId = ns("fig3_width"),
                        label = "width",
                        value = 10),
              textInput(inputId = ns("fig3_height"),
                        label = "height",
                        value = 10),
              selectInput(
                inputId = ns("fig3_format"),label = "format",
                choices = c("jpg","pdf","png","tiff"),
                selected = "pdf",selectize = F
              ),
              downloadButton(ns("fig3_download"),"Download"),
              tags$h3("Status",style = 'color: #008080'),
              hr_main(),
              verbatimTextOutput(ns("object_pos.norm"))
            ),
            tabPanel(
              title = 'Negative',height = '500px',width = "100%",
              icon = icon('minus'),
              tags$h3("Normalized accumulation table",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("neg_norm_tbl")),
              tags$h3("PCA normalized",style = 'color: #008080'),
              hr_main(),
              jqui_resizable(
                uiOutput(ns("pca.neg_plt"),fill = T)
              ),

              textInput(inputId = ns("fig4_width"),
                        label = "width",
                        value = 10),
              textInput(inputId = ns("fig4_height"),
                        label = "height",
                        value = 10),
              selectInput(
                inputId = ns("fig4_format"),label = "format",
                choices = c("jpg","pdf","png","tiff"),
                selected = "pdf",selectize = F
              ),
              downloadButton(ns("fig4_download"),"Download"),
              tags$h3("PCA raw",style = 'color: #008080'),
              hr_main(),
              jqui_resizable(
                uiOutput(ns("pca.neg_plt_raw"),fill = T)
              ),

              textInput(inputId = ns("fig5_width"),
                        label = "width",
                        value = 10),
              textInput(inputId = ns("fig5_height"),
                        label = "height",
                        value = 10),
              selectInput(
                inputId = ns("fig5_format"),label = "format",
                choices = c("jpg","pdf","png","tiff"),
                selected = "pdf",selectize = F
              ),
              downloadButton(ns("fig5_download"),"Download"),
              tags$h3("RSD",style = 'color: #008080'),
              hr_main(),
              jqui_resizable(
                uiOutput(ns("rsd.neg_plt"),fill = T)
              ),
              textInput(inputId = ns("fig6_width"),
                        label = "width",
                        value = 10),
              textInput(inputId = ns("fig6_height"),
                        label = "height",
                        value = 10),
              selectInput(
                inputId = ns("fig6_format"),label = "format",
                choices = c("jpg","pdf","png","tiff"),
                selected = "pdf",selectize = F
              ),
              downloadButton(ns("fig6_download"),"Download"),
              tags$h3("Status",style = 'color: #008080'),
              hr_main(),
              verbatimTextOutput(ns("object_neg.norm"))
            )
          )
        )
      )
    )
  )

}


#' normalization
#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs toggle runjs useShinyjs
#' @importFrom dplyr select all_of left_join pull
#' @importFrom tibble rownames_to_column
#' @importFrom massdataset activate_mass_dataset extract_expression_data extract_sample_info
#' @importFrom masscleaner normalize_data integrate_data
#' @importFrom plotly renderPlotly plotlyOutput
#' @param id module of server
#' @param volumes shinyFiles volumes
#' @param prj_init use project init variables.
#' @param data_clean_rv reactivevalues p2 dataclean
#' @noRd


data_normalize_server <- function(id,volumes,prj_init,data_clean_rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    ### 3.6.5 normalization -----------------------------------------------------
    observeEvent(input$toggleSidebar, {
      shinyjs::toggle(id = "Sidebar")
    })

    p2_norm <- reactiveValues(data = NULL)


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
        fig2_format = as.character(input$fig2_format),
        ##> fig3
        fig3_width = as.numeric(input$fig3_width),
        fig3_height = as.numeric(input$fig3_height),
        fig3_format = as.character(input$fig3_format),
        ##> fig4
        fig4_width = as.numeric(input$fig4_width),
        fig4_height = as.numeric(input$fig4_height),
        fig4_format = as.character(input$fig4_format),
        ##> fig5
        fig5_width = as.numeric(input$fig5_width),
        fig5_height = as.numeric(input$fig5_height),
        fig5_format = as.character(input$fig5_format),
        ##> fig6
        fig6_width = as.numeric(input$fig6_width),
        fig6_height = as.numeric(input$fig6_height),
        fig6_format = as.character(input$fig6_format)
      )
    })

    observe({
      updateSelectInput(session, "pca_col_by",choices = colnames(prj_init$sample_info),selected = colnames(prj_init$sample_info)[3])
    })

    observeEvent(
      input$norm_start,
      {
        if(!is.null(prj_init$object_negative.init) & !is.null(prj_init$object_positive.init) & prj_init$steps == "Normalization"){
          p2_norm$object_neg.impute= prj_init$object_negative.init %>%
            activate_mass_dataset('sample_info') %>%
            dplyr::select('sample_id') %>%
            dplyr::left_join(prj_init$sample_info)
          p2_norm$object_pos.impute = prj_init$object_positive.init %>%
            activate_mass_dataset('sample_info') %>%
            dplyr::select('sample_id') %>%
            dplyr::left_join(prj_init$sample_info)
        } else {
          if(is.null(data_clean_rv$object_pos.impute)) {return()}
          if(is.null(data_clean_rv$object_neg.impute)) {return()}
          p2_norm$object_neg.impute = data_clean_rv$object_neg.impute
          p2_norm$object_pos.impute = data_clean_rv$object_pos.impute
        }

        #> parameters

        temp_method = input$norm_method %>% as.character()
        temp_norm_keep_scale = input$norm_keep_scale %>% as.character()
        temp_norm_optimization = input$norm_optimization %>% as.character()
        temp_norm_pqn_reference = input$norm_pqn_reference %>% as.character()

        temp_norm_begin = input$norm_begin %>% as.numeric()
        temp_norm_end = input$norm_end %>% as.numeric()
        temp_norm_step = input$norm_step %>% as.numeric()
        temp_norm_multiple = input$norm_multiple %>% as.numeric()
        temp_norm_threads = input$norm_threads %>% as.numeric()

        #> acc tbl

        if(temp_norm_keep_scale == "TRUE") {temp_norm_keep_scale = TRUE} else {temp_norm_keep_scale = FALSE}
        if(temp_norm_optimization == "TRUE") {temp_norm_optimization = TRUE} else {temp_norm_optimization = FALSE}

        #> data normalization
        #> positive model
        #>
        pro_step_tbl = c(
          'Positive',
          'Negative',
          'All finish'
        )

        #functions
        withProgress(message = 'Data normalization', value = 0,
                     expr = {
                       for (i in 1:3) {
                         incProgress(1/3,detail = pro_step_tbl[i])
                         if(i == 1) {
                           p2_norm$object_pos.norm <-
                             p2_norm$object_pos.impute %>%
                             normalize_data(
                               method = temp_method,
                               keep_scale = temp_norm_keep_scale,
                               optimization = temp_norm_optimization,
                               pqn_reference = temp_norm_pqn_reference,
                               begin = temp_norm_begin,
                               end = temp_norm_end,
                               step = temp_norm_step,
                               multiple = temp_norm_multiple,
                               threads = temp_norm_threads
                             )
                           if((p2_norm$object_pos.norm %>% extract_sample_info() %>% pull(batch) %>% unique() %>% length()) >1)
                           { p2_norm$object_pos.norm = integrate_data(object = p2_norm$object_pos.norm,method = "qc_mean")}
                         } else if(i == 2){
                           p2_norm$object_neg.norm <-
                             p2_norm$object_neg.impute %>%
                             normalize_data(
                               method = temp_method,
                               keep_scale = temp_norm_keep_scale,
                               optimization = temp_norm_optimization,
                               pqn_reference = temp_norm_pqn_reference,
                               begin = temp_norm_begin,
                               end = temp_norm_end,
                               step = temp_norm_step,
                               multiple = temp_norm_multiple,
                               threads = temp_norm_threads
                             )
                           if((p2_norm$object_neg.norm %>% extract_sample_info() %>% pull(batch) %>% unique() %>% length()) >1)
                           { p2_norm$object_neg.norm = integrate_data(object = p2_norm$object_neg.norm,method = "qc_mean")}
                         } else {Sys.sleep(time = 1.5)}
                       }})
        ##> add rsd info
        ##> export data

        output$pos_norm_tbl = renderDataTable_formated(
          actions = input$norm_start,
          condition1 = p2_norm$object_pos.norm,
          filename.a = "3.6.5.Normalization_Acc_Mat_pos",
          tbl = p2_norm$object_pos.norm %>% extract_expression_data() %>% rownames_to_column("variable_id")
        )
        output$neg_norm_tbl = renderDataTable_formated(
          actions = input$norm_start,condition1 = p2_norm$object_neg.norm,
          filename.a = "3.6.5.Normalization_Acc_Mat_neg",
          tbl = p2_norm$object_neg.norm %>% extract_expression_data() %>% rownames_to_column("variable_id")
        )

        temp_pos_rsd = rsd_plot(obj_old = p2_norm$object_pos.impute,obj_new = p2_norm$object_pos.norm,QC_tag = "QC")
        p2_norm$pos.plt = temp_pos_rsd$plot
        p2_norm$pos.tbl = temp_pos_rsd$rsd_tbl

        temp_neg_rsd = rsd_plot(obj_old = p2_norm$object_neg.impute,obj_new = p2_norm$object_neg.norm,QC_tag = "QC")
        p2_norm$neg.plt = temp_neg_rsd$plot
        p2_norm$neg.tbl = temp_neg_rsd$rsd_tbl

        p2_norm$object_pos.norm <-
          p2_norm$object_pos.norm %>%
          activate_mass_dataset('variable_info') %>%
          left_join(p2_norm$pos.tbl,by = c('variable_id' = 'ID'))


        p2_norm$object_neg.norm <-
          p2_norm$object_neg.norm %>%
          activate_mass_dataset('variable_info') %>%
          left_join(p2_norm$neg.tbl,by = c('variable_id' = 'ID'))
        data_clean_rv$object_pos.norm = p2_norm$object_pos.norm
        data_clean_rv$object_neg.norm = p2_norm$object_neg.norm

        #> save mass object
        save_massobj(
          polarity = 'positive',
          file_path = paste0(prj_init$wd,"/Result/POS/Objects/"),
          stage = 'norm',
          obj = p2_norm$object_pos.norm)

        save_massobj(
          polarity = 'negative',
          file_path = paste0(prj_init$wd,"/Result/NEG/Objects/"),
          stage = 'norm',
          obj = p2_norm$object_neg.norm)


        #> information of mass datasets
        output$object_pos.norm = renderPrint({
          print(p2_norm$object_pos.norm)
        })
        output$object_neg.norm = renderPrint({
          print(p2_norm$object_neg.norm)
        })

      }
    )

    ##> plot
    observeEvent(
      input$norm_plot_start,
      {
        if(is.null(p2_norm$object_pos.norm)) {return()}
        if(is.null(p2_norm$object_neg.norm)) {return()}
        p2_norm$temp_norm_pca_col_by = input$pca_col_by %>% as.character()
        #> mv plot original pos

        output$pca.pos_plt <- renderUI({
          plot_type <- input$norm_plt_format
          if (plot_type) {
            plotlyOutput(outputId = ns("plotly_norm_pca.pos"))
          } else {
            plotOutput(outputId = ns("plot_norm_pca.pos"))
          }
        })
        output$pca.pos_plt_raw <- renderUI({
          plot_type <- input$norm_plt_format
          if (plot_type) {
            plotlyOutput(outputId = ns("plotly_norm_pca2.pos"))
          } else {
            plotOutput(outputId = ns("plot_norm_pca2.pos"))
          }
        })

        output$plot_norm_pca.pos <- renderPlot({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_pos.norm)){return()}
          pca_plot(object = p2_norm$object_pos.norm,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = F)
        })
        output$plot_norm_pca2.pos <- renderPlot({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_pos.norm)){return()}
          pca_plot(object = p2_norm$object_pos.impute,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = F)
        })
        output$plotly_norm_pca.pos <- renderPlotly({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_pos.norm)){return()}
          pca_plot(object = p2_norm$object_pos.norm,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = T)
        })
        output$plotly_norm_pca2.pos <- renderPlotly({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_pos.norm)){return()}
          pca_plot(object = p2_norm$object_pos.impute,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = T)
        })

        #> mv plot original neg
        output$pca.neg_plt <- renderUI({
          plot_type <- input$norm_plt_format

          if (plot_type) {
            plotlyOutput(outputId = ns("plotly_norm_pca.neg"))
          } else {
            plotOutput(outputId = ns("plot_norm_pca.neg"))
          }
        })

        output$pca.neg_plt_raw <- renderUI({
          plot_type <- input$norm_plt_format

          if (plot_type) {
            plotlyOutput(outputId = ns("plotly_norm_pca2.neg"))
          } else {
            plotOutput(outputId = ns("plot_norm_pca2.neg"))
          }
        })


        output$plot_norm_pca.neg <- renderPlot({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_neg.norm)){return()}
          pca_plot(object = p2_norm$object_neg.norm,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = F)
        })

        output$plot_norm_pca2.neg <- renderPlot({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_neg.norm)){return()}
          pca_plot(object = p2_norm$object_neg.impute,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = F)
        })

        output$plotly_norm_pca.neg <- renderPlotly({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_neg.norm)){return()}
          pca_plot(object = p2_norm$object_neg.norm,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = T)
        })

        output$plotly_norm_pca2.neg <- renderPlotly({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_neg.norm)){return()}
          pca_plot(object = p2_norm$object_neg.impute,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = T)
        })

        #> rsd plot
        #> mv plot original pos
        output$rsd.pos_plt <- renderUI({
          plot_type <- input$norm_plt_format
          if (plot_type) {
            plotlyOutput(outputId = ns("plotly_norm_rsd.pos"))
          } else {
            plotOutput(outputId = ns("plot_norm_rsd.pos"))
          }
        })

        output$plot_norm_rsd.pos <- renderPlot({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_pos.norm)){return()}
          p2_norm$pos.plt
        })
        output$plotly_norm_rsd.pos <- renderPlotly({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_pos.norm)){return()}
          p2_norm$pos.plt %>% plotly::ggplotly()
        })

        #> mv plot original neg
        output$rsd.neg_plt <- renderUI({
          plot_type <- input$norm_plt_format
          if (plot_type) {
            plotlyOutput(outputId = ns("plotly_norm_rsd.neg"))
          } else {
            plotOutput(outputId = ns("plot_norm_rsd.neg"))
          }
        })


        output$plot_norm_rsd.neg <- renderPlot({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_neg.norm)){return()}
          p2_norm$neg.plt
        })

        output$plotly_norm_rsd.neg <- renderPlotly({
          if(is.null(input$norm_plot_start)){return()}
          if(is.null(p2_norm$object_neg.norm)){return()}
          p2_norm$neg.plt %>% plotly::ggplotly()
        })
      }
    )



    # download ----------------------------------------------------------------
    ###> fig1 =====
    output$fig1_download = downloadHandler(
      filename = function() {
        paste0("01.pca_normalization-pos.", download_para()$fig1_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # draw condition
        p = pca_plot(object = p2_norm$object_pos.norm,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = F)

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
        paste0("02.pca_raw-pos.", download_para()$fig2_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # draw condition
        p = pca_plot(object = p2_norm$object_pos.impute,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = F)

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
    ###> fig3 =====
    output$fig3_download = downloadHandler(
      filename = function() {
        paste0("03.rsd_plot-pos.", download_para()$fig3_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # draw condition
        p = p2_norm$pos.plt

        # save plot
        ggsave(
          filename = file,
          plot = p,
          width = para_d$fig3_width,
          height = para_d$fig3_height,
          device = para_d$fig3_format
        )
      }
    )
    ###> fig4 ====
    output$fig4_download = downloadHandler(
      filename = function() {
        paste0("01.pca_normalization-neg.", download_para()$fig4_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # draw condition
        p = pca_plot(object = p2_norm$object_neg.norm,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = F)

        # save plot
        ggsave(
          filename = file,
          plot = p,
          width = para_d$fig4_width,
          height = para_d$fig4_height,
          device = para_d$fig4_format
        )
      }
    )
    ###> fig5 ====
    output$fig5_download = downloadHandler(
      filename = function() {
        paste0("02.pca_raw-neg.", download_para()$fig5_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # draw condition
        p = pca_plot(object = p2_norm$object_neg.impute,colby = p2_norm$temp_norm_pca_col_by,center = T,scale = T,removeVar = .1,interactive = F)

        # save plot

        ggsave(
          filename = file,
          plot = p,
          width = para_d$fig5_width,
          height = para_d$fig5_height,
          device = para_d$fig5_format
        )
      }
    )
    ###> fig6 =====
    output$fig6_download = downloadHandler(
      filename = function() {
        paste0("03.rsd_plot-neg.", download_para()$fig6_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # draw condition
        p = p2_norm$neg.plt

        # save plot
        ggsave(
          filename = file,
          plot = p,
          width = para_d$fig6_width,
          height = para_d$fig6_height,
          device = para_d$fig6_format
        )
      }
    )


  })
}

