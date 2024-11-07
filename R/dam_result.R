#' DAM data analysis
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyjqui jqui_resizable
#' @importFrom DT dataTableOutput
#' @importFrom shinyWidgets materialSwitch
#' @noRd


dam_res_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    useShinyjs(),
    title = 'DAM analysis',
    icon = icon("stairs"),
    sidebarLayout(
      div(id = ns('Sidebar'),
          #### sidebarPanel ======
          sidebarPanel(width = 2,
                       tags$h3("Compare group",style = 'color: #008080'),
                       hr_main(),
                       HTML(
                         "<span style='color: #7a8788; font-size: 12px; font-style: italic;'> Automatically group based on 'groups' in sample information.
                       </span>"
                       ),
                       radioButtons(
                         inputId = ns("dam_method_picksample"),
                         label = "Methods",
                         choices = c("Based on group information","Manual selection"),
                         selected = "Based on group information"
                       ),
                       hr_head(),
                       selectInput_div(
                         inputId = ns('dam_left_name_g'),
                         label = "left name (group)",
                         choices = c("case","control"),
                         selected = "case",
                         multiple = F,
                         title = "left vs right. left name"
                       ),
                       selectInput_div(
                         inputId = ns('dam_right_name_g'),
                         label = "right name (group)",
                         choices = c("case","control"),
                         selected = "control",
                         multiple = F,
                         title = "left vs right. right name"
                       ),
                       hr_head(),
                       textInput(
                         inputId = ns('dam_left_name_m'),
                         label = 'left name (manual)',
                         value = "left"
                       ),
                       textInput(
                         inputId = ns('dam_right_name_m'),
                         label = 'right name (manual)',
                         value = "right"
                       ),
                       selectInput_div(
                         inputId = ns('dam_left_sid'),
                         label = "left sample id (manual)",
                         choices = c("S_0001","S_0002"),
                         selected = c("S_0001","S_0002"),
                         multiple = T,
                         title = "left vs right. left sample id"
                       ),
                       selectInput_div(
                         inputId = ns('dam_right_sid'),
                         label = "right sample id (manual)",
                         choices = c("S_0003","S_0004"),
                         selected = c("S_0003","S_0004"),
                         multiple = T,
                         title = "left vs right. right sample id"
                       ),
                       tags$h3("Parameters for DAM",style = 'color: #008080'),
                       hr_main(),
                       radioButtons(
                         inputId = ns("dam_data_form"),
                         label = "Logarithmic transformation",
                         choices = c("TRUE","FALSE"),
                         selected = "FALSE"
                       ),
                       textInput(
                         inputId = ns('dam_pvalue'),
                         label = "pvalue cutoff",
                         value = 0.05
                       ),
                       textInput(
                         inputId = ns('dam_qvalue'),
                         label = "qvalue cutoff",
                         value = 1
                       ),
                       textInput(
                         inputId = ns('dam_VIP'),
                         label = "VIP cutoff",
                         value = 1
                       ),
                       textInput(
                         inputId = ns('dam_log2fc'),
                         label = "log2fc cutoff",
                         value = 0
                       ),
                       selectInput(
                         inputId = ns("dam_method1"),
                         label  = "Univariate test method",
                         choices = c("t-test","wilcox-test"),
                         selected = "t-test"
                       ),
                       radioButtons(
                         inputId = ns('dam_paired'),
                         label = "Paired",
                         choices = c("TRUE","FALSE"),
                         selected = "FALSE"
                       ),
                       selectInput(
                         inputId = ns("dam_method2"),
                         label  = "Multivariate test method",
                         choices = c("opls-da","pls-da"),
                         selected = "opls-da"
                       ),
                       selectInput(
                         inputId = ns("dam_method3"),
                         label  = "p adjust method",
                         choices = c("holm","hochberg","hommel","bonferroni","BH","BY","fdr","none"),
                         selected = "BH"
                       )
          )# sidebarPanel
      ),# div
      #### MainPanel ========
      mainPanel(
        fluidPage(
          actionButton(ns('toggleSidebar'),"Toggle sidebar"),
          actionButton(inputId = ns('dam_initialize'),label = "Wake up",icon = icon("bolt")),
          actionButton(inputId = ns('dam_check'),label = "Check your compaire group",icon = icon("check")),
          br(),
          HTML(
            "<span style='color: #7a8788; font-size: 12px; font-style: italic;'> If you need to start the analysis from the 'Resuming analysis from the unfinished steps (options)', click 'wake up' first to activate the mass_dataset..</span>"
          ),
          hr_bar(),
          materialSwitch(inputId = ns("dam_plt_format"),label = "Interactive plot", status = "primary"),
          tabsetPanel(
            tabPanel(
              title = 'File check',height = '500px',width = "100%",
              icon = icon('check'),
              hr_head(),
              uiOutput(ns("dam_file_check"),fill = T)
            ),
            tabPanel(
              title = 'DAM result',height = '500px',width = "100%",
              icon = icon('table'),
              hr_main(),
              actionButton(inputId = ns("dam_start"),label = "Start DAM analysis",icon = icon("play")),
              tags$h3("Upset plot.",style = 'color: #008080'),
              hr_main(),
              jqui_resizable(
                plotOutput(ns("dam_pos_upset"))
              ),
              tags$h3("DAMs (Differentially Accumulated Metabolites) results.",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("dam_tbl_DAMs")),
              tags$h3("DAMs analysis for all features",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("dam_tbl_all")),
            ),
            tabPanel(
              title = 'Multivariate statistical analysis',height = '500px',width = "100%",
              icon = icon('image'),
              hr_main(),
              actionButton(inputId = ns("dam_pca_show"),label = "Show plot",icon = icon("chart")),
              tags$h3("PCA",style = 'color: #008080'),
              hr_main(),
              fluidRow(
                column(width = 2,
                       tags$h4("parameters",style = "color: #008080"),
                       hr_head(),
                       selectInput(
                         inputId = ns("dam_pca_col_by"),
                         label = "Color by",
                         choices = "group",selected = 'group',multiple = F
                       ),
                       radioButtons(
                         inputId = ns("dam_showloading"),
                         label = "show loading",choices = c("TRUE","FALSE"),
                         selected = "TRUE"
                       )
                ),
                column(width = 5,
                       tags$h4("Score",style = 'color: #008080'),
                       hr_head(),
                       jqui_resizable(
                         uiOutput(ns("dam_pca_score"),fill = T)
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
                       downloadButton(ns("fig1_download"),"Download")
                ),
                column(width = 4,
                       tags$h4("Loading",style = 'color: #008080'),
                       hr_head(),
                       jqui_resizable(
                         uiOutput(ns("dam_pca_loading"),fill = T)
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
              ),
              tags$h3("PLS-DA (OPLS-DA)",style = 'color: #008080'),
              hr_main(),
              fluidRow(
                column(width = 4,
                       tags$h4("Summary",style = 'color: #008080'),
                       hr_head(),
                       verbatimTextOutput(ns("dam_opls_summary"))
                ),
                column(width = 8,
                       tags$h4("Plot auto",style = 'color: #008080'),
                       hr_head(),
                       jqui_resizable(
                         plotOutput(ns("opls_da_plt_auto"))
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
                       downloadButton(ns("fig3_download"),"Download")
                       )
                )
            ),
            tabPanel(
              title = 'Univariate statistical analysis',height = '500px',width = "100%",
              icon = icon('image'),
              tags$h3("Plot",style = 'color: #008080'),
              hr_main(),
              actionButton(inputId = ns("dam_voc_show"),label = "Show plot",icon = icon("chart")),
              fluidRow(
                column(width = 8,
                       tags$h4("Volcano plot",style = 'color: #008080'),
                       hr_head(),
                       jqui_resizable(
                         uiOutput(ns("dam_voc_plot"),fill = T)
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
                       downloadButton(ns("fig4_download"),"Download")
                ),
                column(width = 4,
                       tags$h4("boxplot",style = 'color: #008080'),
                       hr_head(),
                       jqui_resizable(
                         plotOutput(ns("dam_barplot"))
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
                       downloadButton(ns("fig5_download"),"Download")
                )
              ),
              tags$h3("Details",style = 'color: #008080'),
              hr_main(),
              fluidRow(
                column(width = 4,
                       tags$h4("Compound annotation",style = 'color: #008080'),
                       hr_head(),
                       verbatimTextOutput(ns("dam_anno_data"))
                ),
                column(width = 4,
                       tags$h4("MS2 spectra",style = 'color: #008080'),
                       hr_head(),
                       dataTableOutput(ns("dam_ms2_spectra"))
                ),
                column(width = 4,
                       tags$h4("Structure (pubchem)",style = 'color: #008080'),
                       hr_head(),
                       htmlOutput(ns("DAM_structure"))
                )
              )
            )
          )
        )
      )
    )
  )
}


#' DAM analysis
#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs toggle runjs
#' @importFrom ropls plot
#' @importFrom dplyr select pull all_of left_join mutate across where case_when arrange
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom massdataset activate_mass_dataset extract_annotation_table extract_sample_info extract_variable_info extract_expression_data export_mass_dataset4metdna export_ms2_data
#' @importFrom plotly plotlyOutput renderPlotly ggplotly event_register event_data
#' @importFrom writexl write_xlsx
#' @importFrom PCAtools biplot pca plotloadings
#' @importFrom tidyr pivot_longer
#' @importFrom ggstatsplot ggbetweenstats
#' @importFrom ComplexHeatmap make_comb_mat set_size comb_size UpSet comb_degree HeatmapAnnotation anno_barplot rowAnnotation anno_text max_text_width set_name draw
#' @import MDAtoolkits
#' @importFrom circlize colorRamp2
#' @importFrom stringr str_split
#' @importFrom grid gpar
#' @param id module of server
#' @param volumes shinyFiles volumes
#' @param prj_init use project init variables.
#' @param data_clean_rv reactivevalues p2 dataclean
#' @noRd


dam_res_server <- function(id,volumes,prj_init,data_clean_rv,data_download) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    observeEvent(input$toggleSidebar, {
      shinyjs::toggle(id = "Sidebar")
    })

    p3_DAM <- reactiveValues(data = NULL)

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
        fig5_format = as.character(input$fig5_format)
      )
    })

    observeEvent(
      input$dam_initialize,
      {
        if(!is.null(prj_init$object_positive.init) & prj_init$steps == "DAM and rest"){
          p3_DAM$object_merge = prj_init$object_positive.init;
        } else if(!is.null(data_clean_rv$object_merge)){
          p3_DAM$object_merge = data_clean_rv$object_merge
        } else {
          return()
        }
        p3_DAM$Group_method = input$dam_method_picksample %>% as.character()

        p3_DAM$sample_info =  p3_DAM$object_merge %>% extract_sample_info()

        p3_DAM$anno =  p3_DAM$object_merge %>% extract_annotation_table()

        p3_DAM$variable = p3_DAM$object_merge %>% extract_variable_info() %>%
          dplyr::select(1:3)

        p3_DAM$expmat =  p3_DAM$object_merge %>% extract_expression_data() %>%
          rownames_to_column("variable_id")

        p3_DAM$s_group = p3_DAM$sample_info %>% pull(group) %>% unique()

        p3_DAM$s_sample_id =  p3_DAM$sample_info %>% pull(sample_id)

        if(p3_DAM$Group_method == "Based on group information") {
          observe({
            updateSelectInput(session, "dam_left_name_g",choices = p3_DAM$s_group,selected = (p3_DAM$s_group)[1])
            updateSelectInput(session, "dam_right_name_g",choices = p3_DAM$s_group,selected = (p3_DAM$s_group)[2])
          })
        } else if(p3_DAM$Group_method == "Manual selection") {
          observe({
            updateSelectInput(session, "dam_left_sid",choices = p3_DAM$s_sample_id,selected = (p3_DAM$s_sample_id)[1])
            updateSelectInput(session, "dam_right_sid",choices = p3_DAM$s_sample_id,selected = (p3_DAM$s_sample_id)[2])
          })
        }

        p3_DAM$file_check_status = NULL
      }

    )

    observeEvent(
      input$dam_check,
      {
        if(is.null(p3_DAM$sample_info)) {return()}
        if(is.null(p3_DAM$anno)) {return()}
        if(is.null(p3_DAM$variable)) {return()}
        if(is.null(p3_DAM$expmat)) {return()}

        p3_DAM$dam_data_form = input$dam_data_form %>% as.logical()
        p3_DAM$dam_pvalue = input$dam_pvalue %>% as.numeric()
        p3_DAM$dam_qvalue = input$dam_qvalue %>% as.numeric()
        p3_DAM$dam_VIP = input$dam_VIP %>% as.numeric()
        p3_DAM$dam_log2fc = input$dam_log2fc %>% as.numeric()
        p3_DAM$dam_method1 = input$dam_method1 %>% as.character()
        p3_DAM$dam_method2 = input$dam_method2 %>% as.character()
        p3_DAM$dam_method3 = input$dam_method3 %>% as.character()
        p3_DAM$dam_paired = input$dam_paired %>% as.logical()

        if(p3_DAM$Group_method == "Based on group information") {

          p3_DAM$dam_left_name = input$dam_left_name_g %>% as.character()

          p3_DAM$dam_right_name = input$dam_right_name_g %>% as.character()

          p3_DAM$dam_left_sid = p3_DAM$sample_info %>%
            filter(group == p3_DAM$dam_left_name) %>% pull(sample_id)

          p3_DAM$dam_right_sid = p3_DAM$sample_info %>%
            filter(group == p3_DAM$dam_right_name) %>% pull(sample_id)

        } else if(p3_DAM$Group_method == "Manual selection") {

          p3_DAM$dam_left_name = input$dam_left_name_m %>% as.character()

          p3_DAM$dam_right_name = input$dam_right_name_m %>% as.character()

          p3_DAM$dam_left_sid = input$dam_left_sid %>% as.character()

          p3_DAM$dam_right_sid = input$dam_right_sid %>% as.character()
        }



        output$dam_file_check <- renderUI({
          isolate(HTML(paste0(
            '<h4>DAMs Compare group setting</h4><hr>',
            '<p><b>Group A: </b><span style="color: blue;">', as.character(p3_DAM$dam_left_name), '</span></p>',
            '<p><b>Group B: </b><span style="color: blue;">', as.character(p3_DAM$dam_right_name), '</span></p>',
            '<p><b>Sample id (A): </b><span style="color: blue;">', paste(as.character(p3_DAM$dam_left_sid), collapse = " | "), '</span></p>',
            '<p><b>Sample id (B): </b><span style="color: blue;">', paste(as.character(p3_DAM$dam_right_sid), collapse = " | "), '</span></p>',
            '<p><b>Expression data format: </b><span style="color: blue;">', as.character(p3_DAM$dam_data_form), '</span></p>',
            '<h4>DAMs Methods Setting</h4><hr>',
            '<p><b>Univariate: </b><span style="color: blue;">', as.character(p3_DAM$dam_method1), '</span></p>',
            '<p><b>Multivariate: </b><span style="color: blue;">', as.character(p3_DAM$dam_method2), '</span></p>',
            '<p><b>P adjusted: </b><span style="color: blue;">', as.character(p3_DAM$dam_method3), '</span></p>',
            '<p><b>Paired: </b><span style="color: blue;">', as.character(p3_DAM$dam_paired), '</span></p>',
            '<h4>DAMs Cut off Setting</h4><hr>',
            '<p><b>Pvalue: </b><span style="color: blue;">', as.character(p3_DAM$dam_pvalue), '</span></p>',
            '<p><b>Qvalue: </b><span style="color: blue;">', as.character(p3_DAM$dam_qvalue), '</span></p>',
            '<p><b>Log2fc: </b><span style="color: blue;">', as.character(p3_DAM$dam_log2fc), '</span></p>',
            '<p><b>VIP: </b><span style="color: blue;">', as.character(p3_DAM$dam_VIP), '</span></p>'
          )))
        })

        p3_DAM$expmat_paired_mat = p3_DAM$expmat %>% select(variable_id,all_of(c(p3_DAM$dam_left_sid,p3_DAM$dam_right_sid)))
        p3_DAM$sample_info_sim = p3_DAM$sample_info %>% select(sample_id,group)

        p3_DAM$expmat_long <-
          p3_DAM$expmat_paired_mat %>%
          pivot_longer(!variable_id,names_to = "sample_id",values_to = "value") %>%
          left_join(p3_DAM$sample_info_sim)

        p3_DAM$file_check_status <- "finish"
      })

    observeEvent(
      input$dam_start,
      {
        if(is.null(p3_DAM$file_check_status)) {return()}
        if(isTRUE(p3_DAM$dam_data_form)) {
          p3_DAM$expmat = p3_DAM$expmat %>%
            mutate(across(where(is.numeric), ~log2(. + 1)))
        }
        p3_DAM$dam_file_path = paste0(prj_init$wd,"/Result/DAM/",p3_DAM$dam_left_name,"_vs_",p3_DAM$dam_right_name,"/")
        dir.create(path = p3_DAM$dam_file_path,showWarnings = F,recursive = T)

        ## progress
        tags = c("dam analysis","upsetplot")
        dam_steps_tag = c(paste0("Running ",tags," in progress..."),"Finish!")

        dam_steps = length(dam_steps_tag)
        withProgress(message = 'DAM analysis', value = 0,
                     expr = {
                       for (i in 1:(dam_steps)) {
                         incProgress(1/dam_steps,detail = dam_steps_tag[i])
                         if(i == 1) {
                           ##> DAM analysis
                           res_obj = MDAtoolkits::DAM_analysis(
                             x= p3_DAM$expmat,
                             left_index = p3_DAM$dam_left_sid,
                             right_index = p3_DAM$dam_right_sid,
                             left = p3_DAM$dam_left_name,
                             right = p3_DAM$dam_right_name,
                             method = p3_DAM$dam_method1,
                             method2 = p3_DAM$dam_method2,
                             method3 = p3_DAM$dam_method3,
                             paired = p3_DAM$dam_paired
                           )
                           ##> split
                           p3_DAM$res_tbl_all = res_obj$DAM_tbl
                           p3_DAM$oplsda = res_obj$opls_out
                           ##> export
                           writexl::write_xlsx( p3_DAM$res_tbl_all,paste0(p3_DAM$dam_file_path,"DAM_analysis.xlsx"))

                           png(filename = paste0(p3_DAM$dam_file_path,p3_DAM$dam_method2,".png"),width = 3000,height = 2500,res = 300)
                           ropls::plot(res_obj$opls_out)
                           dev.off()
                         } else if (i == 2) {
                           ##> upset plot
                           FvsH_upset = p3_DAM$res_tbl_all %>%
                             mutate(
                               log2fc = case_when(
                                 abs(log2fc) >= 1 ~ 1,
                                 TRUE ~ 0
                               ),
                               pvalue = case_when(
                                 pvalue <= 0.05 ~ 1,
                                 TRUE ~ 0
                               ),
                               FDR = case_when(
                                 FDR <= 0.05 ~ 1,
                                 TRUE ~ 0
                               ),
                               VIP = case_when(
                                 VIP >= 1 ~ 1,
                                 TRUE ~ 0
                               )
                             ) %>%
                             column_to_rownames("CompoundID") %>%
                             as.matrix()

                           m = make_comb_mat(FvsH_upset)
                           ss = set_size(m)
                           cs = comb_size(m)

                           p3_DAM$ht = UpSet(
                             m,
                             set_order = order(ss),
                             comb_order = order(comb_degree(m), -cs),
                             top_annotation = HeatmapAnnotation(
                               "DM Intersections" = anno_barplot(
                                 cs,
                                 ylim = c(0, max(cs)*1.1),
                                 border = FALSE,
                                 gp = gpar(fill = "salmon",alpha = 0.6),
                                 height = unit(4, "cm")
                               ),
                               annotation_name_side = "left",
                               annotation_name_rot = 90
                             ),
                             left_annotation = rowAnnotation(
                               "DM number" = anno_barplot(-ss,
                                                          baseline = 0,
                                                          border = FALSE,
                                                          gp = gpar(fill = "cyan",alpha = 0.6),
                                                          width = unit(4, "cm")
                               ),
                               set_name = anno_text(set_name(m),
                                                    location = 0.5,
                                                    just = "center",
                                                    width = max_text_width(set_name(m)) + unit(4, "mm"))
                             ),
                             right_annotation = NULL,
                             show_row_names = FALSE
                           )


                           output$dam_pos_upset = renderPlot({
                             if(is.null(p3_DAM$ht)){return()}
                             print(draw(p3_DAM$ht))
                           })

                         } else {
                           Sys.sleep(1)
                         }
                       }
                     }
        )

        p3_DAM$res_tbl_filtered =
          p3_DAM$res_tbl_all %>%
          filter(pvalue <= p3_DAM$dam_pvalue) %>%
          filter(FDR <= p3_DAM$dam_qvalue ) %>%
          filter(VIP >= p3_DAM$dam_VIP) %>%
          filter(abs(log2fc) >= p3_DAM$dam_log2fc)

        p3_DAM$dam_pvalue = input$dam_pvalue %>% as.numeric()
        p3_DAM$dam_qvalue = input$dam_qvalue %>% as.numeric()
        p3_DAM$dam_VIP = input$dam_VIP %>% as.numeric()
        p3_DAM$dam_log2fc = input$dam_log2fc %>% as.numeric()


        output$dam_tbl_DAMs = renderDataTable_formated(
          actions = p3_DAM$dam_start,
          condition1 = p3_DAM$res_tbl_filtered,
          filename.a = "3.6.9.DAMs_filtered.csv",
          tbl = p3_DAM$res_tbl_filtered
        )

        output$dam_tbl_all = renderDataTable_formated(
          actions = p3_DAM$dam_start,
          condition1 = p3_DAM$res_tbl_all,
          filename.a = "3.6.9.DAMs_all.csv",
          tbl = p3_DAM$res_tbl_all
        )

        observe({
          updateSelectInput(session, "dam_pca_col_by",choices = colnames(p3_DAM$sample_info),selected = "group")
        })

      }
    )

    observeEvent(input$dam_pca_show,{

      ##> opls-da
      output$dam_opls_summary <- renderPrint({
        if(is.null(input$dam_pca_show)) {return()}
        if(is.null(p3_DAM$oplsda)) {return()}
        p3_DAM$oplsda
      })


      output$opls_da_plt_auto = renderPlot({
        if(is.null(input$dam_pca_show)) {return()}
        if(is.null(p3_DAM$oplsda)) {return()}
        ropls::plot(p3_DAM$oplsda)
      })
      ##> PCA
      dam_showloading = input$dam_showloading %>% as.logical()
      dam_pca_col_by = input$dam_pca_col_by %>% as.character()
      temp_si1 =  p3_DAM$sample_info %>%
        filter(group == p3_DAM$dam_left_name)
      temp_si2 = p3_DAM$sample_info %>%
        filter(group == p3_DAM$dam_right_name)
      temp_sample_info = rbind(temp_si1,temp_si2) %>% column_to_rownames("sample_id")

      temp_expmat <- p3_DAM$expmat %>% column_to_rownames('variable_id') %>%
        filter(rowSums(.) > 0) %>% select(rownames(temp_sample_info))

      temp_pca = PCAtools::pca(
        mat = temp_expmat,
        metadata = temp_sample_info,
        center = T,scale = T,
        removeVar = .1
      )


      ##> PCA biplot
      temp_biplot = PCAtools::biplot(
        pcaobj = temp_pca,
        colby = dam_pca_col_by,legendPosition = "top",
        showLoadings = dam_showloading
      )
      p3_DAM$temp_biplot = temp_biplot


      output$dam_pca_score <- renderUI({
        plot_type <- input$dam_plt_format
        if (plot_type) {
          plotlyOutput(outputId = ns("dam_plotly_pca"))
        } else {
          plotOutput(outputId = ns("dam_plot_pca"))
        }
      })

      output$dam_plot_pca <- renderPlot({
        if(is.null(temp_biplot)){return()}
        temp_biplot
      })


      output$dam_plotly_pca <- renderPlotly({

        if(is.null(temp_biplot)){return()}
        temp_biplot %>% plotly::ggplotly()

      })
      ##> PCA loading plot


      temp_loading_plot = PCAtools::plotloadings(temp_pca)
      p3_DAM$temp_loading_plot = temp_loading_plot

      output$dam_pca_loading <- renderUI({
        plot_type <- input$dam_plt_format
        if (plot_type) {
          plotlyOutput(outputId = ns("dam_plotly_loading"))
        } else {
          plotOutput(outputId = ns("dam_plot_loading"))
        }
      })

      output$dam_plot_loading <- renderPlot({
        if(is.null(temp_loading_plot)){return()}
        temp_loading_plot
      })

      output$dam_plotly_loading <- renderPlotly({
        if(is.null(temp_loading_plot)){return()}
        temp_loading_plot %>% plotly::ggplotly()
      })
    })


# voc plot ----------------------------------------------------------------

    observeEvent(input$dam_voc_show,{
      output$dam_voc_plot <- renderUI({
        plot_type <- input$dam_plt_format
        if (plot_type) {
          plotlyOutput(outputId = ns("voc_plotly"))
        } else {
          plotOutput(outputId = ns("voc_plot"))
        }
      })

      output$voc_plot <- renderPlot({
        if(is.null(p3_DAM$res_tbl_all)){return()}
        MDAtoolkits::DAM_volcano(
          x = p3_DAM$res_tbl_all,
          title = paste0(p3_DAM$dam_left_name,"_vs_",p3_DAM$dam_right_name),
          pval_cut = p3_DAM$dam_pvalue,
          log2fc_cut = p3_DAM$dam_log2fc,
          qval_cut = p3_DAM$dam_qvalue,
          VIP_cut = p3_DAM$dam_VIP
        )
      })

      output$voc_plotly <- renderPlotly({
        if(is.null(p3_DAM$res_tbl_all)){return()}
        p_voc = MDAtoolkits::DAM_volcano_plotly(
          x = p3_DAM$res_tbl_all,
          title = paste0(p3_DAM$dam_left_name,"_vs_",p3_DAM$dam_right_name),
          pval_cut = p3_DAM$dam_pvalue,
          log2fc_cut = p3_DAM$dam_log2fc,
          qval_cut = p3_DAM$dam_qvalue,
          VIP_cut = p3_DAM$dam_VIP
        )
        event_register(p_voc, 'plotly_click')
        p_voc
      })
      d <- event_data("plotly_click")
      d = (d$customdata %>% str_split( pattern = "\\|",n = Inf,simplify = T))[1]
      output$dam_barplot <- renderPlot({
        d <- event_data("plotly_click")
        d = (d$customdata %>% str_split( pattern = "\\|",n = Inf,simplify = T))[1]
        if (is.null(d)) "Hover on a point!" else d
        tbl = p3_DAM$expmat_long %>% filter(variable_id == d)
        if(p3_DAM$dam_method1 == "t-test") {temp_type = "parametric"} else {temp_type = "nonparametric"}
        p3_DAM$boxplot = ggbetweenstats(tbl,group,value,p.adjust.method = p3_DAM$dam_method3,temp_type = "parametric",paired = p3_DAM$dam_paired,
                       title = d)
        p3_DAM$boxplot
      })

      output$dam_anno_data <- renderPrint({
        d <- event_data("plotly_click")
        d = (d$customdata %>% str_split( pattern = "\\|",n = Inf,simplify = T))[1]
        d = p3_DAM$variable %>% left_join(p3_DAM$anno,by = 'variable_id')  %>%
          select(variable_id,mz,rt,Formula,Compound.name,Adduct,KEGG.ID,Lab.ID,Total.score,Level,superclass,class,subclass,parent_levels,description) %>%
          filter(variable_id == d) %>% t()

        if (is.null(d)) "Hover on a point!" else d
      })

      temp_x =  p3_DAM$object_merge@ms2_data
      names(temp_x) = c("neg","pos")
      temp_x_pos = temp_x$pos
      temp_x_neg = temp_x$neg

      output$dam_ms2_spectra <- renderDataTable({
        d <- event_data("plotly_click")
        d = (d$customdata %>% str_split( pattern = "\\|",n = Inf,simplify = T))[1]
        if(d %in% temp_x_pos@variable_id){
          temp_ms2 = temp_x_pos %>% filter(variable_id == d)
          temp_ms2_tbl = temp_ms2@ms2_spectra %>%
            as.data.frame() %>%
            setNames(c("mz","intensity")) %>%
            mutate(intensity = (intensity/max(intensity))*100) %>%
            arrange(desc(intensity))
          temp_ms2_tbl
        } else if(d %in% temp_x_neg@variable_id) {
          temp_ms2 = temp_x_neg %>% filter(variable_id == d)
          temp_ms2_tbl = temp_ms2@ms2_spectra %>%
            as.data.frame() %>%
            setNames(c("mz","intensity")) %>%
            mutate(intensity = (intensity/max(intensity))*100) %>%
            arrange(desc(intensity))
          temp_ms2_tbl
        } else {
          temp_ms2_tbl = data.frame(
            mz = "",intensity = ""
          )
        }
      })

      output$DAM_structure = renderUI({
        tryCatch({
          d <- event_data("plotly_click")
          d = (d$customdata %>% str_split(pattern = "\\|", n = Inf, simplify = TRUE))[1]
          tempx = p3_DAM$anno %>% filter(variable_id == d) %>% pull(Compound.name)
          d = tempx[1]

          # 使用tryCatch来捕获mda_get_cid_fast中的错误
          tempx_inch <- tryCatch({
            MDAtoolkits::mda_get_cid_fast(query = d, probe = 'name', core_num = 1) %>% pull(cid)
          }, error = function(e) {
            NA  # 如果发生错误，返回NA
          })

          # 检查tempx_inch是否为NA
          if (!is.na(tempx_inch)) {
            url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", tempx_inch, "/PNG")
            isolate(HTML(paste0(
              '<img src="', url, '" alt="Image description" style="width:100%; height:auto;">'
            )))
          } else {
            isolate(HTML('<font color="red">"No database were selected!"</font>'))
          }
        }, error = function(e) {
          isolate(HTML('<font color="red">"An error occurred!"</font>'))
        })
      })

    })
    # download ----------------------------------------------------------------
    ###> fig1 =====
    output$fig1_download = downloadHandler(
      filename = function() {
        paste0(p3_DAM$dam_left_name,"_vs_",p3_DAM$dam_right_name,"PCA_score.", download_para()$fig1_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # draw condition
        p = p3_DAM$temp_biplot
        # save plot
        ggsave(
          filename = file,
          plot = p,
          width = para_d$fig1_width,
          height = para_d$fig1_height,
          device = para_d$fig1_format,limitsize = FALSE
        )
      }
    )
    ###> fig2 ====
    output$fig2_download = downloadHandler(
      filename = function() {
        paste0(p3_DAM$dam_left_name,"_vs_",p3_DAM$dam_right_name,"PCA_loading.", download_para()$fig2_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # draw condition
        p = p3_DAM$temp_loading_plot

        # save plot

        ggsave(
          filename = file,
          plot = p,
          width = para_d$fig2_width,
          height = para_d$fig2_height,
          device = para_d$fig2_format,limitsize = FALSE
        )
      }
    )
    ###> fig3 ====
    output$fig3_download = downloadHandler(
      filename = function() {
        paste0(p3_DAM$dam_left_name,"_vs_",p3_DAM$dam_right_name,"OPLS-DA.", download_para()$fig3_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # save plot
        if(para_d$fig2_format == "pdf"){
          pdf(file = file,width = para_d$fi3_width,para_d$fig3_height)
          ropls::plot(p3_DAM$oplsda)
          dev.off()
        } else if(para_d$fig2_format == "tiff") {
          tiff(filename = file,width = para_d$fi3_width,height = para_d$fig3_height,res = 300,units = 'in')
          ropls::plot(p3_DAM$oplsda)
          dev.off()
        } else if(para_d$fig2_format == "png") {
          png(filename = file,width = para_d$fi3_width,height = para_d$fig3_height,res = 300,units = 'in',)
          ropls::plot(p3_DAM$oplsda)
          dev.off()
        } else if(para_d$fig2_format == "jpg") {
          grDevices::jpeg(filename = file,width = para_d$fi3_width,height = para_d$fig3_height,quality = 300,units = 'in')
          ropls::plot(p3_DAM$oplsda)
          dev.off()
        }

      }
    )
    ###> fig4 ====
    output$fig4_download = downloadHandler(
      filename = function() {
        paste0(p3_DAM$dam_left_name,"_vs_",p3_DAM$dam_right_name,"Volcano.", download_para()$fig4_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # save plot
        p = MDAtoolkits::DAM_volcano(
          x = p3_DAM$res_tbl_all,
          title = paste0(p3_DAM$dam_left_name,"_vs_",p3_DAM$dam_right_name),
          pval_cut = p3_DAM$dam_pvalue,
          log2fc_cut = p3_DAM$dam_log2fc,
          qval_cut = p3_DAM$dam_qvalue,
          VIP_cut = p3_DAM$dam_VIP
        )
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
        paste0("variable_boxplot.", download_para()$fig5_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # save plot
        p = p3_DAM$boxplot
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

  })
}

