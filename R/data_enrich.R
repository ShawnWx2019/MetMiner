#' Enrichment
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyjqui jqui_resizable
#' @importFrom shinyWidgets materialSwitch
#' @noRd


data_enrich_ui <- function(id) {
  ns <- NS(id)
  ### Part1.4.2 Remove noisy features -----------------------------------------

  tabPanel(
    useShinyjs(),
    title = 'Enrichment',
    icon = icon("slack"),
    sidebarLayout(
      div(id = ns('Sidebar'),
          #### sidebarPanel ======
          sidebarPanel(width = 2,
                       tags$h3("Parameters",style = 'color: #008080'),
                       hr_main(),
                       radioButtons(
                         inputId = ns('enrich_method'),
                         label = 'Method',
                         choices = c("from new upload annotation file","from mass_dataset"),
                         selected = "from mass_dataset"
                       ),
                       hr_head(),
                       fileInput(
                         inputId = ns('enrich_annotation_file'),
                         label = "upload annotation file (option)",
                         accept = ".csv"
                       ),
                       hr_main(),
                       radioButtons(
                         inputId = ns('enrich_db'),
                         label = 'database',
                         choices = c("ClassyFire","KEGG")
                       ),
                       hr_head(),
                       textInput(
                         inputId = ns('enrich_org'),
                         label = 'Species',
                         value = 'ath'
                       ),
                       hr_main(),
                       HTML(
                         "<span style='color: #7a8788; font-size: 12px; font-style: italic;'>Select the corresponding species code based on </span>
                        <a href='https://www.genome.jp/kegg/catalog/org_list.html' style='color: red; font-size: 12px; font-style: italic;'>https://www.genome.jp/kegg/catalog/org_list.html</a>
                        <span style='color: #7a8788; font-size: 12px; font-style: italic;'> </span>"
                       ),
                       radioButtons(
                         inputId = ns('enrich_plot_type'),
                         label = 'plot type',choices = c("bar","bubble"),
                         selected = 'bar'
                       ),
                       hr_head(),
                       sliderInput(
                         inputId = ns('enrich_pval'),
                         label = 'Pvalue',min = 0,max = 1,step = 0.01,
                         value = 0.05
                       ),
                       sliderInput(
                         inputId = ns('enrich_qval'),
                         label = 'qvalue',min = 0,max = 1,step = 0.01,
                         value = 1
                       ),
                       textInput(
                         inputId = ns('enrich_num'),
                         label = 'Term number',
                         value = 20
                       )
          )# sidebarPanel
      ),# div
      #### MainPanel ========
      mainPanel(
        fluidPage(
          actionButton(ns('toggleSidebar'),"Toggle sidebar"),
          actionButton(inputId = ns('enrich_check'),label = "Start",icon = icon("check")),
          actionButton(inputId = ns('enrich_start'),label = "Start",icon = icon("play")),
          hr_head(),
          materialSwitch(inputId = ns("enrich_plt_format"),label = "Interactive plot", status = "primary"),
          tabsetPanel(
            tabPanel(
              title = 'Database',height = '500px',width = "100%",
              icon = icon('database'),
              tags$h3("Annotation file check",style = 'color: #008080'),
              hr_main(),
              htmlOutput(outputId = ns("class_file_check")),
              tags$h3("Back ground Term to name",style = 'color: #008080'),
              hr_main(),
              dataTableOutput(outputId = ns("class_tbl_t2n")),
              tags$h3("Back ground Term to varaiable id",style = 'color: #008080'),
              hr_main(),
              dataTableOutput(outputId = ns("class_tbl_t2g")),
            ),
            tabPanel(
              title = 'Table',height = '500px',width = "100%",
              icon = icon('table'),
              tags$h3("variable set",style = 'color: #008080'),
              hr_main(),
              HTML(
                "
              Enter a list of metabolite IDs, which can come from differential metabolites or from metabolites in a WGCNA module. Each ID should be on a new line or separated by commas.
              "
              ),
              hr_head(),
              textAreaInput(inputId = ns("enrich_glist"),
                            label = "Input id",
                            value = "",
                            resize = "both"),

              tags$h3("Enrichment table",style = 'color: #008080'),
              hr_main(),
              DT::dataTableOutput(ns("enrich_table")),
            ),
            tabPanel(
              title = 'Plot',height = '500px',width = "100%",
              icon = icon('image'),
              tags$h3("enrichment plot",style = 'color: #008080'),
              hr_main(),
              jqui_resizable(
                uiOutput(ns("enrich_plot"),fill = T)
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
#' @importFrom shinyjs useShinyjs
#' @importFrom dplyr select all_of left_join pull filter rename group_by slice_head ungroup distinct
#' @importFrom tibble rownames_to_column
#' @importFrom massdataset activate_mass_dataset extract_expression_data extract_sample_info extract_annotation_table
#' @importFrom plotly renderPlotly plotlyOutput ggplotly
#' @importFrom clusterProfiler enricher
#' @importFrom purrr map_dfr
#' @importFrom stringr str_detect str_pad str_split
#' @importFrom tidyr pivot_longer separate_longer_delim
#' @importFrom enrichplot dotplot
#' @import MDAtoolkits
#' @param id module of server
#' @param volumes shinyFiles volumes
#' @param prj_init use project init variables.
#' @param data_clean_rv reactivevalues p2 dataclean
#' @noRd


data_enrich_server <- function(id,volumes,prj_init,data_clean_rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns


    #> variable table
    p3_enrich <- reactiveValues(data = NULL)

    #> parameters
    ##> download parameters ================
    download_para = reactive({
      list(
        ##> fig1
        fig1_width = as.numeric(input$fig1_width),
        fig1_height = as.numeric(input$fig1_height),
        fig1_format = as.character(input$fig1_format)
      )
    })

    temp_enrich_anno <- reactive({
      file1 <- input$enrich_annotation_file
      if(is.null(file1)){return()}
      read.csv(file = file1$datapath,
               sep=",",header = T,stringsAsFactors = F)
    })



    #> File check event
    observeEvent(
      input$enrich_check,
      {
        if (!is.null(temp_enrich_anno())) {
          temp_anno_file = temp_enrich_anno()
        } else if(!is.null(prj_init$object_positive.init) & prj_init$steps == "DAM and rest"){
          p3_enrich$object_merge= prj_init$object_positive.init %>%
            activate_mass_dataset('sample_info') %>%
            dplyr::select('sample_id') %>%
            dplyr::left_join(prj_init$sample_info)
        } else if (!is.null(data_clean_rv$object_merge)) {
          p3_enrich$object_merge = data_clean_rv$object_merge
        } else {return()}

        ## para
        temp_method =input$enrich_method %>% as.character()
        p3_enrich$enrich_db = input$enrich_db %>% as.character()
        p3_enrich$enrich_org = input$enrich_org %>% as.character()
        if(temp_method == "from new upload annotation file") {
          temp_upload_file = temp_anno_file
        } else {
          temp_upload_file = p3_enrich$object_merge %>% extract_annotation_table()
        }

        if("KEGG.ID"%in%colnames(temp_upload_file) & "superclass"%in%colnames(temp_upload_file)) {
          n_kegg = temp_upload_file %>%
            dplyr::filter(length(KEGG.ID) > 1) %>% nrow()
          n_class = temp_upload_file %>%
            dplyr::filter(length(superclass) > 1) %>% nrow()
          if(n_kegg ==0 & n_class == 0) {
            output$class_file_check = renderUI({
              isolate(HTML(paste0(
                '<font color = blue> <b>Error: </b> </font> <font color=orange> ','No compounds were annotated in both database, please re-do the classification step.',' </font><br/>'
              )))
            })
            return()
          } else {
            output$class_file_check = renderUI({
              isolate(HTML(paste0(
                '<font color = blue> <b>Total KEGG annotated compounds: </b> </font> <font color=orange> ',n_kegg,' </font><br/>',
                '<font color = blue> <b>Total ClassyFire annotated compounds: </b> </font> <font color=orange> ',n_class,' </font><br/>'
              )))
            })

            if(p3_enrich$enrich_db == "KEGG" & n_kegg > 10) {
              c_id2g <- temp_upload_file %>%
                dplyr::select(variable_id,KEGG.ID) %>%
                dplyr::filter(str_detect(KEGG.ID,"^C")) %>%
                setNames(c('variable_id',"GENE"))
              temp_cid_uniq <- c_id2g %>% pull(GENE) %>% unique()
              p3_enrich$TERM2NAME = get_kegg_pathway_ont(ont = p3_enrich$enrich_org) %>%
                group_by(TERM) %>%
                slice_head(n = 1) %>%
                ungroup() %>%
                distinct()
              ##> t2g -
              p3_enrich$TERM2GENE = purrr::map_dfr(
                .x = 1:length(temp_cid_uniq),.f = function(.x){
                  get_compound_info(compound_id = temp_cid_uniq[.x])
                }
              ) %>% inner_join(data.frame(TERM = p3_enrich$TERM2NAME %>% pull(TERM))) %>%
                left_join(c_id2g,by="GENE") %>%
                select(TERM,variable_id) %>%
                setNames(c("TERM","GENE")) %>%
                distinct()

            } else if(p3_enrich$enrich_db == "ClassyFire" & n_class > 10) {
              db <-
                temp_upload_file %>%
                filter(!is.na(superclass)) %>%
                select(variable_id,Compound.name,superclass,class,subclass,parent_levels) %>%
                separate_longer_delim(cols = parent_levels,delim = " | ") %>% ## seprate parent levels
                pivot_longer(!all_of(c("variable_id","Compound.name")),names_to = "levels",values_to = "NAME") %>%
                filter(!is.na(NAME)) %>%
                distinct()
              p3_enrich$TERM2NAME =
                db %>%
                select(NAME) %>%
                distinct() %>%
                mutate(TERM = paste0("MC:",str_pad(c(1:nrow(.)),5,"left",'0'))) %>% ## add terms
                filter(NAME != "NA") %>%
                select(TERM,NAME) %>%
                distinct()
              p3_enrich$TERM2GENE =
                db %>%
                inner_join(p3_enrich$TERM2NAME) %>%
                select(TERM,variable_id) %>%
                setNames(c("TERM","COMPOUND")) %>%
                distinct()
            }

            output$class_tbl_t2n = renderDataTable_formated(
              actions = input$enrich_check,
              condition1 = p3_enrich$TERM2NAME,
              tbl = p3_enrich$TERM2NAME
            )
            output$class_tbl_t2g = renderDataTable_formated(
              actions = input$enrich_check,
              condition1 = p3_enrich$TERM2GENE,
              tbl = p3_enrich$TERM2GENE
            )
          }

        } else {
          output$class_file_check = renderUI({
            isolate(HTML(paste0(
              '<font color = red> <b>Error: </b> </font> <font color=salmon>','no kegg or classyfire information in your file. Please upload the corrected file', '</font><br/>'
            )

            ))
          })
          return()
        }
      }
    )

    observeEvent(
      input$enrich_start,
      {
        if(is.null(p3_enrich$TERM2NAME)) {return()}
        if(is.null(p3_enrich$TERM2GENE)) {return()}

        p3_enrich$enrich_pval = input$enrich_pval %>% as.numeric()
        p3_enrich$enrich_qval = input$enrich_qval %>% as.numeric()
        p3_enrich$enrich_num = input$enrich_num %>% as.numeric()
        p3_enrich$enrich_plot_type = input$enrich_plot_type %>% as.character()
        p3_enrich$db_method = input$kegg_db_method %>% as.character()
        p3_enrich$enrich_org = input$enrich_org %>% as.character()

        p3_enrich$glist = input$enrich_glist %>% as.character() %>%
          str_split(pattern  = "\\\n|,| ",n = Inf,simplify = F) %>% unlist()

        print(p3_enrich$glist)

        p3_enrich$res = clusterProfiler::enricher(
          gene = p3_enrich$glist,
          pvalueCutoff = p3_enrich$enrich_pval,
          qvalueCutoff = p3_enrich$enrich_qval,
          TERM2GENE = p3_enrich$TERM2GENE,
          TERM2NAME = p3_enrich$TERM2NAME
        )


        p3_enrich$res_tbl = as.data.frame(p3_enrich$res)
        output$enrich_table = renderDataTable_formated(
          actions = input$enrich_start,
          condition1 = p3_enrich$res_tbl,
          tbl = p3_enrich$res_tbl
        )
        print('check point2')
        output$enrich_plot <- renderUI({
          plot_type <- input$enrich_plt_format
          if (plot_type) {
            plotlyOutput(outputId = ns("enrich_p_plotly"))
          } else {
            plotOutput(outputId = ns("enrich_p_plot"))
          }
        })

        output$enrich_p_plot <- renderPlot({
          if(is.null( p3_enrich$res )){return()}
          if(nrow(p3_enrich$res_tbl)==0){return()}
          if(p3_enrich$enrich_plot_type == "bar") {
            barplot(p3_enrich$res, showCategory=p3_enrich$enrich_num)
          } else {
            dotplot(p3_enrich$res, showCategory=p3_enrich$enrich_num)
          }
        })


        output$enrich_p_plotly <- renderPlotly({
          if(is.null( p3_enrich$res )){return()}
          if(nrow(p3_enrich$res_tbl)==0){return()}
          if(p3_enrich$enrich_plot_type == "bar") {
            ggplotly(barplot(p3_enrich$res, showCategory=p3_enrich$enrich_num))
          } else {
            ggplotly(dotplot(p3_enrich$res, showCategory=p3_enrich$enrich_num))
          }
        })
      }
    )

    # download ----------------------------------------------------------------
    ###> fig1 =====
    output$fig1_download = downloadHandler(
      filename = function() {
        paste0("Enrichment.", download_para()$fig1_format)
      },
      content = function(file) {
        # extract parameters
        para_d <- download_para()

        # draw condition
        if(p3_enrich$enrich_plot_type == "bar") {
         p =  barplot(p3_enrich$res, showCategory=p3_enrich$enrich_num)
        } else {
         p =  dotplot(p3_enrich$res, showCategory=p3_enrich$enrich_num)
        }
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


  })
}

