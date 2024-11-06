#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinythemes shinytheme
#' @importFrom shinyjs useShinyjs
#' @importFrom bsicons bs_icon
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    navbarPage(
      theme = shinytheme('spacelab'),
      customLogo(version = "V.1.0.1"),
      homepage_ui("home_id"),
      project_init_ui("project_init_id"),
      navbarMenu(
        title = 'Import file',
        icon = bs_icon("upload"),
        data_import_raw_ui("data_import_raw_id"),
        data_import_tbl_ui("data_import_tbl_id"),
        data_import_massdataset_ui("data_import_massdataset_id")
      ),
      navbarMenu(
        title = 'Data Cleaning',
        icon = bs_icon("filter"),
        data_overview_ui("data_overview_id"),
        data_rm_noise_ui("data_rm_noise_id"),
        data_rm_outlier_ui("data_rm_outlier_id"),
        data_mv_impute_ui("data_mv_impute_id"),
        data_normalize_ui("data_normalize_id")
      ),
      navbarMenu(
        title = 'Downstream Analysis',
        icon = bs_icon("bar-chart-line"),
        annotation_ui("annotation_id"),
        annotation_filter_ui("annotation_filter_id"),
        data_merge_ui("data_merge_id"),
        data_class_ui("data_class_id"),
        dam_res_ui("dam_res_id"),
        data_enrich_ui("data_enrich_id")
      )
    )
  )
}


#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem use_external_css_file  add_resource_path activate_js favicon bundle_resources add_css_file
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )
#  use_external_css_file("styles.css")

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "MetMiner"
    )

    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )

}
