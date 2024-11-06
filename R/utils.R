#' getvolumnes for shinyFiles in different OS.
#'
#' Creates a custom styled horizontal rule.
#' @param style A string indicating the style of the horizontal rule.
#' @return An shinyFiles object for volume.
#' @importFrom shinyFiles getVolumes
#' @noRd

get_volumes <- function() {
  if (Sys.info()["sysname"] == "Windows") {
    return(getVolumes_win())
  } else {
    return(shinyFiles::getVolumes()())
  }
}

#' load_depends
#'
#' load depends of metminer
#' @export
#' @noRd

load_depends_metminer <- function() {
  library(tidyverse)
  library(tidymass)
  library(MDAtoolkits)
}

#' Custom Horizontal Rule
#'
#' Creates a custom styled horizontal rule.
#' @param style A string indicating the style of the horizontal rule.
#' @return An HTML horizontal rule element with custom style.
#' @importFrom shiny tags
#' @noRd

custom_hr <- function(style) {
  tags$hr(style = style)
}

#' Main Horizontal Rule
#'
#' This function creates the main styled horizontal rule used in the UI.
#' @return An HTML horizontal rule element for main style.
#' @noRd

hr_main <- function() {
  custom_hr("border-top: 6px double #008080; border-bottom: 3px solid #008080;")
}

#' Bar Horizontal Rule
#'
#' This function creates a bar styled horizontal rule.
#' @return An HTML horizontal rule element for bar style.
#' @noRd

hr_bar <- function() {
  custom_hr("border-top: 3px dotted #008080;")
}

#' Head Horizontal Rule
#'
#' This function creates a header styled horizontal rule.
#' @return An HTML horizontal rule element for header style.
#' @noRd

hr_head <- function() {
  custom_hr("border: 0; padding-top: 1.5px; background: linear-gradient(to right, transparent, #008080, transparent);")
}


#' ggplot2 themes
#'
#' ggplot2 theme for plots.
#' @return An ggplot2 theme object.
#' @importFrom ggplot2 theme element_rect element_text
#' @noRd
theme1 <- function(){
  temp = theme(
    panel.border = element_rect(linewidth  = 1.5),
    axis.title = element_text(size = 14,color = 'black'),
    axis.text = element_text(size = 12,color = 'black')
  )
  return(temp)
  }

#' website logo
#'
#' MetMiner logo.
#' @return shinyDashboardLogoDIY object.
#' @param version version of metminer
#' @importFrom dashboardthemes shinyDashboardLogoDIY
#' @noRd
#'
customLogo <- function(version) {
  shinyDashboardLogoDIY(
    boldText = "Zhang Lab",
    mainText = "MetMiner",
    textSize = 14,
    badgeText = version,
    badgeTextColor = "white",
    badgeTextSize = 2,
    badgeBackColor = "#40E0D0",
    badgeBorderRadius = 3
  )
}

#' adduct extract
#'
#' to extract all adducts from tidymass.
#' @return A string contains clean adducts.
#' @param string adducts string
#' @importFrom stringr str_replace_all
#' @noRd
#'

re_form_reg = function(string){
  string %>% str_replace_all("\\(","\\\\(") %>% str_replace_all("\\)","\\\\)") %>% str_replace_all("\\-","\\\\-")%>% str_replace_all("\\+","\\\\+")
}

#' reformed textInput
#'
#' to extract all adducts from tidymass.
#' @return A string contains clean adducts.
#' @param inpuitID see `shiny::textInput`
#' @param label see `shiny::textInput`
#' @param value see `shiny::textInput`
#' @param placeholder see `shiny::textInput`
#' @param title hover_text of this input box
#' @importFrom stringr str_replace_all
#' @noRd
#'

textInput_div = function(inputId,label,value,placeholder,title) {
  div(
    textInput(
      inputId = inputId,
      label = label,
      value = value,
      placeholder = placeholder
    ),
    title = title
  )
}

#' reformed selectInput_div
#'
#' to extract all adducts from tidymass.
#' @return A string contains clean adducts.
#' @param inpuitID see `shiny::selectInput`
#' @param label see `shiny::selectInput`
#' @param choices see `shiny::selectInput`
#' @param selected see `shiny::selectInput`
#' @param multiple see `shiny::selectInput`
#' @param title hover_text of this input box
#' @importFrom stringr str_replace_all
#' @noRd
#'

selectInput_div = function(inputId,label,choices,selected,multiple,title) {
  div(
    selectInput(
      inputId = inputId,
      label = label,
      choices = choices,
      selected = selected,
      multiple = multiple,
    ),
    title = title
  )
}

#' fetch kegg pathway via ont
#'
#' Get pathway id via KEGG API.
#' @return A dataframe to get term2name.
#' @param ont species id of KEGG database.
#' @importFrom stringr str_replace str_extract
#' @importFrom dplyr mutate
#' @importFrom RCurl getURL
#' @export
#'

get_kegg_pathway_ont <- function(ont) {
  temp_t2n_url = getURL(paste0("https://rest.kegg.jp/list/pathway/",ont))
  TERM2NAME = read.table(text = temp_t2n_url,sep = "\t",header = F)  |>
    setNames(c("TERM","NAME")) |>
    mutate(TERM = str_replace(TERM,ont,"map")) |>
    mutate(NAME = str_extract(NAME, "^.*?(?= -)"))
  return(TERM2NAME)
}


#' fetch kegg compound id belongs to corresponding pathway
#'
#' match kegg cid 2 pathway.
#' @return A dataframe to get term2gene
#' @param compound_id kegg cid.
#' @importFrom stringr str_extract_all
#' @importFrom RCurl getURL
#' @export
#'

get_compound_info <- function(compound_id) {
  base_url <- "https://rest.kegg.jp/get/cpd:"

  url <- paste0(base_url, compound_id)


  result <- tryCatch({
    response <- RCurl::getURL(url)
    pathways <- stringr::str_extract_all(response, "map[0-9]+")[[1]]

    # check
    if (length(pathways) == 0) {
      stop("No pathways found")
    }

    TERM2GENE <- data.frame(
      TERM = pathways,
      GENE = compound_id
    )

    return(TERM2GENE)
  }, error = function(e) {
    # return NA
    message("Error occurred: ", e$message)
    return(data.frame(
      TERM = NA,
      GENE = compound_id
    ))
  })

  return(result)
}
