#' pie plot for classification
#'
#' Pie plot.
#' @return A ggplot object of missing value.
#' @param x a MDAtoolkits classyfire result.
#' @param tag level of classyfire superclass class subclass
#' @param cut cut-off of small category
#' @param ntop show to n category
#' @importFrom dplyr select group_by summarise n case_when ungroup
#' @importFrom tidyr drop_na
#' @importFrom ggplot2 aes geom_bar coord_polar geom_text position_stack theme_void
#' @importFrom ggforce geom_arc_bar
#' @references see massdataset::show_sample_missing_values
#'
#' @noRd
#' @export
#'

pie_plot <- function(x, tag='superclass',cut = 10,ntop = 15){
  anno_class <-
    x %>%
    select(variable_id,superclass,class,subclass,parent_levels) %>%
    filter(superclass != "NA")

  anno_class_long <-
    anno_class %>%
    pivot_longer(!variable_id,names_to = 'levels',values_to = 'type') %>%
    drop_na() %>%
    filter(type != "NA") %>%
    filter(type != "") %>%
    arrange(levels)

  class_for_plot <-
    anno_class_long %>%
    group_by(levels, type) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    mutate(type = case_when(
      n <= cut ~ "other",
      row_number(desc(n)) >= ntop ~ "other",
      TRUE ~ type
    )) %>%
    group_by(levels, type) %>%
    mutate(n = sum(n)) %>%
    ungroup() %>%
    distinct() %>%
    group_by(levels) %>%
    mutate(
      percentage = n / sum(n) * 100,
      start_angle = cumsum(lag(percentage, default = 0)) / 100 * 2 * pi,
      end_angle = cumsum(percentage) / 100 * 2 * pi,
      angle = (start_angle + end_angle) / 2,
      ypos = 2.2 * cos(angle),
      xpos = 2.2 * sin(angle)
    )


  plt_tbl = class_for_plot %>% filter(levels == tag) %>%
    mutate(type = factor(type,levels = type))

  p =
    ggplot() +
    ggforce::geom_arc_bar(
      data = plt_tbl,
      stat = 'pie',
      mapping = aes(x0 = 0, y0 = 0, r0 = 0, r = 2, amount = n, fill = type)
    ) +
    geom_text(
      data = plt_tbl,
      aes(x = xpos, y = ypos, label = paste0(round(percentage, 1), "%")),
      size = 3, color = "black"  # 标签颜色和大小
    ) +
    coord_fixed() +
    theme_void() +
    labs(fill = tag) +
    theme(legend.position = "right")

  return(p)
}

#' pie plot for classification plotly version
#'
#' Pie plot - plotly.
#' @return A ggplot object of missing value.
#' @param x a MDAtoolkits classyfire result.
#' @param tag level of classyfire superclass class subclass
#' @param cut cut-off of small category
#' @param ntop show to n category
#' @importFrom dplyr select group_by summarise n case_when ungroup
#' @importFrom tidyr drop_na
#' @importFrom ggplot2 aes geom_bar coord_polar geom_text position_stack theme_void
#' @importFrom plotly plot_ly layout
#' @references see massdataset::show_sample_missing_values
#'
#' @noRd
#' @export
#'

pie_plot_plotly <- function(x, tag='superclass',cut = 10,ntop = 15){
  anno_class <-
    x %>%
    select(variable_id,superclass,class,subclass,parent_levels) %>%
    filter(superclass != "NA")

  anno_class_long <-
    anno_class %>%
    pivot_longer(!variable_id,names_to = 'levels',values_to = 'type') %>%
    drop_na() %>%
    filter(type != "NA") %>%
    filter(type != "") %>%
    arrange(levels)

  class_for_plot <-
    anno_class_long %>%
    group_by(levels, type) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    mutate(type = case_when(
      n <= cut ~ "Other",
      row_number(desc(n)) >= ntop ~ "Other",
      TRUE ~ type
    )) %>%
    group_by(levels, type) %>%
    mutate(n = sum(n)) %>%
    ungroup() %>%
    distinct() %>%
    group_by(levels)


  plt_tbl = class_for_plot %>% filter(levels == tag) %>%
    mutate(type = factor(type,levels = type)) %>%
    mutate(percent = (n/sum(n))*100)

  plot_ly(plt_tbl, labels = ~type, values = ~n, type = 'pie', textinfo = 'percent+value',
          textposition = 'inside', insidetextorientation = 'radial', hoverinfo = 'label+percent+value') %>%
    layout(title = paste('Pie Chart of', tag),
           xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

}
