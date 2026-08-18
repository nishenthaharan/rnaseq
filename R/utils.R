`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

abort_user <- function(message, call = NULL) {
  stop(structure(
    list(message = message, call = call),
    class = c("rnaseq_user_error", "error", "condition")
  ))
}

format_count <- function(x) {
  scales::comma(x, accuracy = 1)
}

format_percent <- function(x, accuracy = 0.1) {
  scales::percent(x, accuracy = accuracy)
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "-", x)
  gsub("-+", "-", x)
}

rnaseq_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#14213D"),
      plot.subtitle = ggplot2::element_text(colour = "#526075"),
      axis.title = ggplot2::element_text(face = "bold", colour = "#26354A"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E8EDF3", linewidth = 0.4),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(10, 14, 10, 10)
    )
}

condition_palette <- function(values) {
  values <- unique(as.character(values))
  colours <- scales::hue_pal(l = 55, c = 90)(length(values))
  stats::setNames(colours, values)
}

package_versions <- function(packages) {
  versions <- vapply(
    packages,
    function(package) {
      if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    },
    character(1)
  )
  versions[!is.na(versions)]
}

build_provenance <- function(parameters) {
  list(
    analysed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    r_version = R.version.string,
    platform = R.version$platform,
    parameters = parameters,
    packages = package_versions(c(
      "shiny", "DESeq2", "SummarizedExperiment", "ggplot2", "pheatmap", "DT"
    ))
  )
}
