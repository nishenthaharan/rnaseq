options(shiny.maxRequestSize = 500 * 1024^2)

required_packages <- c(
  "shiny", "bslib", "DT", "DESeq2", "SummarizedExperiment",
  "ggplot2", "ggrepel", "pheatmap", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop(
    "Install missing dependencies before running the app: ",
    paste(missing_packages, collapse = ", "),
    ". Run Rscript install_dependencies.R from the repository root.",
    call. = FALSE
  )
}

invisible(lapply(
  c("utils.R", "io.R", "demo_data.R", "analysis.R", "plots.R"),
  function(file) source(file.path("R", file), local = FALSE)
))

app_theme <- bslib::bs_theme(
  version = 5,
  bg = "#F5F7FA",
  fg = "#182230",
  primary = "#0E7490",
  secondary = "#526075",
  success = "#18864B",
  danger = "#C93448",
  base_font = bslib::font_google("Inter"),
  heading_font = bslib::font_google("Manrope")
)

analysis_card <- function(title, ..., full_screen = TRUE) {
  bslib::card(
    bslib::card_header(title),
    full_screen = full_screen,
    ...
  )
}

ui <- bslib::page_navbar(
  title = shiny::div(
    class = "brand-lockup",
    shiny::span(class = "brand-mark", "R"),
    shiny::div(
      shiny::span(class = "brand-title", "RNASeq Explorer"),
      shiny::span(class = "brand-subtitle", "Differential expression in R")
    )
  ),
  theme = app_theme,
  fillable = FALSE,
  header = shiny::tags$head(
    shiny::tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  sidebar = bslib::sidebar(
    width = 340,
    open = "desktop",
    shiny::h5("Dataset"),
    shiny::checkboxInput("use_demo", "Use the built-in simulated dataset", value = TRUE),
    shiny::conditionalPanel(
      condition = "!input.use_demo",
      shiny::fileInput(
        "counts_file",
        "Raw count matrix",
        accept = c(".csv", ".tsv", ".txt"),
        buttonLabel = "Choose counts"
      ),
      shiny::fileInput(
        "metadata_file",
        "Sample metadata",
        accept = c(".csv", ".tsv", ".txt"),
        buttonLabel = "Choose metadata"
      )
    ),
    shiny::hr(),
    shiny::h5("Comparison"),
    shiny::selectInput("target", "Target condition (numerator)", choices = c("Treated", "Control")),
    shiny::selectInput("reference", "Reference condition", choices = c("Control", "Treated")),
    shiny::checkboxInput("use_batch", "Adjust for batch when available", value = FALSE),
    shiny::hr(),
    shiny::h5("Statistical thresholds"),
    shiny::numericInput("min_count", "Minimum raw count", value = 10, min = 1, step = 1),
    shiny::numericInput("min_samples", "Minimum samples passing count", value = 3, min = 1, step = 1),
    shiny::numericInput("alpha", "Adjusted P-value (FDR)", value = 0.05, min = 0.001, max = 0.5, step = 0.01),
    shiny::numericInput("lfc", "Absolute log2 fold-change", value = 1, min = 0, step = 0.25),
    shiny::actionButton("run_analysis", "Run DESeq2 analysis", class = "btn-primary run-button"),
    shiny::uiOutput("analysis_status")
  ),

  bslib::nav_panel(
    "Overview",
    shiny::div(
      class = "page-intro",
      shiny::h2("RNA-seq differential-expression workspace"),
      shiny::p(
        "Quality control, DESeq2 modelling, effect-size shrinkage, and publication-ready visualisation in one reproducible workflow."
      )
    ),
    shiny::uiOutput("metric_cards"),
    bslib::layout_columns(
      col_widths = c(7, 5),
      analysis_card("Principal component analysis", shiny::plotOutput("overview_pca", height = 420)),
      analysis_card("Analysis specification", shiny::uiOutput("analysis_specification"))
    )
  ),

  bslib::nav_panel(
    "Quality control",
    bslib::layout_columns(
      col_widths = c(6, 6),
      analysis_card("PCA", shiny::plotOutput("pca_plot", height = 480)),
      analysis_card("Library size", shiny::plotOutput("library_plot", height = 480)),
      analysis_card("Sample correlation", shiny::plotOutput("correlation_plot", height = 560)),
      analysis_card("Dispersion estimates", shiny::plotOutput("dispersion_plot", height = 560))
    )
  ),

  bslib::nav_panel(
    "Differential expression",
    bslib::layout_columns(
      col_widths = c(6, 6),
      analysis_card("Volcano plot", shiny::plotOutput("volcano_plot", height = 520)),
      analysis_card("MA plot", shiny::plotOutput("ma_plot", height = 520))
    ),
    analysis_card(
      "Differential-expression results",
      shiny::p(class = "card-note", "Click a column heading to sort. Use the search field to find a gene ID."),
      DT::DTOutput("results_table")
    )
  ),

  bslib::nav_panel(
    "Expression explorer",
    bslib::layout_columns(
      col_widths = c(8, 4),
      analysis_card("Top-gene heatmap", shiny::plotOutput("heatmap_plot", height = 720)),
      shiny::div(
        analysis_card(
          "Heatmap controls",
          shiny::sliderInput("heatmap_genes", "Number of genes", min = 10, max = 60, value = 30, step = 5),
          shiny::p(class = "card-note", "Genes are ranked by adjusted P-value and displayed as row-scaled, variance-stabilised expression.")
        ),
        analysis_card(
          "Single-gene expression",
          shiny::selectizeInput("gene_id", "Gene ID", choices = NULL),
          shiny::plotOutput("gene_plot", height = 400)
        )
      )
    )
  ),

  bslib::nav_panel(
    "Input data",
    bslib::layout_columns(
      col_widths = c(5, 7),
      analysis_card("Sample metadata", DT::DTOutput("metadata_table")),
      analysis_card("Count-matrix preview", DT::DTOutput("count_preview"))
    ),
    analysis_card(
      "Required file format",
      shiny::tags$ul(
        shiny::tags$li("Counts: the first column contains unique gene IDs; remaining columns contain non-negative raw integer counts."),
        shiny::tags$li("Metadata: required columns are sample_id and condition; batch is optional."),
        shiny::tags$li("Count-matrix sample columns and metadata sample_id values must match exactly."),
        shiny::tags$li("TPM, FPKM, CPM, ratios, and pre-normalised expression are not valid DESeq2 inputs.")
      )
    )
  ),

  bslib::nav_panel(
    "Downloads",
    shiny::div(
      class = "page-intro compact",
      shiny::h2("Export reproducible results"),
      shiny::p("Download complete statistics, significant genes, normalised counts, plots, or the full R analysis object.")
    ),
    bslib::layout_column_wrap(
      width = 1 / 3,
      analysis_card(
        "All results",
        shiny::p("Every tested gene with effect size, uncertainty, P-value, FDR, and classification."),
        shiny::downloadButton("download_all", "Download CSV", class = "btn-outline-primary"),
        full_screen = FALSE
      ),
      analysis_card(
        "Significant genes",
        shiny::p("Genes passing both the selected FDR and absolute log2 fold-change thresholds."),
        shiny::downloadButton("download_significant", "Download CSV", class = "btn-outline-primary"),
        full_screen = FALSE
      ),
      analysis_card(
        "Normalised counts",
        shiny::p("DESeq2 size-factor-normalised expression for all genes retained after filtering."),
        shiny::downloadButton("download_normalised", "Download CSV", class = "btn-outline-primary"),
        full_screen = FALSE
      ),
      analysis_card(
        "Publication plots",
        shiny::p("High-resolution PNG files for PCA, library size, volcano, MA, and heatmap views."),
        shiny::downloadButton("download_plots", "Download ZIP", class = "btn-outline-primary"),
        full_screen = FALSE
      ),
      analysis_card(
        "Analysis object",
        shiny::p("Full R object containing DESeq2 results, transformed data, parameters, and provenance."),
        shiny::downloadButton("download_rds", "Download RDS", class = "btn-outline-primary"),
        full_screen = FALSE
      ),
      analysis_card(
        "Reproducibility",
        shiny::verbatimTextOutput("session_details"),
        full_screen = FALSE
      )
    )
  ),

  bslib::nav_spacer(),
  bslib::nav_item(shiny::a("Documentation", href = "https://bioconductor.org/packages/DESeq2", target = "_blank"))
)

server <- function(input, output, session) {
  state <- shiny::reactiveValues(analysis = NULL, error = NULL, started_at = NULL)

  dataset <- shiny::reactive({
    if (isTRUE(input$use_demo)) {
      generate_demo_bundle()
    } else {
      shiny::req(input$counts_file, input$metadata_file)
      load_uploaded_bundle(input$counts_file, input$metadata_file)
    }
  })

  shiny::observe({
    bundle <- tryCatch(dataset(), error = function(error) NULL)
    shiny::req(bundle)
    conditions <- unique(bundle$metadata$condition)
    current_reference <- shiny::isolate(input$reference)
    current_target <- shiny::isolate(input$target)
    reference <- if (current_reference %in% conditions) current_reference else conditions[1L]
    target <- if (current_target %in% setdiff(conditions, reference)) current_target else setdiff(conditions, reference)[1L]
    shiny::updateSelectInput(session, "reference", choices = conditions, selected = reference)
    shiny::updateSelectInput(session, "target", choices = conditions, selected = target)
  })

  shiny::observeEvent(input$run_analysis, {
    state$error <- NULL
    state$analysis <- NULL
    state$started_at <- Sys.time()

    tryCatch(
      shiny::withProgress(message = "Running RNA-seq analysis", value = 0, {
        shiny::incProgress(0.15, detail = "Validating counts and experimental design")
        bundle <- dataset()
        shiny::incProgress(0.25, detail = "Estimating size factors and dispersions")
        analysis <- run_deseq_analysis(
          bundle = bundle,
          reference = input$reference,
          target = input$target,
          min_count = input$min_count,
          min_samples = input$min_samples,
          alpha = input$alpha,
          lfc_threshold = input$lfc,
          use_batch = isTRUE(input$use_batch)
        )
        shiny::incProgress(0.5, detail = "Preparing transformed data and visualisations")
        state$analysis <- analysis
        available_genes <- analysis$result$gene_id
        shiny::updateSelectizeInput(
          session,
          "gene_id",
          choices = available_genes,
          selected = available_genes[1L],
          server = TRUE
        )
      }),
      error = function(error) {
        state$error <- conditionMessage(error)
        shiny::showNotification(state$error, type = "error", duration = 12)
      }
    )
  }, ignoreInit = FALSE, ignoreNULL = FALSE)

  current_analysis <- shiny::reactive({
    message <- state$error %||% "Run the analysis to populate this view."
    shiny::validate(shiny::need(!is.null(state$analysis), message))
    state$analysis
  })

  output$analysis_status <- shiny::renderUI({
    if (!is.null(state$error)) {
      return(shiny::div(class = "status-message status-error", state$error))
    }
    if (is.null(state$analysis)) {
      return(shiny::div(class = "status-message", "Ready to analyse."))
    }
    shiny::div(
      class = "status-message status-success",
      sprintf("Complete: %s", state$analysis$parameters$comparison)
    )
  })

  output$metric_cards <- shiny::renderUI({
    analysis <- current_analysis()
    metrics <- analysis_metrics(analysis)
    subtitles <- c(
      "included in the model",
      "after count filtering",
      sprintf("FDR < %s and |LFC| >= %s", analysis$parameters$alpha, analysis$parameters$absolute_log2_fc),
      "higher in target",
      "lower in target"
    )
    classes <- c("metric-blue", "metric-slate", "metric-purple", "metric-red", "metric-cyan")

    cards <- lapply(seq_len(nrow(metrics)), function(index) {
      shiny::div(
        class = paste("metric-card", classes[index]),
        shiny::span(class = "metric-title", metrics$metric[index]),
        shiny::strong(class = "metric-value", format_count(metrics$value[index])),
        shiny::span(class = "metric-subtitle", subtitles[index])
      )
    })
    do.call(bslib::layout_column_wrap, c(list(width = 1 / 5), cards))
  })

  output$analysis_specification <- shiny::renderUI({
    analysis <- current_analysis()
    parameters <- analysis$parameters
    shiny::tags$dl(
      class = "specification-list",
      shiny::tags$dt("Data source"), shiny::tags$dd(parameters$source),
      shiny::tags$dt("Comparison"), shiny::tags$dd(parameters$comparison),
      shiny::tags$dt("Design"), shiny::tags$dd(parameters$design),
      shiny::tags$dt("Gene filter"), shiny::tags$dd(sprintf(
        ">= %d counts in >= %d samples", parameters$min_count, parameters$min_samples
      )),
      shiny::tags$dt("Decision rule"), shiny::tags$dd(sprintf(
        "FDR < %s and |log2FC| >= %s", parameters$alpha, parameters$absolute_log2_fc
      )),
      shiny::tags$dt("Genes retained"), shiny::tags$dd(sprintf(
        "%s of %s", format_count(parameters$genes_after_filter), format_count(parameters$genes_before_filter)
      )),
      shiny::tags$dt("LFC shrinkage"), shiny::tags$dd(parameters$lfc_shrinkage)
    )
  })

  output$overview_pca <- shiny::renderPlot(plot_pca(current_analysis()), res = 110)
  output$pca_plot <- shiny::renderPlot(plot_pca(current_analysis()), res = 110)
  output$library_plot <- shiny::renderPlot(plot_library_sizes(current_analysis()), res = 110)
  output$correlation_plot <- shiny::renderPlot(plot_sample_correlation(current_analysis()), res = 110)
  output$dispersion_plot <- shiny::renderPlot(plot_dispersion(current_analysis()), res = 110)
  output$volcano_plot <- shiny::renderPlot(plot_volcano(current_analysis()), res = 110)
  output$ma_plot <- shiny::renderPlot(plot_ma(current_analysis()), res = 110)
  output$heatmap_plot <- shiny::renderPlot(
    plot_top_heatmap(current_analysis(), input$heatmap_genes),
    res = 110
  )
  output$gene_plot <- shiny::renderPlot({
    shiny::req(input$gene_id)
    plot_gene_expression(current_analysis(), input$gene_id)
  }, res = 110)

  output$results_table <- DT::renderDT({
    table <- current_analysis()$result
    DT::datatable(
      table,
      rownames = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv"),
        order = list(list(6, "asc"))
      ),
      class = "compact stripe hover"
    ) |>
      DT::formatSignif(c("base_mean", "log2_fold_change", "lfc_se", "statistic"), digits = 4) |>
      DT::formatSignif(c("p_value", "adjusted_p_value"), digits = 3)
  })

  output$metadata_table <- DT::renderDT({
    metadata <- dataset()$metadata
    rownames(metadata) <- NULL
    DT::datatable(metadata, rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  output$count_preview <- DT::renderDT({
    counts <- dataset()$counts
    preview <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
    DT::datatable(
      utils::head(preview, 25L),
      rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 10)
    )
  })

  output$session_details <- shiny::renderText({
    provenance <- current_analysis()$provenance
    package_lines <- paste(names(provenance$packages), provenance$packages, sep = " ", collapse = "\n")
    paste(
      provenance$analysed_at,
      provenance$r_version,
      provenance$platform,
      package_lines,
      sep = "\n"
    )
  })

  comparison_stub <- function(analysis) {
    safe_filename(sprintf("%s-vs-%s", analysis$parameters$target, analysis$parameters$reference))
  }

  output$download_all <- shiny::downloadHandler(
    filename = function() sprintf("deseq2-all-%s.csv", comparison_stub(current_analysis())),
    content = function(file) utils::write.csv(current_analysis()$result, file, row.names = FALSE)
  )

  output$download_significant <- shiny::downloadHandler(
    filename = function() sprintf("deseq2-significant-%s.csv", comparison_stub(current_analysis())),
    content = function(file) utils::write.csv(significant_results(current_analysis()), file, row.names = FALSE)
  )

  output$download_normalised <- shiny::downloadHandler(
    filename = function() sprintf("normalised-counts-%s.csv", comparison_stub(current_analysis())),
    content = function(file) {
      matrix <- current_analysis()$normalised_counts
      table <- data.frame(gene_id = rownames(matrix), matrix, check.names = FALSE)
      utils::write.csv(table, file, row.names = FALSE)
    }
  )

  output$download_rds <- shiny::downloadHandler(
    filename = function() sprintf("rnaseq-analysis-%s.rds", comparison_stub(current_analysis())),
    content = function(file) saveRDS(current_analysis(), file, compress = "xz")
  )

  output$download_plots <- shiny::downloadHandler(
    filename = function() sprintf("rnaseq-plots-%s.zip", comparison_stub(current_analysis())),
    content = function(file) {
      analysis <- current_analysis()
      plot_directory <- tempfile("rnaseq-plots-")
      dir.create(plot_directory, recursive = TRUE)
      on.exit(unlink(plot_directory, recursive = TRUE), add = TRUE)

      ggplot2::ggsave(file.path(plot_directory, "pca.png"), plot_pca(analysis), width = 8, height = 6, dpi = 300)
      ggplot2::ggsave(file.path(plot_directory, "library-size.png"), plot_library_sizes(analysis), width = 8, height = 6, dpi = 300)
      ggplot2::ggsave(file.path(plot_directory, "volcano.png"), plot_volcano(analysis), width = 8, height = 6, dpi = 300)
      ggplot2::ggsave(file.path(plot_directory, "ma-plot.png"), plot_ma(analysis), width = 8, height = 6, dpi = 300)

      grDevices::png(file.path(plot_directory, "top-gene-heatmap.png"), width = 2400, height = 2600, res = 300)
      plot_top_heatmap(analysis, input$heatmap_genes)
      grDevices::dev.off()

      previous_directory <- setwd(plot_directory)
      on.exit(setwd(previous_directory), add = TRUE)
      utils::zip(file, files = list.files(plot_directory), flags = "-j")
    },
    contentType = "application/zip"
  )
}

shiny::shinyApp(ui, server)
