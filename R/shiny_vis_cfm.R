#' Map CFM energy labels to collision energy display names
#' @noRd
.cfm_energy_to_ce_label <- function(energy) {
  dplyr::case_when(
    energy == "energy0" ~ "CE=10",
    energy == "energy1" ~ "CE=20",
    energy == "energy2" ~ "CE=40",
    TRUE ~ as.character(energy)
  )
}

#' Best-scoring fragment assignment per peak (energy + mz)
#' @noRd
.cfm_peak_best_fragment <- function(cfmd) {
  peaks <- cfmd@peak_assignment
  if (is.null(peaks) || !nrow(peaks)) {
    return(peaks[0, , drop = FALSE])
  }
  peaks %>%
    dplyr::group_by(energy, mz) %>%
    dplyr::slice_max(order_by = fragment_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ce_label = factor(
        .cfm_energy_to_ce_label(energy),
        levels = c("CE=10", "CE=20", "CE=40")
      ),
      hover_label = paste0(
        "mz: ", round(mz, 4),
        "\nint: ", signif(intensity, 4),
        "\n", ce_label,
        ifelse(
          is.na(fragment_id),
          "\n(no fragment annotation)",
          paste0("\n", fragment_id)
        )
      )
    )
}

#' Precursor fragment row (highest m/z; fall back to first row)
#' @noRd
.cfm_precursor_fragment <- function(cfmd) {
  frags <- cfmd@fragment_define
  if (is.null(frags) || !nrow(frags)) {
    return(NULL)
  }
  if ("fragment_mz" %in% colnames(frags)) {
    i <- which.max(frags$fragment_mz)
    return(frags[i, , drop = FALSE])
  }
  frags[1, , drop = FALSE]
}

#' Molecular formula from SMILES (NA on failure)
#' @noRd
.cfm_smiles_formula <- function(smiles) {
  if (is.null(smiles) || length(smiles) == 0 || is.na(smiles) || !nzchar(smiles)) {
    return(NA_character_)
  }
  tryCatch({
    sdf <- get_smiles_sdf(smiles)[[1]]
    unname(ChemmineR::MF(sdf, addH = TRUE))
  }, error = function(e) NA_character_)
}

#' Empty visNetwork placeholder
#' @noRd
.shiny_vis_cfm_empty <- function(label = "No fragment selected") {
  visNetwork::visNetwork(
    nodes = data.frame(
      id = 1,
      label = label,
      shape = "text",
      font.size = 18,
      stringsAsFactors = FALSE
    ),
    edges = data.frame()
  )
}

#' Build vertical stem coordinates for a spectrum panel
#' @noRd
.cfm_stem_xy <- function(mz, intensity) {
  n <- length(mz)
  if (!n) {
    return(data.frame(x = numeric(), y = numeric()))
  }
  data.frame(
    x = as.numeric(rbind(mz, mz, rep(NA_real_, n))),
    y = as.numeric(rbind(rep(0, n), intensity, rep(NA_real_, n)))
  )
}

#' Plotly multi-CE spectrum for CFM_data
#'
#' @param cfmd A \code{CFM_data} object.
#' @param source Plotly event source id (default \code{"cfm_spectra"}).
#' @return A plotly object with clickable peaks (customdata = fragment_id).
#' @export
plotly_CFM_spectra <- function(cfmd, source = "cfm_spectra") {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package plotly is required for plotly_CFM_spectra().", call. = FALSE)
  }
  if (!methods::is(cfmd, "CFM_data")) {
    stop("`cfmd` must be a CFM_data object.", call. = FALSE)
  }

  sp.data <- .cfm_peak_best_fragment(cfmd)
  if (!nrow(sp.data)) {
    return(
      plotly::plot_ly(source = source) %>%
        plotly::layout(
          xaxis = list(title = "m/z"),
          yaxis = list(title = "intensity")
        ) %>%
        plotly::event_register("plotly_click")
    )
  }

  ce_colors <- c(
    "CE=10" = "#4DBBD5",
    "CE=20" = "#FF7F0E",
    "CE=40" = "#E64B35"
  )
  ce_levels <- intersect(c("CE=10", "CE=20", "CE=40"), as.character(unique(sp.data$ce_label)))
  if (!length(ce_levels)) {
    ce_levels <- as.character(unique(sp.data$ce_label))
  }
  xmax <- max(sp.data$mz, na.rm = TRUE) * 1.08
  panels <- list()

  for (ce in ce_levels) {
    d <- sp.data %>% dplyr::filter(as.character(ce_label) == ce)
    col <- if (ce %in% names(ce_colors)) unname(ce_colors[[ce]]) else "#555555"
    ymax <- if (nrow(d)) max(d$intensity, na.rm = TRUE) * 1.15 else 1
    stems <- .cfm_stem_xy(d$mz, d$intensity)

    p <- plotly::plot_ly(source = source) %>%
      plotly::add_paths(
        data = stems,
        x = ~x,
        y = ~y,
        line = list(color = col, width = 1.5),
        hoverinfo = "skip",
        showlegend = FALSE
      ) %>%
      plotly::add_markers(
        data = d,
        x = ~mz,
        y = ~intensity,
        customdata = ~fragment_id,
        text = ~hover_label,
        hoverinfo = "text",
        marker = list(color = col, size = 7),
        showlegend = FALSE
      ) %>%
      plotly::layout(
        yaxis = list(title = ce, range = c(0, ymax), showgrid = TRUE),
        xaxis = list(range = c(0, xmax), showgrid = TRUE)
      )
    panels[[ce]] <- p
  }

  plotly::subplot(
    panels,
    nrows = length(panels),
    shareX = TRUE,
    titleY = TRUE,
    margin = 0.04
  ) %>%
    plotly::layout(
      dragmode = "zoom",
      xaxis = list(title = "m/z", range = c(0, xmax))
    ) %>%
    plotly::event_register("plotly_click")
}

#' Shiny app to visualize CFM_data spectra and fragments
#'
#' Layout (left | right):
#' \preformatted{
#' | p1 spectra | p2 precursor molecule |
#' | p1 spectra | p3 fragment molecule  |
#' | p1 spectra | p4 smiles / formula / mz |
#' }
#'
#' Click a peak in the multi-CE spectrum (p1) to update the fragment molecule (p3)
#' and the info panel (p4).
#'
#' @param cfmd A \code{CFM_data} object, typically from
#'   \code{\link{get_CFM_data_from_smiles}} or \code{\link{CFM_annotate_by_predict}}.
#' @param launch.browser Logical; forward to \code{shiny::shinyApp} options
#'   (default \code{TRUE}).
#' @param host Host binding for the app (default \code{"127.0.0.1"}).
#'
#' @return A \code{shiny.appobj} (invisibly when run interactively via the app).
#' @export
#'
#' @examples
#' \dontrun{
#' cfm <- get_CFM_data_from_smiles("NCC(O)=O", compound_id = "glycine")
#' shiny_vis_cfm(cfm)
#' }
shiny_vis_cfm <- function(cfmd,
                          launch.browser = TRUE,
                          host = "127.0.0.1") {
  for (pkg in c("shiny", "plotly", "visNetwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package ", pkg, " is required for shiny_vis_cfm().", call. = FALSE)
    }
  }
  if (!methods::is(cfmd, "CFM_data")) {
    stop("`cfmd` must be a CFM_data object.", call. = FALSE)
  }
  if (is.null(cfmd@peak_assignment) || !nrow(cfmd@peak_assignment)) {
    stop("`cfmd` has no peak_assignment data.", call. = FALSE)
  }

  precursor <- .cfm_precursor_fragment(cfmd)
  precursor_smiles <- if (!is.null(precursor)) precursor$smiles[[1]] else NA_character_
  precursor_mz <- if (!is.null(precursor)) precursor$fragment_mz[[1]] else NA_real_
  precursor_id <- if (!is.null(precursor)) precursor$fragment_id[[1]] else NA_character_
  precursor_formula <- .cfm_smiles_formula(precursor_smiles)

  frags <- cfmd@fragment_define
  peak_lookup <- .cfm_peak_best_fragment(cfmd)

  .ui <- shiny::fluidPage(
    shiny::titlePanel("CFM spectrum viewer"),
    shiny::fluidRow(
      shiny::column(
        width = 7,
        shiny::h4("Multi-CE spectra"),
        plotly::plotlyOutput("p1", height = "720px")
      ),
      shiny::column(
        width = 5,
        shiny::h4("Precursor"),
        visNetwork::visNetworkOutput("p2", height = "220px"),
        shiny::h4("Fragment"),
        visNetwork::visNetworkOutput("p3", height = "220px"),
        shiny::h4("Info"),
        shiny::verbatimTextOutput("p4")
      )
    )
  )

  .server <- function(input, output, session) {
    selected_fid <- shiny::reactiveVal(NA_character_)
    selected_peak <- shiny::reactiveVal(NULL)

    shiny::observeEvent(
      plotly::event_data("plotly_click", source = "cfm_spectra", priority = "event"),
      {
        ed <- plotly::event_data(
          "plotly_click",
          source = "cfm_spectra",
          priority = "event"
        )
        if (is.null(ed) || is.null(ed$x)) {
          return(invisible(NULL))
        }
        mz_click <- as.numeric(ed$x[[1]])
        # Prefer exact peak match on this subplot's y (intensity) when available
        hit <- peak_lookup
        if (nrow(hit)) {
          i <- which.min(abs(hit$mz - mz_click))
          row <- hit[i, , drop = FALSE]
          fid <- row$fragment_id[[1]]
          selected_peak(list(
            mz = row$mz[[1]],
            intensity = row$intensity[[1]],
            ce_label = as.character(row$ce_label[[1]]),
            fragment_id = if (is.na(fid)) NA_character_ else as.character(fid)
          ))
          selected_fid(if (is.na(fid)) NA_character_ else as.character(fid))
        }
      },
      ignoreNULL = TRUE
    )

    output$p1 <- plotly::renderPlotly({
      plotly_CFM_spectra(cfmd, source = "cfm_spectra")
    })

    output$p2 <- visNetwork::renderVisNetwork({
      if (is.na(precursor_smiles) || !nzchar(precursor_smiles)) {
        return(.shiny_vis_cfm_empty("No precursor SMILES"))
      }
      vis_smiles(precursor_smiles, show.formula = FALSE, show_id = FALSE)
    })

    output$p3 <- visNetwork::renderVisNetwork({
      fid <- selected_fid()
      pk <- selected_peak()
      if (is.null(pk)) {
        return(.shiny_vis_cfm_empty("Click a peak"))
      }
      if (is.na(fid) || !nzchar(fid) || !(fid %in% frags$fragment_id)) {
        return(.shiny_vis_cfm_empty("No fragment annotation"))
      }
      smiles <- frags$smiles[match(fid, frags$fragment_id)]
      if (is.na(smiles) || !nzchar(smiles)) {
        return(.shiny_vis_cfm_empty("No SMILES for fragment"))
      }
      vis_smiles(smiles, show.formula = FALSE, show_id = FALSE)
    })

    output$p4 <- shiny::renderText({
      fid <- selected_fid()
      pk <- selected_peak()
      frag_smiles <- NA_character_
      frag_mz <- NA_real_
      frag_formula <- NA_character_
      if (!is.na(fid) && nzchar(fid) && fid %in% frags$fragment_id) {
        row <- frags[match(fid, frags$fragment_id), , drop = FALSE]
        frag_smiles <- row$smiles[[1]]
        frag_mz <- row$fragment_mz[[1]]
        frag_formula <- .cfm_smiles_formula(frag_smiles)
      }

      peak_block <- if (is.null(pk)) {
        "peak:    (none selected)\n"
      } else {
        paste0(
          "peak mz: ", pk$mz, "\n",
          "peak int:", pk$intensity, "\n",
          "CE:      ", pk$ce_label, "\n"
        )
      }

      frag_block <- if (is.null(pk)) {
        "id:      (none selected)\nsmiles:  -\nformula: -\nmz:      -\n"
      } else if (is.na(fid) || !nzchar(fid)) {
        paste0(
          "id:      (no CFM fragment annotation)\n",
          "smiles:  -\n",
          "formula: -\n",
          "mz:      -\n",
          "note:    CFM predicted this peak but did not assign a structure\n"
        )
      } else {
        paste0(
          "id:      ", fid, "\n",
          "smiles:  ", ifelse(is.na(frag_smiles), "-", frag_smiles), "\n",
          "formula: ", ifelse(is.na(frag_formula), "-", frag_formula), "\n",
          "mz:      ", ifelse(is.na(frag_mz), "-", frag_mz), "\n"
        )
      }

      paste0(
        "=== Precursor ===\n",
        "id:      ", precursor_id, "\n",
        "smiles:  ", precursor_smiles, "\n",
        "formula: ", precursor_formula, "\n",
        "mz:      ", precursor_mz, "\n",
        "\n",
        "=== Selected peak ===\n",
        peak_block,
        "\n",
        "=== Fragment ===\n",
        frag_block
      )
    })
  }

  shiny::shinyApp(
    ui = .ui,
    server = .server,
    options = list(host = host, launch.browser = launch.browser)
  )
}
