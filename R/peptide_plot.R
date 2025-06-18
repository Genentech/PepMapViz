#' Create a peptide Plot
#'
#' This function generates a peptide plot using the provided data and allows for customization of the plot layout.
#'
#' @param data A dataframe containing the PSM data or Area data got from peptide_cluster_quantification.
#' @param domain A dataframe containing the domain data with columns including 'domain_start', 'domain_end', and 'domain_type' and condition columns.
#' @param y_axis_vars A list of variables for the donor and type facets.
#' @param x_axis_vars A list of variables for the region facets.
#' @param y_expand A numeric vector of length 2 specifying the expansion for the y-axis. Default is `c(0.1, 0.15)`.
#' @param x_expand A numeric vector of length 2 specifying the expansion for the x-axis. Default is `c(0.6, 0.6)`.
#' @param theme_options A list of additional theme options to customize the plot. Default is an empty list.
#' @param labs_options A list of additional labs options to customize the plot labels. Default is an empty list.
#' @param color_fill_column The name of the column in `data_with_psm` to be used for the fill aesthetic. Default is 'PSM'.
#' @param fill_gradient_options A list of options for `scale_fill_gradient`. Default is an empty list.
#' @param label_size The size of the labels in the plot. Default is 3.
#' @param add_domain A logical value indicating whether to add domain like CDR (Complementarity-Determining Region) to the plot. Default is TRUE.
#' @param domain A dataframe containing the domain data with columns including 'domain_start', 'domain_end', and 'domain_type'.
#' @param domain_start_column The name of the column in `domain` containing the start position of the domain Default is 'domain_start'.
#' @param domain_end_column The name of the column in `domain` containing the end position of the domain Default is 'domain_end'.
#' @param domain_type_column The name of the column in `domain` containing the type of the domain Default is 'domain_type'.
#' @param domain_border_color_column The name of the column in `domain` containing the border color of the domain Default is 'domain_color'.
#' @param domain_fill_color_column The name of the column in `domain` containing the fill color of the domain Default is 'domain_fill_color'.
#' @param add_domain_label Logical; whether to annotate the domain type as text above the domain rectangle. Default is TRUE.
#' @param domain_label_size Numeric; text size for the domain label. Default is 4.
#' @param domain_label_y_column The name of the column in `domain` containing y-axis position for the domain label. Default is 'domain_label_y'.
#' @param domain_label_color Character; color for domain label text. Default is 'black'.
#' @param PTM A logical value indicating whether to include PTM (Post-Translational Modification) data in the plot. Default is FALSE.
#' @param PTM_type_column The name of the column in `data_with_psm` containing the type of the PTM. Default is 'PTM_type'.
#' @param PTM_color A list of colors for the PTM types. Default is NULL.
#' @param add_label A logical value indicating whether to add labels to the plot. Default is TRUE.
#' @param label_column The name of the column in `data_with_psm` containing the labels to be added to the plot. Default is 'Character'.
#' @param label_filter A list of column names and their values to filter the data for the labels. Default is NULL.
#' @param label_y The position of y axis of the label.
#' @param column_order A list of column names and their order for the plot. Default is NULL.
#'
#' @return This function returns a ggplot object representing the PSM plot.
#'
#' @examples
#' data <- data.frame(
#'   Character = c("X", "Y", "Z", "A", "A", "A"),
#'   Position = 1:6,
#'   Condition_1 = rep("Drug1", 6),
#'   Region_2 = rep("Arm_1", 6),
#'   Area = c(0.000000, 0.000000, 0.000000, 6.643856, 6.643856, 6.643856),
#'   Condition_2 = rep("Donor1", 6),
#'   Region_1 = rep("VH", 6),
#'   PTM = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
#'   PTM_type = c(NA, "O", NA, NA, NA, NA)
#' )
#' domain <- data.frame(
#'   domain_type = c("CDR H1", "CDR H2", "CDR H3"),
#'   Region_1 = c("VH", "VH", "VH"),
#'   Region_2 = c("Arm_1", "Arm_1", "Arm_1"),
#'   Condition_1 = c("Drug1", "Drug1", "Drug1"),
#'   domain_start = c(1, 3, 5),
#'   domain_end = c(2, 4, 6),
#'   domain_color = c("#F8766D", "#B79F00", "#00BA38"),
#'   domain_fill_color = c("#F8766D", "#B79F00", "#00BA38"),
#'   domain_label_y = c(1.35, 1.35, 1.35)
#' )
#' x_axis_vars <- c("Region_2", "Region_1")
#' y_axis_vars <- c("Condition_2")
#' PTM_color <- c(
#'   "Ox" = "red",
#'   "Deamid" = "cyan",
#'   "Cam" = "blue",
#'   "Acetyl" = "magenta"
#' )
#' p <- create_peptide_plot(
#'   data,
#'   y_axis_vars,
#'   x_axis_vars,
#'   y_expand = c(0.2, 0.2),
#'   x_expand = c(0.5, 0.5),
#'   theme_options = list(),
#'   labs_options = list(title = "Area Plot", x = "Position", fill = "Area"),
#'   color_fill_column = 'Area',
#'   fill_gradient_options = list(),
#'   label_size = 5,
#'   add_domain = TRUE,
#'   domain = domain,
#'   domain_start_column = "domain_start",
#'   domain_end_column = "domain_end",
#'   domain_type_column = "domain_type",
#'   domain_border_color_column = "domain_color",
#'   domain_fill_color_column = "domain_fill_color",
#'   add_domain_label = TRUE,
#'   domain_label_size = 4,
#'   domain_label_y_column = "domain_label_y",
#'   domain_label_color = "black",
#'   PTM = FALSE,
#'   PTM_type_column = "PTM_type",
#'   PTM_color = PTM_color,
#'   add_label = TRUE,
#'   label_column = "Character",
#'   label_filter = NULL,
#'   label_y = 1,
#'   column_order = list(Region_1 = 'VH')
#' )
#' print(p)
#'
#' @import ggplot2
#' @import DT
#' @import ggforce
#' @import ggnewscale
#' @importFrom grDevices rainbow
#' @importFrom stats setNames
#' @importFrom utils modifyList
#' @importFrom ggh4x facet_nested
#' @importFrom ggplot2 element_blank element_text element_rect element_line
#' @importFrom rlang syms
#'
#' @export
#'
create_peptide_plot <- function(data,
                                y_axis_vars,
                                x_axis_vars,
                                y_expand = c(0.1, 0.15),
                                x_expand = c(0.6, 0.6),
                                theme_options = NULL,
                                labs_options = NULL,
                                color_fill_column,
                                fill_gradient_options = list(),
                                label_size = 3,
                                add_domain = TRUE,
                                domain = NULL,
                                domain_start_column =
                                  "domain_start",
                                domain_end_column = "domain_end",
                                domain_type_column
                                = "domain_type",
                                domain_border_color_column = NULL,
                                domain_fill_color_column = NULL,
                                add_domain_label = TRUE,
                                domain_label_size = 4,
                                domain_label_y_column = NULL,
                                domain_label_color = "black",
                                PTM = FALSE,
                                PTM_type_column =
                                  "PTM_type",
                                PTM_color = NULL,
                                add_label = TRUE,
                                label_column = "Character",
                                label_filter = NULL,
                                label_y = 1.28,
                                column_order = NULL) {
  # Default fill gradient options
  default_fill_gradient_options <- list(
    low = "grey80",
    high = "black",
    space = "Lab",
    na.value = "white",
    guide = "colourbar",
    aesthetics = "fill"
  )

  # Convert column names to quosures
  y_axis_vars <- syms(y_axis_vars)
  x_axis_vars <- syms(x_axis_vars)

  # Merge default and user-provided fill gradient options
  fill_gradient_options <- modifyList(default_fill_gradient_options, fill_gradient_options)

  # Create whole_labels with complete position coverage across facets
  whole_labels <- data
  domain_labels <- domain

  if (!is.null(label_filter)) {
    # Get all facet grouping variables
    facet_vars <- intersect(
      c(sapply(c(y_axis_vars, x_axis_vars), as.character)),
      names(data)
    )

    # Get position-based columns from original data
    pos_cols <- unique(intersect(
      c("Position", label_column, facet_vars),
      names(data)
    ))

    # Create combination grid for non-label_filter columns
    non_label_cols <- setdiff(pos_cols, names(label_filter))
    label_grid <- unique(data[non_label_cols])

    # Process each label_filter column, splitting commas and expanding rows
    for (col in names(label_filter)) {
      # Split comma-separated values into a vector (e.g., "D1,D2" → c("D1", "D2"))
      values <- trimws(strsplit(label_filter[[col]], ",")[[1]])

      # Expand label_grid to include each split value as a new row
      if (length(values) > 0) {
        # Cross-join current label_grid with the split values for this column
        split_df <- data.frame(new_col = values)
        names(split_df) <- col
        label_grid <- merge(label_grid, split_df, by = character())
      }
    }

    # Remove duplicate positions while keeping all columns
    whole_labels <- label_grid

    # Get domain columns that match the filter
    # Get all facet grouping variables
    domain_facet_vars <- intersect(
      c(sapply(c(y_axis_vars, x_axis_vars), as.character)),
      names(domain)
    )

    domain_pos_cols <- unique(intersect(
      c(domain_start_column, domain_end_column, domain_type_column, domain_label_y_column, domain_facet_vars),
      names(domain)
    ))

    # Create combination grid for domain labels
    domain_non_label_cols <- setdiff(domain_pos_cols, names(label_filter))
    domain_label_grid <- unique(domain[domain_non_label_cols])

    for (col in names(label_filter)) {
      values <- trimws(strsplit(label_filter[[col]], ",")[[1]])

      if (length(values) > 0) {
        split_df <- data.frame(new_col = values)
        names(split_df) <- col
        domain_label_grid <- merge(domain_label_grid, split_df, by = character())
      }
    }
    domain_labels <- domain_label_grid
  }

  # Handle NA replacement
  data[[color_fill_column]][data[[color_fill_column]] == 0] <- NA

  # Reorder columns if column_order is provided
  if (!is.null(column_order)) {
    for (col in names(column_order)) {
      order_levels <- unlist(strsplit(column_order[[col]], ","))
      data <- data[data[[col]] %in% order_levels, ]
      data[[col]] <- factor(data[[col]], levels = order_levels)
      if (col %in% colnames(whole_labels)) {
        whole_labels[[col]] <- factor(whole_labels[[col]], levels = order_levels)
      }
      if (col %in% colnames(domain)) {
        domain[[col]] <- factor(domain[[col]], levels = order_levels)
      }
    }
    # Reorder the data frame based on the specified factor levels
    data <- data[do.call(order, data[names(column_order)]), ]
  }

   p <- ggplot(data, aes(x = Position)) +
    geom_raster(data = data,
                aes(
                  x = Position,
                  y = 0.5,
                  fill = !!sym(color_fill_column)
                ),
                interpolate = FALSE) +

    do.call(scale_fill_gradient, fill_gradient_options) +
    theme_minimal()


  if (PTM) {
    PTM_data <- data[data$PTM == TRUE, ]
    unique_ptm_types <- unique(PTM_data[[PTM_type_column]])
    ptm_colors <- if (!is.null(PTM_color))
      PTM_color
    else
      setNames(rainbow(length(unique_ptm_types)), unique_ptm_types)

    p <- p +
      new_scale_fill() +
      geom_rect(
        data = PTM_data,
        inherit.aes = FALSE,
        aes(
          xmin = Position - 0.5,
          xmax = Position + 0.5,
          ymin = 1,
          ymax = 1.2,
          fill = !!sym(PTM_type_column)
        ),
        alpha = 0.5,
        show.legend = TRUE
      ) +
      scale_fill_manual(name = "PTM Types",
                        values = ptm_colors,
                        na.value = "white")
  }


  if (add_domain) {
    # Create new columns for actual color values
    domain$border_color_ <- if (!is.null(domain_border_color_column))
      domain[[domain_border_color_column]]
    else "black"
    domain$fill_color_ <- if (!is.null(domain_fill_color_column))
      domain[[domain_fill_color_column]]
    else "white"

    p <- p +
      geom_rect(
        data = domain,
        inherit.aes = FALSE,
        aes(
          xmin = !!sym(domain_start_column) - 0.5,
          xmax = !!sym(domain_end_column) + 0.5,
          ymin = -Inf,
          ymax = Inf,
          color = I(.data$border_color_),
          fill = I(.data$fill_color_)
        ),
        alpha = 0.2,
        size = 0.1,
        show.legend = TRUE  # Disable automatic legend
      )
  }

  if (add_domain_label) {
    # Calculate label positions using FILTERED domain_labels
    domain_labels$label_x <- (domain_labels[[domain_start_column]] +
                                domain_labels[[domain_end_column]]) / 2

    # Determine Y position
    if (!is.null(domain_label_y_column) && domain_label_y_column %in% names(domain_labels)) {
      domain_labels$label_y <- domain_labels[[domain_label_y_column]]
    } else {
      domain_labels$label_y <- 1.25
    }

    p <- p +
      geom_text(
        data = domain_labels,  # Use filtered domain_labels
        aes(
          x = .data$label_x,
          y = .data$label_y,
          label = !!sym(domain_type_column)
        ),
        size = domain_label_size,
        fontface = "bold",
        color = domain_label_color
      )
  }

  if (add_label) {
    p <- p +
      geom_text(
        data = whole_labels,
        aes(
          x = Position,
          y = 0.1,
          label = !!sym(label_column)
        ),
        size = label_size,
        nudge_y = label_y
      )
  }

  p <- p +
    facet_nested(
      vars(!!!y_axis_vars),
      vars(!!!x_axis_vars),
      switch = "y",
      drop = FALSE,
      scales = "free",
      space = "free"
    ) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    scale_y_continuous(expand = expansion(add = y_expand)) +
    scale_x_continuous(expand = expansion(add = x_expand)) +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      strip.text.y.left = element_text(angle = 0, size = 15),
      strip.text.x = element_text(size = 15, margin = margin(b = 10)),
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(1, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(color = "black", fill = "white", linewidth = 0.1),
      strip.background = element_blank(),
      plot.background = element_rect(color = "black", fill = "white"),
      plot.margin = margin(10, 10, 10, 10),
      legend.box = "horizontal",
      legend.key.size = unit(1, "cm"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )

  # Apply additional theme options if provided
  if (!is.null(theme_options)) {
    p <- p + do.call(theme, theme_options)
  }

  # Apply additional labs options if provided
  if (!is.null(labs_options)) {
    p <- p + do.call(labs, labs_options)
  }

  return(p)
}

