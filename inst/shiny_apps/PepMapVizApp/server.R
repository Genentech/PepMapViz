# inst/shiny-apps/PepMapVizApp/server.R
library(DT)
library(rlang)
library(data.table)
library(mzID)
library(MSnbase)
library(ggplot2)

# Define default PTM table
PTM_table <- data.frame(
  PTM_type = c("Ox", "Deamid","Deamid", "Cam", "Acetyl"),
  PTM_mass = c("15.99", ".98", "0.98", "57.02", "42.01")
)
seq_table <- data.frame(
  Epitope = c("Boco", "Boco"),
  Chain = c("HC", "LC"),
  Region_Sequence = c("QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGEISPFGGRTNYNEKFKSRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARERPLYASDLWGQGTTVTVSSASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSNFGTQTYTCNVDHKPSNTKVDKTVERKCCVECPPCPAPPVAGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTFRVVSVLTVVHQDWLNGKEYKCKVSNKGLPSSIEKTISKTKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPMLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK",
                      "DIQMTQSPSSLSASVGDRVTITCRASQGISSALAWYQQKPGKAPKLLIYSASYRYTGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQRYSLWRTFGQGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC")
)
region_table <- data.frame(
    Epitope = c("Boco", "Boco", "Boco", "Boco", "Boco", "Boco"),
    Chain = c("HC", "HC", "HC", "HC", "LC", "LC"),
    Region = c("VH", "CH1", "CH2", "CH3", "VL", "CL"),
    Region_start = c(1, 119, 229, 338, 1, 108),
    Region_end = c(118, 228, 337, 444, 107, 214)
)
domain_table <- data.frame(
  domain_type = c("VH", "CH1", "CH2", "CH3", "VL", "CL", "CDR H1", "CDR H2", "CDR H3", "CDR L1", "CDR L2", "CDR L3"),
  Chain = c("HC", "HC", "HC", "HC",  "LC", "LC", "HC", "HC", "HC",  "LC", "LC", "LC"),
  Epitope = c("Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco"),
  domain_start = c(1, 119, 229, 338, 1, 108, 26, 50, 97, 24, 50, 89),
  domain_end = c(118, 228, 337, 444, 107, 214, 35, 66, 107,  34, 56, 97),
  domain_color = c("black", "black", "black", "black", "black", "black", "#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"),
  domain_fill_color = c("white", "white", "white", "white", "white", "white", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow"),
  domain_label_y = c(1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4)
)
PTM_color <- data.frame(
  PTM_type = c("Ox", "Deamid", "Cam", "Acetyl"),
  color = c("red", "cyan", "blue", "magenta")
)


function(input, output, session) {
  # Reactive values to store data
  raw_data <- reactiveVal(NULL)
  processed_data <- reactiveVal(NULL)
  strip_data <- reactiveVal(NULL)
  mod_data <- reactiveVal(NULL)

  # Add reactive for metadata
  meta_data <- reactiveVal(NULL)

  # Observe metadata upload
  observeEvent(input$metadata, {
    req(input$metadata)

    tryCatch({
      meta <- data.table::fread(input$metadata$datapath)
      meta_data(meta)
      showNotification("Metadata loaded!", type = "message")
    }, error = function(e) {
      showNotification(paste("Metadata Error:", e$message), type = "error")
    })
  })

  # Read and combine files when "Read File(s)" is clicked
  observeEvent(input$readfile, {
    req(input$inputFolder)

    tryCatch({
      uploaded_files <- input$inputFolder

      withProgress(message = 'Processing files', value = 0, {
        # Initialize progress parameters
        n_files <- nrow(uploaded_files)
        progress_step <- 1/n_files

        # Read and combine files with progress updates
        file_data_list <- lapply(1:n_files, function(i) {
          incProgress(progress_step,
                      detail = paste("Reading", uploaded_files$name[i],
                                     "(", i, "/", n_files, ")"))

          file_path <- uploaded_files$datapath[i]
          file_name <- uploaded_files$name[i]
          ext <- tolower(tools::file_ext(file_name))

          # Read file with format handling
          df <- if (ext %in% c("csv", "tsv", "txt")) {
            data.table::fread(file_path)
          } else if (ext == "mzid") {
            if (!requireNamespace("mzID", quietly = TRUE)) {
              stop("Install 'mzID' with: BiocManager::install('mzID')")
            }
            mzID::flatten(mzID::mzID(file_path))
          } else if (ext == "mztab") {
            if (!requireNamespace("MSnbase", quietly = TRUE)) {
              stop("Install 'MSnbase' with: BiocManager::install('MSnbase')")
            }
            MSnbase::fData(MSnbase::readMzTabData(file_path, "PSM"))
          } else {
            stop("Unsupported file type: ", ext)
          }

          return(df)
        })

        # Combine all files
        incProgress(progress_step, detail = "Combining files")
        combined_df <- data.table::rbindlist(file_data_list, fill = TRUE)

        # Merge metadata if exists
        if (!is.null(meta_data())) {
          incProgress(progress_step, detail = "Merging metadata")
          req(input$merge_column)
          merge_cols <- trimws(unlist(strsplit(input$merge_column, ",")))

          validate(
            need(all(merge_cols %in% names(combined_df)),
                 "Some match columns not found in main data"),
            need(all(merge_cols %in% names(meta_data())),
                 "Some match columns not found in metadata")
          )

          combined_df <- merge(
            x = combined_df,
            y = meta_data(),
            by = merge_cols,
            all.x = TRUE  # Left join behavior
          )
        }

        # Final processing
        incProgress(progress_step, detail = "Finalizing")

        # Store results
        raw_data(combined_df)
        processed_data(combined_df)
      })

      showNotification("Files successfully loaded!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })

  # Strip sequences when "Strip Sequences" is clicked
  observeEvent(input$strip, {
    req(processed_data(), input$datatype, input$strip_column, input$result_column)

    tryCatch({
      # Use withProgress for visual feedback
      withProgress(message = 'Stripping sequences...', value = 0, {
        # Get current data
        current_data <- processed_data()

        # Perform sequence stripping
        stripped_df <- PepMapViz::strip_sequence(
          current_data,
          input$strip_column,
          input$result_column,
          input$datatype
        )

        # Update processed data
        processed_data(stripped_df)
        strip_data(stripped_df)

        incProgress(1, detail = "Done!")
      })

    }, error = function(e) {
      showNotification(paste("Stripping Error:", e$message), type = "error")
    })
  })

  # Create reactive value for editable PTM table
  ptm_data <- reactiveVal(PTM_table)

  observeEvent(input$addRow, {
    req(ptm_data())

    tryCatch({
      # Create new row with same structure as current data
      new_row <- as.list(rep(NA, ncol(ptm_data())))
      names(new_row) <- names(ptm_data())

      # Convert to data.table and bind
      new_row <- data.table::as.data.table(new_row)
      ptm_data(data.table::rbindlist(list(ptm_data(), new_row), fill = TRUE))

    }, error = function(e) {
      showNotification(paste("Error adding row:", e$message), type = "error")
    })
  })

  # Enhanced removeRow function
  observeEvent(input$removeRow, {
    req(ptm_data())

    tryCatch({
      current_data <- ptm_data()

      if (nrow(current_data) == 0) {
        showNotification("No rows to remove", type = "warning")
        return()
      }

      if (!is.null(input$ptmTable_rows_selected)) {
        # Convert to numeric indices if they aren't already
        rows_to_remove <- as.numeric(input$ptmTable_rows_selected)

        # Validate row indices
        if (any(rows_to_remove > nrow(current_data)) || any(rows_to_remove < 1)) {
          showNotification("Invalid row selection", type = "error")
          return()
        }

        ptm_data(current_data[-rows_to_remove, ])
      } else {
        showNotification("No rows selected. Removing last row.", type = "warning")
        ptm_data(current_data[-nrow(current_data), ])
      }

      # Reset row selection
      DT::dataTableProxy("ptmTable") %>% selectRows(NULL)

    }, error = function(e) {
      showNotification(paste("Error removing row:", e$message), type = "error")
    })
  })

  # Handle table edits
  proxy_ptmTable <- DT::dataTableProxy("ptmTable")
  observeEvent(input$ptmTable_cell_edit, {
    info <- input$ptmTable_cell_edit
    info$col <- info$col + 1
    updated_data <- editData(ptm_data(), info)
    ptm_data(updated_data)
    replaceData(proxy_ptmTable, updated_data, resetPaging = FALSE)
  })

  observeEvent(input$extract, {
    req(processed_data(), input$mod_column)

    tryCatch({
      withProgress(message = 'Extracting modifications...', value = 0.5, {

        # Get and validate PTM table
        current_ptm_table <- ptm_data()

        # Convert to proper format
        if (data.table::is.data.table(current_ptm_table)) {
          current_ptm_table <- as.data.frame(current_ptm_table)
        }

        # Validate columns exist in data
        req_cols <- c(input$mod_column)
        if (input$mod_datatype %in% c("MSFragger","mzIdenML","mzTab")) {
          req_cols <- c(req_cols, input$seq_column)
        }

        if (!all(req_cols %in% names(processed_data()))) {
          stop("Required columns not found in data: ",
               paste(setdiff(req_cols, names(processed_data())), collapse = ", "))
        }

        # Prepare parameters
        params <- list(
          data = processed_data(),
          mod_column = input$mod_column,
          type = input$datatype,
          PTM_table = current_ptm_table,
          PTM_annotation = input$PTM_annotation,
          PTM_mass_column = input$PTM_mass_column
        )

        # Add MSFragger column if needed
        if (input$datatype == "MSFragger" && !is.null(input$seq_column)) {
          params$seq_column <- input$seq_column
        }

        # Call modification function
        extracted_data <- do.call(PepMapViz::obtain_mod, params)

        # Update reactive values
        processed_data(extracted_data)
        mod_data(extracted_data)

        incProgress(1, detail = "Done!")
      })

      showNotification("Modifications successfully extracted!", type = "message")

    }, error = function(e) {
      cat("\nERROR:", e$message, "\n")
      traceback()
      showNotification(paste("Extraction Error:", e$message), type = "error")
    })
  })

  # Create reactive value for editable sequence table
  seq_data <- reactiveVal(seq_table)

  observeEvent(input$addseqRow, {
    req(seq_data())

    tryCatch({
      # Create new row with same structure as current data
      new_row <- as.list(rep(NA, ncol(seq_data())))
      names(new_row) <- names(seq_data())

      # Convert to data.table and bind
      new_row <- data.table::as.data.table(new_row)
      seq_data(data.table::rbindlist(list(seq_data(), new_row), fill = TRUE))

    }, error = function(e) {
      showNotification(paste("Error adding row:", e$message), type = "error")
    })
  })

  # Enhanced removeRow function
  observeEvent(input$removeseqRow, {
    req(seq_data())

    tryCatch({
      current_data <- seq_data()

      if (nrow(current_data) == 0) {
        showNotification("No rows to remove", type = "warning")
        return()
      }

      if (!is.null(input$seqTable_rows_selected)) {
        # Convert to numeric indices if they aren't already
        rows_to_remove <- as.numeric(input$seqTable_rows_selected)

        # Validate row indices
        if (any(rows_to_remove > nrow(current_data)) || any(rows_to_remove < 1)) {
          showNotification("Invalid row selection", type = "error")
          return()
        }

        seq_data(current_data[-rows_to_remove, ])
      } else {
        showNotification("No rows selected. Removing last row.", type = "warning")
        seq_data(current_data[-nrow(current_data), ])
      }

      # Reset row selection
      DT::dataTableProxy("SeqTable") %>% selectRows(NULL)

    }, error = function(e) {
      showNotification(paste("Error removing row:", e$message), type = "error")
    })
  })

  # Handle table edits
  proxy_SeqTable <- DT::dataTableProxy("SeqTable")
  observeEvent(input$SeqTable_cell_edit, {
    info <- input$SeqTable_cell_edit
    info$col <- info$col + 1
    new_data <- DT::editData(seq_data(), info)
    seq_data(new_data)
    replaceData(proxy_SeqTable, new_data, resetPaging = FALSE)
  })

  observeEvent(input$uploadseqtable, {
    req(input$uploadseqtable)

    tryCatch({
      # Read uploaded file
      file <- input$uploadseqtable
      ext <- tools::file_ext(file$datapath)

      # Validate file type
      validate(
        need(ext %in% c("csv", "tsv", "txt"),
             "Unsupported file format. Please upload CSV, TSV or TXT file.")
      )

      # Read data based on file type
      uploaded_data <- data.table::fread(file$datapath)

      # Update sequence data
      seq_data(uploaded_data)

      showNotification("Sequence table uploaded successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Upload error:", e$message), type = "error")
    })
  })

  # Reactive value to store matched data
  match_data <- reactiveVal(NULL)

  observeEvent(input$match, {
    req(processed_data(), seq_data(), input$match_column)

    tryCatch({
      withProgress(message = 'Matching sequences...', value = 0.5, {
        # Process multiple column inputs
        match_cols <- trimws(unlist(strsplit(input$match_params, ",")))
        keep_cols <- trimws(unlist(strsplit(input$column_keep, ",")))

        # Validate inputs
        validate(
          need(all(match_cols %in% names(processed_data())),
               "Some match columns not found in main data"),
          need(all(match_cols %in% names(seq_data())),
               "Some match columns not found in sequence table"),
          need(nrow(seq_data()) > 0, "Sequence table is empty")
        )

        # Perform matching
        matched <- PepMapViz::match_and_calculate_positions(
          peptide_data = processed_data(),
          column = input$match_column,
          whole_seq = seq_data(),
          match_columns = match_cols,
          sequence_length = input$length,
          column_keep = keep_cols
        )

        # Store results
        match_data(matched)
        processed_data(matched)

        incProgress(1, detail = "Done!")
      })

      showNotification("Matching completed!", type = "message")

    }, error = function(e) {
      showNotification(paste("Matching Error:", e$message), type = "error")
    })
  })

  # Reactive value to store quantification data
  quant_data <- reactiveVal(NULL)

  observeEvent(input$quantify, {
    req(match_data(), seq_data(), input$quantify_method)

    tryCatch({
      withProgress(message = 'Performing quantification...', value = 0.3, {
        # Process multiple column inputs
        match_quantify_cols <- trimws(unlist(strsplit(input$match_quantify_column, ",")))
        distinct_cols <- trimws(unlist(strsplit(input$distinct_columns, ",")))

        # Validate inputs
        validate(
          need(input$quantify_method %in% c("PSM", "Area"),
               "Invalid quantification method"),
          need(all(match_quantify_cols %in% names(match_data())),
               "Some match columns not found in main data"),
          need(all(distinct_cols %in% names(match_data())),
               "Some distince columns not found in main table")
           )

        # Perform matching
        quant_results <- PepMapViz::peptide_quantification(
          whole_seq = seq_data(),
          matching_result = match_data(),
          matching_columns = match_quantify_cols,
          distinct_columns = distinct_cols,
          quantify_method = input$quantify_method,
          area_column = input$area_column,
          with_PTM = input$with_PTM,
          reps = input$reps
        )

        # Store results
        quant_data(quant_results)
        processed_data(quant_results)  # Update main data if needed

        incProgress(1, detail = "Done!")
      })

      showNotification("Quantification completed successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Matching Error:", e$message), type = "error")
    })
  })

  # Create reactive value for editable PTM table
  region_data <- reactiveVal(region_table)

  observeEvent(input$addRegionRow, {
    req(region_data())

    tryCatch({
      # Create new row with same structure as current data
      new_row <- as.list(rep(NA, ncol(region_data())))
      names(new_row) <- names(region_data())

      # Convert to data.table and bind
      new_row <- data.table::as.data.table(new_row)
      region_data(data.table::rbindlist(list(region_data(), new_row), fill = TRUE))

    }, error = function(e) {
      showNotification(paste("Error adding row:", e$message), type = "error")
    })
  })

  # Enhanced removeRow function
  observeEvent(input$removeRegionRow, {
    req(region_data())

    tryCatch({
      current_data <- region_data()

      if (nrow(current_data) == 0) {
        showNotification("No rows to remove", type = "warning")
        return()
      }

      if (!is.null(input$regionTable_rows_selected)) {
        # Convert to numeric indices if they aren't already
        rows_to_remove <- as.numeric(input$regionTable_rows_selected)

        # Validate row indices
        if (any(rows_to_remove > nrow(current_data)) || any(rows_to_remove < 1)) {
          showNotification("Invalid row selection", type = "error")
          return()
        }

        region_data(current_data[-rows_to_remove, ])
      } else {
        showNotification("No rows selected. Removing last row.", type = "warning")
        region_data(current_data[-nrow(current_data), ])
      }

      # Reset row selection
      DT::dataTableProxy("regionTable") %>% selectRows(NULL)

    }, error = function(e) {
      showNotification(paste("Error removing row:", e$message), type = "error")
    })
  })

  # Handle table edits
  proxy_regionTable <- DT::dataTableProxy("regionTable")
  observeEvent(input$regionTable_cell_edit, {
    info <- input$regionTable_cell_edit
    info$col <- info$col + 1
    updated_data <- editData(region_data(), info)
    region_data(updated_data)
    replaceData(proxy_regionTable, updated_data, resetPaging = FALSE)
  })

  observeEvent(input$region_file, {
    req(input$region_file)

    tryCatch({
      # Read uploaded file
      file <- input$region_file
      ext <- tools::file_ext(file$datapath)

      # Validate file type
      validate(
        need(ext %in% c("csv", "tsv", "txt"),
             "Unsupported file format. Please upload CSV, TSV or TXT file.")
      )

      # Read data based on file type
      uploaded_data <- data.table::fread(file$datapath)

      # Update sequence data
      region_data(uploaded_data)

      showNotification("Region table uploaded successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Upload error:", e$message), type = "error")
    })
  })

  # Reactive value to store matched data
  merge_region_data <- reactiveVal(NULL)

  observeEvent(input$merge_regions, {
    req(processed_data(), region_data())

    tryCatch({
      withProgress(message = "Merging region information...", value = 0.3, {
        # Get current data
        current_data <- processed_data()
        region_info <- region_data()
        match_cols <- trimws(unlist(strsplit(input$match_region_cols, ",")))
        start_col <- input$region_start_col
        end_col <- input$region_end_col
        pos_col <- input$pos_col
        region_name_col <- input$region_col  # From previous position column input

        # Validate required columns exist
         validate(
          need(all(match_cols %in% names(processed_data())),
               "Some match columns not found in main data"),
          need(all(match_cols %in% names(region_data())),
               "Some match columns not found in region data"),
          need(start_col %in% names(region_data()),
               paste("Start column", start_col, "not found in region data")),
          need(end_col %in% names(region_data()),
               paste("End column", end_col, "not found in region data")),
          need(region_name_col %in% names(region_info),
               paste("Region name column", region_name_col, "not found in region data")),
          need(pos_col %in% names(processed_data()),
               paste("End column", pos_col, "not found in main data"))
           )

        # Initialize empty result
        merged_data <- data.frame()

        # Loop through region definitions
        for (i in 1:nrow(region_info)) {
          incProgress(amount = 0.7/nrow(region_info),
                      detail = paste("Processing region", i, "of", nrow(region_info)))

          # Get region parameters
          region_match_vals <- region_info[i, match_cols, drop = FALSE]
          region_start <- region_info[i, start_col]
          region_end <- region_info[i, end_col]
          region_name <- region_info[i, region_name_col]  # Correct column reference

          # Filter and annotate data
          temp <- current_data %>%
            dplyr::inner_join(region_match_vals, by = match_cols) %>%
            dplyr::filter(!!sym(pos_col) >= region_start,
                          !!sym(pos_col) <= region_end) %>%
            dplyr::mutate(Region = region_name)

          merged_data <- dplyr::bind_rows(merged_data, temp)
        }

        # Remove duplicates from overlapping regions
        merged_data <- merged_data %>%
          dplyr::distinct(.keep_all = TRUE)

        # Update processed data with merged results
        processed_data(merged_data)
        merge_region_data(merged_data)

        incProgress(1, detail = "Done!")
      })

      showNotification("Region merge completed!", type = "message")

    }, error = function(e) {
      showNotification(paste("Region merge failed:", e$message), type = "error")
    })
  })

  # Create reactive value for editable sequence table
  domain_data <- reactiveVal(domain_table)

  observeEvent(input$addDomainRow, {
    req(domain_data())

    tryCatch({
      # Create new row with same structure as current data
      new_row <- as.list(rep(NA, ncol(domain_data())))
      names(new_row) <- names(domain_data())

      # Convert to data.table and bind
      new_row <- data.table::as.data.table(new_row)
      domain_data(data.table::rbindlist(list(domain_data(), new_row), fill = TRUE))

    }, error = function(e) {
      showNotification(paste("Error adding row:", e$message), type = "error")
    })
  })

  # Enhanced removeRow function
  observeEvent(input$removeDomainRow, {
    req(domain_data())

    tryCatch({
      current_data <- domain_data()

      if (nrow(current_data) == 0) {
        showNotification("No rows to remove", type = "warning")
        return()
      }

      if (!is.null(input$domainTable_rows_selected)) {
        # Convert to numeric indices if they aren't already
        rows_to_remove <- as.numeric(input$domainTable_rows_selected)

        # Validate row indices
        if (any(rows_to_remove > nrow(current_data)) || any(rows_to_remove < 1)) {
          showNotification("Invalid row selection", type = "error")
          return()
        }

        domain_data(current_data[-rows_to_remove, ])
      } else {
        showNotification("No rows selected. Removing last row.", type = "warning")
        domain_data(current_data[-nrow(current_data), ])
      }

      # Reset row selection
      DT::dataTableProxy("domainTable") %>% selectRows(NULL)

    }, error = function(e) {
      showNotification(paste("Error removing row:", e$message), type = "error")
    })
  })

  # Handle table edits
  proxy_domainTable <- DT::dataTableProxy("domainTable")
  observeEvent(input$domainTable_cell_edit, {
    info <- input$domainTable_cell_edit
    info$col <- info$col + 1
    updated_data <- editData(domain_data(), info)
    domain_data(updated_data)
    replaceData(proxy_domainTable, updated_data, resetPaging = FALSE)
  })

  observeEvent(input$upload_domain, {
    req(input$upload_domain)

    tryCatch({
      # Read uploaded file
      file <- input$upload_domain
      ext <- tools::file_ext(file$datapath)

      # Validate file type
      validate(
        need(ext %in% c("csv", "tsv", "txt"),
             "Unsupported file format. Please upload CSV, TSV or TXT file.")
      )

      # Read data based on file type
      uploaded_data <- data.table::fread(file$datapath)

      # Update sequence data
      domain_data(uploaded_data)

      showNotification("Domain table uploaded successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Upload error:", e$message), type = "error")
    })
  })

  # Create reactive value for editable sequence table
  PTM_color_data <- reactiveVal(PTM_color)

  observeEvent(input$addPTMColorRow, {
    req(PTM_color_data())

    tryCatch({
      # Create new row with same structure as current data
      new_row <- as.list(rep(NA, ncol(PTM_color_data())))
      names(new_row) <- names(PTM_color_data())

      # Convert to data.table and bind
      new_row <- data.table::as.data.table(new_row)
      PTM_color_data(data.table::rbindlist(list(PTM_color_data(), new_row), fill = TRUE))

    }, error = function(e) {
      showNotification(paste("Error adding row:", e$message), type = "error")
    })
  })

  # Enhanced removeRow function
  observeEvent(input$removePTMColorRow, {
    req(PTM_color_data())

    tryCatch({
      current_data <- PTM_color_data()

      if (nrow(current_data) == 0) {
        showNotification("No rows to remove", type = "warning")
        return()
      }

      if (!is.null(input$PTMcolorTable_rows_selected)) {
        # Convert to numeric indices if they aren't already
        rows_to_remove <- as.numeric(input$PTMcolorTable_rows_selected)

        # Validate row indices
        if (any(rows_to_remove > nrow(current_data)) || any(rows_to_remove < 1)) {
          showNotification("Invalid row selection", type = "error")
          return()
        }

        PTM_color_data(current_data[-rows_to_remove, ])
      } else {
        showNotification("No rows selected. Removing last row.", type = "warning")
        PTM_color_data(current_data[-nrow(current_data), ])
      }

      # Reset row selection
      DT::dataTableProxy("PTM_color") %>% selectRows(NULL)

    }, error = function(e) {
      showNotification(paste("Error removing row:", e$message), type = "error")
    })
  })

  # Handle table edits
  proxy_PTMcolorTable <- DT::dataTableProxy("PTMcolorTable")
  observeEvent(input$PTMcolorTable_cell_edit, {
    info <- input$PTMcolorTable_cell_edit
    info$col <- info$col + 1
    updated_data <- editData(PTM_color_data(), info)
    PTM_color_data(updated_data)
    replaceData(proxy_PTMcolorTable, updated_data, resetPaging = FALSE)
  })

  # Reactive value to store column orders
  column_orders <- reactiveValues(orders = list())

  # Update column choices based on data
  observeEvent(processed_data(), {
    req(processed_data())
    updateSelectizeInput(session, "order_columns",
                         choices = names(processed_data()),
                         selected = names(column_orders$orders))
  })

  # Generate dynamic UI inputs for ordering
  output$column_order_inputs <- renderUI({
    req(input$order_columns)
    lapply(input$order_columns, function(col) {
      textInput(
        inputId = paste0("order_", col),
        label = paste("Order for", col),
        value = ifelse(!is.null(column_orders$orders[[col]]),
                       paste(column_orders$orders[[col]], collapse = ","), ""),
        placeholder = "Comma-separated values"
      )
    })
  })

  # Update column orders reactively
  observe({
    req(input$order_columns)
    new_orders <- list()

    for(col in input$order_columns) {
      input_id <- paste0("order_", col)
      if(!is.null(input[[input_id]]) && input[[input_id]] != "") {
        values <- trimws(unlist(strsplit(input[[input_id]], ",")))
        new_orders[[col]] <- values
      }
    }
    column_orders$orders <- new_orders
  })

  # Plot generation
  my_plot <- reactiveVal(NULL)

  observeEvent(input$generate_plot, {
    req(processed_data())

    theme_args <- list(
      legend.position = input$legend_position,
      legend.box = input$legend_box,
      legend.title = element_text(size = input$legend_title_size),
      legend.text = element_text(size = input$legend_text_size),
      axis.text = element_text(size = input$axis_text_size),
      axis.title = element_text(size = input$axis_title_size)
    )

    # Merge with custom theme options
    if(nzchar(input$custom_theme)) {
        custom_theme_list <- eval(parse(text = input$custom_theme))
        theme_args <- modifyList(theme_args, custom_theme_list)
    }

    # Process labs options
    labs_args <- list(
      title = input$plot_title
    )

    if(nzchar(input$custom_labs)) {
        labs_args <- modifyList(labs_args, custom_labs)
    }

    tryCatch({
      # Prepare parameters
      plot_params <- list(
        data = processed_data(),
        y_axis_vars = trimws(unlist(strsplit(input$y_axis, ","))),
        x_axis_vars = trimws(unlist(strsplit(input$x_axis, ","))),
        y_expand = c(input$y_expand_low, input$y_expand_high),
        x_expand = c(input$x_expand_low, input$x_expand_high),
        theme_options = theme_args,
        labs_options = labs_args,
        color_fill_column = input$color_fill_column,
        fill_gradient_options = if (input$fill_max == 0) {
          list(low = input$low_color, high = input$high_color)
        } else {
          list(limits = c(input$fill_min, input$fill_max),
               low = input$low_color, high = input$high_color)
        },
        label_size = input$label_size,
        add_domain = input$add_domain,
        domain = domain_data(),
        domain_start_column = input$domain_start_col,
        domain_end_column = input$domain_end_col,
        domain_type_column = input$domain_type_col,
        domain_border_color_column = input$domain_border_color_column,
        domain_fill_color_column = input$domain_fill_color_column,
        add_domain_label = input$add_domain_label,
        domain_label_size = input$domain_label_size,
        domain_label_y_column = input$domain_label_y_column,
        domain_label_color = input$domain_label_color,
        PTM = input$add_PTM,
        PTM_type_column = input$PTM_type_column,
        PTM_color = setNames(PTM_color_data()$color, PTM_color_data()$PTM_type),
        add_label = input$add_labels,
        label_column = input$label_column,
        label_filter = if (is.null(input$label_filter) || nchar(input$label_filter) == 0) {
          NULL
        } else {
          setNames(strsplit(input$label_filter, "=")[[1]][2],
                   strsplit(input$label_filter, "=")[[1]][1])
        },
        label_y = input$label_y,
        column_order = if (length(column_orders$orders) == 0) NULL else column_orders$orders
      )

      # Generate plot
      p <- do.call(PepMapViz::create_peptide_plot, plot_params)
      my_plot(p)

      # Render plot
      output$peptidePlot <- renderPlot(p)

    }, error = function(e) {
      showNotification(paste("Plot Error:", e$message), type = "error")
    })
  })


  # Render editable PTM table
  output$ptmTable <- renderDT({
    DT::datatable(
      ptm_data(),
      editable = "cell",
      rownames = FALSE,
      options = list(
        dom = 't',
        scrollX = TRUE
      )
    )
  })

  # Render editable PTM table
  output$SeqTable <- renderDT({
    DT::datatable(
      seq_data(),
      editable = "cell",
      rownames = FALSE,
      options = list(
        dom = 't',
        scrollX = TRUE
      )
    )
  })

  # Render editable PTM table
  output$regionTable <- renderDT({
    DT::datatable(
      region_data(),
      editable = "cell",
      rownames = FALSE,
      options = list(
        dom = 't',
        scrollX = TRUE
      )
    )
  })

  # Render editable PTM table
  output$domainTable <- renderDT({
    DT::datatable(
      domain_data(),
      editable = "cell",
      rownames = FALSE,
      options = list(
        dom = 't',
        scrollX = TRUE
      )
    )
  })

  # Render editable PTM table
  output$PTMcolorTable <- renderDT({
    DT::datatable(
      PTM_color_data(),
      editable = "cell",
      rownames = FALSE,
      options = list(
        dom = 't',
        scrollX = TRUE
      )
    )
  })

  # Display data preview
  output$process_dataPreview <- renderDT({
    req(raw_data())
    datatable(raw_data(),
              options = list(
                scrollX = TRUE,
                pageLength = 5,
                autoWidth = TRUE
              ))
  })
  # Display stripped data
  output$strip_dataPreview <- renderDT({
    req(strip_data())
    datatable(strip_data(),
              options = list(
                scrollX = TRUE,
                pageLength = 5,
                autoWidth = TRUE
              ))
  })
  # Add download handler
  output$downloadStrip <- downloadHandler(
    filename = function() {
      paste("stripped_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(strip_data())
      write.csv(strip_data(), file, row.names = FALSE)
    }
  )
  # Display modification data
  output$mod_dataPreview <- renderDT({
    req(mod_data())
    datatable(mod_data(),
              options = list(
                scrollX = TRUE,
                pageLength = 5,
                autoWidth = TRUE
              ))
  })
  # Add download handler
  output$downloadMod <- downloadHandler(
    filename = function() {
      paste("mod_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(mod_data())
      write.csv(mod_data(), file, row.names = FALSE)
    }
  )
  # Render matched results
  output$match_dataPreview <- DT::renderDataTable({
    req(match_data())
    DT::datatable(
      match_data(),
      options = list(
        scrollX = TRUE,
        pageLength = 10
      )
    )
  })
  # Add download handler
  output$downloadMatch <- downloadHandler(
    filename = function() {
      paste("match_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(match_data())
      write.csv(match_data(), file, row.names = FALSE)
    }
  )
  # Display modification data
  output$quantify_dataPreview <- renderDT({
    req(quant_data())
    datatable(quant_data(),
              options = list(
                scrollX = TRUE,
                pageLength = 5,
                autoWidth = TRUE
              ))
  })
  # Add download handler
  output$downloadQuant <- downloadHandler(
    filename = function() {
      paste("quant_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(quant_data())
      write.csv(quant_data(), file, row.names = FALSE)
    }
  )
  # Display modification data
  output$region_merge_dataPreview <- renderDT({
    req(merge_region_data())
    datatable(merge_region_data(),
              options = list(
                scrollX = TRUE,
                pageLength = 5,
                autoWidth = TRUE
              ))
  })
  # Add download handler
  output$downloadRegion <- downloadHandler(
    filename = function() {
      paste("Merge_region_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(merge_region_data())
      write.csv(merge_region_data(), file, row.names = FALSE)
    }
  )
  # Dynamic download handler
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("peptide-plot-", Sys.Date(), ".", input$plot_format, sep = "")
    },
    content = function(file) {
      req(my_plot(), input$plot_width, input$plot_height)

      switch(input$plot_format,
             "pdf" = pdf(file, width = input$plot_width, height = input$plot_height),
             "png" = png(file, width = input$plot_width, height = input$plot_height,
                         units = "in", res = 300)
      )

      print(my_plot())
      dev.off()
    }
  )
}
