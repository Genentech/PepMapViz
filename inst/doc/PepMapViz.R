## ----setup--------------------------------------------------------------------
library(PepMapViz)

## -----------------------------------------------------------------------------
# To access the input files for proteomics results, specify the file path by replacing it with your own directory path.
input_file_folder <- system.file("extdata/example_PEAKS_result", package = "PepMapViz")

# Read the input files
resulting_df <- combine_files_from_folder(input_file_folder)

# Optional. Incorporating metadata into your analysis
meta_data_path <- system.file("extdata/example_PEAKS_metadata", package = "PepMapViz")
meta_data_df <- combine_files_from_folder(meta_data_path)
resulting_df <- merge(
  x = resulting_df,
  y = meta_data_df,
  by = "Source File",
  all.x = TRUE  # Left join behavior
)
head(resulting_df)

## -----------------------------------------------------------------------------
# Strip the sequence
striped_data_peaks <- strip_sequence(resulting_df, "Peptide", "Sequence", "PEAKS")
head(striped_data_peaks)

## -----------------------------------------------------------------------------
# Extract modifications information
PTM_table <- data.frame(PTM_mass = c("15.99", "0.98", "57.02", "42.01"),
                        PTM_type = c("Ox", "Deamid", "Cam", "Acetyl"))
converted_data_peaks <- obtain_mod(
  striped_data_peaks,
  "Peptide",
  "PEAKS",
  seq_column = NULL,
  PTM_table,
  PTM_annotation = TRUE,
  PTM_mass_column = "PTM_mass"
)
head(converted_data_peaks)

## -----------------------------------------------------------------------------
# Match peptide sequence with provided sequence and calculate positions
whole_seq <- data.frame(
  Epitope = c("Boco", "Boco"),
  Chain = c("HC", "LC"),
  Region_Sequence = c("QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGEISPFGGRTNYNEKFKSRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARERPLYASDLWGQGTTVTVSSASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSNFGTQTYTCNVDHKPSNTKVDKTVERKCCVECPPCPAPPVAGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTFRVVSVLTVVHQDWLNGKEYKCKVSNKGLPSSIEKTISKTKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPMLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK", 
                      "DIQMTQSPSSLSASVGDRVTITCRASQGISSALAWYQQKPGKAPKLLIYSASYRYTGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQRYSLWRTFGQGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
  )
)
matching_result <- match_and_calculate_positions(
  converted_data_peaks,
  'Sequence',
  whole_seq,
  match_columns = NULL,
  sequence_length = c(10, 30),
  column_keep = c(
    "PTM_mass",
    "PTM_position",
    "reps",
    "Area",
    "Donor",
    "PTM_type"
  )
)
head(matching_result)

## -----------------------------------------------------------------------------
# Quantify matched peptide sequences by PSM
# Customize the matching_columns and distinct_columns variables to align with your dataset specifics.
matching_columns = c("Chain", "Epitope")
distinct_columns = c("Donor")

data_with_psm <- peptide_quantification(
  whole_seq,
  matching_result,
  matching_columns,
  distinct_columns,
  quantify_method = "PSM",
  with_PTM = TRUE,
  reps = TRUE
)
head(data_with_psm)

## ----fig.width=30, fig.height=6-----------------------------------------------
# Plotting peptide in whole provided sequence
domain <- data.frame(
  domain_type = c("VH", "CH1", "CH2", "CH3", "VL", "CL", "CDR H1", "CDR H2", "CDR H3", "CDR L1", "CDR L2", "CDR L3"),
  Chain = c("HC", "HC", "HC", "HC",  "LC", "LC", "HC", "HC", "HC",  "LC", "LC", "LC"),
  Epitope = c("Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco", "Boco"),
  domain_start = c(1, 119, 229, 338, 1, 108, 26, 50, 97, 24, 50, 89),
  domain_end = c(118, 228, 337, 444, 107, 214, 35, 66, 107,  34, 56, 97),
  domain_color = c("black", "black", "black", "black", "black", "black", "#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"),
  domain_fill_color = c("white", "white", "white", "white", "white", "white", "yellow", "yellow", "yellow", "yellow", "yellow", "yellow"), 
  domain_label_y = c(1.7, 1.7, 1.7, 1.7, 1.7, 1.7, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4)
)
x_axis_vars <- c("Chain")
y_axis_vars <- c("Donor")
column_order <- list(
    Donor = "D1,D2,D3,D4,D5,D6,D7,D8",
    Chain = "HC,LC"
)
PTM_color <- c(
  "Ox" = "red",
  "Deamid" = "cyan",
  "Cam" = "blue",
  "Acetyl" = "magenta"
)
label_filter = list(Donor = "D1")
p_psm <- create_peptide_plot(
  data_with_psm,
  y_axis_vars,
  x_axis_vars,
  y_expand = c(0.2, 0.2),
  x_expand = c(0.5, 0.5),
  theme_options = list(legend.box = "horizontal", legend.position = "bottom"),
  labs_options = list(title = "PSM Plot", x = "Position", fill = "PSM"),
  color_fill_column = 'PSM',
  fill_gradient_options = list(),  # Set the limits for the color scale
  label_size = 1.3,
  add_domain = TRUE,
  domain = domain,
  domain_start_column = "domain_start",
  domain_end_column = "domain_end",
  domain_type_column = "domain_type",
  domain_border_color_column = "domain_color",
  domain_fill_color_column = "domain_fill_color",
  add_domain_label = TRUE,
  domain_label_size = 2,
  domain_label_y_column = "domain_label_y",
  domain_label_color = "black",
  PTM = TRUE,
  PTM_type_column = "PTM_type",
  PTM_color = PTM_color,
  add_label = TRUE,
  label_column = "Character",
  label_filter = label_filter,
  label_y = 1,
  column_order = column_order
)
print(p_psm)


