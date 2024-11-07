## ----setup--------------------------------------------------------------------
library(PepMapViz)

## -----------------------------------------------------------------------------
# Access the input files
input_file_folder <- system.file("extdata", package = "PepMapViz")

# Read the input files
resulting_df <- combine_files_from_folder(input_file_folder)
head(resulting_df)

## -----------------------------------------------------------------------------
# Strip the sequence
striped_data_peaks <- strip_sequence(resulting_df, "Peptide", "Sequence", "PEAKS")
head(striped_data_peaks)

## -----------------------------------------------------------------------------
# Extract modifications information
PTM_table <- data.frame(PTM_mass = c("15.99", ".98", "57.02", "42.01"),
                        PTM_type = c("Ox", "Deamid", "Cam", "Acetyl"))
converted_data_peaks <- obtain_mod(
  striped_data_peaks,
  "Peptide",
  "PEAKS",
  strip_seq_col = NULL,
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
region <- data.frame(
  Epitope = c("Boco", "Boco", "Boco", "Boco", "Boco", "Boco"),
  Chain = c("HC", "HC", "HC", "HC", "LC", "LC"),
  Region = c("VH", "CH1", "CH2", "CH3", "VL", "CL"),
  Region_start = c(1,119,229,338,1,108),
  Region_end = c(118,228,337,444,107,214)
)
result_with_psm <- data.frame()
for (i in 1:nrow(region)) {
  chain <- region$Chain[i]
  region_start <- region$Region_start[i]
  region_end <- region$Region_end[i]
  region_name <- region$Region[i]

  temp <- data_with_psm[data_with_psm$Chain == chain & 
                          data_with_psm$Position >= region_start & 
                          data_with_psm$Position <= region_end, ]
  temp$Region <- region_name

  result_with_psm <- rbind(result_with_psm, temp)
}
  
head(result_with_psm)

## ----fig.width=30, fig.height=6-----------------------------------------------
# Plotting peptide in whole provided sequence
domain <- data.frame(
  domain_type = c("CDR H1", "CDR H2", "CDR H3", "CDR L1", "CDR L2", "CDR L3"),
  Region = c("VH", "VH", "VH",  "VL", "VL", "VL"),
  Epitope = c("Boco", "Boco", "Boco", "Boco", "Boco", "Boco"),
  domain_start = c(26, 50, 97,  24, 50, 89),
  domain_end = c(35, 66, 107,  34, 56, 97)
)
x_axis_vars <- c("Region")
y_axis_vars <- c("Donor")
column_order <- list(
    Donor = "D1,D2,D3,D4,D5,D6,D7,D8",
    Region = "VH,CH1,CH2,CH3,VL,CL"
)
domain_color <- c(
"CDR H1" = "#F8766D",
"CDR H2" = "#B79F00",
"CDR H3" = "#00BA38",
"CDR L1" = "#00BFC4",
"CDR L2" = "#619CFF",
"CDR L3" = "#F564E3"
)
PTM_color <- c(
  "Ox" = "red",
  "Deamid" = "cyan",
  "Cam" = "blue",
  "Acetyl" = "magenta"
)
label_value = list(Donor = "D1")
p_psm <- create_peptide_plot(
  result_with_psm,
  y_axis_vars,
  x_axis_vars,
  y_expand = c(0.2, 0.2),
  x_expand = c(0.5, 0.5),
  theme_options = list(legend.box = "horizontal"),
  labs_options = list(title = "PSM Plot", x = "Position", fill = "PSM"),
  color_fill_column = 'PSM',
  fill_gradient_options = list(limits = c(0, 160)),  # Set the limits for the color scale
  label_size = 1.9,
  add_domain = TRUE,
  domain = domain,
  domain_start_column = "domain_start",
  domain_end_column = "domain_end",
  domain_type_column = "domain_type",
  domain_color = domain_color,
  PTM = TRUE,
  PTM_type_column = "PTM_type",
  PTM_color = PTM_color,
  add_label = TRUE,
  label_column = "Character",
  label_value = label_value,
  column_order = column_order
)
print(p_psm)


