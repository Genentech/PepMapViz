# Increase to 30MB (adjust the value as needed)
options(shiny.maxRequestSize = 1000 * 1024^2)  # 30MB in bytes

# inst/shiny-apps/PepMapVizApp/ui.R
fluidPage(
  titlePanel("PepMapViz - Interactive Peptide Visualization"),

  sidebarLayout(
    sidebarPanel(
      tags$h4(
        "Step1. Select Input Proteomics Files ▼",
        onclick = "toggleVisibility('step1_section', this)",
        style = "cursor: pointer; color: #0E4C92;"
      ),
      div(
        id = "step1_section",
        fileInput("inputFolder", "Select Input Files (Supported formats: CSV, TSV, TXT, mzID, mzTab)", multiple = TRUE,
                  accept = c(".csv", ".tsv", ".txt", ".mzid", ".mztab")),
        fileInput("metadata", "Upload Metadata File to merge (Optional)",
                  multiple = FALSE,
                  accept = c(".csv", ".tsv", ".txt")),
        textInput("merge_column", "Merge Column(s) - comma separated (Optional)",
                  placeholder = "e.g., Source File", value = "Source File"),
        actionButton("readfile", "Read File(s)", class = "btn-primary",
                     style = "margin-bottom: 20px;")
      ),

      tags$h4(
        "Step2. Strip the sequence (Optional) ▼",
        onclick = "toggleVisibility('step2_section', this)",
        style = "cursor: pointer; color: #0E4C92;"
      ),
      div(
        id = "step2_section",
        selectInput("datatype", "Select data type", choices = c("PEAKS","Spectronaut",
                 "MSFragger","Comet","DIANN","Skyline","Maxquant"), multiple = FALSE),
        textInput("strip_column", "Column name for stripping sequences",
                  placeholder = "e.g., Peptide", value = "Peptide"),
        textInput("result_column", "Column name for stripped sequences",
                  placeholder = "e.g., Sequence", value = "Sequence"),
        actionButton("strip", "Strip Sequences", class = "btn-primary",
                     style = "margin-bottom: 20px;"),
      ),

      tags$h4(
        "Step3. Extract Modifications (Optional) ▼",
        onclick = "toggleVisibility('step3_section', this)",
        style = "cursor: pointer; color: #0E4C92;"
      ),
      div(
        id = "step3_section",
        selectInput("mod_datatype", "Select data type", choices = c("PEAKS","Spectronaut","MSFragger",
                    "Comet","DIANN","Skyline","Maxquant","mzIdenML","mzTab"), multiple = FALSE),
        textInput("mod_column", "Column name for modification sequences",
                  placeholder = "e.g., Peptide", value = "Peptide"),
        conditionalPanel(
          condition = "input.mod_datatype === 'MSFragger' || input.mod_datatype === 'mzIdenML' || input.mod_datatype === 'mzTab'",
          textInput("seq_column", "Column name for peptide sequences only for MSFragger, mzIdenML or mzTab")
        ),
        checkboxInput("PTM_annotation", "Annotate PTM with PTM table", value = TRUE),
        h5(strong("PTM table (Editable. Annotation PTM with PTM_type for specified PTM_mass)")),
        DTOutput("ptmTable"),
        fluidRow(
          column(6, actionButton("addRow", "Add Row", icon = icon("plus"))),
          column(6, actionButton("removeRow", "Remove Row", icon = icon("minus")))
        ),
        textInput("PTM_mass_column", "Column name for PTM mass",
                  placeholder = "e.g., PTM_mass", value = "PTM_mass"),
        actionButton("extract", "Extract Modification", class = "btn-primary",
                     style = "margin-bottom: 20px;"),
      ),

      tags$h4(
        "Step4. Match peptide sequences to provided sequence ▼",
        onclick = "toggleVisibility('step4_section', this)",
        style = "cursor: pointer; color: #0E4C92;"
      ),
      div(
        id = "step4_section",
        textInput("match_column", "Column name for peptide sequences to match",
                  placeholder = "e.g., Sequence", value = "Sequence"),
        h5(strong("Sequence table (Editable. Sequence table to be matched on with metadata)")),
        DTOutput("SeqTable"),
        fluidRow(
          column(6, actionButton("addseqRow", "Add Row", icon = icon("plus"))),
          column(6, actionButton("removeseqRow", "Remove Row", icon = icon("minus")))
        ),
        fileInput("uploadseqtable", "Upload Sequence table",
                  multiple = FALSE,
                  accept = c(".csv", ".tsv", ".txt")),
        textInput("match_params", "Column names to match on while matching peptide sequence(comma-separated)",
                  placeholder = "e.g., Molecule"),
        textInput("column_keep", "Column names to keep in result data frame(comma-separated)",
                  placeholder = "e.g., PTM_mass", value = "PTM_mass,PTM_position,reps,Area,Donor,PTM_type"),
        sliderInput(
          "length", "Sequence length range",
          min = 0, max = 100,
          value = c(10, 30)
        ),
        actionButton("match", "Match Sequence", class = "btn-primary",
                     style = "margin-bottom: 20px;"),
      ),
      tags$h4(
        "Step5. Peptide Quantification ▼",
        onclick = "toggleVisibility('step5_section', this)",
        style = "cursor: pointer; color: #0E4C92;"
      ),
      div(
        id = "step5_section",
        selectInput("quantify_method", "Quantification Method:",
                    choices = c("PSM", "Area"),
                    selected = "PSM"),

        # Conditional panel for Area method
        conditionalPanel(
          condition = "input.quantify_method == 'Area'",
          textInput("area_column", "Area Column Name:",
                    placeholder = "e.g., Area")
        ),

        textInput("match_quantify_column", "Matching Columns (comma-separated):",
                  placeholder = "e.g., Chain,Epitope", value = "Chain,Epitope"),
        textInput("distinct_columns", "Distinct Columns (comma-separated):",
                  placeholder = "e.g., Donor", value = "Donor"),
        checkboxInput("with_PTM", "Include PTM Information", value = TRUE),
        checkboxInput("reps", "Include Replicate Information", value = TRUE),
        actionButton("quantify", "Run Quantification", class = "btn-primary")
      ),
      tags$h4(
        "Step6. Merge with Region Information (Optional) ▼",
        onclick = "toggleVisibility('step6_section', this)",
        style = "cursor: pointer; color: #0E4C92;"
      ),
      div(
        id = "step6_section",
        fileInput("region_file", "Upload Region Data (optional)",
                  accept = c(".csv", ".tsv", ".txt")),
        h5(strong("Region table (Editable. Region table to be merged with result)")),
        DTOutput("regionTable"),
        fluidRow(
          column(6, actionButton("addRegionRow", "Add Row", icon = icon("plus"))),
          column(6, actionButton("removeRegionRow", "Remove Row", icon = icon("minus")))
        ),
        textInput("match_region_cols", "Match Column names (comma-separated)", value = "Chain"),
        textInput("region_start_col", "Start Position Column Name in Region table", value = "Region_start"),
        textInput("region_end_col", "End Position Column Name in Region table", value = "Region_end"),
        textInput("region_col", "Region Column Name in Region table", value = "Region"),
        textInput("pos_col", "Position Column Name in processed data", value = "Position"),
        actionButton("merge_regions", "Merge Regions", class = "btn-primary"),
      ),
      tags$h4(
        "Step7. Plot peptides in whole provided sequence ▼",
        onclick = "toggleVisibility('step7_section', this)",
        style = "cursor: pointer; color: #0E4C92;"
      ),
      div(
        id = "step7_section",
        tabsetPanel(
          tabPanel("Plot Settings",
                   textInput("plot_title", "Plot Title:", value = "PSM Plot"),
                   textInput("x_axis", "X-axis Variable(s):", value = "Chain"),
                   textInput("y_axis", "Y-axis Variable(s):", value = "Donor"),
                   textInput("color_fill_column", "Color Fill Column:",
                               value = "PSM"),  # Will be populated server-side
                   fluidRow(
                     column(6, numericInput("fill_min", "Color gradient Min Value:", value = 0)),
                     column(6, numericInput("fill_max", "Color gradient Max Value:", value = 0)),
                     column(6, textInput("low_color", "Color gradient Low Color:", value = "grey80")),
                     column(6, textInput("high_color", "Color gradient High Color:", value = "black"))
                   ),
                   numericInput("label_y", "The position of y axis of the label:", value = 1),
                   h5("Column Ordering"),
                   # Column selection for ordering
                   selectizeInput("order_columns", "Columns to Order:",
                                  choices = NULL, multiple = TRUE,
                                  options = list(placeholder = 'Select columns to define order',
                                                 plugins = list('remove_button'))),

                   # Dynamic UI for selected columns
                   uiOutput("column_order_inputs"),
                   textInput("label_column", "Column name for label:",
                             value = "Character"),
                   textInput("label_filter", "Label Filter (column=value):",
                             value = "Donor=D1"),
                   checkboxInput("add_domain", "Add domain info", value = TRUE),
                   checkboxInput("add_PTM", "Show PTM", value = TRUE),
                   checkboxInput("add_labels", "Show Sequence Labels", value = TRUE),
                   checkboxInput("add_domain_label", "Show Domain Labels", value = TRUE),
                   numericInput("label_size", "Label Size:", value = 1.3, step = 0.1)

          ),

          tabPanel("Domain Settings",
                   h5(strong("Domain definition table (Editable)")),
                   DTOutput("domainTable"),
                   fluidRow(
                     column(6, actionButton("addDomainRow", "Add Row")),
                     column(6, actionButton("removeDomainRow", "Remove Row"))
                   ),
                   fileInput("upload_domain", "Upload Domain Data",
                             accept = c(".csv", ".tsv", ".txt")),
                   fluidRow(
                     column(4, textInput("domain_start_col", "Start Column:", value = "domain_start")),
                     column(4, textInput("domain_end_col", "End Column:", value = "domain_end")),
                     column(4, textInput("domain_type_col", "Type Column:", value = "domain_type"))
                   ),
                   fluidRow(
                     column(4, textInput("domain_border_color_column", "Domain Border Color Column:", value = "domain_color")),
                     column(4, textInput("domain_fill_color_column", "Domain Fill Color Column:", value = "domain_fill_color")),
                     column(4, textInput("domain_label_y_column", "Y axis of domain label Column:", value = "domain_label_y")),
                     column(4, textInput("domain_label_color", "Domain label Color:", value = "black")),
                     column(4, numericInput("domain_label_size", "Domain label Size:", value = 2.5, step = 0.1))
                   )
          ),

          tabPanel("PTM Settings",
                    textInput("PTM_type_column", "PTM type Column:",
                             value = "PTM_type"),
                    h5(strong("PTM color table (Editable)")),
                    DTOutput("PTMcolorTable"),
                    fluidRow(
                     column(6, actionButton("addPTMColorRow", "Add Row")),
                     column(6, actionButton("removePTMColorRow", "Remove Row"))
                    )
          ),

          tabPanel("Advanced Settings",
                   fluidRow(
                     column(6, numericInput("y_expand_low", "Y-axis Lower Expansion:", value = 0.2, step = 0.1)),
                     column(6, numericInput("y_expand_high", "Y-axis Upper Expansion:", value = 0.2, step = 0.1))
                   ),
                   fluidRow(
                     column(6, numericInput("x_expand_low", "X-axis Lower Expansion:", value = 0.5, step = 0.1)),
                     column(6, numericInput("x_expand_high", "X-axis Upper Expansion:", value = 0.5, step = 0.1))
                   ),
                   fluidRow(
                     column(6,
                            selectInput("legend_box", "Legend orientation", choices = c("horizontal", "vertical")),
                     ),
                     column(6,
                            selectInput("legend_position", "Legend Position:",
                                        choices = c("bottom", "right", "left", "top", "none")),
                     ),
                   ),
                   fluidRow(
                     column(6, numericInput("axis_text_size", "Axis Text Size", value = 15)),
                     column(6, numericInput("axis_title_size", "Axis Title Size", value = 15)),
                     column(6, numericInput("legend_title_size", "Legend Title Size", value = 10)),
                     column(6, numericInput("legend_text_size", "Legend Text Size", value = 10))
                   ),
                   h5("Custom Theme Options (Advanced)"),
                   textAreaInput("custom_theme", "theme() parameters as R code:",
                                 placeholder = "e.g., list(plot.margin = margin(20,20,20,20))"),
                   h5("Custom Label Options (Advanced)"),
                   textAreaInput("custom_labs", "labs() parameters as R code:",
                                 placeholder = "e.g., list(fill = 'Intensity', color = 'Domain')"),
          )
        ),
        actionButton("generate_plot", "Generate Plot", class = "btn-primary"),
      ),
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Data Preview",
                 tags$h4(
                   "Raw Data from read file ▼",
                   onclick = "toggleVisibility('raw_section', this)",
                   style = "cursor: pointer;"
                 ),
                 div(
                   id = "raw_section",
                   DT::dataTableOutput("process_dataPreview"),
                 ),
                 tags$h4(
                   "Strip sequence result ▼",
                   onclick = "toggleVisibility('strip_section', this)",
                   style = "cursor: pointer;"
                 ),
                 div(
                   id = "strip_section",
                   DT::dataTableOutput("strip_dataPreview"),
                   downloadButton("downloadStrip", "Download Strip Result",
                                  style = "margin-bottom: 10px;")
                 ),
                 tags$h4(
                   "Extract modification result ▼",
                   onclick = "toggleVisibility('mod_section', this)",
                   style = "cursor: pointer;"
                 ),
                 div(
                   id = "mod_section",
                   DT::dataTableOutput("mod_dataPreview"),
                   downloadButton("downloadMod", "Download Modification Result",
                                  style = "margin-bottom: 10px;")
                 ),
                 tags$h4(
                   "Match Sequence result ▼",
                   onclick = "toggleVisibility('match_section', this)",
                   style = "cursor: pointer;"
                 ),
                 div(
                   id = "match_section",
                   DT::dataTableOutput("match_dataPreview"),
                   downloadButton("downloadMatch", "Download Match Result",
                                  style = "margin-bottom: 10px;")
                 ),
                 tags$h4(
                   "Peptide quantification result ▼",
                   onclick = "toggleVisibility('quantify_section', this)",
                   style = "cursor: pointer;"
                 ),
                 div(
                   id = "quantify_section",
                   DT::dataTableOutput("quantify_dataPreview"),
                   downloadButton("downloadQuant", "Download Peptide Quantification Result",
                                  style = "margin-bottom: 10px;")
                 ),
                 tags$h4(
                   "Region merge result ▼",
                   onclick = "toggleVisibility('region_section', this)",
                   style = "cursor: pointer;"
                 ),
                 div(
                   id = "region_section",
                   DT::dataTableOutput("region_merge_dataPreview"),
                   downloadButton("downloadRegion", "Download Region Merge Result",
                                  style = "margin-bottom: 10px;")
                 )
                 ),
        tabPanel("Visualization",
                 plotOutput("peptidePlot"),
                 selectInput("plot_format", "Select Download File Format:",
                             choices = c("PDF" = "pdf", "PNG" = "png"),
                             selected = "pdf"),
                 numericInput("plot_width", "Plot Width (inches):", value = 30, min = 1, max = 50),
                 numericInput("plot_height", "Plot Height (inches):", value = 6, min = 1, max = 50),
                 downloadButton("downloadPlot", "Download Plot")
                 )
      )
    )
  ),

  # Add JavaScript for toggle functionality
  tags$script(HTML("
    function toggleVisibility(id, element) {
      var section = document.getElementById(id);
      if (section.style.display === 'none') {
        section.style.display = 'block';
        element.innerHTML = element.innerHTML.replace('▼', '▲');
      } else {
        section.style.display = 'none';
        element.innerHTML = element.innerHTML.replace('▲', '▼');
      }
    }
  "))
)
