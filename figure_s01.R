
# setup ------------------------------------------------------------------------

## libraries -------------------------------------------------------------------

library(data.table)
library(DiagrammeR)
library(tidytable)
library(collapse)
library(stringr)
library(here)

## function to open nested tables from all sites and rowbind -------------------

read_grouped_files = function(main_folder, stem) {
  
  site_folders = list.dirs(main_folder, recursive = FALSE, full.names = FALSE)
  
  file_list = list()
  
  for(site in site_folders) {
    
    file_path = file.path(main_folder, site, paste0(stem, "_", site, ".csv"))
    
    if(file.exists(file_path)) {
      file_data         = fread(file_path)
      file_data$site    = site
      file_list[[site]] = file_data
    }
  }
  
  combined = rbindlist(file_list, fill = TRUE)
  
  return(combined)
}

# load and QC flow diagram data ------------------------------------------------

flow_data = read_grouped_files(here(), stem = "figure_s01_flow") 

## exclusion percentages by site -----------------------------------------------

flow_qc = copy(flow_data)
setorder(flow_qc, site, step)

flow_qc[, `:=`(
  pct_excluded_ca = (n_excluded_ca / shift(n_remaining_ca, type = "lag")) * 100,
  pct_excluded_no = (n_excluded_no / shift(n_remaining_no, type = "lag")) * 100
), by = site]

flow_qc = flow_qc[!is.na(n_excluded_ca)]

## mean and SD of exclusion % across sites for each step -----------------------

exclusion_summary = flow_qc[, .(
  mean_pct_ca = mean(pct_excluded_ca, na.rm = TRUE),
  sd_pct_ca   = sd(pct_excluded_ca,   na.rm = TRUE),
  mean_pct_no = mean(pct_excluded_no, na.rm = TRUE),
  sd_pct_no   = sd(pct_excluded_no,   na.rm = TRUE)
), by = step]

## identify outliers -----------------------------------------------------------

flow_qc = merge(flow_qc, exclusion_summary, by = "step")

flow_qc[, `:=`(
  flag_ca = abs(pct_excluded_ca-mean_pct_ca) > 2 * sd_pct_ca,
  flag_no = abs(pct_excluded_no-mean_pct_no) > 2 * sd_pct_no
)]

outliers = flow_qc[flag_ca == TRUE | flag_no == TRUE]

outliers 

if(nrow(outliers) > 0) {
  cat("\n⚠️  WARNING: Outlier exclusion percentages detected:\n\n")
  print(outliers[, .(site, step, 
                     pct_excluded_ca = round(pct_excluded_ca, 1),
                     pct_excluded_no = round(pct_excluded_no, 1))])
  cat("\nReview these before proceeding. Uncomment the stop() line below to halt execution.\n")
} else {cat("\n✓ No outlier exclusion percentages detected\n")}

flow_qc[order(step, site), .(
  site, 
  step, 
  pct_excluded_ca = round(pct_excluded_ca, 1),
  pct_excluded_no = round(pct_excluded_no, 1),
  flag_ca,
  flag_no
)]

# summarize data for plotting --------------------------------------------------

flow_data =
  fgroup_by(flow_data, step) |>
  fsummarize(
    n_remaining_ca = fsum(n_remaining_ca),
    n_excluded_ca  = fsum(n_excluded_ca),
    n_remaining_no = fsum(n_remaining_no),
    n_excluded_no  = fsum(n_excluded_no)
  ) |>
  ftransform(step = str_replace(step, "hitting",   "admission to")) |>
  ftransform(step = str_replace(step, "too early", "before prediction data available")) 
  
step_order = c(
  "Adult inpatient admissions during study period",
  "After excluding patients not admitted through the ED",
  "After excluding patients who were in the ICU before admission to the wards",
  "After excluding encounters with < 6h data",
  "After excluding encounters with outcomes before prediction data available"
)

flow_data = flow_data[match(step_order, flow_data$step), ]

# first node for plotting ------------------------------------------------------

format_n = function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

ca_start = 
  fsubset(flow_data, step == "Adult inpatient admissions during study period") |>
  pull(n_remaining_ca)

no_start = 
  fsubset(flow_data, step == "Adult inpatient admissions during study period") |>
  pull(n_remaining_no)
  
diagram_code = paste0(
  "digraph flowchart {
    node [shape = box, fontname = Helvetica]
    A [label = 'Adult inpatient admissions during study period\\nPatients with Cancer: ", format_n(ca_start), 
    "\\nPatients Without Cancer: ",                                                       format_n(no_start), 
    "']"
)

## add subsequent nodes and connections ----------------------------------------

for(i in 2:nrow(flow_data)) {
  current_node = LETTERS[i]
  prev_node    = LETTERS[i-1]
  excl_node    = paste0("E", i-1)
  excl_label   = gsub("After excluding ", "", flow_data$step[i])
  
  diagram_code = paste0(
    diagram_code, "  ",             excl_node, 
    " [label = 'Excluded ",         excl_label, 
    "\\nPatients With Cancer: ",    format_n(flow_data$n_excluded_ca[i]), 
    "\\nPatients Without Cancer: ", format_n(flow_data$n_excluded_no[i]), "']\n",
    "  ", current_node, " 
    [label = 'Remaining Encounters\\nPatients With Cancer: ", format_n(flow_data$n_remaining_ca[i]), 
    "\\nPatients Without Cancer: ",                           format_n(flow_data$n_remaining_no[i]), "']\n",
    "  ", prev_node, " -> ", excl_node,    "\n",
    "  ", prev_node, " -> ", current_node, "\n"
  )
}

diagram_code = paste0(diagram_code, "}")

grViz(diagram_code)

