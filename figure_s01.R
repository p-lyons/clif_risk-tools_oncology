# figure_s01.R — pooled inclusion flow diagram (coordinating center)
#
# Round-two flow. Each site's figure_s01_flow_<site>.csv now carries SIX rows:
# the starting population, four exclusion steps written by 01_cohort.R, and the
# reconciliation step appended by 02_scores.R ("no calculable score before the
# outcome"). Site files are found recursively under the collection root, so the
# script works whether sites are unpacked as <site>/figure_s01_flow_<site>.csv
# or <site>/upload_to_box_v2/figure_s01_flow_<site>.csv.
#
# Box labels follow the revised eFigure 1 wording: the round-one label
# "outcomes before prediction data available" is replaced by "clinical
# deterioration within 6 hours of ward arrival" (Reviewer 2, comment 19), and
# the reconciliation step is labeled "no score calculable before the outcome"
# (eMethods, cohort exclusions step 5).

# setup ------------------------------------------------------------------------

## libraries -------------------------------------------------------------------

library(data.table)
library(DiagrammeR)
library(tidytable)
library(collapse)
library(stringr)
library(here)

## load all site flow files (recursive) ----------------------------------------

flow_paths = list.files(
  here(),
  pattern    = "^figure_s01_flow_[a-z]+\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

if (length(flow_paths) == 0) stop("No figure_s01_flow_<site>.csv files found.", call. = FALSE)

flow_list = list()

for (fp in flow_paths) {
  file_data      = fread(fp)
  file_data$site = str_remove(str_remove(basename(fp), "^figure_s01_flow_"), "\\.csv$")
  flow_list[[fp]] = file_data
}

flow_data = rbindlist(flow_list, fill = TRUE)

## chronological step index ----------------------------------------------------
# Rows arrive in pipeline order (01_cohort.R writes rows 1-5; 02_scores.R
# appends row 6). An alphabetical sort on the step label does NOT reproduce
# that order, so the lag-based QC below keys on this index, not on the label.

flow_data[, step_idx := seq_len(.N), by = site]

## QC: every site must carry the same six steps --------------------------------

steps_per_site = flow_data[, .(n_steps = .N), by = site]

if (any(steps_per_site$n_steps != 6L)) {
  print(steps_per_site)
  stop("Every site must contribute exactly 6 flow rows (4 exclusions + start + reconciliation).",
       call. = FALSE)
}

# QC: exclusion percentages by site --------------------------------------------

flow_qc = copy(flow_data)
setorder(flow_qc, site, step_idx)

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
  flag_ca = abs(pct_excluded_ca - mean_pct_ca) > 2 * sd_pct_ca,
  flag_no = abs(pct_excluded_no - mean_pct_no) > 2 * sd_pct_no
)]

outliers = flow_qc[flag_ca == TRUE | flag_no == TRUE]

if (nrow(outliers) > 0) {
  cat("\n⚠️  WARNING: Outlier exclusion percentages detected:\n\n")
  print(outliers[, .(site, step,
                     pct_excluded_ca = round(pct_excluded_ca, 1),
                     pct_excluded_no = round(pct_excluded_no, 1))])
  cat("\nReview these before proceeding.\n")
} else {
  cat("\n✓ No outlier exclusion percentages detected\n")
}

flow_qc[order(step_idx, site), .(
  site,
  step,
  pct_excluded_ca = round(pct_excluded_ca, 1),
  pct_excluded_no = round(pct_excluded_no, 1),
  flag_ca,
  flag_no
)]

# summarize data for plotting --------------------------------------------------

flow_data =
  fgroup_by(flow_data, step_idx, step) |>
  fsummarize(
    n_remaining_ca = fsum(n_remaining_ca),
    n_excluded_ca  = fsum(n_excluded_ca),
    n_remaining_no = fsum(n_remaining_no),
    n_excluded_no  = fsum(n_excluded_no)
  ) |>
  ftransform(step = str_replace(step, "hitting", "admission to")) |>
  ftransform(step = str_replace(
    step,
    "outcomes too early",
    "clinical deterioration within 6 hours of ward arrival"
  )) |>
  ftransform(step = str_replace(
    step,
    "no calculable score before the outcome",
    "no score calculable before the outcome"
  ))

step_order = c(
  "Adult inpatient admissions during study period",
  "After excluding patients not admitted through the ED",
  "After excluding patients who were in the ICU before admission to the wards",
  "After excluding encounters with < 6h data",
  "After excluding encounters with clinical deterioration within 6 hours of ward arrival",
  "After excluding encounters with no score calculable before the outcome"
)

if (!all(step_order %in% flow_data$step)) {
  print(setdiff(step_order, flow_data$step))
  stop("Flow steps do not match the expected round-two labels above.", call. = FALSE)
}

flow_data = roworder(flow_data, step_idx)

if (!identical(flow_data$step, step_order)) {
  stop("Chronological row order does not match the expected label order.", call. = FALSE)
}

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

for (i in 2:nrow(flow_data)) {
  current_node = LETTERS[i]
  prev_node    = LETTERS[i - 1]
  excl_node    = paste0("E", i - 1)
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
