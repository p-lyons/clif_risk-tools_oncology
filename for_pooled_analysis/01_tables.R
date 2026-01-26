# 01_tables.R
# Table 1 (flow), Table 2 (characteristics), Table 3 (outcomes)

# setup ------------------------------------------------------------------------

library(flextable)
library(officer)

# output directories -----------------------------------------------------------

fig_dir = here::here("output", "figures")
tbl_dir = here::here("output", "tables")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(tbl_dir)) dir.create(tbl_dir, recursive = TRUE)

# TABLE 1: FLOW DIAGRAM (Figure S1) --------------------------------------------

message("\n== Creating flow diagram ==")

## QC: check exclusion percentages by site -------------------------------------

flow_qc = copy(flow_data_raw)
setorder(flow_qc, site, step)

flow_qc[, `:=`(
  pct_excluded_ca = (n_excluded_ca / shift(n_remaining_ca, type = "lag")) * 100,
  pct_excluded_no = (n_excluded_no / shift(n_remaining_no, type = "lag")) * 100
), by = site]

flow_qc = flow_qc[!is.na(n_excluded_ca)]

exclusion_summary = flow_qc[, .(
  mean_pct_ca = mean(pct_excluded_ca, na.rm = TRUE),
  sd_pct_ca   = sd(pct_excluded_ca,   na.rm = TRUE),
  mean_pct_no = mean(pct_excluded_no, na.rm = TRUE),
  sd_pct_no   = sd(pct_excluded_no,   na.rm = TRUE)
), by = step]

flow_qc = merge(flow_qc, exclusion_summary, by = "step")

flow_qc[, `:=`(
  flag_ca = abs(pct_excluded_ca - mean_pct_ca) > 2 * sd_pct_ca,
  flag_no = abs(pct_excluded_no - mean_pct_no) > 2 * sd_pct_no
)]

outliers = flow_qc[flag_ca == TRUE | flag_no == TRUE]

if (nrow(outliers) > 0) {
  message("  WARNING: Outlier exclusion percentages detected:")
  print(outliers[, .(site, step, 
                     pct_excluded_ca = round(pct_excluded_ca, 1),
                     pct_excluded_no = round(pct_excluded_no, 1))])
} else {
  message("  No outlier exclusion percentages detected")
}

## aggregate flow data ---------------------------------------------------------

flow_data = flow_data_raw |>
  fgroup_by(step) |>
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

## create flow diagram (DiagrammeR) --------------------------------------------

library(DiagrammeR)
library(DiagrammeRsvg)

ca_start = flow_data[step == step_order[1]]$n_remaining_ca
no_start = flow_data[step == step_order[1]]$n_remaining_no

diagram_code = paste0(
  "digraph flowchart {
    node [shape = box, fontname = Helvetica]
    A [label = 'Adult inpatient admissions during study period\\nPatients with Cancer: ", format_n(ca_start),
  "\\nPatients Without Cancer: ", format_n(no_start), "']"
)

for (i in 2:nrow(flow_data)) {
  current_node = LETTERS[i]
  prev_node    = LETTERS[i-1]
  excl_node    = paste0("E", i-1)
  excl_label   = gsub("After excluding ", "", flow_data$step[i])
  
  diagram_code = paste0(
    diagram_code, "  ", excl_node,
    " [label = 'Excluded ", excl_label,
    "\\nPatients With Cancer: ",    format_n(flow_data$n_excluded_ca[i]),
    "\\nPatients Without Cancer: ", format_n(flow_data$n_excluded_no[i]), "']\n",
    "  ", current_node,
    " [label = 'Remaining Encounters\\nPatients With Cancer: ", format_n(flow_data$n_remaining_ca[i]),
    "\\nPatients Without Cancer: ", format_n(flow_data$n_remaining_no[i]), "']\n",
    "  ", prev_node, " -> ", excl_node,    "\n",
    "  ", prev_node, " -> ", current_node, "\n"
  )
}

diagram_code = paste0(diagram_code, "}")

flow_diagram = grViz(diagram_code)

# Save flow diagram
export_svg(flow_diagram) |> charToRaw() |> rsvg::rsvg_pdf(file.path(fig_dir, paste0("figS1_flow_", today, ".pdf")))

# TABLE 2: CHARACTERISTICS -----------------------------------------------------

message("\n== Creating Table 2: Characteristics ==")

## prepare data ----------------------------------------------------------------

cat_data  = copy(cat_data_raw)
cont_data = copy(cont_data_raw)

overall_n = cat_data[var == "female", .(N = sum(n)), by = ca_01]

## clean up category labels ----------------------------------------------------

# Race - title case and clean names
cat_data[var == "race_category", category := fcase(
  category == "american indian or alaska native",           "American Indian or Alaska Native",
  category == "asian",                                      "Asian",
  category == "black or african american",                  "Black or African American",
  category == "native hawaiian or other pacific islander",  "Native Hawaiian or Other Pacific Islander",
  category == "white",                                      "White",
  category == "other/unknown",                              "Other/Unknown",
  default = category
)]

# Ethnicity
cat_data[var == "ethnicity_category", category := fcase(
  category == "hispanic",     "Hispanic",
  category == "non_hispanic", "Non-Hispanic",
  category == "unknown",      "Unknown",
  default = category
)]

# Site - uppercase
cat_data[var == "site", category := toupper(category)]

## continuous variables --------------------------------------------------------

### age ------------------------------------------------------------------------

table_cont_age = cont_data[, .(
  n    = sum(n),
  mean = sum(age_sum) / sum(n),
  sd   = calculate_sd_from_sums(sum(age_sum), sum(age_sumsq), sum(n))
), by = ca_01]

means_age     = table_cont_age$mean
sds_age       = table_cont_age$sd
pooled_sd_age = sqrt((sds_age[1]^2 + sds_age[2]^2) / 2)
smd_age       = abs(means_age[2] - means_age[1]) / pooled_sd_age

table_cont_age[, formatted := paste0(round(mean, 1), " (", round(sd, 1), ")")]
age_wide = dcast(table_cont_age, . ~ ca_01, value.var = "formatted")
age_wide[, `:=`(var = "age", SMD = round(smd_age, 3))]
age_wide[, . := NULL]
setcolorder(age_wide, c("var", "0", "1", "SMD"))
setnames(age_wide, c("0", "1"), c("No Cancer", "Cancer"))

### Van Walraven ---------------------------------------------------------------

table_cont_vw = cont_data[, .(
  n    = sum(n),
  mean = sum(vw_sum) / sum(n),
  sd   = calculate_sd_from_sums(sum(vw_sum), sum(vw_sumsq), sum(n))
), by = ca_01]

means_vw     = table_cont_vw$mean
sds_vw       = table_cont_vw$sd
pooled_sd_vw = sqrt((sds_vw[1]^2 + sds_vw[2]^2) / 2)
smd_vw       = abs(means_vw[2] - means_vw[1]) / pooled_sd_vw

table_cont_vw[, formatted := paste0(round(mean, 1), " (", round(sd, 1), ")")]
vw_wide = dcast(table_cont_vw, . ~ ca_01, value.var = "formatted")
vw_wide[, `:=`(var = "vw", SMD = round(smd_vw, 3))]
vw_wide[, . := NULL]
setcolorder(vw_wide, c("var", "0", "1", "SMD"))
setnames(vw_wide, c("0", "1"), c("No Cancer", "Cancer"))

cont_wide = rbindlist(list(age_wide, vw_wide))

## categorical variables -------------------------------------------------------

### collapse categories --------------------------------------------------------

cat_data[var == "initial_code_status" & category %in% c("presume full", "presumed full", "presume_full"), 
         category := "full"]
cat_data[var == "initial_code_status" & category %in% c("dnar/dni", "dnr/dni", "and"), 
         category := "dnr/dni"]
cat_data[var == "initial_code_status" & category %in% c("other", "special/partial", "dnar"), 
         category := "other"]
cat_data[var == "race_category" & category %in% c("other", "unknown", "two or more races"), 
         category := "other/unknown"]

### aggregate ------------------------------------------------------------------

cat_data  = cat_data[, .(n = sum(n)), by = .(ca_01, var, category, site)]
table_cat = cat_data[!var %in% c("age", "vw"), .(n = sum(n)), by = .(ca_01, var, category)]
table_cat = merge(table_cat, overall_n, by = "ca_01")
table_cat[, pct := (n / N) * 100]

### add site counts ------------------------------------------------------------

site_counts = cat_data[var == "female" & category == "1", .(n = sum(n)), by = .(ca_01, site)]
site_counts = merge(site_counts, overall_n, by = "ca_01")
site_counts[, pct := (n / N) * 100]
site_counts[, `:=`(var = "site", category = site)]
site_counts[, site := NULL]
table_cat = rbindlist(list(table_cat, site_counts), use.names = TRUE)

### calculate SMDs -------------------------------------------------------------

table_cat[, prop := n / N]
props_wide = dcast(table_cat, var + category ~ ca_01, value.var = "prop")
setnames(props_wide, c("0", "1"), c("p0", "p1"))
props_wide[, smd := abs(p1 - p0) / sqrt((p0 * (1 - p0) + p1 * (1 - p1)) / 2)]
table_cat[, formatted := paste0(format_n(n), " (", round(pct, 1), "%)")]

### handle binary variables ----------------------------------------------------

binary_vars        = table_cat[, .(is_binary = all(category %in% c("0", "1"))), by = var][is_binary == TRUE]$var
table_cat_filtered = table_cat[(var %in% binary_vars & category == "1") | !(var %in% binary_vars)]
cat_wide           = dcast(table_cat_filtered, var + category ~ ca_01, value.var = "formatted")
setnames(cat_wide, c("0", "1"), c("No Cancer", "Cancer"))

### merge SMDs -----------------------------------------------------------------

props_filtered = props_wide[(var %in% binary_vars & category == "1") | !(var %in% binary_vars)]
cat_wide = merge(
  cat_wide,
  props_filtered[, .(var, category, SMD = round(smd, 3))],
  by = c("var", "category"),
  all.x = TRUE
)

## variable labels -------------------------------------------------------------

var_labels = data.table(
  var = c(
    "N", "age", "vw", "initial_code_status",
    "female", "race_category", "ethnicity_category",
    "site",
    "los", "wicu", "va", "imv",
    "dead", "hospice", "dead_or_hospice"
  ),
  label = c(
    "N", "Age, years, mean (SD)", "Van Walraven score, mean (SD)", "Admission code status, n (%)",
    "Female sex, n (%)", "Race, n (%)", "Ethnicity, n (%)",
    "Site, n (%)",
    "Length of stay, n (%)", "Ward-ICU transfer, n (%)",
    "Vasoactive medications, n (%)", "Invasive mechanical ventilation, n (%)",
    "Died in hospital, n (%)", "Discharged to hospice, n (%)", "Died or discharged to hospice, n (%)"
  )
)

## create N row ----------------------------------------------------------------

n_row = data.table(
  var         = "N",
  category    = NA_character_,
  `No Cancer` = format_n(overall_n[ca_01 == 0]$N),
  Cancer      = format_n(overall_n[ca_01 == 1]$N),
  SMD         = NA_real_
)

## final assembly --------------------------------------------------------------

all_data         = rbindlist(list(n_row, cont_wide, cat_wide), fill = TRUE)
all_data         = merge(all_data, var_labels, by = "var", all.x = TRUE)
categorical_vars = c("race_category", "ethnicity_category", "initial_code_status", "site", "los")

header_rows = all_data[
  var %in% categorical_vars,
  .(var         = unique(var),
    category    = NA_character_,
    `No Cancer` = NA_character_,
    Cancer      = NA_character_,
    SMD         = NA_real_)
]

header_rows = merge(header_rows, var_labels, by = "var")
all_data    = rbindlist(list(all_data, header_rows), fill = TRUE)

all_data[, display := fcase(
  var %in% binary_vars & category == "1",       label,
  var == "dead_or_hospice",                     label,
  var %in% categorical_vars & is.na(category),  label,
  var %in% c("N", "age", "vw"),                 label,
  !is.na(category) & category != "",            paste0("  ", category),
  default = label
)]

setorder(all_data, var, category, na.last = FALSE)

## build Table 2 ---------------------------------------------------------------

char_vars = c("N", "age", "female", "race_category", "ethnicity_category", "site", "vw", "initial_code_status")

table2 = all_data[var %in% char_vars]
table2 = table2[!(var == "vw" & !is.na(category))]

table2[, sort_key := fcase(
  var == "N",                   1L,
  var == "age",                 2L,
  var == "female",              3L,
  var == "race_category",       4L,
  var == "ethnicity_category",  5L,
  var == "vw",                  6L,
  var == "site",                7L,
  var == "initial_code_status", 8L,
  default = 999L
)]

setorder(table2, sort_key, category, na.last = FALSE)
table2[, sort_key := NULL]
table2 = table2[, .(display, `No Cancer`, Cancer, SMD)]
setnames(table2, "display", "Variable")

message("  Table 2 created with ", nrow(table2), " rows")

ft2 = flextable(table2) |> autofit()
save_as_docx(ft2, path = file.path(tbl_dir, paste0("table2_characteristics_", today, ".docx")))

# TABLE 3: OUTCOMES ------------------------------------------------------------

message("\n== Creating Table 3: Outcomes ==")

## create combined death/hospice variable --------------------------------------

death_hospice = table_cat[var %in% c("dead", "hospice") & category == "1", 
                          .(n = sum(n)), by = .(ca_01, N)]
death_hospice[, pct := (n / N) * 100]
death_hospice[, formatted := paste0(format_n(n), " (", round(pct, 1), "%)")]
death_hospice[, prop := n / N]

dh_smd = death_hospice[, .(p0 = prop[ca_01 == 0], p1 = prop[ca_01 == 1])]
dh_smd[, smd := abs(p1 - p0) / sqrt((p0 * (1 - p0) + p1 * (1 - p1)) / 2)]

dh_row = dcast(death_hospice, . ~ ca_01, value.var = "formatted")
dh_row[, `:=`(var = "dead_or_hospice", category = "1", SMD = round(dh_smd$smd, 3))]
dh_row[, . := NULL]
setcolorder(dh_row, c("var", "category", "0", "1", "SMD"))
setnames(dh_row, c("0", "1"), c("No Cancer", "Cancer"))

cat_wide = rbindlist(list(cat_wide, dh_row), use.names = TRUE)

## chi-square tests ------------------------------------------------------------

outcome_vars     = c("dead_or_hospice", "dead", "hospice", "icu", "va", "imv", "los")
outcome_cat_data = table_cat[var %in% c("dead", "hospice", "icu", "va", "imv", "los")]

chi_results = outcome_cat_data[, {
  tab        = dcast(.SD, category ~ ca_01, value.var = "n")
  tab_matrix = as.matrix(tab[, -1])
  test       = suppressWarnings(chisq.test(tab_matrix))
  .(p_value  = test$p.value)
}, by = var]

# death/hospice chi-square
dh_tab    = death_hospice[, .(category = "1", n), keyby = .(ca_01)]
dh_matrix = matrix(c(dh_tab$n, overall_n$N - dh_tab$n), nrow = 2)
dh_pval   = suppressWarnings(chisq.test(dh_matrix))$p.value
chi_results = rbindlist(list(data.table(var = "dead_or_hospice", p_value = dh_pval), chi_results))

## risk differences ------------------------------------------------------------

binary_outcomes = c("dead", "hospice", "icu", "va", "imv", "dead_or_hospice")

rd_results = table_cat[var %in% binary_outcomes & category == "1", {
  prop_0 = n[ca_01 == 0] / N[ca_01 == 0]
  prop_1 = n[ca_01 == 1] / N[ca_01 == 1]
  rd     = (prop_1 - prop_0) * 100
  .(rd = rd)
}, by = var]

dh_rd = death_hospice[, {
  prop_0 = prop[ca_01 == 0]
  prop_1 = prop[ca_01 == 1]
  rd     = (prop_1 - prop_0) * 100
  .(var = "dead_or_hospice", rd = rd)
}]

rd_results = rbindlist(list(rd_results, dh_rd))

# LOS risk differences
los_rd = table_cat[var == "los", {
  prop_0 = n[ca_01 == 0] / N[ca_01 == 0]
  prop_1 = n[ca_01 == 1] / N[ca_01 == 1]
  rd     = (prop_1 - prop_0) * 100
  .(rd = rd)
}, by = .(var, category)]

rd_all = rbindlist(list(rd_results[, category := "1"], los_rd), fill = TRUE)

## build Table 3 ---------------------------------------------------------------

var_labels_t3 = data.table(
  var   = c("dead_or_hospice", "dead", "hospice", "icu", "va", "imv", "los"),
  label = c(
    "Died or discharged to hospice, n (%)",
    "Died in hospital, n (%)",
    "Discharged to hospice, n (%)",
    "ICU admission, n (%)",
    "Vasoactive medications, n (%)",
    "Invasive mechanical ventilation, n (%)",
    "Length of stay, n (%)"
  )
)

# Rebuild all_data with death/hospice
all_data = rbindlist(list(n_row, cont_wide, cat_wide), fill = TRUE)
all_data = merge(all_data, var_labels, by = "var", all.x = TRUE)

table3 = all_data[var %in% outcome_vars]
table3[var == "los" & category == "0-47h",    category := "< 2 days"]
table3[var == "los" & category == "48h_96h",  category := "2-4 days"]
table3[var == "los" & category == "96h_1wk",  category := "4-7 days"]
table3[var == "los" & category == "1wk_2wk",  category := "7-14 days"]
table3[var == "los" & category == "2wk_plus", category := "> 14 days"]

table3 = merge(table3, chi_results, by = "var", all.x = TRUE)
table3 = merge(table3, rd_all, by = c("var", "category"), all.x = TRUE)
table3 = merge(table3, var_labels_t3, by = "var", all.x = TRUE, suffixes = c("_old", ""))

table3[, display := fcase(
  var %in% c("dead", "hospice", "icu", "va", "imv", "dead_or_hospice") & category == "1", label,
  var == "los" & is.na(category), label,
  !is.na(category) & category != "", paste0("  ", category),
  default = label
)]

table3[, sort_key := fcase(
  var == "dead_or_hospice", 1L,
  var == "dead",            2L,
  var == "hospice",         3L,
  var == "icu",             4L,
  var == "va",              5L,
  var == "imv",             6L,
  var == "los",             7L,
  default = 999L
)]

setorder(table3, sort_key, category, na.last = FALSE)

# Reorder LOS categories
los_header = table3[var == "los" &  is.na(category)]
los_cats   = table3[var == "los" & !is.na(category)]
los_cats   = los_cats[match(c("< 2 days", "2-4 days", "4-7 days", "7-14 days", "> 14 days"), category)]
other_rows = table3[var != "los"]
setorder(other_rows, sort_key, category, na.last = FALSE)
table3 = rbindlist(list(other_rows, los_header, los_cats))

table3[, is_first   := !duplicated(var)]
table3[, p_display  := fifelse(is_first == TRUE, sprintf("%.3f", p_value), NA_character_)]
table3[, rd_display := fifelse(!is.na(rd), sprintf("%.1f", rd), NA_character_)]
table3 = table3[, .(display, `No Cancer`, Cancer, `RD (%)` = rd_display, `P-value` = p_display)]
setnames(table3, "display", "Variable")

message("  Table 3 created with ", nrow(table3), " rows")

ft3 = flextable(table3) |> autofit()
save_as_docx(ft3, path = file.path(tbl_dir, paste0("table3_outcomes_", today, ".docx")))

# exports ----------------------------------------------------------------------

message("\n== Tables complete ==")

# Export objects for figures script
table2_final       = table2
table3_final       = table3
flow_diagram_final = flow_diagram
flow_data_final    = flow_data
