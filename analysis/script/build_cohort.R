library(fs)
library(purrr)
library(here)
purrr::walk(.x = fs::dir_ls(here("R")), .f = source)

bpc_dat_path <- path('data-raw', 'NSCLC', '3.1-consortium')
ca_ind <- readr::read_csv(path(bpc_dat_path, 'cancer_level_dataset_index.csv'))
cohort <- ca_ind

# select just the stuff we need:
cohort %<>%
  select(
    cohort,
    record_id,
    institution,
    ca_seq
  )

# flow_track monitors attrition at each step for us.
flow_track <- flow_record_helper(cohort, "[1] BPC NSCLC v3.1")

cohort %<>% filter(!(institution %in% 'UHN'))
flow_track %<>% flow_record_helper(cohort, "[2] US cases only", .)

cpt <- readr::read_csv(path(
  bpc_dat_path,
  'cancer_panel_test_level_dataset.csv'
))
samp_kras_g12d <- readr::read_rds(path('data', 'genomic', 'samp_kras_g12d.rds'))
records_kras_g12d <- cpt %>%
  filter(cpt_genie_sample_id %in% samp_kras_g12d) %>%
  pull(record_id) %>%
  unique
cohort <- filter(cohort, record_id %in% records_kras_g12d)

flow_track %<>% flow_record_helper(cohort, "[3] KRAS G12D", .)

cohort <- cohort %>%
  filter(record_id %in% (get_dmet_time(ca_ind)$record_id))
flow_track %<>% flow_record_helper(cohort, "Metastatic", .)

lot <- readr::read_rds(path('data', 'lot.rds'))
cohort <- left_join(
  cohort,
  filter(lot, line_of_therapy %in% 2),
  by = 'record_id'
)
cohort %<>% filter(!is.na(line_of_therapy))
flow_track %<>% flow_record_helper(cohort, "[4] Had 2L therapy", .)

# This point I'm making a call to take the first cancer sequence for everyone.
# Empirically there's one case where someone has two cancer and they occur on the same day.
# No reason to mess with join logic for that:
cohort %<>%
  group_by(record_id) %>%
  arrange(ca_seq) %>%
  slice(1) %>%
  ungroup(.)


# second_lines <- lot %>%
#   filter(line_of_therapy %in% 2, record_id %in% cohort$record_id)
#
# reg <- readr::read_csv(path(bpc_dat_path, 'regimen_cancer_level_dataset.csv'))
# reg %<>%
#   mutate(
#     dob_reg_start_int = pmin(
#       drugs_startdt_int_1,
#       drugs_startdt_int_2,
#       drugs_startdt_int_3,
#       drugs_startdt_int_4,
#       drugs_startdt_int_5,
#       na.rm = T
#     )
#   )
# reg_before_2l <- left_join(
#   select(
#     second_lines,
#     record_id,
#     dob_second_line_start_int = dob_reg_start_int
#   ),
#   reg,
#   by = 'record_id'
# )
# reg_before_2l %<>% filter(dob_second_line_start_int > dob_reg_start_int)
# # Add the annotations about drug classes:
# reg_class <- readr::read_rds(path('data', 'reg_class.rds'))
# reg_before_2l <- left_join(
#   reg_before_2l,
#   distinct(select(reg_class, regimen_drugs, matches("^is_"))),
#   by = 'regimen_drugs',
#   relationship = 'many-to-one'
# )
# reg_before_2l_sum <- reg_before_2l %>%
#   group_by(record_id) %>%
#   summarize(
#     across(
#       matches('^is_'),
#       .fns = ~ any(.x, na.rm = T),
#       .names = 'any_{.col}_before_2l'
#     )
#   )
#
# cohort <- left_join(
#   cohort,
#   reg_before_2l_sum,
#   by = 'record_id'
# )
# cohort %<>%
#   replace_na(list(
#     any_is_plat_before_2l = FALSE,
#     any_is_anti_pd1_before_2l = FALSE
#   ))

drug_strings <- readr::read_rds(here('data', 'drug_string_matches.rds'))
cohort %<>%
  mutate(
    is_chemo = str_detect(
      tolower(regimen_drugs),
      paste(drug_strings$chemo, collapse = '|')
    ),
    is_anti_pd1 = str_detect(
      tolower(regimen_drugs),
      paste(drug_strings$anti_pd1, collapse = '|')
    )
  )


# directly add it into the flow tracker because these are now not-cumulative.
flow_track <- cohort %>%
  filter(is_chemo) %>%
  flow_record_helper(., "[5a] 2L contains chemo", flow_track)

flow_track <- cohort %>%
  filter(is_anti_pd1) %>%
  flow_record_helper(., "[5b] 2L contains anti-PD-(L)1", flow_track)

flow_track <- cohort %>%
  filter(is_chemo & is_anti_pd1) %>%
  flow_record_helper(., "[5c] 2L contains both", flow_track)

cohort %<>% filter(any_is_anti_pd1_before_2l)
flow_track %<>% flow_record_helper(cohort, "2L contains anti-PD-L1", .)

readr::write_rds(
  flow_track,
  path('data', 'flow_track.rds')
)
