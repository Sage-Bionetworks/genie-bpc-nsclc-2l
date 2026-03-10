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
flow_track <- flow_record_helper(cohort, "BPC NSCLC v3.1")

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

flow_track %<>% flow_record_helper(cohort, "KRAS G12D", .)

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
flow_track %<>% flow_record_helper(cohort, "Had 2L therapy", .)
