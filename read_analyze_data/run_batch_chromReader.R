# download from https://github.com/mirzahosseini-arash-semmelweis/autoChromReader
# source("../autoChromReader/baseline_flatfit.R")
# source("../autoChromReader/baseline_mile.R")
# source("../autoChromReader/pick_peaks.R")
# source("../autoChromReader/read_chrom.R")
# source("../autoChromReader/peakfinder.R")

path <- "../data/20260811_sampledata"

# Create raw plots for inspection
inspect_chrom_folder(path, baseline_method = "mile")

# Normal first pass that creates lookup table
result <- analyze_chrom_folder_chunked(
  path            = path,
  baseline_method = "mile",
  look_for_t0     = FALSE,
  max_peaks       = 2L,
  n_workers       = 16L
)

saveRDS(
  result,
  "peakfinder_result.rds"
)

# After QC and editing lookup only problematic exports are updated if needed
result <- update_with_lookup(
  result
)

saveRDS(
  result,
  "peakfinder_result_updated.rds"
)

retention <- result$retention
rm(result)

saveRDS(
  retention,
  "retention_data.rds"
)
