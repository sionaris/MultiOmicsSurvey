# Parse all sessionInfo files and consolidate package versions
# Use the most recent version when conflicts exist

library(data.table)

setwd("c:/Users/as3582/OneDrive - University of Cambridge/Desktop/GitHub/MultiOmicsSurvey")

parse_sessioninfo <- function(file_path) {
  lines <- readLines(file_path)
  
  # Find all package_version patterns (e.g., ggplot2_3.5.1)
  pkg_pattern <- "[A-Za-z][A-Za-z0-9.]*_[0-9]+\\.[0-9]+[0-9.\\-]*"
  
  matches <- unlist(regmatches(lines, gregexpr(pkg_pattern, lines, perl = TRUE)))
  
  if (length(matches) == 0) return(NULL)
  
  # Split into package and version
  parts <- strsplit(matches, "_")
  data.table(
    package = sapply(parts, `[`, 1),
    version = sapply(parts, `[`, 2),
    source = basename(file_path)
  )
}

# Get all sessionInfo files
files <- list.files("sessionInfo", pattern = "\\.txt$", full.names = TRUE)
cat("Found", length(files), "sessionInfo files\n\n")

# Parse all files
all_pkgs <- rbindlist(lapply(files, parse_sessioninfo), fill = TRUE)

# For each package, get all versions and pick the highest
pkg_summary <- all_pkgs[, .(
  versions = paste(unique(version), collapse = ", "),
  n_versions = uniqueN(version),
  max_version = max(version)
), by = package][order(package)]

# Show packages with version conflicts
conflicts <- pkg_summary[n_versions > 1]
cat("=== PACKAGES WITH VERSION CONFLICTS ===\n")
cat("Total packages with conflicts:", nrow(conflicts), "\n\n")
print(conflicts[order(-n_versions)], nrows = 100)

# Summary
cat("\n\n=== SUMMARY ===\n")
cat("Total unique packages:", nrow(pkg_summary), "\n")
cat("Packages with single version:", sum(pkg_summary$n_versions == 1), "\n")
cat("Packages with version conflicts:", sum(pkg_summary$n_versions > 1), "\n")

# Save consolidated list
fwrite(
  pkg_summary[, .(package, recommended_version = max_version, all_versions = versions)],
  "sessionInfo/consolidated_packages.csv"
)
cat("\nSaved to sessionInfo/consolidated_packages.csv\n")
