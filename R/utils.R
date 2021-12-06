# TODO: Consider taking in multiple vectors, just a single vector
# Select out all biomarkers matching the chosen experiment; used for all plots
# TODO: export
selectExperiment <- function(biomarkerDt, tissueName="", compoundName="",
                             molecularType="") {
  # Local bindings to satisfy check() and DT syntax
  tissue <- compound <- mDataType <- NULL

  # Validate user input
  if (!data.table::is.data.table(biomarkerDt)) {
    stop("biomarkerDt must be of type data.table.")
  }
  checkmate::assertString(tissueName)
  checkmate::assertString(compoundName)
  checkmate::assertString(molecularType)

  # Select only biomarkers from the chosen experiment
  selectedBiomrks <- biomarkerDt[tissue == tissueName, ]
  selectedBiomrks <- selectedBiomrks[compound == compoundName, ]
  selectedBiomrks <- selectedBiomrks[mDataType == molecularType, ]

  # Check that the resulting df is not empty
  if (nrow(selectedBiomrks) == 0) {
    warning(paste0("There were no experiments with the combination ",
      sprintf("tissue=%s, compound=%s, mDataType=%s.\n", tissueName,
      compoundName, molecularType), "Plots will be empty."))
  }

  return(selectedBiomrks)
}

# [END]
