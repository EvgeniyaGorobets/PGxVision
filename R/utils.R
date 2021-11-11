# Consider taking in multiple vectors, just a single vector
# Select out all biomarkers matching the chosen experiment; used for all plots
selectExperiment <- function(biomarkerDf, experiment) {
  # Validate the experiment vector
  checkmate::assertVector(experiment)
  checkmate::assertNames(names(experiment),
                         permutation.of=c("tissue", "compound", "mDataType"))

  if (is.null(experiment) | !("mDataType" %in% names(experiment) &
                              "tissue" %in% names(experiment) & "compound" %in% names(experiment))) {
    stop("You must choose an experiment to plot. Please provide a vector with
         a tissue, a compound, and a molecular data type (mDataType).")
  }

  # Select only biomarkers from the chosen experiment
  selectedBiomrks <- biomarkerDf[tissue == experiment["tissue"], ]
  selectedBiomrks <- selectedBiomrks[compound == experiment["compound"], ]
  selectedBiomrks <- selectedBiomrks[mDataType == experiment["mDataType"], ]

  # Check that the resulting df is not empty
  if (nrow(selectedBiomrks) == 0) {
    warning(paste0("There were no experiments with the combination ",
                   sprintf("tissue=%s, compound=%s, mDataType=%s.\n",
                           tissue, compound, mDataType),
                   "Manhattan plot will be empty."))
  }

  return(selectedBiomrks)
}
