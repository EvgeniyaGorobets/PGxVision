### INPUT ###

# Helper functions
logicalModelFileUploadBox <- box(
  width=12,
  column(
    width=6,
    fileInput("patientLogicalFile", "Upload Patient Molecular Profile CSV:",
              accept=c("text/csv", ".csv"), buttonLabel="Browse files"
    )
  ),
  column(
    width=6,
    fileInput("logicalModelFile",
              "Upload Logical Model JSON:",
              accept=c("application/json", ".json"), buttonLabel="Browse files"
    )
  )
)

# Return all tab input rows 
logicalTabInputUI = tabItem(
  tabName= "logicalModel",
  h2("Multivariate Logical Model Predictor"),
  fluidRow(logicalModelFileUploadBox),
  fluidRow(
    box(width=12,
        uiOutput("logicalModel")
    )
  )
)

### REACTIVE VALUES AND OBSERVERS ###

# object containing all relevant reactive values for the logical tab
logicalTabCreateRV <- function() {
  return( reactiveValues(logicalModelL = NULL,
                         patientLogicalDf = NULL))
}

# Return all reactive variable observers
logicalTabObservers <- function(input, rv) {
  observeEvent(input$patientLogicalFile, {
    req(input$patientLogicalFile)
    df_ <- data.table::fread(input$patientLogicalFile$datapath)
    rv$patientLogicalDf <- df_
  })
  
  observeEvent(input$logicalModelFile, {
    req(input$logicalModelFile)
    json <- jsonlite::read_json(input$logicalModelFile$datapath)
    rv$logicalModelL <- json
  })
}

### OUTPUT ###

logicalTabOutputUI <- function(input, rv, output) {
  output$logicalModel <- renderUI({
    df_ <- rv$patientLogicalDf
    model <- rv$logicalModelL
    # FIXME:: Remove this work around for mismatching gene names
    remap_genes <- c(SERPINA2="^ATR$",
                     MISP="^C19orf21$|^MISP1$|^Caprice$", CT45A7="^CT45A5$", # CT45A5 is paralog, not equilavent
                     HIST1H1T="^H1-6$")
    for (i in seq_along(remap_genes))
      df_[feature %ilike% remap_genes[i], feature := names(remap_genes)[i]]
    stopifnot(!anyDuplicated(df_$feature))
    # assumes duplicated genes have same threshold
    thresholds <- unlist(model$thresholds[!duplicated(names(model$thresholds))])
    df_sub <- df_[feature %in% names(thresholds), ]
    expression <- unlist(dcast(df_sub, . ~ feature)[, -1])
    expr_gt_thresh <- expression[names(thresholds)] > thresholds
    # clean up formula string for evalution
    form <- model$formula
    for (i in seq_along(remap_genes))
      form <- gsub(remap_genes[i], names(remap_genes)[i], form)
    form <- gsub("\\[", "(", form)
    form <- gsub("\\]", ")", form)
    # error if gene names aren't exclusinvely alphanumeric
    stopifnot(all(vapply(names(expr_gt_thresh),
                         FUN=grepl, FUN.VALUE=logical(1),
                         pattern="^[[:alnum:]]+$")))
    model_match <- with(as.list(expr_gt_thresh), eval(str2lang(form)))
    column(width=12,
           fluidRow(width=12,
                    h3(form)
           ),
           fluidRow(width=12,
                    renderText(paste0("\n\nModel match: ", model_match))
           )
    )
  })
}
