#' Retrieve statistically significant biomarkers per compound per tissue from
#'   PharmacoDB.ca
#'
#' @description NOTE: this function requires valid PharmacoDB database
#'   credentials to work. Not intended for end users.
#'
#' @param file `character(1)` Filepath of database credential file. Read
#'   into the R environment to populate environmental variables for other
#'   parameters. Defaults to file.path("~", ".mysql").
#' @param url `character(1)` Name of enviromental variable with database
#'   url. Defaults to "PHARMACODB".
#' @param port `character(1)` Name of enviromental variable with database
#'   port. Defaults to "PHARMACODB_PORT".
#' @param user `character(1)` Name of environmental variable with database
#'   username. Defaults to "MYSQL_USER".
#' @param password `character(1)` Name of environmental variable with database
#'   user password. Defaults to "MYSQL_PW".
#'
#' @return `data.table` Table of tissue specific biomarkers from PharmacoDB.ca.
#'
#' @import dbplyr
#' @importFrom dplyr select filter rename inner_join left_join tbl collect
#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom RMariaDB MariaDB
#' @importFrom data.table as.data.table
#' @export
fetchPharmacoDbBiomarkers <- function(mysql=file.path("~", ".mysql"),
        url="PHARMACODB", port="PHARMACODB_PORT", user="MYSQL_USER",
        password="MYSQL_PW") {
    readRenviron(mysql)
    CON <- dbConnect(
        MariaDB(),
        dbname="pharmacodb_test",
        username=Sys.getenv(user),
        password=Sys.getenv(password),
        host=Sys.getenv(url),
        port=Sys.getenv(port)
    )
    on.exit(dbDisconnect(CON))  # disconnect after function execution

    # Get tables of interest
    gene <- tbl(CON, "gene")
    gene_annot <- tbl(CON, "gene_annotation")
    compound <- tbl(CON, "compound")
    compound_annot <- tbl(CON, "compound_annotation")
    tissue <- tbl(CON, "tissue")
    biomarker <- tbl(CON, "gene_compound_tissue")

    # Join
    tissue_df <- tissue |>
        rename(tissue="name")

    gene_df <- gene_annot |>
        filter(symbol != "NA") |>
        select(gene_id, symbol) |>
        left_join(gene, by=c(gene_id="id")) |>
        rename(gene_symbol="symbol", ensembl_id="name")


    compound_df <- compound_annot |>
        select(compound_id, inchikey, smiles, pubchem, chembl_id) |>
        left_join(compound, by=c(compound_id="id")) |>
        select(-compound_uid) |>
        rename(compound_name="name")

    biomarker_df <- biomarker |>
        filter(pvalue < 0.05) |>
        select(gene_id, compound_id, tissue_id, estimate, pvalue) |>
        inner_join(gene_df, by="gene_id") |>
        inner_join(compound_df, by="compound_id") |>
        inner_join(tissue_df, by=c("tissue_id"="id")) |>
        select(-gene_id, -compound_id, -tissue_id) |>
        collect()

    return(as.data.table(biomarker_df))
}