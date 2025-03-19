#' @keywords internal
linear_regression <- function(formula = NULL,
                       df = NULL,
                       id_col = NULL,
                       link = NULL,
                       alpha_prior = c(0,1),
                       beta_prior = c(0,1)) {
  model_component("linear_regression", {
    if (is.null(formula)) {
      if (is.null(modeldata$.inputs$formula)) {
        cli::cli_abort("Please provide a formula")
      } else {
        formula <- modeldata$.inputs$formula
      }
    }

    if (is.null(df)) {
      if (is.null(modeldata$.inputs$df)) {
        cli::cli_abort("Please provide a `data.frame` with covariates")
      } else {
        df <- modeldata$.inputs$df
      }
    }

    if (is.null(id_col)) {
      if (is.null(modeldata$.inputs$id_col)) {
        cli::cli_abort("Please provide a column name for the sample ids")
      } else {
        id_col <- modeldata$.inputs$id_col
      }
    }

    if (is.null(link)) {
      if (is.null(modeldata$.inputs$link)) {
        cli::cli_abort("Please provide a link function (identity or log)")
      } else {
        link <- modeldata$.inputs$link
      }
    }

    # check for valid link function
    if (link %in% c("identity")) {
      modeldata$link_type <- 0
    } else if (link %in% c("log", "logarithmic")) {
      modeldata$link_type <- 1
    } else {
      cli::cli_abort("The link function must be either 'identity' or 'log'")
    }

    if (!(id_col %in% names(df))) {
      cli::cli_abort(
        c(paste(
          "The following columns must be present",
          "in the provided `data.frame`:",
          paste(id_col, collapse = ", ")
        ),
        paste("Please adjust the `data.frame` or specify the right column",
              "names via the `_col` arguments of this function.")
        )
      )
    }

    df = as.data.table(df)
    df <- setnames(df, old = id_col, new = "sample_id")

    # check that all modeldata$.metainfo$sample_ids are in the df
    if (!all(modeldata$.metainfo$sample_ids %in% df$sample_id)) {
      cli::cli_abort(
        "The provided `data.frame` does not contain covariates for all sample ids"
      )
    }

    # reduce rows to unique sample ids
    df <- df[!duplicated(df$sample_id), ]

    # sort according to sample ids
    setkey(df, sample_id)

    modeldata$X <- model.matrix(formula, data = df)

    # if first column of X corresponds to an intercept, remove it
    if (all(modeldata$X[, 1] == 1)) {
      modeldata$X <- subset(modeldata$X, select = -1)
    }

    modeldata$K <- ncol(modeldata$X)

    modeldata$alpha_prior <- alpha_prior
    modeldata$beta_prior <- beta_prior

    modeldata$.init$alpha <- 1e-2
    modeldata$.init$beta = rep(1e-4, modeldata$K)

    return(modeldata)
  })
}
