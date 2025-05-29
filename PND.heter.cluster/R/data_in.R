#' data_in
#'
#' A simulated dataset from the 2/1 partially nested design with treatment-incuded clustering
#'
#' @format A data frame with 400 rows and 8 variables:
#' \describe{
#'   \item{Y}{Outcome.}
#'   \item{K}{Cluster assignment in the treatment arm.}
#'   \item{tt}{Treatment assignment. 1 for individuals assigned to the treatment arm. 0 for individuals assigned to the control arm. The control arm is unclustered.}
#'   \item{X_dat.1}{Baseline covariates.}
#'   \item{X_dat.2}{Baseline covariates.}
#'   \item{X_dat.3}{Baseline covariates.}
#'   \item{X_dat.4}{Baseline covariates.}
#'   \item{id}{Individual id.}
#' }
"data_in"
