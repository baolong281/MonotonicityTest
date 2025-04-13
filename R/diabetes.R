#' A Simulated Diabetes Dataset
#'
#' This dataset contains simulated medical measurements for Diabetes and is
#' emulated after data from the Diabetes Prevention Program.
#' Each column represents change in a key metabolic indicators after two years for the placebo group
#' receiving no treatment.
#'
#' @format A data frame with 1000 rows and 4 variables:
#' \describe{
#'   \item{CLDL}{Change in low-density lipoprotein (LDL) cholesterol (mg/dL).}
#'   \item{GLUCOSE}{Change in fasting plasma glucose levels (mg/dL).}
#'   \item{TRIG}{Change in triglyceride levels (mg/dL).}
#'   \item{HBA1C}{Change in hemoglobin A1c levels (\%).}
#' }
#' @examples
#' data("diabetes", package="MonotonicityTest")
#' names(diabetes)
#' @usage data("diabetes", package="MonotonicityTest")
"diabetes"
