#' Stroke prevention therapies for people with atrial fibrillation
#'
#' Contrast-level data.
#'
#' @format A data frame with 49 rows and 14 variables:
#' \describe{
#'   \item{`TE`}{The observed relative treatment effects.}
#'   \item{`seTE`}{The standard error of the observed relative treatment effects.}
#'   \item{`studlab`}{The study labels.}
#'   \item{`treat1`}{The treatment code of the comparator in the observed relative
#'   effect.}
#'   \item{`treat2`}{The treatment code of the intervention in the observed
#'   relative effect.}
#'   \item{`event1`}{The number of participants assigned `treat1` and experienced
#'   a stroke during follow-up.}
#'   \item{`time1`}{The total number of months at risk for participants assigned
#'   `treat1`.}
#'   \item{`event2`}{The number of participants assigned `treat2` and experienced
#'   a stroke during follow-up.}
#'   \item{`time2`}{The total number of months at risk for participants assigned
#'   `treat2`.}
#'   \item{`stroke`}{The proportion of participants in the study who previously
#'   experienced a stroke before randomization.}
#'   \item{`class1_detail`}{The drug class of the treatment in `treat1`.}
#'   \item{`class2_detail`}{The drug class of the treatment in `treat2`.}
#'   \item{`class1_simple`}{`treat1` categorized as "Control" or "Active".}
#'   \item{`class2_simple`}{`treat2` categorized as "Control" or "Active".}
#' }
#'
#' @source Cooper et al. (2009). Addressing between-study heterogeneity and inconsistency in mixed treatment comparisons: Application to stroke prevention treatments in individuals with non-rheumatic atrial fibrillation, Statistics in Medicine 28:1861-1881. [doi:10.1002/sim.3594]. Adopted from gemtc::atrialFibrillation.
#'
#' data(dat_afib)
"dat_afib"
