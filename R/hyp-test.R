# Hypothesis testing

#' Test uncertainty estimates against chi-square distribution
#'
#' @param valdata A data.frame as returned by \code{rt_valdata()}
#' @param debias Remove effect of bias? Defaults to FALSE
#' @param sides 1 for 1-sided, or 2 for 2-sided hypothesis test. If 1, test is on upper tail.
#' @param log.p Report p-value as log-transformed value?
#' @importFrom dplyr group_by mutate summarize
#' @export
rt_hyptest <- function(valdata, debias = FALSE, sides = 2, log.p = FALSE) {
  out <- valdata %>%
    mutate(relerr = pixc_err / sigma_est) %>%
    group_by(variable) %>%
    mutate(meanrelerr = mean(relerr, na.rm = TRUE))

  if (!debias) out[["meanrelerr"]] <- 0

  out <- summarize(out,
                   teststat = sum((relerr - meanrelerr)^2, na.rm = TRUE),
                   df = sum(!is.na(relerr)) - 1)
  logpval <- pchisq(out$teststat, df = out$df, log.p = TRUE, lower.tail = FALSE)
  if (sides == 2) {
    logpval <- log(2) + ifelse(logpval > log(0.5), log(1 - exp(logpval)), logpval)
  }

  out$pval <- if (log.p) logpval else exp(logpval)

  out
}

