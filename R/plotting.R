# Plotting functions

#' Plot a histogram of errors, possibly scaled by uncertainty estimates
#'
#' Currently only implemented for node-level data
#'
#' @param valdata data.frame as returned by code{rt_valdata()}
#' @param center Subtract the mean (bias-correct) the errors?
#' @param scale Scale the errors by 1-sigma uncertainty estimates?
#' @param curve Overlay a standard normal curve? This is the default if
#'   both \code{center} and \code{scale} are \code{TRUE}.
#' @param vars Which variables to plot? Defaults to "all"
#' @param plot if FALSE, return the data used to plot but do not create the plot
#' @param ... Passed to \code{geom_histogram()}
#' @importFrom dplyr group_by ungroup
#' @export

rt_val_hist <- function(valdata, center = FALSE, scale = FALSE,
                        curve = center && scale,
                        vars = "all", plot = TRUE,
                        ...) {
  dots <- list(...)

  if (is.null(dots$bins)) dots$bins <- 15

  # Currently only implemented for nodes

  if (length(vars) > 1 || vars != "all") {
    valdata <- valdata %>%
      dplyr::filter(variable %in% vars)
  }

  plotdf <- valdata %>%
    mutate(err = pixc_err) %>%
    group_by(variable)

  # center and scale plot, if necessary
  if (center) {
    plotdf <- mutate(plotdf, err = err - mean(err, na.rm = TRUE))
  }
  if (scale) {
    plotdf <- mutate(plotdf, err = err / sigma_est)
  }

  plotdf <- ungroup(plotdf)

  if (!plot) return(plotdf)

  out <- plotdf %>%
    ggplot(aes(x = err)) +
    geom_histogram(aes(y = ..density..), bins = dots$bins) +
    facet_wrap(~variable, scales = "free")

  if(curve) {
    out <- out +
      stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "blue")
  }
  out
}

#' Plot predictions and uncertainty bounds across nodes
#'
#' @param valdata data.frame as returned by code{rt_valdata()}
#' @param variable Which variable to plot?
#' @param err_only Only plot the errors (subtract off truth value)?
#' @param ... Passed to \code{rt_valdata()}
#'
#' @importFrom dplyr mutate
#' @export

rt_val_nodeseries <- function(valdata, variable = "wse", err_only = TRUE, ...) {

  valdata <- valdata %>%
    `[`(.$variable == variable, ) %>%
    mutate(value = pixc_val, err = pixc_err)

  if (err_only) {
    valdata$baseval <- 0
  } else {
    valdata$baseval <- valdata$value
  }
  # return(valdata)
  yvar <- ifelse(err_only, "err", "value")

  out <- valdata %>%
    ggplot(aes(x = node_id)) +
    geom_ribbon(aes(ymin = baseval - 1.96 * sigma_est,
                    ymax = baseval + 1.96 * sigma_est), fill = "pink") +
    geom_ribbon(aes(ymin = baseval - sigma_est,
                    ymax = baseval + sigma_est), fill = "#7780ff") +
    geom_point(aes_string(y = yvar)) +
    theme_bw()
  out
}

#' Validation scatterplot
#'
#' @param valdata As returned by \code{rt_valdata()}
#' @param variables Variables to plot
#' @param xvar,yvar x- and y-axis variables
#' @param ci1,ci2 Confidence-interval levels for 2 ribbon geoms.
#'  Either or can be set to FALSE to disable plotting.
#' @param plot if FALSE, just return the data used for plotting
#' @export
rt_val_scatter <- function(valdata, variables = c("wse", "width", "slope"),
                           xvar = "id", yvar = c("value", "err", "relerr"),
                           ci1 = 0.6827, ci2 = 0.95,
                           plot = TRUE) {
  varnames <- c(value = "pixc_val", err = "pixc_err", relerr = "rel_err")
  yvar <- match.arg(yvar)
  yvarname <- yvar # for plot axis label
  yvar <- varnames[yvar]

  idvar <- ifelse(is.null(valdata$node_id), "reach_id", "node_id")
  if (xvar == "id") xvar <- idvar

  plotdata <- valdata %>%
    dplyr::filter(variable %in% variables) %>%
    dplyr::mutate(rel_err = pixc_err / sigma_est)

  # Manually add x and y axis variables based on inputs
  plotdata[["xval"]] <- plotdata[[xvar]]
  plotdata[["yval"]] <- plotdata[[yvar]]
  plotdata[["id"]] <- plotdata[[idvar]]

  if (yvar == "pixc_val") {
    plotdata[["ymiddle"]] <- plotdata[["gdem_val"]]
  } else {plotdata[["ymiddle"]] <- 0}
  if (yvar == "rel_err") {
    plotdata$ysigma <- 1
  } else {plotdata$ysigma <- plotdata$sigma_est}

  # Sigma multipliers for confidence intervals
  if (ci1 > 0) {
    stopifnot(ci1 < 1)
    sigmult1 <- qnorm((1 + ci1) / 2)
    plotdata <- plotdata %>%
      mutate(ci1_lwr = ymiddle - sigmult1 * ysigma,
             ci1_upr = ymiddle + sigmult1 * ysigma)
  }
  if (ci2 > 0) {
    stopifnot(ci2 < 1)
    sigmult2 <- qnorm((1 + ci2) / 2)
    plotdata <- plotdata %>%
      mutate(ci2_lwr = ymiddle - sigmult2 * ysigma,
             ci2_upr = ymiddle + sigmult2 * ysigma)
  }

  # Return the data or build the plot
  if (!plot) return(plotdata)

  out <- ggplot(plotdata, aes(x = xval))

  # Add ribbons
  if (ci1 > 0) {
    out <- out + geom_ribbon(aes(ymin = ci2_lwr, ymax = ci2_upr),
                             fill = "pink")
  }
  if (ci2 > 0) {
    out <- out + geom_ribbon(aes(ymin = ci1_lwr, ymax = ci1_upr),
                             fill = "#7780ff")
  }

  # Add points, wrap variables, label
  out <- out + geom_point(aes(y = yval, text = id)) +
    facet_wrap(~variable, scales = "free") +
    ylab(yvarname) + xlab(xvar)

  out
}


#' Map a given number of nodes' pixcvec locations--gdem versus rivertile
#'
#' @param dir A directory containing rivertile output.
#' @param nodes A vector of node_id's, defaults to worst nodes using \code{badnodes()}
#' @param pcv1,pcv2 Name of pixcvec netcdfs
#' @param pixc1,pixc2 Name of pixel cloud netcdfs
#' @param maxpixels Maximum number of pixels to plot, as a random sample if necessary.
#'
#'
#' @export
val_map_node <- function(dir, nodes = badnodes(rt_valdata(dir)),
                         pcv1 = "pcv.nc", pcv2 = "pcv_gdem.nc",
                         pixc1 = "pixel_cloud.nc", pixc2 = "fake_pixc.nc",
                         maxpixels = 1000) {

  if (!requireNamespace("leaflet", quietly = TRUE))
    stop("The `leaflet` package is required for this function. Please install it.")


  keepvars <- c("reach_index", "node_index", "latitude", "longitude",
                "water_frac", "classification", "pixel_area")
  pcvdata1 <- join_pixc(dir, pcvname = pcv1, pixcname = pixc1)[keepvars] %>%
    dplyr::filter(node_index %in% nodes)
  pcvdata2 <- join_pixc(dir, pcvname = pcv2, pixcname = pixc2)[keepvars] %>%
    dplyr::filter(node_index %in% nodes)

  if (nrow(pcvdata1) > maxpixels) {
    message(sprintf("Subsampling. Change using `maxpixels` argument to > %s.",
                    nrow(pcvdata1)))
    pcvdata1 <- dplyr::sample_n(pcvdata1, maxpixels)
  }
  if (nrow(pcvdata2) > maxpixels) {
    message(sprintf("Subsampling. Change using `maxpixels` argument to > %s.",
                    nrow(pcvdata2)))
    pcvdata2 <- dplyr::sample_n(pcvdata2, maxpixels)
  }

  # color palette
  classes <- c(1, 2, 3, 4, 22, 23, 24)
  classlabs <- c("land", "land_near_water", "water_near_land", "open_water",
                 "land_near_dark_water", "dark_water_edge", "dark_water")
  classcolors <- RColorBrewer::brewer.pal(length(classes), "Set1")
  classpal <- colorFactor(palette = classcolors, domain = classes)

  out <- leaflet::leaflet() %>%
    leaflet::addTiles() %>%
    leaflet::addCircles(data = pcvdata1,
                        lng = ~longitude, lat = ~latitude,
                        radius = ~sqrt(pixel_area / pi),
                        color = ~classpal(classification),
                        stroke = FALSE, fillOpacity = 0.8) %>%
    leaflet::addCircles(data = pcvdata2, lng = ~longitude, lat = ~latitude,
                        radius = ~sqrt(pixel_area / pi), color = "red",
                        stroke = FALSE, fillOpacity = 0.8) %>%
    leaflet::addLegend(colors = classcolors, labels = classlabs)

  out
}



# Area plot ---------------------------------------------------------------



#' Plot the area computation for a node
#'
#' @param pixc_joined A joined pixc data frame, as returned by \code{join_pixc()}
#' @param nodes Vector giving indices of nodes to plot
#' @param node_truth Optional truth data, as returned by \code{rt_read()}
#' @param plot if FALSE, return the plot data but don't construct the ggplot
#' @export
nodearea_plot <- function(pixc_joined, nodes, node_truth = NULL, plot = TRUE) {
  sumrydf <- pixc_joined %>%
    dplyr::filter(node_index %in% nodes) %>%
    dplyr::mutate(water_frac = ifelse(classification < 20, water_frac, 1)) %>%
    group_by(node_index) %>%
    dplyr::arrange(desc(water_frac)) %>%
    mutate(cum_area = cumsum(pixel_area),
           area_lag = dplyr::lag(cum_area, default = 0)) %>%
    ungroup() %>%
    mutate(classification = as.factor(classification))

  if (!is.null(node_truth)) {
    joindf <- node_truth %>%
      transmute(reach_index = reach_id, node_index = node_id,
                true_area = area_total)
    sumrydf <- sumrydf %>%
      left_join(joindf, by = c("node_index", "reach_index"))
  }

  if (!plot) return(sumrydf)

  out <- ggplot(sumrydf)

  out <- out +
    geom_rect(aes(xmin = area_lag, xmax = cum_area,
                    ymin = 0, ymax = water_frac, fill = classification,
                    text = pixel_id)) +
    xlab("Cumulative Pixel Area (m^2)") + ylab("Pixel Water Fraction") +
    facet_wrap(~node_index, scales = "free_x")

  # Add gdem truth
  if (!is.null(node_truth)) {
    out <- out +
      geom_rect(aes(xmin = 0, ymin = 0, xmax = true_area, ymax = 1),
                fill = NA, color = "gray30", linetype = 2)
  }

  out
}


#' Cumulative relative error for reaches
#'
#' Plots cumulative relative error against node index, sorted by
#'  node-level relative error magnitude.
#'
#' @param valdf As returned by \code{rt_valdata()}
#' @param reach_ids Vector of \code{reach_id}s to include, or "all"
#' @param variables Vector of variables to include, or "all"
#' @param desc If TRUE, sort x-axis by descending node-level relative error
#' @param plot if FALSE, return the data used to plot, but do not generate the plot.
#'
#' @export
cumerr_plot <- function(valdf, reach_ids = "all",
                        variables = c("wse", "width", "area_total"),
                        desc = FALSE, plot = TRUE) {

  if (length(reach_ids) == 1 && reach_ids == "all")
    reach_ids <- unique(valdf$reach_id)
  if (length(variables) == 1 && variables == "all")
    variables <- unique(valdf$variable)

  rankfun <- function(x) if (desc) -x else x

  plotdf <- valdf %>%
    dplyr::filter(reach_id %in% reach_ids,
                  variable %in% variables) %>%
    mutate(reach_id = as.factor(reach_id)) %>%
    group_by(variable, reach_id) %>%
    mutate(relerr = pixc_err / sigma_est) %>%
    arrange(rankfun(relerr)) %>%
    mutate(cum_err = cumsum(pixc_err),
           cum_sigma = sqrt(cumsum(sigma_est^2)),
           cum_relerr = cum_err / cum_sigma,
           relerr_rank = rank(rankfun(relerr))) %>%
    ungroup()

  ggplot(plotdf, aes(x = relerr_rank)) +
    geom_line(aes(y = cum_relerr, color = reach_id, group = reach_id)) +
    facet_wrap(~variable, scales = "free")
}


#' Leave-one-out relative error for reaches
#'
#' Plots total relative error against nodes, excluding that node from
#'  the calculation of reach-level relative error.
#'
#' @param valdf As returned by \code{rt_valdata()}
#' @param reach_ids Vector of \code{reach_id}s to include, or "all"
#' @param variables Vector of variables to include, or "all"
#' @param sort How to sort the x-axis? Choices are: \code{id} for node_id,
#'  or \code{relerr} for relative error.
#' @param plot if FALSE, return the data used to plot, but do not generate the plot.
#' @export
looerr_plot <- function(valdf, reach_ids = "all",
                        variables = c("wse", "width", "area_total"),
                        sort = c("id", "relerr"), plot = TRUE) {

  sort <- match.arg(sort)

  if (length(reach_ids) == 1 && reach_ids == "all")
    reach_ids <- unique(valdf$reach_id)
  if (length(variables) == 1 && variables == "all")
    variables <- unique(valdf$variable)


  loosum <- function(x, na.rm = FALSE) sum(x, na.rm = na.rm) - x

  plotdf0 <- valdf %>%
    dplyr::filter(reach_id %in% reach_ids,
                  variable %in% variables) %>%
    mutate(reach_id = as.factor(reach_id),
           relerr = pixc_err / sigma_est,
           ranks = if (sort == "id") rank(node_id)
                     else rank(relerr))

  plotdf <- plotdf0 %>%
    group_by(variable, reach_id) %>%
    arrange(ranks) %>%
    mutate(cum_err = loosum(pixc_err),
           cum_sigma = sqrt(loosum(sigma_est^2)),
           cum_relerr = cum_err / cum_sigma,
           xval = rank(ranks),
           node_ind = node_id - min(node_id)) %>%
    ungroup()

  ggplot(plotdf, aes(x = xval)) +
    geom_line(aes(y = cum_relerr, color = reach_id, group = reach_id)) +
    facet_wrap(~variable, scales = "free")
}
