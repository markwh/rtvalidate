# other rivertile-related functions.
# Starting this script initially as a repository for old,
# inconsistent functions that I might still use but that no longer have
# a place in rivertile package.


#' Join pixcvec to pixel cloud
#'
#' Now replaced (but with different usage) by \code{rivertile::pixc_join()}
#'
#' @param dir Directory containing pixel_cloud netcdfs
#' @param pcvname,pixcname Names of pixcvec and pixel cloud netcdf files
#' @param type What kind of join to perform? Default is "inner"
#' @importFrom fs path
#' @importFrom dplyr inner_join
#' @export
#'
join_pixc <- function(dir, pcvname = "pcv.nc",
                        pixcname = "pixel_cloud.nc",
                        type = c("inner", "outer")) {

    type <- match.arg(type)
    joinfun <- if (type == "inner") dplyr::inner_join else dplyr::full_join

    pcvdf <- pixcvec_read(path(dir, pcvname))
    pixcdf <- pixc_read(path(dir, pixcname))

    outdf <- pixcdf %>%
      joinfun(pcvdf, by = c("azimuth_index", "range_index")) %>%
      dplyr::mutate(pixel_id = 1:nrow(.))
    outattdf <- bind_rows2(list(attr(pcvdf, "atts"),
                                attr(pixcdf, "atts")),
                           addMissing = TRUE)

    out <- structure(outdf, atts = outattdf)
  }
