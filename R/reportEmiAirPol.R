#' Read in GDX and calculate air pollution emissions, used in convGDX2MIF.R for
#' the reporting
#'
#' Read in air pollution emission information from GDX file, information used in
#' convGDX2MIF.R for the reporting
#'
#'
#' @param gdx a GDX object as created by readGDX, or the path to a gdx
#' @param regionSubsetList a list containing regions to create report variables region
#' aggregations. If NULL (default value) only the global region aggregation "GLO" will
#' be created.
#' @param t temporal resolution of the reporting, default:
#' t=c(seq(2005,2060,5),seq(2070,2110,10),2130,2150)
#'
#' @return MAgPIE object - contains the emission variables
#' @author Antoine Levesque, Jerome Hilaire
#' @seealso \code{\link{convGDX2MIF}}
#' @examples
#' \dontrun{
#' reportEmiAirPol(gdx)
#' }
#'
#' @export
#' @importFrom gdx readGDX
#' @importFrom magclass new.magpie mbind getNames<- setNames collapseNames getYears getItems
reportEmiAirPol <- function(gdx, regionSubsetList = NULL, t = c(seq(2005, 2060, 5), seq(2070, 2110, 10), 2130, 2150)) {
  ######### initialisation  ###########
  tmp <- NULL
  out <- NULL

  ######### initialisation  ###########
  airpollutants <- c("so2", "bc", "oc", "CO", "VOC", "NOx", "NH3")

  ######### internal function  ###########
  generateReportingEmiAirPol <- function(pollutant, i_emiAPexsolve = p11_emiAPexsolve, i_emiAPexo = p11_emiAPexo) {
    poll_rep <- toupper(pollutant)
    tmp <- NULL

    # reduce to the pollutant
    emiAPexsolve <- collapseNames(i_emiAPexsolve[, , pollutant])
    emiAPexo <- collapseNames(i_emiAPexo[, , pollutant])
    getSets(emiAPexo) <- getSets(emiAPexsolve)

    # add indprocess to indst
    emiAPexsolve[, , "indst"] <- emiAPexsolve[, , "indst"] + emiAPexsolve[, , "indprocess"]

    # Case distinction to ensure backwarsds compatibility with older REMIND versions
    # The old version will be removed from remind2 with Release 3.5.2
    if ("Waste" %in% getNames(emiAPexo)) {
      ## Old version: Sector "Waste" is part of p11_emiAPexo
      # Replace REMIND sector names by reporting ones
      mapping <- data.frame(
        remind = c("power", "indst", "res", "trans", "solvents", "extraction"),
        reporting = c(
          paste0("Emi|", poll_rep, "|Energy Supply|Electricity (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Energy Demand|Industry (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Energy Demand|Buildings (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Energy Demand|Transport|Ground Trans (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Solvents (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Energy Supply|Extraction (Mt ", poll_rep, "/yr)")
        )
      )

      emiAPexsolve <- setNames(emiAPexsolve[, , mapping$remind], as.character(mapping$reporting))

      tmp <- mbind(
        emiAPexsolve,
        setNames(emiAPexo[, , "AgWasteBurning"], paste0("Emi|", poll_rep, "|Land Use|Agricultural Waste Burning (Mt ", poll_rep, "/yr)")),
        setNames(emiAPexo[, , "Agriculture"], paste0("Emi|", poll_rep, "|Land Use|Agriculture (Mt ", poll_rep, "/yr)")),
        setNames(emiAPexo[, , "ForestBurning"], paste0("Emi|", poll_rep, "|Land Use|Forest Burning (Mt ", poll_rep, "/yr)")),
        setNames(emiAPexo[, , "GrasslandBurning"], paste0("Emi|", poll_rep, "|Land Use|Savannah Burning (Mt ", poll_rep, "/yr)")),
        setNames(emiAPexo[, , "Waste"], paste0("Emi|", poll_rep, "|Waste (Mt ", poll_rep, "/yr)"))
      )
    } else {
      ## New version: Sector "waste" is part of p11_emiAPexsolve
      # Replace REMIND sector names by reporting ones
      mapping <- data.frame(
        remind = c("power", "indst", "res", "trans", "solvents", "extraction", "waste"),
        reporting = c(
          paste0("Emi|", poll_rep, "|Energy Supply|Electricity (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Energy Demand|Industry (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Energy Demand|Buildings (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Energy Demand|Transport|Ground Trans (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Solvents (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Energy Supply|Extraction (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Waste (Mt ", poll_rep, "/yr)")
        )
      )

      emiAPexsolve <- setNames(emiAPexsolve[, , mapping$remind], as.character(mapping$reporting))

      tmp <- mbind(
        emiAPexsolve,
        setNames(emiAPexo[, , "AgWasteBurning"], paste0("Emi|", poll_rep, "|Land Use|Agricultural Waste Burning (Mt ", poll_rep, "/yr)")),
        setNames(emiAPexo[, , "Agriculture"], paste0("Emi|", poll_rep, "|Land Use|Agriculture (Mt ", poll_rep, "/yr)")),
        setNames(emiAPexo[, , "ForestBurning"], paste0("Emi|", poll_rep, "|Land Use|Forest Burning (Mt ", poll_rep, "/yr)")),
        setNames(emiAPexo[, , "GrasslandBurning"], paste0("Emi|", poll_rep, "|Land Use|Savannah Burning (Mt ", poll_rep, "/yr)"))
      )
    }

    # Set NAs to 0
    tmp[is.na(tmp)] <- 0

    return(tmp)
  }

  ####### read in needed data #########
  ## sets
  ttot <- as.numeric(readGDX(gdx, name = c("ttot"), format = "first_found"))
  ## parameter
  p11_emiAPexsolve <- readGDX(
    gdx,
    name = c("p11_emiAPexsolve", "pm_emiAPexsolve"), field = "l", format = "first_found"
  )[, ttot, ]
  p11_emiAPexo <- readGDX(
    gdx,
    name = c("p11_emiAPexo", "pm_emiAPexo"), field = "l", format = "first_found"
  )[, ttot, airpollutants]

  ####### prepare parameter ########################
  getNames(p11_emiAPexsolve) <- gsub("SOx", "so2", getNames(p11_emiAPexsolve))
  getNames(p11_emiAPexsolve) <- gsub("NMVOC", "VOC", getNames(p11_emiAPexsolve))

  ####### calculate reporting parameters ############
  # Loop over air pollutants and call reporting generating function
  out <- do.call("mbind", lapply(airpollutants, generateReportingEmiAirPol))

  # Add global values
  out <- mbind(out, dimSums(out, dim = 1))
  # add other region aggregations
  if (!is.null(regionSubsetList)) {
    out <- mbind(out, calc_regionSubset_sums(out, regionSubsetList))
  }

  # Loop over air pollutants and add some variables
  for (pollutant in airpollutants) {
    poll_rep <- toupper(pollutant)
    # Aggregation: Transport and Energy Supply
    ## Aviation and International Shipping air pollutant emissions are computed in reportExtraEmissions,
    ## and thus neither included in Energy Demand|Transport not in the totals. 
    out <- mbind(
      out,
      setNames(
        dimSums(out[
          , ,
          c(
            paste0("Emi|", poll_rep, "|Energy Demand|Transport|Ground Trans (Mt ", poll_rep, "/yr)")
          )
        ], dim = 3),
        paste0("Emi|", poll_rep, "|Energy Demand|Transport (Mt ", poll_rep, "/yr)")
      ),
      setNames(
        dimSums(out[
          , ,
          c(
            paste0("Emi|", poll_rep, "|Energy Supply|Electricity (Mt ", poll_rep, "/yr)"),
            paste0("Emi|", poll_rep, "|Energy Supply|Extraction (Mt ", poll_rep, "/yr)")
          )
        ], dim = 3),
        paste0("Emi|", poll_rep, "|Energy Supply (Mt ", poll_rep, "/yr)")
      )
    )
    # Aggregation: Energy Demand + Energy Supply, Land Use
    out <- mbind(
      out,
      setNames(
        dimSums(out[, , c(
          paste0("Emi|", poll_rep, "|Energy Demand|Industry (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Energy Demand|Buildings (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Energy Demand|Transport (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Energy Supply (Mt ", poll_rep, "/yr)")
        )], dim = 3),
        paste0("Emi|", poll_rep, "|Energy Supply and Demand (Mt ", poll_rep, "/yr)")
      ),
      setNames(
        dimSums(out[, , c(
          paste0("Emi|", poll_rep, "|Land Use|Savannah Burning (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Land Use|Forest Burning (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Land Use|Agriculture (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Land Use|Agricultural Waste Burning (Mt ", poll_rep, "/yr)")
        )], dim = 3),
        paste0("Emi|", poll_rep, "|Land Use (Mt ", poll_rep, "/yr)")
      )
    )
    # Compute total
    out <- mbind(
      out,
      setNames(
        dimSums(out[, , c(
          paste0("Emi|", poll_rep, "|Energy Supply and Demand (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Solvents (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Land Use (Mt ", poll_rep, "/yr)"),
          paste0("Emi|", poll_rep, "|Waste (Mt ", poll_rep, "/yr)")
        )], dim = 3),
        paste0("Emi|", poll_rep, " (Mt ", poll_rep, "/yr)")
      )
    )
  }

  getSets(out)[3] <- "variable"
  return(out)
}
