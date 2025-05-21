#' Report Extra Emissions
#'
#' Calulate additonal emissions based on CEDS 2020 emissions and
#' REMIND/EDGE-Transport Data
#'
#' @param mif a mif file with reported variables from remind2 and EDGE-Transport
#' @param extraData path to extra data files to be used in the reporting
#' @param gdx a GDX as created by readGDX, or the file name of a gdx, needed to define region subsets
#' @author Gabriel Abrahão, Falk Benke
#' @export
reportExtraEmissions <- function(mif, extraData, gdx) {

  # Define region subsets ----
  regionSubsetList <- toolRegionSubsets(gdx)

  # ADD EU-27 region aggregation if possible
  if ("EUR" %in% names(regionSubsetList)) {
    regionSubsetList <- c(regionSubsetList, list(
      "EU27" = c("ENC", "EWN", "ECS", "ESC", "ECE", "FRA", "DEU", "ESW")
    ))
  }

  # Read in selected mif variables ----
  reportVars <- c(
    "ES|Transport|Bunkers|Freight",
    "ES|Transport|Pass|Aviation",
    "FE|Buildings|Gases|Fossil",
    "FE|Buildings|Gases",
    "FE|Buildings|Liquids",
    "FE|Buildings|Solids"
  )

  report <- quitte::read.quitte(mif, check.duplicates = FALSE) %>%
    deletePlus()

  model <- unique(report$model)[grepl("REMIND", unique(report$model))][1] # Deals with REMIND-MAgPIE mifs
  scenario <- unique(report$scenario)[1]

  report <- report %>%
    filter(.data$variable %in% reportVars, .data$region != "World") %>%
    as.magpie() %>%
    collapseDim()


  # Read in additional reporting files ----

  if (!file.exists(file.path(extraData, "p_emissions4ReportExtraCEDS.cs4r"))) {
    stop("Auxiliary file 'p_emissions4ReportExtraCEDS.cs4r' not found")
  }

  if (!file.exists(file.path(extraData, "p_emissions4ReportExtraIAMC.cs4r"))) {
    stop("Auxiliary file 'p_emissions4ReportExtraIAMC.cs4r' not found")
  }

  cedsceds <- read.magpie(file.path(extraData, "p_emissions4ReportExtraCEDS.cs4r"))
  cedsiamc <- read.magpie(file.path(extraData, "p_emissions4ReportExtraIAMC.cs4r"))

  if (!is.null(regionSubsetList)) {
    cedsceds <- mbind(cedsceds, calc_regionSubset_sums(cedsceds, regionSubsetList))
    cedsiamc <- mbind(cedsiamc, calc_regionSubset_sums(cedsiamc, regionSubsetList))
  }

  # check if regional resolution matches
  if (length(setdiff(getItems(report, dim = 1), getItems(cedsceds, dim = 1))) != 0 ||
      length(setdiff(getItems(cedsceds, dim = 1), getItems(report, dim = 1))) != 0
  ) {
    stop("Regional resolution in 'p_emissions4ReportExtraCEDS.cs4r' and gdx do not match.")
  }

  if (length(setdiff(getItems(report, dim = 1), getItems(cedsiamc, dim = 1))) != 0 ||
      length(setdiff(getItems(cedsiamc, dim = 1), getItems(report, dim = 1))) != 0
  ) {
    stop("Regional resolution in 'p_emissions4ReportExtraIAMC.cs4r' and gdx do not match.")
  }


  # Calculate emissions that are based on emission factors.  ----
  # Derive EFs based on CEDS 2020 emissions and REMIND 2020 activities
  .deriveEF <- function(emirefyear, actreffull, refyear = 2020, convyear = NULL) {

    # EF in the reference year
    ef2020 <- setYears(emirefyear / actreffull[, refyear, ], NULL)

    # Preallocate with actreffull to ensure they will multiply nicely later
    ef <- actreffull
    ef[, , ] <- ef2020

    # If convyear was given, assume linear convergence towards
    # the global emission factor in convyear
    if (!is.null(convyear)) {
      gef <- as.numeric(
        dimSums(emirefyear, dim = 1) / dimSums(actreffull[, refyear, ], dim = 1)
      )
      ef <- convergence(
        ef, gef,
        start_year = 2020, end_year = convyear, type = "linear"
      )
    }
    return(ef)
  }

  MtN_to_ktN2O <- 44 / 28 * 1000 # conversion from MtN to ktN2O

  out <- NULL

  # N2O from international shipping
  # Converge to global EF in 2060
  ef <- .deriveEF(
    dimReduce(cedsiamc[, 2020, "International Shipping.n2o_n"]),
    report[, , "ES|Transport|Bunkers|Freight"],
    refyear = 2020,
    convyear = 2060
  )

  out <- mbind(
    out,
    setNames(
      report[, , "ES|Transport|Bunkers|Freight"] * ef * MtN_to_ktN2O,
      "Emi|N2O|Extra|Transport|Bunkers|Freight (kt N2O/yr)"
    )
  )
  # CH4 from international shipping (should be very small)
  # Converge to global EF in 2060
  ef <- .deriveEF(
    dimReduce(cedsiamc[, 2020, "International Shipping.ch4"]),
    report[, , "ES|Transport|Bunkers|Freight"],
    refyear = 2020,
    convyear = 2060
  )


  out <- mbind(
    out,
    setNames(
      report[, , "ES|Transport|Bunkers|Freight"] * ef,
      "Emi|CH4|Extra|Transport|Bunkers|Freight (Mt CH4/yr)"
    )
  )
  # N2O from domestic+international aviation.
  # Converge to global EF in 2060
  ef <- .deriveEF(
    dimReduce(cedsiamc[, 2020, "Aircraft.n2o_n"]),
    report[, , "ES|Transport|Pass|Aviation"],
    refyear = 2020,
    convyear = 2060
  )

  out <- mbind(
    out,
    setNames(
      report[, , "ES|Transport|Pass|Aviation"] * ef * MtN_to_ktN2O,
      "Emi|N2O|Extra|Transport|Pass|Aviation (kt N2O/yr)"
    )
  )
  # CH4 from residential+commercial, assume most of it is from incomplete biomass/solids burning. Requires CEDS detail
  # Don't assume convergence, as Global South EFs may be more representative of solids burning
  ef <- .deriveEF(
    dimReduce(cedsceds[, 2020, "1A4a_Commercial-institutional.ch4"] + cedsceds[, 2020, "1A4b_Residential.ch4"]),
    report[, , "FE|Buildings|Solids"],
    refyear = 2020,
    convyear = NULL
  )
  out <- mbind(
    out,
    setNames(
      report[, , "FE|Buildings|Solids"] * ef,
      "Emi|CH4|Extra|Buildings|Solids (Mt CH4/yr)"
    )
  )
  # N2O from residential+commercial. Requires CEDS detail, assume it's all from fuel burning byproducts.
  # Common solid fuels tend to have higher N2O EFs than common gaseous and liquid fuels, but here
  # we are implicitly assuming the 2020 mix Solids+Liquids+Gases determines the EF.
  # See https://www.epa.gov/system/files/documents/2024-02/ghg-emission-factors-hub-2024.pdf
  # Don't assume convergence, as Global South EFs may be more representative of solids burning
  tmp <- dimSums(report[, , c("FE|Buildings|Gases", "FE|Buildings|Liquids", "FE|Buildings|Solids")], dim = 3)
  ef <- .deriveEF(
    dimReduce(cedsceds[, 2020, "1A4a_Commercial-institutional.n2o_n"] + cedsceds[, 2020, "1A4b_Residential.n2o_n"]),
    tmp,
    refyear = 2020,
    convyear = NULL
  )
  out <- mbind(
    out,
    setNames(
      tmp * ef,
      "Emi|N2O|Extra|Buildings (kt N2O/yr)"
    )
  )

  # Aggregate to global. Since all variables are emissions, we can just sum them
  out <- mbind(out, setItems(dimSums(out, dim = 1), dim = 1, value = "World"))

  # Set emissions to zero that are not represented but that are required for
  # earth system harmonization
  out <- add_columns(out, addnm = "Emi|BC|Extra|Land Use|Peat Burning (Mt BC/yr)",                  dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CH4|Extra|Land Use|Peat Burning (Mt CH4/yr)",                dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CH4|Extra|Energy Demand|International Aviation (Mt CH4/yr)", dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CH4|Extra|Energy Demand|Domestic Aviation (Mt CH4/yr)",      dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CH4|Extra|Energy Demand|Transport (Mt CH4/yr)",              dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CH4|Extra|Energy Demand|Industry (Mt CH4/yr)",               dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CH4|Extra|Industrial Processes (Mt CH4/yr)",                 dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CH4|Extra|Solvents (Mt CH4/yr)",                             dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CO2|Extra|Land Use|Agriculture (Mt CO2/yr)",                 dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CO2|Extra|Land Use|Agricultural Waste Burning (Mt CO2/yr)",  dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CO2|Extra|Land Use|Forest Burning (Mt CO2/yr)",              dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CO2|Extra|Land Use|Savannah Burning (Mt CO2/yr)",            dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CO2|Extra|Land Use|Peat Burning (Mt CO2/yr)",                dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CO2|Extra|Solvents (Mt CO2/yr)",                             dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|CO|Extra|Land Use|Peat Burning (Mt CO/yr)",                  dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|N2O|Extra|Land Use|Peat Burning (kt N2O/yr)",                dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|N2O|Extra|Energy Demand|Industry (kt N2O/yr)",               dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|N2O|Extra|Solvents (kt N2O/yr)",                             dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|NH3|Extra|Land Use|Peat Burning (Mt NH3/yr)",                dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|NOX|Extra|Land Use|Peat Burning (Mt NOX/yr)",                dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|OC|Extra|Land Use|Peat Burning (Mt OC/yr)",                  dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|SO2|Extra|Land Use|Peat Burning (Mt SO2/yr)",                dim = 3.1, fill = 0)
  out <- add_columns(out, addnm = "Emi|VOC|Extra|Land Use|Peat Burning (Mt VOC/yr)",                dim = 3.1, fill = 0)

  # Convert to quitte and ensure it has the same model and scenario as the original report
  out <- as.quitte(out)
  out$model <- model
  out$scenario <- scenario

  return(out)
}
