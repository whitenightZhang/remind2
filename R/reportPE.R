#' Read in GDX and calculate primary energy, used in convGDX2MIF.R for the
#' reporting
#'
#' Read in primary energy information from GDX file, information used in
#' convGDX2MIF.R for the reporting
#'
#'
#' @param gdx a GDX as created by readGDX, or the file name of a gdx
#' @param regionSubsetList a list containing regions to create report variables region
#' aggregations. If NULL (default value) only the global region aggregation "GLO" will
#' be created.
#' @param t temporal resolution of the reporting, default:
#' t=c(seq(2005,2060,5),seq(2070,2110,10),2130,2150)
#'
#' @author Lavinia Baumstark, Fabrice Lécuyer
#' @examples
#' \dontrun{
#' reportPE(gdx)
#' }
#'
#' @export
#' @importFrom gdx readGDX
#' @importFrom magclass mselect getYears getNames<- mbind setNames

reportPE <- function(gdx, regionSubsetList = NULL, t = c(seq(2005, 2060, 5), seq(2070, 2110, 10), 2130, 2150)) {
  ####### conversion factors ##########
  TWa_2_EJ <- 3600 * 24 * 365 / 1e6
  ####### read in needed data #########
  ## sets
  teCCS    <- readGDX(gdx, "teCCS", format = "first_found") # technologies with carbon capture
  teNoCCS  <- readGDX(gdx, "teNoCCS", format = "first_found") # technologies without CCS
  entySe   <- readGDX(gdx, "entySe", format = "first_found") # secondary energy types
  peFos    <- readGDX(gdx, "peFos", format = "first_found") # primary energy fossil fuels
  peBio    <- readGDX(gdx, "peBio", format = "first_found") # biomass primary energy types
  pe2se    <- readGDX(gdx, "pe2se") # map primary energy carriers to secondary
  pc2te    <- readGDX(gdx, "pc2te", format = "first_found") # prod couple: mapping for own consumption of technologies
  pc2te    <- pc2te[(pc2te$all_enty1 %in% entySe) & (pc2te$all_enty2 %in% entySe), ] # ensure main and couple product are valid entySe

  seLiq    <- intersect(c("seliqfos", "seliqbio"), entySe)
  seGas    <- intersect(c("segafos", "segabio"), entySe)
  seSol    <- intersect(c("sesofos", "sesobio"), entySe)

  ## parameter
  pm_costsPEtradeMp <- readGDX(gdx, "pm_costsPEtradeMp", restore_zeros = FALSE)
  
  ## variables
  demPE  <- readGDX(gdx, name = "vm_demPe", field = "l", restore_zeros = FALSE, format = "first_found") * TWa_2_EJ
  demPE  <- demPE[pe2se]
  prodSE <- readGDX(gdx, name = "vm_prodSe", field = "l", restore_zeros = FALSE, format = "first_found") * TWa_2_EJ
  prodSE <- mselect(prodSE, all_enty1 = entySe)
  fuExtr <- readGDX(gdx, "vm_fuExtr", field = "l", format = "first_found") * TWa_2_EJ
  Mport  <- readGDX(gdx, "vm_Mport", field = "l", format = "first_found") * TWa_2_EJ
  Xport  <- readGDX(gdx, "vm_Xport", field = "l", format = "first_found") * TWa_2_EJ

  ####### calculate minimal temporal resolution #####
  y <- Reduce(intersect, list(getYears(demPE), getYears(prodSE)))
  demPE  <- demPE[, y, ]
  prodSE <- prodSE[, y, ]
  fuExtr <- fuExtr[, y, ]
  Mport  <- Mport[, y, ]
  Xport  <- Xport[, y, ]

  #### adjust regional dimension of prodCouple (own consumption of technologies)
  prodCouple_tmp <- readGDX(gdx, "pm_prodCouple", restore_zeros = FALSE, format = "first_found")
  prodCouple_tmp[is.na(prodCouple_tmp)] <- 0
  prodCouple <- new.magpie(getRegions(prodSE), getYears(prodCouple_tmp), magclass::getNames(prodCouple_tmp), fill = 0)
  prodCouple[getRegions(prodCouple_tmp), , ] <- prodCouple_tmp
  getSets(prodCouple) <- getSets(prodCouple_tmp)
  prodCouple[prodCouple < 0] <- 0 # fix negative values to 0
  
  ####### internal functions for reporting ###########
  get_demPE <- function(PEcarrier, SEcarrier = entySe, te = pe2se$all_te, name = NULL) {
    # Compute the PE of technologies that have SEcarrier as their main product
    PE <- dimSums(mselect(demPE, all_enty = PEcarrier, all_enty1 = SEcarrier, all_te = te), dim = 3)

    # Some technologies output a couple of SE carriers: the main product (all_enty1), and the couple product (all_enty2)
    # Add the couple PE of technologies that have SEcarrier as their couple product (and another carrier as main product)
    pc2te_couple <- pc2te[(pc2te$all_enty %in% PEcarrier) & (pc2te$all_enty2 %in% SEcarrier) & (pc2te$all_te %in% te), ]
    PE <- PE + dimSums(demPE[pc2te_couple] * prodCouple[pc2te_couple] / (1 + prodCouple[pc2te_couple]), dim = 3)
    # Subtract the couple PE of technologies that have SEcarrier as their main product, but also have a couple product
    pc2te_main <- pc2te[(pc2te$all_enty %in% PEcarrier) & (pc2te$all_enty1 %in% SEcarrier) & (pc2te$all_te %in% te), ]
    PE <- PE - dimSums(demPE[pc2te_main] * prodCouple[pc2te_main] / (1 + prodCouple[pc2te_main]), dim = 3)

    if (!is.null(name)) magclass::getNames(PE) <- name
    return(PE)
  }

  get_prodSE <- function(PEcarrier, SEcarrier = entySe, te = pe2se$all_te, name) {
    setNames(dimSums(mselect(prodSE, all_enty = PEcarrier, all_enty1 = SEcarrier, all_te = te), dim = 3), name)
  }

  ####### calculate reporting parameters ############
  out <- mbind(
    # fossil fuels
    get_demPE(peFos,                                                 name = "PE|Fossil (EJ/yr)"),
    get_demPE(peFos, te = teCCS,                                     name = "PE|Fossil|+|w/ CC (EJ/yr)"),
    get_demPE(peFos, te = teNoCCS,                                   name = "PE|Fossil|+|w/o CC (EJ/yr)"),

    get_demPE("pecoal",                                              name = "PE|+|Coal (EJ/yr)"),
    get_demPE("pecoal", te = teCCS,                                  name = "PE|Coal|++|w/ CC (EJ/yr)"),
    get_demPE("pecoal", te = teNoCCS,                                name = "PE|Coal|++|w/o CC (EJ/yr)"),
    get_demPE("pecoal", "seel",                                      name = "PE|Coal|+|Electricity (EJ/yr)"),
    get_demPE("pecoal", "seel", teCCS,                               name = "PE|Coal|Electricity|+|w/ CC (EJ/yr)"),
    get_demPE("pecoal", "seel", teNoCCS,                             name = "PE|Coal|Electricity|+|w/o CC (EJ/yr)"),
    get_demPE("pecoal", seGas,                                       name = "PE|Coal|+|Gases (EJ/yr)"),
    get_demPE("pecoal", seGas, teCCS,                                name = "PE|Coal|Gases|+|w/ CC (EJ/yr)"),
    get_demPE("pecoal", seGas, teNoCCS,                              name = "PE|Coal|Gases|+|w/o CC (EJ/yr)"),
    get_demPE("pecoal", seLiq,                                       name = "PE|Coal|+|Liquids (EJ/yr)"),
    get_demPE("pecoal", seLiq, teCCS,                                name = "PE|Coal|Liquids|+|w/ CC (EJ/yr)"),
    get_demPE("pecoal", seLiq, teNoCCS,                              name = "PE|Coal|Liquids|+|w/o CC (EJ/yr)"),
    get_demPE("pecoal", "seh2",                                      name = "PE|Coal|+|Hydrogen (EJ/yr)"),
    get_demPE("pecoal", "seh2", teCCS,                               name = "PE|Coal|Hydrogen|+|w/ CC (EJ/yr)"),
    get_demPE("pecoal", "seh2", teNoCCS,                             name = "PE|Coal|Hydrogen|+|w/o CC (EJ/yr)"),
    get_demPE("pecoal", "sehe",                                      name = "PE|Coal|+|Heat (EJ/yr)"),
    get_demPE("pecoal", seSol,                                       name = "PE|Coal|+|Solids (EJ/yr)"),

    get_demPE("peoil",                                               name = "PE|+|Oil (EJ/yr)"),
    get_demPE("peoil", te = teCCS,                                   name = "PE|Oil|++|w/ CC (EJ/yr)"),
    get_demPE("peoil", te = teNoCCS,                                 name = "PE|Oil|++|w/o CC (EJ/yr)"),
    get_demPE("peoil", "seel",                                       name = "PE|Oil|+|Electricity (EJ/yr)"),
    get_demPE("peoil", seLiq,                                        name = "PE|Oil|+|Liquids (EJ/yr)"),
    
    get_demPE("pegas",                                               name = "PE|+|Gas (EJ/yr)"),
    get_demPE("pegas", te = teCCS,                                   name = "PE|Gas|++|w/ CC (EJ/yr)"),
    get_demPE("pegas", te = teNoCCS,                                 name = "PE|Gas|++|w/o CC (EJ/yr)"),
    get_demPE("pegas", "seel",                                       name = "PE|Gas|+|Electricity (EJ/yr)"),
    get_demPE("pegas", "seel", teCCS,                                name = "PE|Gas|Electricity|+|w/ CC (EJ/yr)"),
    get_demPE("pegas", "seel", teNoCCS,                              name = "PE|Gas|Electricity|+|w/o CC (EJ/yr)"),
    get_demPE("pegas", seGas,                                        name = "PE|Gas|+|Gases (EJ/yr)"),
    get_demPE("pegas", seLiq,                                        name = "PE|Gas|+|Liquids (EJ/yr)"),
    get_demPE("pegas", seLiq, teCCS,                                 name = "PE|Gas|Liquids|+|w/ CC (EJ/yr)"),
    get_demPE("pegas", seLiq, teNoCCS,                               name = "PE|Gas|Liquids|+|w/o CC (EJ/yr)"),
    get_demPE("pegas", "sehe",                                       name = "PE|Gas|+|Heat (EJ/yr)"),
    get_demPE("pegas", "seh2",                                       name = "PE|Gas|+|Hydrogen (EJ/yr)"),
    get_demPE("pegas", "seh2", teCCS,                                name = "PE|Gas|Hydrogen|+|w/ CC (EJ/yr)"),
    get_demPE("pegas", "seh2", teNoCCS,                              name = "PE|Gas|Hydrogen|+|w/o CC (EJ/yr)"),

    # biomass
    get_demPE(peBio,                                                 name = "PE|+|Biomass (EJ/yr)"),
    get_demPE(peBio, te = teCCS,                                     name = "PE|Biomass|++|w/ CC (EJ/yr)"),
    get_demPE(peBio, te = teNoCCS,                                   name = "PE|Biomass|++|w/o CC (EJ/yr)"),
    get_demPE(c("pebioil", "pebios"),                                name = "PE|Biomass|1st Generation (EJ/yr)"),
    get_demPE(peBio, te = "biotr",                                   name = "PE|Biomass|Traditional (EJ/yr)"),
    get_demPE(peBio, "seel",                                         name = "PE|Biomass|+|Electricity (EJ/yr)"),
    get_demPE(peBio, "seel", teCCS,                                  name = "PE|Biomass|Electricity|+|w/ CC (EJ/yr)"),
    get_demPE(peBio, "seel", teNoCCS,                                name = "PE|Biomass|Electricity|+|w/o CC (EJ/yr)"),
    get_demPE(peBio, seGas,                                          name = "PE|Biomass|+|Gases (EJ/yr)"),
    get_demPE(peBio, seGas, teCCS,                                   name = "PE|Biomass|Gases|+|w/ CC (EJ/yr)"),
    get_demPE(peBio, seGas, teNoCCS,                                 name = "PE|Biomass|Gases|+|w/o CC (EJ/yr)"),
    get_demPE(peBio, "seh2",                                         name = "PE|Biomass|+|Hydrogen (EJ/yr)"),
    get_demPE(peBio, "seh2", teCCS,                                  name = "PE|Biomass|Hydrogen|+|w/ CC (EJ/yr)"),
    get_demPE(peBio, "seh2", teNoCCS,                                name = "PE|Biomass|Hydrogen|+|w/o CC (EJ/yr)"),
    get_demPE(peBio, seLiq,                                          name = "PE|Biomass|+|Liquids (EJ/yr)"),
    get_demPE(peBio, seLiq, teCCS,                                   name = "PE|Biomass|Liquids|+|w/ CC (EJ/yr)"),
    get_demPE(peBio, seLiq, teNoCCS,                                 name = "PE|Biomass|Liquids|+|w/o CC (EJ/yr)"),
    get_demPE("pebiolc", seLiq,                                      name = "PE|Biomass|Liquids|Cellulosic (EJ/yr)"),
    get_demPE("pebiolc", seLiq, teCCS,                               name = "PE|Biomass|Liquids|Cellulosic|+|w/ CC (EJ/yr)"),
    get_demPE("pebiolc", seLiq, teNoCCS,                             name = "PE|Biomass|Liquids|Cellulosic|+|w/o CC (EJ/yr)"),
    get_demPE(c("pebioil", "pebios"), seLiq,                         name = "PE|Biomass|Liquids|Non-Cellulosic (EJ/yr)"),
    get_demPE("pebios", seLiq,                                       name = "PE|Biomass|Liquids|Conventional Ethanol (EJ/yr)"),
    get_demPE(peBio, seLiq, c("bioftrec", "bioftcrec", "biodiesel"), name = "PE|Biomass|Liquids|Biodiesel (EJ/yr)"),
    get_demPE(peBio, seLiq, c("bioftcrec"),                          name = "PE|Biomass|Liquids|Biodiesel|+|w/ CC (EJ/yr)"),
    get_demPE(peBio, seLiq, c("bioftrec", "biodiesel"),              name = "PE|Biomass|Liquids|Biodiesel|+|w/o CC (EJ/yr)"),
    get_demPE(peBio, seSol,                                          name = "PE|Biomass|+|Solids (EJ/yr)"),
    get_demPE(peBio, "sehe",                                         name = "PE|Biomass|+|Heat (EJ/yr)"),

    # renewables and nuclear
    # using prodSE rather than demPE to follow the "direct method" for primary energy:
    # https://ourworldindata.org/energy-substitution-method
    get_prodSE(c("pegeo", "pehyd", "pewin", "pesol"),                name = "PE|Non-Biomass Renewables (EJ/yr)"),
    get_prodSE("pegeo",                                              name = "PE|+|Geothermal (EJ/yr)"),
    get_prodSE("pegeo", "sehe",                                      name = "PE|Geothermal|+|Heat (EJ/yr)"),
    get_prodSE("pegeo", "seel",                                      name = "PE|Geothermal|+|Electricity (EJ/yr)"),
    get_prodSE("peur",                                               name = "PE|+|Nuclear (EJ/yr)"),
    get_prodSE("pehyd",                                              name = "PE|+|Hydro (EJ/yr)"),
    get_prodSE("pewin",                                              name = "PE|+|Wind (EJ/yr)"),
    get_prodSE("pewin", te = "windon",                               name = "PE|Wind|+|Onshore (EJ/yr)"),
    get_prodSE("pewin", te = "windoff",                              name = "PE|Wind|+|Offshore (EJ/yr)"),
    get_prodSE("pesol",                                              name = "PE|+|Solar (EJ/yr)"),
    get_prodSE("pesol", "sehe",                                      name = "PE|Solar|+|Heat (EJ/yr)"),
    get_prodSE("pesol", "seel",                                      name = "PE|Solar|+|Electricity (EJ/yr)")
  )

  # advanced calculation
  out <- mbind(out,
    setNames(dimSums(demPE[, , c(peFos, peBio)], dim = 3)
           + dimSums(prodSE[, , c("pegeo", "pehyd", "pewin", "pesol", "peur")], dim = 3), "PE (EJ/yr)"),
    setNames(dimSums(demPE[, , peBio], dim = 3)
           - dimSums(mselect(demPE, all_enty = peBio, all_te = "biotr"), dim = 3),  "PE|Biomass|Modern (EJ/yr)"),
    setNames(fuExtr[, , "pebiolc.2"], "PE|Biomass|Residues (EJ/yr)"),
    setNames((fuExtr[, , "pebiolc.1"]
              + (1 - pm_costsPEtradeMp[, , "pebiolc"]) * Mport[, , "pebiolc"]
              - Xport[, , "pebiolc"]), "PE|Biomass|Energy Crops (EJ/yr)")
  )

  # add global values
  out <- mbind(out, dimSums(out, dim = 1))
  # add other region aggregations
  if (!is.null(regionSubsetList))
    out <- mbind(out, calc_regionSubset_sums(out, regionSubsetList))

  getSets(out)[3] <- "variable"
  return(out)
}
