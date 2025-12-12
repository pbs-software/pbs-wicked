# Taking cue from Roger Bivand's maptools:
.PBSwickEnv <- new.env(FALSE, parent=globalenv())  # be sure to exportPattern("^\\.PBS") in NAMESPACE

.onAttach <- function(lib, pkg)
{
	pkg_info = utils::sessionInfo( package="PBSwicked" )$otherPkgs$PBSwicked
	if( is.character( pkg_info$Packaged ) )
		pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]
	else
		pkg_date  <- date()

	userguide_path <- system.file( "doc/PBStools-UG.pdf", package = "PBStools" )
	year <- substring(date(),nchar(date())-3,nchar(date()))

	packageStartupMessage("
-----------------------------------------------------------
PBS Wicked ", pkg_info$Version, " -- Copyright (C) 2025-",year," Fisheries and Oceans Canada

Packaged on ", pkg_date, "
Pacific Biological Station, Nanaimo

All available PBS packages can be found at
https://github.com/pbs-software

'Simplicity is the ultimate sophistication' (Leonardo da Vinci)
-----------------------------------------------------------

")
}
.onUnload <- function(libpath) {
	rm(.PBSwickEnv)
}

# No Visible Bindings
# ===================
if(getRversion() >= "2.15.1") utils::globalVariables(names=c(
  "aes",
  "cell_area", "coord_fixed",
  "delta_gamma", "delta_gengamma", "delta_lognormal", "depth",
  "epsilon_st2", "est2", "est_non_rf2", 
  "facet_wrap", "full_join",
  "geom_raster", "get_index", "ggplot", "ggtitle",
  "make_mesh", "mesh", "mnam",
  "omega_s2",
  "predictions", 
  "rawindex",
  "sanity", "scale_fill_gradient2", "scale_fill_viridis_c", "sdmTMB", "sdmTMBcontrol",
  "select", "surv_wUTM", "survey", 
  "theme_set", "tidy", "tweedie",
  "X",
  "Y"
	), package="PBSwicked")


