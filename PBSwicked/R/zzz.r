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
  "aes", "age_reading_id", "ageing_method_desc",
  "cell_area", "coord_fixed",
  "count_bait_only", "count_bent_broken", "count_empty_hooks", "count_non_target_species", "count_target_species",
  "delta_gamma", "delta_gengamma", "delta_lognormal", "depth",
  "epsilon_st2", "est", "est2", "est_non_rf2", 
  "facet_wrap", "fishing_event_id", "full_join",
  "geom_raster", "get_depth", "get_index", "ggplot", "ggtitle",
  "hook_adjust_factor",
  "indReaderSame",
  "make_mesh", "maximum_age", "mesh", "minimum_age", "mnam",
  "omega_s2",
  "PlotOutputFn", "Prec_Age", "predictions", "Primary_Age", "prop_bait_hooks",
  "rawindex", "readerDiff", "Resids", "rMx", "RunFn",
  "sameCount", "sanity", "SaveAll", "scale_fill_gradient2", "scale_fill_viridis_c", "sdmTMB", "sdmTMBcontrol",
  "select", "species", "specimen_age", "surv_wUTM", "survey", "survindex",
  "theme_set", "tidy", "total_hooks", "totalCount", "tweedie",
  "X",
  "Y", "year"
	), package="PBSwicked")

wickiverb = c(
"aick"="add|adjust",
"bick"="back|begin|bring",
"cick"="call|change|choose|cut",
"dick"="detect|display|draw|drop",
"eick"="echo|edit|erase",
"fick"="file|filter|find|fix|flush",
"gick"="get|give|grab",
"hick"="hack|halt|hide|hold",
"iick"="import|include|insert",
"jick"="jibe|join|jot|jump",
"kick"="keep|kill",
"lick"="label|link|list|log",
"mick"="make|match|maximize|mend|merge",
"nick"="name|need|note",
"oick"="omit|open|order",
"pick"="pass|pick|plot|post|print|put",
"qick"="query|quit",
"rick"="read|reduce|replace|run",
"sick"="save|select|set|show|sum",
"tick"="take|tabulate|test|trim|try",
"uick"="undo|unite|unpack|update|upload|use",
"vick"="validate|verify|view|visualize",
"wick"="wait|warn|wipe|wrap|write",
"xick"="exit",## wing it for 'x'
"yick"="yap|yield",
"zick"="zap|zip|zoom"
)

