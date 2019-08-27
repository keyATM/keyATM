.onAttach <- function(...) {
	pkgname <- "keyATM"
  pd <- utils::packageDescription(pkgname);
  packageStartupMessage(pkgname, " v", pd$Version, 
												" successfully loaded.",
												"\n Papers, examples, resources, and other materials are at xxx.xxx");
}
