.onAttach <- function(...) {
  pkgname <- "keyATM"
  pd <- utils::packageDescription(pkgname);
  packageStartupMessage(pkgname, " ", pd$Version, 
                        " successfully loaded.",
                        "\n Papers, examples, resources, and other materials are at",
                        "\n https://keyatm.github.io/keyATM/");
}
