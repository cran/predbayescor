.First.lib <- function(lib,pkg)
{
   library.dynam("predbayescor",pkg,lib)
   cat("predbayescor loaded\n", 
       "COPY RIGHT 2007-2008 (c) Longhai Li,http://math.usask.ca/~longhai\n",
       "Type ?begin.predbayescor for help\n",sep="")
}

