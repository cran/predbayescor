.First.lib <- function(lib,pkg)
{
   library.dynam("predbayescor",pkg,lib)
   cat("predbayescor loaded\n", 
       "COPY RIGHT (c) Longhai Li\n", 
       "http://math.usask.ca/~longhai\n",sep="")
}

