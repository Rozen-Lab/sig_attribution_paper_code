Sys.setenv("LD_LIBRARY_PATH" =
             paste(Sys.getenv("LD_LIBRARY_PATH"),
                   "/home/gmssgr/usr/lib64",
                   sep = ":")
           )
Sys.getenv("LD_LIBRARY_PATH")
