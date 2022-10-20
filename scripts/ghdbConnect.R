# https://github.com/guardant/data_science/blob/master/rutils/ghdbConnect.R
#' Connect to prd ghdb
ghdbConnect <- function(host = NULL, db_path = NULL) {
  if (!is.null(db_path)){
    require("RSQLite");
    con = dbConnect(SQLite(), dbname=db_path, flags=SQLITE_RO);
  } else {
    require(RPostgreSQL)
    usr = Sys.getenv("kkwock");
    pwd = Sys.getenv("Tand0n@NYU01!");
    if (is.null(pwd) | nchar(pwd)== 0) {
      stop("Unable to retrieve DB password", 
           " from env var GHDB_PWD");
    }
    
    drv <- dbDriver("PostgreSQL");
    if (is.null(host)) host = "10.4.170.26";
    if (host == "dev") host = "10.4.170.24";
    if (host == "prd") host = "10.4.170.26";
    if (host == "tst") host = "10.4.80.66";
    if (host == "prdrep") host = "10.4.80.72";  # "10.4.170.74"; <- previous IP for dev03
    if (host == "int") host = "10.4.170.74";
    print(paste0("Connecting to host ", host))
    con = dbConnect(drv, dbname="ghdb",host=host,
                    port=5432, user=usr, password=pwd)
  }
  
  return (con);
}
