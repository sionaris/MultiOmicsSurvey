module purge; \
module load rhel8/cclake/base; \
module load r/4.4.0/gcc/3z5n2tph; \
module load netcdf-c/4.9.2/gcc/intel-oneapi-mpi/izgf3n3g; \
export R_LIBS_USER=$PWD/R/library; \
mkdir -p "$R_LIBS_USER"; \
NC_CONFIG="$(command -v nc-config)"; \
test -x "$NC_CONFIG" || { echo "ERROR: nc-config not found (module load failed?)"; exit 1; }; \
Rscript -e 'lib <- Sys.getenv("R_LIBS_USER"); nc <- Sys.getenv("NC_CONFIG"); install.packages("ncdf4", lib=lib, repos="https://cran.r-project.org", configure.args=paste0("--with-nc-config=", nc)); install.packages("HiClimR", lib=lib, repos="https://cran.r-project.org")'

