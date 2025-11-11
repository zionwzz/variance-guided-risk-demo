# scripts/00_setup_env.R
# Reproducible setup: seed, single-thread BLAS/LAPACK, session info once.

set.seed(12345)

# Make BLAS/LAPACK single-threaded for determinism
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

# Write session info on first run
if (!file.exists("session_info.txt")) {
  sink("session_info.txt"); on.exit(sink(), add = TRUE)
  cat("Session info\n--------------\n")
  print(sessionInfo())
}
