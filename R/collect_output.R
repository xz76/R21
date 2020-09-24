### collect output rds files in given directory

in_dir <- "output"
in_files <- list.files(in_dir, pattern = "\\.rds$")

res <- lapply(seq_along(in_files), function(i) {
    readRDS(file.path(in_dir, in_files[i]))
})

out_dir <- "collected"
if (! dir.exists(out_dir)) dir.create(out_dir)
saveRDS(res, file = file.path(out_dir, "one_sim.rds"))
