# To be used in conjunction with hpc_test.sh
# expected outcome: a file containing a string of numbers 
# appearing in appropriate directory

cat("Starting")
x <- c(1:9)
outpath = "/rds/general/user/dl1220/home/M4RKernel"
cat("Output path:", outpath, "\n")
write.csv(x, file.path(outpath, 'output', 'test.csv'), row.names = F)
cat("Finish")