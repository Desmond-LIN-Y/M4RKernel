To configure your hpc,
log in to your hpc, go to your directory and try running "anaconda-setup"(might be one level above home)
Go to your home directory, place install-dependencies.R, setup_hpc.sh and renv.lock there
  Type in terminal: chmod +x setup_hpc.sh
  Type in terminal: ./setup_hpc.sh
Alternatively, open setup_hpc.sh as a text file and run line by line(attempt if above methods does not work)

Use qsub hpc_test.sh to test your configuration
