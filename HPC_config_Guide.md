Try accessing or search how to connect to SSH on hpc webside
https://imperialcollege.atlassian.net/wiki/spaces/HPC/pages/168624363/Using+SSH


To configure your hpc,\n
log in to your hpc via ssh, go to your directory and try running "anaconda-setup"(might be one level above home)\n
Go to your home directory, place install-dependencies.R, setup_hpc.sh and renv.lock there\n
  Type in terminal: chmod +x setup_hpc.sh\n
  Type in terminal: ./setup_hpc.sh\n
Alternatively, open setup_hpc.sh as a text file and run line by line(attempt if above methods does not work)\n

Use qsub hpc_test.sh to test your configuration\n
