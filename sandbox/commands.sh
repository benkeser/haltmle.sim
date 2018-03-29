#-----------------------------------------
# commands for scp'ing sce and cent over
#-----------------------------------------
cd ~/Dropbox/R/haltmle.sim/sandbox
# random mechanisms
scp  cent.R sce.sh dbenkese@snail.fhcrc.org:~/haltmle.sim
# kang and schafer
scp  cent_ks.R sce_ks.sh dbenkese@snail.fhcrc.org:~/haltmle.sim
# sin wave
scp  cent_sin.R sce_sin.sh dbenkese@snail.fhcrc.org:~/haltmle.sim
# scp mergeOnly.R subMerge.sh dbenkese@snail.fhcrc.org:~/haltmle.sim

ssh dbenkese@snail.fhcrc.org
cd haltmle.sim
# random mechanisms 
scp cent.R sce.sh dbenkese@rhino.fhcrc.org:~/haltmle.sim
# kang and schafer
scp cent_ks.R sce_ks.sh dbenkese@rhino.fhcrc.org:~/haltmle.sim
# sin wave
scp cent_sin.R sce_sin.sh dbenkese@rhino.fhcrc.org:~/haltmle.sim

ssh dbenkese@rhino.fhcrc.org
cd haltmle.sim
# random mechanisms
./sce.sh ./cent.R rdg_n5000_v2
# kang and schafer
./sce_ks.sh ./cent_ks.R ks_mod_v2
# sin wave
./sce_sin.sh ./cent_sin.R sin_redos_v2

#-----------------------
cd ~/Dropbox/R/haltmle.sim/sandbox
scp test.R dbenkese@snail.fhcrc.org:~/haltmle.sim
# scp mergeOnly.R subMerge.sh dbenkese@snail.fhcrc.org:~/haltmle.sim

ssh dbenkese@snail.fhcrc.org
cd haltmle.sim
# scp mergeOnly.R subMerge.sh dbenkese@rhino.fhcrc.org:~/haltmle.sim
scp test.R dbenkese@rhino.fhcrc.org:~/haltmle.sim

ssh dbenkese@rhino.fhcrc.org
cd haltmle.sim
chmod +x test.R
sbatch -N 1 test.R

#-----------------------------------------
# commands for scp'ing mergeOnly
#-----------------------------------------
cd ~/Dropbox/R/haltmle.sim/tests
scp mergeOnly.R subMerge.sh dbenkese@snail.fhcrc.org:~/haltmle.sim

ssh dbenkese@snail.fhcrc.org
cd haltmle.sim
scp mergeOnly.R subMerge.sh dbenkese@rhino.fhcrc.org:~/haltmle.sim

ssh dbenkese@rhino.fhcrc.org
cd haltmle.sim
./subMerge.sh

#-----------------------------------------
# commands to get into rhino and load R
#-----------------------------------------
ssh dbenkese@snail.fhcrc.org
ssh dbenkese@rhino.fhcrc.org
 # enter password
ml R/3.2.2
R

#-----------------------------------------
# scp results from rhino to local machine
#-----------------------------------------
# from rhino
cd ~/haltmle.sim/out
scp allOut_ks.RData dbenkese@snail.fhcrc.org:~/haltmle.sim
	# enter snail password
 	# ctrl + shift + t to open up new term
# scp to snail
cd ~/Dropbox/R/haltmle.sim/results
scp dbenkese@snail.fhcrc.org:~/haltmle.sim/allOut_ks.RData . 

#-----------------------------------------
# misc commands 
#-----------------------------------------
# squeue -u dbenkese
# scancel `seq 51239645 51239655`
