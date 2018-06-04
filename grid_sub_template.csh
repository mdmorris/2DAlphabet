#! /bin/sh
./development/runManySections.py --createCommandFile --envString ENVSTRING --addLog \TEMPDIR/listOfJobs.txt commands.cmd
./runManySections.py --submitCondor commands.cmd
condor_q lcorcodi
