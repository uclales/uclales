 This contains a short description of the way the script llrunon256procs, which automatically resubmits
 a job a X nb of times, should be used. 

 Suppose you want to make a run of 6 hours, and you know that every 2
 hours of simulation take 7 real time hours. Due to the wall limit of 8
 hours you will need do break the run in 3 small runs of 2 hours and you
 want to do that automatically.

 You can use llsubmit llrunon256procs, which requires to have the files
 nsubmit.dat, ntime.dat and nmaxsubmit.dat in your bin directory.

 before you make the llsubmit of llrunon256procs, you have to put in the
 3 files:

 nmaxsubmit.dat - the nb of times you want to resubmit the job after the
 first 2 hours are finished. In our example 2

 nsubmit. dat- contains the counter of the nb of time the job is
 resubmitted, so at the beginning has to be equal to 0.

 ntime.dat -the timemax of the initial simulation, so in this case 7200.

 I'm changing the directories between the bin and the work directory
 because I want to have the outputs in the work directory, so I have to
 have the namelists there too.

  Now about the namelists,  initially you have a namelist with timemax =
  7200. and runtype = INITIAL. But as the next time you want to have a
  namelist with timemax=14400 and runtype = HISTORY, I use a trick in the
  beginning by copying (just the first time, before doing the llsubmit
  llrunon256procs) creating a NAMELIST_test which is the same as NAMELIST,
  but with runtype=HYSTORY. Then at the first loop, the script will change
  the timemax in NAMELIST_test and copy it to NAMELIST, which would make
  you have the NAMELIST you wanted for your second run, with timemax=14400
  and runtype=HISTORY.

  I include here my initial NAMELIST and NAMELIST_test, just in case I
  wasn't clear enough.


NAMELIST


   &model
     nzp =   528
     nxp =   260 
     nyp =   260
     nxpart = .True.
     deltax = 35.
     deltay = 35.
     deltaz = 5.
     dzrat = 1.1
     dzmax = 2500.
     dtlong = 2.
     distim = 100.
     timmax = 7200.
     runtype = "INITIAL"
     igrdtyp = 1
     level = 3
     CCN = 100.e6
     filprf = 'trDct'
     hfilin = 'trDct.rst'
     savg_intvl = 3600.
     ssam_intvl = 60.
     corflg = .true.
     prndtl = -0.3333333
     iradtyp = 4
     strtim= 166.41 
     th00 = 289.  
     umean = -2.
     vmean = -4.
     cntlat=25.
     isfctyp=2  
     zrough = 0.001229
     dthcon = 0.001094
     drtcon = 0.001133
     case_name='trans'
     lsvarflg=.true.
     div=0.00000186
    /


NAMELIST_test

    &model
    nzp =   528
    nxp =   260 
    nyp =   260
    nxpart = .True.
    deltax = 35.
    deltay = 35.
    deltaz = 5.
    dzrat = 1.1
    dzmax = 2500.
    dtlong = 2.
    distim = 100.
    timmax = 7200.
    runtype = "HISTORY"
    igrdtyp = 1
    level = 3
   CCN = 100.e6
   filprf = 'trDct'
   hfilin = 'trDct.rst'
   savg_intvl = 3600.
   ssam_intvl = 60.
  corflg = .true.
  prndtl = -0.3333333
  iradtyp = 4
  strtim= 166.41 
  th00 = 289.  
  umean = -2.
   vmean = -4.
   cntlat=25.
     isfctyp=2  
      zrough = 0.001229
  dthcon = 0.001094
   drtcon = 0.001133
   case_name='trans'
   lsvarflg=.true.
    div=0.00000186
	  /
