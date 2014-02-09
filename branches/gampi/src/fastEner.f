      program fastEnergies
      use fastdef      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! MPI Init Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                 !
      ! rank   = the id of the processor                !
      ! nProc  = total number of processor              !
      ! nProc1 = total number of processor - 1          !
      !                                                 !
      call MPI_Init(ierr)                               ! 
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)    ! 
      call MPI_Comm_size(MPI_COMM_WORLD, nProc, ierr)   ! 
      nProc1 = nProc - 1                                !
      world = MPI_COMM_WORLD
      startTime = MPI_WTIME()                           ! Start timing...
      totAnaTime = 0.d0                                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc Program Input Section cccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	  
       !====================================!
       !  READ INPUT SECTION  (master only) !
       !====================================!
        if (rank .eq. master) then
        write(*,*)""
        tmpStr = all("=",60)
        write(*,*) tmpStr%data(1:tmpStr%length)	             
        
        write(*,*) " ENTER THE CELL FILE FILENAME: "
        read(*,'(a)') cellFile
        write(*,*) " USING CELL FILE NAMED |", 
     x             trim(cellFile), "|"
        write(*,*)""
        
        write(*,*) " ENTER SHAPE OF THE SYSTEM BOX (ORTHO/TRUNC): "
        read(*,'(a)') strBoxType        
        dStrBoxType = findBoxType(strBoxType)
         tmpStr = new(" USING A") + new(" ")
     x     + dStrBoxType + new(" BOX TYPE")
        write(*,*) tmpStr%data(1:tmpStr%length)                  
        write(*,*)""                

        tmpStr = new(" ENTER CALC. PARMS (GR,RDA,MDA,POL,COUL, .....")
        write(*,*) tmpStr%data(1:tmpStr%length)        
        tmpStr = new("        ... P2SH,QDR,IONE,DIONS,AVGS,PFILE,...")
        write(*,*) tmpStr%data(1:tmpStr%length)                
        tmpStr = new("        ... POTL,MSHEL,SOLUT,ICRG,PSF,EFS,....")        
        write(*,*) tmpStr%data(1:tmpStr%length)                
        tmpStr = new("        ... CLUSTER,COVAR,MTRAJ,PRVL,SELECT,..")        
        write(*,*) tmpStr%data(1:tmpStr%length)                
        tmpStr = new("        ... [PROT,SOLV]) ...............:")        
        write(*,*) tmpStr%data(1:tmpStr%length)                
        read(*,'(a)') strCalcParms        
        dStrEnerCalc = findCalcParm(strCalcParms)
         tmpStr = new(" CALCULATION PARMS: ") + dStrEnerCalc
        write(*,*) tmpStr%data(1:tmpStr%length)                  
        write(*,*)"" 

        if (doIndCrg) then
         write(*,*) " ENTER THE F-TRAJECTORY FILENAME: "
         read(*,'(a)') fTrajFile
         write(*,*) " USING F-TRAJ. FILE NAMED |", 
     x             trim(fTrajFile), "|"
         write(*,*)""   

         write(*,*) " ENTER POLARIZ. ATOM DATABASE (PA-DB) FILE: "
         read(*,'(a)') dVosPolDBFile
         write(*,*) " USING PA-DB FILE NAMED |", 
     x             trim(dVosPolDBFile), "|"
         write(*,*)""         
        endif
        
        if (doInduced .or. doAllEner) then
         write(*,*) " ENTER ATOMIC POLARIZ. DATABASE (AP-DB) FILE: "
         read(*,'(a)') atomPolFile
         write(*,*) " USING AP-DB FILE NAMED |", 
     x             trim(atomPolFile), "|"
         write(*,*)""
        endif ! End check for induced-only parameters
        
        write(*,*) " ENTER THE NFO PARAMETER FILENAME: "
        read(*,'(a)') paramsFile
        write(*,*) " USING PARAMETER FILE NAMED |", 
     x             trim(paramsFile), "|"
        write(*,*)""
        
        write(*,*) " ENTER THE AMBER TOPOLOGY OR PSF FILENAME: "
        read(*,'(a)') topFile
        write(*,*) " USING TOPOLOGY/PSF FILE NAMED |", 
     x              trim(topFile), "|"
        write(*,*)""

        write(*,*) " ENTER THE TRAJECTORY FILENAME: "
        read(*,'(a)') datFile
         tmpStr = new(" ") + 
     x    new("USING TRAJECTORY FILE NAMED |") + 
     x   new(datFile) + new("|")
        write(*,*) tmpStr%data(1:tmpStr%length)                  
        write(*,*)""

        write(*,*) " ENTER TRAJ. TYPE (BINPOS/MDCRD/GZMD/DCD): "
        read(*,'(a)') strFileType
        dStrFileType = findFileType(strFileType)
        tmpStr = new(" ") + new("USING THE") + new(" ")
     x   + dStrFileType + new(" FILE TYPE")     
        write(*,*) tmpStr%data(1:tmpStr%length)
        write(*,*)""

        write(*,*) " USING MPI-IO (YES/NO): "
        read(*,'(a)') strMPIIO
        dStrIOType = findIOType(strMPIIO)
        tmpStr = new(" ") + new("USING") + new(" ")
     x   + dStrIOType + new(" FOR READING FILES")     
        write(*,*) tmpStr%data(1:tmpStr%length)
        write(*,*)""

        write(*,*) " ENTER THE NUMBER OF CONFIGS IN MDCRD FILE: "
        read(*,*) nConf
c        tmpStr = new(" ") + new("USING") + new(" ")
c     x   + intToStr(nConf) + new(" CONFIGURATIONS")
c        write(*,*) tmpStr%data(1:tmpStr%length)    
        write(*,*) " USING ", nConf, " CONFIGURATIONS" 
        write(*,*)""
        
        if (doGrCalc) then        
         write(*,*) " ENTER THE ATOMIC INDEX FOR g(r): "
         read(*,*) iProbeAtom        
         tmpStr = new(" ") + new("USING INDEX") + new(" ")
     x   + intToStr(iProbeAtom) + new(" AS A PROBE ATOM")
         write(*,*) tmpStr%data(1:tmpStr%length)     
         write(*,*)""

         write(*,*) " ENTER THE g(r) FILE NAME: "
         read(*,*) strGrFile
         write(*,*) " USING A g(r) FILE NAMED |", 
     x              trim(strGrFile), "|"
         write(*,*)""

         write(*,*) " ENTER THE MIN(r-x) FILE NAME: "
         read(*,*) strMinFile
         write(*,*) " USING A MIN(r-x) FILE NAMED |", 
     x              trim(strMinFile), "|"
         write(*,*)""         
        endif ! End check for g(r) calculation

        if (doCBind) then        

         write(*,*) " ENTER THE AVERAGE TOTAL BINDING ENERGY: "
         read(*,*) epsilonAvg
         write(*,*) " USING AN AVERAGE TOTAL BINDING ENERGY |", 
     x              epsilonAvg, "|"
         write(*,*)""

         write(*,*) " ENTER THE AVERAGE DELTA BINDING ENERGY: "
         read(*,*) DepsilonAvg
         write(*,*) " USING AN AVERAGE DELTA BINDING ENERGY |", 
     x              DepsilonAvg, "|"
         write(*,*)""

         write(*,*) " ENTER THE SIMULATION TEMPERATURE: "
         read(*,*) simTemp
         write(*,*) " USING THE SIMULATION TEMPERATURE |", 
     x              simTemp, "|"
         write(*,*)""
         
        endif ! End check for binding energy cumulants
        
        if (doDIons) then        

         write(*,*) " ENTER THE SOLUTE IONS INPUT FILE NAME: "
         read(*,*) strDionInp
         write(*,*) " USING A SOLUTE IONS INPUT FILE NAMED |", 
     x              trim(strDionInp), "|"
         write(*,*)""

         write(*,*) " ENTER THE SOLUTE IONS OUTPUT FILE NAME: "
         read(*,*) strDionOut
         write(*,*) " USING A SOLUTE IONS OUTPUT FILE NAMED |", 
     x              trim(strDionOut), "|"
         write(*,*)""

        endif ! End check for distribution of ions calculation
        
        if(doRdaCalc) then
        write(*,*) " ENTER THE D-atom DATABASE FILENAME: "
        read(*,'(a)') DAtomDBFile
        write(*,*) " USING D-atom DATABASE FILE NAMED |",
     x               trim(DAtomDBFile), "|"
        write(*,*)""

        write(*,*) " ENTER THE A-atom DATABASE FILENAME: "
        read(*,'(a)') AAtomDBFile
        write(*,*) " USING A-atom DATABASE FILE NAMED |",
     x               trim(AAtomDBFile), "|"
        write(*,*)""        
        endif ! End check for R(D-A) Calculations


        if(doPFile) then
        write(*,*) " ENTER THE NAME IF THE OUTPUT PFILE: "
        read(*,'(a)') outDatFile
        write(*,*) " USING OUTPUT PFILE NAMED |",
     x               trim(outDatFile), "|"
        write(*,*)""
        endif

        if(doSelect) then
        write(*,*) " ENTER THE ATOM SELECTION STRING: "
        read(*,'(a)') strSelect
        write(*,*) " TAGGING THE FOLLOWING ATOM(S)",
     x            " NAMES |",trim(strSelect), "| "
        write(*,*)""
        endif     

        if(doCluster) then
        write(*,*) " ENTER THE CLUSTER RESIDUE'S NAME: "
        read(*,'(a)') strResSel
        write(*,*) " USING THE FOLLOWING RESIDUE FOR CLUSTER",
     x            " ANALYSIS |",trim(strResSel), "| "
        write(*,*)""
        endif

        if(doShellCalc) then
        write(*,*) " ENTER THE SHELL CUTOFF LENGTH (in Ang): "
        read(*,*) rBindCut
        write(*,*) " USING THE SHELL CUTOFF LENGTH |",
     x               rBindCut, "|"
        write(*,*)""
        endif       

        if(doPOTL) then
        write(*,*) " ENTER THE INDEX OF THE PROBE ATOM: "
        read(*,*) iAtomPoint
        write(*,*) " USING THE PROBE ATOM AT INDEX |",
     x               iAtomPoint, "|"
        write(*,*)""
        write(*,*) " ENTER THE NUMBER OF POINTS ON THE LINE: "
        read(*,*) nPts
        write(*,*) " USING THE FOLLOWING NUMBER OF POINTS |",
     x               nPts, "|"
        write(*,*)""
        endif       

        if(doEnerGap) then        
        write(*,*) " ENTER THE dVos ATOM DATABASE FILENAME: "
        read(*,'(a)') dVosDBFile
        write(*,*)" USING dVos DATABASE FILE NAMED |",
     x               trim(dVosDBFile), "|"
        write(*,*)""           
        endif
                
        write(*,*) " USE CUDA OR CPU?: "
        read(*,'(a)') strThreadType        
        dStrThreadType = findThreadType(strThreadType)
        tmpStr = new(" ") + new("USING") + new(" ")
     x   + dStrThreadType + new(" FOR THE CALCULATION")     
        write(*,*) tmpStr%data(1:tmpStr%length)
        write(*,*)""

        if(useCUDA) then
        write(*,*) " ENTER THE NUMBER OF GPUS PER NODE: "
        read(*,*) GPUsPerNode
        write(*,*) " THERE ARE |",
     x               GPUsPerNode, "| GPUs PER NODE"
ccccc REMOVED THIS BUT MAY NEED IT AGAIN cccccc
        !write(*,*)""
ccccc REMOVED THIS BUT MAY NEED IT AGAIN cccccc

        endif       

ccccc REMOVED THIS BUT MAY NEED IT AGAIN cccccc
ccccc   !tmpStr = all("=",60)
ccccc        !write(*,*) tmpStr%data(1:tmpStr%length)
ccccc        !write(*,*)""
ccccc REMOVED THIS BUT MAY NEED IT AGAIN cccccc
        
        !===================================!
        !  DATA INIT SECTION (master only)  !
        !===================================!
         write(*,*) ""		  
         tmpStr = all("=",55)
         write(*,*) tmpStr%data(1:tmpStr%length)	                     
         call readParams(trim(paramsFile))   ! Read internal .nfo file (maybe obsolete?)

         ! verify a selection exists for required routines, and if it doesnt
         call checkSelection() ! then force all atoms to to be                                  

         if (doPSFRead) then
         ! ---- charmm/xplor parameter reading secetion ---- !
          call fillCHARMMData()             ! Read in CHARMM/xplor psf file format
         ! ---- charmm/xplor parameter reading secetion ---- !
         else

         ! ---- amber parameter reading secetion ---- !

         !---- HERE HERE HERE HERE HERE ------ !
         !                                     !
         !     NEED TO ADD ATOM NAMED BASED    !
         !     SELECTION FOR AMBER STYLE       !
         !                                     !
         !---- HERE HERE HERE HERE HERE ------ !

          iTop = getUnit()
          open(unit=iTop, file=trim(topFile), status='old') 
          call readprm(iTop)              ! Read AMBER topology data, from AMBER code
          close(iTop)      
          tmpStr = all("=",55)
          write(*,*) tmpStr%data(1:tmpStr%length)
          write(*,*) ""		  
          call fillAmbData()                ! Read in AMBER toplogy/parameter data
         ! ---- amber parameter reading secetion ---- !
         endif

         if(doEnerGap) then
           call initDVosDbase(dVosDBFile)    ! Create database of atoms in Stokes shift calculation
         endif

         if(doInduced .or. doAllEner) then 
          call initAlphaDbase(atomPolFile) ! Creates database of all atomic polarizabilities
         endif
         if(doIndCrg) then
           call initDVosPolDbase(dVosPolDBFile, fTrajFile)   
         endif         
         if(doRdaCalc) then
          call initDADDbase(DAtomDBFile, AAtomDBFile) ! Creates the databases of all r(D-A) atoms
         endif

         if(.not. doMultiTraj) then
          call readHeader()   ! Read the header of the binpos/gz/mdcrd/dcd traj. file
         endif
         call readCellData()               ! Read cell data (alpha, beta, gamma, box lengths)
         
        endif ! End test for master     

        !==========================================!
        !  DATA INIT SECTION ( master and nodes )  !
        !==========================================!  

        call MPI_BARRIER(world,ierr)    

        !call setThreads()                  ! Tell all nodes how many threads to use  (obsolete?)
        call setMDParm()                    ! Setup MD parameters on all nodes 
        if(doMultiTraj) call setMulti()     ! Setup the multi-traj data and read the 1st traj header
        if(doSelect) call setSelection()    ! Set the indices on all nodes
        call createBoxL()                   ! Create MPI-derived type for boxL type
        call createAtom()                   ! Create MPI-derived type for basic_atom type           
        call setDatParm()                   ! Create MPI-derived type for dat_parm, and bcast it all
        if(doEnerGap) then
         call createDVosDbase()              ! Create MPI-derived type for dVosDbase        
         call sendDVosDbase()                ! Give a copy of dVos database to everyone
        endif
        if(doRdaCalc) then
         call createRDAR()                  ! Create MPI-derived type for the return R(D-A) data
         call createDADbase()               ! Create MPI-derived type for DADbase
         call sendDADbase()                 ! Send a copy of the DADbase to everyone
        endif       
        if(doMdaCalc) then
         call createMDAR()                  ! Create MPI-derived type for the return m(D-A) data        
        endif
        if(doInduced .or. doAllEner) then
          call sendAlphaDbase()             ! Give a copy of the atomic polariz. database to everyone
        endif
        call setQRM()                       ! Give a copy of the qrm database to everyone
        if(doIndCrg) then
         call sendDVosPolDbase()            ! Give a copy of polarizable atoms database to everyone
        endif 
        if(doCBind) then
         call setCBind()                    ! Send out the average delta/total binding energies
        endif
        if(doDIons) then
         call setDIons()                    ! Read and copy over all start indexes for the solute ion dist.
        endif
        if(doQdr) then
         call setAvgFile()                  ! Read and distribute <r> for all atoms
        endif
        if(doShell) then
         call sendAAcidList()               ! Create and distribute the AA-lookup table
        endif

        call allocConfs()                   ! Allocate nodeConf (all) and buffConfs (master only)

        if (doMPIIO .and. 
     x     (.not. doMultiTraj)) then        ! Do all MPI-IO initialization
         call readHeadMPIIO()    
        endif        
        if (doPFile) then
         call writeHeadMPIIO()              ! Write all headers for the PFILE
        endif        

        call outputCTState()                ! Used for simulation prep
        !call outputSolvMol()                ! Used for SolvMol prep
        call MPI_BARRIER(world,ierr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc End Program Input Section cccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc Trajectory Loop Section cccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	    
       if (rank .eq. master) then 
         write(*,*) ""
        tmpStr = all("_",45)
         write(*,*) tmpStr%data(1:tmpStr%length)	                     	                     
         write(*,3226) nConf
3226        format ("  ",i6,
     x              " FRAMES IN TRAJECTORY LOOP")
         print *, ""
         print *, ""
         print *, "    (F)rame     read (F/s)      ana (F/s)"
         print *, "   ---------  -------------  --------------"
       endif        

!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!
        if ((rank .eq. master) .and. doDebug) then 
         print *, "" 
         print *, "" 
         print *, "DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG" 
         print *, "DEBUG                                   DEBUG" 
         print *, "DEBUG      DEBUGGING MODE IS ON AND     DEBUG" 
         print *, "DEBUG   THE MAIN LOOP WILL BE SKIPPED   DEBUG" 
         print *, "DEBUG                                   DEBUG" 
         print *, "DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG" 
         print *, "" 
         print *, "" 
        endif
        if (.not. doDebug) then !!!!!! Only do the main loop if we're not debugging
!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!
!-----------------------------------------------!
! Loop over all configs, with a step            !
! equal to the # of nodes                       !
        do iConf = 1,nConf-nProc+1, nProc
!-----------------------------------------------!               

!       ++++++++++++ read section ++++++++++    !
         if (rank .eq. master) then             ! Start file i/o timing...
          readStartI = MPI_WTIME()      
         endif
         if (isBinary) then
          if (doMPIIO) then
           call readBINMPIIO()                  ! Read data using MPI-IO       
          else
           if (doGZRead) then                    
            call readGZMPI()                    ! Read and distribute the mdcrd.gz data
           elseif (doDCDRead) then              
            call readDCDMPI()                    ! Read and distribute the dcd data
           else ! then ... its a binpos file
            !call readBINMPI()                  ! Read and Distribute the binpos data
            call readBINMMPI()                  ! Read and distribute the modified binpos data
           endif
          endif 
           
         else
          call readCRDMPI()                     ! Read and Distribute the ascii mdcrd data
         endif

         if (rank .eq. master) then             ! End file i/o timing...
          readEndI = MPI_WTIME()   
          dReadI = real(readEndI - readStartI
     x                               ,kind=4)   
         endif
!       ++++++++++ end read section +++++++++   !

         if (rank .eq. master) then             ! Start analysis timing...
          anaStartI = MPI_WTIME()      
         endif
         
         call fillRecip()                       ! Calculate reciprocal vectors for distance reimaging 
         !call reimageMPI()                     ! Reimage on all nodes (obsolete?)
         
         if(doIndCrg) call setPolCharges()      ! Correct the charges on all polarizable atoms 
         if(doGrCalc) call calcGr(iProbeAtom)   ! Calculate a radial distribution function on the probe atom
         if(doDIons) call calcDIon()            ! Calculate the distribution of the solute ion distances
         if(doRdaCalc) call calcRDA()           ! Calculate the mass-weighted distance between donor and acceptors         
         if(doMdaCalc) call calcMDA()           ! Calculate the difference (charge transfer) dipole 
         if(doAvgS) call calcAvgStruct()        ! Calculate the average structure of the protein
         if(doCovar) call calcCovar()           ! Calculate the covariance matrix of a protein selection
         if(doQdr) call calcQdR()               ! Calculate <(qdr)^2> and consequently <(dr)^2>
         if(doShell) call calcEnerShel()        ! Calculate the 1st solvation shell energies/distributions
         if(doShellM) call calcShellM()         ! Calculate the 1st solvation shell dipole moments
         if(doShellP2) call calcShellP2()       ! Calculate the 1st solvation shell <P2>
         if(doCBind) call calcCumuBind()        ! Calculate high order binding energy cumulants
         if(doSolut) call calcEnerSolu()        ! Calculate the total solute-solvent interaction energy
         if(doCoul) call calcEnerCoul()         ! Calculate only coulomb (vert.) energies        
         if(doInduced) call calcEnerInd()       ! Calculate only induced (vert.) energies          
         if(doAllEner) call calcEnerAll()       ! Calculate both coulomb+induced (vert.) energies
         if(doPOTL) call calcPotL()             ! Calculate the potential along a line in the protein
         if(doDQEfield) call calcEfieldShel()   ! Calculate the e-field in water shell created by dq atoms
         if(doVPolCorr) call calcVpolShel()     ! Calculate the V_pol correction
         if(doQPolCorr) call calcQpolShel()     ! Calculate the Q_pol correction
         if(doProtM) call calcProtM()           ! Calculate the total dipole moment of the protein 
         if(doProtV) call calcProtV()           ! Calculate the total volume of the protein via MC
         if(doEFS) call calcEFSHist()           ! Calculate the E-field histogram of the shell
         if(doProtM .and. doShellM) then        ! Calculate the dot product of Mu.M
           call calcProtMuShellM()
         endif
         if(doPFile) call writeBINMPIIO()       ! Write output p-file
         if(doMultiTraj) call multiTraj()       ! Changes the traj pointer when using the multi-traj option 

         if (rank .eq. master) then             ! End analysis timing...
          anaEndI = MPI_WTIME()
          dAnaI = real(anaEndI - anaStartI
     x                               ,kind=4)                              
          totAnaTime = totAnaTime + anaEndI - anaStartI
!          totAnaTime = totAnaTime + dAnaI
         endif        
                  
         if (rank .eq. master) then             ! Periodically report progress to std. output

          thisFrame = iConf+nProc - 1
          dNFrames = real(thisFrame - prevFrame
     x                                ,kind=4)

          write(*,4227) thisFrame,
     x      (dNFrames/dReadI),(dNFrames/dAnaI)
4227        format ("    #",i7,4x,f9.4,4x,f13.6)
           prevFrame = thisFrame
         endif
!-----------------------------------------------!
        enddo
!-----------------------------------------------!

!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!
        endif !!!!!! END OF TEST FOR DEBUGGING
!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!
      ! Close all files used for reading/writing  
      if (doMPIIO) call MPI_FILE_CLOSE(idat, ierr) ! traj. (all nodes)       
      if (doPFile) call MPI_FILE_CLOSE(iOutDat, ierr)            
      if (rank .eq. master) then         
        tmpStr = all("_",45)
         write(*,*) tmpStr%data(1:tmpStr%length)	                        
        print *, ""
        
        ! close the trajectory file (unix-io)
        if (.not. doMPIIO) close(idat)  ! traj. (master only)
        
        ! close all open energy files
        if (doAllEner) then        
         close(iEne)
         close(iInd)
        elseif(doCoul) then
         close(iEne)
        elseif(doInduced) then
         close(iInd)               
        endif ! End check for ener types
        
        if (doSolut) then ! Close up the solute energy file
         close(iSolu)
        endif
        
        if (doShell) then ! Close the 1st shell energy file
         close(iShell)
        endif
        if(doVPolCorr) then
          close(iVpol)
        endif
        if(doQPolCorr) then
          close(iQpol)
        endif
        if (doCBind) then ! Close the binding energy cumulant file
         close(iCBind)
        endif
        if(doDQEfield) then ! Close the dq electric field file
         close(iEfShel)
        endif
        if (doProtM) then ! Close the protein dipole file
          close(iPMu)
        endif
        if (doProtV) then ! Close the protein volume file
          close(iPVol)
        endif
        if (doShellM .and. doProtM) then  ! Close the Mu.M file
          close(iPMuSh)
        endif
        if(doRdaCalc) then  ! Close up the R(D-A) file
         close(iRDA)
        endif      
        if(doMdaCalc) then  ! Close up the m(D-A) file
         close(iMDA)
        endif      
        if(doShellM) then   ! Close the shell dipole output file
         close(iShellM)
        endif
        if(doShellP2) then
         close(iShellP2)    ! Close the shell P2 output file
        endif
        if(doGZRead) then
         call closeGZ()     ! Close the mdcrd.gz file
        endif
        if(doDCDRead) then
         call closeDCD()     ! Close the dcd file
        endif
      endif   ! End check for master     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc End Trajectory Loop Section cccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	        

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc Post-Analysis Section cccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	        

      if(doGrCalc) call grFunc()      ! Converge partial histograms into total g(r)
      if(doDIons)  call dIonsFunc()   ! Converge partial ion distros into one
      if(doShell)  call bindFunc()    ! Converge partial binding energies into one
      if(doAvgS)   call avgSFunc()    ! Converge partial average positions into one
      if(doCovar)  call covarFunc()   ! Converge partial convariance matrices into one
      if(doQdr)    call qdrFunc()     ! Converge partial displacements into a single <(qdr)^2>
      if(doPOTL)   call collectPotL() ! Converge partial potential line data onto master
      if(doProtM)  call prMuFunc()    ! Converge partial P(cos(th),rcut) functions
      if(doEFS)    call efsFunc()     ! Converge partial histograms into total e-field histo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc End Post-Analysis Section cccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	        

      call MPI_BARRIER(world,ierr)     

      call deallocAll()           ! Deallocate all pointers 
      call reportTime()           ! Output the time statistics

      call MPI_BARRIER(world,ierr)     
        
      call MPI_Finalize(ierr)     ! Cleanup MPI processes on all nodes
       
	  stop 'PROGRAM ENDED NORMALLY'	
	  
      end program fastEnergies
