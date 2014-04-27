       module fastdef
        use commonFunctions

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! Amber header files !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !                                                 !	   
       INCLUDE "SIZE.h"                                  !
       INCLUDE "TOP.h"                                   !
       INCLUDE "MISC.h"                                  !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! MPI Params Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !                                                 !
       include "mpif.h"                                  !
       integer :: rank, nProc,ierr                       !
       integer :: proc, nProc1                           !   
       integer :: status(MPI_STATUS_SIZE)                !
       integer :: statusE(MPI_STATUS_SIZE,4)             !
       integer, parameter ::  master = 0                 !
       !integer, parameter ::  world = MPI_COMM_WORLD     !    
       integer ::  world != MPI_COMM_WORLD     !    
       integer            ::  CART_COMM                  !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! For AMBER mdcrd file conversion !!!!!!!!!!!!!!!!
      !                                                  !
       integer iTop, idat, iOutDat                       ! Topology, input, and output files
!       character(len=50) mdFile, pFile                   !
       character(len=80) titleStr                        ! Title string used in the mdcrd file      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! OpenMP Params Section !!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                  !
       integer :: threadID                               !
!       integer :: OMP_GET_THREAD_NUM                    !
       integer :: nThreads                               !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! PFA Derived Types !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       integer, parameter :: PDB_MAXRES=2000             ! The maximum allowed number of prot. residues       
       type dat_parm                                     ! Database that holds all important system info ...
           integer :: numRes                             ! ...
           integer :: numAtoms(PDB_MAXRES)               ! ...
           integer :: numProtAtoms                       ! ...
           integer :: numWaters                          ! ...
           integer :: numIons                            ! ...
       end type dat_parm                                 ! ...
	  type(dat_parm) params                              ! ... and storage for said database                          
                                                         !
       type basic_atom                                   ! Datatype for holding positional data
           real :: x                                     !
           real :: y                                     !
           real :: z                                     !       
       end type basic_atom                               !       
       type(basic_atom), pointer, 
     x      dimension(:) :: buffConfs                    ! All atomic positions for all configurations, i-i+nProc
       type(basic_atom), pointer, 
     x      dimension(:) :: nodeConf                     ! Single configuration of atomic positions on each node
                                                         !            
       type qrm                                          ! Database for holding 
          real :: q                                      ! charges (q) 
          real :: r                                      ! and radii (r)
          real :: m                                      ! and masses (m)
       end type qrm                                      ! 
       type(qrm), pointer, dimension(:) :: qrms          ! Storage for all system charges, radii, and masses
                                                         !
       type ari                                          ! Database for holding 
          integer :: iA                                  ! *A*tomic,
          integer :: iR                                  ! *R*esidue *I*ndex
       end type ari                                      ! 
       type(ari), pointer, dimension(:) :: ariD,ariA     ! Storage for ARI-database for r(D-A) calcs
                                                         !
       type rdaR                                         ! Used to store return data for R(D-A) 
          real :: xD                                     ! calculations
          real :: yD                                     ! ... position vectors
          real :: zD                                     ! ... and
          real :: xA                                     ! ... distance
          real :: yA                                     ! ... 
          real :: zA                                     ! ... 
          real :: rDA                                    ! ...
       end type rdaR                                     ! 
       type(rdaR), pointer, dimension(:) :: rdaRs,rdaRn  ! Storage for a trajectory of returns, and 1 node       
                                                         !
       type mdaR                                         ! Used to store return data for m(D-A) 
          real :: xM                                     ! calculations
          real :: yM                                     ! ... difference dipole vector
          real :: zM                                     ! ... and
          real :: mDA                                    ! ... magnitude
       end type mdaR                                     ! 
       type(mdaR), pointer, dimension(:) :: mdaRs,mdaRn  ! Storage for a trajectory of returns, and 1 node       
                                                         !
                                                         !                                                         
        type dVos_data                                   ! Database filetype for storing
          integer :: atom                                ! atoms involved in the Stokes shift
          integer :: residue                             ! correlation function
          real    :: q                                   !
        end type dVos_data                               !
        type(dVos_data), pointer,
     x      dimension(:) :: dVosDbase                    ! RC-atom Database storage
        type(dVos_data), pointer,
     x      dimension(:) :: dVosPolDbase                 ! Polarizable RC-atom Database storage     
                                                         !
                                                         !
       integer, parameter :: MAX_RESIDUES=20000          ! Maximum number of selected residues
       integer, parameter :: MAX_ATOMS=35                ! Maximum number of atoms per selected residue
       type res_list_data                                ! Derived type for cluster data residue info
          integer :: resid                               ! PDB ID of selected residue
          integer :: numAtomsPerRes                      ! number of atoms per selected residue
          integer :: atomid(MAX_ATOMS)                   ! nodeConf indices of selected atoms
       end type                                          !
                                                         !
                                                         ! ********* MAKE A POINTER ************
       type(res_list_data) :: resList(MAX_RESIDUES)      ! Storage for all selected resList atoms and resids
       type(res_list_data) :: fullResList(MAX_RESIDUES)  ! Storage for all selected residues and corresponding atoms
                                                         ! ********* MAKE A POINTER ************
                                                         !
       type res_atoms                                    ! Derived type for a single cluster residues postions
           integer          :: nResAtoms                 ! number of atoms in the residue
           real             :: atomMass(MAX_ATOMS)       ! masses of res_list atoms
           type(basic_atom) :: atomPos(MAX_ATOMS)        ! x,y,z positions of res_list atom positions
       end type     

       type cluster_postions                             ! Derived type to hold positional data of each cluster
         !integer          :: nRes                        ! number of cluster residues
         type(basic_atom) :: aggCOM                      ! center of mass of the aggregate
         !type(res_atoms)  :: resPos(MAX_RESIDUES)        ! Position storage for all residues in a cluster
       end type      
                                                         !
                                                         !
       integer, parameter :: MAX_PER_CLUSTER=550         ! Maximum number of residues in the cluster
       type agg_list_data                                ! Derived type for cluster data residue info
          integer :: nAgg                                ! number of residues in the aggregation
          integer :: resid(MAX_PER_CLUSTER)                  ! resids of the aggregate
       end type                                          !
                                                         !
       type agg_parm_data                                ! Derived type for cluster parameter
          real :: parm                                   ! (real valued) parameter of the cluster
       end type                                          !

                                                         !
       type boxL                                         ! Used to store boxes 
          real :: x                                      ! in a convenient way to
          real :: y                                      ! pass for every time step
          real :: z                                      ! in in non-constant V simulations
       end type boxL                                     ! 
       type(boxL), pointer, dimension(:) :: boxLs,boxLn  ! Storage for a trajectory of box lengths, and 1 node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! MPI Derived Types !!!!!!!!!!!!!!!!!!!!!!!      
        integer :: PARM_TYPE                      ! MPI-derived type
        integer :: DVOSDB_TYPE                    ! (cont.)
        integer :: ATOM_TYPE                      ! (cont.)        
        integer :: CONF_TYPE                      ! (cont.)  
        integer :: QRM_TYPE                       ! (cont.)        
        integer :: BOXL_TYPE                      ! (cont.)                
        integer :: DADB_TYPE                      ! (cont.)
        integer :: RDAR_TYPE                      ! (cont.)        
        integer :: MDAR_TYPE                      ! (cont.)
        integer :: RESL_TYPE                      ! (cont.)
        integer :: AGG_TYPE                       ! (cont.)
        integer :: AGG_PARM_TYPE                  ! (cont.)
                                                  !
        !integer(kind=MPI_OFFSET_KIND) :: offset   ! Used for file offsets using MPI-IO
        integer*8 :: offset   ! Used for file offsets using MPI-IO
        !integer(kind=MPI_OFFSET_KIND) :: atomOff  ! Used for file offsets using MPI-IO
        integer*8 :: atomOff  ! Used for file offsets using MPI-IO
        !integer(kind=MPI_OFFSET_KIND) :: intOff   ! Used for file offsets using MPI-IO
        integer*8 :: intOff   ! Used for file offsets using MPI-IO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!! Various integer variables !!!!!!!!!        
        integer, parameter :: iEne=56           ! Holds handle for energy file
        integer, parameter :: iSolu=57          ! Holds handle for solute energy file                
        integer, parameter :: iXOfile = 511     ! Holds the handle for the Min(r(X-O)) traj.
        integer, parameter :: iInd=154          ! Holds the handle for the induced energy traj.       
        integer, parameter :: iRDA=155          ! Holds the handle for the r(D-A) output.        
        integer, parameter :: iMDA=156          ! Holds the handle for the m(D-A) output.                
        integer, parameter :: iDIinp=157        ! Holds the handle for the solute dist. input        
        integer, parameter :: iDIout=158        ! Holds the handle for the solute dist. output        
        integer, parameter :: iShell=159        ! Holds the handle for the 1st shell output
        integer, parameter :: iCBind=160        ! Holds the handle for the binding energy cumulant output
        integer, parameter :: iShellM=161       ! Holds the handle for the shell dipole moment output
        integer, parameter :: iShellP2=162      ! Holds the handle for the shell P2 output
        integer, parameter :: iPotL=163         ! Holds the handle for the potential lines (POTL) output       
        integer, parameter :: iPMu=164          ! Holds the handle for the protein dipole output               
        integer, parameter :: iPMuSh=165        ! Holds the handle for the Mu.M output
        integer, parameter :: iEfShel=166       ! Holds the handle water shell efield data
        integer, parameter :: iVpol=167         ! Holds the handle for the Vpol correction data
        integer, parameter :: iQpol=167         ! Holds the handle for the Qpol correction data
        integer, parameter :: iPVol=168         ! Holds the handle for the protein dipole output
        integer, parameter :: iClusterOut=169   ! Holds the handle for the cluster analysis output
        integer, parameter :: iHoleVol=170      ! Holds the handle for the hole volume output
        integer, parameter :: iHoleCOM=171      ! Holds the handle for the hole C.O.M. input
        integer, parameter :: iEnerAll=172      ! Holds the handle for the file with both coul/ind energies
                                                !
        integer            :: nWaters           ! The total number of waters
                                                !        
        integer            :: iWatStart         ! Waters starting atomic index
        integer            :: iWatEnd           ! Waters ending atomic index
        integer            :: iProtStart        ! The proteins starting atomic index
        integer            :: iProtEnd          ! The proteins ending atomic index
        integer            :: iIonStart         ! The ions starting atomic index
        integer            :: iIonEnd           ! The ions ending atomic index
        integer            :: iLipStart         ! The lipids starting atomic index
        integer            :: iLipEnd           ! The lipds ending atomic index
        integer            :: iOthStart         ! Other atoms starting atomic index
        integer            :: iOthEnd           ! Other atoms ending atomic index
                                                !        
    	integer, parameter :: maxSect = 10      ! maximum number of starts/ends
        integer,dimension(maxSect) :: arrWatStart  ! array of starting indices for water
        integer,dimension(maxSect) :: arrWatEnd    ! array of ending indices for water
        integer,dimension(maxSect) :: arrProtStart ! array of starting indices for the protein
        integer,dimension(maxSect) :: arrProtEnd   ! array of ending indices for the protein
        integer,dimension(maxSect) :: arrIonStart  ! array of starting indices for the ions
        integer,dimension(maxSect) :: arrIonEnd    ! array of ending indices for the ions
        integer,dimension(maxSect) :: arrLipStart  ! array of starting indices for the lipids
        integer,dimension(maxSect) :: arrLipEnd    ! array of ending indices for the lipids
        integer,dimension(maxSect) :: arrOthStart  ! array of starting indices for all other atoms
        integer,dimension(maxSect) :: arrOthEnd    ! array of ending indices for all other atoms
        integer            ::  nWatStEnd        ! number of section startings/endings for water
        integer            ::  nProtStEnd       ! number of section startings/endings for the proteins
        integer            ::  nIonStEnd        ! number of section startings/endings for ions
        integer            ::  nLipStEnd        ! number of section startings/endings for lipids
        integer            ::  nOthStEnd        ! number of section startings/endings for other atoms
                                                !        
        integer            :: iProbeAtom        ! Atom for which a radial distrib. func. is calculated
                                                !        
        integer iConf,nConf                     ! Configuration counter, total number
        integer ios                             ! I/O status 
        integer numConfigs                      ! Stores total number of configurations
        integer iNumDB                          ! Total number of dVos database elements
        integer iNumPolDB                       ! Total number of polarizable atoms
        integer iNumFDB                         ! Total number of f-Parameters in the trajectory        
        integer iNumRdaD                        ! Total number of Donor elements in the r(D-A) calcs
        integer iNumRdaA                        ! Total number of Acceptor elements in the r(D-A) calcs        
        integer iSoluteIons                     ! Total number of ionic solute resiudes
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! Various string variables !!!!!!!!!        
        character(len=80) datFile               ! traj file or name of the multi-traj file (input)
        character(len=50) cellFile              ! David's cell filetype (input)        
        character(len=50) fTrajFile             ! Trajectory of f-data
        character(len=50) paramsFile, topFile   ! AMBER param/topology files (input)
        character(len=50) dVosDBFile            ! Database of atoms involved in SS correlation
        character(len=50) dVosPolDBFile         ! Database of atoms polarizable atoms
        character(len=20) strThreadType         ! Type of threads (CPU-only or CUDA)
        character(len=20) strFileType           ! Type of energy calculation (BINPOS or MDCRD)        
        character(len=20) strMPIIO              ! Type of IO-TYPE (MPI-IO or UNIX)
        character(len=20) strBoxType            ! Type of system box (ORTHORHOMBIC OR TRUNCATED)
        character(len=80) strCalcParms          ! Name a string with the energy requests
        character(len=20) strGrFile             ! File name of the g(r) output file
        character(len=20) strBindFile           ! File name of the binding output file (total energy)
        character(len=20) strBindDeFile         ! File name of the binding output file (delta energy)        
        character(len=20) strDionInp            ! File name of the distrib. of solute ions input file        
        character(len=20) strDionOut            ! File name of the distrib. of solute ions output file                
        character(len=20) strMinFile            ! File name of the MIN(r-x) output file    
        character(len=20) DAtomDBFile           ! File name of the donor atoms for the r(D-A) calcs
        character(len=20) AAtomDBFile           ! File name of the donor atoms for the r(D-A) calcs   
        character(len=80) outDatFile            ! Output for parallel-file creation
        character(len=100) strSelect            ! The atomic selection string (from std input)
        character(len=20) strResSel             ! The residue selection string (from std input)
        character(len=30) strClusterOut         ! The output file name for the cluster analysis
        character(len=30) strHoleCOM            ! The (input) file name containing the hole C.O.M. vectors
        character(len=40) exclrFile             ! File name for holding all per-atom exclusions
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
			
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!! Various dstring variables !!!!!!!!!
        type(dstring) :: dFileName              ! dstring of the filename
      	type(dstring) :: tmpStr                 ! used for writing to the screen, files
        type(dstring) :: dStrThreadType           ! used for writing to the screen the Ener. calc. Type
        type(dstring) :: dStrFileType           ! used for writing to the screen the file type
        type(dstring) :: dStrIOType             ! used for writing to the screen the type of IO (MPI-IO or unix)
        type(dstring) :: dStrBoxType            ! used for writing to the screen the shape of the box
        type(dstring) :: dStrEnerCalc           ! used for writing to the screen the requested energy calculations
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! Program decision variables !!!!!!!  
        logical :: sendBoxL                     ! Boolean to tell if all nodes need box lengths
        logical :: isBinary                     ! Boolean to tell if using mdcrd or binpos trajectory
        logical :: doMPIIO                      ! Boolean to tell if using MPI-IO (true) or unix-io (false)
        logical :: doPFile                      ! Boolean to tell if we are creating a parallel-io traj.
        logical :: doProt                       ! Boolean to tell if calculating protein energies
        logical :: doSolv                       ! Boolean to tell if calculating solvent energies
        logical :: doLipid                      ! Boolean to tell if calculating lipid energies
        logical :: doOthers                     ! Boolean to tell if calculating others energies
        logical :: doSolut                      ! Boolean to tell if calculating total solute energies        
        logical :: doShell                      ! Boolean to tell if calculating total 1st shell energies        
        logical :: doShellM                     ! Boolean to tell if calculating total 1st shell dipole moments
        logical :: doShellP2                    ! Boolean to tell if calculating total 1st shell <P2(r.e)>
        logical :: doShellCalc                  ! Boolean to tell if doing ANY shell calculation        
        logical :: doIons                       ! Boolean to tell if calculating ion solvation energies        
        logical :: doDIons                      ! Boolean to tell if calculating solute ion distributions        
        logical :: doCoul                       ! Boolean to tell if only coulomb (protien + solvent) energies
        logical :: doIndCrg                     ! Boolean to tell if need to calculate energies using induced charges        
        logical :: doEnerGap                    ! Boolean to tell if doing an energy gap-style calculation
        logical :: doInduced                    ! Boolean to tell if only calculating induced energies
        logical :: doAllEner                    ! Boolean to tell if calculating both induced + coulomb energies        
        logical :: doGrCalc                     ! Boolean to tell if we are doing a g(r) calculation
        logical :: doRdaCalc                    ! Boolean to tell if we are doing a R(D-A) calculation        
        logical :: doMdaCalc                    ! Boolean to tell if we are doing a m(D-A) calculation                
        logical :: doAvgS                       ! Boolean to tell if we are calculating the avg. structure
        logical :: doCovar                      ! Boolean to tell if we are calculating a proteins covariance matrix
        logical :: doQdr                        ! Boolean to tell if we are calculating <(qdr)^2>
        logical :: doCBind                      ! Boolean to tell if we calculating binding energy cumulants
        logical :: isTruncated                  ! Boolean to tell if box is orthorhombic (true) or truncated (false) 
        logical :: doGZRead                     ! Boolean to tell if we are reading a gzipped mdcrd file
        logical :: doDCDRead                    ! Boolean to tell if we are reading a dcd trajectory file
        logical :: doPSFRead                    ! Boolean to tell if we are reading a psf parameter file
        logical :: doPOTL                       ! Boolean to tell if we are calculating the potential along lines
        logical :: doProtM                      ! Boolean to tell if we are calculating the protein dipole
        logical :: doProtV                      ! Boolean to tell if we are calculating the proteins volume
        logical :: doHoleV                      ! Boolean to tell if we are calculating a cavity volume
        logical :: doGPUInit                    ! Boolean to tell if we need to initialize the GPU ids
        logical :: doDQEfield                   ! Boolean to tell if calculating the electric field due to the vertical energy (probe) charges
        logical :: doEFS                        ! Boolean to tell if calculating the e-field histogram along OH bonds
        logical :: doVPolCorr                   ! Boolean to tell if calculating the Vpol correction
        logical :: doQPolCorr                   ! Boolean to tell if calculating the Qpol correction
        logical :: doSelect                     ! Boolean to tell if setting a selection of atoms
        logical :: doCluster                    ! Boolean to tell if doing cluster analysis (selection required)
        logical :: doMultiTraj                  ! Boolean to tell if we are reading multiple traj files
        logical :: doDebug                      ! Boolean to tell if we are skipping the main loop for debugging
        logical :: doGyr                        ! Boolean to tell if we are doing the radius of gyration calculation
        logical :: doInert                      ! Boolean to tell if we are doing the moment of intertia calculation
        logical :: doExclr                      ! Boolean to tell if we are doing the exclusion calculation
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!! Timing variables !!!!!!!!!!!!!!!!!!!!
        double precision   :: startTime, endTime       ! 
        double precision   :: anaStartI, anaEndI       !
        double precision   :: readStartI, readEndI     !
        real               :: dAnaI,dReadI             !  
        double precision   :: totAnaTime               ! 
        real :: diffTime                               !
        real :: diagTime                               !
        integer :: thisFrame                           !
        integer :: prevFrame                           !
        real    :: dNFrames                            !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! Polarizability variables !!!!!!!!!!!!!!!!!!!!  
        character(len=50) :: atomPolFile                   ! Filename of the atomic polariz.       
        real, pointer, dimension(:) :: alphaDbase          ! Database of atomic polarizabilities
        real, pointer, dimension(:) :: origDeltaQ          ! Normal charge difference for pol. atoms
        real, pointer, dimension(:) :: origQ               ! Normal atomic charges for pol. atoms        
        real, pointer, dimension(:) :: fFacTraj            ! Entire trajectory of f-parameter data
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! g(r) variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	integer, parameter :: maxBin = 50000               ! maximum number of histogram bins
        real,dimension(maxBin) :: hist, hist_i, gr000      ! histogram storage
        real, pointer, dimension(:) :: lowRXO              ! stores partial r(X-O) values
        real, parameter ::  delr = 0.01                    ! dr
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! multi-traj variables !!!!!!!!!!!!!!!!!!!!!!!!
        integer, parameter :: MTMAX = 100                  ! max number of multi-traj files
        type(dstring)      :: mTrajNames(MTMAX)            ! array containing all the traj file names
        integer            :: mTrajLengths(MTMAX)          ! array containing all the file names
        integer            :: nMultiTraj                   ! number of multi-traj names and lengths
        integer            :: idxTraj                      ! working index of our multi-traj array
        integer            :: nFrameTraj                   ! count of current traj length in frames
        integer            :: lastTraj                     ! lastTraj (should be prevFram)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! e-field histogram variables !!!!!!!!!!!!!!!!!
    	integer, parameter :: maxBinEFS = 2000             ! maximum number of histogram bins
        real,dimension(maxBinEFS) :: histEFS, histEFS_i    ! histogram storage
        real, parameter    :: efsShift = 1.0
        real, parameter ::  delEFS = 0.001                 ! EFS bin size
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! P(cos_wat/prot(theta),rCut) variables  !!!!!!!
    	integer, parameter :: maxBinCos = 1000              ! maximum number of histogram bins
        real,dimension(-maxBinCos:maxBinCos) :: 
     x                   histCos,histCos_i,PCos ! histogram storage
        real,dimension(-maxBinCos:maxBinCos) :: histCosCUDA            ! per timestep histo storage on the GPU
c        real, parameter :: dCosShift = 2.0                  ! a variable to shift the storage to a positive value
        real, parameter ::  dCosT = 0.001                  !
        type(basic_atom)  :: extraPmu, extraShellM       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!! CUDA variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer :: GPUsPerNode                               ! The number of GPUs per parallel node
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! distribution of ions variables !!!!!!!!!!!!!!!!
    	integer, parameter :: maxBinIon = 50000              ! maximum number of histogram bins
        integer, pointer, dimension(:) :: starts, resNum     ! index of 1st atom, and residue index of solu. ion dist.
        real,dimension(maxBinIon) :: dI_norm, distI_i, distI ! histogram storage
        real, parameter ::  drI = 0.001                      ! dr
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! Potential line variables !!!!!!!!!!!!!!!!!!!!!!!
          integer :: iAtomPoint                               ! Index of the atom in the potential line of interest
          integer :: nPts                                     ! The number of points in the potential line calc.
          real, pointer, dimension(:) :: arrRMIN              ! Node arrays for storing the min-distance to h2o...
          real, pointer, dimension(:) :: arrDist              ! the distance from the proteins COM
          real, pointer, dimension(:) :: arrPotL              ! potential along a line
          real, pointer, dimension(:) :: arrPotL2             ! squared-potential along a line
          type(basic_atom), pointer, dimension(:) :: arrPts   ! storage of the points along a line
          real, pointer, dimension(:) :: rmin                 !
          real, pointer, dimension(:) :: distPL               !
          real, pointer, dimension(:) :: potl                 !
          real, pointer, dimension(:) :: potl2                !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!! vertical energy charge electric field vars.  !!!!!
       real, pointer, dimension(:) :: efield_dVTraj           ! trajectory of the magnitude of e-field
                                                              ! from all delta q charges
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!! Q/Vpol correction energy vars.                !!!!!
       real, pointer, dimension(:) :: dVpolTraj            ! trajectory of the Vpol(t)
       real, pointer, dimension(:) :: dQpolTraj            ! trajectory of the Vpol(t)
                                                              ! from all nodes
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!! MC volume vars.                          !!!!!
       real, pointer, dimension(:) :: mcVolProtTraj   ! trajectory of V(t) for the protein
       real, pointer, dimension(:) :: mcVolHoleTraj   ! trajectory of V(t) for the hole
       real, pointer, dimension(:) :: mcVolWatTraj   ! trajectory of V(t) for the water shell
       real, pointer, dimension(:) :: ranX                ! cuda random numbers for the x-dimension
       real, pointer, dimension(:) :: ranY                ! cuda random numbers for the y-dimension
       real, pointer, dimension(:) :: ranZ                ! cuda random numbers for the z-dimension
       real, pointer, dimension(:) :: mcEpsProtTraj       ! trajectory of our convergence factor for the protein vol
       real, pointer, dimension(:) :: mcEpsHoleTraj       ! trajectory of our convergence factor for the protein vol
       real, pointer, dimension(:) :: mcEpsWatTraj        ! trajectory of our convergence factor for the water vol
       real                        :: mcVolProt           ! volume of the protein at a given step
       real                        :: mcVolHole           ! volume of the hole at a given step
       real                        :: mcVolWat           ! volume of the water shell at a given step
       real                        :: mcEpsProt           ! convergence factor for the protein at a given step
       real                        :: mcEpsHole           ! convergence factor for the hole at a given step
       real                        :: mcEpsWat           ! convergence factor for the water at a given step
       integer, parameter    :: nAttempts=1000000          ! max number of MC steps
       integer                     :: iAttempt            ! the current MC step
       real                  :: protRadius                ! the radius of the protein
       real                  :: holeRadius                ! the radius of the hole probe
       real                  :: totalRadius               ! the radius of the protein and shell cutoff
       real                  :: totalR2                   ! the radius of the protein and shell cutoff
       type(basic_atom), pointer, dimension(:) :: holeCOMTraj   ! trajectory of the hole C.O.M. vectors
       type(basic_atom)                        :: holeCOM   ! the current hole C.O.M. vector
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! Binding energy variables !!!!!!!!!!!!!!!!!!!!!!!
    	integer, parameter :: maxBinBind = 10000              ! maximum number of histogram bins
    	integer, parameter :: maxBinBindDe = 100000           ! maximum number of histogram bins    	
        real*8, parameter :: enerShift = 10.d0                ! energy shift for making an easy histogram
           !!! Histograms !!                                  !        
        real,dimension(-maxBinBind/2:maxBinBind/2) :: 
     x                       dB_norm, distB_i, distB          ! histo. storage: 1st line-total bindE
        real,dimension(-maxBinBindDe/2:maxBinBindDe/2) ::      
     x                       dBde_norm, distBde_i, distBde    ! histo. 2nd-deltaEpsilon
        real, parameter ::  drB   = 0.001d0                   ! dr, for total binding
        real, parameter ::  drBde = 0.0001d0                  ! dr, for deltaEpsilon
           !!! Total/Delta Binding Averages !!                !
        real            :: epsilonAvg, DepsilonAvg            ! real*4 averages (used for message passing)
        real*8          :: r8epsAvg, r8epsDAvg                ! real*8 averages (used in the calculations)
           !!! Instantaneous Fluctuation moments !!           !        
        real*8          :: deps, dDeps                        ! total/delta binding energy fluctuations
        real*8          :: deps2, dDeps2                      ! squared total/delta binding energy fluctuations
        real*8          :: deps3                              ! cubed total binding energy fluctuations 
           !!! Binding energy Cumulants (nodes) !!            !        
        real*8          :: mdepsdDeps                         ! average/delta binding energies 1st moment
        real*8          :: mdeps2dDeps                        ! mixed average/delta binding energies 2nd moment
        real*8          :: mdeps3dDeps                        ! mixed average/delta binding energies 3rd moment
        real*8          :: mdDeps2                            ! delta binding energy 2nd moment  
        real*8          :: mdDeps2deps                        ! mixed delta/average binding energy 2nd moment
        real*8          :: mdDeps2deps2                       ! mixed delta/average binding energy 2nd moment
           !!! Binding energy Cumulants (master) !!           !        
        real            :: mdepsdDepsM                        ! average/delta binding energies 1st moment
        real            :: mdeps2dDepsM                       ! mixed average/delta binding energies 2nd moment
        real            :: mdeps3dDepsM                       ! mixed average/delta binding energies 3rd moment
        real            :: mdDeps2M                           ! delta binding energy 2nd moment  
        real            :: mdDeps2depsM                       ! mixed delta/average binding energy 2nd moment
        real            :: mdDeps2deps2M                      ! mixed delta/average binding energy 2nd moment  
           !!! Cumulant trajectories (on master) !!!          ! 
        real, pointer, dimension(:) :: depsdDepsTraj          ! trajectory of <(deps)(d(Deps))>
        real, pointer, dimension(:) :: deps2dDepsTraj         ! trajectory of <(deps)^2(d(Deps))>        
        real, pointer, dimension(:) :: deps3dDepsTraj         ! trajectory of <(deps)^3(d(Deps))>
        real, pointer, dimension(:) :: dDeps2Traj             ! trajectory of <(d(Deps))^2>
        real, pointer, dimension(:) :: dDeps2depsTraj         ! trajectory of <(d(Deps))^2(deps)>
        real, pointer, dimension(:) :: dDeps2deps2Traj        ! trajectory of <(d(Deps))^2(deps)^2>
           !!! Shell number trajectories (on master) !!!      !         
        integer, pointer, dimension(:) :: nShellTraj          ! trajectory of the number of 1st shell waters
        integer, pointer, dimension(:) :: nShellLowTraj       ! trajectory of the number of low energy 1st shell waters        
        integer, pointer, dimension(:) :: nShellHighTraj      ! trajectory of the number of high energy 1st shell waters 
        integer :: nShellLow, nShellHigh                      ! low/high energy amino acids (node storage)
           !!! data per amino acid pointers (nodes) !!!       !
        real, pointer, dimension(:) :: ebindTotAA             ! total binding energy
        real, pointer, dimension(:) :: ebindDelAA             ! delta binding energy
        real, pointer, dimension(:) :: nShellLowAA            ! low energy amino acids
        real, pointer, dimension(:) :: nShellHighAA           ! high energy amino acids        
           !!! data per amino acid pointers (master) !!!      !
        real, pointer, dimension(:) :: nHighAA(:)             ! Number of high...
        real, pointer, dimension(:) :: nLowAA(:)              ! ...and low energy amino acids
        real, pointer, dimension(:) :: eTotAA(:)              ! Total and ..
        real, pointer, dimension(:) :: eDelAA(:)              ! ...and partial binding energies per AA           
           !!! amino acid lookup table !!!                    !
        integer, pointer, dimension(:) :: aminoAcidList       ! used for fast look-ups of the amino acid number  
        real*8, parameter :: bindCutoff=-0.011758762169680d0  ! cutoff between low and high energy bound waters
        real              :: rBindCut                         ! Cutoff for all shell calculations
        real*8, parameter :: beta300=557.d0                   ! This is 1/(kBT) at 300K in q^2/Ang (lab units)
        real*8            :: betaT                            ! This is 1/(kBT) at simTemp in q^2/Ang (lab units)       
        real              :: simTemp                          ! The simulation temperature
        real*8            :: expBde                           ! This holds exp(-beta*deps)
        real*8            :: lamYi                            ! This holds the node's individual lambda_Yi
        real*8            :: Yi                               ! This holds the node's individual Y_0i        
        real, pointer, dimension(:) :: lamYiTraj(:)           ! lambda_Yi trajectory of node data
        real, pointer, dimension(:) :: YiTraj(:)              ! <Y_0i> trajectory of node data        
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! Energy variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real, pointer, dimension(:) :: dv0sTraj            ! trajectory of vertical energies (solv.)
        real, pointer, dimension(:) :: dv0iTraj            ! trajectory of vertical energies (ions)        
        real, pointer, dimension(:) :: dv0pTraj            ! trajectory of vertical energies (prot.)       
        real, pointer, dimension(:) :: dv0lTraj            ! trajectory of vertical energies (lipids)             
        real, pointer, dimension(:) :: dv0oTraj            ! trajectory of vertical energies (others)             
        real, pointer, dimension(:) :: dv0shdeTraj         ! trajectory of 1st shell delta energies
        real, pointer, dimension(:) :: dv0shTraj           ! trajectory of 1st shell energies
        real, pointer, dimension(:) :: dv0bshTraj          ! trajectory of 1st shell binding energies
        real, pointer, dimension(:) :: dv0bdeshTraj        ! trajectory of 1st shell binding delta energies
        real, pointer, dimension(:) :: dbe2Traj            ! trajectory of d(eps)^2
        real, pointer, dimension(:) :: dv0cbeTraj          ! trajectory of full cell binding energies
        real, pointer, dimension(:) :: dv0inTraj           ! trajectory of vertical energies (ind.)
        real, pointer, dimension(:) :: dv0inSTraj          ! trajectory of vertical energies (ind., solv.)
        real, pointer, dimension(:) :: dv0inITraj          ! trajectory of vertical energies (ind., ions.)        
        real, pointer, dimension(:) :: dv0inPTraj          ! trajectory of vertical energies (ind., prot.)
        real, pointer, dimension(:) :: dv0inLTraj          ! trajectory of vertical energies (ind., lip.)
        real, pointer, dimension(:) :: dv0inOTraj          ! trajectory of vertical energies (ind., Others)
        real*8 :: v0s_i, v0s,avgV0s                        ! v0s for each ith water, and for each timestep.
        real*8 :: dv0s_i, dv0s                             ! partial and total vertical (solvent) energies
        real*8 :: dv0p_i, dv0p                             ! partial and total vertical (protein) energies
        real*8 :: dv0i_i, dv0i                             ! partial and total vertical (ions) energies
        real*8 :: dv0l_i, dv0l                             ! partial and total vertical (lipids) energies      
        real*8 :: dv0o_i, dv0o                             ! partial and total vertical (other atom) energies   
        real*8 :: dv0sh_i, dv0sh                           ! partial and total 1st shell energies
        real*8 :: dv0shde_i, dv0shde                       ! partial and total 1st shell delta energies
        real*8 :: dv0bsh_i, dv0bsh                         ! partial and total average binding energy
        real*8 :: dv0bshde_i, dv0bshde                     ! partial and total average delta binding energy        
        
        real*8 :: dv0cbe_i, dv0cbe                         ! partial and total full cell binding energy
        real*8 :: dv0inS                                   ! partial and total vertical (solvent, ind.) energies
        real*8 :: dv0inI                                   ! partial and total vertical (ions, ind.) energies        
        real*8 :: dv0inP                                   ! partial and total vertical (protein, ind.  energies
        real*8 :: dv0inL                                   ! partial and total vertical (lipids, ind.) energies
        real*8 :: dv0inO                                   ! partial and total vertical (others, ind.) energies
        real*8 :: dv0in                                    ! the total induced energy        
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!! Position and dqr variables !!!!!!!!!!!!!!!!!!
        real, pointer, dimension(:) :: xAvg                ! partial average positions (x-axis)
        real, pointer, dimension(:) :: yAvg                ! partial average positions (y-axis)
        real, pointer, dimension(:) :: zAvg                ! partial average positions (z-axis)        
        real, pointer, dimension(:) :: xAvgTot             ! total average positions (x-axis)
        real, pointer, dimension(:) :: yAvgTot             ! total average positions (y-axis)
        real, pointer, dimension(:) :: zAvgTot             ! total average positions (z-axis)        
        real, pointer, dimension(:) :: drxAvg              ! partial average fluctuations (x-axis)        
        real, pointer, dimension(:) :: dryAvg              ! partial average fluctuations (y-axis) 
        real, pointer, dimension(:) :: drzAvg              ! partial average fluctuations (z-axis) 
        real, pointer, dimension(:) :: qdrsq               ! partial <(qdr)^2> (but really <(dr)^2>)
        real, pointer, dimension(:) :: drxAvgTot           ! total average fluctuations (x-axis)        
        real, pointer, dimension(:) :: dryAvgTot           ! total average fluctuations (y-axis) 
        real, pointer, dimension(:) :: drzAvgTot           ! total average fluctuations (z-axis) 
        real, pointer, dimension(:) :: qdrsqTot            ! total <(qdr)^2> 
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! Reimaging parameters !!!!!!!!!!!!!!!!!!!!!!!!!!
       !                                                 !   
       double precision ucell(3,3), recip(3,3)           !
       double precision alpha, gamma                     ! beta is defined in TOP.h
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! Ewald parameters !!! (deprecated?) !!!!!!!!!!!!
       !                                                 !       
         real               :: rCut, rCut2               ! cutoff and cutoff^2
         real               :: bX,bY,bZ                  ! x,y,z boxlengths
         real               :: volCell                   ! volume of simulation cell         
         real               :: kappa                     ! convergence parameter    
         integer, parameter :: maxk = 1000               ! maximum number of k-vectors
         real               :: kvec(maxk)                ! k-vectors
         integer, parameter :: kmax = 5                  ! maximum k-vectors in a given direction 
         integer, parameter :: ksqmax = 27               ! maximum k*k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! Selection variables                           !
         integer                :: nSelect               ! the total number of selected atoms
         integer, pointer,
     x             dimension(:) :: selectIndices         ! the atomic indices in "nodeConf" space
         type(basic_atom), pointer,
     x             dimension(:) :: selectAtoms           ! the selected positional data
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! Exclusion analysis variables                  !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         real, pointer,
     x             dimension(:) :: exclrDist             ! the exclusion distance for each system particle
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! Cluster analysis variables                    !
          real                  :: rClusterCut           ! The cutoff for a pair to be called a neighbor
          integer               :: nClusterRes           ! The number of residues selected, (MAX_RES)
          integer               :: nAgg                  ! The number of total aggregates
          integer,parameter     :: MAX_AGG=2000          ! The max size of the aggregates list (aggList)
          integer, pointer, 
     x             dimension(:) :: nAggTraj              ! the trajectory of total aggregates 
          integer, pointer, 
     x             dimension(:) :: neighborCount         ! the local number of neighbors
          integer, pointer, 
     x             dimension(:,:) :: neighborList        ! the full per node array of neighbor matrix
          logical, pointer, 
     x             dimension(:) :: checkList             

          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggEccen              ! eccentricity (shape parm) of each aggregate
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggEccenTraj            ! nProc worth of eccenticities
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggA2                ! a^2 shape parm of each aggregate
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggA2Traj            ! nProc worth of a^2
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggB2                ! b^2 shape parm of each aggregate
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggB2Traj            ! nProc worth of b^2
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggC2                ! c^2 shape parm of each aggregate
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggC2Traj            ! nProc worth of c^2
!!!!!!!!!!!!!!!!!!!!!!!
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggGyr                ! radius of gyration of each aggregate
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggGyrTraj            ! nProc worth of radius of gyrations
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggR1                ! high eig. val. of gyr tensor
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggR1Traj            ! nProc worth of high gyr tensor
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggR2                ! mid eig. val. of gyr tensor
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggR2Traj            ! nProc worth of mid gyr tensor
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggR3                ! low eig. val. of gyr tensor
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggR3Traj            ! nProc worth of low gyr tensor
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggAsp                ! asphericity of gyr tensor
          type(agg_parm_data), pointer, 
     x             dimension(:) :: aggAspTraj            ! nProc worth of asphericity
!!!!!!!!!!!!!!!!!!!!!!!
          type(agg_list_data), pointer, 
     x             dimension(:) :: aggList               ! list of aggregates on each node
          type(agg_list_data), pointer, 
     x             dimension(:) :: aggListTraj           ! full nProc worth of aggregation lists
          type(cluster_postions), pointer,
     x             dimension(:) :: aggPositions          ! list of aggregate positions on each node
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! PCA parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       type(basic_atom), pointer, 
     x             dimension(:) :: proteinCOM            ! storage for the protein in COM
                                                         ! coordinates
       type(basic_atom), pointer, 
     x             dimension(:) :: proteinCRYS           ! storage for the protein crystal
                                                         ! structure
       type(basic_atom), pointer, 
     x             dimension(:) :: proteinROT            ! storage for the rotated protein
                                                         ! structure
       real, pointer, 
     x             dimension(:,:) :: covarNode           ! storage for the nodes own mass-weighted 
                                                         ! covariance matrix
       real, pointer, 
     x             dimension(:,:) :: covarMat            ! storage for the full mass-weighted
                                                         ! covariance matrix
       real, pointer, 
     x             dimension(:) :: vecAvg                ! storage for the mass-weighted
                                                         ! average in vector format (for easy looping)
       real, pointer, 
     x             dimension(:) :: vecConf               ! storage for the mass-weighted
                                                         ! prot. configuration in vector format
                                                         !
       integer                    :: iPCADoF             ! number of degrees of freedom
        type(basic_atom) :: prMuCntr, prAtom             ! prMuCntr is the center at every... These basic_atom types
        type(basic_atom) :: prMuNode, prCrys             ! need to save state through-out
        type(basic_atom) :: prMuHat                      ! the life ofthe program
                                                         !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
       !!! CUDA variables                                !
          logical            :: useCUDA                  ! logical to use CUDA device
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! Operator overloading for position vectors !!!!!
       !                                                 !
        interface operator(*)                            ! 
         module procedure mult_p                         !
        end interface                                    !
                                                         !
        interface operator(+)                            !
         module procedure add_p                          !
      	end interface                                    !
                                                         !		 	
        interface operator(-)                            !
         module procedure diff_p                         !
      	end interface                                    !
                                                         !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
       
         contains

		function add_p(atom1,atom2)
c    *******************************************************************
c    ** find the sum of two position vectors                          **
c    **                                                               **
c    ** principal variables:                                          **
c    ** type(basic_atom)   atom1, atom2  positions of two atoms       **
c    **                                                               **
c    ******************************************************************* 				  
		
          implicit none
          
          ! function type
          type(basic_atom) :: add_p
			
          ! Args
          type(basic_atom), intent(in) :: atom1, atom2

c    *******************************************************************

          add_p%x = atom1%x+atom2%x
          add_p%y = atom1%y+atom2%y
          add_p%z = atom1%z+atom2%z
          
          return
        end function add_p
         
         
		function diff_p(atom1,atom2)
c    *******************************************************************
c    ** find the distance between two position vectors                **
c    **                                                               **
c    ** principal variables:                                          **
c    ** type(basic_atom)   atom1, atom2  positions of two atoms       **
c    **                                                               **
c    ******************************************************************* 				  
		
          implicit none
          
          ! function type
          type(basic_atom) :: diff_p
			
          ! Args
          type(basic_atom), intent(in) :: atom1, atom2

c    *******************************************************************

          diff_p%x = atom1%x-atom2%x
          diff_p%y = atom1%y-atom2%y
          diff_p%z = atom1%z-atom2%z
          
          return
        end function diff_p

		function mult_p(atom1,atom2)
c    *******************************************************************
c    ** find the dot product of two position vectors                  **
c    **                                                               **
c    ** principal variables:                                          **
c    ** type(basic_atom)   atom1, atom2  positions of two atoms       **
c    **                                                               **
c    ******************************************************************* 				  
		
          implicit none
          
          ! function type
          real :: mult_p
			
          ! Args
          type(basic_atom), intent(in) :: atom1, atom2

c    *******************************************************************

          mult_p = atom1%x*atom2%x + 
     x             atom1%y*atom2%y + 
     x     	       atom1%z*atom2%z
          
          return
      end function mult_p         
 
      FUNCTION ran1(idum)
c *******************************************************
c * this the RAN3 routine from the FORRAN recipes book **
c *******************************************************

      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ

      REAL ran1,FAC

      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj

      ran1  =  mj*FAC

      return
      END function ran1

      function findThreadType(strType)
c    *******************************************************************
c    ** determine if using Ewald or Direct energy calculations        **
c    **                                                               **
c    ** principal variables:                                          **
c    ** type(dstring) findThreadType  return dstring data               **
c    ** character(len=*)   strType  string data: "Ewald" or "Direct"  **
c    **                                                               **
c    *******************************************************************		
      implicit none
          
      ! function def
      type(dstring) findThreadType
      
      ! Args
      character(len=*) strType
      integer :: iPos
      
      ! Locals
c    *******************************************************************
      
      
      if (inStr("CUDA",trim(strType),iPos)) then       
       findThreadType = new("CUDA THREADS")   
       useCUDA = .true.
cc#ifdef ENABLE_CUDA
cc       write(*,*) " WARNING !!! : CUDA THREADS",
cc     x            " ARE NOT ENABLE ON TYIS SYSTEM."
cc       write(*,*) " YOU WILL NEED TO RECOMPILE FOR",
cc     x            " FOR CUDA THREADS TO USE THEM"
cc       write(*,*) "   ONLY THE CPU WILL BE USED !!!"
cc       useCUDA = .false.
cc#endif
      else       
       findThreadType = new("ONLY THE CPU")       
       useCUDA = .false.       
      endif
      
      
      return
      end function findThreadType

      function findCalcParm(strCalc)
c    *******************************************************************
c    ** determine the type of energy calculations to be performed.    **
c    **                                                               **
c    ** principal variables:                                          **
c    ** type(dstring) findCalcParm  return dstring data               **
c    ** character(len=*)   strCalc  string data                       **
c    **                                                               **
c    **                                                               **
c    *******************************************************************		
      implicit none
          
      ! function def
      type(dstring) findCalcParm
      
      ! Args
      character(len=*) strCalc
      integer :: iPos
      
      ! Locals
      integer, parameter :: iArr = 10
      integer            :: iNumEner
      integer            :: jj
      logical            :: foundUnk, doEner      
      type(dstring)      :: dArray(iArr), dTmp, dUnk
      type(dstring)      :: dCoulP, dCoulS, dPol, dCoul
      type(dstring)      :: dIndCrg, dGr, dRda, dMda
      type(dstring)      :: dIons, dDIons, dSolut
      type(dstring)      :: dShell, dAvgS, dQdr, dCBind
      type(dstring)      :: dPFILE, dShelM, dShelP2
      type(dstring)      :: dPOTL, dPRMU,dEfldDV,dDVpolC
      type(dstring)      :: dDQpolC, dPSF, dCovar, dEFS
      type(dstring)      :: dProtV, dSel, dMulti, dClust
      type(dstring)      :: dDebug, dGyr, dInert, dOther
      type(dstring)      :: dHoleV, dLipid, dExclR
c    *******************************************************************
      
      ! Set up dstrings used as the test cases
      dCoulP  = new('PROT')
      dCoulS  = new('SOLV')  
      dOther  = new('OTHER')    
      dSolut  = new('SOLUT')      
      dShell  = new('SHELL')
      dShelM  = new('MSHEL')
      dShelP2 = new('P2SH')
      dCBind  = new('CBIND')
      dIons   = new('IONE')
      dLipid  = new('LIPID')
      dDIons  = new('DION')      
      dPol    = new('POL')
      dCoul   = new('COUL')
      dIndCrg = new('ICRG')
      dGr     = new('GR')
      dRda    = new('RDA')
      dMda    = new('MDA')
      dAvgS   = new('AVGS')
      dQdr    = new('QDR')
      dPFILE  = new('PFILE')
      dPOTL   = new('POTL')
      dProtV  = new('PRVL')
      dHoleV  = new('HOLE')
      dPRMU   = new('PRMU')
      dEfldDV = new('EFLD')
      dDVpolC  = new('VPLC')
      dDQpolC = new('QPLC')
      dPSF    = new('PSF')
      dCovar  = new('COVAR')
      dEFS    = new('EFS')
      dSel    = new('SELECT')
      dMulti  = new('MTRAJ')
      dClust  = new('CLUSTER')
      dDebug  = new('DEBUG')
      dGyr    = new('GYR')
      dInert  = new('INERT')
      dExclR  = new('EXCLR')

      ! init all output dstrings
      findCalcParm = new('')
      dUnk         = new('')
      
      ! Init all logical variables
      doProt      = .false.
      doSolv      = .false.
      doOthers    = .false.
      doLipid     = .false.
      doSolut     = .false.
      doShell     = .false.
      doCBind     = .false.
      doIons      = .false.     
      doCoul      = .false.
      doInduced   = .false.
      doIndCrg    = .false.
      doGrCalc    = .false.
      doAllEner   = .false.
      doEnerGap   = .false.
      doRdaCalc   = .false.
      doMdaCalc   = .false.
      doDIons     = .false.       
      foundUnk    = .false.
      doEner      = .false.
      doAvgS      = .false.
      doCovar     = .false.
      doQdr       = .false.
      doPFile     = .false.
      doShellM    = .false.
      doShellP2   = .false.
      doShellCalc = .false.   
      doPOTL      = .false.
      doProtM     = .false.
      doProtV     = .false.
      doHoleV     = .false.
      doDQEfield  = .false.
      doVPolCorr  = .false.
      doQPolCorr  = .false.
      doPSFRead   = .false.
      doEFS       = .false.
      doSelect    = .false.
      doCluster   = .false.
      doMultiTraj = .false.
      doDebug     = .false.
      doGyr       = .false.
      doInert     = .false.
      doExclr     = .false.

      dTmp = new('  INPUT STRING: |') + new(trim(strCalc))
      dTmp = dTmp + new('|')
      call dStrWrite(6,dTmp)       
      call createStrArr(strCalc,dArray,iArr,iNumEner)
      
      ! Loop over all words to see what kind of energies we get
      do jj = 1, iNumEner
       if (indStr(dArray(jj),dCoul,iPos)) then ! Found Coulomb calcs
        doCoul    = .true.      
        doEner    = .true.
        doEnerGap = .true.
        findCalcParm = findCalcParm + new(' COUL')
       elseif (indStr(dArray(jj),dPol,iPos)) then ! Found Induced calcs
        doInduced = .true.       
        doEner    = .true.       
        doEnerGap = .true. 
        findCalcParm = findCalcParm + new(' POL')             
       elseif (indStr(dArray(jj),dCoulS,iPos)) then ! Found solvent calcs
        doSolv    = .true.
        doEner    = .true.     
        doEnerGap = .true.   
        findCalcParm = findCalcParm + new(' SOLV')    
       elseif (indStr(dArray(jj),dSolut,iPos)) then ! Found solute calcs
        doSolut    = .true.
        doEner     = .true.        
        doEnerGap = .true.
        findCalcParm = findCalcParm + new(' SOLUT')          
       elseif (indStr(dArray(jj),dCoulP,iPos)) then ! Found protein calcs
        doProt    = .true.
        doEner    = .true.    
        doEnerGap = .true.    
        findCalcParm = findCalcParm + new(' PROT')   
       elseif (indStr(dArray(jj),dLipid,iPos)) then ! Found lipid calcs
        doLipid   = .true.
        doEner    = .true.    
        doEnerGap = .true.    
        findCalcParm = findCalcParm + new(' LIPID')   
       elseif (indStr(dArray(jj),dOther,iPos)) then ! Found other calcs
        doOthers  = .true.
        doEner    = .true.    
        doEnerGap = .true.    
        findCalcParm = findCalcParm + new(' OTHER')   
       elseif (indStr(dArray(jj),dPOTL,iPos)) then ! Found potential line calcs
        doPOTL    = .true.
        doEnerGap = .true.
        findCalcParm = findCalcParm + new(' POTL')           
       elseif (indStr(dArray(jj),dShell,iPos)) then ! Found shell calcs
        doShell   = .true.
        doShellCalc = .true.
        doEnerGap = .true.
        findCalcParm = findCalcParm + new(' SHELL') 
       elseif (indStr(dArray(jj),dEfldDV,iPos)) then ! Found E-field dV calcs
        doDQEfield   = .true.
        doShellCalc = .true.
        doEnerGap = .true.
        findCalcParm = findCalcParm + new(' EFLD')   
       elseif (indStr(dArray(jj),dEFS,iPos)) then ! Found Efield shell calcs
        doEFS       = .true.
        doShellCalc = .true.
        doEnerGap = .true.  ! NOT SURE IF THIS IS TRUE
        findCalcParm = findCalcParm + new(' EFS')   
       elseif (indStr(dArray(jj),dDVpolC,iPos)) then ! Found Vpol correction calcs
        doVPolCorr   = .true.
        doShellCalc = .true.
        doEnerGap = .true.
        findCalcParm = findCalcParm + new(' VPLC')   
       elseif (indStr(dArray(jj),dDQpolC,iPos)) then ! Found Qpol correction calcs
        doQPolCorr   = .true.
        doShellCalc = .true.
        doEner    = .false.
        doEnerGap = .true.        
        findCalcParm = findCalcParm + new(' QPLC')   
       elseif (indStr(dArray(jj),dShelM,iPos)) then ! Found shell dipole calcs
        doShell = .true.
        doShellM  = .true.
        doShellCalc = .true.        
        findCalcParm = findCalcParm + new(' MSHEL') 
       elseif (indStr(dArray(jj),dSel,iPos)) then ! Found selection
        doSelect = .true.        
        findCalcParm = findCalcParm 
     x               + new(' SELECT') 
       elseif (indStr(dArray(jj),dClust,iPos)) then ! Found cluster analysis
        doCluster = .true.        
        findCalcParm = findCalcParm 
     x               + new(' CLUSTER') 
       elseif (indStr(dArray(jj),dGyr,iPos)) then ! Found radias of gyration analysis
        doGyr = .true.        
        findCalcParm = findCalcParm 
     x               + new(' GYR') 
       elseif (indStr(dArray(jj),dInert,iPos)) then ! Found moment of intertia analysis
        doInert = .true.        
        findCalcParm = findCalcParm 
     x               + new(' INERT')
       elseif (indStr(dArray(jj),dExclR,iPos)) then ! Found exclusion radius calculation
        doExclr = .true.
        doEnerGap = .true.
        findCalcParm = findCalcParm
     x               + new(' EXCLR')
       elseif (indStr(dArray(jj),dDebug,iPos)) then ! Found debugging flag
        doDebug = .true.        
        findCalcParm = findCalcParm 
     x               + new(' DEBUG') 
       elseif (indStr(dArray(jj),dMulti,iPos)) then ! Found multiple trajectory analysis
        doMultiTraj = .true.        
        findCalcParm = findCalcParm 
     x               + new(' MTRAJ')
       elseif (indStr(dArray(jj),dPRMU,iPos)) then ! Found protein/shell dipole calcs
        doProtM  = .true.
        doShellCalc = .true.
        findCalcParm = findCalcParm + new(' PRMU')        
       elseif (indStr(dArray(jj),dProtV,iPos)) then ! Found protein volume calcs
        doProtV  = .true.
        doShellCalc = .true.
        findCalcParm = findCalcParm + new(' PRVL')      
       elseif (indStr(dArray(jj),dHoleV,iPos)) then ! Found hole volume calcs
        doHoleV  = .true.
        findCalcParm = findCalcParm + new(' HOLE')      
       elseif (indStr(dArray(jj),dShelP2,iPos)) then ! Found P2 calcs
        doShellP2  = .true.
        doShellCalc = .true.       
        doEner     = .false.        
        findCalcParm = findCalcParm + new(' P2SH')            
       elseif (indStr(dArray(jj),dCBind,iPos)) then ! Found binding energy cumulant calcs
        doCBind   = .true.
        doShellCalc = .true.        
        doEner    = .true.        
        findCalcParm = findCalcParm + new(' CBIND')         
       elseif (indStr(dArray(jj),dIons,iPos)) then ! Found ion (solvent) calcs
        doIons    = .true.
        doEner    = .true.       
        doEnerGap = .true. 
        findCalcParm = findCalcParm + new(' IONE')               
       elseif (indStr(dArray(jj),dIndCrg,iPos)) then ! Found polarizable charges
        doIndCrg = .true.
        doEner    = .true.       
        doEnerGap = .true. 
        findCalcParm = findCalcParm + new(' ICRG')        
       elseif (indStr(dArray(jj),dGr,iPos)) then ! Found g(r) calculation
        doGrCalc = .true.
        findCalcParm = findCalcParm + new(' GR')    
       elseif (indStr(dArray(jj),dAvgS,iPos)) then ! Found avg. structure calc.
        doAvgS = .true.
        findCalcParm = findCalcParm + new(' AVGS')   
       elseif (indStr(dArray(jj),dCovar,iPos)) then ! Found covariance matrix calc.
        doCovar = .true.
        findCalcParm = findCalcParm + new(' COVAR')   
       elseif (indStr(dArray(jj),dQdr,iPos)) then ! Found <(qdr)^2> calc. -- should remove this option
        doQdr = .true.
        doEnerGap = .true.
        findCalcParm = findCalcParm + new(' QDR')   
       elseif (indStr(dArray(jj),dRda,iPos)) then ! Found r(D-A) calculation
        doRdaCalc = .true.
        findCalcParm = findCalcParm + new(' RDA')         
       elseif (indStr(dArray(jj),dMda,iPos)) then ! Found m(D-A) calculation
        doMdaCalc = .true.
        findCalcParm = findCalcParm + new(' MDA')     
       elseif (indStr(dArray(jj),dDIons,iPos)) then ! Found dist. of solute ions
        doDIons    = .true.
        findCalcParm = findCalcParm + new(' DION')         
       elseif (indStr(dArray(jj),dPFILE,iPos)) then ! Found parallel file output
        doPFile    = .true.
        findCalcParm = findCalcParm + new(' PFILE')        
       elseif (indStr(dArray(jj),dPSF,iPos)) then ! Found psf file input
        doPSFRead    = .true.
        findCalcParm = findCalcParm + new(' PSF')       
       else
        if(dArray(jj)%length .gt. 0) then ! Check for a null string
        foundUnk = .true.
        dUnk = dUnk + new(' |') + dArray(jj)
        dUnk = dUnk + new('|')
        endif
       endif       
     
      enddo ! End loop over iNumEner

      ! Check for default input (prot/solv)
      if((.not. doSolv) .and. (.not. doProt)
     x   .and. (.not. doSolut) 
     x   .and. (.not. doShell) .and. doEner) then
       doSolv = .true.
       doProt = .true.
       findCalcParm = findCalcParm + 
     x                new(' FORCED(SOLV+PROT)')
      endif

      ! Check for forced selection 
      if((.not. doSelect) .and. doCluster) then
        doSelect = .true.        
        findCalcParm = findCalcParm 
     x               + new(' FORCED(SELECT)') 
      endif      
      ! Check for default input (coul/pol)
c       if(doEner .and. (.not. doShell)
c     x    .and. (.not. doSolut)) then     
c        doAllEner = .true.
c        findCalcParm = findCalcParm + new(' COUL+POL')       
c       endif

      ! If we are doing both coul+pol anyway, do them in only 1 routine
       if(doCoul .and. doInduced) then
         doAllEner = .true.
         doCoul = .false.    ! Turn off idividual calcs ...
         doInduced = .false. ! ... same here.
       endif

      ! If we are using exclusions, we must do all energetics
      if (doExclr) then
       !doCoul = .false.    ! Turn off idividual calcs ...
       !doInduced = .false. ! ... same here.
       doAllEner = .false. ! just in case...

       doSolv = .true.
       doProt = .true.
       findCalcParm = findCalcParm + 
     x                new(' FORCED(SOLV+PROT)')

      endif

      ! Report any unknown parameters
      if(foundUnk) then
       dTmp = new('  --> UNKNOWN PARMS:') + dUnk
       call dStrWrite(6,dTmp)      
      endif
      
      return
      end function findCalcParm

      function findBoxType(strBox)
c    *******************************************************************
c    ** determine if the box shape is orthorhombic or a truncated     **
c    ** octahedron                                                    **
c    **                                                               **
c    ** principal variables:                                          **
c    ** type(dstring) findBoxType  return dstring data                **
c    ** character(len=*)   strBox  string data: "Ortho"               **
c    **                                      or "Truncated"           **
c    **                                                               **
c    *******************************************************************		
      implicit none
          
      ! function def
      type(dstring) findBoxType
      
      ! Args
      character(len=*) strBox
      integer :: iPos
      
      ! Locals
c    *******************************************************************
      
      
      if (inStr("ORTHO",trim(strBox),iPos)) then       
       findBoxType = new("ORTHORHOMBIC")   
       isTruncated = .false.
      else              
       findBoxType = new("TRUNCATED")       
       isTruncated = .true.       
      endif      
      
      return
      end function findBoxType
      
      function findFileType(strFile)
c    *******************************************************************
c    ** determine if using Binpos or formatted MDCRD files            **
c    **                                                               **
c    ** principal variables:                                          **
c    ** type(dstring) findFileType  return dstring data               **
c    ** character(len=*)   strType  string data: "BINPOS" or "MDCRD"  **
c    **                                                               **
c    *******************************************************************		
      implicit none
          
      ! function def
      type(dstring) findFileType
      
      ! Args
      character(len=*) strFile

      
      ! Locals
      integer :: iPos      
c    *******************************************************************
      
      
      if (inStr("BINPOS",trim(strFile),iPos)) then       
       findFileType = new("BINPOS")   
       isBinary = .true.
       doGZRead = .false.       
       doDCDRead = .false.
      elseif (inStr("GZMD",trim(strFile),iPos)) then       
       findFileType = new("GZMD")   
       doGZRead = .true.
       isBinary = .true.  
       doDCDRead = .false.
      elseif (inStr("MDCRD",trim(strFile),iPos)) then
       findFileType = new("MDCRD")
       isBinary = .false.       
       doGZRead = .false.              
       doDCDRead = .false.
      else
       findFileType = new("DCD")
       isBinary = .true.       
       doGZRead = .false.              
       doDCDRead = .true.
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! TO DO: MOVE doPSFRead setting to another
!!!        more psf-specific location
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      return
      end function findFileType

      function findIOType(strIO)
c    *******************************************************************
c    ** determine if reading trajectory using MPI or UNIX IO          **
c    **                                                               **
c    ** principal variables:                                          **
c    ** type(dstring) findIOType   return dstring data                **
c    ** character(len=*)   strIO  string data: "YES" is the answer    **
c    **                                                               **
c    *******************************************************************		
      implicit none
          
      ! function def
      type(dstring) findIOType
      
      ! Args
      character(len=*) strIO

      
      ! Locals
      integer :: iPos      
c    *******************************************************************
      
      
      if (inStr("YES",trim(strIO),iPos)) then       
       findIOType = new("MPI-IO")   
       doMPIIO = .true.
      else       
       findIOType = new("UNIX-IO")
       doMPIIO = .false.       
      endif
      
      
      return
      end function findIOType

      subroutine reportTime()
c    *******************************************************************
c    ** reports the time usage to standard output                     **
c    **                                                               **
c    *******************************************************************           
       implicit none
 
      real :: aTime, rTime
c   ********************************************************************
      endTime = MPI_WTIME()
      if (rank .eq. master) then
       write(*,*)""
       write(*,*) "================================================"     
      diffTime = real((endTime - startTime),kind=4)
      aTime = real(totAnaTime,kind=4)
      rTime = diffTime - aTime
      if(doCovar) rTime = rTime - diagTime
c        tmpStr =  new(" Total run time = ")
c     x             + new(" ") 
c     x   + realToStr(diffTime)
c     x             + new(" seconds.")
c        print *,tmpStr%data(1:tmpStr%length)
        print *, " Total run time = ", diffTime,
     x           " seconds"		
c        tmpStr =  new(" Total Analysis time = ")
c     x             + new(" ") 
c     x   + realToStr(aTime)
c     x             + new(" seconds.")
c        print *,tmpStr%data(1:tmpStr%length)		
        print *, " Total analysis time = ", aTime,
     x           " seconds"		

c        tmpStr =  new(" Total Read time = ")
c     x             + new(" ") 
c     x   + realToStr(rTime)
c     x             + new(" seconds.")
c        print *,tmpStr%data(1:tmpStr%length)		
        print *, " Total read time = ", rTime,
     x           " seconds"		
      write(*,*) "================================================"      
      write(*,*)""
      endif
      
      call MPI_BARRIER(world,ierr)            
      
       return
      end subroutine reportTime
      
         subroutine readParams(parmFile)
c    *******************************************************************
c    ** reads internal dat_parm datatype from a file into the memory  **
c    **                                                               **
c    *******************************************************************           
          implicit none
																								
	    	! Args
          character(len=*) parmFile
          !type(dat_parm) params
																								
	       ! Locals
          integer :: iUnitD
          integer myi, myNum
c   ********************************************************************
!          print *, " readParams: 1st line = ", nConf
          
          iUnitD = 854!getUnit()																								
          
!          print *, " readParams: set id = ", nConf
          open(unit=iUnitD, file=parmFile, status='old')
!          print *, " readParams: open = ", nConf          
          read(iUnitD,*) params%numRes 
          read(iUnitD,*) params%numIons
          read(iUnitD,*) params%numWaters
          !pdb_numres = params%numRes
          !pdb_num_H2O = params%numWaters
          !pdb_numNaPlus = params%numIons
!          print *, " readParams: 1 reads = ", nConf          									
																								
          params%numProtAtoms = 0
          do myi = 1, params%numRes
!          print *, " readParams: myi", myi, nConf          
          read(iUnitD,*) params%numAtoms(myi) !myNum 
           params%numProtAtoms = params%numAtoms(myi) + 
     x      params%numProtAtoms
          !pdb_aminos(i)%numAtoms = params%numAtoms(i)
          enddo
          close(iUnitD)
!          print *, " readParams: nconf = ", nConf          
          return
         end subroutine
																							
																							
        subroutine fillAmbData()
c    *******************************************************************
c    ** reads AMBER topology file data to the parameter datatype      **
c    **                                                               **
c    *******************************************************************               
        
	      implicit none
	     
	  INCLUDE "SIZE.h"
          INCLUDE "TOP.h"
          INCLUDE "MISC.h"   
		
		  integer iFirstSod, iNumProt
		  integer nres1
          logical foundSodiums
c    *******************************************************************
          foundSodiums = .false.
		  ! First, find the first sodium
		  do i = 1, natom
		  
		  if(s_eqi_str(igraph(i)(1:3),'Na+') .or. 
     x       s_eqi_str(igraph(i)(1:3),'Cl-')) then
     
           foundSodiums = .true.
           iFirstSod = i
           exit             
           endif
		  
		  enddo

          if(foundSodiums) then		  
		  
          iProtEnd = iFirstSod - 1
          iIonStart = iFirstSod		  
		  do i = iFirstSod, natom
		    if(.not. s_eqi_str(igraph(i)(1:3),'Na+') .and. 
     x       .not. s_eqi_str(igraph(i)(1:3),'Cl-')) exit
		  enddo
		  
           iWatStart = i
           iIonEnd = i - 1           
		  else ! did not found any counter ions
		   ! First, find the first sodium
		  
		   iProtEnd = natom - params%numWaters*3
		   iWatStart = iProtEnd + 1
		   iFirstSod = iProtEnd
		   i = iProtEnd
                   iIonStart = 0
                   iIonEnd = 0
		  endif


           ! This assumes the protein has the starting index
           iProtStart = 1

           ! This assumes the end of the water atoms is the last index
           iWatEnd = natom
           write(*,*) ""

           !tmpStr = all("=",50)
           tmpStr = new("=========================")
     x            + new("=========================")          
           write(*,*) tmpStr%data(1:tmpStr%length)
	             
	   write(*,*) "     Total number of protein atoms: ", 
     x                iProtEnd
	   write(*,*) "      Total number of counter ions: ",
     x                 i - iFirstSod
	   write(*,*) "            Total number of waters: ",
     x               params%numWaters     
	   write(*,*) "       Total number of water atoms: ",
     x               params%numWaters*3
	   write(*,*) "   Total number of atoms in system: ",
     x               natom
           write(*,*) tmpStr%data(1:tmpStr%length)	                      
	   write(*,*) ""
		   
	 nres1 = nres-1
          
         allocate(qrms(natom)) ! allocate qrms on master, must allocate on nodes later
         
         do i = 1, natom         
          qrms(i)%q = chrg(i)/18.222
          qrms(i)%r = rborn(i)
          qrms(i)%m = amass(i)
         enddo

        return
        end subroutine fillAmbData
        
        subroutine fillCHARMMData ()
c    *******************************************************************
c    ** reads CHARMM style PSF file data to the parameter datatype    **
c    **                                                               **
c    *******************************************************************               
        
	  implicit none


          integer :: ii, jj
	     		
          integer iFirstSod, iNumProt
          integer protStart,protEnd, pTot
          integer lipidStart,lipidEnd, lTot
          integer ionStart,ionEnd, iTot
          integer wTot
          integer oTot
          integer testSelect

          !integer, parameter :: iPSFCONF = 5461
          real, pointer :: myMasses(:)
          real, pointer :: myCharges(:)
          real, pointer :: myRadii(:)
          integer, pointer :: mySelect(:)
          character(len=20) :: psfConfName = "psf.conf"
c    *******************************************************************

        call read_psf_conf(maxSect, 
     x       nWatStEnd, arrWatStart, arrWatEnd, 
     x       nProtStEnd, arrProtStart, arrProtEnd, 
     x       nIonStEnd, arrIonStart, arrIonEnd, 
     x       nLipStEnd, arrLipStart, arrLipEnd, 
     x       nOthStEnd, arrOthStart, arrOthEnd, 
     x       natom, 
     x       trim(psfConfName)  // CHAR(0))


         pTot = 0
         do ii = 1, nProtStEnd
          pTot = pTot + 
     x    (arrProtEnd(ii) - arrProtStart(ii)) + 1
c          print *, "p ", ii, arrProtStart(ii),
c     x             arrProtEnd(ii)
         enddo
c         print *, "pTot = ", pTot
        
         iTot = 0
         do ii = 1, nIonStEnd
          iTot = iTot + 
     x    (arrIonEnd(ii) - arrIonStart(ii)) + 1
c          print *, "i ", ii, arrIonStart(ii),
c     x             arrIonEnd(ii)
         enddo
c         print *, "iTot = ", iTot

         lTot = 0
         do ii = 1, nLipStEnd
          lTot = lTot + 
     x    (arrLipEnd(ii) - arrLipStart(ii)) + 1
c          print *, "l ", ii, arrLipStart(ii),
c     x             arrLipEnd(ii)
         enddo
c         print *, "lTot = ", lTot

         wTot = 0
         do ii = 1, nWatStEnd
          wTot = wTot + 
     x    (arrWatEnd(ii) - arrWatStart(ii)) + 1
c          print *, "w ", ii, arrWatStart(ii),
c     x             arrWatEnd(ii)
         enddo
c         print *, "wTot = ", wTot

         oTot = 0
         do ii = 1, nOthStEnd
          oTot = oTot + 
     x    (arrOthEnd(ii) - arrOthStart(ii)) + 1
c          print *, "o ", ii, arrOthStart(ii),
c     x             arrOthEnd(ii)
         enddo
c         print *, "oTot = ", oTot
c         print *, "atoms = ", natom

       ! Set these variables are set for backwards compatibility
       ! to all analysis routines. In the future, will
       ! need to make sure all routines use the loop over
       ! starts/ends
         iProtStart  = arrProtStart(1) 
         iProtEnd    = arrProtEnd(1)
 
         iIonStart   = arrIonStart(1)
         iIonEnd     = arrIonEnd(1)

         iLipStart   = arrLipStart(1)
         iLipEnd     = arrLipEnd(1)

         iWatStart   = arrWatStart(1)
         iWatEnd     = arrWatEnd(1)

         iOthStart   = arrOthStart(1)
         iOthEnd     = arrOthEnd(1)

ccc     ccccccccccccccccccccccccccccccccccccc
ccc     !!!! FOR DEBUG ONLY !!!!!
ccc     ccccccccccccccccccccccccccccccccccccc
!         allocate(myMasses(natom))
!         allocate(myCharges(natom))
!         allocate(myRadii(natom))
!         allocate(mySelect(natom))
!        if ( 0 .eq. 1) then       
ccc     ccccccccccccccccccccccccccccccccccccc
ccc     !!!! END FOR DEBUG ONLY !!!!!
ccc     ccccccccccccccccccccccccccccccccccccc

           write(*,*) ""		  
           !tmpStr = all("=",50)
           tmpStr = new("=========================")
           tmpStr = tmpStr + new("=========================")
           write(*,*) tmpStr%data(1:tmpStr%length)	             
           if(pTot > 1) then
		   write(*,*) "     Total number of protein atoms: ", pTot
           endif
           if(iTot > 1) then
		   write(*,*) "      Total number of counter ions: ",
     x                iTot
           endif
           if(wTot > 1) then
		   write(*,*) "            Total number of waters: ",
     x                params%numWaters     
		   write(*,*) "       Total number of water atoms: ",
     x                params%numWaters*3
           endif
           if(lTot > 1) then
		   write(*,*) "       Total number of lipid atoms: ",
     x                lTot
           endif
           if(oTot > 1) then
		   write(*,*) "       Total number of other atoms: ",
     x                oTot
           endif
		   write(*,*) "   Total number of atoms in system: ",
     x                natom
           write(*,*) tmpStr%data(1:tmpStr%length)	                      
		   write(*,*) ""
		   
        nres = params%numRes

         allocate(myMasses(natom))
         allocate(myCharges(natom))
         allocate(myRadii(natom))
         allocate(mySelect(natom))

         call open_psf(trim(topFile) // CHAR(0)) ! Append a char(0) so c knows the end of the string  

           call read_psf_file_head(myMasses, 
     x                            myCharges, 
     x                              myRadii,
     x                             doSelect,
     x                            doCluster,
     x                              resList,
     x                          fullResList,
     x                            MAX_ATOMS,
     x                          nClusterRes,
     x                             mySelect,
     x                              nSelect,
     x           trim(strSelect) // CHAR(0),
     x           trim(strResSel) // CHAR(0))


         if(doSelect) then
           allocate(selectIndices(nSelect)) ! allocate on master, must allocate on nodes later
           do i = 1, nSelect
            selectIndices(i) = mySelect(i)
           enddo
         endif

         call close_psf()     ! Close the psf file

ccc     ccccccccccccccccccccccccccccccccccccc
ccc     !!!! END FOR DEBUG ONLY !!!!!
ccc     ccccccccccccccccccccccccccccccccccccc
!        endif  ! TO BE REMOVED
ccc     ccccccccccccccccccccccccccccccccccccc
ccc     !!!! END FOR DEBUG ONLY !!!!!
ccc     ccccccccccccccccccccccccccccccccccccc

         allocate(qrms(natom)) ! allocate qrms on master, must allocate on nodes later
         

         testSelect = 0
         do i = 1, natom        
          qrms(i)%q = myCharges(i)
          qrms(i)%r = myRadii(i)
          qrms(i)%m = myMasses(i)

         enddo

         deallocate(myMasses)
         deallocate(myCharges)
         deallocate(myRadii)
         deallocate(mySelect)

        return
        end subroutine fillCHARMMData


        subroutine fillCHARMMDataOLD ()
c    *******************************************************************
c    ** reads CHARMM style PSF file data to the parameter datatype    **
c    **                                                               **
c    **                THIS VERSION IS DEPRECATED.                    **
c    **                                                               **
c    ** Use the fillCHARMMData routine instead of this version, since **
c    ** this one is only here for historical purposes.                **
c    **                                                               **
c    **                THIS VERSION IS DEPRECATED.                    **
c    **                                                               **
c    *******************************************************************               
        
	      implicit none


              integer :: ii, jj
	     		
		  integer iFirstSod, iNumProt
		  integer protStart,protEnd, pTot
		  integer lipidStart,lipidEnd, lTot
		  integer ionStart,ionEnd, iTot
		  integer wTot
		  integer oTot
                  integer testSelect

                  integer, parameter :: iPSFCONF = 5461
                  real, pointer :: myMasses(:)
                  real, pointer :: myCharges(:)
                  real, pointer :: myRadii(:)
                  integer, pointer :: mySelect(:)
c    *******************************************************************

        open(unit=iPSFCONF,file="psf.conf", status='old')
        
         read(iPSFCONF,*) iProtStart, iProtEnd
         read(iPSFCONF,*) iIonStart, iIonEnd
         read(iPSFCONF,*) lipidStart, lipidEnd
         read(iPSFCONF,*) iWatStart, iWatEnd
         read(iPSFCONF,*) iOthStart, iOthEnd
         read(iPSFCONF,*) natom
        
        close(iPSFCONF)

          pTot = (iProtEnd - iProtStart) + 1
          iTot = (iIonEnd - iIonStart) + 1
          lTot = (lipidEnd - lipidStart) + 1
          wTot = (iWatEnd - iWatStart) + 1
          oTot = (iOthEnd - iOthStart) + 1
           
           write(*,*) ""		  
           !tmpStr = all("=",50)
           tmpStr = new("=========================")
           tmpStr = tmpStr + new("=========================")
           write(*,*) tmpStr%data(1:tmpStr%length)	             
           if(pTot > 1) then
		   write(*,*) "     Total number of protein atoms: ", pTot
           endif
           if(iTot > 1) then
		   write(*,*) "      Total number of counter ions: ",
     x                iTot
           endif
           if(wTot > 1) then
		   write(*,*) "            Total number of waters: ",
     x                params%numWaters     
		   write(*,*) "       Total number of water atoms: ",
     x                params%numWaters*3
           endif
           if(lTot > 1) then
		   write(*,*) "       Total number of lipid atoms: ",
     x                lTot
           endif
           if(oTot > 1) then
		   write(*,*) "       Total number of other atoms: ",
     x                oTot
           endif
		   write(*,*) "   Total number of atoms in system: ",
     x                natom
           write(*,*) tmpStr%data(1:tmpStr%length)	                      
		   write(*,*) ""
		   
        nres = params%numRes

         allocate(myMasses(natom))
         allocate(myCharges(natom))
         allocate(myRadii(natom))
         allocate(mySelect(natom))

         call open_psf(trim(topFile) // CHAR(0)) ! Append a char(0) so c knows the end of the string  

c           call read_psf_file_head(myMasses, 
c     x                            myCharges, 
c     x                              myRadii,
c     x                             doSelect,
c     x                            doCluster,
c     x                              resList,
c     x                            MAX_ATOMS,
c     x                          nClusterRes,
c     x                             mySelect,
c     x                              nSelect,
c     x           trim(strSelect) // CHAR(0),
c     x           trim(strResSel) // CHAR(0))

           call read_psf_file_head(myMasses, 
     x                            myCharges, 
     x                              myRadii,
     x                             doSelect,
     x                            doCluster,
     x                              resList,
     x                          fullResList,
     x                            MAX_ATOMS,
     x                          nClusterRes,
     x                             mySelect,
     x                              nSelect,
     x           trim(strSelect) // CHAR(0),
     x           trim(strResSel) // CHAR(0))



cc         print *, "****** ON MASTER ******** "
cc         do ii = 1, 3
cc          do jj = 1, fullResList(ii)%numAtomsPerRes
cc           print *, " fullResList(",ii,")%atomid(",jj,") = ",
cc     x          fullResList(ii)%atomid(jj)
cc          enddo
cc         enddo
cc         print *, "****** ON MASTER ******** "
cc         print *, ""

         if(doSelect) then
           allocate(selectIndices(nSelect)) ! allocate on master, must allocate on nodes later
           do i = 1, nSelect
            selectIndices(i) = mySelect(i)
           enddo
         endif

         call close_psf()     ! Close the psf file

         allocate(qrms(natom)) ! allocate qrms on master, must allocate on nodes later
         

         testSelect = 0
         do i = 1, natom         
          qrms(i)%q = myCharges(i)
          qrms(i)%r = myRadii(i)
          qrms(i)%m = myMasses(i)

         enddo

         deallocate(myMasses)
         deallocate(myCharges)
         deallocate(myRadii)
         deallocate(mySelect)

        return
        end subroutine fillCHARMMDataOLD

       subroutine createBoxL()
c    *******************************************************************
c    ** creates a boxL MPI-datatype                                   **
c    **                                                               **
c    *******************************************************************       
       implicit none
       
       integer :: nTypes
       integer, pointer, dimension(:) :: offsets,types, counts

c    *******************************************************************
       nTypes = 1                   ! # of different variable types 
       allocate(offsets(0:nTypes-1))
       allocate(  types(0:nTypes-1))
       allocate( counts(0:nTypes-1))
       
       ! Give space for 1 x (x,y,z) = 3 reals
       offsets(0) = 0
       types(0)   = MPI_REAL
       counts(0)  = 3

       call MPI_TYPE_STRUCT(nTypes,counts,
     x           offsets,types,BOXL_TYPE,ierr)
       call MPI_TYPE_COMMIT(BOXL_TYPE,ierr)

       deallocate(offsets, types, counts)       
            
          
        return        
       end subroutine createBoxL
       
       subroutine createAtom()
c    *******************************************************************
c    ** creates a basic atom position MPI-datatype                    **
c    **                                                               **
c    *******************************************************************       
       implicit none
       
       integer :: nTypes
       integer, pointer, dimension(:) :: offsets,types, counts

c    *******************************************************************
       nTypes = 1                   ! # of different variable types 
       allocate(offsets(0:nTypes-1))
       allocate(  types(0:nTypes-1))
       allocate( counts(0:nTypes-1))
       
       ! Give space for 1 x (x,y,z) = 3 reals
       offsets(0) = 0
       types(0)   = MPI_REAL
       counts(0)  = 3

       call MPI_TYPE_STRUCT(nTypes,counts,
     x           offsets,types,ATOM_TYPE,ierr)
       call MPI_TYPE_COMMIT(ATOM_TYPE,ierr)

       deallocate(offsets, types, counts)       
            
          
        return        
       end subroutine createAtom
       
       subroutine createConf()
c    *******************************************************************
c    ** creates a MPI-datatype for entire configuration data          **
c    **                                                               **
c    **      NOTE: Creates internal derived type using _createAtom_   **
c    **                                                               **
c    *******************************************************************       
       implicit none
       
       integer :: nTypes
       integer, pointer, dimension(:) :: offsets,types, counts

c    *******************************************************************
       
       !call createAtom() ! Create the atom type first
       
       nTypes = 1                   ! # of different variable types 
       allocate(offsets(0:nTypes-1))
       allocate(  types(0:nTypes-1))
       allocate( counts(0:nTypes-1))
       
       ! Give space for natom x ATOM_TYPE
       offsets(0) = 0
       types(0)   = MPI_REAL !ATOM_TYPE
       counts(0)  = natom*3

       call MPI_TYPE_STRUCT(nTypes,counts,
     x           offsets,types,CONF_TYPE,ierr)
       !call MPI_TYPE_CONTIGUOUS(natom*3,
!     x             MPI_REAL,CONF_TYPE,ierr)
       call MPI_TYPE_COMMIT(CONF_TYPE,ierr)

       deallocate(offsets, types, counts)       
        
          print *, ""
          print *, rank, " Node finished createConf"            
          print *, ""  

        return        
       end subroutine createConf

       subroutine setThreads()
c    *******************************************************************
c    ** sets and b-casts the number of threads to the nodes           **
c    **                                                               **
c    *******************************************************************                      
         implicit none
         
       ! OMP section   
       integer :: iii, nTh
       integer :: OMP_GET_THREAD_NUM, OMP_GET_NUM_PROCS                     
       integer, pointer, dimension(:) :: allThreads
       logical :: useThreads
c    *******************************************************************       
        
c        nThreads = OMP_GET_NUM_PROCS()        
c        call OMP_SET_NUM_THREADS(nThreads)
!        print *, rank, " nThreads ", nThreads
        useThreads = .false.
        call MPI_BARRIER(world,ierr)
        if (rank .eq. master) then
         nTh = nThreads
         allocate(allThreads(nProc1))
         do iii = 1, nProc1
          call MPI_RECV(allThreads(iii), 1, MPI_INTEGER,
     x                  iii,200,world,status,ierr)
          nTh = nTh + allThreads(iii)
          if (allThreads(iii) .gt. 1) useThreads = .true.
         enddo

         write(*,*)""
         write(*,*) "================================================"

         if (useThreads) then
           tmpStr = new("   ") + new("USING") + new(" ")
     x              + intToStr(nTH) +  new(" THREADS ACROSS") 
     x       + new(" ") + intToStr(nProc) + new(" PROCESSORS")
            write(*,*) tmpStr%data(1:tmpStr%length)		              
           write(*,*) ""	           
           write(*,*) "          Node      Thread Count "
           write(*,*) "          ----      ------------ "           
           write(*,2001) rank, nThreads                      
           do iii = 1, nProc1
             write(*,2001) iii, allThreads(iii)        
           enddo
         else ! Not using threads
          write(*,*)  " Multi-threading not available on this system"
         endif ! End test for threads           
         
         write(*,*) "================================================"
         write(*,*)""        
         deallocate(allThreads)
         
        else ! slaves
        
        call MPI_SEND(nThreads, 1, MPI_INTEGER,
     x               master,200,world,ierr)         
        endif
        
c!$OMP PARALLEL PRIVATE(threadID)
c        threadID = OMP_GET_THREAD_NUM()
!        print *,  " Hi From node ", rank, 
!     x   " and thread ", threadID

c!$OMP END PARALLEL        
        

         
2001             format (10x,i3,12x,i3)                    
         return
       end subroutine


        subroutine readExclrList()
c    *******************************************************************
c    ** reads list of per-atom distance exclusions from a file        **
c    ** on master and
c    **                                                               **
c    *******************************************************************
         implicit none

         ! Locals
         integer :: iExclFile
         integer :: ii
         type(dString) :: printStr
c    *******************************************************************

        if (rank .eq. master) then

          iExclFile = getUnit()
         ! open file and read the the exclusion distances
         open(unit=iExclFile, file=trim(exclrFile), status='old')
         do ii = 1, natom
            read(iExclFile,*) exclrDist(ii)
         enddo
         close(iExclFile)

        endif ! End check for master
        
        call MPI_BCAST(exclrDist,natom
     x               ,MPI_REAL,master, world, ierr)        ! BCAST the exclusion distances

        return
        end subroutine



       subroutine readCellData()
c    *******************************************************************
c    ** reads cell data into memory from a file                       **
c    **                                                               **
c    *******************************************************************                      
         implicit none

        ! Locals
        integer :: iUnitC
c    *******************************************************************

        iUnitC = getUnit()
        open(unit=iUnitC,file=trim(cellFile), status='old')
        
        read(iUnitC,*) alpha
        read(iUnitC,*) beta
        read(iUnitC,*) gamma
        read(iUnitC,*) bX
        read(iUnitC,*) bY
        read(iUnitC,*) bZ
        
        close(iUnitC)
        
        return
       end subroutine readCellData

       subroutine setAvgFile()
c    *******************************************************************
c    ** reads avg. positions for all proteins then sends to everyone  **
c    **                                                               **
c    *******************************************************************                      
         implicit none
         
        ! Locals
        integer  :: ii, ij
c    *******************************************************************
            
        ! Allocate the necessities            
        allocate(xAvg(iProtEnd))
        allocate(yAvg(iProtEnd))
        allocate(zAvg(iProtEnd))
        allocate(drxAvg(iProtEnd))
        allocate(dryAvg(iProtEnd))
        allocate(drzAvg(iProtEnd))
        allocate(qdrsq(iProtEnd))
            
        if (rank .eq. master) then        
         ! open file and read the atomic average positions
         open(unit=53, file='avgProt.dat', status='old')
         do ii = 1, iProtEnd
           
            read(53,*) ij, 
     x        xAvg(ii), yAvg(ii), zAvg(ii)
         enddo            
         close(53)        
        endif ! End check for master
        
        call MPI_BCAST(xAvg,iProtEnd 
     x               ,MPI_REAL,master, world, ierr)        ! BCAST all averages  
        call MPI_BCAST(yAvg,iProtEnd 
     x               ,MPI_REAL,master, world, ierr)        ! BCAST all averages
        call MPI_BCAST(zAvg,iProtEnd 
     x               ,MPI_REAL,master, world, ierr)        ! BCAST all averages
     
        return
       end subroutine setAvgFile

       subroutine setSelection()
c    *******************************************************************
c    ** allocates and sets the selection indices on all nodes         **
c    **                                                               **
c    *******************************************************************                      
         implicit none
         
        ! Locals
        integer :: i
c    *******************************************************************
            
        if (rank .ne. master) then        
         ! Allocate the selection on the nodes,
         ! the master already has this allocated
         allocate(selectIndices(nSelect))
        endif ! End check for master

        call MPI_BARRIER(world,ierr)
        
        call MPI_BCAST(selectIndices,nSelect 
     x               ,MPI_INTEGER,master, world, ierr)   ! BCAST selection indices

        return
       end subroutine setSelection
       
       subroutine setMDParm()
c    *******************************************************************
c    ** sets and b-casts all important MD parameters to all nodes     **
c    **                                                               **
c    *******************************************************************                      
         implicit none
         
        ! Locals
        integer,parameter :: nIndices=18
        integer,parameter :: nLogicals=46
        integer,parameter :: nReals=9        
        integer, dimension(nIndices) :: iIndices        
        logical, dimension(nLogicals) :: lDecisions
        real, dimension(nReals) :: rReals
c    *******************************************************************
   
   
       if (rank .eq. master) then ! pack data on master
        iIndices(1) = natom
        iIndices(2) = iWatStart
        iIndices(3) = iProtEnd
        iIndices(4) = iProbeAtom   
        iIndices(5) = nConf
        iIndices(6) = iNumRdaD
        iIndices(7) = iNumRdaA
        iIndices(8) = iNumPolDB
        iIndices(9) = GPUsPerNode
        iIndices(10) = iProtStart
        iIndices(11) = iIonStart
        iIndices(12) = iIonEnd
        iIndices(13) = nSelect
        iIndices(14) = iWatEnd
        iIndices(15) = iOthStart
        iIndices(16) = iOthEnd
        iIndices(17) = iLipStart
        iIndices(18) = iLipEnd
        
        lDecisions(1) = isBinary
        lDecisions(2) = doMPIIO
        lDecisions(3) = doProt
        lDecisions(4) = doSolv
        lDecisions(5) = doCoul
        lDecisions(6) = doInduced
        lDecisions(7) = doIndCrg
        lDecisions(8) = doGrCalc
        lDecisions(9) = doAllEner
        lDecisions(10) = isTruncated
        lDecisions(11) = doRdaCalc        
        lDecisions(12) = doMdaCalc
        lDecisions(13) = doIons 
        lDecisions(14) = doDIons
        lDecisions(15) = doSolut
        lDecisions(16) = doShell
        lDecisions(17) = doAvgS
        lDecisions(18) = doQdr
        lDecisions(19) = doCBind
        lDecisions(20) = doGZRead
        lDecisions(21) = doPFile
        lDecisions(22) = doShellM
        lDecisions(23) = doShellP2
        lDecisions(24) = doPOTL
        lDecisions(25) = doProtM
        lDecisions(26) = useCUDA
        lDecisions(27) = doShellCalc
        lDecisions(28) = doDQEfield
        lDecisions(29) = doVPolCorr              
        lDecisions(30) = doQPolCorr
        lDecisions(31) = doDCDRead
        lDecisions(32) = doPSFRead
        lDecisions(33) = doCovar
        lDecisions(34) = doEFS
        lDecisions(35) = doProtV
        lDecisions(36) = doSelect
        lDecisions(37) = doMultiTraj
        lDecisions(38) = doCluster
        lDecisions(39) = doEnerGap
        lDecisions(40) = doDebug
        lDecisions(41) = doGyr
        lDecisions(42) = doInert
        lDecisions(43) = doOthers
        lDecisions(44) = doHoleV
        lDecisions(45) = doLipid
        lDecisions(46) = doExclr

        rReals(1) = real(alpha,kind=4)
        rReals(2) = real(beta,kind=4)
        rReals(3) = real(gamma,kind=4)        
        rReals(4) = real(bX,kind=4)
        rReals(5) = real(bY,kind=4)
        rReals(6) = real(bZ,kind=4)
        rReals(7) = rBindCut
        rReals(8) = rClusterCut
        rReals(9) = holeRadius
       endif
       
       call MPI_BCAST(iIndices,nIndices
     x               ,MPI_INTEGER,master,world,ierr)       ! BCAST all indices       
       call MPI_BCAST(lDecisions,nLogicals 
     x               ,MPI_LOGICAL,master, world, ierr)     ! BCAST all logical variables     
       call MPI_BCAST(rReals,nReals 
     x               ,MPI_REAL,master, world, ierr)        ! BCAST all real variables     
     
       if (rank .ne. master) then ! re-sort on the nodes
        natom       = iIndices(1)
        iWatStart   = iIndices(2) 
        iProtEnd    = iIndices(3) 
        iProbeAtom  = iIndices(4)
        nConf       = iIndices(5)
        iNumRdaD    = iIndices(6)
        iNumRdaA    = iIndices(7)
        iNumPolDB   = iIndices(8)
        GPUsPerNode = iIndices(9)
        iProtStart  = iIndices(10)
        iIonStart   = iIndices(11)
        iIonEnd     = iIndices(12)
        nSelect     = iIndices(13)
        iWatEnd     = iIndices(14)
        iOthStart = iIndices(15)
        iOthEnd   = iIndices(16)
        iLipStart   = iIndices(17)
        iLipEnd     = iIndices(18)

        
        isBinary    = lDecisions(1)
        doMPIIO     = lDecisions(2)
        doProt      = lDecisions(3)
        doSolv      = lDecisions(4)
        doCoul      = lDecisions(5)
        doInduced   = lDecisions(6)
        doIndCrg    = lDecisions(7)        
        doGrCalc    = lDecisions(8)
        doAllEner   = lDecisions(9)
        isTruncated = lDecisions(10)
        doRdaCalc   = lDecisions(11)        
        doMdaCalc   = lDecisions(12)        
        doIons      = lDecisions(13)
        doDIons     = lDecisions(14)
        doSolut     = lDecisions(15) 
        doShell     = lDecisions(16)
        doAvgS      = lDecisions(17)
        doQdr       = lDecisions(18)
        doCBind     = lDecisions(19)
        doGZRead    = lDecisions(20)
        doPFile     = lDecisions(21)
        doShellM    = lDecisions(22)        
        doShellP2   = lDecisions(23)
        doPOTL      = lDecisions(24)
        doProtM     = lDecisions(25)
        useCUDA     = lDecisions(26)
        doShellCalc = lDecisions(27)
        doDQEfield  = lDecisions(28)
        doVPolCorr  = lDecisions(29)
        doQPolCorr  = lDecisions(30)
        doDCDRead   = lDecisions(31)
        doPSFRead   = lDecisions(32)
        doCovar     = lDecisions(33)
        doEFS       = lDecisions(34)
        doProtV     = lDecisions(35)
        doSelect    = lDecisions(36)
        doMultiTraj = lDecisions(37)
        doCluster   = lDecisions(38)
        doEnerGap   = lDecisions(39)
        doDebug     = lDecisions(40)
        doGyr       = lDecisions(41)
        doInert     = lDecisions(42)
        doOthers    = lDecisions(43)
        doHoleV     = lDecisions(44)
        doLipid     = lDecisions(45)
        doExclr     = lDecisions(46)

        alpha = real(rReals(1),kind=8)
        beta  = real(rReals(2),kind=8)        
        gamma = real(rReals(3),kind=8)        
        bX    = real(rReals(4),kind=8)        
        bY    = real(rReals(5),kind=8)        
        bZ    = real(rReals(6),kind=8)
        rBindCut = rReals(7)
        rClusterCut = rReals(8)
        holeRadius = rReals(9)
       endif
       
       box(1) = bX
       box(2) = bY
       box(3) = bZ
             
       sendBoxL = .true. ! We all need box lengths
c       print *, rank," iWatStart = ", iWatStart
c       print *, rank," iProtEnd = ", iProtEnd
c              doSolv = .true.
c              doProt = .true.
c              doIons = .false.
c              doCoul = .true.
              
         return
       end subroutine setMDParm

       subroutine setStartEnd()
c    *******************************************************************
c    ** sets up the starting and ending indices read from the         **
c    ** psf.conf file                                                 **
c    **                                                               **
c    *******************************************************************               
       implicit none

       integer, dimension(5) :: allNums
c    *******************************************************************

            if (rank .eq. master) then
              allNums(1) = nWatStEnd
              allNums(2) = nProtStEnd
              allNums(3) = nIonStEnd
              allNums(4) = nLipStEnd
              allNums(5) = nOthStEnd
            endif

             call MPI_BCAST(allNums,5
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST all number of all starts/ends

             call MPI_BCAST(arrWatStart,maxSect
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST all number of starting waters
             call MPI_BCAST(arrWatEnd,maxSect
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST all number of ending waters

             call MPI_BCAST(arrProtStart,maxSect
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST all number of starting proteins
             call MPI_BCAST(arrProtEnd,maxSect
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST all number of ending proteins

             call MPI_BCAST(arrIonStart,maxSect
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST all number of starting ions
             call MPI_BCAST(arrIonEnd,maxSect
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST all number of ending ions

             call MPI_BCAST(arrLipStart,maxSect
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST all number of starting lipids
             call MPI_BCAST(arrLipEnd,maxSect
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST all number of ending lipids

             call MPI_BCAST(arrOthStart,maxSect
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST all number of starting other atoms
             call MPI_BCAST(arrOthEnd,maxSect
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST all number of ending other atoms

              ! set all the number of all index types
              nWatStEnd  = allNums(1)
              nProtStEnd = allNums(2)
              nIonStEnd  = allNums(3)
              nLipStEnd  = allNums(4)
              nOthStEnd  = allNums(5)

         return
       end subroutine setStartEnd
              
       subroutine setMulti()
c    *******************************************************************
c    ** sets up the multiple trajectory parameters for all files      **
c    **                                                               **
c    *******************************************************************               
       implicit none

       character(len=100) :: tmpLine
       integer :: ii, ios
       integer :: iArr
       type(dstring)      :: dArray(2)
c    *******************************************************************

       ! Use the nConf line as the number of lines in our file
       nMultiTraj = nConf

       if(rank .eq. master) then
         ! Loop over all lines in the database and count the number of entries
         open(unit=53, file=trim(datFile), status='old')
         nMultiTraj = 0
         do
           read ( 53, '(a)', iostat = ios ) tmpLine
           if ( ios /= 0 ) then
             exit
           else
             nMultiTraj = nMultiTraj + 1
           end if
         enddo
         rewind(53)

         endif

           call MPI_BCAST(nMultiTraj,1
     x      ,MPI_INTEGER,master,world,ierr) ! BCAST nMultiTraj       

           do ii = 1, nMultiTraj           

           if(rank .eq. master) then
            read(53,'(a)') tmpLine
           endif

           call MPI_BCAST(tmpLine,100
     x      ,MPI_CHARACTER,master,world,ierr) ! BCAST the filename       

            call createStrArr(tmpLine,dArray,2,iArr)
            mTrajNames(ii)   = darray(1)
            mTrajLengths(ii) = dStrToInt(darray(2))
 
          enddo           
 
        if(rank .eq. master) close(53)        

        nConf = 0
        do ii = 1, nMultiTraj
         nConf = nConf + mTrajLengths(ii)
        enddo

        datFile = ''
        datFile = mTrajNames(1)%data(1:mTrajNames(1)%length)
        idxTraj = 1
        nFrameTraj = 0
        lastTraj = 0

        ! Open files on all nodes or master
        if(doMPIIO) then
           call readHeadMPIIO()
        else
           if(rank .eq. master) then
             call readHeader()
           endif
        endif

        return
       end subroutine setMulti

       subroutine multiTraj()
c    *******************************************************************
c    ** checks if the next traj file is required, and if so, the      **
c    ** previous traj is closed and the next is opened.               **
c    **                                                               **
c    *******************************************************************               
       implicit none

       integer :: dTraj
       integer :: nowTraj
c    *******************************************************************

          nowTraj = iConf+nProc - 1
          dTraj = nowTraj - lastTraj
          nFrameTraj = nFrameTraj + dTraj

          if( nFrameTraj .eq. mTrajLengths(idxTraj) ) then
            idxTraj = idxTraj + 1
            nFrameTraj = 0

           ! the last file is closed after the main loop
           if(idxTraj .le. nMultiTraj) then

             ! change datFile to the new file name
             datFile = mTrajNames(idxTraj)%data(
     x             1:mTrajNames(idxTraj)%length )

            ! close then open the files on master or the nodes
            if(doMPIIO) then
              call MPI_FILE_CLOSE(idat, ierr) ! mpi-io traj. (all nodes) 
              call readHeadMPIIO()
            else
              if(rank .eq. master) then
               close(idat)          ! unix-io traj. (master only)
               call readHeader()
             endif
            endif

            ! lets make sure to wait for everyone (mainly master)
            !call MPI_BARRIER(world, ierr)

           endif ! check for idxTraj < nMT


          endif ! check if we need to change files
         
         lastTraj = nowTraj
        return
       end subroutine multiTraj
              
       subroutine setDatParm()
c    *******************************************************************
c    ** creates a parameter MPI-datatype and initializes on all nodes **
c    **                                                               **
c    *******************************************************************               
       implicit none

       !Args
!       type(dat_parm) params
              
       ! Locals
       integer :: nTypes
       integer, pointer, dimension(:) :: offsets,types, counts

c    *******************************************************************

       nTypes = 1                   ! # of different variable types 
       allocate(offsets(0:nTypes-1))
       allocate(  types(0:nTypes-1))
       allocate( counts(0:nTypes-1))
       
       ! Give space for (4 + PDB_MAXRES) integers
       offsets(0) = 0
       types(0)   = MPI_INTEGER
       counts(0)  = 4 + PDB_MAXRES

       call MPI_TYPE_STRUCT(nTypes,counts,
     x          offsets,types,PARM_TYPE,ierr)
       call MPI_TYPE_COMMIT(PARM_TYPE,ierr)
       deallocate(offsets, types, counts)

       
       call MPI_BCAST(params,1,PARM_TYPE,master, world, ierr) 

       nWaters = params%numWaters
       
        return        
       end subroutine setDatParm

       subroutine setResList()
c    *******************************************************************
c    ** creates a resList MPI-datatype and initializes on all nodes   **
c    **                                                               **
c    *******************************************************************               
       implicit none

       !Args
              
       ! Locals
       integer :: nTypes
       integer, pointer, dimension(:) :: offsets,types, counts

c    *******************************************************************

       nTypes = 1                   ! # of different variable types 
       allocate(offsets(0:nTypes-1))
       allocate(  types(0:nTypes-1))
       allocate( counts(0:nTypes-1))
       
       ! Give space for ((2 + MAX_ATOMS)*MAX_RESIDUES) integers
       offsets(0) = 0
       types(0)   = MPI_INTEGER
       counts(0)  = (2 + MAX_ATOMS)*MAX_RESIDUES

       ! Create MPI structure for the resList, and bcast
       call MPI_TYPE_STRUCT(nTypes,counts,
     x          offsets,types,RESL_TYPE,ierr)
       call MPI_TYPE_COMMIT(RESL_TYPE,ierr)
       
       call MPI_BCAST(resList,1,RESL_TYPE,
     x               master, world, ierr) 
       
       call MPI_BCAST(fullResList,1,RESL_TYPE,
     x                    master, world, ierr) 


       call MPI_BCAST(nClusterRes,1,MPI_INTEGER,
     x               master, world, ierr) 


       ! Create MPI structure for the aggList
       ! first, Give space for ((1 + MAX_PER_CLUSTER)*nClusterRes) integers
       offsets(0) = 0
       types(0)   = MPI_INTEGER
       counts(0)  = (1 + MAX_PER_CLUSTER)*MAX_AGG

       call MPI_TYPE_STRUCT(nTypes,counts,
     x          offsets,types,AGG_TYPE,ierr)
       call MPI_TYPE_COMMIT(AGG_TYPE,ierr)


       ! Create MPI structure for aggGyr and aggEccen
       ! first, Give space for ((MAX_PER_CLUSTER)*nClusterRes) reals
       offsets(0) = 0
       types(0)   = MPI_REAL
       counts(0)  = MAX_AGG

       call MPI_TYPE_STRUCT(nTypes,counts,
     x          offsets,types,AGG_PARM_TYPE,ierr)
       call MPI_TYPE_COMMIT(AGG_PARM_TYPE,ierr)


       deallocate(offsets, types, counts)

        return        
       end subroutine setResList

       subroutine setCBind()
c    *******************************************************************
c    ** sends out average total/delta binding energy to all nodes     **
c    **                                                               **
c    *******************************************************************               
       implicit none
       
       ! Locals

c    *******************************************************************
       
        ! send out the averages...
        call MPI_BCAST(epsilonAvg,1,MPI_REAL,master, world, ierr) 
        call MPI_BCAST(DepsilonAvg,1,MPI_REAL,master, world, ierr)        
        call MPI_BCAST(simTemp,1,MPI_REAL,master,world,ierr)
        
        ! ... and convert them to real8
        r8epsAvg  = real(epsilonAvg,kind=8)
        r8epsDAvg = real(DepsilonAvg,kind=8)
        betaT = 300.d0*beta300/real(simTemp,kind=8)
        
        return        
       end subroutine setCBind

       subroutine setQRM()
c    *******************************************************************
c    ** creates a qrms MPI-datatype and initializes on all nodes      **
c    **                                                               **
c    *******************************************************************               
       implicit none
       
       ! Locals
       integer :: nTypes
       integer, pointer, dimension(:) :: offsets,types, counts

c    *******************************************************************

       nTypes = 1                   ! # of different variable types 
       allocate(offsets(0:nTypes-1))
       allocate(  types(0:nTypes-1))
       allocate( counts(0:nTypes-1))
       
       ! Give space for 3 reals
       offsets(0) = 0
       types(0)   = MPI_REAL
       counts(0)  = 3

       call MPI_TYPE_STRUCT(nTypes,counts,offsets,types,QRM_TYPE,ierr)
       call MPI_TYPE_COMMIT(QRM_TYPE,ierr)
       deallocate(offsets, types, counts)

       if (rank .ne. master) allocate(qrms(natom)) ! allocate qrms on nodes, master has already allocated
       
       call MPI_BCAST(qrms,natom,QRM_TYPE,master, world, ierr) 
                 
        return        
       end subroutine setQRM

       subroutine initDADDbase(dBaseFileD,dBaseFileA)
c    *******************************************************************
c    ** creates a database for atoms involved in the r(D-A) distance  **
c    ** calculations.  this subroutine assumes the following file     **
c    ** format:                                                       **
c    **                                                               **
c    ** +--------------------+                                        **
c    ** | iAtom_1   iRes_1   |                                        **
c    ** | iAtom_2   iRes_2   |                                        **
c    ** | ... and so on ...  |                                        **
c    ** +--------------------+                                        **
c    **                                                               **
c    ** principal variables:                                          **
c    ** character(len=*)   dBaseFile - the database of the atoms      **
c    **                    involved in the SS correlation function.   **
c    **                                                               **
c    *******************************************************************
       implicit none

       ! Args
       character(len=*) :: dBaseFileD,dBaseFileA

       ! Locals
       integer, parameter :: iFileNum=91
       integer            :: ios
       integer            :: iNumLinesD,iNumLinesA
       integer            :: iLine
       character(len=80)  :: string
c    *******************************************************************

       !!!!!!!!!!!!!!!!!!!!!
       ! DONOR !!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Loop over all lines in the database and count the number of entries
       open(unit=iFileNum, file=trim(dBaseFileD), status='old')
       iNumLinesD = 0
       do
         read ( iFileNum, '(a)', iostat = ios ) string
         if ( ios /= 0 ) then
             exit
         else
             iNumLinesD = iNumLinesD + 1
         end if
       enddo

       rewind(iFileNum)
       iNumRdaD = iNumLinesD        
       ! Allocate and store all donor database entries
       allocate(ariD(iNumRdaD))
       do iLine = 1, iNumRdaD
         read(iFileNum,*) ariD(iLine)%iA,ariD(iLine)%iR
       enddo
       close(iFileNum)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!!!!!!!!!!!!!!!!!
       ! ACCEPTOR !!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Loop over all lines in the database and count the number of entries
       open(unit=iFileNum, file=trim(dBaseFileA), status='old')
       iNumLinesA = 0
       do
         read ( iFileNum, '(a)', iostat = ios ) string
         if ( ios /= 0 ) then
             exit
         else
             iNumLinesA = iNumLinesA + 1
         end if
       enddo

       rewind(iFileNum)
       iNumRdaA = iNumLinesA        
       ! Allocate and store all donor database entries
       allocate(ariA(iNumRdaA))
       do iLine = 1, iNumRdaA
         read(iFileNum,*) ariA(iLine)%iA,ariA(iLine)%iR
       enddo
       close(iFileNum)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       return
       end subroutine initDADDbase

       subroutine initDVosDbase(dBaseFile)
c    *******************************************************************
c    ** creates a database for atoms involved in the Stokes Shift     **
c    ** correlation function.  this subroutine assumes the following  **
c    ** file format:                                                  **
c    **                                                               **
c    ** +-------------------------------+                             **
c    ** | iAtom_1   iRes_1    rCharge_1 |                             **
c    ** | iAtom_2   iRes_2    rCharge_2 |                             **
c    ** | ... and so on ...             |                             **
c    ** +-------------------------------+                             **
c    **                                                               **
c    ** principal variables:                                          **
c    ** character(len=*)   dBaseFile - the database of the atoms      **
c    **                    involved in the SS correlation function.   **
c    **                                                               **
c    *******************************************************************
       implicit none

       ! Args
       character(len=*) :: dBaseFile

       ! Locals
       integer, parameter :: iFileNum=89
       integer            :: ios
       integer            :: iNumLines
       integer            :: iLine
       character(len=80)  :: string
c    *******************************************************************

       ! Loop over all lines in the database and count the number of entries
       open(unit=iFileNum, file=trim(dBaseFile), status='old')
       iNumLines = 0
       do
         read ( iFileNum, '(a)', iostat = ios ) string
         if ( ios /= 0 ) then
             exit
         else
             iNumLines = iNumLines + 1
         end if
       enddo

       rewind(iFileNum)
       iNumDB = iNumLines
       ! Allocate and store all database entries
       allocate(dVosDbase(iNumDB))
       do iLine = 1, iNumLines
         read(iFileNum,*) dVosDbase(iLine)%atom,
     x       dVosDbase(iLine)%residue,dVosDbase(iLine)%q
       enddo
       close(iFileNum)

       return
       end subroutine initDVosDbase

       subroutine initDVosPolDbase(dBaseFile, fFile)
c    *******************************************************************
c    ** creates a database for atoms involved in polarizable          **
c    ** simulations.  this subroutine assumes the following           **
c    ** database file format:                                         **
c    **                                                               **
c    ** +-------------------------------+                             **
c    ** | iAtom_1   iRes_1    dCharge_1 |                             **
c    ** | iAtom_2   iRes_2    dCharge_2 |                             **
c    ** | ... and so on ...             |                             **
c    ** +-------------------------------+                             **
c    ** and f-parameter file format:                                  **
c    **                                                               **
c    ** +-------------------+                                         **
c    ** | fFac_1            |                                         **
c    ** | fFac_2            |                                         **
c    ** | ... and so on ... |                                         **
c    ** +-------------------+                                         **
c    **                                                               **
c    ** principal variables:                                          **
c    ** character(len=*)   dBaseFile - the database of the atoms      **
c    **                    involved as the polarizable center.        **
c    ** character(len=*)   fFile - the trajectory of tje f-parm.      **
c    **                                                               **
c    *******************************************************************
       implicit none

       ! Args
       character(len=*) :: dBaseFile, fFile

       ! Locals
       integer, parameter :: iFileNum=101
       integer            :: ios
       integer            :: iNumLines
       integer            :: iLine
       character(len=80)  :: string
c    *******************************************************************

       ! Loop over all lines in the database and count the number of entries
       open(unit=iFileNum, file=trim(dBaseFile), status='old')
       iNumLines = 0
       do
         read ( iFileNum, '(a)', iostat = ios ) string
         if ( ios /= 0 ) then
             exit
         else
             iNumLines = iNumLines + 1
         end if
       enddo

       rewind(iFileNum)
       iNumPolDB = iNumLines
       ! Allocate and store all database entries
       allocate(dVosPolDbase(iNumPolDB))
       do iLine = 1, iNumPolDB
         read(iFileNum,*) dVosPolDbase(iLine)%atom,
     x       dVosPolDbase(iLine)%residue,dVosPolDbase(iLine)%q
       enddo
       close(iFileNum)
      !----------------------------------------------!
       ! Loop over all lines in the database and count the number of entries
       open(unit=iFileNum, file=trim(fFile), status='old')
       iNumLines = 0
       do
         read ( iFileNum, '(a)', iostat = ios ) string
         if ( ios /= 0 ) then
             exit
         else
             iNumLines = iNumLines + 1
         end if
       enddo

       rewind(iFileNum)
       iNumFDB = iNumLines
       ! Allocate and store all f-parameter entries
       allocate(fFacTraj(iNumFDB))
       do iLine = 1, iNumFDB
         read(iFileNum,*) fFacTraj(iLine)
       enddo
       close(iFileNum)
      
       return
       end subroutine initDVosPolDbase
       
       subroutine initAlphaDbase(dBaseFile)
c    *******************************************************************
c    ** creates an atomic polarizability database for atoms in the    **
c    ** system.  this subroutine assumes the following file format:   **
c    **                                                               **
c    **   NOTE: All polarizabilities are assumed to be in units       **
c    **         of Ang^3.                                             **
c    **                                                               **
c    ** +---------------------+                                       **
c    ** | iAtom_1    rAlpha_1 |                                       **
c    ** | iAtom_2    rAlpha_2 |                                       **
c    ** | ... and so on ...   |                                       **
c    ** +---------------------+                                       **
c    **                                                               **
c    ** principal variables:                                          **
c    ** character(len=*)   dBaseFile - the database of the atoms      **
c    **                                                               **
c    *******************************************************************
       implicit none

       ! Args
       character(len=*) :: dBaseFile

       ! Locals
       integer, parameter :: iFileNum=90
       integer            :: iNumLines
       integer            :: ii, iLine
c    *******************************************************************

       ! Allocate and store all database entries
       allocate(alphaDbase(natom))       
       open(unit=iFileNum, file=trim(dBaseFile), status='old')

       ! Loop over all lines in the database and count the number of entries
       do ii = 1, natom
         read(iFileNum,*) iLine, alphaDbase(ii)
       enddo
       close(iFileNum)

       return
       end subroutine initAlphaDbase

       subroutine createRDAR()
c    *******************************************************************
c    ** creates a RDAR MPI-datatype                                   **
c    **                                                               **
c    *******************************************************************

       implicit none

       integer :: nTypes
       integer, pointer, dimension(:) :: offsets, types, counts
c    *******************************************************************

       nTypes = 1                  ! # of different variable types
       allocate(offsets(0:nTypes-1))
       allocate(  types(0:nTypes-1))
       allocate( counts(0:nTypes-1))

       ! Give space for atom, residue
       offsets(0) = 0
       types(0)   = MPI_REAL
       counts(0)  = 7

       call MPI_TYPE_STRUCT(nTypes,counts,offsets,
     x                 types,RDAR_TYPE,ierr)
       call MPI_TYPE_COMMIT(RDAR_TYPE,ierr)

       deallocate(offsets, types, counts)

        return
       end subroutine createRDAR

       subroutine createMDAR()
c    *******************************************************************
c    ** creates a MDAR MPI-datatype                                   **
c    **                                                               **
c    *******************************************************************

       implicit none

       integer :: nTypes
       integer, pointer, dimension(:) :: offsets, types, counts
c    *******************************************************************

       nTypes = 1                  ! # of different variable types
       allocate(offsets(0:nTypes-1))
       allocate(  types(0:nTypes-1))
       allocate( counts(0:nTypes-1))

       ! Give space for atom, residue
       offsets(0) = 0
       types(0)   = MPI_REAL
       counts(0)  = 4

       call MPI_TYPE_STRUCT(nTypes,counts,offsets,
     x                 types,MDAR_TYPE,ierr)
       call MPI_TYPE_COMMIT(MDAR_TYPE,ierr)

       deallocate(offsets, types, counts)

        return
       end subroutine createMDAR      
       
       subroutine createDADbase()
c    *******************************************************************
c    ** creates a DADbase MPI-datatype                                **
c    **                                                               **
c    *******************************************************************

       implicit none

       integer :: nTypes
       integer, pointer, dimension(:) :: offsets, types, counts
       integer :: mpi_intExt
c    *******************************************************************

       nTypes = 1                  ! # of different variable types
       allocate(offsets(0:nTypes-1))
       allocate(  types(0:nTypes-1))
       allocate( counts(0:nTypes-1))

       ! Give space for atom, residue
       offsets(0) = 0
       types(0)   = MPI_INTEGER
       counts(0)  = 2

       call MPI_TYPE_STRUCT(nTypes,counts,offsets,
     x                 types,DADB_TYPE,ierr)
       call MPI_TYPE_COMMIT(DADB_TYPE,ierr)

       deallocate(offsets, types, counts)

        return
       end subroutine createDADbase
       
       subroutine createDVosDbase()
c    *******************************************************************
c    ** creates a dVosDbase MPI-datatype                              **
c    **                                                               **
c    *******************************************************************

       implicit none

       integer :: nTypes
       integer, pointer, dimension(:) :: offsets, types, counts
       integer :: mpi_intExt
c    *******************************************************************

       nTypes = 2                   ! # of different variable types
       allocate(offsets(0:nTypes-1))
       allocate(  types(0:nTypes-1))
       allocate( counts(0:nTypes-1))

       ! Give space for atom, residue
       offsets(0) = 0
       types(0)   = MPI_INTEGER
       counts(0)  = 2

       ! Give space for q
       call MPI_TYPE_EXTENT(MPI_INTEGER,mpi_intExt,ierr)
       offsets(1) = counts(0)*mpi_intExt
       types(1)   = MPI_REAL
       counts(1)  = 1

       call MPI_TYPE_STRUCT(nTypes,counts,offsets,
     x                 types,DVOSDB_TYPE,ierr)
       call MPI_TYPE_COMMIT(DVOSDB_TYPE,ierr)

       deallocate(offsets, types, counts)

        return
       end subroutine createDVosDbase

       subroutine sendAAcidList()
c    *******************************************************************
c    ** create and send a copy of the lookup table for determining    **
c    ** the amino acid index based on the atomic (nodeConf) index     **
c    **                                                               **
c    *******************************************************************
       implicit none
       
       ! Locals
       integer :: ii, jj,ij, nRes
c    *******************************************************************

          if(rank .eq. master) nRes = params%numRes
          
          call MPI_BCAST(nRes,1,MPI_INTEGER,
     x                     master,world,ierr)  
     
          ! Make sure everybody allocates ...
          allocate(aminoAcidList(iProtEnd))
          allocate(nShellHighAA(nRes))
          allocate(nShellLowAA(nRes))          
          allocate(ebindTotAA(nRes))
          allocate(ebindDelAA(nRes))
          
          if (rank .eq. master) then ! Create on the master          
           ij = 0
           
           do ii = 1, params%numRes
            do jj = 1, params%numAtoms(ii)                     
             ij = ij + 1
             aminoAcidList(ij) = ii           
             enddo
            enddo   
            
          endif ! end check for master
          
          call MPI_BARRIER(world,ierr)
          
          ! ... then send the lookup table to everybody
          call MPI_BCAST(aminoAcidList,iProtEnd,MPI_INTEGER,
     x                     master,world,ierr)                   


        return
       end subroutine sendAAcidList
       
       subroutine sendAlphaDbase()
c    *******************************************************************
c    ** send a copy of the atomic polarizability database from        **
c    ** the master to all nodes                                       **
c    **                                                               **
c    *******************************************************************
       implicit none
       
c    *******************************************************************

          ! Make sure everybody allocates ...
          if(rank .gt. master) allocate(alphaDbase(natom))

          ! ... then send the database
          call MPI_BCAST(alphaDbase,natom,MPI_REAL,
     x                     master,world,ierr)                   


        return
       end subroutine sendAlphaDbase
       
       subroutine sendDVosDbase()
c    *******************************************************************
c    ** send a copy of the dVosDbase from the master to all nodes     **
c    **                                                               **
c    *******************************************************************
       implicit none
       
       ! Locals
        integer :: ii, jj, iCount, iii, iNum
        integer :: shift
        real :: tmpQ
c    *******************************************************************


          if(doPSFRead) then
            shift = arrProtStart(1) - 1
          else
            shift = 0
          endif

          ! Send the total number of elements ...
          call MPI_BCAST(iNumDB,1,MPI_INTEGER,
     x                     master,world,ierr)
          if(rank .gt. master) allocate(dVosDbase(iNumDB))

          ! ... then send the database
          call MPI_BCAST(dVosDbase,iNumDB,DVOSDB_TYPE,
     x                     master,world,ierr)

          
           do iii = 1, iNumDB
            iCount = 0
           
            do ii = 1, dVosDbase(iii)%residue
             do jj = 1, params%numAtoms(ii)

              iCount = iCount + 1             
              if ( jj .eq. dVosDbase(iii)%atom ) then              
               iNum = iCount
              endif

             enddo ! end atom loop
            enddo ! end res loop

            ! vvv WORKING vvv
            !dVosDbase(iii)%atom = iNum ! reset to atomic index in nodeconf space
            ! ^^^ WORKING ^^^       

            ! TODO: INCLUDE A SHIFT PARAMETER THAT PUSHES THE iNum variable to the
            !       PROPER LOCATION IN nodeconf space (e.g. 133657)    
            !dVosDbase(iii)%atom = iNum + 133657! reset to atomic index in nodeconf space
            dVosDbase(iii)%atom = iNum + shift! reset to atomic index in nodeconf space
           enddo ! End database loop


           if(rank .eq. master) then
           write(*,*) ""		  
           write(*,*) "================================================"
  !         tmpStr = new("    ") + new("USING") + new(" ")
  !   x   + intToStr(iNumDB) + new(" ATOMS FOR VERTICAL ENERGIES")
  !          write(*,*) tmpStr%data(1:tmpStr%length)
           write(*,1010) ,iNumDB
1010             format (" ",i4,
     x                 " ATOMS IN VERTICAL ENERGY CALCULATIONS")           
           write(*,*) ""		              
           write(*,*) "          atom index        Charge "
           write(*,*) "          ----------        ------ "           
            do iii = 1, iNumDB
            write(*,1001) iii, dVosDbase(iii)%atom, dVosDbase(iii)%q            
           enddo
           write(*,*) "================================================"    
           write(*,*) ""		             
           endif
           
c1001             format (5x,i3,5x,i5,11x,f6.3)           
1001             format (4x,i4,5x,i6,9x,f6.3)           
        return
       end subroutine sendDVosDbase

       subroutine sendDVosPolDbase()
c    *******************************************************************
c    ** send a copy of the dVosPolDbase from the master to all nodes  **
c    **                                                               **
c    *******************************************************************
       implicit none
       
       ! Locals
        integer :: ii, jj, iCount, iii, iNum
        real :: tmpQ
c    *******************************************************************


          ! Allocate the database
          if(rank .gt. master) then
           allocate(dVosPolDbase(iNumPolDB))
          endif

          ! ... then send the database
          call MPI_BCAST(dVosPolDbase,iNumPolDB,DVOSDB_TYPE,
     x                     master,world,ierr)

          ! allocate space for the original deltaQ charges
           allocate(origDeltaQ(iNumPolDB))
           allocate(origQ(iNumPolDB))
           
           ! Reset to become atomic indexes
           do iii = 1, iNumPolDB
            iCount = 0
           
            do ii = 1, dVosPolDbase(iii)%residue
             do jj = 1, params%numAtoms(ii)

              iCount = iCount + 1             
              if ( jj .eq. dVosPolDbase(iii)%atom ) then              
               iNum = iCount
              endif

             enddo ! end atom loop
            enddo ! end res loop
           
            dVosPolDbase(iii)%atom = iNum ! reset to atomic index
           enddo ! End database loop

           ! Make a copy of the original q's and deltaQ's
           do ii = 1, iNumPolDB
           origQ(ii) = qrms(dVosPolDbase(ii)%atom)%q
            do jj = 1, iNumDB
             if( dVosPolDbase(ii)%atom .eq. 
     x           dVosDbase(jj)%atom ) then
              origDeltaQ(ii) = dVosDbase(jj)%q
             endif
            enddo
           enddo
           
           if(rank .eq. master) then
           write(*,*) ""		  
           write(*,*) "================================================"
           tmpStr = new("    ") + new("USING") + new(" ")
     x   + intToStr(iNumPolDB) + new(" ATOMS FOR POLARIZABLE ANALYSIS")
            write(*,*) tmpStr%data(1:tmpStr%length)
           write(*,*) ""		              
           write(*,*) "          atom index        Charge "
           write(*,*) "          ----------        ------ "           
            do iii = 1, iNumPolDB
            write(*,801) iii, dVosPolDbase(iii)%atom, 
     x                   dVosPolDbase(iii)%q            
           enddo
           write(*,*) "================================================"    
           write(*,*) ""		             
           endif
           
801             format (5x,i3,5x,i5,11x,f6.3)           

        return
       end subroutine sendDVosPolDbase
       
       subroutine sendDADbase()
c    *******************************************************************
c    ** send a copy of the donor-acceptor from the master             **
c    ** to all nodes                                                  **
c    **                                                               **
c    *******************************************************************
       implicit none
       
       ! Locals
        integer :: ii, jj, iCount, iii, iNum
        integer :: iMax
        real :: tmpQ
c    *******************************************************************

          ! We already know the numbers of both from setMDParm ...
          if(rank .gt. master) then
           allocate(ariD(iNumRdaD))           
           allocate(ariA(iNumRdaA))          
          endif

          ! ... then send the donor acceptor databases
          call MPI_BCAST(ariA,iNumRdaA,DADB_TYPE,
     x                     master,world,ierr)
          call MPI_BCAST(ariD,iNumRdaD,DADB_TYPE,
     x                     master,world,ierr)

          ! reorganize the donors
           do iii = 1, iNumRdaD
            iCount = 0
           
            do ii = 1, ariD(iii)%iR
             do jj = 1, params%numAtoms(ii)

              iCount = iCount + 1             
              if ( jj .eq. ariD(iii)%iA ) then              
               iNum = iCount
              endif

             enddo ! end atom loop
            enddo ! end res loop
           
            ariD(iii)%iA = iNum ! reset to atomic index
           enddo ! End database loop
           
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          ! reorganize the acceptors
           do iii = 1, iNumRdaA
            iCount = 0
           
            do ii = 1, ariA(iii)%iR
             do jj = 1, params%numAtoms(ii)

              iCount = iCount + 1             
              if ( jj .eq. ariA(iii)%iA ) then              
               iNum = iCount
              endif

             enddo ! end atom loop
            enddo ! end res loop
           
            ariA(iii)%iA = iNum ! reset to atomic index
           enddo ! End database loop
           
           if(rank .eq. master) then
           
           ! Find the maximum extension
           if (iNumRdaA .gt. iNumRdaD ) then
            iMax = iNumRdaA
           elseif (iNumRdaD .gt. iNumRdaA ) then
            iMax = iNumRdaD
           else
            iMax = iNumRdaD                     
           endif
           
           write(*,*) ""		  
           write(*,*) "================================================"
           tmpStr = new("    ") + new("R(D-A) :: USING") + new(" ")
     x   + intToStr(iNumRdaD) + new(" ATOMS AS DONOR ATOMS")
            write(*,*) tmpStr%data(1:tmpStr%length)
           tmpStr = new("               ") + new("... and") + new(" ")
     x   + intToStr(iNumRdaA) + new(" AS ACCEPTORS")
            write(*,*) tmpStr%data(1:tmpStr%length)
            
           write(*,*) ""		              
           write(*,*) "          Donor index       Acceptor index "
           write(*,*) "          ------------      -------------- "           
            do iii = 1, iMax
            
           if ((iii .le. iNumRdaD) .and. 
     x         (iii .le. iNumRdaA)) then
            write(*,1114) iii, ariD(iii)%iA,  ariA(iii)%iA
            
           elseif(iii .le. iNumRdaD) then
            write(*,1113) iii, ariD(iii)%iA
           else
            write(*,1115) iii, ariA(iii)%iA            
           endif
           enddo
           write(*,*) "================================================"    
           write(*,*) ""		             
           endif
           
1113             format (5x,i3,6x,i5)  ! on the left    
1114             format (5x,i3,6x,i5,14x,i5) ! both
1115             format (5x,i3,25x,i5) ! on the right

        return
       end subroutine sendDADbase
       
       subroutine searchDbase(iAtom,iRes,dvosQ,foundAtom)
c    *******************************************************************
c    ** search the dVosDbase to tell if an atom is involved in the    **
c    ** stokes shift correlation function calculations.               **
c    **                                                               **
c    ** principal variables:                                          **
c    ** integer    iAtom   atom's index inside the residue            **
c    ** integer    iRes    residue's index inside the protein         **
c    ** real       dvosQ   (output) atomic q in other redox state     **
c    ** logical    foundAtom  (output) found atom in dVosDbase        **
c    **                                                               **
c    *******************************************************************
       implicit none

       ! Args
       integer :: iAtom, iRes
       real    :: dvosQ
       logical :: foundAtom
c    *******************************************************************

        foundAtom = .false.
        dvosQ = 0.d0
        do i = 1, iNumDB
         if(dVosDbase(i)%atom .eq. iAtom .and.
     x      dVosDbase(i)%residue .eq. iRes) then
            dvosQ = dVosDbase(i)%q
            foundAtom = .true.
            exit
         endif
        enddo

        return
       end subroutine searchDbase
        
       subroutine searchDBAtom(iAtom,dvosQ,foundAtom)
c    *******************************************************************
c    ** search the dVosDbase to tell if an atom is involved in the    **
c    ** stokes shift correlation function calculations.               **
c    **                                                               **
c    ** principal variables:                                          **
c    ** integer    iAtom   atom's index in nodeconf space             **
c    ** real       dvosQ   (output) atomic q in other redox state     **
c    ** logical    foundAtom  (output) found atom in dVosDbase        **
c    **                                                               **
c    *******************************************************************
       implicit none

       ! Args
       integer :: iAtom
       real    :: dvosQ
       logical :: foundAtom
c    *******************************************************************

        foundAtom = .false.
        dvosQ = 0.d0
        do i = 1, iNumDB
         if(dVosDbase(i)%atom .eq. iAtom) then
            dvosQ = dVosDbase(i)%q
            foundAtom = .true.
            exit
         endif
        enddo

        return
       end subroutine searchDBAtom
       
		subroutine allocConfs()
c    *******************************************************************
c    ** allocates nProc configurations on master 1 on the nodes       **
c    **                                                               **
c    ******************************************************************* 		
			implicit none
c    ******************************************************************* 			             
	     if (rank .eq. master) then
              allocate(buffConfs(natom*nProc))
              allocate(boxLs(nProc))
              if(doGrCalc) allocate(lowRXO(nProc))              
              if (doAllEner .or. doExclr) then
               if(doIons) then 
                 allocate(dv0inITraj(nProc))
                 allocate(dv0iTraj(nProc))  
               endif 

               if(doLipid) then 
                 allocate(dv0inLTraj(nProc))
                 allocate(dv0lTraj(nProc))  
               endif 

               if(doOthers) then 
                 allocate(dv0inOTraj(nProc))
                 allocate(dv0oTraj(nProc))  
               endif 

              if(doSolv .and. doProt) then
               allocate(dv0sTraj(nProc))
               allocate(dv0inSTraj(nProc))
               allocate(dv0pTraj(nProc))
               allocate(dv0inPTraj(nProc))               
              elseif(doSolv) then
               allocate(dv0sTraj(nProc))
               allocate(dv0inSTraj(nProc))
              else ! protein only              
               allocate(dv0pTraj(nProc))
               allocate(dv0inPTraj(nProc))
              endif ! end check for solv/prot
             
             elseif(doCoul .or. doSolut) then
              if(doSolv .and. doProt 
     x           .and. doIons) then             
               allocate(dv0sTraj(nProc))
               allocate(dv0iTraj(nProc))               
               allocate(dv0pTraj(nProc))               
              elseif(doSolv .and. doProt) then
               allocate(dv0sTraj(nProc))
               allocate(dv0pTraj(nProc))               
              elseif(doIons .and. doProt) then
               allocate(dv0iTraj(nProc))
               allocate(dv0pTraj(nProc))               
              elseif(doSolv .and. doIons) then
               allocate(dv0sTraj(nProc))
               allocate(dv0iTraj(nProc))                              
              elseif(doIons) then !ions only
               allocate(dv0iTraj(nProc))               
              elseif(doSolv) then !solvent only
               allocate(dv0sTraj(nProc))                              
              else !protein only
               allocate(dv0pTraj(nProc))               
              endif ! end check for solv/prot/ions     
          
              if(doOthers) then 
                allocate(dv0oTraj(nProc))
              endif ! end check for others

              if(doLipid) then 
                allocate(dv0lTraj(nProc))
              endif ! end check for others

             elseif(doInduced) then
              if(doSolv .and. doProt) then                     
               allocate(dv0inSTraj(nProc))  
               allocate(dv0inPTraj(nProc))                 
              elseif(doSolv) then
               allocate(dv0inSTraj(nProc))                
              else !protein only
               allocate(dv0inPTraj(nProc))  
              endif ! end check for solv/prot  
              if(doIons) then 
                allocate(dv0inITraj(nProc))
              endif ! end check for others
              if(doOthers) then 
                allocate(dv0inOTraj(nProc))
              endif ! end check for others
              if(doLipid) then 
                allocate(dv0inLTraj(nProc))
              endif ! end check for others
            
             endif ! End check for ener types
             
              if(doShell) then
               allocate(dv0bshTraj(nProc))
               allocate(dv0shTraj(nProc))
               allocate(dv0bdeshTraj(nProc))
               allocate(dv0shdeTraj(nProc))    
               allocate(nShellLowTraj(nProc))
               allocate(nShellHighTraj(nProc))               
        
              endif
              
             if(doCBind) then
!               allocate(depsdDepsTraj(nProc))
!               allocate(deps2dDepsTraj(nProc))
!               allocate(deps3dDepsTraj(nProc))
!               allocate(dDeps2Traj(nProc))
!               allocate(dDeps2depsTraj(nProc))
!               allocate(dDeps2deps2Traj(nProc))
               allocate(lamYiTraj(nProc))
               allocate(YiTraj(nProc))               
c               allocate(nShellTraj(nProc))
             endif

c             if(doShellM) then
c               allocate(nShellTraj(nProc))              
c             endif

              if(doProtV) then
                allocate(mcVolProtTraj(nProc))
                allocate(mcVolWatTraj(nProc))
                allocate(mcEpsProtTraj(nProc))
                allocate(mcEpsWatTraj(nProc))
              endif         

              if(doHoleV) then
                allocate(mcVolHoleTraj(nProc))
                allocate(mcEpsHoleTraj(nProc))
                allocate(holeCOMTraj(nConf))
              endif         


              if(doCluster) then
                allocate(aggListTraj(MAX_AGG*nProc))
                allocate(nAggTraj(nProc))

               if(doGyr) then
                 allocate(aggGyrTraj(MAX_AGG*nProc))

                 allocate(aggAspTraj(MAX_AGG*nProc))

                 allocate(aggR1Traj(MAX_AGG*nProc))
                 allocate(aggR2Traj(MAX_AGG*nProc))
                 allocate(aggR3Traj(MAX_AGG*nProc))
               endif

               if(doInert) then
                 allocate(aggEccenTraj(MAX_AGG*nProc))

                 allocate(aggA2Traj(MAX_AGG*nProc))
                 allocate(aggB2Traj(MAX_AGG*nProc))
                 allocate(aggC2Traj(MAX_AGG*nProc))
               endif

              endif

              if(doShellCalc) then
                allocate(nShellTraj(nProc))              
              endif
c             if(doShellM) then
c               allocate(nShellTraj(nProc))              
c             endif

c             if(doProtM) then
c               allocate(nShellTraj(nProc))              
c             endif

c             if(doShellP2) then
c               allocate(nShellTraj(nProc))              
c             endif
             if(doVPolCorr) allocate(dVpolTraj(nProc))
             if(doQPolCorr) allocate(dQpolTraj(nProc))
             if(doDQEfield) allocate(efield_dVTraj(nProc))
             if(doRdaCalc) allocate(rdaRs(nProc))
             if(doMdaCalc) allocate(mdaRs(nProc))
             if(doAvgS) then             

              if(doSelect) then
               allocate(xAvgTot(nSelect))
               allocate(yAvgTot(nSelect))
               allocate(zAvgTot(nSelect))
              else
               allocate(xAvgTot(iProtEnd))
               allocate(yAvgTot(iProtEnd))
               allocate(zAvgTot(iProtEnd))
              endif

             endif     
             if(doQdr) then
               allocate(drxAvgTot(iProtEnd))
               allocate(dryAvgTot(iProtEnd))
               allocate(drzAvgTot(iProtEnd))
               allocate(qdrsqTot(iProtEnd))
             endif
		    endif ! End check for master

            !!!!! below are allocaitons for ALL CPUS, including master !!!!


            if (doProtV .or. doHoleV) then
               allocate(ranX(nAttempts))
               allocate(ranY(nAttempts))
               allocate(ranZ(nAttempts))
            endif

              if(doCluster) then
                allocate(aggList(MAX_AGG))
                allocate(checkList(nClusterRes))
                allocate(neighborCount(nClusterRes))
                allocate(neighborList(nClusterRes,MAX_PER_CLUSTER))

                if(doGyr) then
                  allocate(aggGyr(MAX_AGG))

                  allocate(aggAsp(MAX_AGG))

                  allocate(aggR1(MAX_AGG))
                  allocate(aggR2(MAX_AGG))
                  allocate(aggR3(MAX_AGG))
                endif

                if(doInert) then
                  allocate(aggEccen(MAX_AGG))

                  allocate(aggA2(MAX_AGG))
                  allocate(aggB2(MAX_AGG))
                  allocate(aggC2(MAX_AGG))
                endif

              endif

            if(doExclr) then
               allocate(exclrDist(natom))
            endif
		    
            allocate(nodeConf(natom))
            allocate(boxLn(1))
            if(doRdaCalc) allocate(rdaRn(1))
            if(doMdaCalc) allocate(mdaRn(1))
            if(doAvgS .or. doCovar) then  ! Note: doQdr has its allocation in setAvgFile           
             if (doSelect) then
               allocate(xAvg(nSelect))
               allocate(yAvg(nSelect))
               allocate(zAvg(nSelect))
              else
               allocate(xAvg(iProtEnd))
               allocate(yAvg(iProtEnd))
               allocate(zAvg(iProtEnd))
             endif
            endif ! check for avg/covar   

             if(useCUDA) then
#ifdef ENABLE_CUDA
              ! allocate space on the GPU for program-wide storage

              doGPUInit = .true.
              if(doCoul) then
               call alloc_gpu_dvos(rank, iWatStart, natom,
     x                   GPUsPerNode,doGPUInit)
              endif

              if(doAvgS .or. doCovar) then
               call alloc_gpu_pca(rank, iProtStart,iProtEnd,
     x                qrms,doSelect, selectIndices, nSelect,
     x                                GPUsPerNode,doGPUInit)
              endif

              if(doShellCalc .or. doAvgS .or. doCovar) then
               call alloc_gpu_shell(rank, iWatStart, iProtEnd,
     x                     natom,iNumDB,nAttempts, nSelect,
     x                           GPUsPerNode,doGPUInit)
              endif

#endif
             endif ! End check for GPU settings
		    return
		 end subroutine allocConfs		

             subroutine deallocAll()
c    *******************************************************************
c    ** deallocates everything on master the nodes                    **
c    **                                                               **
c    ******************************************************************* 		
             implicit none
c    ******************************************************************* 		    
             if (rank .eq. master) then
              deallocate(buffConfs)
              deallocate(boxLs)              
              if(doGrCalc) deallocate(lowRXO)
             if (doAllEner .or. doExclr) then

               if(doIons) then 
                 deallocate(dv0inITraj)
                 deallocate(dv0iTraj)  
               endif 

               if(doLipid) then 
                 deallocate(dv0inLTraj)
                 deallocate(dv0lTraj)  
               endif 

               if(doOthers) then 
                 deallocate(dv0inOTraj)
                 deallocate(dv0oTraj)  
               endif 


              if(doSolv .and. doProt) then
               deallocate(dv0sTraj)
               deallocate(dv0inSTraj)
               deallocate(dv0pTraj)
               deallocate(dv0inPTraj)               
              elseif(doSolv) then
               deallocate(dv0sTraj)
               deallocate(dv0inSTraj)
              else !protein only
               deallocate(dv0pTraj)
               deallocate(dv0inPTraj)
              endif ! end check for solv/prot
             
             elseif(doCoul .or. doSolut) then
              if(doSolv .and. doProt 
     x                .and. doIons) then             
               deallocate(dv0sTraj)
               deallocate(dv0pTraj)
               deallocate(dv0iTraj)   
              elseif(doSolv .and. doProt) then             
               deallocate(dv0sTraj)
               deallocate(dv0pTraj)            
              elseif(doSolv .and. doIons) then             
               deallocate(dv0sTraj)
               deallocate(dv0iTraj)               
              elseif(doIons .and. doProt) then             
               deallocate(dv0iTraj)
               deallocate(dv0pTraj)               
              elseif(doSolv) then
               deallocate(dv0sTraj)
              elseif(doIons) then
               deallocate(dv0iTraj)              
              else ! protein only
               deallocate(dv0pTraj)
              endif ! end check for solv/prot/ions
              if(doOthers) then
               deallocate(dv0oTraj)
              endif              
              if(doLipid) then
               deallocate(dv0lTraj)
              endif              

             elseif(doInduced) then
              if(doSolv .and. doProt) then                     
               deallocate(dv0inPTraj)                
               deallocate(dv0inSTraj)
              elseif(doSolv) then
               deallocate(dv0inSTraj)              
              else !protein only
               deallocate(dv0inPTraj)  
              endif ! end check for solv/prot     
              if(doOthers) then
               deallocate(dv0inOTraj)
              endif              
              if(doLipid) then
               deallocate(dv0inLTraj)
              endif                       
              if(doIons) then
               deallocate(dv0inITraj)
              endif                       
             endif ! End check for ener types              
             
              if(doShell) then
                deallocate(dv0bshTraj)
                deallocate(dv0shTraj)
                deallocate(dv0bdeshTraj)
                deallocate(dv0shdeTraj)
                deallocate(nShellLowTraj)
                deallocate(nShellHighTraj)                  
              endif 
              
             if(doCBind) then
!               deallocate(depsdDepsTraj)
!               deallocate(deps2dDepsTraj)
!               deallocate(deps3dDepsTraj)
!               deallocate(dDeps2Traj)
!               deallocate(dDeps2depsTraj)
!               deallocate(dDeps2deps2Traj)
               allocate(lamYiTraj(nProc))
               allocate(YiTraj(nProc))               
c               deallocate(nShellTraj)
             endif

c             if(doShellM) then
c               deallocate(nShellTraj) 
c             endif

c             if(doProtM) then
c               deallocate(nShellTraj) 
c             endif

c             if(doShellP2) then
c               deallocate(nShellTraj) 
c             endif

              if(doProtV) then
                deallocate(mcVolProtTraj)
                deallocate(mcVolWatTraj)
                deallocate(mcEpsProtTraj)
                deallocate(mcEpsWatTraj)
              endif         

              if(doHoleV) then
                deallocate(mcVolHoleTraj)
                deallocate(mcEpsHoleTraj)
                deallocate(holeCOMTraj)
              endif         

             if(doShellCalc) then
               deallocate(nShellTraj) 
             endif

              if(doCluster) then
                deallocate(aggListTraj)
                deallocate(nAggTraj)

               if(doGyr) then
                 deallocate(aggGyrTraj)

                 deallocate(aggAspTraj)

                 deallocate(aggR1Traj)
                 deallocate(aggR2Traj)
                 deallocate(aggR3Traj)
               endif

               if(doInert) then
                 deallocate(aggEccenTraj)

                 deallocate(aggA2Traj)
                 deallocate(aggB2Traj)
                 deallocate(aggC2Traj)
               endif

              endif !end cluster
             
             if(doRdaCalc) deallocate(rdaRs)
             if(doMdaCalc) deallocate(mdaRs)             
             if(doIndCrg) deallocate(fFacTraj)
             if(doDQEfield) deallocate(efield_dVTraj)
             if(doVPolCorr) deallocate(dVpolTraj)
             if(doQPolCorr) deallocate(dQpolTraj)
             if(doAvgS) then             
              deallocate(xAvgTot)
              deallocate(yAvgTot)
              deallocate(zAvgTot)
             endif             
             if(doQdr) then
               deallocate(drxAvgTot)
               deallocate(dryAvgTot)
               deallocate(drzAvgTot)
               deallocate(qdrsqTot)
             endif             
           endif ! End check for master

            if(doExclr) deallocate(exclrDist)

            deallocate(nodeConf)
            deallocate(qrms)
            if(doEnerGap) then
             deallocate(dVosDbase)
            endif
            if(doInduced .or. doAllEner .or.
     x         doExclr) then
             deallocate(alphaDbase)
            endif
            
            if(doRdaCalc) then 
             deallocate(ariA) 
             deallocate(ariD)             
             deallocate(rdaRn)             
            endif            
            
            if(doMdaCalc) deallocate(mdaRn)

            if (doShell) then
             deallocate(aminoAcidList)
             deallocate(nShellHighAA)
             deallocate(nShellLowAA)          
             deallocate(ebindTotAA)
             deallocate(ebindDelAA)             
            endif
            
            if(doDIons) then
             deallocate(starts)
             deallocate(resNum)
            endif

            if(doAvgS .or. 
     x          doQdr .or. doCovar) then

              deallocate(xAvg)
              deallocate(yAvg)
              deallocate(zAvg)

              deallocate(proteinCOM)
              deallocate(proteinCRYS)

             endif          

            if(doAvgS .or. doQdr) then
              deallocate(proteinROT)
            endif

            if(doSelect) then
             deallocate(selectIndices)
            endif

            if(doCovar) then
              deallocate(vecAvg)
              deallocate(vecConf)
              deallocate(covarNode)
              deallocate(covarMat)
            endif

            if(doQdr) then
              deallocate(drxAvg)
              deallocate(dryAvg)
              deallocate(drzAvg)
              deallocate(qdrsq)
            endif        
            
            if(doIndCrg) then
             deallocate(origDeltaQ)
             deallocate(origQ)
             deallocate(dVosPolDbase)
            endif            

              if(doCluster) then
                deallocate(aggList)
                deallocate(checkList)
                deallocate(neighborCount)
                deallocate(neighborList)

                if(doGyr .or. doInert) then
                  deallocate(aggPositions)
                endif

                if(doGyr) then
                  deallocate(aggGyr)

                  deallocate(aggAsp)

                  deallocate(aggR1)
                  deallocate(aggR2)
                  deallocate(aggR3)
                endif

                if(doInert) then
                  deallocate(aggEccen)

                  deallocate(aggA2)
                  deallocate(aggB2)
                  deallocate(aggC2)
                endif

              endif

            if (doProtV .or. doHoleV) then
               deallocate(ranX)
               deallocate(ranY)
               deallocate(ranZ)
            endif
            
            deallocate(boxLn)

             if(useCUDA) then
#ifdef ENABLE_CUDA
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! deallocate space on the GPU
              if(doCoul) then
               call dealloc_gpu_dvos()
              endif

              if(doAvgS .or. doCovar) then
               call dealloc_gpu_pca()
              endif

              if(doShellCalc 
     x      .or. doAvgS .or. doCovar) then
               call dealloc_gpu_shell()
              endif
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
             endif ! End check for GPU settings
		    return
		 end subroutine deallocAll
		 
		subroutine readCRDMPI()
c    *******************************************************************
c    ** reads nProc configurations on master and distributes via      **
c    ** mpi_scatter 1 config to each node                             **
c    **                                                               **
c    ******************************************************************* 		
			implicit none
					    
		    ! Locals
		    integer :: ii, jj, iCnt,ir
            real :: box(3)
c    ******************************************************************* 

!+=============================================================+!
!           Read configs on master and scatter to nodeConf
!+=============================================================+!
            !!!!!!!!!! READ ACTUAL RECORD !!!!!!!!!!!!!!!!!!!!!!!!
            if (rank .eq. master) then    ! CHECK FOR MASTER !!!!!
     	    jj = 1
     	    iCnt = 0
     	    ir = 0
     	    do ii = 1, nProc !!!!!!!!!!!!!!!!!!!!!!!!     	    
     	    
             read(idat,'(10F8.3)') 
     x           (buffConfs(ir+i)%x,buffConfs(ir+i)%y,buffConfs(ir+i)%z,
     x           i=1,natom)
             read(idat,'(10F8.3)') box(1), box(2), box(3)  !!!!!!!!!!!!!!     

             ir = ii*natom
             
            enddo                        !!!!!!!!!!!!!!!!!!!!!!!!
            endif                        ! END CHECK FOR MASTER !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            ! Distribute a configuration to each node
            call MPI_SCATTER(buffConfs,natom,ATOM_TYPE,
     x                       nodeConf,natom,ATOM_TYPE,
     x                       master,world,ierr)     
            !if (sendBoxL) then              
              call MPI_BCAST(box,3,MPI_REAL,
     x                    master,world,ierr)              
     
             bX = box(1)
             bY = box(2)
             bZ = box(3)  

             boxLn(1)%x = bX !box(1)
             boxLn(1)%y = bY !box(2)
             boxLn(1)%z = bZ !box(3)        
             
            ! sendBoxL = .false.
            !endif
            
		    return
		 end subroutine readCRDMPI	
		 
		subroutine writeHeadMPIIO()
c    *******************************************************************
c    ** reads the header on all nodes                                 **
c    **                                                               **
c    ******************************************************************* 		
			implicit none
			
			! Locals

           integer :: mpi_charExt, mpi_realExt, mpi_intExt			
c			character :: mstr(4)
           character(len=4)  :: mstr
c    *******************************************************************
       call MPI_BCAST(outDatFile,80
     x      ,MPI_CHARACTER,master,world,ierr) ! BCAST the filename       

c ********* UNCOMMENT THE FOLLOWING SECTION TO ENABLE MPI-IO ****************
c      print *, rank, " OPENING OUTPUT PFILE: ", outDatFile
           
       call MPI_FILE_OPEN(world,outDatFile,
     x      MPI_MODE_CREATE + MPI_MODE_RDWR, 
     x      MPI_INFO_NULL,        
     x      iOutDat,ierr)

           
       offset = 0
       mstr = '1234'

       call MPI_FILE_WRITE_AT_ALL(iOutDat,
     x      offset,
     x      mstr,4,
     x      MPI_CHARACTER,     
     x      MPI_STATUS_IGNORE,ierr)

c POSSIBLY ADD ERROR CHECKS HERE TO MAKE SURE EVERYONE HAS IT RIGHT? 
          if(rank .eq. master) then                   
          write( *,*) " WROTE MAGIC NUMBER: |",mstr,"|" !char(mstr)
        endif

       call MPI_TYPE_EXTENT(MPI_CHARACTER,mpi_charExt,ierr)
       call MPI_TYPE_EXTENT(MPI_INTEGER,mpi_intExt,ierr)
       call MPI_TYPE_EXTENT(MPI_REAL,mpi_realExt,ierr)       
c ********* UNCOMMENT THE PREVIOUS SECTION TO ENABLE MPI-IO ****************

       atomOff = mpi_realExt*3
       intOff = mpi_intExt

       offset = 4*mpi_charExt +
     x    rank * (intOff+(1+natom)*atomOff)

       !call MPI_FILE_SYNC(iOutDat,ierr)
       call MPI_BARRIER(world,ierr)

        return
      end subroutine writeHeadMPIIO
      
      subroutine writeBINMPIIO()
c    *******************************************************************
c    ** reads nProc configurations on master and distributes via      **
c    ** mpi_scatter 1 config to each node                             **
c    **                                                               **
c    ******************************************************************* 		
			implicit none
					    
		    ! Locals
            integer :: ii, jj, iCnt,ir, natm
            real :: box(3)
c    ********************************************************************
!+=============================================================+!
!           Read configs on all nodes using explicit offsets
!+=============================================================+!
            

c ********* UNCOMMENT THE FOLLOWING SECTION TO ENABLE MPI-IO ****************
c
c    NOTE: After each MPI write, the offset must be updated to the correct
c          file position by simply adding the offset of what was written  
c          to the offset variable.  And, after writing the whole configuration,
c          the offset is pushed down another nProc1, not just nProc.
c
c****************************************************************************      

            call MPI_FILE_WRITE_AT_ALL(iOutDat,
     x      offset,
     x      natom,1,
     x      MPI_INTEGER,     
     x      MPI_STATUS_IGNORE,ierr)

            offset = offset + intOff
            
            call MPI_FILE_WRITE_AT_ALL(iOutDat,
     x      offset,
     x      nodeConf,natom,
     x      ATOM_TYPE,     
     x      MPI_STATUS_IGNORE,ierr)

            offset = offset + atomOff*natom ! For box(x,y,z)
            call MPI_FILE_WRITE_AT_ALL(iOutDat,
     x      offset,
     x      boxLn,1,
     x      BOXL_TYPE,     
     x      MPI_STATUS_IGNORE,ierr)            
            
c ********* UNCOMMENT THE PREVIOUS SECTION TO ENABLE MPI-IO ****************
            offset = offset + atomOff +
     x         nProc1*(atomOff*(natom+1)+intOff)          

            !call MPI_FILE_SYNC(iOutDat,ierr)
            call MPI_BARRIER(world,ierr)

         
		    return
		 end subroutine writeBINMPIIO	      
		 
		subroutine readHeadMPIIO()
c    *******************************************************************
c    ** reads the header on all nodes                                 **
c    **                                                               **
c    ******************************************************************* 		
			implicit none
			
			! Locals

           integer :: mpi_charExt, mpi_realExt, mpi_intExt			
!			character :: mstr(4)
			character(len=4) :: mstr
c    *******************************************************************
       call MPI_BCAST(datFile,80
     x      ,MPI_CHARACTER,master,world,ierr) ! BCAST the filename       

c ********* UNCOMMENT THE FOLLOWING SECTION TO ENABLE MPI-IO ****************
       call MPI_FILE_OPEN(world,trim(datFile),
     x      MPI_MODE_RDONLY, 
     x      MPI_INFO_NULL,        
     x      idat,ierr)

       offset = 0

       call MPI_FILE_READ_AT_ALL(idat,
     x      offset,
     x      mstr,4,
     x      MPI_CHARACTER,     
     x      MPI_STATUS_IGNORE,ierr)
     
c POSSIBLY ADD ERROR CHECKS HERE TO MAKE SURE EVERYONE HAS IT RIGHT? 
        if (rank .eq. master) then
          write( *,*) " READ MAGIC NUMBER: |",mstr,"|" !char(mstr)
        endif

       call MPI_TYPE_EXTENT(MPI_CHARACTER,mpi_charExt,ierr)
       call MPI_TYPE_EXTENT(MPI_INTEGER,mpi_intExt,ierr)
       call MPI_TYPE_EXTENT(MPI_REAL,mpi_realExt,ierr)       
c ********* UNCOMMENT THE PREVIOUS SECTION TO ENABLE MPI-IO ****************

       atomOff = mpi_realExt*3
       intOff = mpi_intExt

       offset = 4*mpi_charExt +
     x    rank * (intOff+(1+natom)*atomOff)

       call MPI_BARRIER(world,ierr)

        return
      end subroutine readHeadMPIIO
        
		subroutine readBINMPIIO()
c    *******************************************************************
c    ** reads nProc configurations on master and distributes via      **
c    ** mpi_scatter 1 config to each node                             **
c    **                                                               **
c    ******************************************************************* 		
			implicit none
					    
		    ! Locals
		    integer :: ii, jj, iCnt,ir, natm
            real :: box(3)
c    ********************************************************************
!+=============================================================+!
!           Read configs on all nodes using explicit offsets
!+=============================================================+!
            

c ********* UNCOMMENT THE FOLLOWING SECTION TO ENABLE MPI-IO ****************
c
c    NOTE: After each MPI read, the offset must be updated to the correct
c          file position simply adding the offset of what was read to the 
c          offset variable.  And, after reading the whole configuration,
c          the offset is pushed down another nProc1, not just nProc.
c
c****************************************************************************
            call MPI_FILE_READ_AT_ALL(idat,
     x      offset,
     x      natom,1,
     x      MPI_INTEGER,     
     x      MPI_STATUS_IGNORE,ierr)

            offset = offset + intOff
            
            call MPI_FILE_READ_AT_ALL(idat,
     x      offset,
     x      nodeConf,natom,
     x      ATOM_TYPE,     
     x      MPI_STATUS_IGNORE,ierr)

            offset = offset + atomOff*natom     
            
            call MPI_FILE_READ_AT_ALL(idat,
     x      offset,
     x      boxLn,1,
     x      BOXL_TYPE,     
     x      MPI_STATUS_IGNORE,ierr)
                   
c ********* UNCOMMENT THE PREVIOUS SECTION TO ENABLE MPI-IO ****************
            offset = offset + atomOff +           
     x         nProc1*(atomOff*(natom+1)+intOff)            
     
             bX = boxLn(1)%x !box(1)
             bY = boxLn(1)%y!box(2)
             bZ = boxLn(1)%z!box(3)               
             box(1) = bX
             box(2) = bY
             box(3) = bZ
             
            call MPI_BARRIER(world,ierr)
       
		    return
		 end subroutine readBINMPIIO	

		subroutine readBINMMPI()
c    *******************************************************************
c    ** reads nProc configurations on master and distributes via      **
c    ** mpi_scatter 1 config to each node using a modified binpos     **
c    ** file type.                                                    **
c    **                                                               **
c    ******************************************************************* 		
			implicit none
					    
		    ! Locals
		    integer :: ii, jj, iCnt,ir, natm
            real :: box(3)
c    ******************************************************************* 

!+=============================================================+!
!           Read configs on master and scatter to nodeConf
!+=============================================================+!
            !!!!!!!!!! READ ACTUAL RECORD !!!!!!!!!!!!!!!!!!!!!!!!
            if (rank .eq. master) then    ! CHECK FOR MASTER !!!!!
     	    jj = 1
     	    iCnt = 0
     	    ir = 0
     	    do ii = 1, nProc !!!!!!!!!!!!!!!!!!!!!!!!     	    

             read(idat) natm  !!!!!!!!!!!!!!     	    
             read(idat) 
     x           (buffConfs(ir+i)%x,buffConfs(ir+i)%y,buffConfs(ir+i)%z,
     x           i=1,natom)
!             read(idat) box(1), box(2), box(3)  !!!!!!!!!!!!!
             read(idat) boxLs(ii)%x,boxLs(ii)%y,boxLs(ii)%z

             ir = ii*natom
            
            enddo
            
                                         !!!!!!!!!!!!!!!!!!!!!!!!
            endif                        ! END CHECK FOR MASTER !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            ! Distribute a configuration and a box to each node
            call MPI_SCATTER(buffConfs,natom,ATOM_TYPE,
     x                       nodeConf,natom,ATOM_TYPE,
     x                       master,world,ierr)    
            call MPI_SCATTER(boxLs,1,BOXL_TYPE,
     x                       boxLn,1,BOXL_TYPE,
     x                       master,world,ierr)          
            !if (sendBoxL) then
!              call MPI_BCAST(box,3,MPI_REAL,
!     x                    master,world,ierr)              
     
!             bX = box(1)
!             bY = box(2)
!             bZ = box(3)  

             !bX = 52.955!boxLn(1)%x !box(1)
             !bY = 65.338!boxLn(1)%y!box(2)
             !bZ = 55.402!boxLn(1)%z!box(3)  
             bX = boxLn(1)%x !box(1)
             bY = boxLn(1)%y!box(2)
             bZ = boxLn(1)%z!box(3)               
             box(1) = bX
             box(2) = bY
             box(3) = bZ
            ! sendBoxL = .false.
            !endif
            
		    return
		 end subroutine readBINMMPI

		subroutine readGZMPI()
c    *******************************************************************
c    ** reads nProc configurations on master and distributes via      **
c    ** mpi_scatter 1 config to each node using a mdcrd.gz            **
c    ** file type.                                                    **
c    **                                                               **
c    ******************************************************************* 		
			implicit none
					    
		    ! Locals
		    integer :: ii, jj, iCnt,ir, natm
            real :: box(3)
c    ******************************************************************* 

!+=============================================================+!
!           Read configs on master and scatter to nodeConf
!+=============================================================+!
            !!!!!!!!!! READ ACTUAL RECORD !!!!!!!!!!!!!!!!!!!!!!!!
            if (rank .eq. master) then    ! CHECK FOR MASTER !!!!!            
              call read_gz(nProc, natom, buffConfs, boxLs)
            
                                         !!!!!!!!!!!!!!!!!!!!!!!!
            endif                        ! END CHECK FOR MASTER !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call MPI_BARRIER(world,ierr)
            ! Distribute a configuration and a box to each node
            call MPI_SCATTER(buffConfs,natom,ATOM_TYPE,
     x                       nodeConf,natom,ATOM_TYPE,
     x                       master,world,ierr)    
            call MPI_SCATTER(boxLs,1,BOXL_TYPE,
     x                       boxLn,1,BOXL_TYPE,
     x                       master,world,ierr)          
             bX = boxLn(1)%x !box(1)
             bY = boxLn(1)%y!box(2)
             bZ = boxLn(1)%z!box(3)               


c            print *, rank, "node...boxes: ", 
c     x         bX, bY, bZ

c            print *, rank, "node...buff: ", 
c     x         buffConfs(1)%x,buffConfs(1)%y,buffConfs(1)%z

             box(1) = bX
             box(2) = bY
             box(3) = bZ
            
		    return
		 end subroutine readGZMPI

             subroutine readDCDMPI()
c    *******************************************************************
c    ** reads nProc configurations on master and distributes via      **
c    ** mpi_scatter 1 config to each node using a dcd file type.      **
c    **                                                               **
c    ******************************************************************* 		
               implicit none
					    
		    ! Locals
		    integer :: ii, jj, iCnt,ir, natm
            real :: box(3)
c    ******************************************************************* 

!+=============================================================+!
!           Read configs on master and scatter to nodeConf
!+=============================================================+!
            !!!!!!!!!! READ ACTUAL RECORD !!!!!!!!!!!!!!!!!!!!!!!!
            if (rank .eq. master) then    ! CHECK FOR MASTER !!!!!            
              call read_dcd(nProc, natom, buffConfs, boxLs)
            
                                         !!!!!!!!!!!!!!!!!!!!!!!!
            endif                        ! END CHECK FOR MASTER !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call MPI_BARRIER(world,ierr)
            ! Distribute a configuration and a box to each node
            call MPI_SCATTER(buffConfs,natom,ATOM_TYPE,
     x                       nodeConf,natom,ATOM_TYPE,
     x                       master,world,ierr)    
            call MPI_SCATTER(boxLs,1,BOXL_TYPE,
     x                       boxLn,1,BOXL_TYPE,
     x                       master,world,ierr)          
             bX = boxLn(1)%x !box(1)
             bY = boxLn(1)%y!box(2)
             bZ = boxLn(1)%z!box(3)               


c            print *, rank, "node: ", rank, " boxes(x,y,z) = ",
c     x         bX, bY, bZ

c            print *, rank, "node ", rank, " buff(x,y,z) = ",
c     x         nodeConf(1)%x,nodeConf(1)%y,nodeConf(1)%z

             box(1) = bX
             box(2) = bY
             box(3) = bZ
            
		    return
		 end subroutine readDCDMPI
		 
		subroutine readBINMPI()
c    *******************************************************************
c    ** reads nProc configurations on master and distributes via      **
c    ** mpi_scatter 1 config to each node                             **
c    **                                                               **
c    ******************************************************************* 		
			implicit none
					    
		    ! Locals
		    integer :: ii, jj, iCnt,ir, natm
            real :: box(3)
c    ******************************************************************* 

!+=============================================================+!
!           Read configs on master and scatter to nodeConf
!+=============================================================+!
            !!!!!!!!!! READ ACTUAL RECORD !!!!!!!!!!!!!!!!!!!!!!!!
            if (rank .eq. master) then    ! CHECK FOR MASTER !!!!!
     	    jj = 1
     	    iCnt = 0
     	    ir = 0
     	    do ii = 1, nProc !!!!!!!!!!!!!!!!!!!!!!!!     	    

             read(idat) natm  !!!!!!!!!!!!!!     	    
             read(idat) 
     x           (buffConfs(ir+i)%x,buffConfs(ir+i)%y,buffConfs(ir+i)%z,
     x           i=1,natom)


             ir = ii*natom
            
            enddo
            
                                         !!!!!!!!!!!!!!!!!!!!!!!!
            endif                        ! END CHECK FOR MASTER !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            ! Distribute a configuration to each node
            call MPI_SCATTER(buffConfs,natom,ATOM_TYPE,
     x                       nodeConf,natom,ATOM_TYPE,
     x                       master,world,ierr)     
            if (sendBoxL) then
              call MPI_BCAST(box,3,MPI_REAL,
     x                    master,world,ierr)              
     
             bX = box(1)
             bY = box(2)
             bZ = box(3)  
             
             sendBoxL = .false.
            endif
            
		    return
		 end subroutine readBINMPI	

        subroutine closeGZ()
c    *******************************************************************
c    ** a fortran wrapper for closing the pre-opened mdcrd.gz file    **
c    **                                                               **
c    *******************************************************************          
         implicit none

         ! Locals         
c    *******************************************************************

         call close_gz()
           
         return
        end subroutine closeGZ

        subroutine closeDCD()
c    *******************************************************************
c    ** a fortran wrapper for closing the pre-opened binary DCD file  **
c    **                                                               **
c    *******************************************************************          
         implicit none

         ! Locals         
c    *******************************************************************

         call close_dcd()
           
         return
        end subroutine closeDCD

        subroutine setPolCharges()
c    *******************************************************************
c    ** sends out current f-parameters and resets the charges on      **
c    ** a polarizable center.                                         **
c    **                                                               **
c    *******************************************************************          
         implicit none

         ! Locals         
         integer :: ii, jj, iPolAtom
         real    :: fFac, rPolQ
c    *******************************************************************
         ! Give out the f-parameter corresponding to the current configuration
         call MPI_SCATTER(fFacTraj(iConf),1,MPI_REAL,
     x                    fFac,1,MPI_REAL,
     x                    master,world,ierr)  

         ! Loop over all polarizable atoms and setup the right deltaQ charge
           do ii = 1, iNumPolDB
           
            iPolAtom = dVosPolDbase(ii)%atom
            rPolQ    = dVosPolDbase(ii)%q
            
            qrms(iPolAtom)%q = origQ(ii) + 
     x      fFac*rPolQ
     
            do jj = 1, iNumDB
             if( iPolAtom .eq. 
     x           dVosDbase(jj)%atom ) then

c ????? POSSSIBLE CORRECT ?????       
!              dVosDbase(jj)%q = origDeltaQ(ii) - 
!     x        fFac*rPolQ
c ????? POSSSIBLE CORRECT ?????       
              dVosDbase(jj)%q = origDeltaQ(ii) - 
     x        fFac*rPolQ

             endif
            enddo
           enddo
           
         return
        end subroutine setPolCharges

        subroutine checkSelection()
c    *******************************************************************
c    **                                                               **
c    *******************************************************************          
         implicit none

         ! Locals         
         integer :: irc, iAtm
         real    :: rcQ1, delQ, rcQ2
c    *******************************************************************


         if(.not. doSelect) then
            ! Here we need one big boolean for all
            ! required functions that force the usage
            ! of atom selections    
            if(doCovar .or. doAvgS) then

             write(*,*) ""
             write(*,*) "  !!!!!!!!!!!!!!!!!!!!",
     x                  "!!!!!!!!!!!!!!!!!!!"
             write(*,*) "  !!! WARNING WARNING ",
     x                  "WARNING WARNING !!!"
             write(*,*) "  !!!                 ",
     x                  "                !!!"
            write(*,*) "  !!!  ASSUMING AN ALL PROTEIN ATOM   !!!"
            write (*,*)"  !!! SELECTION: GOODBYE FREE MEMORY  !!!"
             write(*,*) "  !!!                 ",
     x                  "                !!!"
             write(*,*) "  !!! WARNING WARNING ",
     x                  "WARNING WARNING !!!"
             write(*,*) "  !!!!!!!!!!!!!!!!!!!!",
     x                  "!!!!!!!!!!!!!!!!!!!"


            end if

         end if
         
         return
        end subroutine checkSelection
        
        subroutine outputCTState()
c    *******************************************************************
c    ** makes an input file for the othe CT state from the dv0s       **
c    ** database and the current charges. Used (only?) for preparing  **
c    ** a simulation.                                                 **
c    ** NOTE: To be careful, it forces master only, and synchs.       **
c    **                                                               **
c    *******************************************************************          
         implicit none

         ! Locals         
         integer :: irc, iAtm
         real    :: rcQ1, delQ, rcQ2
c    *******************************************************************

         if(rank .eq. master) then
         open(unit=58,file='afterCT.dat', status='unknown')
         do irc = 1, iNumDB  ! Loop over all RC atoms     
          iAtm = dVosDbase(irc)%atom          
           rcQ1 = qrms(iAtm)%q ! Initial state
           delQ = dVosDbase(irc)%q ! delta state
           rcQ2 = rcQ1 + delQ ! Final state
         write(58,*) iAtm, rcQ1, delQ, rcQ2           
         enddo                  
         close(58)
         
         endif ! End check for master
         
         call MPI_BARRIER(world,ierr)
         
         return
        end subroutine outputCTState

        subroutine outputSolvMol()
c    *******************************************************************
c    ** makes the three required input files for solvMol:             **
c    **    x,y,z,r file  **
c    **    x,y,z,q/2 file                                             **
c    **    x,y,z,deltaq file                                          **
c    ** NOTE: To be careful, it forces master only, and synchs.       **
c    **                                                               **
c    *******************************************************************          
         implicit none

         ! Locals         
         integer :: irc, iAtm, ipr
         real    :: rcQ1, delQ, rcQ2, rcMid
         real    :: qq, rr, mm
c    *******************************************************************

         if(rank .eq. master) then
         open(unit=58,file='diff.q1-dq-q2-qmid.dat', status='unknown')
         do irc = 1, iNumDB  ! Loop over all RC atoms     
          iAtm = dVosDbase(irc)%atom          
           rcQ1 = qrms(iAtm)%q ! Initial state
           delQ = dVosDbase(irc)%q ! delta state
           rcQ2 = rcQ1 + delQ ! Final state
           rcMid = (rcQ1+rcQ2)/2. !mid-point state
         write(58,*) iAtm, rcQ1, delQ, rcQ2, rcMid          
         enddo                  
         close(58)
         
         open(unit=58,file='r-q-m.dat', status='unknown')
        do ipr = 1, iProtEnd ! Loop over all protein atoms
         rr = qrms(ipr)%r
         qq = qrms(ipr)%q
         mm = qrms(ipr)%m         
         write(58,*) ipr, rr, qq, mm
         enddo                  
         close(58)
         
         endif ! End check for master
         
         call MPI_BARRIER(world,ierr)
         
         return
        end subroutine outputSolvMol
        
        subroutine readHeader()
c    *******************************************************************
c    ** reads header from the configuration data file                 **
c    **                                                               **
c    *******************************************************************          
         implicit none

         ! Locals         
!         character :: mstr(4)
          character(len=4) :: mstr
          character(len=80) :: psfFile          
c    *******************************************************************
         if (isBinary) then
!!!!!!!!!!!!!! IF THE FILE IS BINARY !!!!!!!!!!!!!!!!!!!!
          if ((.not. doMPIIO) .and. (.not. doGZRead) 
     x                        .and. (.not. doDCDRead)) then          
           ! Open the binpos file on master

#ifdef ENABLE_PGI
           ! use below for pgi comiles
           open(unit=idat,file=trim(datFile),
     x         status='old', form='binary')
#endif

#ifdef ENABLE_GFORT
           ! use the one below for gfortran compiles  
           open(unit=idat,file=trim(datFile),
     x          status='old', access='stream')  
#endif

#ifdef ENABLE_IFORT
           !wtf is access=transparent?... is this for ifort?
            open(unit=idat,file=trim(datFile),
     x          status='old',form='unformatted',access='transparent')  
#endif
           !read(idat) titleStr(1:4) ! Read the magic number 
           read(idat) mstr
           write( *,*) " MAGIC NUMBER: |",mstr,"|" !char(mstr)                      
          !!!!! READ A GZMD file
          elseif((.not. doMPIIO) .and. doGZRead) then          
            call open_gz(trim(datFile) // CHAR(0)) ! Append a char(0) so c knows the end of the string   
          !!!!! READ A DCD file
          elseif((.not. doMPIIO) .and. doDCDRead) then
            write(*,*) ""		  
            !tmpStr = all("=",55)
            tmpStr = new("=========================")
            tmpStr = tmpStr + new("=========================")
            tmpStr = tmpStr + new("=====")
            write(*,*) tmpStr%data(1:tmpStr%length)
            tmpStr = new(" ") + new("DCD HEADER INFO:")
            write(*,*) tmpStr%data(1:tmpStr%length)
c            write(*,*) ""		  
            call open_dcd(trim(datFile) // CHAR(0)) ! Append a char(0) so c knows the end of the string  
!            tmpStr = all("=",55)
            tmpStr = new("=========================")
            tmpStr = tmpStr + new("=========================")
            tmpStr = tmpStr + new("=====")
            write(*,*) tmpStr%data(1:tmpStr%length)
            write(*,*) ""		  

          endif
!!!!!!!!!!!!!! ELSE THE FILE IS NOT BINARY !!!!!!!!!!!!!!!!!!!!          
         else         
          ! Open the mdcrd file on master
          open(unit=idat,file=trim(datFile),
     x         status='old', form='formatted')  
          read(idat,'(a)') titleStr ! Read the title
         endif
		 
         return         
        end subroutine readHeader

           subroutine fillRecip()
c    *******************************************************************
c    ** gets reciprocal box info from geometric parameters            **
c    **                                                               **
c    *******************************************************************     
           implicit none
           
           ! Locals
           integer jj
           double precision u12(3), u23(3), u31(3)
           double precision factor, a, b, c, volume
c    *******************************************************************
            factor = 3.14159265358979323846d0/180.d0
            a = bX
            b = bY
            c = bZ

            ucell(1,1) = a
            ucell(2,1) = 0.d0
            ucell(3,1) = 0.d0
            ucell(1,2) = b*cos(factor*gamma)
            ucell(2,2) = b*sin(factor*gamma)
            ucell(3,2) = 0.d0
            ucell(1,3) = c*cos(factor*beta)
            ucell(2,3) = (b*c*cos(factor*alpha)-ucell(1,3)*ucell(1,2))/
     x                   ucell(2,2)
            ucell(3,3) = sqrt( c*c - ucell(1,3)*ucell(1,3) -
     x                    ucell(2,3)*ucell(2,3) )


C     now get reciprocal vectors
       call cross(ucell(1,2),ucell(1,3),u23)
       call cross(ucell(1,3),ucell(1,1),u31)
       call cross(ucell(1,1),ucell(1,2),u12)
       call dot(ucell(1,1),u23,volume)
       volCell = volume
       do jj = 1,3
        recip(jj,1) = u23(jj)/volume
        recip(jj,2) = u31(jj)/volume
        recip(jj,3) = u12(jj)/volume
       end do

       return
      end subroutine fillRecip
           
      subroutine dot(v1,v2,result)
c    *******************************************************************
c    ** calculates the scalar product on double precision vectors     **
c    **                                                               **
c    *******************************************************************  
       implicit none

       double precision v1(3),v2(3),result
c    *******************************************************************       
      result = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
      return
      end subroutine dot



      subroutine cross(v1,v2,v12)
c    *******************************************************************
c    ** calculates the cross product on double precision vectors      **
c    **                                                               **
c    *******************************************************************  
       implicit none

       double precision v1(3),v2(3),v12(3)
c    *******************************************************************       
       v12(1) = v1(2)*v2(3)-v1(3)*v2(2)
       v12(2) = v1(3)*v2(1)-v1(1)*v2(3)
       v12(3) = v1(1)*v2(2)-v1(2)*v2(1)
      return
      end subroutine cross
        
        subroutine reimageMPI()
c    *******************************************************************
c    ** reimages configuration data on all nodes                      **
c    **                                                               **
c    *******************************************************************  
           implicit none

           
           ! Locals
           integer :: ii
           real    :: rX, rY, rZ

      ! OMP section         
       integer :: OMP_GET_THREAD_NUM           
c    *******************************************************************

!Everyone reimages their own copy of the protein
!$OMP PARALLEL DO 
!$OMP&DEFAULT(PRIVATE), 
!$OMP&SHARED(nodeConf,bX,bY,bZ,natom)
            do ii = 1, natom
              
              rX = nodeConf(ii)%x
              rY = nodeConf(ii)%y
              rZ = nodeConf(ii)%z

!$OMP ATOMIC                                                         
              rX = rX - bX * ANINT(rX/bX)
!$OMP ATOMIC                                                                       
              rY = rY - bY * ANINT(rY/bY)
!$OMP ATOMIC                                                                       
              rZ = rZ - bZ * ANINT(rZ/bZ)

              nodeConf(ii)%x = rX
              nodeConf(ii)%y = rY
              nodeConf(ii)%z = rZ                            

            enddo
!$OMP END PARALLEL DO

             return
          end subroutine reimageMPI		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! Analysis Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcRDA()
c    *******************************************************************
c    ** calculates the positions and distance between the center of   **
c    ** of the donor and acceptor atomic databases.                   **
c    **                                                               **
c    *******************************************************************            
        implicit none
 
        ! Locals        
        integer :: iid,iia, iWrite
        integer :: iiatm
        real    :: massTot, mi
        real    :: ddX, ddY, ddZ
c    ******************************************************************* 

        
        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then 
          open(unit=iRDA,file='rDA.dat', status='unknown')                   
         endif
        endif 

        rdaRn(1)%xD = 0.d0 
        rdaRn(1)%yD = 0.d0 
        rdaRn(1)%zD = 0.d0 

        rdaRn(1)%xA = 0.d0 
        rdaRn(1)%yA = 0.d0 
        rdaRn(1)%zA = 0.d0       
       
        massTot = 0.d0       
       do iid = 1, iNumRdaD
        iiatm =  ariD(iid)%iA ! index in nodeConf space
        mi = qrms(iiatm)%m
        massTot = massTot + mi
        
        rdaRn(1)%xD = rdaRn(1)%xD + mi*nodeConf(iiatm)%x
        rdaRn(1)%yD = rdaRn(1)%yD + mi*nodeConf(iiatm)%y
        rdaRn(1)%zD = rdaRn(1)%zD + mi*nodeConf(iiatm)%z
      
       enddo
       rdaRn(1)%xD = rdaRn(1)%xD/massTot
       rdaRn(1)%yD = rdaRn(1)%yD/massTot
       rdaRn(1)%zD = rdaRn(1)%zD/massTot


       massTot = 0.d0       
       do iia = 1, iNumRdaA
        iiatm =  ariA(iia)%iA ! index in nodeConf space
        mi = qrms(iiatm)%m
        massTot = massTot + mi
        
        rdaRn(1)%xA = rdaRn(1)%xA + mi*nodeConf(iiatm)%x
        rdaRn(1)%yA = rdaRn(1)%yA + mi*nodeConf(iiatm)%y
        rdaRn(1)%zA = rdaRn(1)%zA + mi*nodeConf(iiatm)%z
      
       enddo
       rdaRn(1)%xA = rdaRn(1)%xA/massTot
       rdaRn(1)%yA = rdaRn(1)%yA/massTot
       rdaRn(1)%zA = rdaRn(1)%zA/massTot
       
       ddX = rdaRn(1)%xD - rdaRn(1)%xA
       ddY = rdaRn(1)%yD - rdaRn(1)%yA
       ddZ = rdaRn(1)%zD - rdaRn(1)%zA
       
        rdaRn(1)%rDA = sqrt(ddX**2+ddY**2+ddZ**2)

        call MPI_GATHER(rdaRn,1,RDAR_TYPE,
     x                 rdaRs,1,RDAR_TYPE,
     x                 master,world,ierr)     
       
           if(rank .eq. master)  then
           
            do iWrite = 1, nProc
              write(iRDA,3331) iConf+iWrite-1, 
     x         rdaRs(iWrite)%xD, rdaRs(iWrite)%yD,
     x         rdaRs(iWrite)%zD, rdaRs(iWrite)%xA,
     x         rdaRs(iWrite)%yA, rdaRs(iWrite)%zA,
     x         rdaRs(iWrite)%rDA   
            enddo
           endif       

        !call MPI_BARRIER(world,ierr)
        
3331             format (3x,i6,5x,7(F6.3,4x))

        return
      end subroutine calcRDA

      subroutine calcMDA()
c    *******************************************************************
c    ** calculates the x,y,z components and magnitude of the differ-  **
c    ** ence dipole moment.                                           **
c    **                                                               **
c    *******************************************************************            
        implicit none
 
        ! Locals        
        integer :: irc, iAtm, iWrite
        real    :: qTot, qi
        real    :: ddX, ddY, ddZ
c    ******************************************************************* 

        
        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then 
          open(unit=iMDA,file='mDA.dat', status='unknown')   
   
         endif
        endif 

        mdaRn(1)%xM = 0.d0 
        mdaRn(1)%yM = 0.d0 
        mdaRn(1)%zM = 0.d0 

       ! Find the difference dipole moment
        do irc = 1, iNumDB
         iAtm = dVosDbase(irc)%atom
         qi = dVosDbase(irc)%q

         mdaRn(1)%xM = mdaRn(1)%xM + 
     x               (nodeConf(iAtm)%x)*qi
         mdaRn(1)%yM = mdaRn(1)%yM + 
     x               (nodeConf(iAtm)%y)*qi
         mdaRn(1)%zM = mdaRn(1)%zM + 
     x               (nodeConf(iAtm)%z)*qi      

        enddo         
       
        ddX = mdaRn(1)%xM
        ddY = mdaRn(1)%yM
        ddZ = mdaRn(1)%zM
       
        mdaRn(1)%mDA = sqrt(ddX**2+ddY**2+ddZ**2)

        call MPI_GATHER(mdaRn,1,MDAR_TYPE,
     x                  mdaRs,1,MDAR_TYPE,
     x                  master,world,ierr)     
       
           if(rank .eq. master)  then
           
            do iWrite = 1, nProc
              write(iMDA,*) iConf+iWrite-1, 
     x         mdaRs(iWrite)%xM, mdaRs(iWrite)%yM,
     x         mdaRs(iWrite)%zM, mdaRs(iWrite)%mDA  
            enddo
           endif       
     
        !call MPI_BARRIER(world,ierr)
        
3332             format (3x,i6,5x,4(F6.3,4x))

        return
      end subroutine calcMDA
      
      subroutine calcGr(iXAtom)
c    *******************************************************************
c    ** calculates radial distribution function for X-O, where O is   **
c    ** from solvent water molecules and X is any atom.               **
c    **                                                               **
c    ** principal variables:                                          **
c    ** integer   iXAtom     the index of the probe (X) atom          **
c    **                                                               **
c    *******************************************************************            
        implicit none
 
        ! Args
        integer :: iXAtom
   

        ! Locals        
        integer :: iw, iWrite
        integer :: xInt, yInt, zInt
        integer :: ibin
        integer :: ix, iy, iz
        real  :: atomX, atomY, atomZ
        real  :: oxyX, oxyY, oxyZ
        real  :: dX, dY, dZ, W 
        real  :: rIJ, testLow, rIJ2
        real*8 :: dfrx, dfry, dfrz, delr2
c    ******************************************************************* 

        
        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then 
          open(unit=iXOfile, 
     x        file=trim(strMinFile), status='unknown')                   
         endif
        endif 
         
!        print *, rank, " allocated lows..."                  
!        call MPI_BARRIER(world,ierr)      
        testLow = 88888888.d0        
        
        if (rank .eq. master) lowRXO(:)  = 88888888.d0          

          ! determine the position of the Cu atom
          atomX = nodeConf(iXAtom)%x
          atomY = nodeConf(iXAtom)%y
          atomZ = nodeConf(iXAtom)%z

 ! Loop over waters and build the partial histogram on each node
!$OMP PARALLEL
!$OMP&DEFAULT(shared) 
!$OMP&REDUCTION(+:hist_i)
!$OMP&PRIVATE(iw,oxyX,oxyY,oxyZ,dX,dY,dZ)
!$OMP&PRIVATE(rIJ,ibin,xInt,yInt,zInt)
!$OMP DO SCHEDULE(static)
          do iw = iWatStart, natom, 3
           
           oxyX = nodeConf(iw)%x
           oxyY = nodeConf(iw)%y
           oxyZ = nodeConf(iw)%z

           if(isTruncated) then
            dX = atomX-oxyX
            dY = atomY-oxyY
            dZ = atomZ-oxyZ

            !xInt = ANINT((atomX-oxyX)/bX)
            !yInt = ANINT((atomY-oxyY)/bY)
            !zInt = ANINT((atomZ-oxyZ)/bZ)
!$OMP END CRITICAL                
 
            !dX = (atomX-oxyX) - xInt*bX
            !dY = (atomY-oxyY) - yInt*bY
            !dZ = (atomZ-oxyZ) - zInt*bZ
            
           dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
           dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
           dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

           dfrx = dfrx - ANINT(dfrx)
           dfry = dfry - ANINT(dfry)
           dfrz = dfrz - ANINT(dfrz)

          
          rIJ2 = dX**2+dY**2+dZ**2
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ2) then 
              rIJ2 = delr2
            endif

            enddo
           enddo
          enddo
          
          rIJ = sqrt(rIJ2)
          
           else ! orthorhombic
!$OMP CRITICAL
            !xInt = ANINT((atomX-oxyX)/bX)
            !yInt = ANINT((atomY-oxyY)/bY)
            !zInt = ANINT((atomZ-oxyZ)/bZ)
!$OMP END CRITICAL                
           dX = atomX-oxyX
           dY = atomY-oxyY
           dZ = atomZ-oxyZ
           if (dX .lt. 0.d0) dX = -dX
           if (dY .lt. 0.d0) dY = -dY
           if (dZ .lt. 0.d0) dZ = -dZ
           
           do while (dX .gt. bX)
             dX = dX - bX
           enddo           
           do while (dY .gt. bY)
             dY = dY - bY
           enddo
           do while (dZ .gt. bZ)
             dZ = dZ - bZ
           enddo
           
           W = bX - dX
           if(W .lt. dX) dX = W
           W = bY - dY
           if(W .lt. dY) dY = W
           W = bZ - dZ
           if(W .lt. dZ) dZ = W
           
          !  dX = dX - xInt*bX
          !  dY = dY - yInt*bY
          !  dZ = dZ - zInt*bZ           
!           dX = (cuX-oxyX) - bX * ANINT((cuX-oxyX)/bX)
!           dY = (cuY-oxyY) - bY * ANINT((cuY-oxyY)/bY)
!           dZ = (cuZ-oxyZ) - bZ * ANINT((cuZ-oxyZ)/bZ)
         
           rIJ = sqrt(dX**2+dY**2+dZ**2 ) 
           
!$OMP ATOMIC         
          !  rIJ = sqrt(dX**2+dY**2+dZ**2 )  ! Reimaged X-O distance
            
           endif ! End test for truncated box
            
            if (rIJ .lt. testLow) then
!$OMP ATOMIC                   
             testLow = rIJ             
            endif


            !Bin the distance into the partial histrogram on the node
!$OMP ATOMIC            
           ibin = int ( rIJ/delr ) + 1

           if(ibin .lt. maxBin) then
            hist_i(ibin) = hist_i(ibin) + 1.d0
           endif   

          enddo ! End water loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

         if(rank .eq. master) lowRXO(:) = 0.d0
              
           call MPI_GATHER(testLow,1,MPI_REAL,
     x                      lowRXO,1,MPI_REAL,
     x                       master,world,ierr)     
     
           if(rank .eq. master)  then
           
            do iWrite = 1, nProc
              write(iXOfile, *) iConf+iWrite-1, 
     x         lowRXO(iWrite)
            enddo
           endif

        return
      end subroutine      

      subroutine calcQdR()
c    *******************************************************************
c    ** sums the current possition to the average which will be       **
c    ** processed to give average positions.                          **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    *******************************************************************            
        implicit none
 
        ! Locals        
        integer :: iip
        integer :: iatm    
        real    :: xxx, yyy, zzz
        real    :: ddX, ddY, ddZ
        real    :: qdr2, qdrx, qdry, qdrz
c    ******************************************************************* 

           ! Initialize the avgerages
           if (iConf .eq. 1 ) then
            drxAvg(:) = 0.d0
            dryAvg(:) = 0.d0
            drzAvg(:) = 0.d0
            qdrsq(:)  = 0.d0
           endif 
       
           ! Loop over all protein atoms and collect the avg. struct.       
           do iip = 1, iProtEnd

!**** add code here to remove trans./rot. movement **********
           xxx  = nodeConf(iip)%x
           yyy  = nodeConf(iip)%y
           zzz  = nodeConf(iip)%z
!************************************************************

           ddX = xxx - xAvg(iip)
           ddY = yyy - yAvg(iip)
           ddZ = zzz - zAvg(iip)
           
           drxAvg(iip) = drxAvg(iip) + ddX
           dryAvg(iip) = dryAvg(iip) + ddY
           drzAvg(iip) = drzAvg(iip) + ddZ

           !qdrx = q*ddX
           !qdry = q*ddY
           !qdrz = q*ddZ

           qdr2 = ddX**2+ddY*2+ddZ*2 !qdrx**2 + qdry**2+ qdrz**2
           qdrsq(iip) = qdrsq(iip) + qdr2 ! needs q^2 term
           
           enddo ! loop over atoms
           
        return
      end subroutine calcQdR

      subroutine calcAvgStruct()
c    *******************************************************************
c    ** sums the current possition to the average which will be       **
c    ** processed to give average positions.                          **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    *******************************************************************            
        implicit none
 
        ! Locals        
        logical :: goodMatrix
        logical :: doOutput
        logical :: setVecs
        integer :: nP
        integer :: iip, jjp, kkp
        integer :: ilow, massOrQ
        integer :: iDUMMY, laINFO
        integer, parameter :: iWORK=48
        real    :: xxx, yyy, zzz
        real    :: xi,yi,zi
        real    :: totM, prM
        real    :: rTmp
        real*8  :: r11, r12, r13
        real*8  :: r21, r22, r23
        real*8  :: r31, r32, r33
        real*8  :: one1,one2,one3
        real*8  :: two1,two2,two3
        real*8  :: three1,three2,three3
        real*8  :: kMat_tmp(3,3)
        real*8  :: kMat(3,3), det
        real*8  :: bMat(3,3), aMat(3,3)
        real*8  :: rMat(3,3)
        real*8  :: eig_valR(3), eig_valI(3)
        real*8  :: lowVal, bConst
        real*8  :: laWORK(iWORK), laDUMMY
c    ******************************************************************* 
            nP = (iProtEnd - iProtStart) + 1
           ! Initialize the avgerages
           if (iConf .eq. 1 ) then
            xAvg(:) = 0.d0
            yAvg(:) = 0.d0
            zAvg(:) = 0.d0


            if(doSelect) then
             allocate(proteinCOM(nSelect))
             allocate(proteinCRYS(nSelect))
             allocate(proteinROT(nSelect))
            else
             allocate(proteinCOM(nP))
             allocate(proteinCRYS(nP))
             allocate(proteinROT(nP))
            endif


            if(rank .eq. master) then
             !!! IMPORTANT !!!
             !!! --> should get this structure from a file

            if(doSelect) then
              do iip=1, nSelect
               jjp = selectIndices(iip)
               proteinCRYS(iip) = nodeConf(jjp)
              enddo
            else ! no selection
               jjp = 0
              do iip=iProtStart,iProtEnd
               jjp = jjp + 1
               proteinCRYS(jjp) = nodeConf(iip)
              enddo
            endif ! check for select

            endif ! check for master

            if(doSelect) then
             call MPI_BCAST(proteinCRYS,nSelect,ATOM_TYPE,
     x                     master,world,ierr)  
            else
             call MPI_BCAST(proteinCRYS,nP,ATOM_TYPE,
     x                     master,world,ierr)  
            endif

            call MPI_BARRIER(world,ierr)

           if (useCUDA) then
#ifdef ENABLE_CUDA

           doOutput = .false.
           massOrQ = 0 !set to 0 to calculate the center of mass vector
            
           if(doSelect) then
            ! send the crystal structure to the GPU
            call set_prot_sel_gpu(proteinCRYS)

            ! Calculate the crystal protein's COM vector
            call get_prmu_sel_gpu(qrms,
     x                   0.00000, 0.00000, 0.000000,  
     x              prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x                 totM,massOrQ,selectIndices, 
     x                     doOutput,proteinCOM)

            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM

            massOrQ = 2 !set to 2 to calculate the COM coordinates
            doOutput = .true.
            call get_prmu_sel_gpu(qrms,
     x           prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x              prMuHat%x,prMuHat%y,prMuHat%z,
     x                 totM,massOrQ,selectIndices, 
     x                     doOutput,proteinCOM)


           else

            ! send the crystal structure to the GPU
            call set_prot_gpu(proteinCRYS) 
            ! Calculate the crystal protein's COM vector
            call get_prmu_gpu(qrms,
     x                   0.00000, 0.00000, 0.000000,  
     x              prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x              totM,massOrQ,doOutput,proteinCOM)
            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM

            massOrQ = 2 !set to 2 to calculate the COM coordinates
            doOutput = .true.
           ! now calculate the proteins coordinates in the COM frame
            call get_prmu_gpu(qrms,
     x           prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x        prMuHat%x,prMuHat%y,prMuHat%z,totM,massOrQ,
     x         doOutput,proteinCOM)

            endif ! end check for selection

            prMuCntr%x = prMuHat%x/totM
            prMuCntr%y = prMuHat%y/totM
            prMuCntr%z = prMuHat%z/totM

            setVecs = .false.
           if (doSelect) then
           do iip = 1, nSelect  ! Move all crystal atoms to their COM position
            proteinCRYS(iip) = proteinCOM(iip)
!proteinCRYS(iip) -
!     x                         prMuCntr
           enddo ! End loop over all protein atoms 
            ! Store the crystal structure in COM frame on the GPU

            call store_crys_sel_gpu(proteinCRYS,qrms,
     x                        setVecs,selectIndices) 

           else ! not a selection

           do iip = 1, nP  ! Move all crystal atoms to their COM position
            proteinCRYS(iip) = proteinCOM(iip)
!            proteinCRYS(iip) = proteinCRYS(iip) -
!     x                         prMuCntr
           enddo ! End loop over all protein atoms 
            ! Store the crystal structure in COM frame on the GPU
            call store_crys_gpu(proteinCRYS,qrms,setVecs)

           endif  ! End check for selection
#endif
           else ! a non CUDA run
           prMuCntr%x = 0.d0
           prMuCntr%y = 0.d0
           prMuCntr%z = 0.d0
           totM = 0.d0
           jjp = 0

           if (doSelect) then
            do iip = 1, nSelect  ! Loop over protein atoms to find the COM coordinates
             jjp = selectIndices(iip)
             prAtom = proteinCRYS(iip)
            
             prM = qrms(jjp)%m
             
             totM = totM + prM           

             prMuCntr%x = prMuCntr%x + prM*prAtom%x
             prMuCntr%y = prMuCntr%y + prM*prAtom%y
             prMuCntr%z = prMuCntr%z + prM*prAtom%z
             
            enddo ! End loop over all protein atoms 
            
            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM

            do iip = 1, nSelect  ! Move all crystal atoms to their COM position
             prAtom = proteinCRYS(iip) ! set the position            
             proteinCRYS(iip) = prAtom - prMuCntr
            enddo ! End loop over all protein atoms 


           else ! no selection
            do iip = iProtStart, iProtEnd  ! Loop over protein atoms to find the COM coordinates
             jjp = jjp + 1
             prAtom = proteinCRYS(jjp)
            
             prM = qrms(iip)%m
             
             totM = totM + prM           

             prMuCntr%x = prMuCntr%x + prM*prAtom%x
             prMuCntr%y = prMuCntr%y + prM*prAtom%y
             prMuCntr%z = prMuCntr%z + prM*prAtom%z
             
            enddo ! End loop over all protein atoms 
            
            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM

            do iip = 1, nP  ! Move all crystal atoms to their COM position
             prAtom = proteinCRYS(iip) ! set the position            
             proteinCRYS(iip) = prAtom - prMuCntr
            enddo ! End loop over all protein atoms 

           endif ! end check for selection
           endif ! check for CUDA

           endif !check for 1st conf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! END OF 1st CONF !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Init all matrices
            do iip = 1, 3
             do jjp = 1, 3
              kMat(iip,jjp) = 0.d0
              kMat_tmp(iip,jjp) = 0.d0
              bMat(iip,jjp) = 0.d0
             aMat(iip,jjp) = 0.d0
             rMat(iip,jjp) = 0.d0
             enddo
            enddo

           if (useCUDA) then
#ifdef ENABLE_CUDA

ccccccccccc ROTATION AND TRANSLATION REMOVAL cccccccccccccc 
           doOutput = .false.
           massOrQ = 0 !set to 0 to calculate the center of mass vector

           if(doSelect) then

            ! send the crystal structure to the GPU
            ! but, use proteinCOM as a vessel since
            ! set_prot_sel_gpu() cant take nodeConf
            ! as input
            do iip = 1, nSelect
             jjp = selectIndices(iip)             
             proteinCOM(iip) = nodeConf(jjp)
            enddo          
            call set_prot_sel_gpu(proteinCOM)

            ! Calculate the crystal protein's COM vector
            call get_prmu_sel_gpu(qrms,
     x                   0.00000, 0.00000, 0.000000,  
     x              prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x                 totM,massOrQ,selectIndices, 
     x                     doOutput,proteinCOM)
            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM
            massOrQ = 2 !set to 2 to calculate the COM coordinates
            doOutput = .true.
            call get_prmu_sel_gpu(qrms,
     x           prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x              prMuHat%x,prMuHat%y,prMuHat%z,
     x                 totM,massOrQ,selectIndices, 
     x                     doOutput,proteinCOM)

           else  ! no selection

            ! send the crystal structure to the GPU
            call set_prot_gpu(nodeConf)
           ! Calculate the configuration protein's COM vector
            call get_prmu_gpu(qrms,
     x                 0.000000, 0.000000, 0.000000,
     x              prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x              totM,massOrQ,doOutput,proteinCOM)
            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM
            massOrQ = 2 !set to 2 to calculate the COM coordinates
            doOutput = .true.
           ! now calculate the proteins coordinates in the COM frame
            call get_prmu_gpu(qrms,
     x           prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x        prMuHat%x,prMuHat%y,prMuHat%z,totM,massOrQ,
     x         doOutput,proteinCOM)
           endif ! end check for selection

           if(doSelect) then
            ! now calculate the proteins coordinates in the COM frame
            call calc_kmat_sel_gpu(proteinCOM, 
     x                       one1, one2, one3, 
     x                       two1, two2, two3, 
     x                 three1, three2, three3)
           else  ! no selection

            ! now calculate the proteins coordinates in the COM frame  
            call calc_kmat_gpu(proteinCOM,
     x                        one1, one2, one3, 
     x                        two1, two2, two3, 
     x                        three1, three2, three3)

           endif ! end of check for the selection 

           kMat_tmp(1,1) = one1
           kMat_tmp(1,2) = one2
           kMat_tmp(1,3) = one3

           kMat_tmp(2,1) = two1
           kMat_tmp(2,2) = two2
           kMat_tmp(2,3) = two3

           kMat_tmp(3,1) = three1
           kMat_tmp(3,2) = three2
           kMat_tmp(3,3) = three3
ccccccccccc END ROTATION AND TRANSLATION REMOVAL cccccccccccccc 
#endif
           else  ! a non-cuda run

           prMuCntr%x = 0.d0
           prMuCntr%y = 0.d0
           prMuCntr%z = 0.d0
           totM = 0.d0

           if(doSelect) then
            do iip = 1, nSelect  ! Loop over sel atoms to find the COM coordinates

             jjp = selectIndices(iip)

             prAtom = nodeConf(jjp)            
             prM = qrms(jjp)%m
            
             totM = totM + prM           

             prMuCntr%x = prMuCntr%x + prM*prAtom%x
             prMuCntr%y = prMuCntr%y + prM*prAtom%y
             prMuCntr%z = prMuCntr%z + prM*prAtom%z
            
            enddo ! End loop over all sel atoms 


           else ! no selection

            do iip = iProtStart, iProtEnd  ! Loop over protein atoms to find the COM coordinates
             prAtom = nodeConf(iip)            
             prM = qrms(iip)%m
            
             totM = totM + prM           

             prMuCntr%x = prMuCntr%x + prM*prAtom%x
             prMuCntr%y = prMuCntr%y + prM*prAtom%y
             prMuCntr%z = prMuCntr%z + prM*prAtom%z
            
            enddo ! End loop over all protein atoms 

           endif !end selection 
            
            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM

           if(doSelect) then
            do iip = 1, nSelect
             jjp = selectIndices(iip)
             prAtom = nodeConf(jjp) ! set the position            
             proteinCOM(iip) = prAtom - prMuCntr !!! move to the COM
            enddo ! End loop over all select atoms 

            do iip = 1, nSelect
             jjp = selectIndices(iip)

             prAtom = proteinCOM(iip)
             prCrys = proteinCRYS(iip)
             prM = qrms(jjp)%m

             kMat_tmp(1,1) = kMat_tmp(1,1) + 
     x       real(prM * prCrys%x * prAtom%x,kind=8)
             kMat_tmp(1,2) = kMat_tmp(1,2) + 
     x       real(prM * prCrys%x * prAtom%y,kind=8)
             kMat_tmp(1,3) = kMat_tmp(1,3) + 
     x       real(prM * prCrys%x * prAtom%z,kind=8)
             kMat_tmp(2,1) = kMat_tmp(2,1) + 
     x       real(prM * prCrys%y * prAtom%x,kind=8)
             kMat_tmp(2,2) = kMat_tmp(2,2) + 
     x       real(prM * prCrys%y * prAtom%y,kind=8)
             kMat_tmp(2,3) = kMat_tmp(2,3) + 
     x       real(prM * prCrys%y * prAtom%z,kind=8)
             kMat_tmp(3,1) = kMat_tmp(3,1) + 
     x       real(prM * prCrys%z * prAtom%x,kind=8)
             kMat_tmp(3,2) = kMat_tmp(3,2) + 
     x       real(prM * prCrys%z * prAtom%y,kind=8)
             kMat_tmp(3,3) = kMat_tmp(3,3) + 
     x       real(prM * prCrys%z * prAtom%z,kind=8)
            end do
           else ! no selection

            jjp = 0
            do iip = iProtStart, iProtEnd
             jjp = jjp + 1
             prAtom = nodeConf(iip) ! set the position            
             proteinCOM(jjp) = prAtom - prMuCntr  ! Move to the COM
            enddo ! End loop over all protein atoms 

            jjp = 0
            do iip = iProtStart, iProtEnd
             jjp = jjp + 1

             prAtom = proteinCOM(jjp)
             prCrys = proteinCRYS(jjp)
             prM = qrms(iip)%m

             kMat_tmp(1,1) = kMat_tmp(1,1) + 
     x       real(prM * prCrys%x * prAtom%x,kind=8)
             kMat_tmp(1,2) = kMat_tmp(1,2) + 
     x       real(prM * prCrys%x * prAtom%y,kind=8)
             kMat_tmp(1,3) = kMat_tmp(1,3) + 
     x       real(prM * prCrys%x * prAtom%z,kind=8)
             kMat_tmp(2,1) = kMat_tmp(2,1) + 
     x       real(prM * prCrys%y * prAtom%x,kind=8)
             kMat_tmp(2,2) = kMat_tmp(2,2) + 
     x       real(prM * prCrys%y * prAtom%y,kind=8)
             kMat_tmp(2,3) = kMat_tmp(2,3) + 
     x       real(prM * prCrys%y * prAtom%z,kind=8)
             kMat_tmp(3,1) = kMat_tmp(3,1) + 
     x       real(prM * prCrys%z * prAtom%x,kind=8)
             kMat_tmp(3,2) = kMat_tmp(3,2) + 
     x       real(prM * prCrys%z * prAtom%y,kind=8)
             kMat_tmp(3,3) = kMat_tmp(3,3) + 
     x       real(prM * prCrys%z * prAtom%z,kind=8)
            end do

         endif ! end check for selection

       endif ! End check for CUDA

      ! create a positive definitive kMat by multiplying kMat_tmp by
      ! by its transpose
           do iip=1,3
             do jjp=1,3
              do kkp=1,3
               kMat(iip,jjp) = kMat(iip,jjp) + 
     x         kMat_tmp(kkp,iip)*kMat_tmp(kkp,jjp)
              enddo
             enddo
           enddo

      ! diagonalize kMat and store the eigenvectors in
      ! aMat(i=1,3 ; j) for each j-th eigenvector.
      ! 'eig_valR' has the eigenvalues.   
         iDUMMY = 1
         call dgeev('N','V',3,kMat,
     x         3,eig_valR,eig_valI,
     x              laDUMMY,iDUMMY,
     x        aMat,3,laWORK,iWORK,laINFO)

      ! make sure the determinant is non-zero
       det = kMat(1,1)*kMat(2,2)*kMat(3,3) - 
     x       kMat(1,1)*kMat(2,3)*kMat(3,2) -  
     x       kMat(1,2)*kMat(2,1)*kMat(3,3) + 
     x       kMat(1,2)*kMat(2,3)*kMat(3,1) +  
     x       kMat(1,3)*kMat(2,1)*kMat(3,2) - 
     x       kMat(1,3)*kMat(2,2)*kMat(3,1)
   
        if (abs(det) < 1.d-5) then

         write (*,*) "POSSIBLE ERROR IN KABSCH MATRIX!: ", 
     x       abs(det)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        lowVal = 1.d20
        ! find lowest magnitude eigenvalue
        do iip=1,3
         if (eig_valR(iip) < lowVal) then
            ilow = iip
            lowVal  = eig_valR(iip)
         end if
        end do
   
       ! create the B vectors fromt the eigenvectors
       do jjp=1,3
        bConst = 
     x  1.d0/sqrt(abs(eig_valR(jjp)))

        do iip=1,3
         do kkp=1,3
           bMat(iip,jjp) = bMat(iip,jjp) + 
     x               bConst*aMat(kkp,jjp)*
     x                  kMat_tmp(iip,kkp)

         enddo
        enddo
       enddo
   
      ! now iteratively calculate the rotation matrix
      goodMatrix = .false.

       do while (.not. goodMatrix)
        rMat(:,:) = 0.d0  ! array assignment

         do iip=1,3
          do jjp=1,3
           do kkp=1,3
           rMat(iip,jjp) = rMat(iip,jjp) + 
     x             bMat(iip,kkp)*aMat(jjp,kkp)
           enddo
          enddo
         enddo

        ! check if the determinant is negative
        det =  rMat(1,1)*rMat(2,2)*rMat(3,3) -  
     x        rMat(1,1)*rMat(2,3)*rMat(3,2) -  
     x        rMat(1,2)*rMat(2,1)*rMat(3,3) +  
     x        rMat(1,2)*rMat(2,3)*rMat(3,1) +  
     x        rMat(1,3)*rMat(2,1)*rMat(3,2) -  
     x        rMat(1,3)*rMat(2,2)*rMat(3,1)
   
        if (abs(det) < 1.d-10) then
          write (*,*) "SMALL rMat DETERMINANT!: ", 
     x       abs(det)
        end if
   
        ! if the determinant is negative then invert to
        ! obtain a new rotation matrix
        if (det < 0) then
          do iip=1,3
             bMat(iip,ilow) = -bMat(iip,ilow)
          enddo
          goodMatrix = .false.
        else 
          goodMatrix = .true.
        endif

      end do ! test for goodMatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           if (useCUDA) then
#ifdef ENABLE_CUDA
           if(doSelect) then
            call accum_avg_sel_gpu(proteinCOM,
     x             rMat(1,1), rMat(1,2), rMat(1,3), 
     x             rMat(2,1), rMat(2,2), rMat(2,3), 
     x             rMat(3,1), rMat(3,2), rMat(3,3))

           else ! no selection
            call accum_avg_gpu(proteinCOM,
     x             rMat(1,1), rMat(1,2), rMat(1,3), 
     x             rMat(2,1), rMat(2,2), rMat(2,3), 
     x             rMat(3,1), rMat(3,2), rMat(3,3))

           endif ! check for select
#endif
           else  ! a non-cuda run            

            r11 = rMat(1,1)
            r12 = rMat(1,2)
            r13 = rMat(1,3)

            r21 = rMat(2,1)
            r22 = rMat(2,2)
            r23 = rMat(2,3)

            r31 = rMat(3,1)
            r32 = rMat(3,2)
            r33 = rMat(3,3)


           if(doSelect) then
            ! Loop over all protein atoms and collect the avg. struct.       
            do iip = 1, nSelect

            xxx  = proteinCOM(iip)%x
            yyy  = proteinCOM(iip)%y
            zzz  = proteinCOM(iip)%z

            ! rotate the protein atom thats beens trans. to its COM
            xi = real(r11*xxx + 
     x            r12*yyy + 
     x            r13*zzz,kind=4)
            yi = real(r21*xxx + 
     x            r22*yyy + 
     x            r23*zzz,kind=4)
            zi = real(r31*xxx + 
     x            r32*yyy + 
     x            r33*zzz,kind=4)

            ! then accumulate the average
            xAvg(iip) = xAvg(iip) + xxx
            yAvg(iip) = yAvg(iip) + yyy
            zAvg(iip) = zAvg(iip) + zzz
           
            enddo ! loop over select atoms

           else ! no selection
            ! Loop over all protein atoms and collect the avg. struct.       
            do iip = 1, nP

            xxx  = proteinCOM(iip)%x
            yyy  = proteinCOM(iip)%y
            zzz  = proteinCOM(iip)%z

            ! rotate the protein atom thats been trans. to its COM
            xi = real(r11*xxx + 
     x            r12*yyy + 
     x            r13*zzz,kind=4)
            yi = real(r21*xxx + 
     x            r22*yyy + 
     x            r23*zzz,kind=4)
            zi = real(r31*xxx + 
     x            r32*yyy + 
     x            r33*zzz,kind=4)

            xAvg(iip) = xAvg(iip) + xxx
            yAvg(iip) = yAvg(iip) + yyy
            zAvg(iip) = zAvg(iip) + zzz
           
            enddo ! loop over atoms
          endif ! check for selection
         endif ! End check for CUDA
          
         return
      end subroutine calcAvgStruct

      subroutine calcCovar()
c    *******************************************************************
c    ** uses the average file to calculate the wass-weighted          **
c    ** covariance matrix for a selection of (protein) atoms.         **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    *******************************************************************            
        implicit none
 
        ! Locals        
        logical :: goodMatrix
        logical :: doOutput
        logical :: setVecs
        integer :: nP
        integer :: iip, jjp, kkp
        integer :: ilow, massOrQ
        integer :: iDUMMY, laINFO
        integer, parameter :: iWORK=48
        real    :: xxx, yyy, zzz
        real    :: xi,yi,zi
        real    :: totM, prM
        real    :: rTmp, sqrtM
        real*8  :: r11, r12, r13
        real*8  :: r21, r22, r23
        real*8  :: r31, r32, r33
        real*8  :: one1,one2,one3
        real*8  :: two1,two2,two3
        real*8  :: three1,three2,three3
        real*8  :: kMat_tmp(3,3)
        real*8  :: kMat(3,3), det
        real*8  :: bMat(3,3), aMat(3,3)
        real*8  :: rMat(3,3)
        real*8  :: eig_valR(3), eig_valI(3)
        real*8  :: lowVal, bConst
        real*8  :: laWORK(iWORK), laDUMMY
c    ******************************************************************* 
            nP = (iProtEnd - iProtStart) + 1
!            print *, "ABOUT TO INIT MATRICES..."
            ! Init all matrices
            do iip = 1, 3
             do jjp = 1, 3
              kMat(iip,jjp) = 0.d0
              kMat_tmp(iip,jjp) = 0.d0
              bMat(iip,jjp) = 0.d0
             aMat(iip,jjp) = 0.d0
             rMat(iip,jjp) = 0.d0
             enddo
            enddo
!            print *, "ABOUT TO DO 1st CONFIG STUFF..."
!!!!!!!!!!!!!!!!!!!!!!!! FIRST CONFIGURATION ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! FIRST CONFIGURATION ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! FIRST CONFIGURATION ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!
           if (iConf .eq. 1 ) then
            xAvg(:) = 0.d0
            yAvg(:) = 0.d0
            zAvg(:) = 0.d0


            if(doSelect) then
             iPCADoF = 3*nSelect
             allocate(proteinCOM(nSelect))
             allocate(proteinCRYS(nSelect))
!             allocate(proteinROT(nSelect))
            else
             iPCADoF = 3*nP
             allocate(proteinCOM(nP))
             allocate(proteinCRYS(nP))
!             allocate(proteinROT(nP))
            endif

            allocate(covarNode(iPCADoF,iPCADoF))
            allocate(covarMat(iPCADoF,iPCADoF))
            allocate(vecAvg(iPCADoF))
            allocate(vecConf(iPCADoF))

            covarMat(:,:) = 0.d0
            covarNode(:,:) = 0.d0

           if(rank .eq. master) then

           open(unit=1233, file='avgProt.dat', 
     x        status='old')
             ! Read the avgerage file

            if(doSelect) then
              do iip = 1, nSelect
               read(1233,*) jjp, 
     x         proteinCRYS(iip)%x, 
     x         proteinCRYS(iip)%y, 
     x         proteinCRYS(iip)%z
              enddo            

            else ! no selection
              kkp = 0
              do iip = iProtStart, iProtEnd
               kkp = kkp + 1        
               read(1233,*) jjp, 
     x         proteinCRYS(kkp)%x, 
     x         proteinCRYS(kkp)%y, 
     x         proteinCRYS(kkp)%z
              enddo            
            endif ! end check for selection

            close(1233)
           end if ! check for master
           call MPI_BARRIER(world,ierr)

          if(doSelect) then
            call MPI_BCAST(proteinCRYS,nSelect,
     x                              ATOM_TYPE,
     x                     master,world,ierr)  
          else ! no select
            call MPI_BCAST(proteinCRYS,nP,ATOM_TYPE,
     x                     master,world,ierr)  
          endif ! check for select

           if (useCUDA) then
#ifdef ENABLE_CUDA
             setVecs = .true.
            if(doSelect) then
             call store_crys_sel_gpu(proteinCRYS,qrms,
     x                        setVecs,selectIndices) 
             call init_covar_sel_gpu()              
            else ! no selection
             call store_crys_gpu(proteinCRYS,qrms,setVecs)
             call init_covar_gpu()
            endif ! no selection

#endif
           else ! a non CUDA run

            if(doSelect) then

             ! put the average coordinates into an easily loop-able vector
              jjp = 0
              do iip = 1, nSelect
               kkp = selectIndices(iip)
               sqrtM = sqrt(qrms(kkp)%m)

               jjp = jjp + 1
               vecAvg(jjp) = proteinCRYS(iip)%x*
     x                                    sqrtM 
               jjp = jjp + 1
               vecAvg(jjp) = proteinCRYS(iip)%y*
     x                                    sqrtM 
               jjp = jjp + 1
               vecAvg(jjp) = proteinCRYS(iip)%z*
     x                                    sqrtM 
!              if (iip .eq. iProtEnd) then
!                 print *, "***** AVG IN FORTRAN: *****"
!                 print *, "      sqrt(m)     = ", 
!     x                               sqrtM
!                 print *, "      crys(x,y,z) = ", 
!     x                   proteinCRYS(kkp)%x,
!     x                   proteinCRYS(kkp)%y,
!     x                   proteinCRYS(kkp)%z
!                 print *, "     vec avg(k,k-1,k-2) = ", 
!     x                   vecAvg(jjp),
!     x                   vecAvg(jjp-1),
!     x                   vecAvg(jjp-2)
!                 print *, "***** AVG IN FORTRAN: *****"
!              endif
             enddo ! end loop over putting avgs into a vector

            else ! no selection

             ! put the average coordinates into an easily loop-able vector
              jjp = 0
              kkp = 0
              do iip = iProtStart, iProtEnd
               kkp = kkp + 1
               sqrtM = sqrt(qrms(iip)%m)
               jjp = jjp + 1
               vecAvg(jjp) = proteinCRYS(kkp)%x*
     x                                    sqrtM 
               jjp = jjp + 1
               vecAvg(jjp) = proteinCRYS(kkp)%y*
     x                                    sqrtM 
               jjp = jjp + 1
               vecAvg(jjp) = proteinCRYS(kkp)%z*
     x                                    sqrtM 
!              if (iip .eq. iProtEnd) then
!                 print *, "***** AVG IN FORTRAN: *****"
!                 print *, "      sqrt(m)     = ", 
!     x                               sqrtM
!                 print *, "      crys(x,y,z) = ", 
!     x                   proteinCRYS(kkp)%x,
!     x                   proteinCRYS(kkp)%y,
!     x                   proteinCRYS(kkp)%z
!                 print *, "     vec avg(k,k-1,k-2) = ", 
!     x                   vecAvg(jjp),
!     x                   vecAvg(jjp-1),
!     x                   vecAvg(jjp-2)
!                 print *, "***** AVG IN FORTRAN: *****"
!              endif
             enddo ! end loop over putting avgs into a vector
            endif ! check for selection

           endif ! check for CUDA

           endif !check for 1st conf

!!!!!!!!!!!!!!!!!!!!!!!! END FIRST CONFIGURATION ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! END FIRST CONFIGURATION ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! END FIRST CONFIGURATION ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!


           if (useCUDA) then
#ifdef ENABLE_CUDA
ccccccccccc ROTATION AND TRANSLATION REMOVAL cccccccccccccc 
           doOutput = .false.
           massOrQ = 0 !set to 0 to calculate the center of mass vector

           if(doSelect) then

            ! send the crystal structure to the GPU
            ! but, use proteinCOM as a vessel since
            ! set_prot_sel_gpu() cant take nodeConf
            ! as input
            do iip = 1, nSelect
             jjp = selectIndices(iip)             
             proteinCOM(iip) = nodeConf(jjp)
            enddo          
            call set_prot_sel_gpu(proteinCOM)

            ! Calculate the crystal protein's COM vector
            call get_prmu_sel_gpu(qrms,
     x                   0.00000, 0.00000, 0.000000,  
     x              prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x                 totM,massOrQ,selectIndices, 
     x                     doOutput,proteinCOM)
            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM
            massOrQ = 2 !set to 2 to calculate the COM coordinates
            doOutput = .true.
            call get_prmu_sel_gpu(qrms,
     x           prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x              prMuHat%x,prMuHat%y,prMuHat%z,
     x                 totM,massOrQ,selectIndices, 
     x                     doOutput,proteinCOM)

           else  ! no selection

            ! send the crystal structure to the GPU
            call set_prot_gpu(nodeConf)
           ! Calculate the configuration protein's COM vector
            call get_prmu_gpu(qrms,
     x                 0.000000, 0.000000, 0.000000,
     x              prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x              totM,massOrQ,doOutput,proteinCOM)
            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM
            massOrQ = 2 !set to 2 to calculate the COM coordinates
            doOutput = .true.
           ! now calculate the proteins coordinates in the COM frame
            call get_prmu_gpu(qrms,
     x           prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x        prMuHat%x,prMuHat%y,prMuHat%z,totM,massOrQ,
     x         doOutput,proteinCOM)
           endif ! end check for selection

           if(doSelect) then
            ! now calculate the proteins coordinates in the COM frame
!            do iip = 1, nSelect  ! Move all protein atoms to their COM position
!             proteinROT(iip) = proteinCOM(iip)
!            enddo ! End loop over all protein atoms 

!            call calc_kmat_sel_gpu(proteinROT, 
            call calc_kmat_sel_gpu(proteinCOM, 
     x                       one1, one2, one3, 
     x                       two1, two2, two3, 
     x                 three1, three2, three3)
           else  ! no selection

            ! now calculate the proteins coordinates in the COM frame  

!           do iip = 1, nP ! Move all protein atoms to their COM position
!            proteinROT(iip) = proteinCOM(iip)
!           enddo ! End loop over all protein atoms 


!            call calc_kmat_gpu(proteinROT,
            call calc_kmat_gpu(proteinCOM,
     x                        one1, one2, one3, 
     x                        two1, two2, two3, 
     x                        three1, three2, three3)

           endif ! end of check for the selection 

           kMat_tmp(1,1) = one1
           kMat_tmp(1,2) = one2
           kMat_tmp(1,3) = one3

           kMat_tmp(2,1) = two1
           kMat_tmp(2,2) = two2
           kMat_tmp(2,3) = two3

           kMat_tmp(3,1) = three1
           kMat_tmp(3,2) = three2
           kMat_tmp(3,3) = three3

ccccccccccc END ROTATION AND TRANSLATION REMOVAL cccccccccccccc 

#endif
           else  ! a non-cuda run

           prMuCntr%x = 0.d0
           prMuCntr%y = 0.d0
           prMuCntr%z = 0.d0
           totM = 0.d0

           if(doSelect) then
            do iip = 1, nSelect  ! Loop over sel atoms to find the COM coordinates

             jjp = selectIndices(iip)

             prAtom = nodeConf(jjp)            
             prM = qrms(jjp)%m
            
             totM = totM + prM           

             prMuCntr%x = prMuCntr%x + prM*prAtom%x
             prMuCntr%y = prMuCntr%y + prM*prAtom%y
             prMuCntr%z = prMuCntr%z + prM*prAtom%z
            
            enddo ! End loop over all sel atoms 


           else ! no selection

            do iip = iProtStart, iProtEnd  ! Loop over protein atoms to find the COM coordinates
             prAtom = nodeConf(iip)            
             prM = qrms(iip)%m
            
             totM = totM + prM           

             prMuCntr%x = prMuCntr%x + prM*prAtom%x
             prMuCntr%y = prMuCntr%y + prM*prAtom%y
             prMuCntr%z = prMuCntr%z + prM*prAtom%z
            
            enddo ! End loop over all protein atoms 

           endif !end selection 
            
            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM
           

           if(doSelect) then
            do iip = 1, nSelect
             jjp = selectIndices(iip)

             prAtom = nodeConf(jjp) ! set the position            
             proteinCOM(iip) = prAtom - prMuCntr
            enddo ! End loop over all select atoms 

            do iip = 1, nSelect
             jjp = selectIndices(iip)

             prAtom = proteinCOM(iip)
             prCrys = proteinCRYS(iip)
             prM = qrms(jjp)%m

             kMat_tmp(1,1) = kMat_tmp(1,1) + 
     x       real(prM * prCrys%x * prAtom%x,kind=8)
             kMat_tmp(1,2) = kMat_tmp(1,2) + 
     x       real(prM * prCrys%x * prAtom%y,kind=8)
             kMat_tmp(1,3) = kMat_tmp(1,3) + 
     x       real(prM * prCrys%x * prAtom%z,kind=8)
             kMat_tmp(2,1) = kMat_tmp(2,1) + 
     x       real(prM * prCrys%y * prAtom%x,kind=8)
             kMat_tmp(2,2) = kMat_tmp(2,2) + 
     x       real(prM * prCrys%y * prAtom%y,kind=8)
             kMat_tmp(2,3) = kMat_tmp(2,3) + 
     x       real(prM * prCrys%y * prAtom%z,kind=8)
             kMat_tmp(3,1) = kMat_tmp(3,1) + 
     x       real(prM * prCrys%z * prAtom%x,kind=8)
             kMat_tmp(3,2) = kMat_tmp(3,2) + 
     x       real(prM * prCrys%z * prAtom%y,kind=8)
             kMat_tmp(3,3) = kMat_tmp(3,3) + 
     x       real(prM * prCrys%z * prAtom%z,kind=8)
            end do
           else ! no selection

            jjp = 0
            do iip = iProtStart, iProtEnd
             jjp = jjp + 1
             prAtom = nodeConf(iip) ! set the position            
             proteinCOM(jjp) = prAtom - prMuCntr
            enddo ! End loop over all protein atoms 

            jjp = 0
            do iip = iProtStart, iProtEnd
             jjp = jjp + 1

             prAtom = proteinCOM(jjp)
             prCrys = proteinCRYS(jjp)
             prM = qrms(iip)%m

             kMat_tmp(1,1) = kMat_tmp(1,1) + 
     x       real(prM * prCrys%x * prAtom%x,kind=8)
             kMat_tmp(1,2) = kMat_tmp(1,2) + 
     x       real(prM * prCrys%x * prAtom%y,kind=8)
             kMat_tmp(1,3) = kMat_tmp(1,3) + 
     x       real(prM * prCrys%x * prAtom%z,kind=8)
             kMat_tmp(2,1) = kMat_tmp(2,1) + 
     x       real(prM * prCrys%y * prAtom%x,kind=8)
             kMat_tmp(2,2) = kMat_tmp(2,2) + 
     x       real(prM * prCrys%y * prAtom%y,kind=8)
             kMat_tmp(2,3) = kMat_tmp(2,3) + 
     x       real(prM * prCrys%y * prAtom%z,kind=8)
             kMat_tmp(3,1) = kMat_tmp(3,1) + 
     x       real(prM * prCrys%z * prAtom%x,kind=8)
             kMat_tmp(3,2) = kMat_tmp(3,2) + 
     x       real(prM * prCrys%z * prAtom%y,kind=8)
             kMat_tmp(3,3) = kMat_tmp(3,3) + 
     x       real(prM * prCrys%z * prAtom%z,kind=8)
            end do

         endif ! end check for selection

       endif ! End check for CUDA
!       print *, "CREATING POS-DEF Kmat ON THE CPU..."
      ! create a positive definitive kMat by multiplying kMat_tmp by
      ! by its transpose
           do iip=1,3
             do jjp=1,3
              do kkp=1,3
               kMat(iip,jjp) = kMat(iip,jjp) + 
     x         kMat_tmp(kkp,iip)*kMat_tmp(kkp,jjp)
              enddo
             enddo
           enddo

      ! diagonalize kMat and store the eigenvectors in
      ! aMat(i=1,3 ; j) for each j-th eigenvector.
      ! 'eig_valR' has the eigenvalues.   
         iDUMMY = 1
         call dgeev('N','V',3,kMat,
     x         3,eig_valR,eig_valI,
     x              laDUMMY,iDUMMY,
     x        aMat,3,laWORK,iWORK,laINFO)

      ! make sure the determinant is non-zero
       det = kMat(1,1)*kMat(2,2)*kMat(3,3) - 
     x       kMat(1,1)*kMat(2,3)*kMat(3,2) -  
     x       kMat(1,2)*kMat(2,1)*kMat(3,3) + 
     x       kMat(1,2)*kMat(2,3)*kMat(3,1) +  
     x       kMat(1,3)*kMat(2,1)*kMat(3,2) - 
     x       kMat(1,3)*kMat(2,2)*kMat(3,1)
   
        if (abs(det) < 1.d-5) then

         write (*,*) "POSSIBLE ERROR IN KABSCH MATRIX!: ", 
     x       abs(det)
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        lowVal = 9999999999999999.0
        ! find lowest magnitude eigenvalue
        do iip=1,3
         if (eig_valR(iip) < lowVal) then
            ilow = iip
            lowVal  = eig_valR(iip)
         end if
        end do
   
       ! create the B vectors fromt the eigenvectors
       do jjp=1,3
        bConst = 
     x  1.d0/sqrt(abs(eig_valR(jjp)))

        do iip=1,3
         do kkp=1,3
           bMat(iip,jjp) = bMat(iip,jjp) + 
     x     bConst*aMat(kkp,jjp)*
     x     kMat_tmp(iip,kkp)

         enddo
        enddo
       enddo
   
      ! now iteratively calculate the rotation matrix
      goodMatrix = .false.

       do while (.not. goodMatrix)
        rMat(:,:) = 0.d0  ! array assignment

         do iip=1,3
          do jjp=1,3
           do kkp=1,3
           rMat(iip,jjp) = rMat(iip,jjp) + 
     x             bMat(iip,kkp)*aMat(jjp,kkp)
           enddo
          enddo
         enddo

        ! check if the determinant is negative
       det = rMat(1,1)*rMat(2,2)*rMat(3,3) -  
     x       rMat(1,1)*rMat(2,3)*rMat(3,2) -  
     x       rMat(1,2)*rMat(2,1)*rMat(3,3) +  
     x       rMat(1,2)*rMat(2,3)*rMat(3,1) +  
     x       rMat(1,3)*rMat(2,1)*rMat(3,2) -  
     x       rMat(1,3)*rMat(2,2)*rMat(3,1)
   
        if (abs(det) < 1.d-10) then
          write (*,*) "SMALL rMat DETERMINANT!: ", 
     x       abs(det)
        end if
   
        ! if the determinant is negative then invert to
        ! obtain a new rotation matrix
        if (det < 0 ) then
          do iip=1,3
             bMat(iip,ilow) = -bMat(iip,ilow)
          enddo
          goodMatrix = .false.
        else 
          goodMatrix = .true.
        endif

      end do ! test for goodMatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           if (useCUDA) then
#ifdef ENABLE_CUDA
           if(doSelect) then
            call accum_covar_sel_gpu(proteinCOM,qrms,
     x                                  selectIndices,
     x            rMat(1,1), rMat(1,2), rMat(1,3), 
     x            rMat(2,1), rMat(2,2), rMat(2,3), 
     x            rMat(3,1), rMat(3,2), rMat(3,3))
           else ! no selection
            call accum_covar_gpu(proteinCOM,qrms,
     x            rMat(1,1), rMat(1,2), rMat(1,3), 
     x            rMat(2,1), rMat(2,2), rMat(2,3), 
     x            rMat(3,1), rMat(3,2), rMat(3,3))
           endif

#endif
           else  ! a non-cuda run            

            ! Loop over all protein atoms and collect the avg. struct.
            r11 = rMat(1,1)
            r12 = rMat(1,2)
            r13 = rMat(1,3)

            r21 = rMat(2,1)
            r22 = rMat(2,2)
            r23 = rMat(2,3)

            r31 = rMat(3,1)
            r32 = rMat(3,2)
            r33 = rMat(3,3)

            if(doSelect) then
             jjp = 0
             do iip = 1, nSelect

              kkp = selectIndices(iip)
              sqrtM = sqrt(qrms(kkp)%m)

              xi  = proteinCOM(iip)%x
              yi  = proteinCOM(iip)%y
              zi  = proteinCOM(iip)%z
            
              ! rotate the protein atom
              xxx = real(r11*xi + 
     x             r12*yi + 
     x             r13*zi,kind=4)
              yyy = real(r21*xi + 
     x            r22*yi + 
     x            r23*zi,kind=4)
              zzz = real(r31*xi + 
     x            r32*yi + 
     x            r33*zi,kind=4)

              ! put the protein coordinates into an easily loop-able vector
              jjp = jjp + 1
              vecConf(jjp) = xxx*
     x                       sqrtM 
              jjp = jjp + 1
              vecConf(jjp) = yyy*
     x                       sqrtM 
              jjp = jjp + 1
              vecConf(jjp) = zzz*
     x                       sqrtM 
             enddo ! loop over selected atoms

            else ! no selection
             jjp = 0       
             kkp = 0
             do iip = iProtStart, iProtEnd
              sqrtM = sqrt(qrms(iip)%m)

              kkp = kkp + 1
              xi  = proteinCOM(kkp)%x
              yi  = proteinCOM(kkp)%y
              zi  = proteinCOM(kkp)%z

              ! rotate the protein atom in COM coords
              xxx = real(r11*xi + 
     x             r12*yi + 
     x             r13*zi,kind=4)
              yyy = real(r21*xi + 
     x            r22*yi + 
     x            r23*zi,kind=4)
              zzz = real(r31*xi + 
     x            r32*yi + 
     x            r33*zi,kind=4)

              ! put the protein coordinates into an easily loop-able vector
              jjp = jjp + 1
              vecConf(jjp) = xxx*
     x                       sqrtM 
              jjp = jjp + 1
              vecConf(jjp) = yyy*
     x                       sqrtM 
              jjp = jjp + 1
              vecConf(jjp) = zzz*
     x                       sqrtM 
             enddo ! loop over atoms
         endif ! end check for selction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do iip = 1, iPCADoF
           do jjp = 1, iPCADoF

            if ((iip .eq. iPCADoF) .and.
     x          (jjp .eq. iPCADoF)) then

         xxx = vecConf(iip) - vecAvg(iip)
         yyy = vecConf(jjp) - vecAvg(jjp)

!         print *, "CPUCPUCPUCPUCPUCPUCPUCPUCPUCPUCPUCPU"
!         print *, "CPUCPU: COVARIANCE MATRIX :CPUCPU"
!         print *, "  C(dof,dof),old = ",
!     x       covarNode(iPCADoF,iPCADoF)
!         print *, "         new val = ",
!     x       xxx*yyy
!         print *, "         vec(i) = ", vecConf(iip),
!     x           "         avg(i) = ", vecAvg(iip)
!         print *, "         vec(j) = ", vecConf(jjp),
!     x            "         avg(j) = ", vecAvg(jjp)

            endif

            covarNode(iip,jjp) = covarNode(iip,jjp) + 
     x      ( ( vecConf(iip) - vecAvg(iip) ) *
     x      ( vecConf(jjp) - vecAvg(jjp) ) )

!            if (iip .ne. jjp) covarNode(jjp,iip) = 
!     x                   covarNode(iip,jjp)
           enddo  ! loop over jjp/iPCADoF
          enddo  ! loop over iip/iPCADoF

         endif ! End check for CUDA
          
         return
      end subroutine calcCovar

      subroutine calcDIon()
c    *******************************************************************
c    ** calculates the distribution of ions and the average of the    **
c    ** ion-ion distances.                                            **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    *******************************************************************            
        implicit none
 
        ! Locals        
        integer :: iion, jion
        integer :: ibin        
        real    :: distIJ
c    ******************************************************************* 
       
           ! Loop over all solute ion pairs and collect the partial distributions          
           do iion = 2, iSoluteIons
            do jion = 1, iSoluteIons-1
              distIJ = calcCMDist(iion,jion)
            enddo ! loop over jions           
           enddo ! loop over iions
           
           ibin = int ( distIJ/drI ) + 1

           if(ibin .lt. maxBin) then
            distI_i(ibin) = distI_i(ibin) + 1.d0
           endif   
              
        return
      end subroutine calcDIon

      function calcCMDist(res1,res2)
c    *******************************************************************
c    ** find the center of mass distance between two residues         **
c    **                                                               **
c    ** principal variables:                                          **
c    ** integer      res1,res2  indexes of the two residues           **
c    **                                                               **
c    ******************************************************************* 				  
		
          implicit none
          
          ! function type
          real :: calcCMDist
			
          ! Args
          integer, intent(in) :: res1,res2

          ! Locals
          integer :: ii, rr1, rr2
          integer :: rStart, rEnd
          real    :: cmx1, cmy1, cmz1 
          real    :: cmx2, cmy2, cmz2          
          real    :: massTot, mi
          real    :: ddX, ddY, ddZ
          real    :: dXYZ
c    *******************************************************************

             ! Find center of mass for residue 1
             massTot = 0.d0                   
             cmx1 = 0.d0
             cmy1 = 0.d0
             cmz1 = 0.d0             
             rStart = starts(res1)
             rr1  = resNum(res1)
             rEnd = rStart + params%numAtoms(rr1) - 1
             do ii = rStart, rEnd ! Loop over res1 atoms
             
              mi = qrms(ii)%m             
              massTot = massTot + mi
              
              cmx1 = cmx1 + mi*nodeConf(ii)%x
              cmy1 = cmy1 + mi*nodeConf(ii)%y
              cmz1 = cmz1 + mi*nodeConf(ii)%z
                      
             enddo ! end atom loop
             
             cmx1 = cmx1/massTot
             cmy1 = cmy1/massTot
             cmz1 = cmz1/massTot
             
             ! Find center of mass for residue 2
             massTot = 0.d0                   
             cmx2 = 0.d0             
             cmy2 = 0.d0             
             cmz2 = 0.d0             
             rStart = starts(res2)
             rr2  = resNum(res2)
             rEnd = rStart + params%numAtoms(rr2) - 1
             do ii = rStart, rEnd ! Loop over res1 atoms
             
              mi = qrms(ii)%m             
              massTot = massTot + mi
              
              cmx2 = cmx2 + mi*nodeConf(ii)%x
              cmy2 = cmy2 + mi*nodeConf(ii)%y
              cmz2 = cmz2 + mi*nodeConf(ii)%z
                      
             enddo ! end atom loop
             
             cmx2 = cmx2/massTot
             cmy2 = cmy2/massTot
             cmz2 = cmz2/massTot
          
             ddX = cmx1 - cmx2
             ddY = cmy1 - cmy2
             ddZ = cmz1 - cmz2
          
             dXYZ =sqrt(ddX**2+ddY**2+ddZ**2)
             calcCMDist = dXYZ
          return
        end function calcCMDist

      subroutine setDIons()
c    *******************************************************************
c    ** find the center of mass distance between two residues         **
c    **                                                               **
c    ** principal variables:                                          **
c    ** integer      res1,res2  indexes of the two residues           **
c    **                                                               **
c    ******************************************************************* 				  		
          implicit none
          

          ! Locals
          integer :: ii, jj, ios
          character (len=100) :: str1
c    *******************************************************************

            if (rank .eq. master) then
             iSoluteIons = 0
             open(unit=iDIinp, file=trim(strDionInp), status='old')
             do
              read (iDIinp, '(a)', iostat = ios ) str1
              if ( ios /= 0 ) then
                exit
               else
                iSoluteIons = iSoluteIons + 1
              end if
             enddo
             allocate(starts(iSoluteIons))
             allocate(resNum(iSoluteIons))             
             rewind(iDIinp)             

             do ii = 1, iSoluteIons
              read(iDIinp, *) resNum(ii)
             enddo
             close(iDIinp)
             
             do ii = 1, iSoluteIons ! Loop over all ions
              starts(ii) = 0
              do jj = 1, resNum(ii) - 1 
                starts(ii) = starts(ii) + params%numAtoms(jj) 
              enddo ! jj loop         
              starts(ii) = starts(ii) + 1
             enddo ! ii loop

            endif ! End check for master
            
            call MPI_BARRIER(world,ierr)

           call MPI_BCAST(iSoluteIons,1
     x               ,MPI_INTEGER,master,world,ierr)       ! BCAST all the number of solute ions      
           if (rank .ne. master) then
             allocate(starts(iSoluteIons))
             allocate(resNum(iSoluteIons))
           endif
           call MPI_BCAST(starts,iSoluteIons
     x               ,MPI_INTEGER,master,world,ierr)       ! BCAST all indices of solute ions      
           call MPI_BCAST(resNum,iSoluteIons
     x               ,MPI_INTEGER,master,world,ierr)       ! BCAST all indices of solute ions      
     
           distI_i(:) = 0.d0  ! initialize the individual distributions  
          return
        end subroutine setDIons

      subroutine dIonsFunc()
c ************************************************************
c This subroutine calculates the solute ion distribution  
c function on master from partial histograms
c ************************************************************
        implicit none
        integer :: ii
        
c ************************************************************
        
        distI(:) = 0.d0 ! Initialize the distributions
        
        call MPI_REDUCE(distI_i,distI,maxBinIon,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
     
        if (rank .eq. master) then
        do ii = 1, maxBinIon
            dI_norm(ii) = distI(ii)/real(nConf)
        enddo   

        ! open file and write out the distribution of solute ions
        open(unit=iDIout, file=trim(strDionOut), status='unknown')
        do ii = 1, maxBinIon
           write(iDIout,*) drI*ii, distI(ii), dI_norm(ii)
        enddo            
        close(iDIout)

        endif ! end check for master
        
        call MPI_BARRIER(world,ierr)
        
        return
      end subroutine dIonsFunc

      subroutine covarFunc()
c ************************************************************
c This subroutine calculates the covariance matrix 
c on master from partial matrices on the nodes
c ************************************************************
        implicit none
       
        ! Locals
        integer :: ii,jj 
        integer :: kk,ll

        integer :: iDofBy3
        integer :: nDofSquared     
        real, pointer, dimension(:,:) :: covarOut
        real*8, pointer, dimension(:,:) :: 
     x                                   covarDiag
        real    :: norm        
        real    :: cm, cm1, cm2
        real    :: cmTrace
        real    :: diagTrace
        integer :: iDUMMY, laINFO
        integer, parameter :: iWORK=50000
        real*8, pointer,dimension(:,:)  :: aMat
        real*8, pointer, dimension(:)  :: eig_valR
        real*8, pointer, dimension(:)  :: eig_valI
        real*8  :: laWORK(iWORK), laDUMMY
        double precision   :: sDiagTime,eDiagTime
c ************************************************************

        iDofBy3 = iPCADoF/3        
        norm = 0.010d0/real(nConf)  !!! converts from ang^2 to nm^2, as well
        cmTrace = 0.d0

        ! Initialize the total covariance matrix ...
        if (rank .eq. master) then

          allocate( covarOut(iDofBy3,iDofBy3) )
          allocate( covarDiag(iDofBy3,iDofBy3) )
          allocate( aMat(iDofBy3,iDofBy3)  )  
          allocate( eig_valR(iDofBy3) )
          allocate( eig_valI(iDofBy3) )
          covarMat(:,:) = 0.d0
          covarOut(:,:) = 0.d0
          eig_valR(:) = 0.d0
          eig_valI(:) = 0.d0
        endif

       if (useCUDA) then
#ifdef ENABLE_CUDA
           if(doSelect) then
             call get_covar_sel_gpu(covarNode)
           else
             call get_covar_gpu(covarNode)
           endif

             call MPI_BARRIER(world,ierr)
#endif
       endif        

       nDofSquared = iPCADoF*iPCADoF

       ! ... and collect them
        call MPI_REDUCE(covarNode,covarMat,nDofSquared,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
     
        if (rank .eq. master) then

        ! open file and write the atomic average positions
        open(unit=103, file='covar.dat', status='unknown')
        kk = 0
        do ii = 1, iPCADoF,3
         kk = kk + 1
         ll = 0
         do jj = 1, iPCADoF,3
           ll = ll + 1

           ! this line below will hopefully use a fused-multiply add cycle
           ! and in the summation below use only need 2 adds intead of three
!           cm = covarMat(ii,jj)*norm
!           cm1 = covarMat(ii+1,jj+1)*norm + cm
!           cm2 = covarMat(ii+2,jj+2)*norm

           cm = covarMat(ii,jj)
           cm1 = covarMat(ii+1,jj+1) + cm
           cm2 = covarMat(ii+2,jj+2)

           covarOut(kk,ll) = covarOut(kk,ll) +
     x     cm1 + cm2          
          enddo
        enddo            

        do ii = 1, iDofBy3
          do jj = 1, iDofBy3
          covarOut(ii,jj) = covarOut(ii,jj)*norm
          covarDiag(ii,jj) = real(covarOut(ii,jj),
     x                             kind=8)
          enddo
        enddo

        print *, " :: Diagonalizing Cov(dx_i,dx_j) :: "
      ! diagonalize kMat and store the eigenvectors in
      ! aMat(i=1,3 ; j) for each j-th eigenvector.
      ! 'eig_valR' has the eigenvalues.   
         iDUMMY = 1
         sDiagTime = MPI_WTIME()
         call dgeev('N','V',
     x            iDofBy3,covarDiag,
     x         iDofBy3,eig_valR,eig_valI,
     x              laDUMMY,iDUMMY,
     x        aMat,iDofBy3,laWORK,iWORK,laINFO)
         eDiagTime = MPI_WTIME()
        print *, " :: Diagonalizing COMPLETE! :: "
        diagTime = real( (eDiagTime - sDiagTime)
     x            ,kind=4)
        print *, " :: Time(diag) = ",
     x          diagTime
        print *, ""
        diagTrace = 0.d0
        do ii = 1, iDofBy3
          cmTrace = cmTrace + covarOut(ii,ii)
          diagTrace = diagTrace + eig_valR(ii)
        enddo
        print *, " Tr(covariance matrix) = ", cmTrace
        print *, " Tr(diag[covariance matrix]) = ", 
     x                                       diagTrace
        print *, " iDofBy3 = ", iDofBy3
        do ii = 1,iDofBy3

         do jj = 1,iDofBy3
          write(103,*) ii,jj,covarOut(ii,jj)
!          write(103,*) ii,jj,covarDiag(ii,jj)
         enddo
!          write(103,9922) (covarOut(ii,jj),jj=1, iDofBy3)
!          write(103,*) (covarOut(ii,jj),jj=1, iDofBy3)
!          write(103,*) (covarOut(ii,jj),jj=1, iPCADoF,3)
!          write(103,*) covarOut(ii,ii)
        enddo
        close(103)
9922   format (10F13.6)
        deallocate( covarOut ) 
        deallocate( covarDiag ) 
        deallocate( eig_valR )
        deallocate( eig_valI )

        endif ! end check for master


        
        call MPI_BARRIER(world,ierr)
        
        return
      end subroutine covarFunc

      subroutine avgSFunc()
c ************************************************************
c This subroutine calculates the average structure of each  
c atom on master from partial sums on the nodes
c ************************************************************
        implicit none
       
        ! Locals
        integer :: ii      
c ************************************************************

        ! Initialize the averages ...
        if (rank .eq. master) then
        xAvgTot(:) = 0.d0
        yAvgTot(:) = 0.d0
        zAvgTot(:) = 0.d0
        endif

       if (useCUDA) then
#ifdef ENABLE_CUDA
          if(doSelect) then
           call get_avg_sel_gpu(xAvg,yAvg,zAvg)
          else
           call get_avg_gpu(xAvg,yAvg,zAvg)
          endif
#endif
       endif        


        if(doSelect) then

        ! ... and collect them
         call MPI_REDUCE(xAvg,xAvgTot,nSelect,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
         call MPI_REDUCE(yAvg,yAvgTot,nSelect,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
         call MPI_REDUCE(zAvg,zAvgTot,nSelect,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
     
        if (rank .eq. master) then
         ! open file and write the atomic average positions
          open(unit=103, file='avgProt.dat', status='unknown')
          do ii = 1, nSelect
          xAvgTot(ii) = xAvgTot(ii)/real(nConf)
          yAvgTot(ii) = yAvgTot(ii)/real(nConf)
          zAvgTot(ii) = zAvgTot(ii)/real(nConf)           
           write(103,*) ii, 
     x        xAvgTot(ii), yAvgTot(ii), zAvgTot(ii)
          enddo            
          close(103)
         endif ! end check for master

        else ! no selection

        ! ... and collect them
         call MPI_REDUCE(xAvg,xAvgTot,iProtEnd,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
         call MPI_REDUCE(yAvg,yAvgTot,iProtEnd,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
         call MPI_REDUCE(zAvg,zAvgTot,iProtEnd,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
     
        if (rank .eq. master) then
         ! open file and write the atomic average positions
          open(unit=103, file='avgProt.dat', status='unknown')
          do ii = 1, iProtEnd
          xAvgTot(ii) = xAvgTot(ii)/real(nConf)
          yAvgTot(ii) = yAvgTot(ii)/real(nConf)
          zAvgTot(ii) = zAvgTot(ii)/real(nConf)           
           write(103,*) ii, 
     x        xAvgTot(ii), yAvgTot(ii), zAvgTot(ii)
          enddo            
          close(103)
         endif ! end check for master

        endif ! end check for selection

        call MPI_BARRIER(world,ierr)
        
        return
      end subroutine avgSFunc
      
      subroutine qdrFunc()
c ************************************************************
c This subroutine calculates <(qdr)^2>, where dr is the     
c fluctuation of the atomic position and q is the charge.
c ************************************************************
        implicit none
        integer :: ii, iatm
        real    :: qq, ddr2
        real    :: ddx, ddy, ddz
c ************************************************************
           
        ! Initialize the fluctuations ...
        if (rank .eq. master) then
        drxAvgTot(:) = 0.d0
        dryAvgTot(:) = 0.d0
        drzAvgTot(:) = 0.d0 
        qdrsqTot(:)  = 0.d0
        endif        

       ! ... and collect them
        call MPI_REDUCE(drxAvg,drxAvgTot,iProtEnd,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
        call MPI_REDUCE(dryAvg,dryAvgTot,iProtEnd,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
        call MPI_REDUCE(drzAvg,drzAvgTot,iProtEnd,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
        call MPI_REDUCE(qdrsq,qdrsqTot,iProtEnd,
     x          MPI_REAL,MPI_SUM,master,world,ierr)     
     
        if (rank .eq. master) then
!         drxAvgTot(:) = drxAvgTot(:)/real(nConf)
!         dryAvgTot(:) = dryAvgTot(:)/real(nConf)
!         drzAvgTot(:) = drzAvgTot(:)/real(nConf) 
         
        ! open file and write the fluctuation data
        open(unit=53, file='qdr.dat', status='unknown')
        do ii = 1, iProtEnd
         qq = qrms(ii)%q
         ddr2 = qdrsqTot(ii)/real(nConf)         
         qdrsqTot(ii) = (qq*qq*qdrsqTot(ii))/real(nConf) 
         !ddx = drxAvgTot(ii)
         !ddy = dryAvgTot(ii)
         !ddz = drzAvgTot(ii)         
         !ddr2 = ddx**2+ ddy**2+ddz**2

         
           write(53,*) ii, qdrsqTot(ii), qq, ddr2
!     x        drxAvgTot(ii), dryAvgTot(ii), drzAvgTot(ii)
        enddo            
        close(53)

        endif ! end check for master
        
        call MPI_BARRIER(world,ierr)
        
        return
      end subroutine qdrFunc
      
      subroutine bindFunc()
c ************************************************************
c This subroutine calculates the water binding energy  
c  distribution function on master from partial histograms
c ************************************************************
        implicit none
        integer :: ii
        real    :: rWatTot
        integer :: nRes
        
c ************************************************************
        
!         nRes = params%numRes
          if(rank .eq. master) nRes = params%numRes
          
          call MPI_BCAST(nRes,1,MPI_INTEGER,
     x                     master,world,ierr)  

         
        if(rank .eq. master) then
        
         allocate(nHighAA(nRes))         
         allocate(nLowAA(nRes))         
         allocate(eTotAA(nRes))         
         allocate(eDelAA(nRes))    
         
         distB(:) = 0.d0 ! Initialize the distributions
         distBde(:) = 0.d0 ! Initialize the distributions

     
         ! number of high/low energy waters on all amino acids         
         nHighAA(:) = 0.d0
         nLowAA(:) = 0.d0
         
         ! total binding energies on all amino acids
         eTotAA(:) = 0.d0
         eDelAA(:) = 0.d0         
        endif

c      Per-amino acid based Stuff        
        do ii = 1, nRes
         rWatTot = nShellLowAA(ii) + nShellHighAA(ii)
         ebindTotAA(ii) = ebindTotAA(ii)/rWatTot
         ebindDelAA(ii) = ebindDelAA(ii)/rWatTot         
        enddo
        
        call MPI_REDUCE(ebindTotAA,eTotAA,nRes,
     x          MPI_REAL,MPI_SUM,master,world,ierr) ! total binding energy on each AA

        call MPI_REDUCE(ebindDelAA,eDelAA,nRes,
     x          MPI_REAL,MPI_SUM,master,world,ierr) ! delta binding energy on each AA
     
        call MPI_REDUCE(nShellHighAA,nHighAA,nRes,
     x          MPI_REAL,MPI_SUM,master,world,ierr) ! number of high energy waters on each AA

        call MPI_REDUCE(nShellLowAA,nLowAA,nRes,
     x          MPI_REAL,MPI_SUM,master,world,ierr) ! number of low energy waters on each AA


c      Distribution Stuff
        call MPI_REDUCE(distB_i,distB,maxBinBind,
     x          MPI_REAL,MPI_SUM,master,world,ierr) ! full

        call MPI_REDUCE(distBde_i,distBde,maxBinBindDe,
     x          MPI_REAL,MPI_SUM,master,world,ierr) ! delta energy
     
        if (rank .eq. master) then

c      Per-amino acid based Stuff            
        do ii = 1, nRes
         eTotAA(ii) = eTotAA(ii)/real(nConf)
         eDelAA(ii) = eDelAA(ii)/real(nConf)
         nHighAA(ii) = nHighAA(ii)/real(nConf)
         nLowAA(ii) = nLowAA(ii)/real(nConf)
        enddo
        
       ! open file and write out the AA info
        open(unit=53, file='enerBindAA.dat', status='unknown')       
        do ii = 1, nRes
           write(53,*) ii,
     x                eTotAA(ii), eDelAA(ii),
     x                nLowAA(ii), nHighAA(ii)     
        enddo            
        close(53)

         deallocate(nHighAA)         
         deallocate(nLowAA)         
         deallocate(eTotAA)         
         deallocate(eDelAA)   
         
        !normalize the delta energy distros.
        do ii = -maxBinBindDe/2, maxBinBindDe/2
            dBde_norm(ii) = distBde(ii)/real(nConf)
        enddo   

        !normalize the full energy distros.        
        do ii = -maxBinBind/2, maxBinBind/2
            dB_norm(ii) = distB(ii)/real(nConf)
        enddo           

        ! open file and write out the distribution of solute ions (delta)
        open(unit=53, file=trim(strBindDeFile), status='unknown')       
        do ii = -maxBinBindDe/2, maxBinBindDe/2
           write(53,*) drBde*ii,!(drBde*ii - enerShift), 
     x                    distBde(ii), dBde_norm(ii)
        enddo            
        close(53)

        ! open file and write out the distribution of solute ions (total)
        open(unit=53, file=trim(strBindFile), status='unknown')
        do ii = -maxBinBind/2, maxBinBind/2
           write(53,*) drB*ii,!(drB*ii - enerShift), 
     x                    distB(ii), dB_norm(ii)
        enddo            
        close(53)

        endif ! end check for master
        
        call MPI_BARRIER(world,ierr)
        
        return
      end subroutine bindFunc

      subroutine grFunc()
c ************************************************************
c This subroutine calculates the radial distribution function 
c from a given histogram
c ************************************************************
        implicit none
        integer :: ii
        real    :: nideal, rl, ru, pi, const
        real    :: totVol
        
c ************************************************************
        pi      = 4.*atan( 1.0 )
        const   = 4.0 * pi / 3.0

        !!!!WARNINGWARNINGWARNINGWARNINGWARNING!!!
        !!! MUST BE CHANGED -- HARDCODED FOR PC !!
        !!!!WARNINGWARNINGWARNINGWARNINGWARNING!!!
        totVol = (bX*bY*bZ - const*16.**3)     !!! 
        !!!!WARNINGWARNINGWARNINGWARNINGWARNING!!!
        
        if (rank .eq. master) then
         close(iXOfile)
        endif 
        
          call MPI_REDUCE(hist_i,hist,maxBin,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
     
        if (rank .eq. master) then
        do ii = 1, maxBin
            rl     = real (ii - 1) * delr
            ru     = real (ii) * delr
            nideal = const * (ru**3 - rl**3)*nWaters/totVol
            gr000(ii) = hist(ii)/real(nConf)/nideal
        enddo   

        ! open file and write out gr(Cu-O)
        open(unit=53, file=trim(strGrFile), status='unknown')
        !open(unit=87, file='grCu-O.raw', status='unknown')                   
        do ii = 1, maxBin
           write(53,*) delr*ii, gr000(ii)
        !   write(87,*) delr*ii, hist(ii)                      
        enddo            
        close(53)
        !close(87)
        endif ! end check for master
        
        call MPI_BARRIER(world,ierr)
        
        return
      end subroutine grFunc      

      subroutine efsFunc()
c ************************************************************
c This subroutine calculates the vibrational OH spectrum 
c from a set of e-field projections along the OH bonds
c ************************************************************
        implicit none
        integer :: ii
        
c ************************************************************
        
        
        call MPI_REDUCE(histEFS_i,histEFS,maxBinEFS,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
     
        if (rank .eq. master) then

        ! open file and write out P(efs)
        open(unit=1980, file='hisEFS.dat', status='unknown')                   
        do ii = 1, maxBinEFS
           write(1980,*) delEFS*ii - efsShift,
     x                    histEFS(ii)
        enddo            
        close(1980)

        endif ! end check for master
        
        call MPI_BARRIER(world,ierr)
        
        return
      end subroutine efsFunc     

      subroutine prMuFunc()
c ************************************************************
c This subroutine calculates the P(cos(theta),rcut) function 
c from a given histogram
c ************************************************************
        implicit none
        integer :: ii
        real    :: nideal, rl, ru, pi, const
        real    :: cosIntegral
        
c ************************************************************
        cosIntegral = 0.d0        
        
        call MPI_REDUCE(histCos_i,histCos,maxBinCos,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
     
        if (rank .eq. master) then
        do ii = -maxBinCos, maxBinCos
          cosIntegral = cosIntegral + 
     x                  histCos(ii)*dCosT   
        enddo   

        do ii = -maxBinCos, maxBinCos
          Pcos(ii) = histCos(ii)/real(cosIntegral,kind=4)
        enddo
        ! open file and write out P(cos(theta),rCut)
        open(unit=85, file='Pcos_pw.dat', status='unknown')
        !open(unit=87, file='grCu-O.raw', status='unknown')                   
        do ii = -maxBinCos, maxBinCos
           write(85,*) dCosT*ii, Pcos(ii), histCos(ii)
        !   write(87,*) delr*ii, hist(ii)                      
        enddo            
        close(85)
        !close(87)
        endif ! end check for master
        
        call MPI_BARRIER(world,ierr)
        
        return
      end subroutine prMuFunc      

      subroutine calcEnerSolu()
c    *******************************************************************
c    ** calculates solvent (+ ions) energetics at every timestep      **
c    ** using all solute (protein) atoms.                             **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: iw, irc, jj, iWrite
        real    :: ddX, ddY, ddZ, ddQ        
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then 
          open(unit=iSolu,file="enerSolu.dat",status="unknown")               
         endif
        endif 
        
        dv0s = 0.d0
        dv0i = 0.d0
c---------------------------------------------------------c         
         do irc = 1, iProtEnd      
c---------------------------------------------------------c
            ddX = nodeConf(irc)%x
            ddY = nodeConf(irc)%y
            ddZ = nodeConf(irc)%z
            ddQ = qrms(irc)%q
            
          if (doSolv) then
           if (useCUDA) then
#ifdef ENABLE_CUDA
             call dvos_gpu(ddX,ddY,ddZ,ddQ,
     x                     bX,bY,bZ,nodeConf,qrms,dv0s_i)
#endif
           else
             call mu0s_direct(dv0s_i,
     x            ddX, ddY, ddZ, ddQ)                  
           endif ! end check for cuda           
            dv0s = dv0s + dv0s_i                  
          endif ! end test for solvent calculations

         if (doIons) then
           call mu0i_direct(dv0i_i,
     x          ddX, ddY, ddZ, ddQ)
     
           dv0i = dv0i + dv0i_i                      
         endif ! end test for ion energy calculations
         
c---------------------------------------------------------c          
         enddo ! end loop over all protein atoms
c---------------------------------------------------------c         

        if (doSolv) then         
         call MPI_GATHER(real(dv0s,kind=4),1,MPI_REAL,
     x                     dv0sTraj,1,MPI_REAL,
     x                     master,world,ierr)     
        endif ! End test for solvent calculation

        if (doIons) then                
          call MPI_GATHER(real(dv0i,kind=4),1,MPI_REAL,
     x                     dv0iTraj,1,MPI_REAL,
     x                     master,world,ierr)      
        endif ! End test for protein calc  
                   
          if(rank .eq. master)  then
           
            do iWrite = 1, nProc
              if (doSolv .and. doIons) then
               write(iSolu, *) iConf+iWrite-1, 
     x         dv0sTraj(iWrite), dv0iTraj(iWrite)
              elseif(doIons) then
               write(iSolu, *) iConf+iWrite-1, 
     x         dv0iTraj(iWrite)     
              else ! solvent only
               write(iSolu, *) iConf+iWrite-1, 
     x         dv0sTraj(iWrite)
              endif ! end check for doProt/doSolv
            enddo
           endif ! End check for master       

        call MPI_BARRIER(world,ierr)
       
        return
      end subroutine calcEnerSolu

      subroutine calcEnerShel()
c    *******************************************************************
c    ** calculates 1st shell (total and vertical energies) at every   **
c    ** timestep.  It also makes a histogram of the results.          **
c    **                                                               **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: iw, iw2, irc, jj, iWrite
        integer :: iAtm, irc2
        integer :: iHigh, iLow
        integer :: nShell, ibin
        integer :: lShell(25000)  ! list of first atoms of water
        integer :: ix, iy, iz        
        integer :: xInt, yInt, zInt
        integer :: iwCnt, iLowAA        
        real    :: ddX, ddY, ddZ, ddQ
        real*8 :: prX, prY, prZ, prQ, prR
        real    :: dX, dY, dZ, myQ         
        real*8  :: dQ, rX, rY, rZ, rIJ
        real*8  :: pQ, pX, pY, pZ
        real*8  :: wQ, wX, wY, wZ        
        real*8  :: dfrx, dfry, dfrz, delr2   
        real*8  :: eShift, lowRIJ        
        logical :: firstShell, doAtom
        type(basic_atom) :: dvAtoms(iNumDB)
        real    :: dVosQ(iNumDB), dVFull, dVos_only
c        real    :: ttQ
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then 
          strBindFile = 'distBind.dat' !!!!! <-- Should be from input
          strBindDeFile = 'distBindDe.dat' !!!!! <-- Should be from input          
          open(unit=iShell,file="enerShel.dat",status="unknown")               
         endif

        endif ! End test for 1st configuration 
        
         nShell = 0        

         if(useCuda) then 
#ifdef ENABLE_CUDA
            call get_shell_gpu(bX,bY,bZ,
     x                     nodeConf,qrms,
     x                    rBindCut,lShell,nShell,.false.)

c           ttQ = 0.d0
            ! Loop over all RC atoms, and fill up the dvAtoms array
          do irc = 1, iNumDB 
             iatm = dVosDbase(irc)%atom
             dvAtoms(irc) = nodeConf(iatm)
             dVosQ(irc) = dVosDbase(irc)%q

c             ttQ = ttQ + dVosQ(irc)
          enddo                  

c           print *, "IN FORTRAN ::: total dQ = ", ttQ, 
c     x              " for ", iNumDB, " atoms"

         ! Calculate the total binding/shell energy
         ! and the dV0s shell energy
            call get_dv0s_shl_gpu(bX,bY,bZ,
     x                            nodeConf,qrms,
     x                            dvAtoms, dVosQ, 
     x                            lShell,nShell,
     x                            dVFull, dVos_only)      

            dv0sh = real(dVFull,kind=8)     
            dv0shde = real(dVos_only,kind=8)
#endif /* CUDA PART */
         else
         
         do iw = iWatStart, natom, 3    ! Loop over all water-oxygen atoms
               
          wX = real(nodeConf(iw)%x,kind=8)                     
          wY = real(nodeConf(iw)%y,kind=8)
          wZ = real(nodeConf(iw)%z,kind=8)          
         
          ! Find the 1st solvation shell waters
          firstShell = isFirstShell(nodeConf(iw)%x,
     x                              nodeConf(iw)%y,
     x                              nodeConf(iw)%z)
     
          ! Store the index of the oxygen of the 1st shell
          if(firstShell) then
           nShell = nShell + 1
           lShell(nShell) = iw         
          endif
          
         enddo ! End the loop over all waters     

         dv0sh_i = 0.d0
         dv0sh   = 0.d0
         dv0shde_i = 0.d0
         dv0shde   = 0.d0    
         dv0bsh  = 0.d0
         dv0bshde = 0.d0
         nShellLow = 0
         nShellHigh = 0
          do jj = 1, nShell ! Loop over all 1st shell waters          
          iw = lShell(jj)
          
          dv0bsh_i  = 0.d0          
          dv0bshde_i  = 0.d0      
          lowRIJ = 9999999999.99
          iwCnt = 0
           do iw2 = iw, iw+2 ! Loop over all 3 atoms
            iwCnt = iwCnt + 1
            wX = real(nodeConf(iw2)%x,kind=8)
            wY = real(nodeConf(iw2)%y,kind=8)
            wZ = real(nodeConf(iw2)%z,kind=8)
            wQ = real(qrms(iw2)%q,kind=8)
            
           do irc = 1, iProtEnd  ! Loop over protein atoms
            prX = real(nodeConf(irc)%x,kind=8)                     
            prY = real(nodeConf(irc)%y,kind=8)
            prZ = real(nodeConf(irc)%z,kind=8)
            prQ = real(qrms(irc)%q,kind=8)            

            !*******************         
         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then 
          dX = wX - prX
          dY = wY - prY
          dZ = wZ - prZ

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          rIJ = dX**2 + dY**2 + dZ**2
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
            endif

            enddo
           enddo
          enddo
          rIJ = sqrt(rIJ)
         else ! its an orthorhombic box
         
          xInt = ANINT((wX-prX)/bX)
          yInt = ANINT((wY-prY)/bY)
          zInt = ANINT((wZ-prZ)/bZ)

          dX = (wX-prX) - xInt*bX
          dY = (wY-prY) - yInt*bY
          dZ = (wZ-prZ) - zInt*bZ         
          rIJ = sqrt(dX*dX+dY*dY+dZ*dZ)

         endif ! End check for truncation            
            !*******************
            

            if(iwCnt .eq. 1) then             
             if(rIJ .lt. lowRIJ) then
               lowRIJ = rIJ ! Make the lowest
               iLowAA = aminoAcidList(irc) ! use the aa-lookup table to find the amino acid
             endif
            endif

           ! Do a calculation for the energy of the solute
           dv0bsh_i = dv0bsh_i + wQ*prQ/rIJ 
           
           ! Make sure you do a separate deltaE calc. for rc atoms
            call searchDBAtom(irc,myQ,doAtom)
             if (doAtom) then

             prQ = real(myQ,kind=8)
             dv0bshde_i = dv0bshde_i + wQ*prQ/rIJ                          
            endif

           
            enddo ! End loop over protein atoms   
      
           enddo ! End loop over all 3 water atoms

          ! Determine if the energy is high or low...          
          if (dv0bsh_i .lt. bindCutoff) then ! ... and bin the results
           nShellLow = nShellLow + 1
           nShellLowAA(iLowAA) = nShellLowAA(iLowAA) + 1.d0
          else
           nShellHigh = nShellHigh + 1
           nShellHighAA(iLowAA) = nShellHighAA(iLowAA) + 1.d0
          endif          
           
          ! Accumulate the full/delta energy to the closest AA
          ebindTotAA(iLowAA) = ebindTotAA(iLowAA) + dv0bsh_i
          ebindDelAA(iLowAA) = ebindDelAA(iLowAA) + dv0bshde_i
          
          ! Bin the full and delta binding energies
          dv0sh = dv0sh + dv0bsh_i
          dv0shde = dv0shde + dv0bshde_i
          
           ! Bin the full binding energy distros
           ibin = int(dv0bsh_i/drB ) + 1
           if(abs(ibin) .le. maxBinBind/2) then
            distB_i(ibin) = distB_i(ibin) + 1.d0
           endif

          ! Bin the delta binding energy distros
           ibin = int(dv0bshde_i/drBde) + 1
           if(abs(ibin) .le. maxBinBindDe/2) then
            distBde_i(ibin) = distBde_i(ibin) + 1.d0
           endif
        
          enddo ! End loop over 1st shell waters         
        endif ! End check for CUDA
 
          dv0bsh = dv0sh/real(nShell)
          dv0bshde = dv0shde/real(nShell)           
          
         call MPI_GATHER(real(dv0bsh,kind=4),1,MPI_REAL,
     x                     dv0bshTraj,1,MPI_REAL,
     x                     master,world,ierr) 
     
         call MPI_GATHER(real(dv0sh,kind=4),1,MPI_REAL,
     x                     dv0shTraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0bshde,kind=4),1,MPI_REAL,
     x                     dv0bdeshTraj,1,MPI_REAL,
     x                     master,world,ierr) 
     
         call MPI_GATHER(real(dv0shde,kind=4),1,MPI_REAL,
     x                     dv0shdeTraj,1,MPI_REAL,
     x                     master,world,ierr)  

         call MPI_GATHER(nShell,1,MPI_INTEGER,
     x                   nShellTraj,1,MPI_INTEGER,
     x                   master,world,ierr) 

c         call MPI_GATHER(nShellLow,1,MPI_INTEGER,
c     x                   nShellLowTraj,1,MPI_INTEGER,
c     x                   master,world,ierr) 

c         call MPI_GATHER(nShellHigh,1,MPI_INTEGER,
c     x                   nShellHighTraj,1,MPI_INTEGER,
c     x                   master,world,ierr) 
     
          if(rank .eq. master)  then           
            do iWrite = 1, nProc
c               iHigh = int(nShellHighTraj(iWrite))
c               iLow = int(nShellLowTraj(iWrite))
c               iHigh = nShellHighTraj(iWrite)
c               iLow = nShellLowTraj(iWrite)

               write(iShell, 8888) iConf+iWrite-1,                
     x         dv0bshTraj(iWrite), dv0shTraj(iWrite), 
     x         dv0bdeshTraj(iWrite), dv0shdeTraj(iWrite),
     x         nShellTraj(iWrite) 
c     x         iLow, iHigh
c     x         nShellLowTraj(iWrite),nShellHighTraj(iWrite)
            enddo
           endif ! End check for master       

!801             format (5x,i3,5x,i5,11x,f6.3)           
c8888             format (i6,2x,4(f15.12,2x),2(i4,2x))
8888             format (i6,2x,4(f15.12,2x),i6)                                        
        return
      end subroutine calcEnerShel


      subroutine calcVpolShel()
c    *******************************************************************
c    ** calculates 1st shell V_pol correction for all shell waters    **
c    ** at every timestep.                                            **
c    **                                                               **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: iw, iw2, irc, jj, iWrite
        integer :: iAtm, irc2
        integer :: iHigh, iLow
        integer :: nShell, ibin
        integer :: lShell(25000)  ! list of first atoms of water
        integer :: ix, iy, iz        
        integer :: xInt, yInt, zInt
        integer :: iwCnt, iLowAA        
        real    :: ddX, ddY, ddZ, ddQ
        real*8 :: prX, prY, prZ, prQ, prR
        real    :: dX, dY, dZ, myQ         
        real*8  :: dQ, rX, rY, rZ, rIJ
        real*8  :: pQ, pX, pY, pZ
        real*8  :: wQ, wX, wY, wZ        
        real*8  :: dfrx, dfry, dfrz, delr2
        real*8  :: pi, dVpol_i
        real*8  :: dVpol,rij2,twoPiByThreeV    
        real    :: twoPiByThreeVTraj(nProc)    
        logical :: firstShell, doAtom
        type(basic_atom) :: dvAtoms(iNumDB)
        real    :: dVosQ(iNumDB), dVpol32
c    ******************************************************************* 

        if (iConf .eq. 1 ) then

         if ( rank .eq. master) then 
          open(unit=iVpol,file="enerVpol.dat",status="unknown")               
         endif

        endif ! End test for 1st configuration 

         pi = 4.*atan( 1.0 )        
         nShell = 0        
         ! find the 2pi/3V constant
         twoPiByThreeV = (-2.d0*pi)/(3.d0*real(volCell,kind=8))

         if(useCuda) then 
#ifdef ENABLE_CUDA
            call get_shell_gpu(bX,bY,bZ,
     x                     nodeConf,qrms,
     x                    rBindCut,lShell,nShell,.false.)

            ! Loop over all RC atoms, and fill up the dvAtoms array
          do irc = 1, iNumDB 
             iatm = dVosDbase(irc)%atom
             dvAtoms(irc) = nodeConf(iatm)
             dVosQ(irc) = dVosDbase(irc)%q
          enddo                  

         ! Calculate -Vpol*3V/2Pi
            call get_vpol_shl_gpu(bX,bY,bZ,
     x                            nodeConf,qrms,
     x                            dvAtoms, dVosQ, 
     x                            lShell,nShell,
     x                            dVpol32)      

          ! and then mulitply the correction factor
            dVpol = real(dVpol32,kind=8)
            dVpol = twoPiByThreeV * dVpol
     
#endif /* CUDA PART */
         else
         
         do iw = iWatStart, natom, 3    ! Loop over all water-oxygen atoms
               
          wX = real(nodeConf(iw)%x,kind=8)                     
          wY = real(nodeConf(iw)%y,kind=8)
          wZ = real(nodeConf(iw)%z,kind=8)          
         
          ! Find the 1st solvation shell waters
          firstShell = isFirstShell(nodeConf(iw)%x,
     x                              nodeConf(iw)%y,
     x                              nodeConf(iw)%z)
     
          ! Store the index of the oxygen of the 1st shell
          if(firstShell) then
           nShell = nShell + 1
           lShell(nShell) = iw         
          endif
          
         enddo ! End the loop over all waters     

   
          dVpol_i  = 0.d0          
          do jj = 1, nShell ! Loop over all 1st shell waters          
          iw = lShell(jj)
          
           do iw2 = iw, iw+2 ! Loop over all 3 atoms
            iwCnt = iwCnt + 1
            wX = real(nodeConf(iw2)%x,kind=8)
            wY = real(nodeConf(iw2)%y,kind=8)
            wZ = real(nodeConf(iw2)%z,kind=8)
            wQ = real(qrms(iw2)%q,kind=8)
            
           do irc = 1, iProtEnd  ! Loop over protein atoms
            prX = real(nodeConf(irc)%x,kind=8)                     
            prY = real(nodeConf(irc)%y,kind=8)
            prZ = real(nodeConf(irc)%z,kind=8)
            prQ = real(qrms(irc)%q,kind=8)            

            !*******************         
         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then 
          dX = wX - prX
          dY = wY - prY
          dZ = wZ - prZ

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          rIJ = dX**2 + dY**2 + dZ**2
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
              rij2 = delr2
            endif

            enddo
           enddo
          enddo
          rIJ = sqrt(rIJ)
         else ! its an orthorhombic box
         
          xInt = ANINT((wX-prX)/bX)
          yInt = ANINT((wY-prY)/bY)
          zInt = ANINT((wZ-prZ)/bZ)

          dX = (wX-prX) - xInt*bX
          dY = (wY-prY) - yInt*bY
          dZ = (wZ-prZ) - zInt*bZ

          rij2 = dX*dX+dY*dY+dZ*dZ         
          rIJ = sqrt(rij2)

         endif ! End check for truncation            
            !*******************
            
           ! Do a calculation for Vpol of the deltaQ atoms only           
            call searchDBAtom(irc,myQ,doAtom)
            if (doAtom) then
             prQ = real(myQ,kind=8)
             dVpol_i = dVpol_i + prQ*wQ*rij2                          
            endif

           
            enddo ! End loop over protein atoms   
      
           enddo ! End loop over all 3 water atoms
                   
          enddo ! End loop over 1st shell waters  

          !dVpol = dVpol_i
          dVpol = twoPiByThreeV * dVpol_i
       
        endif ! End check for CUDA
           
         call MPI_GATHER(real(dVpol,kind=4),1,MPI_REAL,
     x                     dVpolTraj,1,MPI_REAL,
     x                     master,world,ierr) 
     
         call MPI_GATHER(real(twoPiByThreeV,kind=4),1,
     x                     MPI_REAL,
     x                twoPiByThreeVTraj,1,MPI_REAL,
     x                     master,world,ierr) 

         call MPI_GATHER(nShell,1,MPI_INTEGER,
     x                   nShellTraj,1,MPI_INTEGER,
     x                   master,world,ierr) 

     
          if(rank .eq. master)  then           
            do iWrite = 1, nProc
               write(iVpol, 8989) iConf+iWrite-1,                
     x         dVpolTraj(iWrite), 
     x    twoPiByThreeVTraj(iWrite), nShellTraj(iWrite) 
            enddo
           endif ! End check for master       

8989             format (i6,2x,2(f15.12,2x),i6)                                        
        return
      end subroutine calcVpolShel

      subroutine calcQpolShel()
c    *******************************************************************
c    ** calculates 1st shell Q_pol correction for all shell waters    **
c    ** at every timestep.                                            **
c    **                                                               **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: iw, iw2, jj, iWrite
        integer :: nShell
        integer,parameter :: maxShellN=30000
        integer :: lShell(maxShellN)  ! list of first atoms of water
        integer :: ix, iy, iz        
        integer :: xInt, yInt, zInt     
        integer :: iii
        real    :: rIJ
        real    :: Qpol32
        real    :: qS(3)
        real*8  :: Qpol, Qpol_i
        logical :: firstShell, doAtom, keepFinal
        type(basic_atom)  :: initWater(3), finalWater(3), dWat
        type(basic_atom)  :: watM, watOxyR
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then          
          open(unit=iQpol,file="Qpol.dat",status="unknown")               
         endif
        endif 
        
         nShell = 0        
         
         if(useCuda) then
#ifdef ENABLE_CUDA
            call get_shell_gpu(bX,bY,bZ,
     x                     nodeConf,qrms,
     x                    rBindCut,lShell,nShell,.false.)

         ! Calculate Q(t)
            call get_qpol_shl_gpu(bX,bY,bZ,
     x                            nodeConf,qrms,
     x                            lShell,nShell,
     x                            Qpol32)      

          ! and then mulitply the correction factor
            Qpol = real(Qpol32,kind=8)

#endif
         else
!!!!!!!!!!!!!! start cpu-only part !!!!!!!!!!
            do iw = iWatStart, natom, 3    ! Loop over all water-oxygen atoms
                        
             ! Find the 1st solvation shell waters
             firstShell = isFirstShell(nodeConf(iw)%x,
     x                              nodeConf(iw)%y,
     x                              nodeConf(iw)%z)
     
             ! Store the index of the oxygen of the 1st shell
             if(firstShell) then
              nShell = nShell + 1
              lShell(nShell) = iw         
             endif
          
           enddo ! End the loop over all waters     

           ! Initialize vars
           Qpol = 0.d0
          do jj = 1, nShell ! Loop over all 1st shell waters          
           iw = lShell(jj) ! Get shell water index
                    
           iii = 0
           do iw2 = iw, iw+2 ! Loop over all 3 atoms
           
           iii = iii + 1           
           initWater(iii) = nodeConf(iw2) ! Store the initial coordinates of h2o atoms
           qS(iii) = qrms(iw2)%q          ! Store the charge of the h2o atoms
           
           call imageCoords(initWater(iii)%x,
     x          initWater(iii)%y,initWater(iii)%z,
     x          finalWater(iii)%x,finalWater(iii)%y,
     x     finalWater(iii)%z)
                                     
      
           enddo ! End loop over all 3 water atoms

           keepFinal = .true.
           do iii = 2, 3 ! Loop over both hydrogens
           
            dWat = finalWater(1) - finalWater(iii)                                          
            rIJ = sqrt(dWat%x*dWat%x + dWat%y*dWat%y + dWat%z*dWat%z)
 
            if(rIJ .lt. 2.5d0) keepFinal = .false.
           
           enddo ! End loop over both hydros

           
          if(keepFinal) then
            watOxyR = finalWater(1)
          else
            watOxyR = initWater(1)
          endif
         
          ! Calculate dipole moment
           watM%x = 0.d0
           watM%y = 0.d0
           watM%z = 0.d0

           do iii = 1, 3
            if(keepFinal) then

              watM%x = watM%x + qS(iii)*finalWater(iii)%x
              watM%y = watM%y + qS(iii)*finalWater(iii)%y
              watM%z = watM%z + qS(iii)*finalWater(iii)%z              
              
            else
           
              watM%x = watM%x + qS(iii)*initWater(iii)%x
              watM%y = watM%y + qS(iii)*initWater(iii)%y
              watM%z = watM%z + qS(iii)*initWater(iii)%z              
          
            endif
          
           enddo ! End loop of all 3 water atoms
           
           ! now calculate Q((t) = Rm.mu(t)
           Qpol_i = real(watOxyR*watM,kind=8)
           Qpol = Qpol + Qpol_i
          enddo ! End loop over 1st shell waters

       endif  ! End check for useCuda

         call MPI_BARRIER(world,ierr)
         ! gather all output data
         call MPI_GATHER(real(Qpol,kind=4),1,MPI_REAL,
     x                   dQpolTraj,1,MPI_REAL,
     x                   master,world,ierr) 


         call MPI_GATHER(nShell,1,MPI_INTEGER,
     x                   nShellTraj,1,MPI_INTEGER,
     x                   master,world,ierr) 

     
          if(rank .eq. master)  then           
            do iWrite = 1, nProc
               
               write(iQpol, 7678) iConf+iWrite-1, 
     x         dQpolTraj(iWrite), nShellTraj(iWrite)
     
            enddo 
           endif ! End check for master       
          
7678             format (i6,2x,E14.6,2x,i6)                    
        return

      end subroutine calcQpolShel


      subroutine calcEfieldShel()
c    *******************************************************************
c    ** calculates 1st shell electric field due to all energy gap     **
c    ** charges (aka - any probe charges) at every timestep.          **
c    **                                                               **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: iw, iw2, irc, jj, iWrite
        integer :: iAtm, irc2
        integer :: nShell, ibin
        integer :: lShell(25000)  ! list of first atoms of water
        integer :: ix, iy, iz        
        integer :: xInt, yInt, zInt
        integer :: iwCnt, iLowAA        
        real    :: ddX, ddY, ddZ, ddQ
        real*8 :: prX, prY, prZ, prQ, prR
        real    :: dX, dY, dZ, myQ         
        real*8  :: dQ, rX, rY, rZ, rIJ
        real*8  :: pQ, pX, pY, pZ
        real*8  :: wQ, wX, wY, wZ 
        real*8  :: efX0,efY0,efZ0
        real*8  :: efX,efY,efZ
        real*8  :: efMag2, efMag, efield_dV       
        real*8  :: dfrx, dfry, dfrz, delr2   
        real*8  :: rcDQ, r3dq, rinv, rinv3       
        real*8  :: lowDX, lowDY, lowDZ           
        logical :: firstShell, doAtom
        type(basic_atom) :: dvAtoms(iNumDB)
        real    :: dVosQ(iNumDB), dVEfield
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then        
          open(unit=iEfShel,file="efieldShel.dat",status="unknown")               
         endif

        endif ! End test for 1st configuration 
        
         nShell = 0        

         if(useCuda) then 
#ifdef ENABLE_CUDA
            call get_shell_gpu(bX,bY,bZ,
     x                     nodeConf,qrms,
     x                    rBindCut,lShell,nShell,.false.)

c           ttQ = 0.d0
            ! Loop over all RC atoms, and fill up the dvAtoms array
          do irc = 1, iNumDB 
             iatm = dVosDbase(irc)%atom
             dvAtoms(irc) = nodeConf(iatm)
             dVosQ(irc) = dVosDbase(irc)%q

c             ttQ = ttQ + dVosQ(irc)
          enddo                  

c           print *, "IN FORTRAN ::: total dQ = ", ttQ, 
c     x              " for ", iNumDB, " atoms"

         ! Calculate the total binding/shell energy
         ! and the dV0s shell energy

            call get_efld_shl_gpu(bX,bY,bZ,
     x                            nodeConf,qrms,
     x                            dvAtoms, dVosQ, 
     x                            lShell,nShell,
     x                            dVEfield)      

        efield_dV = real(dVEfield,kind=8)

#endif /* CUDA PART */
         else
         
         do iw = iWatStart, natom, 3    ! Loop over all water-oxygen atoms
               
          wX = real(nodeConf(iw)%x,kind=8)                     
          wY = real(nodeConf(iw)%y,kind=8)
          wZ = real(nodeConf(iw)%z,kind=8)          
         
          ! Find the 1st solvation shell waters
          firstShell = isFirstShell(nodeConf(iw)%x,
     x                              nodeConf(iw)%y,
     x                              nodeConf(iw)%z)
     
          ! Store the index of the oxygen of the 1st shell
          if(firstShell) then
           nShell = nShell + 1
           lShell(nShell) = iw         
          endif
          
         enddo ! End the loop over all waters     

         efield_dV = 0.d0
          do jj = 1, nShell ! Loop over all 1st shell waters          
          iw = lShell(jj)
          
          ! get the oxygen's position
          wX = real(nodeConf(iw)%x,kind=8)
          wY = real(nodeConf(iw)%y,kind=8)
          wZ = real(nodeConf(iw)%z,kind=8)      

         efX0 = 0.d0
         efY0 = 0.d0
         efZ0 = 0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do irc = 1, iNumDB  ! Loop over all vertical energy charges
          iAtm = dVosDbase(irc)%atom
          dQ = real(dVosDbase(irc)%q,kind=8)
          rX = real(nodeConf(iAtm)%x,kind=8)                     
          rY = real(nodeConf(iAtm)%y,kind=8)
          rZ = real(nodeConf(iAtm)%z,kind=8)

         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then 
          dX = wX - rX
          dY = wY - rY
          dZ = wZ - rZ

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          rIJ = 1000000.0
         
          !rIJ = dX**2 + dY**2 + dZ**2
          
          !rIJ = 1000000.0
          lowDX = dX
          lowDY = dY              
          lowDZ = dZ   
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
              lowDX = dX
              lowDY = dY              
              lowDZ = dZ              
            endif

            enddo
           enddo
          enddo
          
           dX = lowDX
           dY = lowDY
           dZ = lowDZ
           
         else ! its an orthorhombic box
         
          xInt = ANINT((wX-rX)/bX)
          yInt = ANINT((wY-rY)/bY)
          zInt = ANINT((wZ-rZ)/bZ)

          dX = (wX-rX) - xInt*bX
          dY = (wY-rY) - yInt*bY
          dZ = (wZ-rZ) - zInt*bZ         
          rIJ = sqrt(dX*dX+dY*dY+dZ*dZ)

         endif ! end check for orthorhombicity
         
         !rinv = 1.d0/sqrt(rIJ)
         rinv = 1.d0/rIJ
         rinv3 = rinv**3

         r3dq = dQ*rinv3

         efX = dX*r3dq
         efY = dY*r3dq
         efZ = dZ*r3dq

         efX0 = efX + efX0
         efY0 = efY + efY0
         efZ0 = efZ + efZ0
        
         enddo ! end loop over all vertical energy charges

         efMag2 = (efX0*efX0) + (efY0*efY0) + (efZ0*efZ0)
         efMag = sqrt(efMag2)

         efield_dV = efield_dV + efMag
          enddo ! End loop over 1st shell waters         

         efield_dV = efield_dV/real(nShell)  

        endif ! End check for CUDA
 
       
          
         call MPI_GATHER(real(efield_dV,kind=4),1,MPI_REAL,
     x                     efield_dVTraj,1,MPI_REAL,
     x                     master,world,ierr) 
     
         call MPI_GATHER(nShell,1,MPI_INTEGER,
     x                   nShellTraj,1,MPI_INTEGER,
     x                   master,world,ierr) 

     
          if(rank .eq. master)  then           
            do iWrite = 1, nProc


               write(iEfShel, 4657) iConf+iWrite-1,                
     x         efield_dVTraj(iWrite), nShellTraj(iWrite) 
            enddo
           endif ! End check for master       

!801             format (5x,i3,5x,i5,11x,f6.3)           
c8888             format (i6,2x,4(f15.12,2x),2(i4,2x))
4657             format (i6,2x,f15.12,2x,i6)                                        
        return
      end subroutine calcEfieldShel

      subroutine calcEFSHist()
c    *******************************************************************
c    ** calculates a histogram of the electric field of all atoms     **
c    ** within a given shell along a given O-H bond of a water        **
c    ** within that shell.                                            **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: iw, iw2, kk, jj, ii, iWrite
        integer :: nShell, niShell
        integer :: n
        integer,parameter :: maxShellN=30000
        integer :: lShell(maxShellN)      ! list of shell waters
        integer :: liShell(maxShellN/10)  ! list of shell ions
        integer :: ix, iy, iz        
        integer :: xInt, yInt, zInt     
        integer :: iii, ibin
        real    :: rIJ, watSq, rinv, rinv3
        real    :: mMag, qS(3), rtmp, r3dq, dQ
        real    :: eOH1, eOH2
        real    :: eHigh1, eLow1, eHigh2, eLow2
        logical :: firstShell, doAtom, keepFinal
        type(basic_atom)  :: initWater(3), finalWater(3), dWat
        type(basic_atom)  :: nodeShellM, masterShellM(nProc)
        type(basic_atom)  :: diH1, diH2, probeAtom
        type(basic_atom)  :: dfH1, dfH2, uOH1, uOH2
        type(basic_atom)  :: eH1, eH2, eiH1, eiH2  
        type(basic_atom)  :: probeH1, probeH2
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         do ii=1,maxBinEFS
           histEFS_i(ii) = 0.d0
         enddo

c         if ( rank .eq. master) then          
c          open(unit=iShellM,file="efs.dat",status="unknown")               
c         endif
        endif 
        
         nShell = 0        
         niShell = 0
         if(useCuda) then
#ifdef ENABLE_CUDA
            call get_shell_gpu(bX,bY,bZ,
     x                     nodeConf,qrms,
     x                    rBindCut,lShell,nShell,.false.)

#endif
         else


            !do iw = iWatStart, natom, 3    ! Loop over all water-oxygen atoms
            do iw = iWatStart, iWatEnd, 3    ! Loop over all water-oxygen atoms                        
             ! Find the 1st solvation shell waters
             firstShell = isFirstShell(nodeConf(iw)%x,
     x                              nodeConf(iw)%y,
     x                              nodeConf(iw)%z)
     
             ! Store the index of the oxygen of the 1st shell
             if(firstShell) then
              nShell = nShell + 1
              lShell(nShell) = iw         
             endif
          
           enddo ! End the loop over all waters     
       endif  ! End check for useCuda           

       do ii = iIonStart, iIonEnd    ! Loop over all ions
                        
             ! Find the 1st solvation shell ion
             firstShell = isFirstShell(nodeConf(ii)%x,
     x                              nodeConf(ii)%y,
     x                              nodeConf(ii)%z)
     
             ! Store the index of the oxygen of the 1st shell
             if(firstShell) then
              niShell = niShell + 1
              liShell(nShell) = iw         
             endif
          
        enddo ! End the loop over all waters 

          do ii = 1, nShell ! Loop over all 1st shell waters          
           ! init each hydros efield
           eH1%x = 0.d0
           eH1%y = 0.d0
           eH1%z = 0.d0

           eH2%x = 0.d0
           eH2%y = 0.d0
           eH2%z = 0.d0

           iw = lShell(ii) ! Get shell water index
           probeH1 = nodeConf(iw+1)
           probeH2 = nodeConf(iw+2)
           
!         --- here we need the O as the origin
!             of the rOH vector
           uOH1 = probeH1 - nodeConf(iw)
           uOH2 = probeH2 - nodeConf(iw)


           rTmp = uOH1 * uOH1
           rTmp = sqrt(rTmp)
           uOH1%x = uOH1%x/rTmp
           uOH1%y = uOH1%y/rTmp
           uOH1%z = uOH1%z/rTmp

           rTmp = uOH2 * uOH2
           rTmp = sqrt(rTmp)
           uOH2%x = uOH2%x/rTmp
           uOH2%y = uOH2%y/rTmp
           uOH2%z = uOH2%z/rTmp

           ! loop over all protein atoms
             do jj = iProtStart, iProtEnd
           ! get the imaged coordinates
               dQ = qrms(jj)%q
               probeAtom = nodeConf(jj)
               diH1 = probeAtom - probeH1
               diH2 = probeAtom - probeH2

           call imageCoords(diH1%x,diH1%y,diH1%z,
     x                      dfH1%x,dfH1%y,dfH1%z)

           call imageCoords(diH2%x,diH2%y,diH2%z,
     x                      dfH2%x,dfH2%y,dfH2%z)

           ! calculate the field from the probe H1
           rTmp = dfH1 * dfH1
           rTmp = sqrt(rTmp)
           rinv = 1.d0/rTmp
           rinv3 = rinv**3
           r3dq = dQ*rinv3

           eiH1%x = dfH1%x*r3dq
           eiH1%y = dfH1%y*r3dq
           eiH1%z = dfH1%z*r3dq

           ! calculate the field from the probe H2
           rTmp = dfH2 * dfH2
           rTmp = sqrt(rTmp)
           rinv = 1.d0/rTmp
           rinv3 = rinv**3
           r3dq = dQ*rinv3
           eiH2%x = dfH2%x*r3dq
           eiH2%y = dfH2%y*r3dq
           eiH2%z = dfH2%z*r3dq

            ! acucumulate the fields
            eH1 = eH1 + eiH1
            eH2 = eH2 + eiH2

             enddo
           ! end loop over all protein atoms                    

           ! loop over all ions
             do jj = iIonStart, iIonEnd
           ! get the imaged coordinates
               dQ = qrms(jj)%q
               probeAtom = nodeConf(jj)
               diH1 = probeAtom - probeH1
               diH2 = probeAtom - probeH2

           call imageCoords(diH1%x,diH1%y,diH1%z,
     x                      dfH1%x,dfH1%y,dfH1%z)

           call imageCoords(diH2%x,diH2%y,diH2%z,
     x                      dfH2%x,dfH2%y,dfH2%z)

           ! calculate the field from the probe H1
           rTmp = dfH1 * dfH1
           rTmp = sqrt(rTmp)
           rinv = 1.d0/rTmp
           rinv3 = rinv**3
           r3dq = dQ*rinv3
           eiH1%x = dfH1%x*r3dq
           eiH1%y = dfH1%y*r3dq
           eiH1%z = dfH1%z*r3dq

           ! calculate the field from the probe H2
           rTmp = dfH2 * dfH2
           rTmp = sqrt(rTmp)
           rinv = 1.d0/rTmp
           rinv3 = rinv**3
           r3dq = dQ*rinv3
           eiH2%x = dfH2%x*r3dq
           eiH2%y = dfH2%y*r3dq
           eiH2%z = dfH2%z*r3dq

            ! acucumulate the fields
            eH1 = eH1 + eiH1
            eH2 = eH2 + eiH2

             enddo
           ! end loop over all ions

           ! loop over all other atoms
             do jj = iOthStart, iOthEnd
           ! get the imaged coordinates
               dQ = qrms(jj)%q
               probeAtom = nodeConf(jj)
               diH1 = probeAtom - probeH1
               diH2 = probeAtom - probeH2

           call imageCoords(diH1%x,diH1%y,diH1%z,
     x                      dfH1%x,dfH1%y,dfH1%z)

           call imageCoords(diH2%x,diH2%y,diH2%z,
     x                      dfH2%x,dfH2%y,dfH2%z)

           ! calculate the field from the probe H1
           rTmp = dfH1 * dfH1
           rTmp = sqrt(rTmp)
           rinv = 1.d0/rTmp
           rinv3 = rinv**3
           r3dq = dQ*rinv3
           eiH1%x = dfH1%x*r3dq
           eiH1%y = dfH1%y*r3dq
           eiH1%z = dfH1%z*r3dq

           ! calculate the field from the probe H2
           rTmp = dfH2 * dfH2
           rTmp = sqrt(rTmp)
           rinv = 1.d0/rTmp
           rinv3 = rinv**3
           r3dq = dQ*rinv3
           eiH2%x = dfH2%x*r3dq
           eiH2%y = dfH2%y*r3dq
           eiH2%z = dfH2%z*r3dq

            ! acucumulate the fields
            eH1 = eH1 + eiH1
            eH2 = eH2 + eiH2

             enddo
           ! end loop over all other atoms

           ! loop over all remaining shell waters
            do jj = 1, nShell ! Loop over all 1st shell waters
             if (jj .ne. ii ) then
              iw = lShell(jj) 

             do kk = 0, 2 ! loop over all 3 water atoms              
           ! get the imaged coordinates
               dQ = qrms(iw+kk)%q
               probeAtom = nodeConf(iw+kk)
               diH1 = probeAtom - probeH1
               diH2 = probeAtom - probeH2

           call imageCoords(diH1%x,diH1%y,diH1%z,
     x                      dfH1%x,dfH1%y,dfH1%z)

           call imageCoords(diH2%x,diH2%y,diH2%z,
     x                     dfH2%x,dfH2%y,dfH2%z)

           ! calculate the field from the probe H1
           rTmp = dfH1 * dfH1
           rTmp = sqrt(rTmp)
           rinv = 1.d0/rTmp
           rinv3 = rinv**3
           r3dq = dQ*rinv3
           eiH1%x = dfH1%x*r3dq
           eiH1%y = dfH1%y*r3dq
           eiH1%z = dfH1%z*r3dq

           ! calculate the field from the probe H2
           rTmp = dfH2 * dfH2
           rTmp = sqrt(rTmp)
           rinv = 1.d0/rTmp
           rinv3 = rinv**3
           r3dq = dQ*rinv3
           eiH2%x = dfH2%x*r3dq
           eiH2%y = dfH2%y*r3dq
           eiH2%z = dfH2%z*r3dq

            ! acucumulate the fields
            eH1 = eH1 + eiH1
            eH2 = eH2 + eiH2
             enddo ! End loop over 3 wat atoms
             endif
            enddo  ! end loop over shell waters
          

           ! Find the projections along both O-H bonds
           eOH1 = eH1 * uOH1
           eOH2 = eH2 * uOH2

c           print *, "node : ", rank 
c           print *, " delEFS = ", delEFS
c           print *, " eOH1   = ", eOH1
c           print *, " eOH2   = ", eOH2

           ! and make a histogram of the results
           ibin = int((eOH1+efsShift)/delEFS) + 1
           if(ibin .lt. maxBinEFS .and. 
     x        ibin .gt. 0) then
            histEFS_i(ibin) = histEFS_i(ibin) + 1.d0
           endif   

           ibin = int((eOH2 + efsShift)/delEFS) + 1
           if(ibin .lt. maxBinEFS .and. 
     x        ibin .gt. 0) then
            histEFS_i(ibin) = histEFS_i(ibin) + 1.d0
           endif   

          enddo ! End loop over 1st shell waters

         call MPI_BARRIER(world,ierr)
c         ! Send the total dipople moment as an atom type..yes, lazy
c         call MPI_GATHER(nodeShellM,1,ATOM_TYPE,
c     x                   masterShellM,1,ATOM_TYPE,
c     x                   master,world,ierr) 


c         call MPI_GATHER(nShell,1,MPI_INTEGER,
c     x                   nShellTraj,1,MPI_INTEGER,
c     x                   master,world,ierr) 

     
c          if(rank .eq. master)  then           
c            do iWrite = 1, nProc
            
c               mMag = masterShellM(iWrite)%x*masterShellM(iWrite)%x+
c     x                masterShellM(iWrite)%y*masterShellM(iWrite)%y+         
c     x                masterShellM(iWrite)%z*masterShellM(iWrite)%z
     
c               mMag = sqrt(mMag)
               
c               write(iShellM, 7778) iConf+iWrite-1, 
c     x         mMag, masterShellM(iWrite)%x, masterShellM(iWrite)%y, 
c     x         masterShellM(iWrite)%z, nShellTraj(iWrite)
     
c            enddo 



c           if (iConf .eq. 1 ) then
         
c            open(unit=8557,file="dipole-shell.dat",status="unknown") 
c             do jj = 1, nShell
c              write(8557, 7558) watOpos(jj)%x, watOpos(jj)%y, 
c     x         watOpos(jj)%z,
c     x         watEHat(jj)%x, watEHat(jj)%y, watEHat(jj)%z
c             enddo
c            close(8557)
c           endif  ! end  check for 1st conf

c           endif ! End check for master       
          
c7778             format (i6,2x,4(E14.6,2x),i6)
c7558             format (6(E14.6,2x))                                        
        return
      end subroutine calcEFSHist

      subroutine calcCumuBind()
c    *******************************************************************
c    ** calculates cumulants of 1st shell energies (total and         **
c    **  vertical) at every timestep.                                 **
c    **                                                               **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: iw, iw2, irc, jj, iWrite
        integer :: iAtm, irc2
        integer :: nShell, ibin
        integer :: lShell(25000)  ! list of first atoms of water
        integer :: ix, iy, iz        
        integer :: xInt, yInt, zInt        
        real    :: ddX, ddY, ddZ, ddQ
        real*8 :: prX, prY, prZ, prQ, prR
        real    :: dX, dY, dZ, myQ         
        real*8  :: dQ, rX, rY, rZ, rIJ
        real*8  :: pQ, pX, pY, pZ
        real*8  :: wQ, wX, wY, wZ        
        real*8  :: dfrx, dfry, dfrz, delr2   
        real*8  :: eShift, rNShell
        logical :: firstShell, doAtom
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then          
          open(unit=iCBind,file="enerCumuBind.dat",status="unknown")               
         endif
        endif 
        
         nShell = 0        
         
         if(useCuda) then
#ifdef ENABLE_CUDA 
            call get_shell_gpu(bX,bY,bZ,
     x                     nodeConf,qrms,
     x                    lShell,nShell,.false.)
#endif
         else

         do iw = iWatStart, natom, 3    ! Loop over all water-oxygen atoms
               
          wX = real(nodeConf(iw)%x,kind=8)                     
          wY = real(nodeConf(iw)%y,kind=8)
          wZ = real(nodeConf(iw)%z,kind=8)          
         
          ! Find the 1st solvation shell waters
          firstShell = isFirstShell(nodeConf(iw)%x,
     x                              nodeConf(iw)%y,
     x                              nodeConf(iw)%z)
     
          ! Store the index of the oxygen of the 1st shell
          if(firstShell) then
           nShell = nShell + 1
           lShell(nShell) = iw         
          endif
          
        enddo ! End the loop over all waters     
        endif ! CHeck for CUDA

         dv0sh_i = 0.d0
         dv0sh   = 0.d0
         dv0shde_i = 0.d0
         dv0shde   = 0.d0    
         dv0bsh  = 0.d0
         dv0bshde = 0.d0
         
         ! Set all node's moments to 0
!         mdepsdDeps = 0.d0
!         mdeps2dDeps = 0.d0
!         mdeps3dDeps = 0.d0
!         mdDeps2 = 0.d0
!         mdDeps2deps = 0.d0
!         mdDeps2deps2 = 0.d0
          lamYi = 0.d0
          Yi = 0.d0
         rNShell = real(nShell,kind=8)
         
          do jj = 1, nShell ! Loop over all 1st shell waters          
          iw = lShell(jj)
          
          dv0bsh_i  = 0.d0          
          dv0bshde_i  = 0.d0              
           do iw2 = iw, iw+2 ! Loop over all 3 atoms
            wX = real(nodeConf(iw2)%x,kind=8)
            wY = real(nodeConf(iw2)%y,kind=8)
            wZ = real(nodeConf(iw2)%z,kind=8)
            wQ = real(qrms(iw2)%q,kind=8)
            
           do irc = 1, iProtEnd  ! Loop over protein atoms
            prX = real(nodeConf(irc)%x,kind=8)                     
            prY = real(nodeConf(irc)%y,kind=8)
            prZ = real(nodeConf(irc)%z,kind=8)
            prQ = real(qrms(irc)%q,kind=8)            

            !*******************         
         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then 
          dX = wX - prX
          dY = wY - prY
          dZ = wZ - prZ

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          rIJ = dX**2 + dY**2 + dZ**2
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
            endif

            enddo
           enddo
          enddo
          rIJ = sqrt(rIJ)
         else ! its an orthorhombic box
         
          xInt = ANINT((wX-prX)/bX)
          yInt = ANINT((wY-prY)/bY)
          zInt = ANINT((wZ-prZ)/bZ)

          dX = (wX-prX) - xInt*bX
          dY = (wY-prY) - yInt*bY
          dZ = (wZ-prZ) - zInt*bZ         
          rIJ = sqrt(dX*dX+dY*dY+dZ*dZ)

         endif ! End check for truncation            
            !*******************
           dv0bsh_i = dv0bsh_i + wQ*prQ/rIJ 
           
           ! Make sure you do a separate deltaE calc. for rc atoms
            call searchDBAtom(irc,myQ,doAtom)
             if (doAtom) then

             prQ = real(myQ,kind=8)
             dv0bshde_i = dv0bshde_i + wQ*prQ/rIJ                          
            endif

           
            enddo ! End loop over protein atoms   
      
           enddo ! End loop over all 3 water atoms
          
           ! Calculate the instantaneous fluctuations
           deps = dv0bsh_i - r8epsAvg
           dDeps = dv0bshde_i - r8epsDAvg
           
!           deps2 = deps*deps  ! 2nd moments
           dDeps2 = dDeps*dDeps
           

           ! Calculate exponent and reor., and average         
           expBde = -betaT*abs(deps)
           expBde = exp(expBde)
           lamYi = lamYi + dDeps2*expBde
           
           Yi = Yi + dDeps*(expBde - 1)

!           deps3 = deps2*deps  ! 3nd moments
           
           ! Calculate the moments on the nodes
!           mdepsdDeps = mdepsdDeps + deps*dDeps
!           mdeps2dDeps = mdeps2dDeps + deps2*dDeps
!           mdeps3dDeps = mdeps3dDeps + deps3*dDeps
           
!           mdDeps2 = mdDeps2 + dDeps2
!           mdDeps2deps = mdDeps2deps + dDeps2*deps
!           mdDeps2deps2 = mdDeps2deps2 + dDeps2*deps2                   
        
          enddo ! End loop over 1st shell waters
            

        ! Average the moments on each node according to the number of waters
            lamYi = lamYi/rNShell
            Yi = Yi/rNShell + r8epsDAvg
!           mdepsdDeps = mdepsdDeps/rNShell
!           mdeps2dDeps = mdeps2dDeps/rNShell
!           mdeps3dDeps = mdeps3dDeps/rNShell
           
!           mdDeps2 = mdDeps2/rNShell
!           mdDeps2deps = mdDeps2deps/rNShell
!           mdDeps2deps2 = mdDeps2deps2/rNShell
        
         call MPI_GATHER(real(lamYi,kind=4),1,MPI_REAL,
     x                   lamYiTraj,1,MPI_REAL,
     x                   master,world,ierr) 

         call MPI_GATHER(real(Yi,kind=4),1,MPI_REAL,
     x                   YiTraj,1,MPI_REAL,
     x                   master,world,ierr) 
     
!         call MPI_GATHER(real(mdepsdDeps,kind=4),1,MPI_REAL,
!     x                   depsdDepsTraj,1,MPI_REAL,
!     x                   master,world,ierr) 

!         call MPI_GATHER(real(mdeps2dDeps,kind=4),1,MPI_REAL,
!     x                   deps2dDepsTraj,1,MPI_REAL,
!     x                   master,world,ierr) 

!         call MPI_GATHER(real(mdeps3dDeps,kind=4),1,MPI_REAL,
!     x                   deps3dDepsTraj,1,MPI_REAL,
!     x                   master,world,ierr) 

!         call MPI_GATHER(real(mdDeps2,kind=4),1,MPI_REAL,
!     x                   dDeps2Traj,1,MPI_REAL,
!     x                   master,world,ierr) 
     
!         call MPI_GATHER(real(mdDeps2deps,kind=4),1,MPI_REAL,
!     x                   dDeps2depsTraj,1,MPI_REAL,
!     x                   master,world,ierr) 
     
!         call MPI_GATHER(real(mdDeps2deps2,kind=4),1,MPI_REAL,
!     x                   dDeps2deps2Traj,1,MPI_REAL,
!     x                   master,world,ierr) 

         call MPI_GATHER(nShell,1,MPI_INTEGER,
     x                   nShellTraj,1,MPI_INTEGER,
     x                   master,world,ierr) 

     
          if(rank .eq. master)  then           
            do iWrite = 1, nProc
               write(iCBind, 8889) iConf+iWrite-1, 
     x          lamYiTraj(iWrite), YiTraj(iWrite),
     x          nShellTraj(iWrite)               
!     x         depsdDepsTraj(iWrite), deps2dDepsTraj(iWrite), 
!     x         deps3dDepsTraj(iWrite), dDeps2Traj(iWrite), 
!     x         dDeps2depsTraj(iWrite), dDeps2deps2Traj(iWrite),      
            enddo
           endif ! End check for master       
          
8889             format (i6,2x,2(E14.6,2x),i4)                    
        return
      end subroutine calcCumuBind      

      subroutine collectPotL()
c    *******************************************************************
c    ** sums POTL data onto master and writes to a file               **
c    **                                                               **
c    **                                                               **
c    *******************************************************************            
        implicit none

        ! Locals
        integer :: iPts

c    *******************************************************************        
       if(rank .eq. master) then
        allocate(rmin(nPts))
        allocate(potl(nPts))
        allocate(potl2(nPts))
        allocate(distPL(nPts))

        do iPts = 1, nPts
          rmin(iPts) = 0.d0
          potl(iPts) = 0.d0
          potl2(iPts) = 0.d0
          distPL(iPts) = 0.d0
        enddo
        
       endif
       
          call MPI_REDUCE(arrRMIN,rmin,nPts,
     x          MPI_REAL,MPI_SUM,master,world,ierr)
     
          call MPI_REDUCE(arrPotL,potl,nPts,
     x          MPI_REAL,MPI_SUM,master,world,ierr)

          call MPI_REDUCE(arrPotL2,potl2,nPts,
     x          MPI_REAL,MPI_SUM,master,world,ierr)

          call MPI_REDUCE(arrDist,distPL,nPts,
     x          MPI_REAL,MPI_SUM,master,world,ierr)


       if(rank .eq. master) then
         open(unit=iPotL,file="potL.dat",status="unknown")

         do iPts = 1, nPts
           rmin(iPts) = rmin(iPts)/real(nConf)
           distPL(iPts) = distPL(iPts)/real(nConf)
           potl(iPts) = potl(iPts)/real(nConf)
           potl2(iPts) = potl2(iPts)/real(nConf)
           
           write(iPotL,*) distPL(iPts), 
     x           rmin(iPts), potl(iPts), potl2(iPts)
         enddo
         
         close(iPotL)

         deallocate(rmin)
         deallocate(potl)
         deallocate(potl2)
         deallocate(distPL)

         endif
         
         call MPI_BARRIER(world,ierr)
         
         deallocate(arrRMIN)
         deallocate(arrPotL)
         deallocate(arrPotL2)
         deallocate(arrDist)
         deallocate(arrPts)

        !!!! NORM: arrDist, arrPotL, arrPotL2, arrRMIN
        return
      end subroutine collectPotL
      
      subroutine calcPotL()
c    *******************************************************************
c    ** calculates the moments of the potential on a line from the    **
c    ** center of mass of the protein to the atom at index iAtomPoint **
c    **                                                               **
c    **   NOTE: This routines assumes the trajectory has been         **
c    **         reimaged so the origin is at the COM of the protein   **
c    **                                                               **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: iw, iPts
        real :: rXpt, rYpt, rZpt, dPt, lam, rIJ
        real :: rMin, sumPot, qw, sumPot2, pot_i
        type(basic_atom)  :: rWat, rDiff, rPt
c    ******************************************************************* 
        if (iConf .eq. 1 ) then  ! Alloc and Normalize

        call MPI_BCAST(nPts,1 
     x               ,MPI_INTEGER,master, world, ierr)      

        call MPI_BCAST(iAtomPoint,1 
     x               ,MPI_INTEGER,master, world, ierr)      
        
        allocate(arrRMIN(nPts))
        allocate(arrPotL(nPts))
        allocate(arrPotL2(nPts))
        allocate(arrDist(nPts))
        allocate(arrPts(nPts))
        
        do iPts = 1, nPts
          arrRMIN(iPts) = 0.d0
          arrPotL(iPts) = 0.d0
          arrPotL2(iPts) = 0.d0
          arrDist(iPts) = 0.d0
        enddo

        
        endif 
                  
        rPt = nodeConf(iAtomPoint)  ! Get position
        rXpt = rPt%x
        rYpt = rPt%y
        rZpt = rPt%z
        dPt = 1./real(nPts)  !!! MUST SEND nPts, iAtomPoint        
        lam = -dPt/2.
        
        do iPts =1, nPts
        lam = lam + dPt
        arrPts(iPts)%x = rXpt*lam
        arrPts(iPts)%y = rYpt*lam
        arrPts(iPts)%z = rZpt*lam  
        
        rDiff = arrPts(iPts)
        rIJ = rDiff*rDiff
        rIJ = Sqrt(rIJ)
        !call imageCoords(xxI,yyI,zzI,xxF,yyF,zzF)
        arrDist(iPts) = arrDist(iPts) + rIJ
          
        !if(rank .eq. master) then
        !  print *, "rIJ = ", rIJ, "lam = ", lam
        !endif
        enddo
        
        do iPts =1, nPts
        rMin = 99999999.99
        sumPot = 0.d0
        !sumPot2 = 0.d0
        rPt = arrPts(iPts)
         do iw = iWatStart, natom  ! Loop over all waters and sum up the potential
           rWat = nodeConf(iw) ! Store the initial coordinates of h2o atoms
           qw = qrms(iw)%q          ! Store the charge of the h2o atoms
         
           rDiff = rWat - rPt
           
           rIJ = rDiff*rDiff
           rIJ = Sqrt(rIJ)
           
           pot_i = qw/rIJ
           sumPot = sumPot + pot_i
           
           
           if(rIJ .lt. rMin) then
            rMin = rIJ
           endif
           
         enddo ! End loop over all water atoms
         
         sumPot2 = sumPot*sumPot
         
         arrRMIN(iPts) = arrRMIN(iPts) + rMin
         arrPotL(iPts) = arrPotL(iPts) + sumPot
         arrPotL2(iPts) = arrPotL2(iPts) + sumPot2
         
        enddo  ! End loop over all potential line points
        
        return
      end subroutine calcPotL  

      subroutine calcShellM()
c    *******************************************************************
c    ** calculates the total dipole moment of 1st shell waters        **
c    **                                                               **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: iw, iw2, jj, iWrite
        integer :: nShell
        integer,parameter :: maxShellN=30000
        integer :: lShell(maxShellN)  ! list of first atoms of water
        integer :: ix, iy, iz        
        integer :: xInt, yInt, zInt     
        integer :: iii
        real    :: rIJ, watSq
        real    :: mMag, qS(3)
        logical :: firstShell, doAtom, keepFinal
        type(basic_atom)  :: initWater(3), finalWater(3), dWat
        type(basic_atom)  :: nodeShellM, masterShellM(nProc)
        type(basic_atom)  :: watOpos(maxShellN)
        type(basic_atom)  :: watEHat(maxShellN)
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then          
          open(unit=iShellM,file="shellM.dat",status="unknown")               
         endif
        endif 
        
         nShell = 0        

         if(useCuda) then
#ifdef ENABLE_CUDA
            call get_shell_gpu(bX,bY,bZ,
     x                     nodeConf,qrms,
     x                    rBindCut,lShell,nShell,.false.)
#endif
         else

            do iw = iWatStart, natom, 3    ! Loop over all water-oxygen atoms
                        
             ! Find the 1st solvation shell waters
             firstShell = isFirstShell(nodeConf(iw)%x,
     x                              nodeConf(iw)%y,
     x                              nodeConf(iw)%z)
     
             ! Store the index of the oxygen of the 1st shell
             if(firstShell) then
              nShell = nShell + 1
              lShell(nShell) = iw         
             endif
          
           enddo ! End the loop over all waters     
       endif  ! End check for useCuda



           ! Initialize dipole moments
           nodeShellM%x = 0.d0
           nodeShellM%y = 0.d0
           nodeShellM%z = 0.d0

          do jj = 1, nShell ! Loop over all 1st shell waters          
           iw = lShell(jj) ! Get shell water index
                    
           iii = 0
           do iw2 = iw, iw+2 ! Loop over all 3 atoms
           
           iii = iii + 1           
           initWater(iii) = nodeConf(iw2) ! Store the initial coordinates of h2o atoms
           qS(iii) = qrms(iw2)%q          ! Store the charge of the h2o atoms
           
           call imageCoords(initWater(iii)%x,
     x          initWater(iii)%y,initWater(iii)%z,
     x          finalWater(iii)%x,finalWater(iii)%y,
     x     finalWater(iii)%z)
                                     
      
           enddo ! End loop over all 3 water atoms

           keepFinal = .true.
           do iii = 2, 3 ! Loop over both hydrogens
           
            dWat = finalWater(1) - finalWater(iii)                                          
            rIJ = sqrt(dWat%x*dWat%x + dWat%y*dWat%y + dWat%z*dWat%z)
 
            if(rIJ .lt. 2.5d0) keepFinal = .false.
           
           enddo ! End loop over both hydros

           watEHat(jj)%x = 0.d0
           watEHat(jj)%y = 0.d0
           watEHat(jj)%z = 0.d0
          ! Calculate dipole moments
           do iii = 1, 3
            if(keepFinal) then

              nodeShellM%x = nodeShellM%x + qS(iii)*finalWater(iii)%x
              nodeShellM%y = nodeShellM%y + qS(iii)*finalWater(iii)%y
              nodeShellM%z = nodeShellM%z + qS(iii)*finalWater(iii)%z   
              if (iConf .eq. 1 ) then
               if (iii .eq. 1) then
                watOpos(jj) = finalWater(iii)
               endif
                watEHat(jj)%x = watEHat(jj)%x + 
     x                          qS(iii)*finalWater(iii)%x
                watEHat(jj)%y = watEHat(jj)%y + 
     x                           qS(iii)*finalWater(iii)%y
                watEHat(jj)%z = watEHat(jj)%z + 
     x                           qS(iii)*finalWater(iii)%z
              endif            
              
            else
           
              nodeShellM%x = nodeShellM%x + qS(iii)*initWater(iii)%x
              nodeShellM%y = nodeShellM%y + qS(iii)*initWater(iii)%y
              nodeShellM%z = nodeShellM%z + qS(iii)*initWater(iii)%z              
              if (iConf .eq. 1 ) then
               if ( iii .eq. 1 ) then
                watOpos(jj) = initWater(iii)
               endif
                watEHat(jj)%x = watEHat(jj)%x + 
     x                           qS(iii)*initWater(iii)%x
                watEHat(jj)%y = watEHat(jj)%y + 
     x                           qS(iii)*initWater(iii)%y
                watEHat(jj)%z = watEHat(jj)%z + 
     x                           qS(iii)*initWater(iii)%z
              endif              
            endif
          
           enddo ! End loop of all 3 water atoms
           watSq = watEHat(jj)*watEHat(jj)
           watSq = sqrt(watSq)
           watEHat(jj)%x = watEHat(jj)%x/watSq
           watEHat(jj)%y = watEHat(jj)%y/watSq
           watEHat(jj)%z = watEHat(jj)%z/watSq
          enddo ! End loop over 1st shell waters

         extraShellM = nodeShellM  ! keep the current shell dipole for further calcs


         call MPI_BARRIER(world,ierr)
         ! Send the total dipople moment as an atom type..yes, lazy

         call MPI_GATHER(nodeShellM,1,ATOM_TYPE,
     x                   masterShellM,1,ATOM_TYPE,
     x                   master,world,ierr) 

         call MPI_GATHER(nShell,1,MPI_INTEGER,
     x                   nShellTraj,1,MPI_INTEGER,
     x                   master,world,ierr) 

          if(rank .eq. master)  then           
            do iWrite = 1, nProc
            
             mMag = masterShellM(iWrite)%x*masterShellM(iWrite)%x+
     x              masterShellM(iWrite)%y*masterShellM(iWrite)%y+         
     x              masterShellM(iWrite)%z*masterShellM(iWrite)%z
     
               mMag = sqrt(mMag)
               
               write(iShellM, 7778) iConf+iWrite-1, 
     x       mMag, masterShellM(iWrite)%x, masterShellM(iWrite)%y, 
     x       masterShellM(iWrite)%z, nShellTraj(iWrite)
     
            enddo 

c           if (iConf .eq. 1 ) then
         
c            open(unit=8557,file="dipole-shell.dat",status="unknown") 
c             do jj = 1, nShell
c              write(8557, 7558) watOpos(jj)%x, watOpos(jj)%y, 
c     x         watOpos(jj)%z,
c     x         watEHat(jj)%x, watEHat(jj)%y, watEHat(jj)%z
c             enddo
c            close(8557)
c           endif  ! end  check for 1st conf

           endif ! End check for master       

         call MPI_BARRIER(world,ierr)          
7778             format (i6,2x,4(E14.6,2x),i6)
7558             format (6(E14.6,2x))                                        
        return
      end subroutine calcShellM  

      subroutine calcProtM()
c    *******************************************************************
c    ** calculates Mu, dipole vector of the protein with respect      **
c    ** to the center of charge.                                      **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: ipr, iWrite
        integer :: ii,jj,iw2, iw, iCosBin
        integer :: massOrQ
        integer :: nShell
        integer, parameter :: nMaxShell = 30000
        integer :: lShell(nMaxShell), iBins(nMaxShell) ! list of first and bins atoms of water
        real    :: totM, totQ, prM, prQ, prMuLen, prMuLen2
        real    :: wwX, wwY, wwZ, wwQ
        real    :: eLen2, eLen
        logical :: firstShell, doOutput
        type(basic_atom)  :: prMuCntr, prAtom,eHat
        type(basic_atom)  :: prMuNode, watMu, prMuHat
        type(basic_atom)  :: initWater(3), finalWater(3), dWat
        type(basic_atom)  :: nodePmu, masterPmu(nProc)  !
        type(basic_atom)  :: tmpAtom  ! for temp storage
        real  :: cosProtWatNode, cosProtWat
        real  :: cosProtWatTraj(nProc)
c    ******************************************************************* 
         
        if (iConf .eq. 1 ) then
         histCos_i(:) = 0.0d0
         if ( rank .eq. master) then
          open(unit=iPMu,file="proteinMu.dat",status="unknown")               
         endif
        endif 

         if(useCuda) then
#ifdef ENABLE_CUDA 
            call get_shell_gpu(bX,bY,bZ,
     x                         nodeConf,qrms,
     x                         rBindCut,lShell,nShell,.false.)


           ! Calculate the protein's diople
           prAtom%x = 0.d0
           prAtom%y = 0.d0
           prAtom%z = 0.d0
           massOrQ = 0 !set to 0 to calculate the center of mass vector
           doOutput = .false. ! we dont want all the positions copied back

            call get_prmu_gpu(qrms,
     x                         prAtom%x,prAtom%y,prAtom%z, 
     x                         prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x                         totM,massOrQ,doOutput,tmpAtom)

            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM

            massOrQ = 1 !set to 1 to calculate the dipole moment

           ! now calculate the protein diople with respect to the center of mass
            call get_prmu_gpu(qrms,
     x                         prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x               prMuHat%x,prMuHat%y,prMuHat%z,totQ,
     x                      massOrQ,doOutput,tmpAtom)

           ! save the protein dipole moment, before its normalized
           nodePmu = prMuHat
           extraPmu = prMuHat
           !prMuNode = prMuCntr

           ! normalize the protein vector
           eLen2 = prMuHat*prMuHat
           eLen = sqrt(eLen2)
           prMuHat%x = prMuHat%x/eLen
           prMuHat%y = prMuHat%y/eLen
           prMuHat%z = prMuHat%z/eLen

           ! now calculate the array of water-dipole/protein-dipole cosigns
           call get_shl_cos_gpu(nodeConf,qrms,
     x                          prMuHat%x,prMuHat%y,prMuHat%z, 
     x                          dCosT,lShell,nShell,
     x                          cosProtWatNode,iBins)

           ! Now complete the partial histogram !!! SHOULD MOVE TO THE GPU!!!
            do i=1, nShell
              iCosBin = iBins(i)
              if(iCosBin .lt. maxBinCos) then
              histCos_i(iCosBin) = histCos_i(iCosBin) + 1.d0
              endif  

              
            enddo ! end loop over 1st shell waters

#endif /* CUDA CODE */
         else ! we are only doing non-CUDA, CPU code

           prMuCntr%x = 0.d0
           prMuCntr%y = 0.d0
           prMuCntr%z = 0.d0
           totM = 0.d0           
           prMuNode%x = 0.d0
           prMuNode%y = 0.d0
           prMuNode%z = 0.d0           
           
           do ipr = 1, iProtEnd  ! Loop over protein atoms to find the center of charge

            prAtom = nodeConf(ipr)
            
c           call imageCoords(prImag%x,prImag%y,prImag%z,
c     x                      prAtom%x,prAtom%y,prAtom%z)
            prM = qrms(ipr)%m
            
            totM = totM + prM           

            prMuCntr%x = prMuCntr%x + prM*prAtom%x
            prMuCntr%y = prMuCntr%y + prM*prAtom%y
            prMuCntr%z = prMuCntr%z + prM*prAtom%z
            
           enddo ! End loop over all protein atoms 
            
            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM

           
           do ipr = 1, iProtEnd  ! Loop over protein atoms to find the protein dipole

            prAtom = nodeConf(ipr) ! set the position
            !prAtom = nodeConf(ipr)

            prQ = qrms(ipr)%q
            
            prMuNode%x = prMuNode%x + prQ*(prAtom%x - prMuCntr%x)
            prMuNode%y = prMuNode%y + prQ*(prAtom%y - prMuCntr%y)
            prMuNode%z = prMuNode%z + prQ*(prAtom%z - prMuCntr%z)
            
           enddo ! End loop over all protein atoms 

           ! Find the unit vector of the protein dipole moment
           prMuLen2 = prMuNode*prMuNode
           prMuLen = sqrt(prMuLen2)

           extraPmu = prMuNode  ! Keep the current protein dipole for further calcs
          
           ! save the protein dipole moment, before its normalized
           nodePmu = prMuNode

           prMuHat%x = prMuNode%x/prMuLen
           prMuHat%y = prMuNode%y/prMuLen
           prMuHat%z = prMuNode%z/prMuLen


c          print *, "about to find the 1st shell atoms..."          

          nShell = 0
          do iw = iWatStart, natom, 3    ! Loop over all water-oxygen atoms
                        
           ! Find the 1st solvation shell waters
           firstShell = isFirstShell(nodeConf(iw)%x,
     x                              nodeConf(iw)%y,
     x                              nodeConf(iw)%z)
     
           ! Store the index of the oxygen of the 1st shell
           if(firstShell) then
            nShell = nShell + 1
            lShell(nShell) = iw         
           endif
          
          enddo ! End the loop over all waters     

c          print *, "about to loop over the 1st shell atoms..."          

           cosProtWatNode = 0.d0
          do jj = 1, nShell ! Loop over all 1st shell waters          
           iw = lShell(jj) ! Get shell water index
                    
           eHat%x = 0.d0
           eHat%y = 0.d0
           eHat%z = 0.d0 
           
           do iw2 = iw, iw+2 ! Loop over all 3 atoms           
           
           wwX = nodeConf(iw2)%x
           wwY = nodeConf(iw2)%y
           wwZ = nodeConf(iw2)%z
           wwQ = qrms(iw2)%q
           
c            call imageCoords(wwX,wwY,wwZ,
c     x                       dXn, dYn, dZn)
c
           watMu%x = wwQ*wwX
           watMu%y = wwQ*wwY
           watMu%z = wwQ*wwZ
c           watMu%x = wwQ*dXn
c           watMu%y = wwQ*dXn
c           watMu%z = wwQ*dZn

           eHat = eHat + watMu
c           eHat%x = eHat%x + wwQ*wwX
c           eHat%y = eHat%y + wwQ*wwY
c           eHat%z = eHat%z + wwQ*wwZ           

          enddo ! loop over water atoms

           eLen2 = eHat*eHat
           eLen = sqrt(eLen2)

                    
           eHat%x = eHat%x/eLen
           eHat%y = eHat%y/eLen
           eHat%z = eHat%z/eLen
          
           cosProtWat = prMuHat*eHat

           if(cosProtWat < -1.0) cosProtWat = -1.0 
           if(cosProtWat > 1.0) cosProtWat = 1.0 

           cosProtWatNode = cosProtWatNode + cosProtWat

           ! find the bin numbers
           iCosBin = int (cosProtWat/dCosT ) + 1

           if(iCosBin .lt. maxBinCos) then
            histCos_i(iCosBin) = histCos_i(iCosBin) + 1.d0
           endif  

        enddo ! loop over 1shell waters            
        ! normalize cosProtWatNode
        cosProtWatNode = cosProtWatNode/real(nShell,kind=4)

       endif ! END CHECK FOR CUDA

         call MPI_BARRIER(world,ierr)


c          print *, "about to gather the 1st data..."          

         call MPI_GATHER(cosProtWatNode,1,MPI_REAL,
     x                   cosProtWatTraj,1,MPI_REAL,
     x                   master,world,ierr) 


         call MPI_GATHER(nShell,1,MPI_INTEGER,
     x                   nShellTraj,1,MPI_INTEGER,
     x                   master,world,ierr) 

         call MPI_GATHER(nodePmu,1,ATOM_TYPE,
     x                   masterPmu,1,ATOM_TYPE,
     x                   master,world,ierr) 

     
          if(rank .eq. master)  then           
            do iWrite = 1, nProc
                           
               write(iPMu, 7779) iConf+iWrite-1, 
     x         cosProtWatTraj(iWrite),
     x         masterPmu(iWrite)%x,
     x         masterPmu(iWrite)%y,
     x         masterPmu(iWrite)%z,
     x         nShellTraj(iWrite)
     
            enddo 
           endif ! End check for master       
          
7779             format (i6,2x,4(E14.6,2x),i6)                    
        return
      end subroutine calcProtM


      subroutine calcProtV()
c    *******************************************************************
c    ** calculates the volume of a protein using MC integration       **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: ipr, iWrite
        integer :: ii,jj,iw2, iw, idum
        integer :: iOxy, iH1, iH2
        integer :: massOrQ
        integer :: nShell
        integer, parameter :: nMaxShell = 30000
        integer :: lShell(nMaxShell), iBins(nMaxShell) ! list of first and bins atoms of water
        real    :: totM, totQ, prM, prQ, prMuLen, prMuLen2
        real    :: prR, prR2
        real    :: dxyz2, x1, y1, z1
        real    :: dX, dY, dZ
        real    :: sum_test, sum_prot_vol, sum_wat_vol
        real    :: prot_e, prot_ee
        real    :: wat_e, wat_ee
        real, parameter :: pi = 3.14159d0
        real    :: testR2, sphereVol
        real    :: oxyR, oxyR2
        real    :: h1R, h1R2
        real    :: h2R, h2R2
        real    :: h1Val,h2Val,oVal,wVal
        logical :: firstShell, doOutput
        logical :: insideOxy, insideH1, insideH2
        logical :: outsideProt
        logical :: outsideOxy, outsideH1
        type(basic_atom)  :: prMuCntr, prAtom, tmpAtom
        type(basic_atom)  :: oxyAtom, h1Atom, h2Atom
        type(basic_atom)  :: dXYZ, testXYZ  ! for temp storage
c    ******************************************************************* 
         
        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then
          open(unit=iPVol,file="proteinVol.dat",status="unknown")               
         endif
        endif 

         if(useCuda) then
#ifdef ENABLE_CUDA 
            call get_shell_gpu(bX,bY,bZ,
     x                         nodeConf,qrms,
     x                         rBindCut,lShell,nShell,.true.)


           ! Calculate the protein's diople
           prAtom%x = 0.d0
           prAtom%y = 0.d0
           prAtom%z = 0.d0
           massOrQ = 0 !set to 0 to calculate the center of mass vector
           doOutput = .false. ! we dont want all the positions copied back

            call get_prmu_gpu(qrms,
     x                         prAtom%x,prAtom%y,prAtom%z, 
     x                         prMuCntr%x,prMuCntr%y,prMuCntr%z,
     x                         totM,massOrQ,doOutput,tmpAtom)

            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM

            protRadius = 0.d0
           do ipr = iProtStart, iProtEnd  ! Loop over protein atoms to find the proteins effective radius

            prAtom = nodeConf(ipr)

            dXYZ = prAtom - prMuCntr
            dxyz2 = dXYZ * dXYZ

            if (dxyz2 > protRadius) protRadius = dxyz2
            
           enddo ! End loop over all protein atoms 

           protRadius = sqrt(protRadius)
    
           ! add an extra 3.5A to make sure all water atoms are included
           totalRadius = protRadius + rBindCut + 3.5d0

            ! generate random numbers for the volume calculation
            ! NOTE: this should move to the GPU, seriously.
           do iAttempt = 1, nAttempts
            ranX(iAttempt) = ran1(idum)
            ranY(iAttempt) = ran1(idum)
            ranZ(iAttempt) = ran1(idum)
           enddo

            call get_vol_shl_gpu(bX,bY,bZ,
     x                           nodeConf,qrms,
     x                            lShell,nShell,
     x                            totalRadius,
     x                            ranX, ranY, ranZ,            
     x                  prMuCntr%x, prMuCntr%y, prMuCntr%z,
     x                        mcVolProt, mcVolWat,
     x                        mcEpsProt, mcEpsWat)


#endif /* CUDA CODE */
         else ! we are only doing non-CUDA, CPU code

           prMuCntr%x = 0.d0
           prMuCntr%y = 0.d0
           prMuCntr%z = 0.d0
           totM = 0.d0           
           nShell = 0

           ! first, find all of the waters in the shell
           ! and this time, all atoms of the water must 
           ! be within
          do iw = iWatStart, natom,3    ! Loop over all waters

           insideOxy = isFirstShell(nodeConf(iw)%x,
     x                              nodeConf(iw)%y,
     x                              nodeConf(iw)%z)

           insideH1 = isFirstShell(nodeConf(iw+1)%x,
     x                             nodeConf(iw+1)%y,
     x                             nodeConf(iw+1)%z)

           insideH2 = isFirstShell(nodeConf(iw+2)%x,
     x                             nodeConf(iw+2)%y,
     x                             nodeConf(iw+2)%z)
     
           ! Store the index of the oxygen of the 1st shell
           if(insideOxy .and. 
     x         insideH1 .and. insideH2) then
            nShell = nShell + 1
            lShell(nShell) = iw         
           endif
          
          enddo ! End the loop over all waters     
           
           do ipr = 1, iProtEnd  ! Loop over protein atoms to find the center of mass

            prAtom = nodeConf(ipr)            
            prM = qrms(ipr)%m
            
            totM = totM + prM           

            prMuCntr%x = prMuCntr%x + prM*prAtom%x
            prMuCntr%y = prMuCntr%y + prM*prAtom%y
            prMuCntr%z = prMuCntr%z + prM*prAtom%z
            
           enddo ! End loop over all protein atoms 
            
            prMuCntr%x = prMuCntr%x/totM
            prMuCntr%y = prMuCntr%y/totM
            prMuCntr%z = prMuCntr%z/totM

            protRadius = 0.d0
           do ipr = iProtStart, iProtEnd  ! Loop over protein atoms to find the proteins effective radius

            prAtom = nodeConf(ipr)

            dXYZ = prAtom - prMuCntr
            dxyz2 = dXYZ * dXYZ

            if (dxyz2 > protRadius) protRadius = dxyz2
            
           enddo ! End loop over all protein atoms 

           protRadius = sqrt(protRadius)
    
           ! add an extra 3.5A to make sure all water atoms are included
           totalRadius = protRadius + rBindCut + 3.5d0
           totalR2 = totalRadius * totalRadius

           sphereVol = (4.0d0/3.0d0)*pi* (totalRadius**3)

           sum_prot_vol  = 0.0d0
           sum_wat_vol  = 0.0d0
           sum_test = 0.0d0
           prot_e = 1.0d0
           wat_e = 1.0d0
          ! now do MC to find our volume
           do iAttempt = 1, nAttempts


           ! get test point
           x1 = (2.*ran1(idum) - 1.)
     x              *totalRadius
           y1 = (2.*ran1(idum) - 1.)
     x              *totalRadius
           z1 = (2.*ran1(idum) - 1.)
     x              *totalRadius


           testXYZ%x = x1
           testXYZ%y = y1
           testXYZ%z = z1

           testR2 = testXYZ * testXYZ
                     
            ! we're inside the sphere
            if (testR2 < totalR2) then
              sum_test = sum_test + 1
              outsideProt = .true.

              ! check if out point falls inside one of our protein atom spheres
              do ipr = iProtStart, iProtEnd 

               prAtom = nodeConf(ipr)

               prR    = qrms(ipr)%r                           
               prR2   = prR * prR

               dXYZ = prAtom - prMuCntr
               dXYZ = dXYZ - testXYZ
               dxyz2 = dXYZ * dXYZ

               ! yay, we're in a protein sphere
               if (dxyz2 < prR2) then
                sum_prot_vol = sum_prot_vol + 1.0d0
                outsideProt = .false.
               endif

              enddo ! End loop over all protein atoms 
            endif ! end check for being in the large sphere


            ! if our test point is outside the protein
            ! yet inside our sphere, then it must be 
            ! in the first shell volume
            if (outsideProt) then

              h1Val = 0.d0
              h2Val = 0.d0
              oVal  = 0.d0

              do jj = 1, nShell ! Loop over all 1st shell waters          
                 iOxy = lShell(jj) ! Get shell water index
                 iH1 = iOxy + 1
                 iH2 = iOxy + 2

                 oxyAtom = nodeConf(iOxy)
                 oxyR    = qrms(iOxy)%r                           
                 oxyR2   = oxyR * oxyR

                 dXYZ = oxyAtom - prMuCntr
            call imageCoords(dXYZ%x,dXYZ%y,dXYZ%z,
     x                       oxyAtom%x,oxyAtom%y,oxyAtom%z)
                 dXYZ = oxyAtom - testXYZ
!                 dXYZ = dXYZ - testXYZ
                 dxyz2 = dXYZ * dXYZ

!        real    :: h1Val,h2Val,oVal

               ! yay, we're in a oxgen sphere
               if (dxyz2 < oxyR2) then
                oVal = oVal + 1.0d0
!                sum_wat_vol = sum_wat_vol + 1.0d0
               endif
!               else  ! now we're either in H1 or H2 or nowhere

                 h1Atom = nodeConf(iH1)
                 h1R    = qrms(iH1)%r                           
                 h1R2   = h1R * h1R

                 dXYZ = h1Atom - prMuCntr
            call imageCoords(dXYZ%x,dXYZ%y,dXYZ%z,
     x                  h1Atom%x,h1Atom%y,h1Atom%z)
                 dXYZ = h1Atom - testXYZ
!                 dXYZ = dXYZ - testXYZ
                 dxyz2 = dXYZ * dXYZ

                 if (dxyz2 < h1R2) then
                  h1Val = h1Val + 1.0d0
!                  sum_wat_vol = sum_wat_vol + 1.0d0
                 endif

!                 else  ! now, were either in H2 or nowhere

                  h2Atom = nodeConf(iH2)
                  h2R    = qrms(iH2)%r   
                  h2R2   = h2R * h2R

                  dXYZ = h2Atom - prMuCntr
            call imageCoords(dXYZ%x,dXYZ%y,dXYZ%z,
     x                  h2Atom%x,h2Atom%y,h2Atom%z)
                  dXYZ = h2Atom - testXYZ
!                  dXYZ = dXYZ - testXYZ
                  dxyz2 = dXYZ * dXYZ

                  if (dxyz2 < h2R2) then
                   h2Val = h2Val + 1.0d0
!                   sum_wat_vol = sum_wat_vol + 1.0d0
                  endif

!                  endif ! end check for dxyz2 < h2R2
!                 endif ! end check for dxyz2 < h1R2
!               endif  ! end check for dxyz2 < oxyR2
              enddo ! end loop over shell waters

              wVal = oVal + h1Val + h2Val

              if(wVal .gt. 1.0) wVal = 1.0

              sum_wat_vol = sum_wat_vol + wVal

            endif


           if(mod(iAttempt,10000) .eq. 0) then
             prot_ee = sum_prot_vol/sum_test
             mcEpsProt = abs((prot_e - prot_ee)/prot_ee)
             prot_e = prot_ee

             wat_ee = sum_wat_vol/sum_test
             mcEpsWat = abs((wat_e - wat_ee)/wat_ee)
             wat_e = wat_ee
           endif

           enddo ! end attempts


           ! now, the volume of the protein and the water
           mcVolProt = (sum_prot_vol/sum_test) * sphereVol      
           mcVolWat = (sum_wat_vol/sum_test) * sphereVol

       endif ! END CHECK FOR CUDA

         call MPI_BARRIER(world,ierr)


c          print *, "about to gather the 1st data..."          

         call MPI_GATHER(mcVolProt,1,MPI_REAL,
     x                   mcVolProtTraj,1,MPI_REAL,
     x                   master,world,ierr) 

         call MPI_GATHER(mcVolWat,1,MPI_REAL,
     x                   mcVolWatTraj,1,MPI_REAL,
     x                   master,world,ierr) 

         call MPI_GATHER(mcEpsProt,1,MPI_REAL,
     x                   mcEpsProtTraj,1,MPI_REAL,
     x                   master,world,ierr) 

         call MPI_GATHER(mcEpsWat,1,MPI_REAL,
     x                   mcEpsWatTraj,1,MPI_REAL,
     x                   master,world,ierr) 

     
          if(rank .eq. master)  then           
            do iWrite = 1, nProc
                           
               write(iPVol, 5559) iConf+iWrite-1, 
     x         mcVolProtTraj(iWrite),
     x         mcEpsProtTraj(iWrite),
     x         mcVolWatTraj(iWrite),
     x         mcEpsWatTraj(iWrite)
     
            enddo 
           endif ! End check for master       
          
5559             format (i6,4(2x,E14.6))  
        return
      end subroutine calcProtV

      subroutine calcHoleV()
c    *******************************************************************
c    ** calculates the volume of a protein cavity using MC            **
c    ** integration                                                   **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: ipr, iWrite
        integer :: ii,jj,iw2, iw, idum
!        integer :: iOxy, iH1, iH2
!        integer :: massOrQ
        integer :: thisI
c        real    :: totM, totQ, prM, prQ, prMuLen, prMuLen2
        real    :: prR, prR2
        real    :: dxyz2, x1, y1, z1
        real    :: dX, dY, dZ
        real    :: sum_test, sum_hole_vol
        real    :: hole_e, hole_ee
        real, parameter :: pi = 3.14159d0
        real    :: testR2, sphereVol
c        real    :: oxyR, oxyR2
c        real    :: h1R, h1R2
c        real    :: h2R, h2R2
c        real    :: h1Val,h2Val,oVal,wVal
c        logical :: firstShell, doOutput
c        logical :: insideOxy, insideH1, insideH2
        logical :: outsideProt
c        logical :: outsideOxy, outsideH1
c        type(basic_atom)  :: prMuCntr, prAtom, tmpAtom
c        type(basic_atom)  :: oxyAtom, h1Atom, h2Atom
        type(basic_atom)  :: dXYZ, testXYZ  ! for temp storage
c    ******************************************************************* 
         
        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then
          open(unit=iHoleVol,file="holeVol.dat",status="unknown")

          open(unit=iHoleCOM,file=trim(strHoleCOM),status="old")
          do ii = 1, nConf
            read(iHoleCOM, *) thisI, holeCOMTraj(ii)%x,
     x           holeCOMTraj(ii)%y, holeCOMTraj(ii)%z
          enddo
          close(iHoleCOM)

         endif

         ! make sure everyone waits while master reads the C.O.M. data
         call MPI_BARRIER(world,ierr)
        endif 

       ! Distribute another set of hole C.O.M. to all processors  
        call MPI_SCATTER(holeCOMTraj(iConf),1,ATOM_TYPE,
     x                   holeCOM,1,ATOM_TYPE,
     x                   master,world,ierr)    


         if(useCuda) then
#ifdef ENABLE_CUDA 
c           vvvv lShell, nShell is undefined here vvvv
c            call get_shell_gpu(bX,bY,bZ,
c     x                         nodeConf,qrms,
c     x                         rBindCut,lShell,nShell,.true.)
c           ^^^^ lShell, nShell is undefined here ^^^^

           ! Calculate the protein's diople
c           prAtom%x = 0.d0
c           prAtom%y = 0.d0
c           prAtom%z = 0.d0
c           massOrQ = 0 !set to 0 to calculate the center of mass vector
c           doOutput = .false. ! we dont want all the positions copied back

c            call get_prmu_gpu(qrms,
c     x                         prAtom%x,prAtom%y,prAtom%z, 
c     x                         prMuCntr%x,prMuCntr%y,prMuCntr%z,
c     x                         totM,massOrQ,doOutput,tmpAtom)

c            prMuCntr%x = prMuCntr%x/totM
c            prMuCntr%y = prMuCntr%y/totM
c            prMuCntr%z = prMuCntr%z/totM

c            protRadius = 0.d0
c           do ipr = iProtStart, iProtEnd  ! Loop over protein atoms to find the proteins effective radius

c            prAtom = nodeConf(ipr)

c            dXYZ = prAtom - prMuCntr
c            dxyz2 = dXYZ * dXYZ

c            if (dxyz2 > protRadius) protRadius = dxyz2
            
c           enddo ! End loop over all protein atoms 

c           protRadius = sqrt(protRadius)
    
           ! add an extra 3.5A to make sure all water atoms are included
c           totalRadius = protRadius + rBindCut + 3.5d0

            ! generate random numbers for the volume calculation
            ! NOTE: this should move to the GPU, seriously.
c           do iAttempt = 1, nAttempts
c            ranX(iAttempt) = ran1(idum)
c            ranY(iAttempt) = ran1(idum)
c            ranZ(iAttempt) = ran1(idum)
c           enddo

c           vvvv lShell, nShell is undefined here vvvv
c            call get_vol_shl_gpu(bX,bY,bZ,
c     x                           nodeConf,qrms,
c     x                            lShell,nShell,
c     x                            totalRadius,
c     x                            ranX, ranY, ranZ,            
c     x                  prMuCntr%x, prMuCntr%y, prMuCntr%z,
c     x                        mcVolProt, mcVolWat,
c     x                        mcEpsProt, mcEpsWat)
c           ^^^^ lShell, nShell is undefined here ^^^^

#endif /* CUDA CODE */
         else ! we are only doing non-CUDA, CPU code

!           prMuCntr%x = 0.d0
!           prMuCntr%y = 0.d0
!           prMuCntr%z = 0.d0
!           totM = 0.d0           
           
!           do ipr = 1, iProtEnd  ! Loop over protein atoms to find the center of mass

!            prAtom = nodeConf(ipr)            
!            prM = qrms(ipr)%m
            
!            totM = totM + prM           

!            prMuCntr%x = prMuCntr%x + prM*prAtom%x
!            prMuCntr%y = prMuCntr%y + prM*prAtom%y
!            prMuCntr%z = prMuCntr%z + prM*prAtom%z
            
!           enddo ! End loop over all protein atoms 
            
!            prMuCntr%x = prMuCntr%x/totM
!            prMuCntr%y = prMuCntr%y/totM
!            prMuCntr%z = prMuCntr%z/totM

!            protRadius = 0.d0
!           do ipr = iProtStart, iProtEnd  ! Loop over protein atoms to find the proteins effective radius

!            prAtom = nodeConf(ipr)

!            dXYZ = prAtom - prMuCntr
!            dxyz2 = dXYZ * dXYZ

!            if (dxyz2 > protRadius) protRadius = dxyz2
            
!           enddo ! End loop over all protein atoms 

!           protRadius = sqrt(protRadius)
    
           ! add an extra 3.5A to make sure all water atoms are included
!           totalRadius = protRadius + rBindCut + 3.5d0
!           totalR2 = totalRadius * totalRadius

!holeRadius

           totalR2 = holeRadius * holeRadius
           sphereVol = (4.0d0/3.0d0)*pi* 
     x                 (holeRadius**3)

           sum_hole_vol  = 0.0d0
           sum_test = 0.0d0
           hole_e = 1.0d0
!           wat_e = 1.0d0
          ! now do MC to find our volume
           do iAttempt = 1, nAttempts


           ! get test point
           x1 = (2.*ran1(idum) - 1.)
     x              *holeRadius
           y1 = (2.*ran1(idum) - 1.)
     x              *holeRadius
           z1 = (2.*ran1(idum) - 1.)
     x              *holeRadius

           testXYZ%x = x1
           testXYZ%y = y1
           testXYZ%z = z1


           testR2 = testXYZ * testXYZ
                     
            ! we're inside the sphere
            if (testR2 < totalR2) then
              sum_test = sum_test + 1.d0
              outsideProt = .true.

              ! shift the test point to the C.O.M. center
!              testXYZ = testXYZ - holeCOM

              ! check if out point falls outside all protein atom spheres
              do ipr = iProtStart, iProtEnd 

               prAtom = nodeConf(ipr)

               prR    = qrms(ipr)%r                           
               prR2   = prR * prR

               dXYZ = prAtom - holeCOM
               dXYZ = dXYZ - testXYZ
               dxyz2 = dXYZ * dXYZ

               ! we're inside a protein sphere
               if (dxyz2 < prR2) then
                outsideProt = .false.
                exit
               endif

              enddo ! End loop over all protein atoms 

              if(outsideProt) then
               sum_hole_vol = sum_hole_vol + 1.0d0
              endif

            endif ! end check for being in the probe sphere


           if(mod(iAttempt,10000) .eq. 0) then
             hole_ee = sum_hole_vol/sum_test
             mcEpsHole = abs((hole_e - hole_ee)/hole_ee)
             hole_e = hole_ee
           endif

           enddo ! end attempts

           ! now, the volume of the hole
           mcVolHole = (sum_hole_vol/sum_test) * sphereVol      

       endif ! END CHECK FOR CUDA

         call MPI_BARRIER(world,ierr)


c          print *, "about to gather the 1st data..."          

         call MPI_GATHER(mcVolHole,1,MPI_REAL,
     x                   mcVolHoleTraj,1,MPI_REAL,
     x                   master,world,ierr) 

         call MPI_GATHER(mcEpsHole,1,MPI_REAL,
     x                   mcEpsHoleTraj,1,MPI_REAL,
     x                   master,world,ierr) 
     
          if(rank .eq. master)  then           
            do iWrite = 1, nProc
                           
               write(iHoleVol, 5559) iConf+iWrite-1, 
     x         mcVolHoleTraj(iWrite),
     x         mcEpsHoleTraj(iWrite)
     
            enddo 
           endif ! End check for master       
          
5559             format (i6,2(2x,E14.6))  
        return
      end subroutine calcHoleV

      subroutine calcProtMuShellM()
c    *******************************************************************
c    ** calculates Mu.M, where Mu is dipole vector of the protein     **
c    ** with respect to the center of charge, and M is dipole moment  **
c    ** of the shell of 1st waters.                                   **
c    **                                                               **
c    ** NOTE: This subroutine assumes that both calcShellM and        **
c    **       calcProtM have already been called.                     **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        type(basic_atom) :: eHat
        real  :: nodePmuShellM, masterPmuShellM(nProc)  
        real  :: eLen2, eLen
        integer :: iWrite        !
c    ******************************************************************* 
         
        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then
          open(unit=iPMuSh,file="muDotShellM.dat",status="unknown")               
         endif
        endif 

c        type(basic_atom)  :: extraPmu, extraShellM       

         eLen2 = extraPmu*extraPmu
         eLen = sqrt(eLen2)
         extraPmu%x = extraPmu%x/eLen
         extraPmu%y = extraPmu%y/eLen
         extraPmu%z = extraPmu%z/eLen

         eLen2 = extraShellM*extraShellM
         eLen = sqrt(eLen2)
         extraShellM%x = extraShellM%x/eLen
         extraShellM%y = extraShellM%y/eLen
         extraShellM%z = extraShellM%z/eLen

         ! Find the dot product of the protein dipole and the 1st shell water dipole
         nodePmuShellM = extraPmu*extraShellM

         call MPI_BARRIER(world,ierr)


c          print *, "about to gather the 1st data..."          

         call MPI_GATHER(nodePmuShellM,1,MPI_REAL,
     x                   masterPmuShellM,1,MPI_REAL,
     x                   master,world,ierr) 

     
          if(rank .eq. master)  then           
            do iWrite = 1, nProc
                           
               write(iPMuSh, 6678) iConf+iWrite-1, 
     x         masterPmuShellM(iWrite)
     
            enddo 
           endif ! End check for master       
          
6678             format (i6,2x,E14.6)                    
        return
      end subroutine calcProtMuShellM
     
      subroutine calcShellP2()
c    *******************************************************************
c    ** calculates <P2(r.e)> where e is the unit dipole vector of all **
c    ** 1st shell waters, and r is the unit vector (normal to the     **
c    ** surface) of the closest atom to a given water.                **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: iw, iw2, jj, iWrite
        integer :: nShell
        integer :: lShell(25000)  ! list of first atoms of water
        integer :: ix, iy, iz, irc          
        real    :: rIJ, rIJ2, lowestRij
        real    :: eLen, eLen2, wwX, wwY, wwZ, wwQ
        real    :: p2arg, p2Traj(nProc), p2node
        real    :: p1Traj(nProc), p1node
        logical :: firstShell
        type(basic_atom)  :: dRW, dRWi, eHat, rHat
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then          
          open(unit=iShellP2,file="shellP2.dat",status="unknown")               
         endif
        endif 
        
         nShell = 0        
         
         if(useCuda) then
#ifdef ENABLE_CUDA 
            call get_shell_gpu(bX,bY,bZ,
     x                     nodeConf,qrms,
     x                    rBindCut,lShell,nShell,.false.)

            call get_p2_gpu(bX,bY,bZ,
     x                     nodeConf,qrms,
     x                    lShell,
     x                   nShell,p1node,p2node)
#endif
         else

          do iw = iWatStart, natom, 3    ! Loop over all water-oxygen atoms
                        
           ! Find the 1st solvation shell waters
           firstShell = isFirstShell(nodeConf(iw)%x,
     x                              nodeConf(iw)%y,
     x                              nodeConf(iw)%z)
     
           ! Store the index of the oxygen of the 1st shell
           if(firstShell) then
            nShell = nShell + 1
            lShell(nShell) = iw         
           endif
          
          enddo ! End the loop over all waters     

           p2node = 0.d0
           p1node = 0.d0
          do jj = 1, nShell ! Loop over all 1st shell waters          
           iw = lShell(jj) ! Get shell water index
                    
           eHat%x = 0.d0
           eHat%y = 0.d0
           eHat%z = 0.d0 
           
           do iw2 = iw, iw+2 ! Loop over all 3 atoms           
           
           wwX = nodeConf(iw2)%x
           wwY = nodeConf(iw2)%y
           wwZ = nodeConf(iw2)%z
           wwQ = qrms(iw2)%q
           
           eHat%x = eHat%x + wwQ*wwX
           eHat%y = eHat%y + wwQ*wwY

           eHat%z = eHat%z + wwQ*wwZ           

           if (iw2 .eq. (iw+2)) then
           eLen2 = eHat*eHat
           eLen = sqrt(eLen2)
          
           eHat%x = eHat%x/eLen
           eHat%y = eHat%y/eLen
           eHat%z = eHat%z/eLen
          
            lowestRij = 9999999.d0
            do irc = 1, iProtEnd  ! Loop over protein atoms to find the closest

            dRW = nodeConf(irc) - nodeConf(iw)
            call imageCoords(dRW%x,dRW%y,dRW%z,
     x                       dRWi%x,dRWi%y,dRWi%z)
            
            rIJ2 = dRWi*dRWi
            rIJ = sqrt(rIJ2)
            
            if (rIJ .lt. lowestRij) then
              lowestRij = rIJ
              rHat%x = dRWi%x/rIJ
              rHat%y = dRWi%y/rIJ
              rHat%z = dRWi%z/rIJ
            endif
            
            enddo ! End loop over all protein atoms       
           endif  ! End If iw2=iw+2 (ie. the h2 atom)
         enddo ! End loop over all 3 water atoms
                    
              
          p2arg = eHat*rHat
          
          p1node = p1node + p2arg
          p2node = p2node + 3*(p2arg**2) - 1.d0
       enddo ! End loop over 1st shell waters

        endif ! Check for CUDA
          
          p1node = p1node/real(nShell)
          p2node = p2node/real(2*nShell)  ! the '2' is for P2 = (3x^2-1)/2
          
         call MPI_BARRIER(world,ierr)

         call MPI_GATHER(p2node,1,MPI_REAL,
     x                   p2Traj,1,MPI_REAL,
     x                   master,world,ierr) 

         call MPI_GATHER(p1node,1,MPI_REAL,
     x                   p1Traj,1,MPI_REAL,
     x                   master,world,ierr) 

         call MPI_GATHER(nShell,1,MPI_INTEGER,
     x                   nShellTraj,1,MPI_INTEGER,
     x                   master,world,ierr) 

     
          if(rank .eq. master)  then           
            do iWrite = 1, nProc
                           
               write(iShellP2, 7878) iConf+iWrite-1, 
     x         p2Traj(iWrite), p1Traj(iWrite), 
     x         nShellTraj(iWrite)
     
            enddo 
           endif ! End check for master       
          
7878             format (i6,2x,2(E14.6,2x),i6)                    
        return
      end subroutine calcShellP2  


      function checkIJ(resI,resJ)
c    *******************************************************************
c    ** finds all clusters of selected residues within a pre-defined  **
c    **  cutoff (rClusterCut)                                         **
c    **                                                               **
c    *******************************************************************            
        implicit none

        ! Function def
        logical :: checkIJ

        ! Arguments
        type(res_list_data) :: resI,resJ

        ! Locals
        logical :: testIJ 
        type(basic_atom) :: posI, posJ, dIJ, dfIJ
        integer :: ii, jj, idx, jdx
        integer :: iTotal, jTotal
        real    :: dR2, dR
c    ******************************************************************* 
     
       testIJ = .false.

       iTotal = resI%numAtomsPerRes
       jTotal = resJ%numAtomsPerRes

       iloop: do ii = 1, iTotal
         idx = resI%atomid(ii)

         posI = nodeConf(idx)

        jloop: do jj = 1, jTotal 
         jdx = resJ%atomid(jj)
         posJ = nodeConf(jdx)

         dIJ = posI - posJ

         call imageCoords(dIJ%x,dIJ%y,dIJ%z,
     x               dfIJ%x,dfIJ%y,dfIJ%z)
        
         dR2 = dfIJ*dfIJ
         dR = sqrt(dR2)

         if(dR .le. rClusterCut) then
           testIJ = .true.          
           exit iloop 
         endif

        end do jloop
       end do iloop

       checkIJ = testIJ

       return
      end function checkIJ

      recursive subroutine addNeighborList(cumRes,feederRes)
c    *******************************************************************
c    ** adds all the neighbors of one residue (feederRes) to another  **
c    ** residue (feederRes) of the neighborList matrix                **
c    **                                                               **
c    *******************************************************************            
        implicit none

        ! Arguments
        integer :: cumRes,feederRes

        ! Locals
        integer :: ii, jj, kk, idx, jdx
        integer :: jNc, iNc, iRes, jRes
        integer :: tmpIdx
        logical :: foundRes, notIn
c    ******************************************************************* 
     
       idx = cumRes
       jdx = feederRes

       iNc = neighborCount(idx)
       jNc = neighborCount(jdx)

       jloop: do jj = 1, jNc
         jRes = neighborList(jdx,jj)
         foundRes = .false.

         iloop: do ii = 1, iNc
          iRes = neighborList(idx,ii)
         
          if(iRes .eq. jRes) then
           foundRes = .true.
           exit iloop
          endif

         end do iloop ! ii loop

         if ((.not. foundRes) 
     x     .and. (idx .ne. jRes)) then

            notIn = .true.
            do kk = 1, neighborCount(idx)
              if(jRes .eq. neighborList(idx,kk)) then
                notIn = .false.
              endif
                
            enddo ! end check for notin

            if(notIn) then
             neighborCount(idx) = neighborCount(idx) + 1
             tmpIdx = neighborCount(idx)
             neighborList(idx,tmpIdx) = jRes          
             call addNeighborList(idx, jRes) !!! new addition
             checklist(jRes) = .false. !!! new addition
            endif

         endif

       end do jloop ! jj loop


       return
      end subroutine addNeighborList

      subroutine findClusters()
c    *******************************************************************
c    ** finds all clusters of selected residues within a pre-defined  **
c    **  cutoff (rClusterCut)                                         **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: ii, jj, iWrite
        integer :: neighbor
        logical :: foundIJCluster
        real    :: avgNagg
        real    :: nClustRes
        real    :: n_in_micelle
        real    :: n_micelle
        real    :: n_agg
        real    :: n_tot
        real    :: n_cluster
        real    :: avgNaggTraj(nProc)
        real    :: avgMicelle
        real    :: avgMicelleTraj(nProc)
        real    :: avg
        real    :: avgTraj(nProc)
        integer :: nAgg2, nCluster, myCount
        type(dstring) :: clusterOut
        type(res_list_data) :: resI, resJ
        character(len=250) :: aggString
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then          
          open(unit=iClusterOut,file=trim(strClusterOut),
     x      form="formatted",status="unknown")   

           write(iClusterOut,1124) ,rClusterCut
1124             format ("# ANALYZED USING A CUTOFF OF ",F5.2)

           clusterOut = new("# WARNING:  USING type") + 
     x                  new(" ") + new(strSelect) + new(";") +
     x                  new(" ") + new("MAKE SURE")
           write(iClusterOut,'(a)') clusterOut%data(1:clusterOut%length)

           clusterOut = new("# THESE ARE THE RIGHT TYPES!!!!")
           write(iClusterOut,'(a)') clusterOut%data(1:clusterOut%length)

         endif
        endif 

!========================================================!
         ! initialize the variables for each step
         neighborCount(:) = 0
         checkList(:) = .true.

         ! double loop over all pairs, minus the diagonal
         do ii = 2, nClusterRes
           do jj = 1, ii - 1
               resI = resList(ii)
               resJ = resList(jj)

               foundIJCluster = checkIJ(resI,resJ)

               if(foundIJCluster) then
                 neighborCount(ii) = neighborCount(ii) + 1
                 neighborCount(jj) = neighborCount(jj) + 1

                 neighborList(ii,neighborCount(ii)) = jj
                 neighborList(jj,neighborCount(jj)) = ii

               endif 
           enddo ! end jj
         enddo ! end ii
!========================================================!

!========================================================!
         ! double loop over selected residues and neighbors
         do ii = 1, nClusterRes

           ! only bother if we are not  monomoers
          if(checkList(ii)) then
            if(neighborCount(ii) .gt. 0) then
  
              do jj = 1, neighborCount(ii)

               neighbor = neighborList(ii,jj)

               checklist(neighbor) = .false. !!! new spot
               call addNeighborList(ii, neighbor)
               ! checklist(neighbor) = .false. !!! orig spot

              enddo ! end jj
            endif ! end check for neighbors .gt. 0)
           endif ! end check for checklist
         enddo ! end ii
!========================================================!


!========================================================!
         n_micelle = 0.d0
         n_tot = 0.d0
         n_in_micelle = 0.d0

         ! build aggList for resturning to master
         ! double loop over selected residues and neighbors
         nAgg = 0            
         nClustRes = 0.d0
         do ii = 1, nClusterRes          
          if( checkList(ii) ) then
          nCluster = neighborCount(ii) + 1

          nAgg = nAgg + 1
          aggList(nAgg)%nAgg = nCluster

          do jj = 1, neighborCount(ii)          
             aggList(nAgg)%resid(jj + 1) = 
     x        resList(neighborList(ii,jj))%resid - 1
          enddo
          aggList(nAgg)%resid(1) = resList(ii)%resid - 1 ! self term

          n_tot = n_tot + real(nCluster,kind=4)
          ! only count micelles if they are above a 5-mer cutoff
           if( neighborCount(ii) .ge. 5) then
            
            n_micelle = n_micelle + 1.d0
            n_in_micelle = n_in_micelle + 
     x           real(nCluster,kind=4)

           endif ! end check for minimum cutoff

          endif ! end check for still being on the checklist

         enddo ! end ii
!========================================================!


         n_micelle = real(nAgg,kind=4)
         n_cluster =  real(nAgg,kind=4)

         avg = n_tot/n_cluster
         avgMicelle = n_in_micelle/n_micelle

!========================================================!
         if(doGyr .or. doInert) call aggImage() ! Wrap all positions to be in the box as the cluster
         if(doGyr) call calcGyrCluster()        ! Calculate the radius of gyration for all aggregates
         if(doInert) call calcInertCluster()    ! Calculate the moment of intertia tensor for all aggregates

!========================================================!
         ! Gather all cluster data on the master node
         call MPI_GATHER(nAgg,1,MPI_INTEGER,
     x               nAggTraj,1,MPI_INTEGER,
     x               master,world,ierr)

         call MPI_GATHER(aggList,1,AGG_TYPE,
     x              aggListTraj,1,AGG_TYPE,
     x              master,world,ierr)

         call MPI_GATHER(avg,1,MPI_REAL,
     x               avgTraj,1,MPI_REAL,
     x               master,world,ierr)

         call MPI_GATHER(avgMicelle,1,MPI_REAL,
     x               avgMicelleTraj,1,MPI_REAL,
     x               master,world,ierr)

         if(doGyr) then
           call MPI_GATHER(aggGyr,1,AGG_PARM_TYPE,
     x                 aggGyrTraj,1,AGG_PARM_TYPE,
     x               master,world,ierr)

           call MPI_GATHER(aggAsp,1,AGG_PARM_TYPE,
     x                 aggAspTraj,1,AGG_PARM_TYPE,
     x               master,world,ierr)

           call MPI_GATHER(aggR1,1,AGG_PARM_TYPE,
     x                 aggR1Traj,1,AGG_PARM_TYPE,
     x               master,world,ierr)
           call MPI_GATHER(aggR2,1,AGG_PARM_TYPE,
     x                 aggR2Traj,1,AGG_PARM_TYPE,
     x               master,world,ierr)
           call MPI_GATHER(aggR3,1,AGG_PARM_TYPE,
     x                 aggR3Traj,1,AGG_PARM_TYPE,
     x               master,world,ierr)
!                  deallocate(aggEccen)
         endif

         if(doInert) then
           call MPI_GATHER(aggEccen,1,AGG_PARM_TYPE,
     x                 aggEccenTraj,1,AGG_PARM_TYPE,
     x               master,world,ierr)

           call MPI_GATHER(aggA2,1,AGG_PARM_TYPE,
     x                 aggA2Traj,1,AGG_PARM_TYPE,
     x               master,world,ierr)
           call MPI_GATHER(aggB2,1,AGG_PARM_TYPE,
     x                 aggB2Traj,1,AGG_PARM_TYPE,
     x               master,world,ierr)
           call MPI_GATHER(aggC2,1,AGG_PARM_TYPE,
     x                 aggC2Traj,1,AGG_PARM_TYPE,
     x               master,world,ierr)

!                  deallocate(aggEccen)
         endif

!========================================================!

          if(rank .eq. master)  then

!           jj = 0
            do iWrite = 1, nProc
            ! --++++++++++++++++++++++++--
!             jj = jj + 1

             write(iClusterOut,'(a)',advance="no") "FRAME"
             write(iClusterOut,'(i5)',advance="no")
     x           (iConf+iWrite-2) ! 2, instead of 1, to make FRAME 0-based
             write(iClusterOut,'(a)',advance="yes") ""

              do ii = 1, nAggTraj(iWrite) 

               write(iClusterOut,'(a)',advance="no") "AGG"
               write(iClusterOut,'(i5)',advance="no") 
     x         aggListTraj((iWrite - 1)*MAX_AGG + ii)%nAgg
               write(iClusterOut,'(a)',advance="no") " "

               do jj = 1, 
     x         aggListTraj((iWrite - 1)*MAX_AGG + ii)%nAgg

                write(iClusterOut,'(i5)',advance="no") 
     x                  aggListTraj((iWrite - 1)*MAX_AGG + ii)
     x                             %resid(jj)
               enddo ! end loop over all jj

              if(doGyr .or. doInert) then
               write(iClusterOut,'(a)',advance="no") " # "
              endif

              if(doGyr) then
c               write(iClusterOut,'(f12.6)',advance="no") 
               write(iClusterOut,'(E10.4)',advance="no") 
     x         aggGyrTraj((iWrite - 1)*MAX_AGG + ii)%parm
               write(iClusterOut,'(a)',advance="no") " "

               write(iClusterOut,'(E10.4)',advance="no") 
     x         aggR1Traj((iWrite - 1)*MAX_AGG + ii)%parm
               write(iClusterOut,'(a)',advance="no") " "
               write(iClusterOut,'(E10.4)',advance="no") 
     x         aggR2Traj((iWrite - 1)*MAX_AGG + ii)%parm
               write(iClusterOut,'(a)',advance="no") " "
               write(iClusterOut,'(E10.4)',advance="no") 
     x         aggR3Traj((iWrite - 1)*MAX_AGG + ii)%parm
               write(iClusterOut,'(a)',advance="no") " "

               write(iClusterOut,'(E10.4)',advance="no") 
     x         aggAspTraj((iWrite - 1)*MAX_AGG + ii)%parm
               write(iClusterOut,'(a)',advance="no") " "

              endif

              if(doInert) then
               write(iClusterOut,'(E10.4)',advance="no") 
     x         aggA2Traj((iWrite - 1)*MAX_AGG + ii)%parm
               write(iClusterOut,'(a)',advance="no") " "
               write(iClusterOut,'(E10.4)',advance="no") 
     x         aggB2Traj((iWrite - 1)*MAX_AGG + ii)%parm
               write(iClusterOut,'(a)',advance="no") " "
               write(iClusterOut,'(E10.4)',advance="no") 
     x         aggC2Traj((iWrite - 1)*MAX_AGG + ii)%parm
               write(iClusterOut,'(a)',advance="no") " "

c               write(iClusterOut,'(f8.6)',advance="no") 
               write(iClusterOut,'(E10.4)',advance="no") 
     x         aggEccenTraj((iWrite - 1)*MAX_AGG + ii)%parm
               write(iClusterOut,'(a)',advance="no") " "
              endif


            write(iClusterOut,'(a)',advance="yes") ""
            enddo ! loop over all ii

              write(iClusterOut,1126) avgTraj(iWrite)
              write(iClusterOut,1127) avgMicelleTraj(iWrite)

             enddo ! loop over all processor data
            ! --++++++++++++++++++++++++--
          endif ! check for master

!        call MPI_BARRIER(world,ierr)

c1126             format ("AVG ",F5.2)
c1127             format ("AVGMICELLE ",F5.2)

1126             format ("AVG ",E10.4)
1127             format ("AVGMICELLE ",E10.4)
        return
      end subroutine findClusters

      function findResInList(pdbID)
c    *******************************************************************
c    ** finds a residue in the list by index of the resList. the      **
c    ** index is returned.                                            **
c    **                                                               **
c    **  principle variabes                                           **
c    **      integer :: pdbID - pdb ID to search for in the resList   **
c    **                                                               **
c    *******************************************************************            
        implicit none      

        ! func def
        integer :: findResInList

        ! Args
        integer :: pdbID
        
        ! Locals
        integer :: ii, resID
c    ******************************************************************* 

        resID = 0
        do ii = 1, nClusterRes
          if(resList(ii)%resID .eq. pdbID) then
           resID = ii
           exit 
          endif 
        enddo

        findResInList = resID

        return
      end function findResInList


      subroutine aggImage()
c    *******************************************************************
c    ** reimages the positions of individual monomers to the position **
c    ** of the COM of their corresponding aggregate.                  **
c    **                                                               **
c    *******************************************************************            
        implicit none      

        ! Args
        
        ! Locals
        integer :: ii,jj,kk
        integer :: idx, kdx
        integer :: nAggInMic, nAtomsInMon
        integer :: aggID, aggTest, monoID, monoTest
        real    :: mTot, mi
        type(basic_atom) :: atomCOM, atomPOS, aggCOM
        type(basic_atom) :: atomMPos 
        type(basic_atom) :: dIJ, dfIJ
c    ******************************************************************* 

        if(iConf .eq. 1) then
          allocate(aggPositions(MAX_AGG))
        endif
     
        
        do ii = 1, nAgg

          nAggInMic = aggList(ii)%nAgg

          !aggPositions(ii)%nRes = nAggInMic

          aggTest = aggList(ii)%resid(1)+1
          aggID = findResInList(aggTest)   

          ! old version with resList
          !idx = resList(aggID)%atomid(1)
          idx = fullResList(aggID)%atomid(1)

          atomCOM = nodeConf(idx)

          mTot = 0.d0

          aggCOM%x = 0.d0
          aggCOM%y = 0.d0
          aggCOM%z = 0.d0

          do jj = 1, nAggInMic


           ! get the id of the individual monomer
           monoTest = aggList(ii)%resid(jj) + 1
           monoID   = findResInList(monoTest)

           ! old version in resList
           !nAtomsInMon = resList(monoID)%numAtomsPerRes
           nAtomsInMon = fullResList(monoID)%numAtomsPerRes

           do kk = 1, nAtomsInMon

           ! old version with resList
           !kdx = resList(monoID)%atomid(kk)
           kdx = fullResList(monoID)%atomid(kk)

           atomPOS = nodeConf(kdx)

           mi = qrms(kdx)%m
           mTot = mTot + mi

           dIJ = atomPOS - atomCOM
           call imageCoords(dIJ%x,dIJ%y,dIJ%z,
     x                      dfIJ%x,dfIJ%y,dfIJ%z)                     

           !aggPositions(ii)%resPos(jj)%atomPos(kk) = dfIJ
           !aggPositions(ii)%resPos(jj)%atomMass(kk) = mi

           atomMPos%x = dfIJ%x * mi
           atomMPos%y = dfIJ%y * mi
           atomMPos%z = dfIJ%z * mi

           aggCOM = aggCOM + atomMPos

           ! get positions 

           enddo ! end kk (atoms in mon)

          enddo ! end jj (monomers in agg)

          aggCOM%x = (1.d0/mTot)*aggCOM%x
          aggCOM%y = (1.d0/mTot)*aggCOM%y
          aggCOM%z = (1.d0/mTot)*aggCOM%z

          aggPositions(ii)%aggCOM = aggCOM + atomCOM

        enddo! end ii (aggregates)

        return
      end subroutine aggImage

      subroutine calcGyrCluster()
c    *******************************************************************
c    ** calculates the radius of gyration for ALL clusters found by   **
c    ** using the residue selection.                                  **
c    **                                                               **
c    *******************************************************************            
        implicit none      

        ! Args
        
        ! Locals
        integer :: ii,jj,kk
        integer :: idx, kdx
        integer :: lowIDx,midIDx,highIDx
        integer :: nAggInMic, nAtomsInMon
        integer :: aggID, aggTest, monoID, monoTest
        real    :: mTot, mi, dfIJ2
        integer :: iDUMMY, laINFO
        integer, parameter :: iWORK=48
        real*8  :: gyrMat(3,3)
        real*8  :: eig_Vec(3,3)
        real*8  :: eig_valR(3), eig_valI(3)
        real*8  :: laWORK(iWORK), laDUMMY
        real*8  :: lowR2, midR2, highR2
        real*8  :: oneByN
        real*8  :: nFloat
        real    :: s11,s22,s33
        real    :: s12,s23,s13
        type(basic_atom) :: atomCOM, atomPOS, aggCOM
        type(basic_atom) :: atomMPos 
        type(basic_atom) :: dIJ, dfIJ
c    ******************************************************************* 

        !if(iConf .eq. 1) then
        !  allocate(aggGyr(MAX_AGG))
        !endif
     
        
        do ii = 1, nAgg

          nAggInMic = aggList(ii)%nAgg
          aggCOM = aggPositions(ii)%aggCOM

          aggGyr(ii)%parm = 0.d0
          mTot = 0.d0
          gyrMat(:,:) = 0.d0 
          nFloat = 0.d0

          do jj = 1, nAggInMic


           ! get the id of the individual monomer
           monoTest = aggList(ii)%resid(jj) + 1
           monoID   = findResInList(monoTest)

           ! old version with the resList
           !nAtomsInMon = resList(monoID)%numAtomsPerRes
           nAtomsInMon = fullResList(monoID)%numAtomsPerRes

           do kk = 1, nAtomsInMon

           nFloat = nFloat + 1.d0

           ! old version with the resList
           !kdx = resList(monoID)%atomid(kk)
           kdx = fullResList(monoID)%atomid(kk)

           atomPOS = nodeConf(kdx)

!           mi = qrms(kdx)%m
!           mTot = mTot + mi

           dIJ = atomPOS - aggCOM
           call imageCoords(dIJ%x,dIJ%y,dIJ%z,
     x                      dfIJ%x,dfIJ%y,dfIJ%z)                     

!           dfIJ2 = dfIJ * dfIJ

           ! calculate the gyration tensor
           s11 = dfIJ%x*dfIJ%x
           gyrMat(1,1) = gyrMat(1,1) +
     x      real(s11,kind=8)

           s22 = dfIJ%y*dfIJ%y
           gyrMat(2,2) = gyrMat(2,2) +
     x      real(s22,kind=8)

           s33 = dfIJ%z*dfIJ%z
           gyrMat(3,3) = gyrMat(3,3) +
     x      real(s33,kind=8)

           s12 = dfIJ%x*dfIJ%y
           gyrMat(1,2) = gyrMat(1,2) +
     x      real(s12,kind=8)

           s13 = dfIJ%x*dfIJ%z
           gyrMat(1,3) = gyrMat(1,3) +
     x      real(s13,kind=8)

           s23 = dfIJ%y*dfIJ%z
           gyrMat(2,3) = gyrMat(2,3) +
     x      real(s23,kind=8)

!           aggGyr(ii)%parm = aggGyr(ii)%parm +
!     x       dfIJ2 * mi

           enddo ! end kk (atoms in mon)

          enddo ! end jj (monomers in agg)

          ! normalize the gyration tensor
          oneByN = 1.d0/nFloat
          gyrMat(1,1) = oneByN*gyrMat(1,1)
          gyrMat(2,2) = oneByN*gyrMat(2,2)
          gyrMat(3,3) = oneByN*gyrMat(3,3)
          gyrMat(1,2) = oneByN*gyrMat(1,2)
          gyrMat(1,3) = oneByN*gyrMat(1,3)
          gyrMat(2,3) = oneByN*gyrMat(2,3)

          ! account for symmetry in the gyration tensor
          gyrMat(2,1) = gyrMat(1,2)
          gyrMat(3,1) = gyrMat(1,3)
          gyrMat(3,2) = gyrMat(2,3)

         ! diagonalize intertia matrix and store the eigenvectors in
         ! aMat(i=1,3 ; j) for each j-th eigenvector.
         ! 'eig_valR' has the eigenvalues.   
         iDUMMY = 1
         call dgeev('N','V',3,gyrMat,
     x         3,eig_valR,eig_valI,
     x              laDUMMY,iDUMMY,
     x        eig_Vec,3,laWORK,iWORK,laINFO)


          lowR2 = 999999999999999999.d0
          highR2 = -999999999999999999.d0

          do jj = 1, 3

!           avgI = avgI + eig_valR(jj)

           if(eig_valR(jj) .lt. lowR2) then 
             lowR2 = eig_valR(jj)
             lowIDx = jj
           endif

           if(eig_valR(jj) .gt. highR2) then 
             highR2 = eig_valR(jj)
             highIDx = jj
           endif

          enddo


          do jj = 1, 3
            if( (jj .ne. highIDx) .and. 
     x          (jj .ne. lowIDx) ) then
             midR2 = eig_valR(jj)           
            endif
          enddo

          aggAsp(ii)%parm = 
     x    real((lowR2 - 0.5d0*(highR2 + midR2))
     x                       ,kind=4)

          aggR1(ii)%parm = highR2
          aggR2(ii)%parm = midR2
          aggR3(ii)%parm = LowR2

          aggGyr(ii)%parm = sqrt( highR2 + 
     x          midR2 + lowR2)

!          aggGyr(ii)%parm = (1.d0/mTot)* 
!     x                   aggGyr(ii)%parm


        enddo! end ii (aggregates)

        return
      end subroutine calcGyrCluster

      subroutine calcRDens()
c    *******************************************************************
c    ** calculates the radial density profile centered about the COM  **
c    ** of the individual aggregates.                                 **
c    **                                                               **
c    ** TO DO: Remove all extraneous variables and actually write     **
c    **        the code.                                              **
c    *******************************************************************            
        implicit none      

        ! Args
        
        ! Locals
        integer :: ii,jj,kk
        integer :: idx, kdx
        integer :: nAggInMic, nAtomsInMon
        integer :: aggID, aggTest, monoID, monoTest
        integer :: lowIDx,midIDx,highIDx
        integer :: iDUMMY, laINFO
        integer, parameter :: iWORK=48
        real    :: mTot, mi, dfIJ2
        real*8  :: gyrMat(3,3)
        real*8  :: eig_Vec(3,3)
        real*8  :: eig_valR(3), eig_valI(3)
        real*8  :: laWORK(iWORK), laDUMMY
        real*8  :: lowR2, midR2, highR2
        real*8  :: oneByN
        real*8  :: nFloat
        real    :: s11,s22,s33
        real    :: s12,s23,s13
        type(basic_atom) :: atomCOM, atomPOS, aggCOM
        type(basic_atom) :: atomMPos 
        type(basic_atom) :: dIJ, dfIJ
c    ******************************************************************* 

        !if(iConf .eq. 1) then
        !  allocate(aggGyr(MAX_AGG))
        !endif
     
        
        do ii = 1, nAgg

          nAggInMic = aggList(ii)%nAgg
          aggCOM = aggPositions(ii)%aggCOM

          aggGyr(ii)%parm = 0.d0
          mTot = 0.d0
          gyrMat(:,:) = 0.d0 
          nFloat = 0.d0

          do jj = 1, nAggInMic


           ! get the id of the individual monomer
           monoTest = aggList(ii)%resid(jj) + 1
           monoID   = findResInList(monoTest)

           ! old version with the resList
           !nAtomsInMon = resList(monoID)%numAtomsPerRes
           nAtomsInMon = fullResList(monoID)%numAtomsPerRes

           do kk = 1, nAtomsInMon

           nFloat = nFloat + 1.d0

           ! old version with the resList
           !kdx = resList(monoID)%atomid(kk)
           kdx = fullResList(monoID)%atomid(kk)

           atomPOS = nodeConf(kdx)

!           mi = qrms(kdx)%m
!           mTot = mTot + mi

           dIJ = atomPOS - aggCOM
           call imageCoords(dIJ%x,dIJ%y,dIJ%z,
     x                      dfIJ%x,dfIJ%y,dfIJ%z)                     

!           dfIJ2 = dfIJ * dfIJ

           ! calculate the gyration tensor
           s11 = dfIJ%x*dfIJ%x
           gyrMat(1,1) = gyrMat(1,1) +
     x      real(s11,kind=8)

           s22 = dfIJ%y*dfIJ%y
           gyrMat(2,2) = gyrMat(2,2) +
     x      real(s22,kind=8)

           s33 = dfIJ%z*dfIJ%z
           gyrMat(3,3) = gyrMat(3,3) +
     x      real(s33,kind=8)

           s12 = dfIJ%x*dfIJ%y
           gyrMat(1,2) = gyrMat(1,2) +
     x      real(s12,kind=8)

           s13 = dfIJ%x*dfIJ%z
           gyrMat(1,3) = gyrMat(1,3) +
     x      real(s13,kind=8)

           s23 = dfIJ%y*dfIJ%z
           gyrMat(2,3) = gyrMat(2,3) +
     x      real(s23,kind=8)

!           aggGyr(ii)%parm = aggGyr(ii)%parm +
!     x       dfIJ2 * mi

           enddo ! end kk (atoms in mon)

          enddo ! end jj (monomers in agg)

          ! normalize the gyration tensor
          oneByN = 1.d0/nFloat
          gyrMat(1,1) = oneByN*gyrMat(1,1)
          gyrMat(2,2) = oneByN*gyrMat(2,2)
          gyrMat(3,3) = oneByN*gyrMat(3,3)
          gyrMat(1,2) = oneByN*gyrMat(1,2)
          gyrMat(1,3) = oneByN*gyrMat(1,3)
          gyrMat(2,3) = oneByN*gyrMat(2,3)

          ! account for symmetry in the gyration tensor
          gyrMat(2,1) = gyrMat(1,2)
          gyrMat(3,1) = gyrMat(1,3)
          gyrMat(3,2) = gyrMat(2,3)

         ! diagonalize intertia matrix and store the eigenvectors in
         ! aMat(i=1,3 ; j) for each j-th eigenvector.
         ! 'eig_valR' has the eigenvalues.   
         iDUMMY = 1
         call dgeev('N','V',3,gyrMat,
     x         3,eig_valR,eig_valI,
     x              laDUMMY,iDUMMY,
     x        eig_Vec,3,laWORK,iWORK,laINFO)


          lowR2 = 999999999999999999.d0
          highR2 = -999999999999999999.d0

          do jj = 1, 3

!           avgI = avgI + eig_valR(jj)

           if(eig_valR(jj) .lt. lowR2) then 
             lowR2 = eig_valR(jj)
             lowIDx = jj
           endif

           if(eig_valR(jj) .gt. highR2) then 
             highR2 = eig_valR(jj)
             highIDx = jj
           endif

          enddo


          do jj = 1, 3
            if( (jj .ne. highIDx) .and. 
     x          (jj .ne. lowIDx) ) then
             midR2 = eig_valR(jj)           
            endif
          enddo

          aggAsp(ii)%parm = 
     x    real((lowR2 - 0.5d0*(highR2 + midR2))
     x                       ,kind=4)

          aggR1(ii)%parm = highR2
          aggR2(ii)%parm = midR2
          aggR3(ii)%parm = LowR2

          aggGyr(ii)%parm = sqrt( highR2 + 
     x          midR2 + lowR2)

!          aggGyr(ii)%parm = (1.d0/mTot)* 
!     x                   aggGyr(ii)%parm


        enddo! end ii (aggregates)

        return
      end subroutine calcRDens

      subroutine calcInertCluster()
c    *******************************************************************
c    ** calculates the diagonalized moment of intertia tensor from    **
c    ** ALL clusters found by using the residue selection.            **
c    **                                                               **
c    *******************************************************************            
        implicit none      

        ! Args
        
        ! Locals
        integer :: ii,jj,kk
        integer :: idx, kdx
        integer :: lowIDx,midIDx,highIDx
        integer :: nAggInMic, nAtomsInMon
        integer :: aggID, aggTest, monoID, monoTest
        integer :: iDUMMY, laINFO
        integer, parameter :: iWORK=48
        real    :: mTot, mi, dfIJ2
        real*8  :: inertMat(3,3)
        real*8  :: eig_Vec(3,3)
        real*8  :: eig_valR(3), eig_valI(3)
        real*8  :: laWORK(iWORK), laDUMMY
        real*8  :: lowI, midI, highI, avgI
        real    :: i11,i12,i13
        real    :: i22,i23, i33
        type(basic_atom) :: atomCOM, atomPOS, aggCOM
        type(basic_atom) :: atomMPos 
        type(basic_atom) :: dIJ, dfIJ
c    ******************************************************************* 
        !if(iConf .eq. 1) then
        !  allocate(aggEccen(MAX_AGG))
        !endif
     
        
        do ii = 1, nAgg

          nAggInMic = aggList(ii)%nAgg
          aggCOM = aggPositions(ii)%aggCOM

          inertMat(:,:) = 0.d0 

          do jj = 1, nAggInMic


           ! get the id of the individual monomer
           monoTest = aggList(ii)%resid(jj) + 1
           monoID   = findResInList(monoTest)

           ! old version with resList
           !nAtomsInMon = resList(monoID)%numAtomsPerRes
           nAtomsInMon = fullResList(monoID)%numAtomsPerRes

           do kk = 1, nAtomsInMon

           ! old version with resList
           !kdx = resList(monoID)%atomid(kk)
           kdx = fullResList(monoID)%atomid(kk)

           atomPOS = nodeConf(kdx)

           mi = qrms(kdx)%m
           mTot = mTot + mi

           dIJ = atomPOS - aggCOM
           call imageCoords(dIJ%x,dIJ%y,dIJ%z,
     x                      dfIJ%x,dfIJ%y,dfIJ%z)                     


           ! calculate the moment of interia tensor
           i11 = mi*(dfIJ%y*dfIJ%y + dfIJ%z*dfIJ%z)
           inertMat(1,1) = inertMat(1,1) +
     x      real(i11,kind=8)

           i22 = mi*(dfIJ%x*dfIJ%x + dfIJ%z*dfIJ%z)
           inertMat(2,2) = inertMat(2,2) +
     x      real(i22,kind=8)

           i33 = mi*(dfIJ%x*dfIJ%x + dfIJ%y*dfIJ%y)
           inertMat(3,3) = inertMat(3,3) +
     x      real(i33,kind=8)

           ! swap sign later on 1,2
           i12 = mi*(dfIJ%x*dfIJ%y)
           inertMat(1,2) = inertMat(1,2) +
     x      real(i12,kind=8)

           ! swap sign later on 1,3
           i13 = mi*(dfIJ%x*dfIJ%z)
           inertMat(1,3) = inertMat(1,3) +
     x      real(i13,kind=8)

           ! swap sign later on 2,3
           i23 = mi*(dfIJ%y*dfIJ%z)
           inertMat(2,3) = inertMat(2,3) +
     x      real(i23,kind=8)

           enddo ! end kk (atoms in mon)

          enddo ! end jj (monomers in agg)


          ! change sign to account for previous signage
          inertMat(1,2) = -1.d0*inertMat(1,2)
          inertMat(1,3) = -1.d0*inertMat(1,3)
          inertMat(2,3) = -1.d0*inertMat(2,3)

          ! account for symmetry
          inertMat(2,1) = inertMat(1,2)
          inertMat(3,1) = inertMat(1,3)
          inertMat(3,2) = inertMat(2,3)

         ! diagonalize intertia matrix and store the eigenvectors in
         ! aMat(i=1,3 ; j) for each j-th eigenvector.
         ! 'eig_valR' has the eigenvalues.   
         iDUMMY = 1
         call dgeev('N','V',3,inertMat,
     x         3,eig_valR,eig_valI,
     x              laDUMMY,iDUMMY,
     x        eig_Vec,3,laWORK,iWORK,laINFO)


          lowI = 999999999999999999.d0
          highI = -999999999999999999.d0
          avgI = 0.d0
          do jj = 1, 3

           avgI = avgI + eig_valR(jj)

           if(eig_valR(jj) .lt. lowI) then 
             lowI = eig_valR(jj)
             lowIDx = jj
           endif

           if(eig_valR(jj) .gt. highI) then 
             highI = eig_valR(jj)
             highIDx = jj
           endif

          enddo


          do jj = 1, 3
            if( (jj .ne. highIDx) .and. 
     x          (jj .ne. lowIDx) ) then
             midI = eig_valR(jj)           
            endif
          enddo

          avgI = avgI*(0.333333333333333d0)

          aggEccen(ii)%parm = 
     x          real((1.d0 - lowI/avgI)
     x                       ,kind=4)

          aggA2(ii)%parm = HighI
          aggB2(ii)%parm = midI
          aggC2(ii)%parm = LowI
c          aa = (5.0d0/2.0d0/mTot) * 
c     x                (highI - lowI + midI)

c          if(aa < 0.00000000d0) then
c             aa = 0.000000000
c          endif
c          aggA2(ii)%parm = 
c     x          real( sqrt(aa),kind=4)

c          bb = (5.0d0/2.0d0/mTot) * 
c     x                (highI - midI + lowI)
c          if(bb < 0.00000000d0) then
c             bb = 0.000000000
c          endif
c          aggB2(ii)%parm = 
c     x          real( sqrt(bb),kind=4)

c          cc = (5.0d0/2.0d0/mTot) * 
c     x                (lowI - highI + midI)
c          if(cc < 0.00000000d0) then
c             cc = 0.000000000
c          endif
c          aggC2(ii)%parm = 
c     x          real( sqrt(cc),kind=4)

        enddo! end ii (aggregates)

        return
      end subroutine calcInertCluster

      
      subroutine imageCoords(xxI,yyI,zzI,xxF,yyF,zzF)
c    *******************************************************************
c    ** returns reimaged coordinates for either a truncated or        **
c    ** orthorhombic cell type.                                       **
c    **                                                               **
c    ** principal variables:                                          **
c    **          real    :: xxI,yyI,zzI - initial x,y,z coords.       **
c    **          real    :: xxF,yyF,zzF - final x,y,z coords.         **
c    **                                                               **
c    *******************************************************************            
        implicit none      

        ! Args
        real    :: xxI,yyI,zzI,xxF,yyF,zzF
        
        ! Locals
        integer :: ix, iy, iz               
        integer :: xInt, yInt, zInt         
        real    :: dX, dY, dZ       
        real*8  :: dfrx, dfry, dfrz, delr2, rIJ    
        real*8  :: lowDX, lowDY, lowDZ                
c    ******************************************************************* 

         ! Reimage the box in the for truncated or orthorhombic boxes
         dX = xxI
         dY = yyI
         dZ = zzI
         
         if (isTruncated) then          

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          rIJ = dX**2 + dY**2 + dZ**2
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
              lowDX = dX
              lowDY = dY              
              lowDZ = dZ              
            endif

            enddo
           enddo
          enddo
          
           dX = lowDX
           dY = lowDY
           dZ = lowDZ                      
           
         else ! its an orthorhombic box
          
          xInt = ANINT(dX/bX)
          yInt = ANINT(dY/bY)
          zInt = ANINT(dZ/bZ)

          dX = dX - xInt*bX
          dY = dY - yInt*bY
          dZ = dZ - zInt*bZ         

         endif ! End check for truncation
         
         xxF = dX
         yyF = dY
         zzF = dZ
         
       return      
      end subroutine imageCoords
      
      function cellBindE(wwX,wwY,wwZ,wwQ,iWNum)
c    *******************************************************************
c    ** determines the full cell binding energy depending for a given **
c    ** atom in a binding water.                                      **
c    **                                                               **
c    ** principal variables:                                          **
c    **          real    :: wwX, wwY, wwZ - x,y,z coords. of h2o      **
c    **          integer :: iWNum - the index of the water position   **
c    **                                                               **
c
c
c      IS THIS CODE EVEN USED????
c
c    *******************************************************************            
        implicit none
       
        ! Function def
        real :: cellBindE

        ! Args
        real    :: wwX, wwY, wwZ, wwQ
        integer :: iWNum
        
        ! Locals
        integer :: iww
        integer :: ix, iy, iz        
        integer :: xInt, yInt, zInt        
        real*8  :: dfrx, dfry, dfrz, delr2        
        real*8  :: wX2, wY2, wZ2, wQ2, wQ   
        real*8  :: rIJ, dX, dY, dZ
        real*8  :: wX, wY, wZ
c    ******************************************************************* 
         cellBindE = 0.d0
         wX = real(wwX,kind=8)
         wY = real(wwY,kind=8)
         wZ = real(wwZ,kind=8)
         wQ = real(wwQ,kind=8)
         
         do iww = iWatStart, natom  ! Loop over protein atoms
         
         if (iWNum .ne. iww) then  
          wX2 = real(nodeConf(iww)%x,kind=8)
          wY2 = real(nodeConf(iww)%y,kind=8)
          wZ2 = real(nodeConf(iww)%z,kind=8)
          wQ2 = real(qrms(iww)%q,kind=8)

         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then 
          dX = wX - wX2
          dY = wY - wY2
          dZ = wZ - wZ2

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          rIJ = dX**2 + dY**2 + dZ**2
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
            endif

            enddo
           enddo
          enddo
          rIJ = sqrt(rIJ)
         else ! its an orthorhombic box
          
          xInt = ANINT((wX-wX2)/bX)
          yInt = ANINT((wY-wY2)/bY)
          zInt = ANINT((wZ-wZ2)/bZ)

          dX = (wX-wX2) - xInt*bX
          dY = (wY-wY2) - yInt*bY
          dZ = (wZ-wZ2) - zInt*bZ         
          rIJ = sqrt(dX*dX+dY*dY+dZ*dZ)

         endif ! End check for truncation

         cellBindE = cellBindE + wQ2/rIJ     
        ! enddo ! End loop over water atoms    
                
        endif ! End test for finding other waters
        
        enddo ! End loop over protein (really?) atoms
        
        cellBindE = cellBindE*wQ
      return
      end function cellBindE      

      function isFirstShell(wwX,wwY,wwZ)
c    *******************************************************************
c    ** determines if a water molecule (oxygen) is within the first   **
c    ** solvation shell of a protein.                                 **
c    **                                                               **
c    ** principal variables:                                          **
c    **          real    :: wwX, wwY, wwZ - x,y,z coords. of h2o      **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Function def
        logical :: isFirstShell

        ! Args
        real    :: wwX, wwY, wwZ
        
        ! Locals
        integer :: irc
        integer :: ix, iy, iz        
        integer :: xInt, yInt, zInt        
        real*8  :: dfrx, dfry, dfrz, delr2        
        real*8  :: prX, prY,prZ , prR   
        real*8  :: rIJ, dX, dY, dZ
        real*8  :: wX, wY, wZ, myCut
c    ******************************************************************* 
         isFirstShell = .false.
         
         wX = real(wwX,kind=8)
         wY = real(wwY,kind=8)
         wZ = real(wwZ,kind=8)
         
         do irc = 1, iProtEnd  ! Loop over protein atoms
          prX = real(nodeConf(irc)%x,kind=8)
          prY = real(nodeConf(irc)%y,kind=8)
          prZ = real(nodeConf(irc)%z,kind=8)
          prR = real(qrms(irc)%r,kind=8)
          myCut = rBindCut !2.87d0        
          !myCut = prR + rBindCut !2.87d0        
          !myCut = rBindCut

         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then 
          dX = wX - prX
          dY = wY - prY
          dZ = wZ - prZ

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          rIJ = dX**2 + dY**2 + dZ**2
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
            endif

            enddo
           enddo
          enddo
          rIJ = sqrt(rIJ)
         else ! its an orthorhombic box
         
          xInt = ANINT((wX-prX)/bX)
          yInt = ANINT((wY-prY)/bY)
          zInt = ANINT((wZ-prZ)/bZ)

          dX = (wX-prX) - xInt*bX
          dY = (wY-prY) - yInt*bY
          dZ = (wZ-prZ) - zInt*bZ         
          rIJ = sqrt(dX*dX+dY*dY+dZ*dZ)

         endif ! End check for truncation

         ! Find if the water is in the proteins 1st solvation shell         
         if(rIJ .le. myCut) then
          isFirstShell = .true.          
          exit         
         endif          
         enddo ! End loop over protein atoms       

         ! check that we havent already found the first shell water
         if(.not. isFirstShell) then
         do irc = iOthStart, iOthEnd ! loop over other atoms
          prX = real(nodeConf(irc)%x,kind=8)
          prY = real(nodeConf(irc)%y,kind=8)
          prZ = real(nodeConf(irc)%z,kind=8)
          prR = real(qrms(irc)%r,kind=8)
          myCut = rBindCut !2.87d0        
          !myCut = prR + rBindCut !2.87d0        
          !myCut = rBindCut

         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then 
          dX = wX - prX
          dY = wY - prY
          dZ = wZ - prZ

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          rIJ = dX**2 + dY**2 + dZ**2
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
            endif

            enddo
           enddo
          enddo
          rIJ = sqrt(rIJ)
         else ! its an orthorhombic box
         
          xInt = ANINT((wX-prX)/bX)
          yInt = ANINT((wY-prY)/bY)
          zInt = ANINT((wZ-prZ)/bZ)

          dX = (wX-prX) - xInt*bX
          dY = (wY-prY) - yInt*bY
          dZ = (wZ-prZ) - zInt*bZ         
          rIJ = sqrt(dX*dX+dY*dY+dZ*dZ)

         endif ! End check for truncation

         ! Find if the water is in the proteins 1st solvation shell         
         if(rIJ .le. myCut) then
          isFirstShell = .true.          
          exit         
         endif          
         enddo ! End loop over other atoms
         endif ! end check that we havent already found a first shell atom

      return
      end function isFirstShell
      
      subroutine calcEnerCoul()
c    *******************************************************************
c    ** calculates solvent (+ protein) energetics at every timestep   **
c    ** using a database of atoms that are involved in the Stokes     **
c    ** Shift correlation function.                                   **
c    **                                                               **
c    ** principal variables:                                          **
c    ** logical    doProt  (input) calculate protein energetics       **
c    **                                                               **
c    *******************************************************************            
        implicit none
       
        ! Locals
        integer :: iw, irc, jj, iWrite, iatm 
        real    :: ddX, ddY, ddZ, ddQ
        real    :: dV_i     ! the 32-bit version of dv0s_i
c    ******************************************************************* 
       doIons = .false.
        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then 
          open(unit=iEne,file="enerCoul.dat",status="unknown")               
         endif
        endif 
        
        dv0s = 0.d0
        dv0p = 0.d0
        dv0i = 0.d0
        dv0o = 0.d0
c---------------------------------------------------------c         
         do irc = 1, iNumDB      
c---------------------------------------------------------c
            iatm = dVosDbase(irc)%atom
            ddX = nodeConf(iatm)%x
            ddY = nodeConf(iatm)%y
            ddZ = nodeConf(iatm)%z
            ddQ = dVosDbase(irc)%q
            
          if (doSolv) then          
           if (useCUDA) then
#ifdef ENABLE_CUDA
             call dvos_gpu(ddX,ddY,ddZ,ddQ,
     x                     bX,bY,bZ,nodeConf,qrms,dV_i)
             dv0s_i = real(dV_i,kind=8) ! convert our value to the 64-bit version
#endif
           !  call mu0s_direct(dv0s_i,
           !x            ddX, ddY, ddZ, ddQ) ! MUST BE CHANGED FOR EWALD CALC  
           else
             call mu0s_direct(dv0s_i,
     x            ddX, ddY, ddZ, ddQ)                  
           endif ! end check for ewald           
            dv0s = dv0s + dv0s_i                  
          endif ! end test for solvent calculations

         if (0 .eq. 1) then
           call mu0i_direct(dv0i_i,
     x          ddX, ddY, ddZ, ddQ)
     
           dv0i = dv0i + dv0i_i                      
         endif ! end test for ion energy calculations
         
         if (doProt) then
           call mu0p_direct(dv0p_i,
     x          ddX, ddY, ddZ, ddQ)
     
           dv0p = dv0p + dv0p_i                 
         endif ! End test for protein calc        

         if (doOthers) then
           call mu0o_direct(dv0o_i,
     x          ddX, ddY, ddZ, ddQ)
     
           dv0o = dv0o + dv0o_i                 
         endif ! End test for others calc    
c---------------------------------------------------------c          
         enddo ! end loop over reaction center atoms
c---------------------------------------------------------c         

        if (doSolv) then         
         call MPI_GATHER(real(dv0s,kind=4),1,MPI_REAL,
     x                     dv0sTraj,1,MPI_REAL,
     x                     master,world,ierr)     
        endif ! End test for solvent calculation

        if (0 .eq. 1) then                
          call MPI_GATHER(real(dv0i,kind=4),1,MPI_REAL,
     x                     dv0iTraj,1,MPI_REAL,
     x                     master,world,ierr)      
        endif ! End test for protein calc  
        
        if (doProt) then                
          call MPI_GATHER(real(dv0p,kind=4),1,MPI_REAL,
     x                     dv0pTraj,1,MPI_REAL,
     x                     master,world,ierr)      
        endif ! End test for protein calc           

        if (doOthers) then                
          call MPI_GATHER(real(dv0o,kind=4),1,MPI_REAL,
     x                     dv0oTraj,1,MPI_REAL,
     x                     master,world,ierr)      
        endif ! End test for protein calc           

           
          if(rank .eq. master)  then
           
            !*** NEED TO MAKE THIS MORE GENERAL, NO IF STATEMENTS ***
            do iWrite = 1, nProc
              if (doProt .and. doSolv .and. doIons) then
               write(iEne, *) iConf+iWrite-1, 
     x         dv0sTraj(iWrite), dv0pTraj(iWrite), 
     x         dv0iTraj(iWrite)
              elseif(doSolv .and. doProt 
     x              .and. doOthers) then
               write(iEne, *) iConf+iWrite-1, 
     x         dv0sTraj(iWrite), dv0pTraj(iWrite), 
     x         dv0oTraj(iWrite)
              elseif(doSolv .and. doProt) then
               write(iEne, *) iConf+iWrite-1, 
     x         dv0sTraj(iWrite), dv0pTraj(iWrite)
              elseif(doSolv .and. doIons) then
               write(iEne, *) iConf+iWrite-1, 
     x         dv0sTraj(iWrite), dv0iTraj(iWrite)
              elseif(doIons .and. doProt) then
               write(iEne, *) iConf+iWrite-1, 
     x         dv0pTraj(iWrite), dv0iTraj(iWrite)
              elseif(doIons) then
               write(iEne, *) iConf+iWrite-1, 
     x         dv0iTraj(iWrite)     
              elseif(doProt) then
               write(iEne, *) iConf+iWrite-1, 
     x         dv0pTraj(iWrite)     
              else ! solvent only
               write(iEne, *) iConf+iWrite-1, 
     x         dv0sTraj(iWrite)
              endif ! end check for doProt/doSolv
            enddo
           endif ! End check for master       

        return
      end subroutine calcEnerCoul  

      subroutine calcEnerInd()
c    *******************************************************************
c    ** calculates induction energetics at every timestep             **
c    ** using a database of atoms that are involved in the Stokes     **
c    ** Shift correlation function and all waters.                    **
c    **                                                               **
c    *******************************************************************            
        implicit none
        
        ! Locals
        integer :: iw, irc, iWrite
        integer :: iAtm, iiP, iIon
        integer :: ix, iy, iz        
        integer, parameter :: XX=1
        integer, parameter :: YY=2
        integer, parameter :: ZZ=3
        integer :: xInt, yInt, zInt
        real    :: dX, dY, dZ         
        real*8  :: Efield1(3),Efield2(3)
        real*8  :: dQ, rX, rY, rZ, rIJ
        real*8  :: pQ, pX, pY, pZ
        real*8  :: lowDX, lowDY, lowDZ
        real*8  :: ionQ, ionX, ionY, ionZ        
        real*8  :: wQ, wX, wY, wZ        
        real*8  :: rcX, rcY, rcZ, rcQ1, rcQ2
        real*8  :: rinv, rinv3, E1sq,E2sq        
        real*8  :: dfrx, dfry, dfrz, delr2
        real*8  :: alphaP, alphaI, alphaW 
        logical :: skipAtom
        real    :: myQ
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then 
          open(unit=iInd,file="enerInd.dat",status="unknown")               
         endif
        endif 
        
        dv0in = 0.d0         
        dv0inS = 0.d0
        dv0inP = 0.d0
        
c -----------------------------------------------------------------
c ---------------  Non-RC protein Atom Section --------------------
c -----------------------------------------------------------------
        if(doProt) then
        do iip = 1, iProtEnd ! Loop over all non-RC, protein atoms
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        
        
         !pQ = real(qrms(iip)%q,kind=8)
         pX = real(nodeConf(iip)%x,kind=8)                     
         pY = real(nodeConf(iip)%y,kind=8)
         pZ = real(nodeConf(iip)%z,kind=8)
         alphaP = alphaDbase(iip)
         
         do irc = 1, iNumDB  ! Loop over all RC atoms     
          iAtm = dVosDbase(irc)%atom          
          call searchDBAtom(iip,myQ,skipAtom)
          if (.not. skipAtom) then
          
          !if (iip .ne. iAtm) then ! Make sure we are on different atoms          
           rcQ1 = real(qrms(iAtm)%q,kind=8)
           rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)
           rcX = real(nodeConf(iAtm)%x,kind=8)                     
           rcY = real(nodeConf(iAtm)%y,kind=8)
           rcZ = real(nodeConf(iAtm)%z,kind=8)

           dX = pX - rcX
           dY = pY - rcY
           dZ = pZ - rcZ
          
           rinv = 1.d0/sqrt(dX*dX + dY*dY + dZ*dZ)
           rinv3 = rinv**3
           
           Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
           Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
           Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
           Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
           Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
           Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          
     
          endif
          
         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         !dv0in = dv0in - 0.5d0*alphaP*(E2sq-E1sq)     
        dv0inP = dv0inP - 0.5d0*alphaP*(E2sq-E1sq) 
        enddo ! End the loop over non-RC, protein atoms
        
         call MPI_GATHER(real(dv0inP,kind=4),1,MPI_REAL,
     x                     dv0inPTraj,1,MPI_REAL,
     x                     master,world,ierr)     
        
       endif ! End check for doProt
c -----------------------------------------------------------------
c ---------------  Counter Ion Section ----------------------------
c -----------------------------------------------------------------
       if(0 .eq. 1) then ! always off, for now...
c -----------------------------------------------------------------       
        do iIon = iProtEnd+1, iWatStart-1 ! Loop over all counter ions
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        
        
         ionQ = real(qrms(iIon)%q,kind=8)
         ionX = real(nodeConf(iIon)%x,kind=8)                     
         ionY = real(nodeConf(iIon)%y,kind=8)
         ionZ = real(nodeConf(iIon)%z,kind=8)
         alphaI = alphaDbase(iIon)
         
         do irc = 1, iNumDB  ! Loop over RC atoms
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)
          rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)
          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then 
          dX = ionX - rcX
          dY = ionY - rcY
          dZ = ionZ - rcZ

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          rIJ = 1000000.0
          !rIJ = dX**2 + dY**2 + dZ**2
          
          lowDX = dX
          lowDY = dY              
          lowDZ = dZ   
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
              lowDX = dX
              lowDY = dY              
              lowDZ = dZ              
            endif

            enddo
           enddo
          enddo
          
           dX = lowDX
           dY = lowDY
           dZ = lowDZ

          else  ! its orthorhombic

           xInt = ANINT((ionX-rcX)/bX)
           yInt = ANINT((ionY-rcY)/bY)
           zInt = ANINT((ionZ-rcZ)/bZ)

           dX = (ionX-rcX) - xInt*bX
           dY = (ionY-rcY) - yInt*bY
           dZ = (ionZ-rcZ) - zInt*bZ         
           rIJ = sqrt(dX*dX+dY*dY+dZ*dZ)
          
          endif ! End test for truncated box           
          
           rinv = 1.d0/sqrt(rIJ)
           rinv3 = rinv**3
           
           Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
           Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
           Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
           Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
           Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
           Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          
     
          
         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
        dv0in = dv0in - 0.5d0*alphaI*(E2sq-E1sq)               
        enddo ! End the loop over counter ions   
c -----------------------------------------------------------------          
       endif ! end always off
c -----------------------------------------------------------------       
        
c -----------------------------------------------------------------
c ---------------  Water Atom Section -----------------------------
c -----------------------------------------------------------------
         if(doSolv) then
         do iw = iWatStart, natom        
          Efield1(:) = 0.d0
          Efield2(:) = 0.d0         
          wX = real(nodeConf(iw)%x,kind=8)                     
          wY = real(nodeConf(iw)%y,kind=8)
          wZ = real(nodeConf(iw)%z,kind=8)          
          alphaW = alphaDbase(iw)
         
         do irc = 1, iNumDB  ! Loop over all counter ions   
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)
          rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)
          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then 
          dX = wX - rcX
          dY = wY - rcY
          dZ = wZ - rcZ

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          rIJ = 1000000.0
         
          !rIJ = dX**2 + dY**2 + dZ**2
          
          !rIJ = 1000000.0
          lowDX = dX
          lowDY = dY              
          lowDZ = dZ   
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
              lowDX = dX
              lowDY = dY              
              lowDZ = dZ              
            endif

            enddo
           enddo
          enddo
          
           dX = lowDX
           dY = lowDY
           dZ = lowDZ
           
         else ! its an orthorhombic box
         
          xInt = ANINT((wX-rcX)/bX)
          yInt = ANINT((wY-rcY)/bY)
          zInt = ANINT((wZ-rcZ)/bZ)

          dX = (wX-rcX) - xInt*bX
          dY = (wY-rcY) - yInt*bY
          dZ = (wZ-rcZ) - zInt*bZ         
          rIJ = sqrt(dX*dX+dY*dY+dZ*dZ)

         endif
         
           rinv = 1.d0/sqrt(rIJ)
           rinv3 = rinv**3
           
           Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
           Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
           Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
           Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
           Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
           Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          
     
          
         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         !dv0in = dv0in - 0.5d0*alphaW*(E2sq-E1sq)     
         dv0inS = dv0inS - 0.5d0*alphaW*(E2sq-E1sq)
        enddo ! End the loop over all waters     
        
         call MPI_GATHER(real(dv0inS,kind=4),1,MPI_REAL,
     x                     dv0inSTraj,1,MPI_REAL,
     x                     master,world,ierr)     
        
        endif ! Check for doSolv
c -----------------------------------------------------------------
c ---------------  End Calculations  ------------------------------
c -----------------------------------------------------------------
         
c         call MPI_GATHER(real(dv0in,kind=4),1,MPI_REAL,
c     x                     dv0inTraj,1,MPI_REAL,
c     x                     master,world,ierr)     
     
          
          if(rank .eq. master)  then
           
            do iWrite = 1, nProc
c               write(iInd, *) iConf+iWrite-1, 
c     x         dv0inTraj(iWrite)
              if (doProt .and. doSolv) then
               write(iInd, *) iConf+iWrite-1, 
     x         dv0inSTraj(iWrite), dv0inPTraj(iWrite)
              elseif(doProt) then
               write(iInd, *) iConf+iWrite-1, 
     x         dv0inPTraj(iWrite)
              else ! Just solvent
               write(iInd, *) iConf+iWrite-1, 
     x         dv0inSTraj(iWrite)              
              endif !end check for solv/prot

            enddo
           endif ! End check for master       
        
        return
      end subroutine calcEnerInd        

      subroutine calcEnerExclr()
c    *******************************************************************
c    ** calculates induction + coulomb energetics at every timestep   **
c    ** using a database of atoms that are involved in the Stokes     **
c    ** Shift correlation function and all other atoms.               **
c    **                                                               **
c    ** This function is identical to calcEnerAll, but includes       **
c    ** per-atom exclusions. i.e., if the atom pair distance is not   **
c    ** BEYOND a given exclusion distance, the atom pair energy is    **
c    ** excluded.                                                     **
c    **                                                               **
c    *******************************************************************
        implicit none
        
        ! Locals
        integer :: iw, irc, iWrite
        integer :: iAtm, iiP, iIon
        integer :: iLip, iOth
        integer :: ix, iy, iz        
        integer :: pidx, widx, lidx
        integer :: oidx, iidx
        integer, parameter :: XX=1
        integer, parameter :: YY=2
        integer, parameter :: ZZ=3
        integer :: xInt, yInt, zInt
        real    :: dXi, dYi, dZi
        real    :: dX, dY, dZ         
        real*8  :: Efield1(3),Efield2(3)
        real*8  :: dQ, rX, rY, rZ, rIJ
        real*8  :: pQ, pX, pY, pZ
        real*8  :: lowDX, lowDY, lowDZ
        real*8  :: ionQ, ionX, ionY, ionZ        
        real*8  :: lipQ, lipX, lipY, lipZ
        real*8  :: othQ, othX, othY, othZ        
        real*8  :: wQ, wX, wY, wZ, delQ        
        real*8  :: rcX, rcY, rcZ, rcQ1, rcQ2
        real*8  :: rinv, rinv3, E1sq,E2sq        
        real*8  :: dfrx, dfry, dfrz, delr2
        real*8  :: rExcl
        real*8  :: alphaP, alphaI, alphaW 
        logical :: skipAtom
        real    :: myQ
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then 
          open(unit=iEnerAll,file="enerAll.dat",status="unknown")

          ! write the file header
          write(iEnerAll, '(a)',advance="no") "# frame "
          if (doProt) then
            write(iEnerAll, '(a)',advance="no") 
     x      "Prot(C)  Prot(Ind)  "
          endif 
          if (doSolv) then
            write(iEnerAll, '(a)',advance="no") 
     x      "Wat(C)  Wat(Ind)  "
          endif 
          if (doIons) then
            write(iEnerAll, '(a)',advance="no") 
     x      "Ion(C)  Ion(Ind)  "
          endif 
          if (doLipid) then
            write(iEnerAll, '(a)',advance="no") 
     x       "Lip(C)  Lip(Ind)  "
          endif 
          if (doOthers) then
            write(iEnerAll, '(a)',advance="no") 
     x        "Oth(C)  Oth(Ind)  "
          endif 
          write(iEnerAll,'(a)',advance="yes") ""



         endif
        endif 
        
        dv0in = 0.d0         
        dv0s = 0.d0   ! solvent coulomb
        dv0p = 0.d0   ! protein Coulomb
        dv0l = 0.d0   ! lipid Coulomb
        dv0i = 0.d0   ! ions Coulomb
        dv0o = 0.d0   ! Other atoms Coulomb
        dv0inS = 0.d0 ! Solvent (water) induction
        dv0inP = 0.d0 ! Protein induction
        dv0inI = 0.d0 ! Ion induction
        dv0inL = 0.d0 ! Lipid induction
        dv0inO = 0.d0 ! Other atoms induction

c -----------------------------------------------------------------
c ---------------  Non-RC protein Atom Section --------------------
c -----------------------------------------------------------------    
       if (doProt) then 

        do pidx = 1, nProtStEnd
         iProtStart = arrProtStart(pidx)
         iProtEnd   = arrProtEnd(pidx)

        do iip = iProtStart, iProtEnd ! Loop over all non-RC, protein atoms -- new version
        !do iip = 1, iProtEnd ! Loop over all non-RC, protein atoms -- OLD Version
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        

         rExcl = exclrDist(iip) ! get the exclusion distance for the protein atom

         pQ = real(qrms(iip)%q,kind=8)
         pX = real(nodeConf(iip)%x,kind=8)                     
         pY = real(nodeConf(iip)%y,kind=8)
         pZ = real(nodeConf(iip)%z,kind=8)
         alphaP = alphaDbase(iip)
         
         do irc = 1, iNumDB  ! Loop over all RC atoms     
          iAtm = dVosDbase(irc)%atom          
          call searchDBAtom(iip,myQ,skipAtom)
          if (.not. skipAtom) then
          
          !if (iip .ne. iAtm) then ! Make sure we are on different atoms          
           rcQ1 = real(qrms(iAtm)%q,kind=8)
           delQ = real(dVosDbase(irc)%q,kind=8)
           !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)
           rcQ2 = rcQ1 + delQ
           rcX = real(nodeConf(iAtm)%x,kind=8)                     
           rcY = real(nodeConf(iAtm)%y,kind=8)
           rcZ = real(nodeConf(iAtm)%z,kind=8)


           ! coordinates are NOT reimaged for the protein
           ! because interactions with the proteins images
           ! are not physically relavant
           dX = pX - rcX
           dY = pY - rcY
           dZ = pZ - rcZ

           rIJ = sqrt(dX*dX + dY*dY + dZ*dZ)

           ! Only add if not within exclusion distance
           if(rIJ > rExcl) then

            rinv = 1.d0/rIJ
            rinv3 = rinv**3
           
            Efield1(XX) = Efield1(XX) +
     x                 dX*rcQ1*rinv3
            Efield1(YY) = Efield1(YY) +
     x                 dY*rcQ1*rinv3
            Efield1(ZZ) = Efield1(ZZ) +
     x                 dZ*rcQ1*rinv3          
     
            Efield2(XX) = Efield2(XX) +
     x                 dX*rcQ2*rinv3
            Efield2(YY) = Efield2(YY) +
     x                 dY*rcQ2*rinv3
            Efield2(ZZ) = Efield2(ZZ) +
     x                 dZ*rcQ2*rinv3          

           dv0p = dv0p + pQ*delQ*rinv ! Calculate the columb(prot) part
           endif ! end check for exclusion
          endif ! end check for skip Atom
          
         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inP = dv0inP - 0.5d0*alphaP*(E2sq-E1sq)     
          
        enddo ! End the loop over non-RC, protein atoms
        enddo ! End the loop over the protein starts/ends
         
         call MPI_GATHER(real(dv0inP,kind=4),1,MPI_REAL,
     x                     dv0inPTraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0p,kind=4),1,MPI_REAL,
     x                     dv0pTraj,1,MPI_REAL,
     x                     master,world,ierr) 

       endif ! End check for doProt 
        
c -----------------------------------------------------------------
c ---------------  Counter Ion Section ----------------------------
c -----------------------------------------------------------------

       if (doIons) then 
c -----------------------------------------------------------------   
        do iidx = 1, nIonStEnd
         iIonStart = arrIonStart(iidx)
         iIonEnd   = arrIonEnd(iidx)

        do iIon = iIonStart, iIonEnd ! Loop over all counter ions
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        

         rExcl = exclrDist(iIon) ! get the exclusion distance for the ion

         ionQ = real(qrms(iIon)%q,kind=8)
         ionX = real(nodeConf(iIon)%x,kind=8)                     
         ionY = real(nodeConf(iIon)%y,kind=8)
         ionZ = real(nodeConf(iIon)%z,kind=8)
         alphaI = alphaDbase(iIon)

         
         do irc = 1, iNumDB  ! Loop over RC atoms
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)

          delQ =  real(dVosDbase(irc)%q,kind=8) ! new
          !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)          
          rcQ2 = rcQ1 + delQ  !new

          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

          dXi = ionX - rcX
          dYi = ionY - rcY
          dZi = ionZ - rcZ

          call imageCoords(dXi,dYi,dZi,dX,dY,dZ)

          rIJ = sqrt(dX*dX + dY*dY + dZ*dZ)

          ! Only add if not within exclusion distance
          if(rIJ > rExcl) then

          
           rinv = 1.d0/rIJ
           rinv3 = rinv**3
           
           Efield1(XX) = Efield1(XX) +
     x                 dX*rcQ1*rinv3
           Efield1(YY) = Efield1(YY) +
     x                 dY*rcQ1*rinv3
           Efield1(ZZ) = Efield1(ZZ) +
     x                 dZ*rcQ1*rinv3          
     
           Efield2(XX) = Efield2(XX) +
     x                 dX*rcQ2*rinv3
           Efield2(YY) = Efield2(YY) +
     x                 dY*rcQ2*rinv3
           Efield2(ZZ) = Efield2(ZZ) +
     x                 dZ*rcQ2*rinv3          
     
          dv0i = dv0i +
     x          ionQ*delQ*rinv ! Calculate the columb(ion) part
          endif ! End check for exclusion

         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inI = dv0inI - 
     x            0.5d0*alphaI*(E2sq-E1sq)     
          
        enddo ! End the loop over counter ions  
        enddo ! End the loop over the protein starts/ends 
c -----------------------------------------------------------------

         call MPI_GATHER(real(dv0inI,kind=4),1,MPI_REAL,
     x                     dv0inITraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0i,kind=4),1,MPI_REAL,
     x                     dv0iTraj,1,MPI_REAL,
     x                     master,world,ierr) 

       endif ! end check for do ions
c -----------------------------------------------------------------       


c ###############################################################
c -----------------------------------------------------------------
c ---------------  Lipid Section ----------------------------
c -----------------------------------------------------------------

       if (doLipid) then 
c -----------------------------------------------------------------
        do lidx = 1, nLipStEnd
         iLipStart = arrLipStart(lidx)
         iLipEnd   = arrLipEnd(lidx)

        do iLip = iLipStart, iLipEnd ! Loop over all counter lipids
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        

         rExcl = exclrDist(iLip) ! get the exclusion distance for the lipid atom

         lipQ = real(qrms(iLip)%q,kind=8)
         lipX = real(nodeConf(iLip)%x,kind=8)                     
         lipY = real(nodeConf(iLip)%y,kind=8)
         lipZ = real(nodeConf(iLip)%z,kind=8)
         alphaI = alphaDbase(iLip)
         
         do irc = 1, iNumDB  ! Loop over RC atoms
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)

          delQ =  real(dVosDbase(irc)%q,kind=8) ! new
          !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)          
          rcQ2 = rcQ1 + delQ  !new


          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

          dXi = lipX - rcX
          dYi = lipY - rcY
          dZi = lipZ - rcZ

          call imageCoords(dXi,dYi,dZi,dX,dY,dZ)

          rIJ = sqrt(dX*dX + dY*dY + dZ*dZ)

          ! Only add if not within exclusion distance
          if(rIJ > rExcl) then

           rinv = 1.d0/rIJ
           rinv3 = rinv**3
           
           Efield1(XX) = Efield1(XX) +
     x                 dX*rcQ1*rinv3
           Efield1(YY) = Efield1(YY) +
     x                 dY*rcQ1*rinv3
           Efield1(ZZ) = Efield1(ZZ) +
     x                 dZ*rcQ1*rinv3          
     
           Efield2(XX) = Efield2(XX) +
     x                 dX*rcQ2*rinv3
           Efield2(YY) = Efield2(YY) +
     x                 dY*rcQ2*rinv3
           Efield2(ZZ) = Efield2(ZZ) +
     x                 dZ*rcQ2*rinv3          
     
          dv0l = dv0l +
     x          lipQ*delQ*rinv ! Calculate the columb(ion) part
          endif ! end check for exclusion atom

         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inL = dv0inL - 0.5d0*alphaI*(E2sq-E1sq)     
          
        enddo ! End the loop over lipids
        enddo ! End the loop over the protein starts/ends
c -----------------------------------------------------------------
c         if(iConf .eq. 1) then

c          print *, "dv0l= ", dv0l
c          print *, "lipQ= ", lipQ
c         endif

         call MPI_GATHER(real(dv0inL,kind=4),1,MPI_REAL,
     x                     dv0inLTraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0l,kind=4),1,MPI_REAL,
     x                     dv0lTraj,1,MPI_REAL,
     x                     master,world,ierr) 


       endif ! end check for do lipids
c ----------------------------------------------------------------- 

c -----------------------------------------------------------------
c ---------------  Other Atom Section ----------------------------
c -----------------------------------------------------------------

       if (doOthers) then 
c -----------------------------------------------------------------       
        do oidx = 1, nOthStEnd
         iOthStart = arrOthStart(oidx)
         iOthEnd   = arrOthEnd(oidx)

        do iOth = iOthStart, iOthEnd ! Loop over all counter lipids
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        

         rExcl = exclrDist(iOth) ! get the exclusion distance for the other atom

         othQ = real(qrms(iOth)%q,kind=8)
         othX = real(nodeConf(iOth)%x,kind=8)                     
         othY = real(nodeConf(iOth)%y,kind=8)
         othZ = real(nodeConf(iOth)%z,kind=8)
         alphaI = alphaDbase(iOth)
         
         do irc = 1, iNumDB  ! Loop over RC atoms
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)

          delQ =  real(dVosDbase(irc)%q,kind=8) ! new
          !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)          
          rcQ2 = rcQ1 + delQ  !new


          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

          dXi = othX - rcX
          dYi = othY - rcY
          dZi = othZ - rcZ

          call imageCoords(dXi,dYi,dZi,dX,dY,dZ)

          rIJ = sqrt(dX*dX + dY*dY + dZ*dZ)

          ! Only add if not within exclusion distance
          if(rIJ > rExcl) then

           rinv = 1.d0/rIJ
           rinv3 = rinv**3
           
           Efield1(XX) = Efield1(XX) +
     x                 dX*rcQ1*rinv3
           Efield1(YY) = Efield1(YY) +
     x                 dY*rcQ1*rinv3
           Efield1(ZZ) = Efield1(ZZ) +
     x                 dZ*rcQ1*rinv3          
     
           Efield2(XX) = Efield2(XX) +
     x                 dX*rcQ2*rinv3
           Efield2(YY) = Efield2(YY) +
     x                 dY*rcQ2*rinv3
           Efield2(ZZ) = Efield2(ZZ) +
     x                 dZ*rcQ2*rinv3          
     
          dv0o = dv0o +
     x          othQ*delQ*rinv ! Calculate the columb(others) part
          endif ! End check for exclusion atom

         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inO = dv0inO - 0.5d0*alphaI*(E2sq-E1sq)     
          
        enddo ! End the loop over other atoms
        enddo ! End the loop over the protein starts/ends
c -----------------------------------------------------------------

         call MPI_GATHER(real(dv0inO,kind=4),1,MPI_REAL,
     x                     dv0inOTraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0o,kind=4),1,MPI_REAL,
     x                     dv0oTraj,1,MPI_REAL,
     x                     master,world,ierr) 

       endif ! end check for do other atoms
c ----------------------------------------------------------------- 

  
c ###############################################################        
c -----------------------------------------------------------------
c ---------------  Water Atom Section -----------------------------
c -----------------------------------------------------------------    
        if(doSolv) then 
        do widx = 1, nWatStEnd
         iWatStart = arrWatStart(widx)
         iWatEnd   = arrWatEnd(widx)

         do iw = iWatStart, iWatEnd
         !do iw = iWatStart, natom        
          Efield1(:) = 0.d0
          Efield2(:) = 0.d0    

          rExcl = exclrDist(iw) ! get the exclusion distance for the water atom

          wQ = real(qrms(iw)%q,kind=8)
          wX = real(nodeConf(iw)%x,kind=8)                     
          wY = real(nodeConf(iw)%y,kind=8)
          wZ = real(nodeConf(iw)%z,kind=8)          
          alphaW = alphaDbase(iw)


         do irc = 1, iNumDB  ! Loop over all database reaction center atoms
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)
          delQ =  real(dVosDbase(irc)%q,kind=8)
          !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)          
          rcQ2 = rcQ1 + delQ
          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

          dXi = wX - rcX
          dYi = wY - rcY
          dZi = wZ - rcZ

          call imageCoords(dXi,dYi,dZi,dX,dY,dZ)

          rIJ = sqrt(dX*dX + dY*dY + dZ*dZ)

          ! Only add if not within exclusion distance
          if(rIJ > rExcl) then

          rinv = 1.d0/rIJ
          rinv3 = rinv**3
           
          Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
          Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
          Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
          Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
          Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
          Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          
     
         dv0s = dv0s + wQ*delQ*rinv ! Calculate the columb(water) part        
         endif ! end check for exclusion atom

        enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inS = dv0inS - 0.5d0*alphaW*(E2sq-E1sq)     
          
        enddo ! End the loop over all waters     
       enddo ! End loop over water starts/ends

         call MPI_GATHER(real(dv0inS,kind=4),1,MPI_REAL,
     x                     dv0inSTraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0s,kind=4),1,MPI_REAL,
     x                     dv0sTraj,1,MPI_REAL,
     x                     master,world,ierr) 
     
       endif ! End check for doSolv        
c -----------------------------------------------------------------
c ---------------  End Calculations  ------------------------------
c -----------------------------------------------------------------                        
          if(rank .eq. master)  then
           
            do iWrite = 1, nProc

          write(iEnerAll, '(i6)',advance="no") 
     x    iConf+iWrite-1
          write(iEnerAll, '(a)',advance="no") " "

          if (doProt) then
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0pTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0inPTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
          endif 

          if (doSolv) then
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0sTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0inSTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
          endif 

          if (doIons) then
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0iTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0inITraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
          endif 

          if (doLipid) then
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0lTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0inLTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
          endif 

          if (doOthers) then
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0oTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0inOTraj(iWrite)
          endif 

          write(iEnerAll,'(a)',advance="yes") ""
            
            enddo ! end loop over iwrite
           endif ! End check for master       
        
        return
      end subroutine calcEnerExclr

      subroutine calcEnerAll()
c    *******************************************************************
c    ** calculates induction + coulomb energetics at every timestep   **
c    ** using a database of atoms that are involved in the Stokes     **
c    ** Shift correlation function and all other atoms.               **
c    **                                                               **
c    *******************************************************************            
        implicit none
        
        ! Locals
        integer :: iw, irc, iWrite
        integer :: iAtm, iiP, iIon
        integer :: iLip, iOth
        integer :: ix, iy, iz        
        integer :: pidx, widx, lidx
        integer :: oidx, iidx
        integer, parameter :: XX=1
        integer, parameter :: YY=2
        integer, parameter :: ZZ=3
        integer :: xInt, yInt, zInt
        real    :: dXi, dYi, dZi
        real    :: dX, dY, dZ         
        real*8  :: Efield1(3),Efield2(3)
        real*8  :: dQ, rX, rY, rZ, rIJ
        real*8  :: pQ, pX, pY, pZ
        real*8  :: lowDX, lowDY, lowDZ
        real*8  :: ionQ, ionX, ionY, ionZ        
        real*8  :: lipQ, lipX, lipY, lipZ
        real*8  :: othQ, othX, othY, othZ        
        real*8  :: wQ, wX, wY, wZ, delQ        
        real*8  :: rcX, rcY, rcZ, rcQ1, rcQ2
        real*8  :: rinv, rinv3, E1sq,E2sq        
        real*8  :: dfrx, dfry, dfrz, delr2
        real*8  :: alphaP, alphaI, alphaW 
        logical :: skipAtom
        real    :: myQ
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then 
c          open(unit=iInd,file="enerInd.dat",status="unknown") 
c          open(unit=iEne,file="enerCoul.dat",status="unknown")                         
          open(unit=iEnerAll,file="enerAll.dat",status="unknown")

          ! write the file header
          write(iEnerAll, '(a)',advance="no") "# frame "
          if (doProt) then
            write(iEnerAll, '(a)',advance="no") 
     x      "Prot(C)  Prot(Ind)  "
          endif 
          if (doSolv) then
            write(iEnerAll, '(a)',advance="no") 
     x      "Wat(C)  Wat(Ind)  "
          endif 
          if (doIons) then
            write(iEnerAll, '(a)',advance="no") 
     x      "Ion(C)  Ion(Ind)  "
          endif 
          if (doLipid) then
            write(iEnerAll, '(a)',advance="no") 
     x       "Lip(C)  Lip(Ind)  "
          endif 
          if (doOthers) then
            write(iEnerAll, '(a)',advance="no") 
     x        "Oth(C)  Oth(Ind)  "
          endif 
          write(iEnerAll,'(a)',advance="yes") ""



         endif
        endif 
        
        dv0in = 0.d0         
        dv0s = 0.d0   ! solvent coulomb
        dv0p = 0.d0   ! protein Coulomb
        dv0l = 0.d0   ! lipid Coulomb
        dv0i = 0.d0   ! ions Coulomb
        dv0o = 0.d0   ! Other atoms Coulomb
        dv0inS = 0.d0 ! Solvent (water) induction
        dv0inP = 0.d0 ! Protein induction
        dv0inI = 0.d0 ! Ion induction
        dv0inL = 0.d0 ! Lipid induction
        dv0inO = 0.d0 ! Other atoms induction

c        integer            :: iWatStart         ! Waters starting atomic index
c        integer            :: iWatEnd           ! Waters ending atomic index
c        integer            :: iProtStart        ! The proteins starting atomic index
c        integer            :: iProtEnd          ! The proteins ending atomic index
c        integer            :: iIonStart         ! The ions starting atomic index
c        integer            :: iIonEnd           ! The ions ending atomic index
c        integer            :: iOtherStart       ! Other atoms starting atomic index
c        integer            :: iOtherEnd         ! Other atoms ending atomic index
                                                !        
c    	integer, parameter :: maxSect = 10      ! maximum number of starts/ends
c        integer,dimension(maxSect) :: arrWatStart  ! array of starting indices for water
c        integer,dimension(maxSect) :: arrWatEnd    ! array of ending indices for water
c        integer,dimension(maxSect) :: arrProtStart ! array of starting indices for the protein
c        integer,dimension(maxSect) :: arrProtEnd   ! array of ending indices for the protein
c        integer,dimension(maxSect) :: arrIonStart  ! array of starting indices for the ions
c        integer,dimension(maxSect) :: arrIonEnd    ! array of ending indices for the ions
c        integer,dimension(maxSect) :: arrLipStart  ! array of starting indices for the lipids
c        integer,dimension(maxSect) :: arrLipEnd    ! array of ending indices for the lipids
c        integer,dimension(maxSect) :: arrOthStart  ! array of starting indices for all other atoms
c        integer,dimension(maxSect) :: arrOthEnd    ! array of ending indices for all other atoms
c        integer            ::  nWatStEnd        ! number of section startings/endings for water
c        integer            ::  nProtStEnd       ! number of section startings/endings for the proteins
c        integer            ::  nIonStEnd        ! number of section startings/endings for ions
c        integer            ::  nLipStEnd        ! number of section startings/endings for lipids
c        integer            ::  nOthStEnd        ! number of section startings/endings for other atoms        
c -----------------------------------------------------------------
c ---------------  Non-RC protein Atom Section --------------------
c -----------------------------------------------------------------    
       if (doProt) then 

        do pidx = 1, nProtStEnd
         iProtStart = arrProtStart(pidx)
         iProtEnd   = arrProtEnd(pidx)

        do iip = iProtStart, iProtEnd ! Loop over all non-RC, protein atoms -- new version
        !do iip = 1, iProtEnd ! Loop over all non-RC, protein atoms -- OLD Version
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        
        
         pQ = real(qrms(iip)%q,kind=8)
         pX = real(nodeConf(iip)%x,kind=8)                     
         pY = real(nodeConf(iip)%y,kind=8)
         pZ = real(nodeConf(iip)%z,kind=8)
         alphaP = alphaDbase(iip)
         
         do irc = 1, iNumDB  ! Loop over all RC atoms     
          iAtm = dVosDbase(irc)%atom          
          call searchDBAtom(iip,myQ,skipAtom)
          if (.not. skipAtom) then
          
          !if (iip .ne. iAtm) then ! Make sure we are on different atoms          
           rcQ1 = real(qrms(iAtm)%q,kind=8)
           delQ = real(dVosDbase(irc)%q,kind=8)
           !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)
           rcQ2 = rcQ1 + delQ
           rcX = real(nodeConf(iAtm)%x,kind=8)                     
           rcY = real(nodeConf(iAtm)%y,kind=8)
           rcZ = real(nodeConf(iAtm)%z,kind=8)


           ! coordinates are NOT reimaged for the protein
           ! because interactions with the proteins images
           ! are not physically relavant
           dX = pX - rcX
           dY = pY - rcY
           dZ = pZ - rcZ
          
           rinv = 1.d0/sqrt(dX*dX + dY*dY + dZ*dZ)
           rinv3 = rinv**3
           
           Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
           Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
           Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
           Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
           Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
           Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          

          dv0p = dv0p + pQ*delQ*rinv ! Calculate the columb(prot) part
          
          endif
          
         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inP = dv0inP - 0.5d0*alphaP*(E2sq-E1sq)     
          
        enddo ! End the loop over non-RC, protein atoms
        enddo ! End the loop over the protein starts/ends
         
         call MPI_GATHER(real(dv0inP,kind=4),1,MPI_REAL,
     x                     dv0inPTraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0p,kind=4),1,MPI_REAL,
     x                     dv0pTraj,1,MPI_REAL,
     x                     master,world,ierr) 

       endif ! End check for doProt 
        
c -----------------------------------------------------------------
c ---------------  Counter Ion Section ----------------------------
c -----------------------------------------------------------------

       if (doIons) then 
c -----------------------------------------------------------------   
        do iidx = 1, nIonStEnd
         iIonStart = arrIonStart(iidx)
         iIonEnd   = arrIonEnd(iidx)

        do iIon = iIonStart, iIonEnd ! Loop over all counter ions
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        
        
         ionQ = real(qrms(iIon)%q,kind=8)
         ionX = real(nodeConf(iIon)%x,kind=8)                     
         ionY = real(nodeConf(iIon)%y,kind=8)
         ionZ = real(nodeConf(iIon)%z,kind=8)
         alphaI = alphaDbase(iIon)

         
         do irc = 1, iNumDB  ! Loop over RC atoms
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)

          delQ =  real(dVosDbase(irc)%q,kind=8) ! new
          !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)          
          rcQ2 = rcQ1 + delQ  !new

          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

          dXi = ionX - rcX
          dYi = ionY - rcY
          dZi = ionZ - rcZ

          call imageCoords(dXi,dYi,dZi,dX,dY,dZ)
          rIJ = dX*dX+dY*dY+dZ*dZ
          
          rinv = 1.d0/sqrt(rIJ)
          rinv3 = rinv**3
           
          Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
          Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
          Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
          Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
          Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
          Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          
     
         dv0i = dv0i + 
     x          ionQ*delQ*rinv ! Calculate the columb(ion) part
          
         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inI = dv0inI - 
     x            0.5d0*alphaI*(E2sq-E1sq)     
          
        enddo ! End the loop over counter ions  
        enddo ! End the loop over the protein starts/ends 
c -----------------------------------------------------------------

         call MPI_GATHER(real(dv0inI,kind=4),1,MPI_REAL,
     x                     dv0inITraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0i,kind=4),1,MPI_REAL,
     x                     dv0iTraj,1,MPI_REAL,
     x                     master,world,ierr) 

       endif ! end check for do ions
c -----------------------------------------------------------------       


c ###############################################################
c -----------------------------------------------------------------
c ---------------  Lipid Section ----------------------------
c -----------------------------------------------------------------

       if (doLipid) then 
c -----------------------------------------------------------------
        do lidx = 1, nLipStEnd
         iLipStart = arrLipStart(lidx)
         iLipEnd   = arrLipEnd(lidx)

        do iLip = iLipStart, iLipEnd ! Loop over all counter lipids
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        
        
         lipQ = real(qrms(iLip)%q,kind=8)
         lipX = real(nodeConf(iLip)%x,kind=8)                     
         lipY = real(nodeConf(iLip)%y,kind=8)
         lipZ = real(nodeConf(iLip)%z,kind=8)
         alphaI = alphaDbase(iLip)
         
         do irc = 1, iNumDB  ! Loop over RC atoms
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)

          delQ =  real(dVosDbase(irc)%q,kind=8) ! new
          !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)          
          rcQ2 = rcQ1 + delQ  !new


          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

          dXi = lipX - rcX
          dYi = lipY - rcY
          dZi = lipZ - rcZ

          call imageCoords(dXi,dYi,dZi,dX,dY,dZ)

          rIJ = dX*dX+dY*dY+dZ*dZ
          
          rinv = 1.d0/sqrt(rIJ)
          rinv3 = rinv**3
           
          Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
          Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
          Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
          Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
          Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
          Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          
     
         dv0l = dv0l + 
     x          lipQ*delQ*rinv ! Calculate the columb(ion) part
          
         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inL = dv0inL - 0.5d0*alphaI*(E2sq-E1sq)     
          
        enddo ! End the loop over lipids
        enddo ! End the loop over the protein starts/ends
c -----------------------------------------------------------------
c         if(iConf .eq. 1) then

c          print *, "dv0l= ", dv0l
c          print *, "lipQ= ", lipQ
c         endif

         call MPI_GATHER(real(dv0inL,kind=4),1,MPI_REAL,
     x                     dv0inLTraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0l,kind=4),1,MPI_REAL,
     x                     dv0lTraj,1,MPI_REAL,
     x                     master,world,ierr) 


       endif ! end check for do lipids
c ----------------------------------------------------------------- 

c -----------------------------------------------------------------
c ---------------  Other Atom Section ----------------------------
c -----------------------------------------------------------------

       if (doOthers) then 
c -----------------------------------------------------------------       
        do oidx = 1, nOthStEnd
         iOthStart = arrOthStart(oidx)
         iOthEnd   = arrOthEnd(oidx)

        do iOth = iOthStart, iOthEnd ! Loop over all counter lipids
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        
        
         othQ = real(qrms(iOth)%q,kind=8)
         othX = real(nodeConf(iOth)%x,kind=8)                     
         othY = real(nodeConf(iOth)%y,kind=8)
         othZ = real(nodeConf(iOth)%z,kind=8)
         alphaI = alphaDbase(iOth)
         
         do irc = 1, iNumDB  ! Loop over RC atoms
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)

          delQ =  real(dVosDbase(irc)%q,kind=8) ! new
          !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)          
          rcQ2 = rcQ1 + delQ  !new


          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

          dXi = othX - rcX
          dYi = othY - rcY
          dZi = othZ - rcZ

          call imageCoords(dXi,dYi,dZi,dX,dY,dZ)
          rIJ = dX*dX+dY*dY+dZ*dZ
          
          rinv = 1.d0/sqrt(rIJ)
          rinv3 = rinv**3
           
          Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
          Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
          Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
          Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
          Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
          Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          
     
         dv0o = dv0o + 
     x          othQ*delQ*rinv ! Calculate the columb(others) part
          
         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inO = dv0inO - 0.5d0*alphaI*(E2sq-E1sq)     
          
        enddo ! End the loop over other atoms
        enddo ! End the loop over the protein starts/ends
c -----------------------------------------------------------------

         call MPI_GATHER(real(dv0inO,kind=4),1,MPI_REAL,
     x                     dv0inOTraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0o,kind=4),1,MPI_REAL,
     x                     dv0oTraj,1,MPI_REAL,
     x                     master,world,ierr) 

       endif ! end check for do other atoms
c ----------------------------------------------------------------- 

  
c ###############################################################        
c -----------------------------------------------------------------
c ---------------  Water Atom Section -----------------------------
c -----------------------------------------------------------------    
        if(doSolv) then 
        do widx = 1, nWatStEnd
         iWatStart = arrWatStart(widx)
         iWatEnd   = arrWatEnd(widx)

         do iw = iWatStart, iWatEnd
         !do iw = iWatStart, natom        
          Efield1(:) = 0.d0
          Efield2(:) = 0.d0    
          wQ = real(qrms(iw)%q,kind=8)
          wX = real(nodeConf(iw)%x,kind=8)                     
          wY = real(nodeConf(iw)%y,kind=8)
          wZ = real(nodeConf(iw)%z,kind=8)          
          alphaW = alphaDbase(iw)
         
         do irc = 1, iNumDB  ! Loop over all database reaction center atoms
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)
          delQ =  real(dVosDbase(irc)%q,kind=8)
          !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)          
          rcQ2 = rcQ1 + delQ
          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

          dXi = wX - rcX
          dYi = wY - rcY
          dZi = wZ - rcZ

          call imageCoords(dXi,dYi,dZi,dX,dY,dZ)
          rIJ = dX*dX+dY*dY+dZ*dZ
         
          rinv = 1.d0/sqrt(rIJ)
          rinv3 = rinv**3
           
          Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
          Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
          Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
          Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
          Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
          Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          
     
         dv0s = dv0s + wQ*delQ*rinv ! Calculate the columb(water) part          
        enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inS = dv0inS - 0.5d0*alphaW*(E2sq-E1sq)     
          
        enddo ! End the loop over all waters     
       enddo ! End loop over water starts/ends

         call MPI_GATHER(real(dv0inS,kind=4),1,MPI_REAL,
     x                     dv0inSTraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0s,kind=4),1,MPI_REAL,
     x                     dv0sTraj,1,MPI_REAL,
     x                     master,world,ierr) 
     
       endif ! End check for doSolv        
c -----------------------------------------------------------------
c ---------------  End Calculations  ------------------------------
c -----------------------------------------------------------------                        
          if(rank .eq. master)  then
           
            do iWrite = 1, nProc

          write(iEnerAll, '(i6)',advance="no") 
     x    iConf+iWrite-1
          write(iEnerAll, '(a)',advance="no") " "

          if (doProt) then
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0pTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0inPTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
          endif 

          if (doSolv) then
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0sTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0inSTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
          endif 

          if (doIons) then
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0iTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0inITraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
          endif 

          if (doLipid) then
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0lTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0inLTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
          endif 

          if (doOthers) then
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0oTraj(iWrite)
            write(iEnerAll, '(a)',advance="no") " "
            write(iEnerAll, '(E11.4)',advance="no")  
     x            dv0inOTraj(iWrite)
          endif 

          write(iEnerAll,'(a)',advance="yes") ""
            
            enddo ! end loop over iwrite
           endif ! End check for master       
        
        return
      end subroutine calcEnerAll

      subroutine calcEnerAllOLD()
c    *******************************************************************
c    ** calculates induction + coulomb energetics at every timestep   **
c    ** using a database of atoms that are involved in the Stokes     **
c    ** Shift correlation function and all other atoms.               **
c    **                                                               **
c    **                   DEPRECATED VERSION                          **
c    **                                                               **
c    **  The current version has calculations for the protein, ions,  **
c    **  water, lipids, and other atoms, including the support        **
c    **  for multi-section psf.conf files.                            **
c    **                                                               **
c    **                   DEPRECATED VERSION                          **
c    **                                                               **
c    *******************************************************************            
        implicit none
        
        ! Locals
        integer :: iw, irc, iWrite
        integer :: iAtm, iiP, iIon
        integer :: ix, iy, iz        
        integer, parameter :: XX=1
        integer, parameter :: YY=2
        integer, parameter :: ZZ=3
        integer :: xInt, yInt, zInt
        real    :: dX, dY, dZ         
        real*8  :: Efield1(3),Efield2(3)
        real*8  :: dQ, rX, rY, rZ, rIJ
        real*8  :: pQ, pX, pY, pZ
        real*8  :: lowDX, lowDY, lowDZ
        real*8  :: ionQ, ionX, ionY, ionZ        
        real*8  :: wQ, wX, wY, wZ, delQ        
        real*8  :: rcX, rcY, rcZ, rcQ1, rcQ2
        real*8  :: rinv, rinv3, E1sq,E2sq        
        real*8  :: dfrx, dfry, dfrz, delr2
        real*8  :: alphaP, alphaI, alphaW 
        logical :: skipAtom
        real    :: myQ
c    ******************************************************************* 

        if (iConf .eq. 1 ) then
         if ( rank .eq. master) then 
          open(unit=iInd,file="enerInd.dat",status="unknown") 
          open(unit=iEne,file="enerCoul.dat",status="unknown")                         
         endif
        endif 
        
        dv0in = 0.d0         
        dv0s = 0.d0
        dv0p = 0.d0
        dv0inS = 0.d0
        dv0inP = 0.d0
        
c -----------------------------------------------------------------
c ---------------  Non-RC protein Atom Section --------------------
c -----------------------------------------------------------------    
       if (doProt) then 
        do iip = 1, iProtEnd ! Loop over all non-RC, protein atoms
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        
        
         pQ = real(qrms(iip)%q,kind=8)
         pX = real(nodeConf(iip)%x,kind=8)                     
         pY = real(nodeConf(iip)%y,kind=8)
         pZ = real(nodeConf(iip)%z,kind=8)
         alphaP = alphaDbase(iip)
         
         do irc = 1, iNumDB  ! Loop over all RC atoms     
          iAtm = dVosDbase(irc)%atom          
          call searchDBAtom(iip,myQ,skipAtom)
          if (.not. skipAtom) then
          
          !if (iip .ne. iAtm) then ! Make sure we are on different atoms          
           rcQ1 = real(qrms(iAtm)%q,kind=8)
           delQ = real(dVosDbase(irc)%q,kind=8)
           !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)
           rcQ2 = rcQ1 + delQ
           rcX = real(nodeConf(iAtm)%x,kind=8)                     
           rcY = real(nodeConf(iAtm)%y,kind=8)
           rcZ = real(nodeConf(iAtm)%z,kind=8)

           dX = pX - rcX
           dY = pY - rcY
           dZ = pZ - rcZ
          
           rinv = 1.d0/sqrt(dX*dX + dY*dY + dZ*dZ)
           rinv3 = rinv**3
           
           Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
           Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
           Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
           Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
           Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
           Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          

          dv0p = dv0p + pQ*delQ*rinv ! Calculate the columb(prot) part
          
          endif
          
         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inP = dv0inP - 0.5d0*alphaP*(E2sq-E1sq)     
          
        enddo ! End the loop over non-RC, protein atoms
        
         call MPI_GATHER(real(dv0inP,kind=4),1,MPI_REAL,
     x                     dv0inPTraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0p,kind=4),1,MPI_REAL,
     x                     dv0pTraj,1,MPI_REAL,
     x                     master,world,ierr) 

       endif ! End check for doProt
        
c -----------------------------------------------------------------
c ---------------  Counter Ion Section ----------------------------
c -----------------------------------------------------------------
       if(0 .eq. 1) then ! always off, for now...
c -----------------------------------------------------------------       
        do iIon = iProtEnd+1, iWatStart-1 ! Loop over all counter ions
         Efield1(:) = 0.d0
         Efield2(:) = 0.d0        
        
         ionQ = real(qrms(iIon)%q,kind=8)
         ionX = real(nodeConf(iIon)%x,kind=8)                     
         ionY = real(nodeConf(iIon)%y,kind=8)
         ionZ = real(nodeConf(iIon)%z,kind=8)
         alphaI = alphaDbase(iIon)
         
         do irc = 1, iNumDB  ! Loop over RC atoms
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)
          rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)
          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then 
          dX = ionX - rcX
          dY = ionY - rcY
          dZ = ionZ - rcZ

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          rIJ = dX**2 + dY**2 + dZ**2
          !rIJ = 1000000.0
          lowDX = dX
          lowDY = dY              
          lowDZ = dZ   
              
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
              lowDX = dX
              lowDY = dY              
              lowDZ = dZ              
            endif

            enddo
           enddo
          enddo
          
           dX = lowDX
           dY = lowDY
           dZ = lowDZ

          else  ! its orthorhombic

           xInt = ANINT((ionX-rcX)/bX)
           yInt = ANINT((ionY-rcY)/bY)
           zInt = ANINT((ionZ-rcZ)/bZ)

           dX = (ionX-rcX) - xInt*bX
           dY = (ionY-rcY) - yInt*bY
           dZ = (ionZ-rcZ) - zInt*bZ         
           rIJ = dX*dX+dY*dY+dZ*dZ
          
          endif ! End test for truncated box           
          
           rinv = 1.d0/sqrt(rIJ)
           rinv3 = rinv**3
           
           Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
           Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
           Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
           Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
           Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
           Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          
     
          
         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0in = dv0in - 0.5d0*alphaI*(E2sq-E1sq)     
          
        enddo ! End the loop over counter ions   
c -----------------------------------------------------------------
       endif ! end always off
c -----------------------------------------------------------------       

        
c -----------------------------------------------------------------
c ---------------  Water Atom Section -----------------------------
c -----------------------------------------------------------------    
        if(doSolv) then 
         do iw = iWatStart, natom        
          Efield1(:) = 0.d0
          Efield2(:) = 0.d0    
          wQ = real(qrms(iw)%q,kind=8)
          wX = real(nodeConf(iw)%x,kind=8)                     
          wY = real(nodeConf(iw)%y,kind=8)
          wZ = real(nodeConf(iw)%z,kind=8)          
          alphaW = alphaDbase(iw)
         
         do irc = 1, iNumDB  ! Loop over all database reaction center atoms
          iAtm = dVosDbase(irc)%atom
          rcQ1 = real(qrms(iAtm)%q,kind=8)
          delQ =  real(dVosDbase(irc)%q,kind=8)
          !rcQ2 = rcQ1 + real(dVosDbase(irc)%q,kind=8)          
          rcQ2 = rcQ1 + delQ
          rcX = real(nodeConf(iAtm)%x,kind=8)                     
          rcY = real(nodeConf(iAtm)%y,kind=8)
          rcZ = real(nodeConf(iAtm)%z,kind=8)

         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then 
          dX = wX - rcX
          dY = wY - rcY
          dZ = wZ - rcZ

          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          !rIJ = dX**2 + dY**2 + dZ**2 ! working? ...maybe
          rIJ = 1000000.0
          lowDX = dX
          lowDY = dY              
          lowDZ = dZ   
              
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             delr2 = dX*dX + dY*dY + dZ*dZ

            if(delr2.lt.rIJ) then 
              rIJ = delr2
              lowDX = dX
              lowDY = dY              
              lowDZ = dZ              
            endif

            enddo
           enddo
          enddo
          
           dX = lowDX
           dY = lowDY
           dZ = lowDZ
           
         else ! its an orthorhombic box
         
          xInt = ANINT((wX-rcX)/bX)
          yInt = ANINT((wY-rcY)/bY)
          zInt = ANINT((wZ-rcZ)/bZ)

          dX = (wX-rcX) - xInt*bX
          dY = (wY-rcY) - yInt*bY
          dZ = (wZ-rcZ) - zInt*bZ         
          rIJ = dX*dX+dY*dY+dZ*dZ

         endif
         
           rinv = 1.d0/sqrt(rIJ)
           rinv3 = rinv**3
           
           Efield1(XX) = Efield1(XX) + 
     x                 dX*rcQ1*rinv3
           Efield1(YY) = Efield1(YY) + 
     x                 dY*rcQ1*rinv3
           Efield1(ZZ) = Efield1(ZZ) + 
     x                 dZ*rcQ1*rinv3          
     
           Efield2(XX) = Efield2(XX) + 
     x                 dX*rcQ2*rinv3
           Efield2(YY) = Efield2(YY) + 
     x                 dY*rcQ2*rinv3
           Efield2(ZZ) = Efield2(ZZ) + 
     x                 dZ*rcQ2*rinv3          
     
          dv0s = dv0s + wQ*delQ*rinv ! Calculate the columb(water) part          
         enddo ! End loop over RC atoms
         
          E2sq = Efield2(XX)*Efield2(XX) + 
     x           Efield2(YY)*Efield2(YY) + 
     x           Efield2(ZZ)*Efield2(ZZ)

          E1sq = Efield1(XX)*Efield1(XX) + 
     x           Efield1(YY)*Efield1(YY) + 
     x           Efield1(ZZ)*Efield1(ZZ)
     
         dv0inS = dv0inS - 0.5d0*alphaW*(E2sq-E1sq)     
          
        enddo ! End the loop over all waters     

         call MPI_GATHER(real(dv0inS,kind=4),1,MPI_REAL,
     x                     dv0inSTraj,1,MPI_REAL,
     x                     master,world,ierr)     

         call MPI_GATHER(real(dv0s,kind=4),1,MPI_REAL,
     x                     dv0sTraj,1,MPI_REAL,
     x                     master,world,ierr) 
     
       endif ! End check for doSolv        
c -----------------------------------------------------------------
c ---------------  End Calculations  ------------------------------
c -----------------------------------------------------------------                        
          if(rank .eq. master)  then
           
            do iWrite = 1, nProc
              if (doProt .and. doSolv) then
               write(iEne, *) iConf+iWrite-1, 
     x         dv0sTraj(iWrite), dv0pTraj(iWrite)   
               write(iInd, *) iConf+iWrite-1, 
     x         dv0inSTraj(iWrite), dv0inPTraj(iWrite)
              elseif(doProt) then
               write(iEne, *) iConf+iWrite-1, 
     x         dv0pTraj(iWrite)   
               write(iInd, *) iConf+iWrite-1, 
     x         dv0inPTraj(iWrite)
              else ! Just solvent
               write(iEne, *) iConf+iWrite-1, 
     x         dv0sTraj(iWrite)
               write(iInd, *) iConf+iWrite-1, 
     x         dv0inSTraj(iWrite)              
              endif
            
            enddo
           endif ! End check for master       
        
        return
      end subroutine calcEnerAllOLD
      
      
      subroutine mu0s_direct(dMu0s_i,ddX, ddY, ddZ,dQ)
c    *******************************************************************
c    ** calculates potential energy of solvation directly for         **
c    ** all water molecules and a given reaction center atom.         **
c    **                                                               **
c    ** principal variables:                                          **
c    ** real*8    dMu0s_i   potential energy difference (for SS-TCF)  **
c    ** real   ddX       the x position                               **
c    ** real   ddY       the y position                               **
c    ** real   ddZ       the z position                               **
c    ** real   dQ        the charge on the protein atom               **
c    **                                                               **
c    *******************************************************************
       implicit none
       
        ! Args
        real*8 :: dMu0s_i
        !integer :: iProt
        real    :: ddX, ddY, ddZ, dQ

        ! Locals
        integer :: ii, iw
        integer :: xInt, yInt, zInt
        integer ix, iy, iz
        real*8 :: delr2m, delr2
        real*8 :: wQ, wX, wY, wZ
        real*8 :: prX, prY, prZ, prQ
        real*8 :: rIJ, dX, dY, dZ, rIJ2, lowRIJ2
        real*8 :: dfrx, dfry, dfrz        
c    ******************************************************************* 

        dMu0s_i = 0.d0

        prX = real(ddX,kind=8) !real(nodeConf(iProt)%x,kind=8)
        prY = real(ddY,kind=8) !real(nodeConf(iProt)%y,kind=8)
        prZ = real(ddZ,kind=8) !real(nodeConf(iProt)%z,kind=8)
        prQ = real(dQ,kind=8)
                     
!$OMP PARALLEL
!$OMP&DEFAULT(shared) 
!$OMP&REDUCTION(+:dMu0s_i)
!$OMP&PRIVATE(ii,wX,wY,wZ,wQ,dX,dY,dZ)
!$OMP&PRIVATE(rIJ,xInt,yInt,zInt)
!$OMP DO SCHEDULE(static)
        !print *, " iWatStart = ", iWatStart
        !print *, " natom = ", natom  
        !print *, rank, "bXYZ = ", bX, bY, bZ        


        !do ii = iWatStart, natom  ! ORIGINAL WITH natom being the last atom that is searched
        do ii = iWatStart, iWatEnd
        ! Store original water position               

         wX = real(nodeConf(ii)%x,kind=8)
         wY = real(nodeConf(ii)%y,kind=8)
         wZ = real(nodeConf(ii)%z,kind=8)
         wQ = real(qrms(ii)%q,kind=8)

         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then         
          dX = prX-wX
          dY = prY-wY
          dZ = prZ-wZ
           
          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          lowRIJ2 = 10000000.00
          !lowRIJ2 = dX**2+dY**2+dZ**2
          
          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             rIJ2 = dX*dX + dY*dY + dZ*dZ

            if(rIJ2.lt.lowRIJ2) then
             lowRIJ2 = rIJ2
             
            endif
            
            enddo
           enddo
          enddo
         
!$OMP ATOMIC                  
          rIJ = sqrt(lowRIJ2)
         
         else ! orthorhombic box

!$OMP CRITICAL
c           xInt = ANINT((prX-wX)/bX)
c           yInt = ANINT((prY-wY)/bY)
c           zInt = ANINT((prZ-wZ)/bZ)
!$OMP END CRITICAL                

c           dX = (prX-wX) - xInt*bX
c           dY = (prY-wY) - yInt*bY
c           dZ = (prZ-wZ) - zInt*bZ         
c           rIJ = sqrt(dX**2+dY**2+dZ*2)
           
         dX = (prX-wX) - bX * ANINT((prX-wX)/bX)
         dY = (prY-wY) - bY * ANINT((prY-wY)/bY)
         dZ = (prZ-wZ) - bZ * ANINT((prZ-wZ)/bZ)
         
         rIJ = sqrt(dX**2+dY**2+dZ**2 )
         
         endif ! end test for a truncated box
         

c         if(ii .eq. iWatStart) then
c           print *, "FORTRAN: O:: prQ * wQ / rIJ = ", (prQ * wQ / rIJ)
c           print *, "FORTRAN: O:: dX,dY,dZ  = ", dX,dY,dZ
c           print *, "FORTRAN: O:: pX,pY,pZ  = ", prX,prY,prZ
c           print *, "FORTRAN: O:: wX,wY,wZ  = ", wX,wY,wZ
c           print *, "FORTRAN: O:: prQ = ", prQ, 
c     x              " wQ  = ", wQ, " rIJ = ", rIJ
c         endif

c         if(ii .eq. (iWatStart + 1)) then
c           print *, "FORTRAN: H1:: prQ * wQ / rIJ = ", (prQ * wQ / rIJ)
c           print *, "FORTRAN: H1:: dX,dY,dZ  = ", dX,dY,dZ
c           print *, "FORTRAN: H1:: pX,pY,pZ  = ", prX,prY,prZ
c           print *, "FORTRAN: H1:: wX,wY,wZ  = ", wX,wY,wZ
c           print *, "FORTRAN: H1:: prQ = ", prQ, 
c     x              " wQ  = ", wQ, " rIJ = ", rIJ
c         endif

c         if(ii .eq. (iWatStart + 2)) then
c           print *, "FORTRAN: H2:: prQ * wQ / rIJ = ", (prQ * wQ / rIJ)
c           print *, "FORTRAN: H2:: dX,dY,dZ  = ", dX,dY,dZ
c           print *, "FORTRAN: H2:: pX,pY,pZ  = ", prX,prY,prZ
c           print *, "FORTRAN: H2:: wX,wY,wZ  = ", wX,wY,wZ
c           print *, "FORTRAN: H2:: prQ = ", prQ, 
c     x              " wQ  = ", wQ, " rIJ = ", rIJ
c         endif
         
         dMu0s_i = dMu0s_i + wQ / rIJ

       enddo ! end loop over atoms in molecule       
!$OMP END DO NOWAIT
!$OMP END PARALLEL 
        dMu0s_i = dMu0s_i*prQ
        
        return
      end subroutine mu0s_direct

      subroutine mu0i_direct(dMu0i_i,ddX, ddY, ddZ,dQ)
c    *******************************************************************
c    ** calculates potential energy of solvation directly for         **
c    ** all counter ions and a given reaction center atom.            **
c    **                                                               **
c    ** principal variables:                                          **
c    ** real*8    dMu0i_i   potential energy difference (for SS-TCF)  **
c    ** real   ddX       the x position                               **
c    ** real   ddY       the y position                               **
c    ** real   ddZ       the z position                               **
c    ** real   dQ        the charge on the protein atom               **
c    **                                                               **
c    *******************************************************************
       implicit none
       
        ! Args
        real*8 :: dMu0i_i
        !integer :: iProt
        real    :: ddX, ddY, ddZ,dQ

        ! Locals
        integer :: ii, iw
        integer :: xInt, yInt, zInt
        integer :: ix, iy, iz
        integer :: iS, iE
        real*8 :: delr2m, delr2
        real*8 :: iQ, iiX, iiY, iiZ
        real*8 :: prX, prY, prZ, prQ
        real*8 :: rIJ, dX, dY, dZ, rIJ2, lowRIJ2
        real*8 :: dfrx, dfry, dfrz        
c    ******************************************************************* 

        dMu0i_i = 0.d0

        prX = real(ddX,kind=8) !real(nodeConf(iProt)%x,kind=8)
        prY = real(ddY,kind=8) !real(nodeConf(iProt)%y,kind=8)
        prZ = real(ddZ,kind=8) !real(nodeConf(iProt)%z,kind=8)
        prQ = real(dQ,kind=8)
        
        iS = iProtEnd+1
        iE = iWatStart - 1 
        do ii = iS, iE
        ! Store original ion position               

         iiX = real(nodeConf(ii)%x,kind=8)
         iiY = real(nodeConf(ii)%y,kind=8)
         iiZ = real(nodeConf(ii)%z,kind=8)
         iQ = real(qrms(ii)%q,kind=8)

         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then         
          dX = prX-iiX
          dY = prY-iiY
          dZ = prZ-iiZ
           
          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          lowRIJ2 = dX**2+dY**2+dZ**2

          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             rIJ2 = dX*dX + dY*dY + dZ*dZ

            if(rIJ2.lt.lowRIJ2) then
             lowRIJ2 = rIJ2
            endif
            
            enddo
           enddo
          enddo
         
          rIJ = sqrt(lowRIJ2)
         
         else ! orthorhombic box

c           xInt = ANINT((prX-iiX)/bX)
c           yInt = ANINT((prY-iiY)/bY)
c           zInt = ANINT((prZ-iiZ)/bZ)

c           dX = (prX-iiX) - xInt*bX
c           dY = (prY-iiY) - yInt*bY
c           dZ = (prZ-iiZ) - zInt*bZ         
c           rIJ = sqrt(dX**2+dY**2+dZ*2)
           
c         dX = (prX-iiX) - bX * ANINT((prX-iiX)/bX)
c         dY = (prY-iiY) - bY * ANINT((prY-iiY)/bY)
c         dZ = (prZ-iiZ) - bZ * ANINT((prZ-iiZ)/bZ)
         
c         rIJ = sqrt(dX**2+dY**2+dZ**2 )
         
           xInt = ANINT((prX-iiX)/bX)
           yInt = ANINT((prY-iiY)/bY)
           zInt = ANINT((prZ-iiZ)/bZ)

           dX = (prX-iiX) - xInt*bX
           dY = (prY-iiY) - yInt*bY
           dZ = (prZ-iiZ) - zInt*bZ         
           rIJ = sqrt(dX*dX+dY*dY+dZ*dZ)         
         
         endif ! end test for a truncated box
         
         
         dMu0i_i = dMu0i_i + iQ / rIJ

       enddo ! end loop over all counter ions       

       dMu0i_i = dMu0i_i*prQ
        
        return
      end subroutine mu0i_direct
      
      subroutine mu0p_direct(dMu0p_i,ddX, ddY, ddZ,dQ)
c    *******************************************************************
c    ** calculates protein's portion of the potential energy          **
c    ** given a particular protein atom and the other non-RC atoms.   **
c    **                                                               **
c    ** principal variables:                                          **
c    ** real*8    dMu0p_i   protein vertical energy                   **
c    ** real   ddX       the x position                               **
c    ** real   ddY       the y position                               **
c    ** real   ddZ       the z position                               **
c    ** real   dQ        the charge on the protein atom               **
c    **                                                               **
c    *******************************************************************
       implicit none
       
        ! Args
        real*8 :: dMu0p_i
        !integer :: iProt
        real    :: ddX, ddY, ddZ,dQ

        ! Locals
        integer :: ii
        real*8 :: pQ, pX, pY, pZ
        real*8 :: prX, prY, prZ, prQ
        real*8 :: rIJ, dX, dY, dZ
        real   :: myQ
        logical :: skipAtom
c    ******************************************************************* 

        dMu0p_i = 0.d0

        prX = real(ddX,kind=8) !real(nodeConf(iProt)%x,kind=8)
        prY = real(ddY,kind=8) !real(nodeConf(iProt)%y,kind=8)
        prZ = real(ddZ,kind=8) !real(nodeConf(iProt)%z,kind=8)
        prQ = real(dQ,kind=8)          

!$OMP PARALLEL
!$OMP&DEFAULT(shared) 
!$OMP&REDUCTION(+:dMu0p_i)
!$OMP&PRIVATE(ii,pX,pY,pZ,pQ,dX,dY,dZ,rIJ)
!$OMP DO SCHEDULE(static)
        do ii = 1, iProtEnd
         
        ! Skip any atom in the reaction center         
         call searchDBAtom(ii,myQ,skipAtom)
         if (.not. skipAtom) then
         !if (iProt .ne. ii) then
         
          pX = real(nodeConf(ii)%x,kind=8)
          pY = real(nodeConf(ii)%y,kind=8)
          pZ = real(nodeConf(ii)%z,kind=8)
          pQ = real(qrms(ii)%q,kind=8)

          dX = prX-pX
          dY = prY-pY
          dZ = prZ-pZ

!$OMP ATOMIC         
          rIJ = sqrt(dX*dX+dY*dY+dZ*dZ)         
          dMu0p_i = dMu0p_i + pQ / rIJ
          
         endif ! Chek to see if the atom is skipped
         
        enddo ! Loop over protein atoms
!$OMP END DO NOWAIT
!$OMP END PARALLEL 
        dMu0p_i = dMu0p_i*prQ
        
        return
      end subroutine mu0p_direct

      subroutine mu0o_direct(dMu0o_i,ddX, ddY, ddZ,dQ)
c    *******************************************************************
c    ** calculates potential energy of solvation directly for         **
c    ** all other atoms and a given reaction center atom.             **
c    **                                                               **
c    ** principal variables:                                          **
c    ** real*8    dMu0o_i   potential energy difference (for SS-TCF)  **
c    ** real   ddX       the x position                               **
c    ** real   ddY       the y position                               **
c    ** real   ddZ       the z position                               **
c    ** real   dQ        the charge on the protein atom               **
c    **                                                               **
c    *******************************************************************
       implicit none
       
        ! Args
        real*8 :: dMu0o_i
        !integer :: iProt
        real    :: ddX, ddY, ddZ,dQ

        ! Locals
        integer :: ii, iw
        integer :: xInt, yInt, zInt
        integer :: ix, iy, iz
        integer :: iS, iE
        real*8 :: delr2m, delr2
        real*8 :: iQ, iiX, iiY, iiZ
        real*8 :: prX, prY, prZ, prQ
        real*8 :: rIJ, dX, dY, dZ, rIJ2, lowRIJ2
        real*8 :: dfrx, dfry, dfrz        
c    ******************************************************************* 

        dMu0o_i = 0.d0

        prX = real(ddX,kind=8) !real(nodeConf(iProt)%x,kind=8)
        prY = real(ddY,kind=8) !real(nodeConf(iProt)%y,kind=8)
        prZ = real(ddZ,kind=8) !real(nodeConf(iProt)%z,kind=8)
        prQ = real(dQ,kind=8)
        
        do ii = iOthStart, iOthEnd
        ! Store original ion position               

         iiX = real(nodeConf(ii)%x,kind=8)
         iiY = real(nodeConf(ii)%y,kind=8)
         iiZ = real(nodeConf(ii)%z,kind=8)
         iQ = real(qrms(ii)%q,kind=8)

         ! Reimage the box in the for truncated or orthorhombic boxes
         if (isTruncated) then         
          dX = prX-iiX
          dY = prY-iiY
          dZ = prZ-iiZ
           
          dfrx = dX*recip(1,1)+dY*recip(2,1)+dZ*recip(3,1)
          dfry = dX*recip(1,2)+dY*recip(2,2)+dZ*recip(3,2)
          dfrz = dX*recip(1,3)+dY*recip(2,3)+dZ*recip(3,3)

          dfrx = dfrx - ANINT(dfrx)
          dfry = dfry - ANINT(dfry)
          dfrz = dfrz - ANINT(dfrz)

          lowRIJ2 = dX**2+dY**2+dZ**2

          do ix=-1,1
           do iy=-1,1
            do iz=-1,1

             dX = (dfrx+ix)*ucell(1,1) + 
     x            (dfry+iy)*ucell(1,2) + 
     x            (dfrz+iz)*ucell(1,3)
             dY = (dfrx+ix)*ucell(2,1) + 
     x            (dfry+iy)*ucell(2,2) + 
     x            (dfrz+iz)*ucell(2,3)
             dZ = (dfrx+ix)*ucell(3,1) + 
     x            (dfry+iy)*ucell(3,2) + 
     x            (dfrz+iz)*ucell(3,3)

             rIJ2 = dX*dX + dY*dY + dZ*dZ

            if(rIJ2.lt.lowRIJ2) then
             lowRIJ2 = rIJ2
            endif
            
            enddo
           enddo
          enddo
         
          rIJ = sqrt(lowRIJ2)
         
         else ! orthorhombic box

c           xInt = ANINT((prX-iiX)/bX)
c           yInt = ANINT((prY-iiY)/bY)
c           zInt = ANINT((prZ-iiZ)/bZ)

c           dX = (prX-iiX) - xInt*bX
c           dY = (prY-iiY) - yInt*bY
c           dZ = (prZ-iiZ) - zInt*bZ         
c           rIJ = sqrt(dX**2+dY**2+dZ*2)
           
c         dX = (prX-iiX) - bX * ANINT((prX-iiX)/bX)
c         dY = (prY-iiY) - bY * ANINT((prY-iiY)/bY)
c         dZ = (prZ-iiZ) - bZ * ANINT((prZ-iiZ)/bZ)
         
c         rIJ = sqrt(dX**2+dY**2+dZ**2 )
         
           xInt = ANINT((prX-iiX)/bX)
           yInt = ANINT((prY-iiY)/bY)
           zInt = ANINT((prZ-iiZ)/bZ)

           dX = (prX-iiX) - xInt*bX
           dY = (prY-iiY) - yInt*bY
           dZ = (prZ-iiZ) - zInt*bZ         
           rIJ = sqrt(dX*dX+dY*dY+dZ*dZ)         
         
         endif ! end test for a truncated box
         
         
         dMu0o_i = dMu0o_i + iQ / rIJ

       enddo ! end loop over all counter ions       

       dMu0o_i = dMu0o_i*prQ
        
        return
      end subroutine mu0o_direct
      

      subroutine forceMPIStop()
c    *******************************************************************
c    ** this is a debugging routine to help put to the program to bed **
c    ** immediately without memory leaks, aside from deallocation.    **
c    **                                                               **
c    ** principal variables:                                          **
c    *******************************************************************
       implicit none

        call MPI_Finalize(ierr)     ! Cleanup MPI processes on all nodes     
        stop '** PROGRAM PUT TO REST GENTLY **'

        return
      end subroutine forceMPIStop

      end module fastdef
