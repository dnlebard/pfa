      character*80 ititl
      character*4 igraph(maxnatom),labres(maxnres),isymbl(maxnatom)
      character*4 itree(maxnatom)

      integer natom,ntypes,nbonh,mbona,ntheth,mtheta
      integer nphih,mphia,nhparm,nparm,next,nres,nbona,ntheta
      integer nphia,numbnd,numang,nptra,natyp,nphb,ifpert,nbper
      integer ngper,ndper,mbper,mgper,mdper,ifbox,nmxrs,ifcap

      integer iac(maxnatom),numex(maxnatom),ico(maxntypes*maxntypes)
      integer ipres(maxnres),ibh(maxbnd),jbh(maxbnd),icbh(maxbnd) 
      integer ib(maxbnd),jb(maxbnd),icb(maxbnd) 
      integer ith(maxang),jth(maxang),kth(maxang),icth(maxang) 
      integer it(maxang),jt(maxang),kt(maxang),ict(maxang) 
      integer iph(maxdih),jph(maxdih),kph(maxdih),lph(maxdih),
     x        icph(maxdih) 
      integer ip(maxdih),jp(maxdih),kp(maxdih),lp(maxdih),icp(maxdih) 
      integer natex(maxnext),join(maxnatom),irotat(maxnatom)
      integer iptres,nspm,nspsol,nsp(maxnspm),natcap,ipol,numextra
      integer ibper(maxbnd),jbper(maxbnd),icbper(maxbnd)
      integer itper(maxang),jtper(maxang),ktper(maxang),ictper(maxang)
      integer ipper(maxdih),jpper(maxdih),kpper(maxdih),lpper(maxdih)
      integer icpper(maxdih)
      integer iaper(maxnatom),iacper(maxnatom)
      integer inegp(maxdih),jnegp(maxdih),knegp(maxdih),lnegp(maxdih)
      integer natomo
      real*8  cgper(maxnatom),almper(maxnatom)
      real*8  ew_alpha,ew_beta,ew_gamma,time
      logical ewald_box,havetime
      character*4 labre2(maxnres),igrper(maxnatom),ismper(maxnatom)
      integer inegh(maxdih),jnegh(maxdih),knegh(maxdih),lnegh(maxdih)
      integer ineg(maxdih),jneg(maxdih),kneg(maxdih),lneg(maxdih)
      real*8 chrg(maxnatom),amass(maxnatom),rkk(maxbnd),req(maxbnd)
      real*8 tk(maxang),teq(maxang),pk(maxdih),pn(maxdih),phase(maxdih)
      real*8 solty(maxnatyp),cn1(maxntypes*(maxntypes+1)/2)
      real*8 cn2(maxntypes*(maxntypes+1)/2),asol(maxnphb),bsol(maxnphb)
      real*8 hbcut(maxnphb),beta,box(3),cutcap,xcap,ycap,zcap 
      real*8 rborn(maxnatom),fs(maxnatom)
      common /c1/ ititl,igraph,labres,isymbl,itree
      common /cp1/ labre2,igrper,ismper
      common /i1/ natom,ntypes,nbonh,mbona,ntheth,mtheta
      common /i2/ nphih,mphia,nhparm,nparm,next,nres,nbona,ntheta
      common /i3/ nphia,numbnd,numang,nptra,natyp,nphb,ifpert,nbper
      common /i4/ ngper,ndper,mbper,mgper,mdper,ifbox,nmxrs,ifcap
      common /i5/ iac,numex,ico,ipres,ibh,jbh,icbh
      common /i6/ ib,jb,icb,ith,jth,kth,icth
      common /i7/ it,jt,kt,ict,iph,jph,kph,lph,icph
      common /i8/ ip,jp,kp,lp,icp,natex,join,irotat
      common /i9/ iptres,nspm,nspsol,nsp,natcap,ipol
      common /i10/ ineg,jneg,kneg,lneg,inegh,jnegh,knegh,lnegh
      common /ip1/ ibper,jbper,icbper,itper,jtper,ktper,ictper
      common /ip2/ ipper,jpper,kpper,lpper,icpper,iaper,iacper
      common /ip3/ inegp,jnegp,knegp,lnegp,numextra
      common /r1/ chrg,amass,rkk,req,tk,teq,pk,pn,phase,rborn,fs
      common /r2/ solty,cn1,cn2,asol,bsol,hbcut
      common /r3/ beta,box,cutcap,xcap,ycap,zcap
      common /r4/ ew_alpha,ew_beta,ew_gamma,time
      common /rp1/ cgper,almper
      common /orig/ natomo
      common /l1/ ewald_box,havetime
