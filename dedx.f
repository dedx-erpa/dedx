c ******************************************************************** c
c                                                                      c
c                            D E V S D X                               c
C  This program computes ion stopping power using Lindhard's stopping
c  power theory.
c  (1) The spherically averaged charge density rho(r) of the target 
c      atom is calculated self-consistently with free electron sea:
c      Muffin-Tin atomic model.
c  (2) Relativistic Dirac-Fock equation is solved.                     c
c                                                                      c
c ******************************************************************** c

      program devsdx
      implicit real*8 (a-h,o-z)
      parameter (nrx=850, npx=1500, nex=100, nre=10)
      parameter (qea = 0.75, qei = 0.01)
      
      character string*80, odir*80, floss*80, floss0*80
      
      dimension r(nrx), dr(nrx),   rho(nrx), rhol(nrx), sum(nrx)
      dimension ee(npx), elog(npx), dedx(npx), dedxt(50,50,200)
      dimension rhof(nrx),rho0(nrx),rhoe(nrx),fdedx(npx),rhofl(npx)
      dimension rdedx(nre,nrx), cdedx(nre,nrx), rage(npx)
      dimension a(10), fits(nex), temp(50), deni(50)

      common /lhtab/eione(50),rl(30),cl(50,30),cl0(50,30),
     +     neee,nrho,xte,xte0
      common / hlbdy / elow, ehig, epeak
      common /rrang/range(50,50,200)

      namelist /dedxinp/zzp, qmass, ztg, amass, mep, emin, emax,
     +     mloss, mout

c                                                                      c
c ************************* START THE EXECUTION ********************** c
c                                                                      c

      Pi = 3.1415926
c     proton mass in g
      pmass = 1.672623100e-24
c     Bohr radius in (cm)
      a0 = 5.29177249e-9

c      write(*,201)
c 201  format(//5x,53(h*)
c     :     /5x,'*                     D E V S D X                   *'
c     :     /5x,'*  a code for computing ion stopping power using an *'
c     :     /5x,'*         unified self-consistent field model       *'
c     :     /5x,'*                                                   *'
c     :     /5x,'*  To run the code, you need copy the corresponding *'
c     :     /5x,'*  electron density distribution function data to   *'
c     :     /5x,'*  input data file (rho.function), and then input   *'
c     :     /5x,'*  parameters through namelist file (dedx.inp)      *'/)

      read(5,*) odir
      odir = trim(adjustl(odir))
      len1 = len(trim(odir))
      read(5,*) floss
      floss = trim(adjustl(floss))
      len2 = len(trim(floss))
      read(5,*) floss0
      floss0 = trim(adjustl(floss0))
      len0 = len(trim(floss0))
      if (floss0(1:len0) .eq. 'NONE') then
         floss0(1:len2) = floss(1:len2)
         len0 = len2
      endif
c      write(*,*) odir
      open(unit=1, file=odir(1:len1)//'/dedx.inp', status='old')
      open(unit=2, file=odir(1:len1)//'/rho.functions', status='old')     
c      open(unit=11,file=odir(1:len1)//'/lh.dat', status='unknown')       
c      open(unit=12,file=odir(1:len1)//'/lh.tab', status='unknown')      
c      open(unit=13,file=odir(1:len1)//'/dedx.dat', status='unknown')
      open(unit=14,file=odir(1:len1)//'/dedx.dat', status='unknown')
      if (mout .gt. 0) then
         open(unit=15,file=odir(1:len1)//'/rdedx.dat',status='unknown')
      endif
c
c ... namelist input      
c     
      open(unit=17, file=odir(1:len1)//'temp', status='unknown')
1999  read(1,9999, end=2999) string
      if(string(1:1).eq.'*') goto 1999 
      write(17, 9999) string(1:20)
      goto 1999

 2999 rewind (17)
      read(17, dedxinp)
      close(17, status='delete')
 9999 format(a)

      
      acoef = 6.023d20/amass
      bcoef = acoef*1d-21
c     write(*,*) mep, emin, emax
c
c ... projectile energy mesh
c
      emin = max(emin,1.0e-4)
      mep  = min(mep,100)

      ee(1) = emin
      mre = mep/nre
      if(mep.gt.1) then
         de = log10(emax/emin)/(mep-1)
         do 2 j=2,mep
            ee(j) = ee(j-1)*10**de
 2       continue
      endif
      do j=1,mep
         elog(j) = dlog(ee(j))
      enddo
      he = elog(2)-elog(1)
c
c ... read in temperature and density points
c

      read(2,*) kzz,nmt, nmd
      do 4 it=1,nmt
      do 3 id=1,nmd
         read(2,*) temp(it), deni(id)      
         read(2,*) nrd,h,rs
         read(2,*) (r(i),  i=1,nrd)
         read(2,*) (dr(i), i=1,nrd)
         read(2,*) (rho0(i),i=1,nrd)
         read(2,*) (rhof(i),i=1,nrd)
         read(2,*) (rhoe(i),i=1,nrd)
 3    continue
 4    continue
      rewind (2)
      read(2,*) kzz, nmt, nmd
      if(abs(kzz-ztg) .gt. 0.1) stop '!!! rho.function incorrect !!!'      
c 916  format(*)
     

c      write(13,1001) 
c 1001 format(130(h*)
c     : /'*'
c     : /'*',40x,'ION STOPPING POWER DATA TABLE'
c     : /'*'
c     : /'*',2x, 
c     :  'This table provides fitting coefficients for stopping '
c     :  'of ion projectile in targets of verious temperature and'
c     : /'*',2x,
c     :  'density conditions. Stopping formulae using tabulated '
c     :  'coefficients are listed in the following:'
c     : /'*',30x,
c     :  ' for        0 < E < 0.1*E0:  dE/dX = A1 * E**A2'
c     : /'*',30x,
c     :  ' for   0.1*E0 < E < 10.*E0:  dE/dX = S1*S2/(S1+S2)'
c     : /'*',30x,                       
c     :  ' where                          S1 = A3 * E**A4'
c     : /'*',30x,                       
c     :  '                                S2 = (A5/E)*ln(1+A6/E+A7*E)'
c     : /'*',30x,
c     :  ' for            E > 10*E0:   dE/dX = (A8/E)*ln(1+A9/E+A10*E)'
c     : /'*'/130(h*)/)            
      
c      write(13,1002) int(ztg)
c 1002 format(2x,'Z_target     = ',i3)
c      write(13,1003) int(zzp)
c 1003 format(2x,'Z_projectile = ',i3)
c      write(13,1004) nmt
c 1004 format(2x,'# of T grid  = ',i3)      
c      write(13,1005) nmd
c 1005 format(2x,'# of N grid  = ',i3)      

      write(14,2000) 1
 2000 format('#   nzt = ', i3)
      write(14,2001) ztg
 2001 format('#    zt = ', e10.4)
      write(14,2002) 1.0
 2002 format('#    wt = ', e10.4)
      write(14,2003) zzp
 2003 format('#    zp = ', e10.4)
c 2001 format('#'/'#   stopping power for (Zp) ', i2,' in (Zt) ',i2) 

      rh1 = 0d0
c
c ... loop over termperatures
c ---------------------------
      do 502 it=1,nmt
         ttt = temp(it)
c         write(13,1006) ttt
c 1006    format(/2x,'T = ', f8.1, ' eV')
c         write(13,1007)
c 1007    format(2x,'rho(g/cc)',4x,'E0(MeV)',
c     :          6x,'A1',8x,'A2',9x,'A3',9x,'A4',9x,'A5'
c     :          9x,'A6',9x,'A7',9x,'A8',9x,'A9',9x,'A10')
c
c .... compute Lindhard's stopping number
c ----------------------------------------
         tt0 = ttt
         ttt = max(ttt,0.025)

         rewind (12)
         icb = 1
         if (mloss .ge. 100) then
            icb = 0
            mloss = mod(mloss, 100)
         endif
         if (mloss .gt. 0) then
            xte0 = -1d31
            if (floss0(1:len0) .ne. floss(1:len2)) then
               open(unit=3, file=floss0(1:len0), status='old')
               call rloss(mloss)
               close(3)
               do i = 1, neee
                  do j = 1, nrho
                     cl0(i,j) = cl(i,j)
                  enddo
               enddo
               xte0 = xte
            endif            
            open(unit=3, file=floss(1:len2), status='old')
            call rloss(mloss)
            close(3)
            if (xte0 .lt. -1d30) xte0 = xte
            xt = log(ttt)
            if (xte0 .lt. xte) then
               if (xt .le. xte0) then
                  do i = 1, neee
                     do j = 1, nrho
                        cl(i,j) = cl0(i,j)
                     enddo
                  enddo     
               else if (xt .lt. xte) then
                  xf = (xt-xte0)/(xte-xte0)                  
                  do i = 1, neee
                     do j = 1, nrho
                        cl(i,j) = cl0(i,j)*(1-xf)+cl(i,j)*xf
                     enddo
                  enddo
               endif
            else if (xte0 .gt. xte) then
               if (xt .ge. xte0) then
                  do i = 1, neee
                     do j = 1, nrho
                        cl(i,j) = cl0(i,j)
                     enddo
                  enddo     
               else if (xt .gt. xte0) then
                  xf = (xt-xte0)/(xte-xte0)                  
                  do i = 1, neee
                     do j = 1, nrho
                        cl(i,j) = cl0(i,j)*(1-xf)+cl(i,j)*xf
                     enddo
                  enddo
               endif
            endif            
         else if (mloss .eq. 0) then
            call lindh(ttt)
         endif
c         do k = 1, neee
c            eion = eione(k)
c            do j = 1, nrho
c               write(11,3201) rl(j), cl(k,j)
c            enddo
c            write(11,3203)
      
c            if(k.eq.1) then
c               write(12,3300) ttt
c               write(12,3301)
c               write(12,3302) (eione(i),i=1,neee)
c               write(12,3303)
c               write(12,3202) (rl(i),i=1,nrho)
c            endif
c            write(12,3204) eion
c            write(12,3205) k
c            write(12,3202) (cl(k,i),i=1,nrho)
c         enddo               
c
c .... do loop over densities
c ----------------------------
         do 501 id=1,nmd
c
c ... readin radial mesh and electron charge dnesity function
c ------------------------------------------------------------
            read(2,*) t, d            
            read(2,*) nrd,h,rs
            read(2,*) (r(i),  i=1,nrd)
            read(2,*) (dr(i), i=1,nrd)
            read(2,*) (rho0(i),i=1,nrd)
            read(2,*) (rhof(i),i=1,nrd)
            read(2,*) (rhoe(i),i=1,nrd)
            chkt = abs(t-temp(it))
            chkd = abs(d-deni(id))
            if(chkt.gt.1.0e-3 .or. chkd.gt.1.0e-3) 
     :      stop 'input data not match !!!'

            do 5 i=1,nrd
               if(rho0(i).lt.0) rho0(i) = abs(rho0(i))
               if(r(i).ge.rs) then
                  nrd = i
                  rho0(i) = rho0(i-1)*(r(i)/r(i-1))**2
                  goto 6
               endif
 5          continue
 6          continue
 
c
c .... computing Zbar 
c ----------------------------------
            
            do 1101 i=1,nrd
               rho(i) = rhof(i)*dr(i)*r(i)
 1101       continue
            call intgrl(nrd,+h,rho,sum,rh1)
            zbar=cubint(rs,nrd,r,sum,nrd)
            write(14,1102) rs
 1102       format('#    rs = ',e15.8)
c            write(*,*) 'zb=', zbar, nrd, r(nrd), rs, rhof(nrd)
            do 100 ie=1,mep
               eion = ee(ie)
               vion = sqrt(2*eion*1.6021917e-6/pmass)
c
c ... effective charge for projectile ion
c  ---------------------------------------
               if(zzp.gt.2) then
                  beta = vion/3.0e10
                  factor = 137.04*beta/(zzp)**0.69
                  zeff = zzp*(1-1.034*dexpw(-factor))
               else
                  zeff = zzp
               endif
c
c ... compute rho(r)*L(rho,v)
c     ---------------------------
               qe0 = rhoe(nrd)/(4*pi*r(nrd)**2)
               ir0 = nrd
               ir1 = nrd
               zzr0 = 0.75
               zzr = rho0(1)*r(1)/3
               zz0 = rho0(1)
               qm = rho0(1)/(4*pi*r(1)**2)
               do k=2,nrd
                  zz1 = rho0(k)
                  q = zz1/(4*pi*r(k)**2)
                  if (q .gt. qm) qm = q
                  zzr0 = zzr
                  zzr = zzr + 0.5*(zz1+zz0)*(r(k)-r(k-1))
                  zz0 = zz1
                  if (zzr .ge. 0.75*ztg) then
                     ir0 = k-1
                     ir1 = k                     
                     exit
                  endif
               enddo
               rms = r(ir0) + (0.75*ztg-zzr0)*(r(ir1)-r(ir0))/(zzr-zzr0)
               qs = 6.748333e24*0.75*ztg/(4*pi*rms**3/3)
               zf = rho0(nrd)/(4*pi*r(nrd)**2)/qm
               zf = min(1.0, zf)
               zf = (1-zf)**4
               do 10 k=1,nrd
                  rho(k) = rho0(k)
                  fpr = 4*pi*r(k)**2
                  q = 6.748333e24*(rho(k))/fpr
                  if (icb .gt. 0) then
                     qe1 = qe0*(r(k)/r(nrd))**qei
                     qe = rhoe(k)/fpr-qe1
                     if (qe .gt. 0) then
                        qe = 6.748333e24*qea*qe
                     else
                        qe = 0d0
                     endif
                  else
                     qe = 0d0
                     qe1 = 0d0
                  endif
                  qf = 6.748333e24*rhof(k)/fpr
                  xcl = vlhfit(zeff, eion, q+qe, qe, qs, zf, ttt, mloss)
                  ycl = vlhfit(zeff, eion, qf, 0d0, qs, zf, ttt, mloss)
                  rhol(k) = rho(k)*xcl*r(k)*dr(k)
                  rhofl(k) = rhof(k)*ycl*r(k)*dr(k)
 10            continue

c
c ... do integration
c     ------------------
               call intgrl(nrd,+h,rhol,sum,rh1)
               sumitg = cubint(rs,nrd,r,sum,nrd)
               dedx(ie) = 4.58284e17*(zeff/vion)**2 * sumitg
               if (mod(ie-1,mre) .eq. 0) then
                  ire = 1 + (ie-1)/mre
                  do iir=1,nrd 
                     rdedx(ire,iir) = sum(iir)/sumitg
                     cdedx(ire,iir) = rhol(iir)*h
                  enddo
               endif
               fdedx(ie) = ee(ie)/(dedx(ie)*bcoef)
c               call intgrl(nrd,+h,rhofl,sum,rh1)
c               sumitg = cubint(rs,nrd,r,sum,nrd)
c               fdedx(ie) = 4.58284e17*(zeff/vion)**2 * sumitg            

 100        continue

            rh2 = 2*fdedx(1)
            call intgrl(mep,+he,fdedx,rage,rh2)

c     
c ... do 10 parameter fitting
c     ---------------------------
c            write(*,*) 'dexfit:', it, id
c            call dexfit(it,id,zzp, ztg, mep, ee, dedx,
c     :                                                a)

            do 30 i=1,mep
               eion = ee(i)
               if(eion.le.elow) then
                  fits(i) = a(1)*eion**a(2)
               endif
               if(eion.gt.elow .and.eion.lt.ehig) then
                  sl = a(3)*eion**a(4)
                  sh = (a(5)/eion)*log(1.0+a(6)/eion+a(7)*eion)
                  fits(i) = sl*sh/(sl+sh)
               endif
               if(eion.ge.ehig) then
                  fits(i) = (a(8)/eion)*log(1.0+a(9)/eion+a(10)*eion)
               endif
               fits(i) = 0.0
 30         continue

c
c .... output (1): to table for BUCKY
c -----------------------------------

c            write(13, 901) deni(id), epeak,  
c     :      a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), a(9), a(10)
c 901        format(1x,1p,e9.3,2x,11(1p,e10.3,1x),0p)

c
c .... output (2): compare data and fitting result
c ------------------------------------------------
            write(14,2110) temp(it)
            write(14,2111) deni(id)
            write(14,2112) zbar
            write(14,2113) mep
            write(14,2114) mloss
 2110       format('#    Te = ', e15.8)
 2111       format('#   rho = ', e15.8)
 2112       format('#  zbar = ', e15.8)
 2113       format('#   mep = ', i4)
 2114       format('# mloss = ', i3)
c 902        format('#'/'#  T = ',1p,e10.4,2x,'rho = ',1p,e10.4,0p,
c     :              '  zbar = ',f5.2)
            write(14,2115)
 2115       format('#',3x,'E/AMU (MeV)',12x,'dEdX',11x,'range') 
            do 31 i=1,mep
               dedxt(it,id,i) = dedx(i)
               write(14,903) ee(i), dedx(i), rage(i)
 31         continue
c            write(14,904)
 903        format(e15.8,1x,e15.8,1x,e15.8)
c 904        format('A')

            if (mout .gt. 0) then
               do ire=1,nre
                  do iir=1,nrd
                     write(15,905) temp(it),deni(id),ee(1+(ire-1)*mre),
     +                    cdedx(ire,iir),rdedx(ire,iir)
                  enddo
               enddo
            endif
 905        format(2x,5(1p,e15.8,3x))

c            call getrange(kzz,it,id,mep,ee,elog,dedx)

 501     continue
 502  continue

c      do 550 id=1,nmd
c      do 503 it=1,nmt
c         write(15,2101) temp(it), deni(id)
c         do 504 ie=1,mep
c            rag = (range(it,id,mep)-range(it,id,ie))/acoef
c            write(15,2002) ee(ie), rag, dedxt(it,id,ie)
c 504     continue
c         write(15,2003)
c 503  continue
c      write(15,2102)
c 550  continue

c      do 560 it=1,nmt
c      do 505 id=1,nmd
c         write(16,2004) deni(id), temp(it)
c         do 506 ie=1,mep
c            rag = (range(it,id,mep)-range(it,id,ie))/acoef
c            write(16,2002) ee(ie), rag, dedxt(it,id,ie)
c 506     continue
c         write(16,2003)
c 505  continue
c      write(16,2103)
c 560  continue

c 2101 format('###  temperature = ', 1p,e10.3,' density = ',1p,e9.3)
c 2002 format(3x,1p,e10.4,3x,1p,e10.4,3x,1p,e10.4)
c 2003 format('A')
c 2004 format('###  density = ', 1p,e10.3, ' T = ',1p,e9.3)
c 2102 format('### density cut boundary ###')
c 2103 format('### temperature cut boundary ###')

      
c 3201 format(5x,1pe12.5,3x,1pe12.5,3x,1pe12.5)
c 3202 format(5x,5(1pe12.5,1x), 4(5x,5(1pe12.5,1x)),
c     :     5x,4(1pe12.5,1x),1pe12.5)
c 3204 format('c     energy=', 1pe9.3) 
c 3203 format('A')
c 3205 format(6x,'data (vlhtab(',i2,',i),i=1,30)')
c 3300 format(5x,1p,e12.5)
c 3301 format('c     energy mesh (MeV)')
c 3302 format(5x,5(1pe12.5,1x),/8(5x,5(1pe12.5,1x)/),
c     :     5x,4(1pe12.5,1x),1pe12.5)
c     3303 format('c      log10(rho) mesh (#/cm**3)')
      close(14)
      if (mout .gt. 0) then
         close(15)
      endif
      end






c
c ******************************************************************* c
c                                                                     c
c                          L H F I T 
c   This routine is used for fitting Lindhard's stopping number in
c   the density range of 1.0e21 -- 1.0e33 (cm-3) and the energy 
c   range of 100 -- 5.0e-4 (MeV).
c
      function vlhfit(ze,e0,r0,re,rf,zf,te,md)
      implicit real*8 (a-h,o-z)
      common /lhtab/etab(50),rtab(30),vlhtab(50,30),vlhtab0(50,30),
     +     neee,nrho,xte,xte0

      ex = e0
      yb = 0d0
      rx = log10(max(1d-10,r0))
      if (r0 .gt. 0 .and.
     +     (md .lt. 0 .or. (md .ge. 10 .and. md .lt. 30))) then
         v2 = (2*e0*1e6/27.2114/1823)
         v = sqrt(v2)
         v3 = v*v2
         v4 = v2*v2
         xn = r0/6.748333e24
         xf = rf/6.748333e24
         wp = 3.545*sqrt(xn)
         wf = 3.545*sqrt(xf)
         xkf = 3.094*xf**(0.3333333333)
         y = 2*v2/wp
         if (y < 1d6) then
            y = 1/(1-exp(-1d0/y))
         endif
         y2 = y*y
         y8 = sqrt(1+8/y2)
         y28 = 0.5*y2*(1-y8)
         yx = 4*log(0.25*(1+y8))
c         xi = 0.0
         xi = re/r0
         xi2 = xi*xi
         yb0 = (5.-xi2/2.)*log(y)
         yb1 = (1.-xi2/4.)*(1+y28+yx)
         yb = (yb0+yb1)/3.0
         yb = yb * ze*wp*3.14159/v3
         ef = xkf*xkf*27.21
         vp = 4.0*(rf/1d25)**(-0.025)*sqrt(1+te/ef)
         vp4 = vp**4
         vf = (v4/(vp4+v4))*zf
         yb = yb * vf
         yb = max(yb, 0.0)
         if (md .ge. 20 .or. md .eq. -2) then
            y = ze/v
            y2 = y*y
            y4 = y2*y2
            vp4 = (2*vp)**4
            vf = (v4/(vp4+v4))*zf
            yc = -y2*1.202*vf
            yb = yb + yc
         endif
         if (md .lt. 0) then
            vlhfit = yb
            return
         endif
      endif
c
c .... (1) looking up the energy index
c
      if(ex.ge.etab(1)) then
         inde1 = 1
         inde2 = 1
         goto 10
      endif
      if(ex.le.etab(neee)) then
         inde1 = neee
         inde2 = neee
         goto 10
      endif
      do 1 ie=2,neee-1
         if(ex.eq.etab(ie)) then
            inde1 = ie
            inde2 = ie
            goto 10
         endif
 1    continue
      do 2 ie=1,neee-1
         if(ex.lt.etab(ie) .and. ex.gt.etab(ie+1)) then
            inde1 = ie
            inde2 = ie + 1
            goto 10
         endif
 2    continue
c
c ... looking up the density index
c
 10   if(rx.le.rtab(1)) then
         indr1 = 1
         indr2 = 1
         goto 20
      endif
      if(rx.ge.rtab(nrho)) then
         indr1 = nrho
         indr2 = nrho
         goto 20
      endif
      do 11 ir=2,nrho-1
         if(rx.eq.rtab(ir)) then
            indr1 = ir
            indr2 = ir
            goto 20
         endif
 11   continue
      do 12 ir=1,nrho-1
         if(rx.gt.rtab(ir) .and. rx.lt.rtab(ir+1)) then
            indr1 = ir
            indr2 = ir + 1
            goto 20
         endif
 12   continue

 20   continue
c
c .... if the point is right at one of the  nodes
c
      if(inde1.eq.inde2 .and. indr1.eq.indr2) then
         vlhfit = 10**vlhtab(inde1,indr1)+yb
         return
      endif

c
c .... do linear interolation 
c
      if(inde1.eq.inde2) then
         r1 = rtab(indr1)
         r2 = rtab(indr2)
         v1 = vlhtab(inde1,indr1)
         v2 = vlhtab(inde1,indr2)
         vlhfit = v1 + (v2-v1)/(r2-r1) * (rx-r1)         
         vlhfit = 10**vlhfit
c         vlhfit = cubint(rx, indr1, rtab, vlhtab(inde1,1:30), 30)
      else
c         v1 = cubint(rx, indr1, rtab, vlhtab(inde1, 1:30), 30)
c         v2 = cubint(rx, indr1, rtab, vlhtab(inde2, 1:30), 30)
         e1 = log10(etab(inde1))
         e2 = log10(etab(inde2))
         ex = log10(ex)
c         vlhfit = v1 + (v2-v1)/(e2-e1)*(ex-e1)
         r1 = rtab(indr1)
         r2 = rtab(indr2)
         u1 = vlhtab(inde1,indr1)
         u2 = vlhtab(inde2,indr1)
         u3 = vlhtab(inde1,indr2)
         u4 = vlhtab(inde2,indr2)
         v1 = u1 + (u2-u1)/(e2-e1) * (ex-e1)
         if(indr1.eq.indr2) then
            vlhfit = v1
            vlhfit = 10**vlhfit
         else
            v2 = u3 + (u4-u3)/(e2-e1) * (ex-e1)
            vlhfit = v1 + (v2-v1)/(r2-r1) * (rx-r1)
            vlhfit = 10**vlhfit
         endif
      endif
      vlhfit = vlhfit + yb
      return
      end
         
         
           

c ******************************************************************** c
c

      SUBROUTINE INTGRL(N,H,Y,F,FIN)
 
C       F(X) = F(X1) + INTEGRAL(X1,X;Y(T))           H>0
C       F(X) = INTEGRAL(X,X2;Y(T)) + F(X2)           H<0
C    SEE ABRAMOWITZ: P 886, (25.4.12); P 915, TABLE (25.3)
C    THIS SUBROUTINE TAKES ADVANTAGE OF TH VECTOR PROPERTIES OF THE
C    CRAY 1 COMPUTER. TWO STATEMENTS BELOW PREFIXED BY CDIR$ ARE
C    COMPILER DIRECTIVES. THE EQUIVALENT SUBROUTINE IN GOLIATH MAY BE
C    SUBSTITUTED FOR THIS ONE.
 
      PARAMETER  (NRX=850,IX=35,JX=6*IX,LX=1,NX=NRX+2)
      implicit real*8 (a-h,o-z)
      REAL*8  Y(N),F(N),G(NX)
 
      IF (N .LT. 1)  GO TO 7
      N1 = N - 1
      H24 = ABS(H)/24.0
      IF (H .LT. 0.0)  GO TO 4
      IF (N .LT. 4)  GO TO 3
      G(2) = H24*( 9.0*Y(1)+19.0*Y(2)-5.0*Y(3)+Y(4) )
      DO  1  K=3,N1
 1    G(K) = H24*( -Y(K-2)+13.0*(Y(K-1)+Y(K))-Y(K+1) )
      G(N) = H24*( Y(N-3)-5.0*Y(N-2)+19.0*Y(N-1)+9.0*Y(N) )
      NB = 4*((N-1)/4) + 1
      NC = NB + 1
      F(1) = FIN
      DO  21  K=5,NB,4
 21   F(K) = F(K-4) + G(K-3) + G(K-2) + G(K-1) + G(K)
CDIR$  IVDEP
      DO  22  K=5,NB,4
      F(K-3) = F(K-4) + G(K-3)
      F(K-2) = F(K-4) + G(K-3) + G(K-2)
 22   F(K-1) = F(K-4) + G(K-3) + G(K-2) + G(K-1)
      IF (NB .EQ. N)  GO TO 24
      DO  23  K=NC,N
 23   F(K) = F(K-1) + G(K)
 24   CONTINUE
      RETURN
 3    F(1) = FIN
      IF (N .EQ. 1)  RETURN
      HA = ABS(H)
      F(2) = FIN + HA*(Y(1)+Y(2))/2.0
      IF (N .EQ. 2) RETURN
      F(2) = F(2) - HA*(Y(1)-2.0*Y(2)+Y(N-2))/12.0
      F(3) = FIN + HA*(Y(1)+4.0*Y(2)+Y(3))/3.0
      RETURN
 
 4    IF (N .LT. 4)  GO TO 6
      NP = N + 1
      G(N-1) = H24*( 9.0*Y(N)+19.0*Y(N-1)-5.0*Y(N-2)+Y(N-3) )
      DO  5  K=3,N1
      J = NP - K
 5    G(J) = H24*( -Y(J+2)+13.0*(Y(J+1)+Y(J))-Y(J-1) )
      G(1) = H24*( Y(4)-5.0*Y(3)+19.0*Y(2)+9.0*Y(1) )
      NB = 4*((N-1)/4) + 1
      NC = NB + 1
      F(N) = FIN
      DO  31  K=5,NB,4
      J = NP - K
 31   F(J) = F(J+4) + G(J+3) + G(J+2) + G(J+1) + G(J)
CDIR$  IVDEP
      DO  32  K=5,NB,4
      J = NP - K
      F(J+3) = F(J+4) + G(J+3)
      F(J+2) = F(J+4) + G(J+3) + G(J+2)
 32   F(J+1) = F(J+4) + G(J+3) + G(J+2) + G(J+1)
      IF (NB .EQ. N)  GO TO 34
      DO  33  K=NC,N
      J = NP - K
 33   F(J) = F(J+1) + G(J)
 34   CONTINUE
      RETURN
 6    F(N) = FIN
      IF (N .EQ. 1)  RETURN
      HA = ABS(H)
      F(N-1) = FIN + HA*(Y(N)+Y(N-1))/2.0
      IF (N .EQ. 2) RETURN
      F(N-1) = FIN - HA*(Y(N)-2.0*Y(N-1)+Y(N-2))/12.0
      F(N-2) = F(N-1) + HA*(Y(N)+4.0*Y(N-1)+Y(N-2))/3.0
      RETURN
 
 7    PRINT 8,  N
      STOP  7
 8    FORMAT('       STOP IN INTGRL. N<1 .  N =',I3)
      END
                          SUBROUTINE ONE(Y,R,DR,G,S)
      implicit real*8 (a-h,o-z)
      REAL*8  F(2),R(2),DR(2),Y(2)
 
      F(1) = Y(1)/DR(1)/R(1)
      F(2) = Y(2)/DR(2)/R(2)
      IF (F(1)*F(2) .LE. 0.0)  GO TO 1
      S = DLOG(F(2)/F(1))/DLOG(R(2)/R(1)) + 1.0
      G = F(1)*R(1)/S
      RETURN
 
 1    S = 1.0
      G = 0.0
      RETURN
      END



      subroutine rloss(md0)
      implicit real*8 (a-h,o-z)

      common /lhtab/eione(50),rl(30),cl(50,30),cl0(50,30),
     +     neee,nrho,xte,xte0

      md = mod(md0,10)
      if (md .eq. 5) md = md0
      
      read(3, 101) nttt
      read(3, 102) nrho
      read(3, 103) neee
c      write(*,*) nttt, nrho, neee
      do i = 1, nttt
         do j = 1, nrho
            read(3, 104) ii, ttt
            if (i .eq. 1 .and. j .eq. 1) then
               if (ttt .gt. 0) then
                  xte = log(ttt)
               else
                  xte = -1d2
               endif
            endif
            read(3, 105) jj, rhoj
            rl(j) = log10(rhoj)
            do k = 1, neee
               read(3, *) ttti, rhoj, ek, cl1, cl2, cl3, cl4, cl5, cl6
c               read(3, *) ttti, rhoj, ek, cl1, cl2, cl3
               kk = neee-k+1
               eione(kk) = ek               
               if (md .eq. 1) then
                  cl(kk,j) = cl1
               else if (md .eq. 2) then
                  cl(kk,j) = cl2
               else if (md .eq. 3) then
                  cl(kk,j) = cl3
               else if (md .eq. 4) then
                  cl(kk,j) = cl4
               else if (md .eq. 5) then
                  cl(kk,j) = cl5
               else if (md .eq. 15) then
                  cl(kk,j) = cl1+cl5
               else if (md .eq. 25) then
                  cl(kk,j) = cl2+cl5
               else if (md .eq. 35) then
                  cl(kk,j) = cl3+cl5
               else if (md .eq. 45) then
                  cl(kk,j) = cl4+cl5
               endif
            enddo
         enddo
      enddo
      do k = 1, neee
         do j = 1, nrho
            cl(k,j) = log10(max(1d-50,cl(k,j)))
         enddo
      enddo
      
 101  format(7x, I2)
 102  format(7x, I2)
 103  format(7x, I2)     
 104  format(6x, I2, 1x, 1pe12.5)
 105  format(6x, I2, 1x, 1pe12.5) 
 106  format(5(1pe12.5,1x),1pe12.5)
      end
      
c ******************************************************************* c
c                                                                     c
c                                                                     c
c                       L I N D H
c  This routine evaluates Lindhard's stopping number of ions with
c  any velocity in dense plasmas at any temperatures
c

      subroutine lindh(ttt)
      implicit real*8 (a-h,o-z)
      common /lhtab/eione(50),rl(30),cl(50,30),cl0(50,30),
     +     neee,nrho,xte,xte0
      dimension  rho(30)


      pi = 3.1415926
c     electron mass in g
      emass = 9.109389700e-28
c     electron charge in cgs unit
      e  = 4.8065e-10
c     proton mass in g
      pmass = 1.672623100e-24
c     Plank constant h/2pi erg*s
      hbar = 1.0545887e-27
c     Bohr velocity (cm/s)
      v0   = e**2/hbar
c     Bohr radius (cm)
      a0   = 5.291772e-9


      neee = 50            
      eione(1) = 100
      eione(2) = 75
      eione(3) = 50
      eione(4) = 25
      eione(5) = 10
      de = -1.0/9
      j = 5
 900  do 901 i=2,46
         j = j+1
         eione(j) = eione(j-1)*10**(de)
 901  continue

c
c ... constructing density mesh
c     
      rho1 = 1.0e19
      rho2 = 1.0e29
      nrho = 30

      rho(1) = rho1
      if(nrho.gt.1) then
         drho = log(rho2/rho1)/(nrho-1)
         do 2 k=2,nrho
            rho(k) = rho(k-1)*dexpw(drho)
 2       continue
      endif


      do 100 iee=1,neee
         eion = eione(iee)         
         do 50 irho=1,nrho
            rrho = rho(irho)

c
c ... define some variables
c --------------------------
c     inter-electron distance rs (in a0)
            rs = 1./(4./3*pi*rrho)**(1./3.)/a0
c     fermi energy in (2Ry)
            ef = 1.84/rs**2
c     ion velocity in (cm/s)
            vion = sqrt(2*eion*1.6021917e-6/pmass)
c     fermi velocity in (cm/s)
            vf = (hbar/emass) * (3*pi**2*rrho)**(1./3.)
c     reduce temperature
            te = ttt/(27.2116*ef)

c                                                                      c
c ************************ generating table  ************************* c
c     c
            if (ttt .le. 1e-10) then
               call lind0(eion, rrho, vlh2)
            else
               call lindt(eion,rrho,ttt,vlh2)
            endif
            rl(irho)     = log10(rho(irho))
            cl(iee,irho) = vlh2
            cl(iee,irho) = log10(max(1d-50,vlh2))

 50   continue

 100  continue      
      
      return
      end


c ************************************************************************
c 
c                     L I N D H A R D
c    This routine evaluates zero temperature Lindhard's stopping number
c
c    reference: J. Appl. Phys. 50, 5579, 1979.
c ************************************************************************

      subroutine lind0(eion,rho,cl)
      implicit real*8 (a-h,o-z)
      
      dimension uu(1000), up(100), wup(100)
      dimension aa(1000), ap(100), wap(100)


c
c ... constants
c
      pi = 3.1415926
c     electron mass in g
      emass = 9.109389700e-28
c     proton mass in g
      pmass = 1.672623100e-24
c     Plank constant h/2pi erg*s
      hbar = 1.0545887e-27
c     Bohr velocity (cm/s)
      v0   = 2.18769066e8
c     electron charge in cgs unit
      e  = 4.8065e-10
c
c ... Fermi velocity
c
      vf = (hbar/emass) * (3*pi**2*rho)**(1./3.)
      x2  = v0/pi/vf


c
c ... incident ion velocity (cm/s)
c                          
      vion = sqrt(2*eion*1.6021917e-6/pmass)

c 
c ... plasma frequency
c
      wp = sqrt(4*pi*rho*e**2/emass)



c            
c ... dividing the first integration range 0 to v/vf into 
c     ten equal logritham sub-ranges
c

         vratio = vion/vf
         xx = x2


c
c .... high velocity cases
c
         if(vratio.gt.5) then
            v1 = log(2*emass*vion**2/(hbar*wp))
            sumu = pi*xx/6*v1
            cl = 6/pi/xx*sumu 
            return
         endif

c 
c .... bridging point
c
         vbp = 5.0*vf
         v1 = log(2*emass*vbp**2/(hbar*wp))
         sumu = pi*xx/6*v1     
         cl1  = 6/pi/xx*sumu      
         
         icheck = 1
         ratio1 = vratio
         vratio = 5
c
c .... for general cases
c

 1234    if(vratio .gt. 1) then            
            uu(1) = 0
            uu(2) = 1
            uu(3) = vratio
            nu = 3
            np = 50
         else
            nu = 2
            uu(1) = 0
            uu(2) = vratio
            np = 100
         endif
         
c
c ... outer integral over u
c
         sumu = 0.0
         do 50 iu=1,nu-1
            u1 = uu(iu)
            u2 = uu(iu+1)
            
            call nodewt(u1,u2,up,wup,np)
            
            sumu1 = 0.0
            do 40 ju=1,np
               u = up(ju)
c
c ... searching for the root of function g
c
               if(u.gt.1) call groot(u, xx,a0)
               
c
c ... inner integral over a
c     the whole range -u to infinity is divided into 3 parts:
c     (1)-u to -1; (2) -1 to 1; (3) 1 to infinity
c
               if(u.gt.1) then
                  naa = 4
                  aa(1) = -u
                  aa(2) = -1
                  aa(3) = 0
                  aa(4) = 1
               else
                  naa = 3
                  aa(1) = -u
                  aa(2) = 0
                  aa(3) = 1
               endif
               
               suma = 0.0
               
               do 20 ia=1,naa-1
                  a1 = aa(ia)
                  a2 = aa(ia+1)
                  
                  if(abs(a1).gt.1.0 .or. abs(a2).gt.1.0) then
                     if(abs(a0).gt.1) then
                        suma=suma+(u+a0)**3*ff0(xx,u,a0)
                     endif
                  else
                     call nodewt(a1,a2,ap,wap,100)                     
                     suma1 = 0.0
                     do 10 j=1,100
                        a = ap(j)
                        suma1 = suma1+ wap(j)*(u+a)**3 * ff(xx,u,a)
 10                  continue
                     suma = suma + suma1
                  endif
 20            continue
               
               sumu1 = sumu1 + wup(ju)*u*suma
               
 40         continue
            
            sumu = sumu + sumu1
            
 50      continue
         

         cl = 6/pi/xx*sumu
         if(icheck.eq.1) then
            cg = (cl/cl1-1.)/vbp**2
            vratio = ratio1
            icheck = 2
            goto 1234         
         else
            cl = cl/(1+cg*vion**2)
         endif

      return
      end

                  
c
c ********************************************************************* c
c              This routine define function F

      function ff(x2,u,a)
      implicit real*8 (a-h,o-z)

      pi = 3.1415926

      ff = x2*f2(u,a)/(g(a,u,x2)**2 + (x2*f2(u,a))**2)

      return
      end


      function f1(u,a)
      implicit real*8 (a-h,o-z)

c     
c ... define f1
c      
      z = u + a
      al1 = abs((z-u+1)/(z-u-1))
      al2 = abs((z+u+1)/(z+u-1))
      f1 = 0.5 + 1./(8*z) * (1-(z-u)**2) * log(al1)
     :         + 1./(8*z) * (1-(z+u)**2) * log(al2)


      return
      end

      function f2(u,a)
      implicit real*8 (a-h,o-z)

      pi = 3.1415926
c
c ... define f2
c
      if((2*u+a) .lt. 1) then
         f2 = 0.5 * pi * u
      endif
      if(abs(a).lt.1 .and. (2*u+a).gt.1) then
         f2 = pi*(1-a**2)/(8*(u+a))
      endif
      if(abs(a).gt.1) then
         f2 = 0.
      endif

      return
      end


      function g(a,u,x2)
      implicit real*8 (a-h,o-z)
c
c ... define g
c
      g = (u+a)**2 + x2*f1(u,a)

      return
      end

      function dgda(a,u,x2)
      implicit real*8 (a-h,o-z)


      z = u+a
      zuu1 = z-u+1
      zub1 = z-u-1
      zuu2 = z+u+1
      zub2 = z+u-1

      dgda1 = 2*z
      pm = 1.
      if(a.lt.-1.0 .or. a.gt.1.0) pm=-1
      dgda2 = -1./(8*z**2) * (1-(z-u)**2) * log(abs(zuu1/zub1))
     :        -2./(8*z)*(z-u) * log(abs(zuu1/zub1))
     :        +1./(8*z) * (1-(z-u)**2) 
     :        *abs(zub1/zuu1) * (pm*2)/zub1**2

      pm = 1
      b = a + 2*u
      if(b.lt.-1.0 .or. b.gt.1.0) pm=-1
      dgda3 = -1./(8*z**2) * (1-(z+u)**2) * log(abs(zuu2/zub2))
     :        -2./(8*z)*(z+u) * log(abs(zuu2/zub2))
     :        +1./(8*z) * (1-(z+u)**2) 
     :        *abs(zub2/zuu2) * (pm*2)/zub2**2
      
      dgda = dgda1 + x2*(dgda2+dgda3)

      return
      end



      function ff0(x2,u,a0)
      implicit real*8 (a-h,o-z)

      pi = 3.1415926

      ff0 = pi/abs(dgda(a0,u,x2))

      return
      end
      



c
c ***************************************************************** c
c      searching for the root of the g(u,a) for a given u value


      subroutine groot(u,xx, rtsafe)
      implicit real*8 (a-h,o-z)

      x1 = 0
      x2 = -u

      fl = g(x1, u, xx)
      fh = g(x2, u, xx)

      if(fl*fh.gt.0)  stop 'searching a0 fail !'

      if(fl.eq.0.)then
        rtsafe=x1
        return
      else if(fh.eq.0.)then
        rtsafe=x2
        return
      else if(fl.lt.0.)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
      rtsafe=.5*(x1+x2)
      dxold=abs(x2-x1)
      dx=dxold
      f = g(rtsafe,u, xx)
      df= dgda(rtsafe,u, xx)
      do 11 j=1,100
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0. .or. abs(2*f)
     : .gt.abs(dxold*df) ) then
          dxold=dx
          dx=0.5*(xh-xl)
          rtsafe=xl+dx
          if(xl.eq.rtsafe)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp.eq.rtsafe)return
        endif
        if(abs(dx).lt.1.0e-7) return
        f = g(rtsafe,u,xx)
        df = dgda(rtsafe,u,xx)
        if(f.lt.0.) then
          xl=rtsafe
        else
          xh=rtsafe
        endif
11    continue
      pause 'rtsafe exceeding maximum iterations'
      return
      end





c ************************************************************
c This data block contains the abscissas and weight factors
c  for gaussian integration
c

      subroutine nodewt(x1,x2,x,w,np)
      implicit real*8 (a-h,o-z)
      common/nodes/ x30(15), w30(15), x50(25), w50(25), 
     :              x80(40), w80(40), x100(50), w100(50)
      dimension x(np), w(np)

      
      bma = (x2-x1)/2.
      bpa = (x2+x1)/2.
      if(np.eq.30) then
         do 1 i=1,15
            x(i)    = -bma*x30(16-i) + bpa
            x(15+i) =  bma*x30(i)    + bpa
            w(i)    =  bma*w30(16-i)
            w(15+i) =  bma*w30(i)
 1       continue
         return
       endif

      if(np.eq.50) then
         do 2 i=1,25
            x(i)    = -bma*x50(26-i) + bpa
            x(25+i) =  bma*x50(i)    + bpa
            w(i)    =  bma*w50(26-i)
            w(25+i) =  bma*w50(i)
 2       continue
         return
       endif

      if(np.eq.80) then
         do 3 i=1,40
            x(i)    = -bma*x80(41-i) + bpa
            x(40+i) =  bma*x80(i)    + bpa
            w(i)    =  bma*w80(41-i)
            w(40+i) =  bma*w80(i)
 3          continue
         return
       endif

      if(np.eq.100) then
         do 4 i=1,50
            x(i)    = -bma*x100(51-i) + bpa
            x(50+i) =  bma*x100(i)    + bpa
            w(i)    =  bma*w100(51-i)
            w(50+i) =  bma*w100(i)
 4       continue
       endif
       
       return
       end

      block data
      implicit real*8 (a-h,o-z)
      common/nodes/ x30(15), w30(15), x50(25), w50(25), 
     :              x80(40), w80(40), x100(50), w100(50)

      data x30/
     :   .0514718426,  .1538699136,  .2546369262,  .3527047255,
     :    .4470337695,  .5366241481,  .6205261830,  .6978504948,
     :    .7677774321,  .8295657624,  .8825605358,  .9262000474,
     :    .9600218650,  .9836681233,  .9968934841/
      data w30/
     :    .1028526529,  .1017623897,  .0995934206,  .0963687372,
     :    .0921225222,  .0868997872,  .0807558952,  .0737559747,
     :    .0659742299,  .0574931562,  .0484026728,  .0387991926,
     :    .0287847079,  .0184664683,  .0079681925/
      data x50/
     :    .0310983383,  .0931747016,  .1548905900,  .2160072369,
     :    .2762881938,  .3355002454,  .3934143119,  .4498063350,
     :    .5044581449,  .5571583045,  .6077029272,  .6558964657,
     :    .7015524687,  .7444943022,  .7845558329,  .8215820709,
     :    .8554297694,  .8859679795,  .9130785567,  .9366566189,
     :    .9566109552,  .9728643851,  .9853540840,  .9940319694,
     :    .9988664044/
      data w50/
     :    .0621766167,  .0619360674,  .0614558996,  .0607379708,
     :    .0597850587,  .0586008498,  .0571899256,  .0555577448,
     :    .0537106219,  .0516557031,  .0494009384,  .0469550513,
     :    .0443275043,  .0415284631,  .0385687566,  .0354598356,
     :    .0322137282,  .0288429936,  .0253606736,  .0217802432,
     :    .0181155607,  .0143808228,  .0105905484,  .0067597992,
     :    .0029086226/
      data x80/
     :    .0195113833,  .0585044372,  .0974083984,  .1361640228,
     :    .1747122918,  .2129945029,  .2509523584,  .2885280549,
     :    .3256643707,  .3623047535,  .3983934059,  .4338753708,
     :    .4686966152,  .5028041119,  .5361459209,  .5686712681,
     :    .6003306228,  .6310757730,  .6608598990,  .6896376443,
     :    .7173651854,  .7440002976,  .7695024201,  .7938327175,
     :    .8169541387,  .8388314736,  .8594314067,  .8787225677,
     :    .8966755794,  .9132631026,  .9284598772,  .9422427613,
     :    .9545907663,  .9654850890,  .9749091406,  .9828485727,
     :    .9892913025,  .9942275410,  .9976498644,  .9995538227/
      data w80/  
     :    .0390178137,  .0389583960,  .0388396511,  .0386617598,
     :    .0384249930,  .0381297113,  .0377763644,  .0373654902,
     :    .0368977146,  .0363737499,  .0357943940,  .0351605290,
     :    .0344731205,  .0337332150,  .0329419394,  .0321004987,
     :    .0312101742,  .0302723218,  .0292883696,  .0282598161,
     :    .0271882275,  .0260752358,  .0249225358,  .0237318829,
     :    .0225050902,  .0212440261,  .0199506109,  .0186268142,
     :    .0172746521,  .0158961836,  .0144935080,  .0130687616,
     :    .0116241141,  .0101617660,  .0086839453,  .0071929048,
     :    .0056909225,  .0041803131,  .0026635336,  .0011449500/
      data x100/
     :    .0156289844,  .0468716824,  .0780685828,  .1091892036,
     :    .1402031372,  .1710800805,  .2017898641,  .2323024818,
     :    .2625881204,  .2926171880,  .3223603439,  .3517885264,
     :    .3808729816,  .4095852917,  .4378974022,  .4657816498,
     :    .4932107892,  .5201580199,  .5465970121,  .5725019326,
     :    .5978474702,  .6226088602,  .6467619085,  .6702830156,
     :    .6931491994,  .7153381176,  .7368280898,  .7575981185,
     :    .7776279096,  .7968978924,  .8153892383,  .8330838799,
     :    .8499645279,  .8660146885,  .8812186794,  .8955616450,
     :    .9090295710,  .9216092981,  .9332885350,  .9440558701,
     :    .9539007829,  .9628136543,  .9707857758,  .9778093585,
     :    .9838775407,  .9889843952,  .9931249370,  .9962951347,
     :    .9984919506,  .9997137268/
      data w100/
     :    .0312554235,  .0312248843,  .0311638357,  .0310723374,
     :    .0309504789,  .0307983790,  .0306161866,  .0304040795,
     :    .0301622651,  .0298909796,  .0295904881,  .0292610841,
     :    .0289030896,  .0285168543,  .0281027557,  .0276611982,
     :    .0271926134,  .0266974592,  .0261762192,  .0256294029,
     :    .0250575445,  .0244612027,  .0238409603,  .0231974232,
     :    .0225312203,  .0218430024,  .0211334421,  .0204032326,
     :    .0196530875,  .0188837396,  .0180959407,  .0172904606,
     :    .0164680862,  .0156296211,  .0147758845,  .0139077107,
     :    .0130259479,  .0121314577,  .0112251140,  .0103078026,
     :    .0093804197,  .0084438715,  .0074990733,  .0065469485,
     :    .0055884280,  .0046244501,  .0036559612,  .0026839254,
     :    .0017093927,  .0007346345/



       end


c ********************************************************************* c
c                                                                       c
c                       T L I N D H A R D
c
c   This routine calculates temperature dependent Lindhard's stopping
c   number.
c
c                  reference: Phys. Rev. A26, 665(1982).
c                             J. De. Phys. Vol 46, 71(1985)
c ********************************************************************* c

 
      subroutine lindt(eion, rho, 
     :                           t, vlh)
      implicit real*8 (a-h,o-z)
      dimension vii(10), vlhp(10)

c 
c .... constants
c
      zero = 0.0
      f12  = 1./2.
      f32  = 3./2.
      f52  = 5./2.
      pi = 3.1415926
c     electron mass in g
      emass = 9.109389700e-28
c     electron charge in cgs unit
      e  = 4.8065e-10
c     proton mass in g
      pmass = 1.672623100e-24
c     Plank constant h/2pi erg*s
      hbar = 1.0545887e-27
c     Bohr velocity (cm/s)
      v0   = e**2/hbar
c     Bohr radius (cm)
      a0   = 5.291772e-9

c
c ... define osme variables
c --------------------------
c     inter-electron distance rs (in a0)
      rs = 1./(4./3*pi*rho)**(1./3.)/a0
c     X2
      x2 = 0.5211*rs/pi
c     fermi energy in (2Ry)
      ef = 1.84/rs**2
c     ion velocity in (cm/s)
      vion = sqrt(2*eion*1.6021917e-6/pmass)
c     fermi velocity in (cm/s)
      vf = (hbar/emass) * (3*pi**2*rho)**(1./3.)
c     reduce temperature
      te = t/(27.2116*ef)
c     plasma frequency
      wp = sqrt(4*pi*rho*e**2/emass)


c
c ... chemical potential
c
      call findmu(te,alpha)
      amu = alpha*t/27.2116
c     reduce chemical potential
      gamma = amu/ef

c
c  determine the cut point between the 'high' and low': Vint
c -----------------------------------------------------------

      ve2 = vf**2 * te * fd(f32,alpha)/fd(f12,alpha)
      cnt = hbar*wp/emass
      vint= sqrt(1.5*(ve2 + cnt))
      vratio = vion/vint

      call highv(te,alpha,vint,vf,emass,hbar,wp, vlhh)
      call lowv (te,gamma,alpha,vint,vf,x2, vlhl)             

      cg = (vlhl/vlhh - 1.)/vint**2
c
c  for high velocity  limit
c -------------------------
      if(vion.gt.1.5*vint) then
         call highv(te,alpha,vion,vf,emass,hbar,wp, vlh)
         goto 1234
       endif

c
c ... for low velocity case
c
        
        if(vion.le.0.75*vint) then
           call lowv(te,gamma,alpha,vion,vf,x2,vlh)       
           vlh = vlh/(1+cg*vion**2)
           goto 1234
        endif


c
c  for general cases, do interpolation
c -------------------------------------

        vii(1) =  0.25*vint
        vii(2) =  0.50*vint
        vii(3) =  0.75*vint
        vii(4) =  1.50*vint
        vii(5) =  2.00*vint
        vii(6) =  2.50*vint

        do 1 i=1,3
           vionp = vii(i)
           call lowv(te,gamma,alpha,vionp,vf,x2,vlh)
           vlh = vlh/(1+cg*vionp**2)
           vlhp(i) = log10(vlh)
           vii(i) = log10(vii(i))           
 1      continue

        do 2 i=4,6
 3         vionp = vii(i)
           call highv(te,alpha,vionp,vf,emass,hbar,wp, vlh)       
           vlhp(i) = log10(vlh)
           vii(i) = log10(vii(i)) 
 2      continue

        vpoit = log10(vion)
        call itpolt(6,vii,vlhp,vpoit,vlh)
        vlh = 10**(vlh)

 1234  continue
      

       return
       end


       function fff(u,z,x2,gamma,te)
       implicit real*8 (a-h,o-z)

       fz1 = z**2 + x2*ft1(u,z,te,gamma)
       fz2 = x2*ft2(u,z,te,gamma)

       fff = z**3*fz2/(fz1**2 + fz2**2) 

       return
       end




      function ft2(u,z,te,gamma)
      implicit real*8 (a-h,o-z)

      pi = 3.1415926

      pp = u+z
      pm = u-z
      
      ex1 = (gamma-pp**2)/te
      ex2 = (gamma-pm**2)/te
      if(ex1 .le. 10.) then
         eup = log(1+dexpw(ex1))
      else
         eup = ex1
      endif
      if(ex2 .le. 10.) then
         elw = log(1+dexpw(ex2))
      else
         elw = ex2
      endif

      ft2 = -pi*te/(8*z) * (eup-elw)

      return
      end

      function ft1(u,z,te,gamma)
      implicit real*8 (a-h,o-z)

      common/vft0v/sumft0
      dimension q0m(100), w0m(100)

      an(x) = 0.7071*dsqrt( r + sqrt(r**2+(2*x+1)**2*pi**2*te**2) )
      bn(x) = 0.7071*dsqrt(-r + sqrt(r**2+(2*x+1)**2*pi**2*te**2) )


      pi = 3.1415926      
      r  = gamma
 
c
c ... do summation for the infinit convergent series
c

      sumv1 = 0.0
      sumv2 = 0.0
c ... n=0 point
      n = 0
      x = n 
      rn = an(x)**2 + bn(x)**2
      pp = ((u+z) + an(x))/bn(x)
      pm = ((u+z) - an(x))/bn(x)
      qp = ((u-z) + an(x))/bn(x)
      qm = ((u-z) - an(x))/bn(x)         
      v1 = bn(x)/rn
      v2 = atan(pp) + atan(pm) - atan(qp) - atan(qm)
      sumv1 = sumv1 + v1
      sumv2 = sumv2 + v2
c
c .... for V1
c -----------
c
c ... integrate over 0 to infinit
      v11 = 0
      v21 = 0
      q1 = 0
      q2 = 10
      do 5 in=1,10
         call nodewt(q1,q2,q0m,w0m,30)
      do 1 i=1,30
         x = q0m(i)
         rn = an(x)**2 + bn(x)**2
         v11 = v11 + w0m(i)*bn(x)/rn

         pp = ((u+z) + an(x))/bn(x)
         pm = ((u+z) - an(x))/bn(x)
         qp = ((u-z) + an(x))/bn(x)
         qm = ((u-z) - an(x))/bn(x)                  
         atgs = atan(pp) + atan(pm) - atan(qp) - atan(qm)
         v21 = v21 + w0m(i)*atgs
 1    continue
      q1 = q2 
      q2 = 10*q2
 5    continue

c
c ... function value at 0 
      v12 = 0.5*v1
c
c ... derivative value at 0      
      f0 = v1
      dx = 1.0e-3
      rx = an(dx)**2 + bn(dx)**2
      fx = bn(dx)/rx
      v13 = 1./12 * (fx-f0)/dx

      sumv1 = sumv1 + v11 - v12 + v13

c .... for V2
c ------------

      v22 = 0.5*v2

      f0 = v2
      dx = 1.0e-3
      pp = ((u+z) + an(dx))/bn(dx)
      pm = ((u+z) - an(dx))/bn(dx)
      qp = ((u-z) + an(dx))/bn(dx)
      qm = ((u-z) - an(dx))/bn(dx)                  
      fx = atan(pp) + atan(pm) - atan(qp) - atan(qm)
      v23 = 1./12 * (fx-f0)/dx

      sumv2 = 0.25/z*(sumv2 + v21 - v22 + v23)

      sum = pi*te*(sumv1 - sumv2)

 111  format('A')
 112  format(10(1pe12.4,2x))


      ft1 = sumft0 + sum

      return
      end



      subroutine coeft0(te,gamma)
      implicit real*8 (a-h,o-z)

      common/vft0v/sumft0

      x = -0.5
      y = gamma/te
      sumft0 = sqrt(te)/2*fd(x,y)
      
       return
       end

      function ft0(te,gamma)
      implicit real*8 (a-h,o-z)

      dimension qn(100), qw(100)
      
      q(x)   = (x*te+gamma)

      np = 30

      sum = 0.0
      q1 = 0
      xk = 1.0
      q2 = q(xk)
      if(q2.gt.0) then
         q2 = sqrt(q2)
      else
         q2 =1       
      endif
 1    call nodewt(q1,q2,qn,qw,np)
      sum1 = 0.0
      do 2 i=1,np
         sum1 = sum1 + qw(i)*an0(qn(i),gamma,te)
 2    continue
      sum = sum + sum1
      
      ratio = abs(sum1/sum)
      if(ratio.lt.1.0e-2) goto 3
      q1 = q2
      xk = 10*xk
      q2 = q(xk)
      if(q2.gt.0) then
         q2 = sqrt(q2)
      else
         q2 = 5*q1
      endif
      goto 1

 3    ft0 = sum
      return
      end

      function an0(x,gamma,te)
      implicit real*8 (a-h,o-z)

      y = (x**2 - gamma)/te
      if(y.gt.10) then
         an0 = dexpw(-y)
      else
         an0 = 1./(1+dexpw(y))
      endif
      
      return
      end

c
c ***************************************************************** c
c      searching for the root of the 
c
c           2/3*Te**(-2./3) - fd(0.5,alpha) = 0
c      for 'alpha'


      subroutine findmu(zz,rtsafe)
      implicit real*8 (a-h,o-z)

      itry = 1
      dmu = 1.0e-3
      x1 = -10
      x2 = 10

 123  fl = fzz(zz,x1)
      fh = fzz(zz,x2)
      
      if(fl*fh.gt.0) then
         x1=10*x1
         x2=10*x2
         itry=itry+1         
         if(itry.gt.10)  stop 'searching a0 fail at findmu!'
         goto 123
      endif


      if(fl.eq.0.)then
        rtsafe=x1
        return
      else if(fh.eq.0.)then
        rtsafe=x2
        return
      else if(fl.lt.0.)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
      rtsafe=.5*(x1+x2)
      dxold=abs(x2-x1)
      dx=dxold
      f = fzz(zz,rtsafe)
      rp =rtsafe+dmu
      fp= fzz(zz,rp)
      df= (fp-f)/dmu
      do 11 j=1,100
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0. .or. abs(2*f)
     : .gt.abs(dxold*df) ) then
          dxold=dx
          dx=0.5*(xh-xl)
          rtsafe=xl+dx
          if(xl.eq.rtsafe)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp.eq.rtsafe)return
        endif
        if(abs(dx).lt.1.0e-5) return
        f = fzz(zz,rtsafe)
        rp =rtsafe+dmu
        fp= fzz(zz,rp)
        df= (fp-f)/dmu
        if(f.lt.0.) then
          xl=rtsafe
        else
          xh=rtsafe
        endif
11    continue
      pause 'rtsafe exceeding maximum iterations at findmu'
      return
      end


      function fzz(zz,x)
      implicit real*8 (a-h,o-z)
      
      half = 0.5
      fzz = 2./3./zz**(1.5) - fd(half,x)
      
      return
      end



      function fd (xnu, alpha)
      implicit real*8 (a-h,o-z)
c w. fullerton, imsl houston, apr 1983 edition.
c
c fermi-dirac function of order xnu and argument alpha, where
c xnu is an integer multiple of 0.5 between -0.5 and 4.0
c inclusive.  this routine is accurate everywhere to approx-
c imately 10 digits.
c
      dimension fdscs(21,10), fdmcs(33,10), fdbcs(26,10),
     1  nscs(10), nmcs(10), nbcs(10), alfmax(10)
c
c  ixer, jxer, kxer flags for xerr returns
      data ixer,jxer,kxer/-1,2,2/
c
c series for fm5s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   5.75e-11
c                                         log weighted error  10.24
c                               significant figures required  10.32
c                                    decimal places required  10.90
c
      data fdscs(  1, 1) /   2.1468096512e0 /
      data fdscs(  2, 1) /   -.5099961735 5e0 /
      data fdscs(  3, 1) /    .1357694225 8e0 /
      data fdscs(  4, 1) /   -.0379066840 1e0 /
      data fdscs(  5, 1) /    .0108700231 1e0 /
      data fdscs(  6, 1) /   -.0031712624 6e0 /
      data fdscs(  7, 1) /    .0009364562 1e0 /
      data fdscs(  8, 1) /   -.0002790240 1e0 /
      data fdscs(  9, 1) /    .0000837159 4e0 /
      data fdscs( 10, 1) /   -.0000252565 8e0 /
      data fdscs( 11, 1) /    .0000076541 8e0 /
      data fdscs( 12, 1) /   -.0000023283 6e0 /
      data fdscs( 13, 1) /    .0000007105 2e0 /
      data fdscs( 14, 1) /   -.0000002174 1e0 /
      data fdscs( 15, 1) /    .0000000666 8e0 /
      data fdscs( 16, 1) /   -.0000000204 9e0 /
      data fdscs( 17, 1) /    .0000000063 1e0 /
      data fdscs( 18, 1) /   -.0000000019 4e0 /
      data fdscs( 19, 1) /    .0000000006 0e0 /
      data fdscs( 20, 1) /   -.0000000001 8e0 /
      data fdscs( 21, 1) /    .0000000000 57e0 /
c
c series for fm5m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   5.27e-11
c                                         log weighted error  10.28
c                               significant figures required  10.55
c                                    decimal places required  10.94
c
      data fdmcs(  1, 1) /   3.7374179482 1e0 /
      data fdmcs(  2, 1) /    .0702840192 74e0 /
      data fdmcs(  3, 1) /    .0049804498 35e0 /
      data fdmcs(  4, 1) /   -.0102460511 86e0 /
      data fdmcs(  5, 1) /    .0046134524 88e0 /
      data fdmcs(  6, 1) /   -.0015084055 45e0 /
      data fdmcs(  7, 1) /    .0004453628 58e0 /
      data fdmcs(  8, 1) /   -.0001323433 48e0 /
      data fdmcs(  9, 1) /    .0000410236 76e0 /
      data fdmcs( 10, 1) /   -.0000130973 69e0 /
      data fdmcs( 11, 1) /    .0000042172 38e0 /
      data fdmcs( 12, 1) /   -.0000013563 51e0 /
      data fdmcs( 13, 1) /    .0000004356 20e0 /
      data fdmcs( 14, 1) /   -.0000001400 58e0 /
      data fdmcs( 15, 1) /    .0000000451 47e0 /
      data fdmcs( 16, 1) /   -.0000000145 91e0 /
      data fdmcs( 17, 1) /    .0000000047 25e0 /
      data fdmcs( 18, 1) /   -.0000000015 33e0 /
      data fdmcs( 19, 1) /    .0000000004 97e0 /
      data fdmcs( 20, 1) /   -.0000000001 61e0 /
      data fdmcs( 21, 1) /    .0000000000 527e0 /
c
c series for fm5b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   2.77e-11
c                                         log weighted error  10.56
c                               significant figures required  10.85
c                                    decimal places required  11.26
c
      data fdbcs(  1, 1) /   3.9538986547 2e0 /
      data fdbcs(  2, 1) /   -.0314341746 06e0 /
      data fdbcs(  3, 1) /   -.0086537614 36e0 /
      data fdbcs(  4, 1) /   -.0000437949 13e0 /
      data fdbcs(  5, 1) /    .0003173634 01e0 /
      data fdbcs(  6, 1) /    .0000743852 22e0 /
      data fdbcs(  7, 1) /   -.0000283890 50e0 /
      data fdbcs(  8, 1) /   -.0000089391 79e0 /
      data fdbcs(  9, 1) /    .0000045170 94e0 /
      data fdbcs( 10, 1) /    .0000007855 81e0 /
      data fdbcs( 11, 1) /   -.0000009131 13e0 /
      data fdbcs( 12, 1) /    .0000000758 35e0 /
      data fdbcs( 13, 1) /    .0000001572 30e0 /
      data fdbcs( 14, 1) /   -.0000000671 69e0 /
      data fdbcs( 15, 1) /   -.0000000100 31e0 /
      data fdbcs( 16, 1) /    .0000000194 28e0 /
      data fdbcs( 17, 1) /   -.0000000059 47e0 /
      data fdbcs( 18, 1) /   -.0000000020 95e0 /
      data fdbcs( 19, 1) /    .0000000025 45e0 /
      data fdbcs( 20, 1) /   -.0000000007 47e0 /
      data fdbcs( 21, 1) /   -.0000000002 91e0 /
      data fdbcs( 22, 1) /    .0000000003 71e0 /
      data fdbcs( 23, 1) /   -.0000000001 31e0 /
      data fdbcs( 24, 1) /   -.0000000000 290e0 /
      data fdbcs( 25, 1) /    .0000000000 573e0 /
      data fdbcs( 26, 1) /   -.0000000000 277e0 /
c
c series for f00s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   3.18e-11
c                                         log weighted error  10.50
c                               significant figures required  10.36
c                                    decimal places required  11.15
c
      data fdscs(  1, 2) /   1.3659874054e0 /
      data fdscs(  2, 2) /   -.2438973101 9e0 /
      data fdscs(  3, 2) /    .0547680108 5e0 /
      data fdscs(  4, 2) /   -.0135159352 3e0 /
      data fdscs(  5, 2) /    .0035158670 3e0 /
      data fdscs(  6, 2) /   -.0009461111 9e0 /
      data fdscs(  7, 2) /    .0002607200 1e0 /
      data fdscs(  8, 2) /   -.0000731250 4e0 /
      data fdscs(  9, 2) /    .0000207911 1e0 /
      data fdscs( 10, 2) /   -.0000059759 5e0 /
      data fdscs( 11, 2) /    .0000017329 5e0 /
      data fdscs( 12, 2) /   -.0000005062 5e0 /
      data fdscs( 13, 2) /    .0000001488 2e0 /
      data fdscs( 14, 2) /   -.0000000439 8e0 /
      data fdscs( 15, 2) /    .0000000130 5e0 /
      data fdscs( 16, 2) /   -.0000000038 9e0 /
      data fdscs( 17, 2) /    .0000000011 6e0 /
      data fdscs( 18, 2) /   -.0000000003 4e0 /
      data fdscs( 19, 2) /    .0000000001 0e0 /
      data fdscs( 20, 2) /   -.0000000000 31e0 /
c
c series for f00m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   6.63e-11
c                                         log weighted error  10.18
c                               significant figures required  10.22
c                                    decimal places required  10.85
c
      data fdmcs(  1, 2) /   2.1728598421e0 /
      data fdmcs(  2, 2) /   -.1263072405 5e0 /
      data fdmcs(  3, 2) /    .0627040952 7e0 /
      data fdmcs(  4, 2) /   -.0248015923 3e0 /
      data fdmcs(  5, 2) /    .0086910218 0e0 /
      data fdmcs(  6, 2) /   -.0028969250 7e0 /
      data fdmcs(  7, 2) /    .0009558659 0e0 /
      data fdmcs(  8, 2) /   -.0003167553 1e0 /
      data fdmcs(  9, 2) /    .0001054756 2e0 /
      data fdmcs( 10, 2) /   -.0000351869 0e0 /
      data fdmcs( 11, 2) /    .0000117376 0e0 /
      data fdmcs( 12, 2) /   -.0000039134 7e0 /
      data fdmcs( 13, 2) /    .0000013044 3e0 /
      data fdmcs( 14, 2) /   -.0000004347 7e0 /
      data fdmcs( 15, 2) /    .0000001449 1e0 /
      data fdmcs( 16, 2) /   -.0000000483 0e0 /
      data fdmcs( 17, 2) /    .0000000161 0e0 /
      data fdmcs( 18, 2) /   -.0000000053 6e0 /
      data fdmcs( 19, 2) /    .0000000017 8e0 /
      data fdmcs( 20, 2) /   -.0000000005 9e0 /
      data fdmcs( 21, 2) /    .0000000001 9e0 /
      data fdmcs( 22, 2) /   -.0000000000 66e0 /
c
c series for f00b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   2.19e-11
c                                         log weighted error  10.66
c                               significant figures required  10.66
c                                    decimal places required  11.33
c
      data fdbcs(  1, 2) /   2.0021325464 9e0 /
      data fdbcs(  2, 2) /    .0018258034 26e0 /
      data fdbcs(  3, 2) /    .0011225105 54e0 /
      data fdbcs(  4, 2) /    .0004555439 28e0 /
      data fdbcs(  5, 2) /    .0000880185 91e0 /
      data fdbcs(  6, 2) /   -.0000137060 73e0 /
      data fdbcs(  7, 2) /   -.0000090015 63e0 /
      data fdbcs(  8, 2) /    .0000013512 84e0 /
      data fdbcs(  9, 2) /    .0000010449 37e0 /
      data fdbcs( 10, 2) /   -.0000003164 95e0 /
      data fdbcs( 11, 2) /   -.0000001071 28e0 /
      data fdbcs( 12, 2) /    .0000000783 47e0 /
      data fdbcs( 13, 2) /   -.0000000017 22e0 /
      data fdbcs( 14, 2) /   -.0000000148 37e0 /
      data fdbcs( 15, 2) /    .0000000055 82e0 /
      data fdbcs( 16, 2) /    .0000000010 85e0 /
      data fdbcs( 17, 2) /   -.0000000017 29e0 /
      data fdbcs( 18, 2) /    .0000000005 16e0 /
      data fdbcs( 19, 2) /    .0000000001 80e0 /
      data fdbcs( 20, 2) /   -.0000000002 22e0 /
      data fdbcs( 21, 2) /    .0000000000 691e0 /
      data fdbcs( 22, 2) /    .0000000000 219e0 /
c
c series for f05s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   2.87e-11
c                                         log weighted error  10.54
c                               significant figures required  10.38
c                                    decimal places required  11.18
c
      data fdscs(  1, 3) /   1.3284029124e0 /
      data fdscs(  2, 3) /   -.1783849741 0e0 /
      data fdscs(  3, 3) /    .0338871501 3e0 /
      data fdscs(  4, 3) /   -.0074116985 5e0 /
      data fdscs(  5, 3) /    .0017528135 4e0 /
      data fdscs(  6, 3) /   -.0004358562 3e0 /
      data fdscs(  7, 3) /    .0001122556 1e0 /
      data fdscs(  8, 3) /   -.0000296748 6e0 /
      data fdscs(  9, 3) /    .0000080041 3e0 /
      data fdscs( 10, 3) /   -.0000021938 6e0 /
      data fdscs( 11, 3) /    .0000006092 5e0 /
      data fdscs( 12, 3) /   -.0000001710 4e0 /
      data fdscs( 13, 3) /    .0000000484 6e0 /
      data fdscs( 14, 3) /   -.0000000138 4e0 /
      data fdscs( 15, 3) /    .0000000039 8e0 /
      data fdscs( 16, 3) /   -.0000000011 5e0 /
      data fdscs( 17, 3) /    .0000000003 3e0 /
      data fdscs( 18, 3) /   -.0000000000 97e0 /
      data fdscs( 19, 3) /    .0000000000 28e0 /
c
c series for f05m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   2.90e-11
c                                         log weighted error  10.54
c                               significant figures required  10.52
c                                    decimal places required  11.23
c
      data fdmcs(  1, 3) /   1.8324909204 1e0 /
      data fdmcs(  2, 3) /   -.2831869649 90e0 /
      data fdmcs(  3, 3) /    .1228778370 71e0 /
      data fdmcs(  4, 3) /   -.0473020885 11e0 /
      data fdmcs(  5, 3) /    .0172527816 16e0 /
      data fdmcs(  6, 3) /   -.0061561464 96e0 /
      data fdmcs(  7, 3) /    .0021769305 63e0 /
      data fdmcs(  8, 3) /   -.0007657618 44e0 /
      data fdmcs(  9, 3) /    .0002681487 29e0 /
      data fdmcs( 10, 3) /   -.0000935040 12e0 /
      data fdmcs( 11, 3) /    .0000324836 07e0 /
      data fdmcs( 12, 3) /   -.0000112489 79e0 /
      data fdmcs( 13, 3) /    .0000038849 42e0 /
      data fdmcs( 14, 3) /   -.0000013385 73e0 /
      data fdmcs( 15, 3) /    .0000004602 70e0 /
      data fdmcs( 16, 3) /   -.0000001579 79e0 /
      data fdmcs( 17, 3) /    .0000000541 36e0 /
      data fdmcs( 18, 3) /   -.0000000185 24e0 /
      data fdmcs( 19, 3) /    .0000000063 30e0 /
      data fdmcs( 20, 3) /   -.0000000021 60e0 /
      data fdmcs( 21, 3) /    .0000000007 36e0 /
      data fdmcs( 22, 3) /   -.0000000002 50e0 /
      data fdmcs( 23, 3) /    .0000000000 854e0 /
      data fdmcs( 24, 3) /   -.0000000000 290e0 /
c
c series for f05b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   2.38e-11
c                                         log weighted error  10.62
c                               significant figures required  10.46
c                                    decimal places required  11.27
c
      data fdbcs(  1, 3) /   1.3737552319e0 /
      data fdbcs(  2, 3) /    .0271865732 7e0 /
      data fdbcs(  3, 3) /    .0071381409 2e0 /
      data fdbcs(  4, 3) /    .0001635361 7e0 /
      data fdbcs(  5, 3) /   -.0000117682 0e0 /
      data fdbcs(  6, 3) /   -.0000141643 9e0 /
      data fdbcs(  7, 3) /   -.0000002415 6e0 /
      data fdbcs(  8, 3) /    .0000012719 9e0 /
      data fdbcs(  9, 3) /   -.0000000281 9e0 /
      data fdbcs( 10, 3) /   -.0000001602 3e0 /
      data fdbcs( 11, 3) /    .0000000306 4e0 /
      data fdbcs( 12, 3) /    .0000000187 1e0 /
      data fdbcs( 13, 3) /   -.0000000101 5e0 /
      data fdbcs( 14, 3) /   -.0000000003 8e0 /
      data fdbcs( 15, 3) /    .0000000021 2e0 /
      data fdbcs( 16, 3) /   -.0000000007 1e0 /
      data fdbcs( 17, 3) /   -.0000000001 7e0 /
      data fdbcs( 18, 3) /    .0000000002 3e0 /
      data fdbcs( 19, 3) /   -.0000000000 70e0 /
      data fdbcs( 20, 3) /   -.0000000000 23e0 /
c
c series for f10s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   3.46e-11
c                                         log weighted error  10.46
c                               significant figures required  10.38
c                                    decimal places required  11.09
c
      data fdscs(  1, 4) /   1.6098847156e0 /
      data fdscs(  2, 4) /   -.1624060718 4e0 /
      data fdscs(  3, 4) /    .0261468226 3e0 /
      data fdscs(  4, 4) /   -.0050782609 1e0 /
      data fdscs(  5, 4) /    .0010937471 2e0 /
      data fdscs(  6, 4) /   -.0002516893 5e0 /
      data fdscs(  7, 4) /    .0000606609 8e0 /
      data fdscs(  8, 4) /   -.0000151303 0e0 /
      data fdscs(  9, 4) /    .0000038751 7e0 /
      data fdscs( 10, 4) /   -.0000010136 9e0 /
      data fdscs( 11, 4) /    .0000002697 7e0 /
      data fdscs( 12, 4) /   -.0000000728 3e0 /
      data fdscs( 13, 4) /    .0000000199 0e0 /
      data fdscs( 14, 4) /   -.0000000054 9e0 /
      data fdscs( 15, 4) /    .0000000015 3e0 /
      data fdscs( 16, 4) /   -.0000000004 3e0 /
      data fdscs( 17, 4) /    .0000000001 2e0 /
      data fdscs( 18, 4) /   -.0000000000 34e0 /
c
c series for f10m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   3.92e-11
c                                         log weighted error  10.41
c                               significant figures required  10.44
c                                    decimal places required  11.11
c
      data fdmcs(  1, 4) /   1.8823916245e0 /
      data fdmcs(  2, 4) /   -.4963882506 6e0 /
      data fdmcs(  3, 4) /    .2216981727 0e0 /
      data fdmcs(  4, 4) /   -.0903287142 0e0 /
      data fdmcs(  5, 4) /    .0352537623 9e0 /
      data fdmcs(  6, 4) /   -.0134370871 3e0 /
      data fdmcs(  7, 4) /    .0050410084 6e0 /
      data fdmcs(  8, 4) /   -.0018681458 3e0 /
      data fdmcs(  9, 4) /    .0006853986 5e0 /
      data fdmcs( 10, 4) /   -.0002493647 5e0 /
      data fdmcs( 11, 4) /    .0000900867 6e0 /
      data fdmcs( 12, 4) /   -.0000323503 7e0 /
      data fdmcs( 13, 4) /    .0000115572 5e0 /
      data fdmcs( 14, 4) /   -.0000041103 4e0 /
      data fdmcs( 15, 4) /    .0000014560 9e0 /
      data fdmcs( 16, 4) /   -.0000005140 2e0 /
      data fdmcs( 17, 4) /    .0000001808 9e0 /
      data fdmcs( 18, 4) /   -.0000000634 8e0 /
      data fdmcs( 19, 4) /    .0000000222 2e0 /
      data fdmcs( 20, 4) /   -.0000000077 6e0 /
      data fdmcs( 21, 4) /    .0000000027 0e0 /
      data fdmcs( 22, 4) /   -.0000000009 4e0 /
      data fdmcs( 23, 4) /    .0000000003 2e0 /
      data fdmcs( 24, 4) /   -.0000000001 1e0 /
      data fdmcs( 25, 4) /    .0000000000 39e0 /
c
c series for f10b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   3.24e-11
c                                         log weighted error  10.49
c                               significant figures required  10.22
c                                    decimal places required  11.12
c
      data fdbcs(  1, 4) /   1.0766097168e0 /
      data fdbcs(  2, 4) /    .0509708968 7e0 /
      data fdbcs(  3, 4) /    .0125669010 0e0 /
      data fdbcs(  4, 4) /   -.0001333933 2e0 /
      data fdbcs(  5, 4) /   -.0000390220 0e0 /
      data fdbcs(  6, 4) /   -.0000033829 5e0 /
      data fdbcs(  7, 4) /    .0000018566 7e0 /
      data fdbcs(  8, 4) /    .0000003251 2e0 /
      data fdbcs(  9, 4) /   -.0000001931 9e0 /
      data fdbcs( 10, 4) /   -.0000000183 8e0 /
      data fdbcs( 11, 4) /    .0000000281 7e0 /
      data fdbcs( 12, 4) /   -.0000000030 6e0 /
      data fdbcs( 13, 4) /   -.0000000037 4e0 /
      data fdbcs( 14, 4) /    .0000000016 2e0 /
      data fdbcs( 15, 4) /    .0000000001 6e0 /
      data fdbcs( 16, 4) /   -.0000000003 7e0 /
      data fdbcs( 17, 4) /    .0000000001 1e0 /
      data fdbcs( 18, 4) /    .0000000000 32e0 /
c
c series for f15s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   5.17e-11
c                                         log weighted error  10.29
c                               significant figures required  10.34
c                                    decimal places required  10.90
c
      data fdscs(  1, 5) /   2.2601818297e0 /
      data fdscs(  2, 5) /   -.1708722765 0e0 /
      data fdscs(  3, 5) /    .0233363666 4e0 /
      data fdscs(  4, 5) /   -.0040304134 4e0 /
      data fdscs(  5, 5) /    .0007916285 5e0 /
      data fdscs(  6, 5) /   -.0001687845 2e0 /
      data fdscs(  7, 5) /    .0000381078 6e0 /
      data fdscs(  8, 5) /   -.0000089765 6e0 /
      data fdscs(  9, 5) /    .0000021848 5e0 /
      data fdscs( 10, 5) /   -.0000005458 3e0 /
      data fdscs( 11, 5) /    .0000001393 0e0 /
      data fdscs( 12, 5) /   -.0000000361 8e0 /
      data fdscs( 13, 5) /    .0000000095 4e0 /
      data fdscs( 14, 5) /   -.0000000025 4e0 /
      data fdscs( 15, 5) /    .0000000006 8e0 /
      data fdscs( 16, 5) /   -.0000000001 8e0 /
      data fdscs( 17, 5) /    .0000000000 51e0 /
c
c series for f15m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   5.45e-11
c                                         log weighted error  10.26
c                               significant figures required  10.43
c                                    decimal places required  10.97
c
      data fdmcs(  1, 5) /   2.2310484506e0 /
      data fdmcs(  2, 5) /   -.8435086752 5e0 /
      data fdmcs(  3, 5) /    .4030118348 2e0 /
      data fdmcs(  4, 5) /   -.1765850337 0e0 /
      data fdmcs(  5, 5) /    .0738231089 9e0 /
      data fdmcs(  6, 5) /   -.0299237390 3e0 /
      data fdmcs(  7, 5) /    .0118550148 7e0 /
      data fdmcs(  8, 5) /   -.0046128514 7e0 /
      data fdmcs(  9, 5) /    .0017688631 8e0 /
      data fdmcs( 10, 5) /   -.0006701511 6e0 /
      data fdmcs( 11, 5) /    .0002513292 1e0 /
      data fdmcs( 12, 5) /   -.0000934452 7e0 /
      data fdmcs( 13, 5) /    .0000344851 9e0 /
      data fdmcs( 14, 5) /   -.0000126440 2e0 /
      data fdmcs( 15, 5) /    .0000046095 2e0 /
      data fdmcs( 16, 5) /   -.0000016719 6e0 /
      data fdmcs( 17, 5) /    .0000006037 1e0 /
      data fdmcs( 18, 5) /   -.0000002171 0e0 /
      data fdmcs( 19, 5) /    .0000000777 9e0 /
      data fdmcs( 20, 5) /   -.0000000277 8e0 /
      data fdmcs( 21, 5) /    .0000000098 9e0 /
      data fdmcs( 22, 5) /   -.0000000035 1e0 /
      data fdmcs( 23, 5) /    .0000000012 4e0 /
      data fdmcs( 24, 5) /   -.0000000004 3e0 /
      data fdmcs( 25, 5) /    .0000000001 5e0 /
      data fdmcs( 26, 5) /   -.0000000000 54e0 /
c
c series for f15b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   2.21e-11
c                                         log weighted error  10.66
c                               significant figures required  10.32
c                                    decimal places required  11.28
c
      data fdbcs(  1, 5) /    .9138458031 3e0 /
      data fdbcs(  2, 5) /    .0756461485 3e0 /
      data fdbcs(  3, 5) /    .0185325720 6e0 /
      data fdbcs(  4, 5) /   -.0002173856 5e0 /
      data fdbcs(  5, 5) /   -.0000237328 8e0 /
      data fdbcs(  6, 5) /    .0000042673 3e0 /
      data fdbcs(  7, 5) /    .0000012018 7e0 /
      data fdbcs(  8, 5) /   -.0000002038 1e0 /
      data fdbcs(  9, 5) /   -.0000000983 9e0 /
      data fdbcs( 10, 5) /    .0000000298 3e0 /
      data fdbcs( 11, 5) /    .0000000073 9e0 /
      data fdbcs( 12, 5) /   -.0000000054 3e0 /
      data fdbcs( 13, 5) /    .0000000001 9e0 /
      data fdbcs( 14, 5) /    .0000000008 2e0 /
      data fdbcs( 15, 5) /   -.0000000002 9e0 /
      data fdbcs( 16, 5) /   -.0000000000 48e0 /
      data fdbcs( 17, 5) /    .0000000000 77e0 /
      data fdbcs( 18, 5) /   -.0000000000 22e0 /
c
c series for f20s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   2.47e-11
c                                         log weighted error  10.61
c                               significant figures required  10.86
c                                    decimal places required  11.22
c
      data fdscs(  1, 6) /   3.5445815749 6e0 /
      data fdscs(  2, 6) /   -.2001497509 36e0 /
      data fdscs(  3, 6) /    .0231937129 11e0 /
      data fdscs(  4, 6) /   -.0035654858 18e0 /
      data fdscs(  5, 6) /    .0006393090 63e0 /
      data fdscs(  6, 6) /   -.0001264180 87e0 /
      data fdscs(  7, 6) /    .0000267615 70e0 /
      data fdscs(  8, 6) /   -.0000059580 71e0 /
      data fdscs(  9, 6) /    .0000013790 87e0 /
      data fdscs( 10, 6) /   -.0000003292 55e0 /
      data fdscs( 11, 6) /    .0000000806 22e0 /
      data fdscs( 12, 6) /   -.0000000201 61e0 /
      data fdscs( 13, 6) /    .0000000051 32e0 /
      data fdscs( 14, 6) /   -.0000000013 26e0 /
      data fdscs( 15, 6) /    .0000000003 47e0 /
      data fdscs( 16, 6) /   -.0000000000 921e0 /
      data fdscs( 17, 6) /    .0000000000 246e0 /
c
c series for f20m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   2.78e-11
c                                         log weighted error  10.56
c                               significant figures required  10.91
c                                    decimal places required  11.28
c
      data fdmcs(  1, 6) /   2.9777839001 0e0 /
      data fdmcs(  2, 6) /  -1.4577413536 2e0 /
      data fdmcs(  3, 6) /    .7528269512 35e0 /
      data fdmcs(  4, 6) /   -.3549647428 05e0 /
      data fdmcs(  5, 6) /    .1584014924 88e0 /
      data fdmcs(  6, 6) /   -.0680073485 74e0 /
      data fdmcs(  7, 6) /    .0283569566 67e0 /
      data fdmcs(  8, 6) /   -.0115545568 43e0 /
      data fdmcs(  9, 6) /    .0046209871 95e0 /
      data fdmcs( 10, 6) /   -.0018197198 20e0 /
      data fdmcs( 11, 6) /    .0007073370 31e0 /
      data fdmcs( 12, 6) /   -.0002719114 95e0 /
      data fdmcs( 13, 6) /    .0001035295 30e0 /
      data fdmcs( 14, 6) /   -.0000390900 34e0 /
      data fdmcs( 15, 6) /    .0000146509 87e0 /
      data fdmcs( 16, 6) /   -.0000054554 02e0 /
      data fdmcs( 17, 6) /    .0000020195 19e0 /
      data fdmcs( 18, 6) /   -.0000007436 80e0 /
      data fdmcs( 19, 6) /    .0000002725 59e0 /
      data fdmcs( 20, 6) /   -.0000000994 63e0 /
      data fdmcs( 21, 6) /    .0000000361 53e0 /
      data fdmcs( 22, 6) /   -.0000000130 94e0 /
      data fdmcs( 23, 6) /    .0000000047 26e0 /
      data fdmcs( 24, 6) /   -.0000000017 01e0 /
      data fdmcs( 25, 6) /    .0000000006 10e0 /
      data fdmcs( 26, 6) /   -.0000000002 18e0 /
      data fdmcs( 27, 6) /    .0000000000 780e0 /
      data fdmcs( 28, 6) /   -.0000000000 277e0 /
c
c series for f20b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   1.41e-11
c                                         log weighted error  10.85
c                               significant figures required  10.48
c                                    decimal places required  11.47
c
      data fdbcs(  1, 6) /    .8211121274 8e0 /
      data fdbcs(  2, 6) /    .1030146856 2e0 /
      data fdbcs(  3, 6) /    .0258442758 9e0 /
      data fdbcs(  4, 6) /    .0000739478 7e0 /
      data fdbcs(  5, 6) /    .0000269632 3e0 /
      data fdbcs(  6, 6) /    .0000055393 5e0 /
      data fdbcs(  7, 6) /   -.0000000666 0e0 /
      data fdbcs(  8, 6) /   -.0000002863 2e0 /
      data fdbcs(  9, 6) /    .0000000098 7e0 /
      data fdbcs( 10, 6) /    .0000000250 2e0 /
      data fdbcs( 11, 6) /   -.0000000043 8e0 /
      data fdbcs( 12, 6) /   -.0000000022 7e0 /
      data fdbcs( 13, 6) /    .0000000011 2e0 /
      data fdbcs( 14, 6) /    .0000000000 40e0 /
      data fdbcs( 15, 6) /   -.0000000001 9e0 /
      data fdbcs( 16, 6) /    .0000000000 60e0 /
      data fdbcs( 17, 6) /    .0000000000 14e0 /
c
c series for f25s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   4.95e-11
c                                         log weighted error  10.31
c                               significant figures required  10.79
c                                    decimal places required  10.91
c
      data fdscs(  1, 7) /   6.0776352655 5e0 /
      data fdscs(  2, 7) /   -.2553089335 90e0 /
      data fdscs(  3, 7) /    .0250962593 28e0 /
      data fdscs(  4, 7) /   -.0034359138 78e0 /
      data fdscs(  5, 7) /    .0005628501 70e0 /
      data fdscs(  6, 7) /   -.0001033045 44e0 /
      data fdscs(  7, 7) /    .0000205192 58e0 /
      data fdscs(  8, 7) /   -.0000043206 22e0 /
      data fdscs(  9, 7) /    .0000009516 32e0 /
      data fdscs( 10, 7) /   -.0000002172 43e0 /
      data fdscs( 11, 7) /    .0000000510 64e0 /
      data fdscs( 12, 7) /   -.0000000122 98e0 /
      data fdscs( 13, 7) /    .0000000030 23e0 /
      data fdscs( 14, 7) /   -.0000000007 56e0 /
      data fdscs( 15, 7) /    .0000000001 92e0 /
      data fdscs( 16, 7) /   -.0000000000 495e0 /
c
c series for f25m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   4.15e-11
c                                         log weighted error  10.38
c                               significant figures required  10.96
c                                    decimal places required  11.11
c
      data fdmcs(  1, 7) /   4.4112295532 9e0 /
      data fdmcs(  2, 7) /  -2.6064119236 9e0 /
      data fdmcs(  3, 7) /   1.4541147675 7e0 /
      data fdmcs(  4, 7) /   -.7348922045 75e0 /
      data fdmcs(  5, 7) /    .3485252773 58e0 /
      data fdmcs(  6, 7) /   -.1579050700 24e0 /
      data fdmcs(  7, 7) /    .0690927288 72e0 /
      data fdmcs(  8, 7) /   -.0294110574 69e0 /
      data fdmcs(  9, 7) /    .0122427694 58e0 /
      data fdmcs( 10, 7) /   -.0050026362 88e0 /
      data fdmcs( 11, 7) /    .0020124730 44e0 /
      data fdmcs( 12, 7) /   -.0007988327 90e0 /
      data fdmcs( 13, 7) /    .0003134423 09e0 /
      data fdmcs( 14, 7) /   -.0001217494 74e0 /
      data fdmcs( 15, 7) /    .0000468708 54e0 /
      data fdmcs( 16, 7) /   -.0000179017 70e0 /
      data fdmcs( 17, 7) /    .0000067890 45e0 /
      data fdmcs( 18, 7) /   -.0000025582 83e0 /
      data fdmcs( 19, 7) /    .0000009584 71e0 /
      data fdmcs( 20, 7) /   -.0000003572 13e0 /
      data fdmcs( 21, 7) /    .0000001324 92e0 /
      data fdmcs( 22, 7) /   -.0000000489 26e0 /
      data fdmcs( 23, 7) /    .0000000179 94e0 /
      data fdmcs( 24, 7) /   -.0000000065 93e0 /
      data fdmcs( 25, 7) /    .0000000024 07e0 /
      data fdmcs( 26, 7) /   -.0000000008 76e0 /
      data fdmcs( 27, 7) /    .0000000003 17e0 /
      data fdmcs( 28, 7) /   -.0000000001 15e0 /
      data fdmcs( 29, 7) /    .0000000000 415e0 /
c
c series for f25b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   4.96e-11
c                                         log weighted error  10.30
c                               significant figures required   9.92
c                                    decimal places required  10.91
c
      data fdbcs(  1, 7) /    .7721541990 3e0 /
      data fdbcs(  2, 7) /    .1348979022 5e0 /
      data fdbcs(  3, 7) /    .0353565117 1e0 /
      data fdbcs(  4, 7) /    .0009476728 1e0 /
      data fdbcs(  5, 7) /    .0001277291 1e0 /
      data fdbcs(  6, 7) /    .0000007218 8e0 /
      data fdbcs(  7, 7) /   -.0000009206 8e0 /
      data fdbcs(  8, 7) /   -.0000001154 7e0 /
      data fdbcs(  9, 7) /    .0000000580 6e0 /
      data fdbcs( 10, 7) /    .0000000047 4e0 /
      data fdbcs( 11, 7) /   -.0000000060 8e0 /
      data fdbcs( 12, 7) /    .0000000005 1e0 /
      data fdbcs( 13, 7) /    .0000000006 5e0 /
      data fdbcs( 14, 7) /   -.0000000002 4e0 /
      data fdbcs( 15, 7) /   -.0000000000 28e0 /
      data fdbcs( 16, 7) /    .0000000000 49e0 /
c
c series for f30s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   2.87e-11
c                                         log weighted error  10.54
c                               significant figures required  11.29
c                                    decimal places required  11.14
c
      data fdscs(  1, 8) /  11.2341939777e0 /
      data fdscs(  2, 8) /   -.3495789894 02e0 /
      data fdscs(  3, 8) /    .0291275872 60e0 /
      data fdscs(  4, 8) /   -.0035525827 96e0 /
      data fdscs(  5, 8) /    .0005319821 79e0 /
      data fdscs(  6, 8) /   -.0000906823 61e0 /
      data fdscs(  7, 8) /    .0000169110 38e0 /
      data fdscs(  8, 8) /   -.0000033697 23e0 /
      data fdscs(  9, 8) /    .0000007066 16e0 /
      data fdscs( 10, 8) /   -.0000001543 15e0 /
      data fdscs( 11, 8) /    .0000000348 35e0 /
      data fdscs( 12, 8) /   -.0000000080 83e0 /
      data fdscs( 13, 8) /    .0000000019 20e0 /
      data fdscs( 14, 8) /   -.0000000004 65e0 /
      data fdscs( 15, 8) /    .0000000001 14e0 /
      data fdscs( 16, 8) /   -.0000000000 287e0 /
c
c series for f30m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   2.33e-11
c                                         log weighted error  10.63
c                               significant figures required  11.47
c                                    decimal places required  11.38
c
      data fdmcs(  1, 8) /   7.1719763668 6e0 /
      data fdmcs(  2, 8) /  -4.8571185185 4e0 /
      data fdmcs(  3, 8) /   2.9113197795 3e0 /
      data fdmcs(  4, 8) /  -1.5684642300 9e0 /
      data fdmcs(  5, 8) /    .7869998941 23e0 /
      data fdmcs(  6, 8) /   -.3749544936 90e0 /
      data fdmcs(  7, 8) /    .1716880079 87e0 /
      data fdmcs(  8, 8) /   -.0761760305 76e0 /
      data fdmcs(  9, 8) /    .0329421355 00e0 /
      data fdmcs( 10, 8) /   -.0139450242 81e0 /
      data fdmcs( 11, 8) /    .0057976754 88e0 /
      data fdmcs( 12, 8) /   -.0023734227 03e0 /
      data fdmcs( 13, 8) /    .0009586830 99e0 /
      data fdmcs( 14, 8) /   -.0003827164 22e0 /
      data fdmcs( 15, 8) /    .0001512084 34e0 /
      data fdmcs( 16, 8) /   -.0000591925 75e0 /
      data fdmcs( 17, 8) /    .0000229809 46e0 /
      data fdmcs( 18, 8) /   -.0000088559 00e0 /
      data fdmcs( 19, 8) /    .0000033897 35e0 /
      data fdmcs( 20, 8) /   -.0000012895 26e0 /
      data fdmcs( 21, 8) /    .0000004878 14e0 /
      data fdmcs( 22, 8) /   -.0000001835 85e0 /
      data fdmcs( 23, 8) /    .0000000687 64e0 /
      data fdmcs( 24, 8) /   -.0000000256 43e0 /
      data fdmcs( 25, 8) /    .0000000095 24e0 /
      data fdmcs( 26, 8) /   -.0000000035 23e0 /
      data fdmcs( 27, 8) /    .0000000012 99e0 /
      data fdmcs( 28, 8) /   -.0000000004 77e0 /
      data fdmcs( 29, 8) /    .0000000001 74e0 /
      data fdmcs( 30, 8) /   -.0000000000 638e0 /
      data fdmcs( 31, 8) /    .0000000000 232e0 /
c
c series for f30b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   5.48e-11
c                                         log weighted error  10.26
c                               significant figures required   9.88
c                                    decimal places required  10.85
c
      data fdbcs(  1, 8) /    .7554309639 2e0 /
      data fdbcs(  2, 8) /    .1734863060 3e0 /
      data fdbcs(  3, 8) /    .0481579481 5e0 /
      data fdbcs(  4, 8) /    .0027149874 6e0 /
      data fdbcs(  5, 8) /    .0003217541 9e0 /
      data fdbcs(  6, 8) /   -.0000071412 9e0 /
      data fdbcs(  7, 8) /   -.0000009676 6e0 /
      data fdbcs(  8, 8) /    .0000001160 1e0 /
      data fdbcs(  9, 8) /    .0000000450 4e0 /
      data fdbcs( 10, 8) /   -.0000000103 6e0 /
      data fdbcs( 11, 8) /   -.0000000025 9e0 /
      data fdbcs( 12, 8) /    .0000000014 6e0 /
      data fdbcs( 13, 8) /   -.0000000000 04e0 /
      data fdbcs( 14, 8) /   -.0000000001 8e0 /
      data fdbcs( 15, 8) /    .0000000000 54e0 /
c
c series for f35s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   7.30e-11
c                                         log weighted error  10.14
c                               significant figures required  11.18
c                                    decimal places required  10.72
c
      data fdscs(  1, 9) /  22.1653046970e0 /
      data fdscs(  2, 9) /   -.5086526539 69e0 /
      data fdscs(  3, 9) /    .0358871327 23e0 /
      data fdscs(  4, 9) /   -.0038993959 73e0 /
      data fdscs(  5, 9) /    .0005339699 07e0 /
      data fdscs(  6, 9) /   -.0000845770 08e0 /
      data fdscs(  7, 9) /    .0000148157 47e0 /
      data fdscs(  8, 9) /   -.0000027951 08e0 /
      data fdscs(  9, 9) /    .0000005582 82e0 /
      data fdscs( 10, 9) /   -.0000001166 84e0 /
      data fdscs( 11, 9) /    .0000000253 06e0 /
      data fdscs( 12, 9) /   -.0000000056 60e0 /
      data fdscs( 13, 9) /    .0000000012 99e0 /
      data fdscs( 14, 9) /   -.0000000003 05e0 /
      data fdscs( 15, 9) /    .0000000000 730e0 /
c
c series for f35m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   3.73e-11
c                                         log weighted error  10.43
c                               significant figures required  11.56
c                                    decimal places required  11.18
c
      data fdmcs(  1, 9) /  12.6701567503 6e0 /
      data fdmcs(  2, 9) /  -9.4636920164 00e0 /
      data fdmcs(  3, 9) /   6.0480654287 69e0 /
      data fdmcs(  4, 9) /  -3.4530339209 92e0 /
      data fdmcs(  5, 9) /   1.8250457226 69e0 /
      data fdmcs(  6, 9) /   -.9112870648 186e0 /
      data fdmcs(  7, 9) /    .4354953493 280e0 /
      data fdmcs(  8, 9) /   -.2009635658 884e0 /
      data fdmcs(  9, 9) /    .0901210173 526e0 /
      data fdmcs( 10, 9) /   -.0394612435 160e0 /
      data fdmcs( 11, 9) /    .0169328410 948e0 /
      data fdmcs( 12, 9) /   -.0071407017 340e0 /
      data fdmcs( 13, 9) /    .0029661522 591e0 /
      data fdmcs( 14, 9) /   -.0012158829 523e0 /
      data fdmcs( 15, 9) /    .0004926051 670e0 /
      data fdmcs( 16, 9) /   -.0001975006 123e0 /
      data fdmcs( 17, 9) /    .0000784453 353e0 /
      data fdmcs( 18, 9) /   -.0000308953 181e0 /
      data fdmcs( 19, 9) /    .0000120749 876e0 /
      data fdmcs( 20, 9) /   -.0000046864 594e0 /
      data fdmcs( 21, 9) /    .0000018072 775e0 /
      data fdmcs( 22, 9) /   -.0000006928 714e0 /
      data fdmcs( 23, 9) /    .0000002641 967e0 /
      data fdmcs( 24, 9) /   -.0000001002 365e0 /
      data fdmcs( 25, 9) /    .0000000378 535e0 /
      data fdmcs( 26, 9) /   -.0000000142 333e0 /
      data fdmcs( 27, 9) /    .0000000053 303e0 /
      data fdmcs( 28, 9) /   -.0000000019 887e0 /
      data fdmcs( 29, 9) /    .0000000007 393e0 /
      data fdmcs( 30, 9) /   -.0000000002 739e0 /
      data fdmcs( 31, 9) /    .0000000001 011e0 /
      data fdmcs( 32, 9) /   -.0000000000 372e0 /
c
c series for f35b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   5.48e-11
c                                         log weighted error  10.26
c                               significant figures required   9.91
c                                    decimal places required  10.85
c
      data fdbcs(  1, 9) /    .7665702804 4e0 /
      data fdbcs(  2, 9) /    .2216671315 0e0 /
      data fdbcs(  3, 9) /    .0657600143 6e0 /
      data fdbcs(  4, 9) /    .0058625466 3e0 /
      data fdbcs(  5, 9) /    .0006969431 7e0 /
      data fdbcs(  6, 9) /   -.0000101109 1e0 /
      data fdbcs(  7, 9) /   -.0000000660 4e0 /
      data fdbcs(  8, 9) /    .0000002536 7e0 /
      data fdbcs(  9, 9) /   -.0000000021 4e0 /
      data fdbcs( 10, 9) /   -.0000000132 6e0 /
      data fdbcs( 11, 9) /    .0000000015 0e0 /
      data fdbcs( 12, 9) /    .0000000009 5e0 /
      data fdbcs( 13, 9) /   -.0000000003 4e0 /
      data fdbcs( 14, 9) /   -.0000000000 30e0 /
      data fdbcs( 15, 9) /    .0000000000 54e0 /
c
c series for f40s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   4.91e-11
c                                         log weighted error  10.31
c                               significant figures required  11.67
c                                    decimal places required  10.90
c
      data fdscs(  1,10) /  46.3350918683 9e0 /
      data fdscs(  2,10) /   -.7807026145 628e0 /
      data fdscs(  3,10) /    .0465789224 743e0 /
      data fdscs(  4,10) /   -.0045080435 976e0 /
      data fdscs(  5,10) /    .0005646381 627e0 /
      data fdscs(  6,10) /   -.0000831331 627e0 /
      data fdscs(  7,10) /    .0000136850 757e0 /
      data fdscs(  8,10) /   -.0000024454 133e0 /
      data fdscs(  9,10) /    .0000004654 205e0 /
      data fdscs( 10,10) /   -.0000000931 320e0 /
      data fdscs( 11,10) /    .0000000194 128e0 /
      data fdscs( 12,10) /   -.0000000041 863e0 /
      data fdscs( 13,10) /    .0000000009 290e0 /
      data fdscs( 14,10) /   -.0000000002 113e0 /
      data fdscs( 15,10) /    .0000000000 491e0 /
c
c series for f40m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   6.13e-11
c                                         log weighted error  10.21
c                               significant figures required  11.66
c                                    decimal places required  10.97
c
      data fdmcs(  1,10) /  24.0980879457 2e0 /
      data fdmcs(  2,10) / -19.2973238247 9e0 /
      data fdmcs(  3,10) /  13.0411335433 2e0 /
      data fdmcs(  4,10) /  -7.8442177530 69e0 /
      data fdmcs(  5,10) /   4.3484777309 21e0 /
      data fdmcs(  6,10) /  -2.2682065310 65e0 /
      data fdmcs(  7,10) /   1.1283901506 72e0 /
      data fdmcs(  8,10) /   -.5404257776 187e0 /
      data fdmcs(  9,10) /    .2508755478 873e0 /
      data fdmcs( 10,10) /   -.1134573260 212e0 /
      data fdmcs( 11,10) /    .0501830996 299e0 /
      data fdmcs( 12,10) /   -.0217756382 802e0 /
      data fdmcs( 13,10) /    .0092927921 972e0 /
      data fdmcs( 14,10) /   -.0039080375 664e0 /
      data fdmcs( 15,10) /    .0016223028 722e0 /
      data fdmcs( 16,10) /   -.0006656886 564e0 /
      data fdmcs( 17,10) /    .0002703262 001e0 /
      data fdmcs( 18,10) /   -.0001087475 762e0 /
      data fdmcs( 19,10) /    .0000433751 765e0 /
      data fdmcs( 20,10) /   -.0000171663 866e0 /
      data fdmcs( 21,10) /    .0000067455 409e0 /
      data fdmcs( 22,10) /   -.0000026333 256e0 /
      data fdmcs( 23,10) /    .0000010217 923e0 /
      data fdmcs( 24,10) /   -.0000003942 636e0 /
      data fdmcs( 25,10) /    .0000001513 391e0 /
      data fdmcs( 26,10) /   -.0000000578 112e0 /
      data fdmcs( 27,10) /    .0000000219 841e0 /
      data fdmcs( 28,10) /   -.0000000083 247e0 /
      data fdmcs( 29,10) /    .0000000031 398e0 /
      data fdmcs( 30,10) /   -.0000000011 798e0 /
      data fdmcs( 31,10) /    .0000000004 418e0 /
      data fdmcs( 32,10) /   -.0000000001 648e0 /
      data fdmcs( 33,10) /    .0000000000 613e0 /
c
c series for f40b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   1.68e-11
c                                         log weighted error  10.77
c                               significant figures required  10.47
c                                    decimal places required  11.36
c
      data fdbcs(  1,10) /    .8056894147 6e0 /
      data fdbcs(  2,10) /    .2834447403 7e0 /
      data fdbcs(  3,10) /    .0903522185 7e0 /
      data fdbcs(  4,10) /    .0111606016 1e0 /
      data fdbcs(  5,10) /    .0014164744 6e0 /
      data fdbcs(  6,10) /    .0000100892 5e0 /
      data fdbcs(  7,10) /    .0000022449 5e0 /
      data fdbcs(  8,10) /    .0000001741 4e0 /
      data fdbcs(  9,10) /   -.0000000486 2e0 /
      data fdbcs( 10,10) /   -.0000000054 1e0 /
      data fdbcs( 11,10) /    .0000000035 1e0 /
      data fdbcs( 12,10) /   -.0000000000 82e0 /
      data fdbcs( 13,10) /   -.0000000003 1e0 /
      data fdbcs( 14,10) /    .0000000000 81e0 /
      data fdbcs( 15,10) /    .0000000000 16e0 /
c
      data nscs / 21, 20, 19, 18, 17, 17, 16, 16, 15, 15 /
      data nmcs / 21, 22, 24, 25, 26, 28, 29, 31, 32, 33 /
      data nbcs / 26, 22, 20, 18, 18, 17, 16, 15, 15, 15 /
c
      data exp1, alfsml, alfbig, alfmax / 13*0.0 /
c
      if (exp1.ne.0.0) go to 20
      exp1 = 2.0/exp(1.0)
c
      alfsml = dlog (r1mach(1)/0.8863)
      alfbig = dexp (dmin1 (-dlog(r1mach(1)),dlog(r1mach(2)))-log(2.0))
c
      alfmax(1) = r1mach(2)
      alfmax(2) = r1mach(2)
      do 10 ndx=3,10
        xk = (ndx-2)*0.5
        alfmax(ndx) = r1mach(2)**(1.0/(xk+1.0)) * 0.99
 10   continue
c
 20   ndx = (xnu+1.0)*2.0 + 0.01
      fd = 0.0
      if (alpha.lt.alfsml) return
c
      if (ndx.ne.2) go to 30
c
c calculate the fermi-dirac function for xnu = 0
c
      if (alpha.lt.10.0) fd = alnrel (exp(alpha))
      if (alpha.ge.10.0) fd = alpha + alnrel (exp(-alpha))
      return
c
 30   if (alpha.ge.1.0) go to 40
      expalf = exp (alpha)
      fd = expalf * csevl (expalf*exp1-1.0, fdscs(1,ndx), nscs(ndx))
      return
c
 40   if (alpha.ge.4.0) go to 50
      fd = alpha**(xnu+1.0) * csevl ((alpha-2.5)/1.5, fdmcs(1,ndx),
     1  nmcs(ndx))
      return
c
 50   alfi = 0.0
      if (alpha.lt.alfbig) alfi = 1.0/alpha
      fd = alpha**(xnu+1.0) * csevl ((alfi-0.125)*8.0, fdbcs(1,ndx),
     1  nbcs(ndx))
      return
c
      end
 

      function alnrel(x)
      implicit real*8 (a-h,o-z)
c
c ****  description
c
c     alnrel(x) evaluates ln(1+x) accurately in the sense of relative
c     error when x is very small.  this routine must be used to
c     maintain relative error accuracy whenever x is small and
c     accurately known.
c
c series for alnr       on the interval -3.75000d-01 to  3.75000d-01
c                                        with weighted error   1.93e-17
c                                         log weighted error  16.72
c                               significant figures required  16.44
c                                    decimal places required  17.40


      dimension alnrcs(23)
      data alnrcs( 1) /   1.0378693562 743770e0 /
      data alnrcs( 2) /   -.1336430150 4908918e0 /
      data alnrcs( 3) /    .0194082491 35520563e0 /
      data alnrcs( 4) /   -.0030107551 12753577e0 /
      data alnrcs( 5) /    .0004869461 47971548e0 /
      data alnrcs( 6) /   -.0000810548 81893175e0 /
      data alnrcs( 7) /    .0000137788 47799559e0 /
      data alnrcs( 8) /   -.0000023802 21089435e0 /
      data alnrcs( 9) /    .0000004164 04162138e0 /
      data alnrcs(10) /   -.0000000735 95828378e0 /
      data alnrcs(11) /    .0000000131 17611876e0 /
      data alnrcs(12) /   -.0000000023 54670931e0 /
      data alnrcs(13) /    .0000000004 25227732e0 /
      data alnrcs(14) /   -.0000000000 77190894e0 /
      data alnrcs(15) /    .0000000000 14075746e0 /
      data alnrcs(16) /   -.0000000000 02576907e0 /
      data alnrcs(17) /    .0000000000 00473424e0 /
      data alnrcs(18) /   -.0000000000 00087249e0 /
      data alnrcs(19) /    .0000000000 00016124e0 /
      data alnrcs(20) /   -.0000000000 00002987e0 /
      data alnrcs(21) /    .0000000000 00000554e0 /
      data alnrcs(22) /   -.0000000000 00000103e0 /
      data alnrcs(23) /    .0000000000 00000019e0 /
      data nlnrel, xmin /0, 0./
c
c **** first executable statement  alnrel
c
      if (nlnrel.ne.0) go to 10
      nlnrel = inits (alnrcs, 23, 0.1*r1mach(3))
      xmin = -1.0 + sqrt(r1mach(4))
c
 10   continue
 
c
      if (abs(x).le.0.375) alnrel = x*(1. -
     1  x*csevl (x/.375, alnrcs, nlnrel))
      if (abs(x).gt.0.375) alnrel = dlog (1.0+x)
c
      return
      end



      function csevl(x,cs,n)
      implicit real*8 (a-h,o-z)
c
c ....description:
c     evaluate the n-term chebyshev series cs at x.  adapted from
c     r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973). also see fox
c     and parker, chebyshev polynomials in numerical analysis, oxford press
c     page 56.
c
c     input arguments --
c     x    value at which the series is to be evaluated.
c     cs   array of n terms of a chebyshev series.  in eval-
c          uating cs, only half the first coefficient is summed.
c     n    number of terms in array cs.
c     
       dimension cs(1)

       b1=0.
       b0=0.
       twox=2.*x
       do 10 i=1,n
       b2=b1
       b1=b0
       ni=n+1-i
       b0=twox*b1-b2+cs(ni)
 10    continue
c
       csevl = 0.5 * (b0-b2)
c
       return
      end



      function r1mach(i)
      implicit real*8 (a-h,o-z)
      real rmach(5)

c  single-precision machine constants
c
c  machine constants for the cray 1, cray x-mp
c  smallest useable real number
c      data rmach(1) / 200034000000000000000b /
c  largest useable real number
c      data rmach(2) / 577767777777777777776b /
c  smallest eps such that 1.+eps .ne. 1.
c      data rmach(3) / 377224000000000000000b /
c  2.*rmach(3)
c      data rmach(4) / 377234000000000000000b /
c  log10(2)
c      data rmach(5) / 377774642023241175720b /


c
c for  HP 730
c

      DATA RMACH(1) / Z'00800000' /
      DATA RMACH(2) / Z'7F7FFFFF' /
      DATA RMACH(3) / Z'33800000' /
      DATA RMACH(4) / Z'34000000' /
      DATA RMACH(5) / Z'3E9A209B' /
 

      r1mach = rmach(i)

      return
      end





      function inits(os,nos,eta)
      implicit real*8 (a-h,o-z)
      dimension os(nos)
c 
c .... initialize,orthogonal series,special function
c      initializes an orthogonal series so that it defines the
c     number of terms to carry in the series to meet a specified
c     error.
c
c     input arguments --
c     os     array of nos coefficients in an orthogonal series.
c     nos    number of coefficients in os.
c     eta    requested accuracy of series.


      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
 10   continue

 20   continue
      inits = i
c
      return
      end




c **************************************************************** c
c                                                                  c
c                                                                  c
c                        L O W V
c This routine calculates temperature depensent stopping number
c at low velocity limit.
c
      subroutine lowv(te,gamma,alpha,vion,vf,x2,vlh)
      implicit real*8 (a-h,o-z)
      dimension qz(100), wz(100)

      pi = 3.1415926
c
c ... compute g1
c
      half = -0.5
      g1 = -1./3.*ggg1(half,alpha)
      f1 = fd(half,alpha)

      nz = 100
      z1 = 0
      z2 = 1
      
      sumz = 0.
 1    call nodewt(z1,z2,qz,wz,nz)
      sumz1 = 0
      do 2 iz=1,nz
         g0=sqrt(te)/2*(f1 + 2*qz(iz)**2/te*g1)
         sumz1 = sumz1 + wz(iz)*u0fu(qz(iz),te,gamma,g0,x2)
 2    continue
      sumz = sumz + sumz1
      
      ratio = abs(sumz1/sumz)
      if(ratio.lt.1.0e-5) goto 5

      z1 = z2
      z2 = 5*z1
      goto 1
 
 5    vlh = 2./pi * (vion/vf)**3 * sumz

      return
      end

      function u0fu(z,te,gamma,g0,x2)
      implicit real*8 (a-h,o-z)

c
c ... compute df2/du
c
      pi = 3.1415926
      exz = (z**2 - gamma)/te
      if(exz.gt.30) then
         f2p = pi/2 * dexpw(-exz)
      else
         f2p = pi/2 /(1+dexpw(exz))
      endif
      u0fu = z**3 * f2p / (z**2+x2*g0)**2 

      return
      end



      subroutine highv(te,alpha,vion,vf,emass,hbar,wp, vlh)
      implicit real*8 (a-h,o-z)


      f12  = 1./2.
      f32  = 3./2.
      f52  = 5./2.

      v1 = log(2*emass*vion**2/(hbar*wp))
      v2 = te * fd(f32,alpha)/fd(f12,alpha)*(vf/vion)**2
      v3 = te**2/2 * fd(f52,alpha)/fd(f12,alpha)
     :   * (vf/vion)**4

      vlh = v1 - v2 - v3

      return
      end




      subroutine itpolt(n,x,y,xvalue,yvalue)
      implicit real*8 (a-h,o-z)	
      dimension x(10),y(10),u(10),v(10),z(10),w(10),b(10)

      call zspl3(n,x,y,w,b,u,v,z)
	
      yvalue = spl3(n,x,y,z,xvalue)
	
      return
      end



      subroutine zspl3(n,t,y,h,b,u,v,z)
      implicit real*8 (a-h,o-z)	
      dimension t(n),y(n),h(n),b(n),u(n),v(n),z(n)

      do 2 i = 1,n-1
         h(i) = t(i+1) - t(i)
         b(i) = (y(i+1) - y(i))/h(i)
 2    continue
      u(2) = 2.0*(h(1) + h(2))
      v(2) = 6.0*(b(2) - b(1))
      do 3 i = 3,n-1
         u(i) = 2.0*(h(i)+h(i-1)) - h(i-1)**2/u(i-1)
         v(i) = 6.0*(b(i)-b(i-1)) - h(i-1)*v(i-1)/u(i-1)
 3    continue
      z(n) = 0.0
      do 4 i = n-1,2,-1
         z(i) = (v(i)-h(i)*z(i+1))/u(i)
 4    continue
      z(1) = 0.0
      return
      end
      
      
      function spl3(n,t,y,z,x)
      implicit real*8 (a-h,o-z)		
      dimension t(n),y(n),z(n)
      
      do 2 i= n-1,2,-1
	 diff = x - t(i)
	 if(diff.ge.0.0) goto 3
 2    continue
      i=1
      diff = x - t(1)
 3    h = t(i+1) - t(i)
      b = (y(i+1)-y(i))/h - h*(z(i+1)+2.0*z(i))/6.0
      p = 0.5*z(i) + diff*(z(i+1)-z(i))/(6.0*h)
      p = b+ diff*p
      spl3 = y(i) + diff*p
      
      return
      end


c
c ********************************************************************c
c                                                                     c
c                           G1
c      This routine is used for calculating G1(alpha) integral


      function ggg1(sigma,alpha)
      implicit real*8 (a-h,o-z)
      dimension qn(100), qwn(100)
      
      if(alpha.gt.15) then
         q1 = -alpha
         if(q1.lt.-20) q1=-20
         q2 = 20 
         call nodewt(q1,q2,qn,qwn,50)
         sum = 0
         do 12 i=1,50
            sum = sum + qwn(i)*gg1(sigma,alpha,qn(i))
 12      continue
      else
         q1 = 0
         q2 = 30
         call nodewt(q1,q2,qn,qwn,50)
         sum = 0
         do 22 i=1,50
            sum = sum + qwn(i)*gg2(sigma,alpha,qn(i))
 22      continue
      endif
               
      ggg1 = sum

      return
      end


      
      function gg1(sigma,alpha,x)
      implicit real*8 (a-h,o-z)

      gg1 = (alpha+x)**sigma*dexpw(-x) / (1+dexpw(-x))**2

      return
      end

      
      function gg2(sigma,alpha,x)
      implicit real*8 (a-h,o-z)

      gg2 = x**sigma*dexpw(alpha-x) / (1+dexpw(alpha-x))**2

      return
      end



      function dexpw(x)
      implicit real*8 (a-h,o-z)
      
      y = x
      if(y.lt.-200) y=-200
      dexpw = dexp(y)

      return
      end

      function cubint (xj,m,x,y,nx)
 
      IMPLICIT  REAL*8 (a-h,o-z)
c  cubic lagrange interpolation to find y(xi). n is a location in the
c  table of x near xi. usually xi is between x(n-1) and x(n).
      dimension x(nx), y(nx)
      if (m .lt. 0) then
         do ie = 1,nx
            if (x(ie) .ge. xj) goto 5
         enddo
 5       continue
         m = ie-1
      endif
      n=m
      if (m.gt.nx) go to 10
      if (m.lt.2) go to 20
      if (m.eq.nx) n=m-1
      if (m.eq.2) n=m+1
      xi=xj
      amm=(xi-x(n-1))*(xi-x(n))*(xi-x(n+1))
      am=(xi-x(n-2))*(xi-x(n))*(xi-x(n+1))
      a=(xi-x(n-2))*(xi-x(n-1))*(xi-x(n+1))
      ap=(xi-x(n-2))*(xi-x(n-1))*(xi-x(n))
      xi=x(n-2)
      bmm=(xi-x(n-1))*(xi-x(n))*(xi-x(n+1))
      xi=x(n-1)
      if (bmm.eq.0D0) go to 10
      bm=(xi-x(n-2))*(xi-x(n))*(xi-x(n+1))
      xi=x(n)
      if (bm.eq.0D0) go to 10
      b=(xi-x(n-2))*(xi-x(n-1))*(xi-x(n+1))
      xi=x(n+1)
      if (b.eq.0D0) go to 10
      bp=(xi-x(n-2))*(xi-x(n-1))*(xi-x(n))
      if (bp.eq.0D0) go to 10
      cubint=y(n-2)*amm/bmm+y(n-1)*am/bm+y(n)*a/b+y(n+1)*ap/bp
      return

 10   CUBINT=Y(NX)
     :      + (y(nx-1)-y(nx))*(xj-x(nx))/(x(nx-1)-x(nx))
c     PRINT 5, N,X(N-2),X(N-1),X(N),X(N+1)
      RETURN
 20   CUBINT=Y(1) 
     :      + (y(2) - y(1))*(xj-x(1))/(x(2)-x(1))
c     PRINT 5, N,X(N-2),X(N-1),X(N),X(N+1)
      RETURN
                    entry xlint(xj,m,x,y,nx)
c      linear interpolation to find y(xj).
      n=m
      if (m.gt.nx) go to 10
      if (m.lt.2) go to 20
      b=x(n-1)-x(n)
      if (b.eq.0D0) go to 10
      cubint=(y(n-1)*(xj-x(n))+y(n)*(x(n-1)-xj))/b

      xlint = cubint

      return
c
   30 format (' trouble in cubint. set to y(n).',i5,1p4e15.7)
      end




      subroutine dexfit(iitt,iidd,zpj, ztg, nme, eamu, dedx,
     :                                                      a)
      implicit real*8 (a-h,o-z)
      parameter (mxep=100)

      dimension eamu(mxep), dedx(mxep), w(mxep), u(10), e(10), a(10)

      common / oldat / olda(10)
      common / fitss / np, x(mxep), y(mxep)
      common / hlkey / key
      common / hlbdy / elow, ehigh, epeak

c --------------------------------------------------------------------
c .... cutoff boundary of 'low' and 'high' energy (in MeV/amu)
c      low energies: below 30 keV/amu; high energies: above 1 MeV/amu
 
c      data elow,ehigh/0.01, 1.00/
c ---------------------------------------------------------------------      

c
c ... find the peak of dE/dX
      epeak = 0
      dpeak = 0
      do 1 i=1,nme
         if(dedx(i).gt.dpeak) then
            epeak = eamu(i)
            dpeak = dedx(i)
         endif
 1    continue
      elow  = min(0.1,0.1*epeak)
      ehigh = max(0.5,10.*epeak) 

c
c ... first of all, fit low energy data     

      ilow = 0
      do 3 i=1,nme
         if(eamu(i).le.elow) then
            ilow = ilow + 1
            x(ilow) = eamu(i)
            y(ilow) = dedx(i)
         endif
 3    continue
      if(ilow.lt.1) then
         a(1) = 0
         a(2) = 0
         goto 1002
      endif

c ... fit the data
      key = 1
      np = ilow
      n = 2
      do 4 i=1,10
         e(i) = 1.e-3
 4    continue
      ef = 0.
      escale = 0.07 / e(1)
      iprint = 1
      icon = 1
      maxit = 70
      lunbtm = 18
         
      u(1) = 50
      u(2) = 0.5
      
      call botm ( u,e,n,ef,escale,iprint,icon,maxit,w,lunbtm )
      
      a(1) = u(1)
      a(2) = u(2)

c     
c ... now, fit data between (elow, ehigh)     

 1002 ihi = 0
      do 5 i=1,nme
         if(eamu(i).gt.elow .and. eamu(i).lt.ehigh) then
            ihi = ihi + 1
            x(ihi) = eamu(i)
            y(ihi) = dedx(i)
         endif
 5    continue
      if(ihi .lt.1) then
         a(3) = 0
         a(4) = 0
         a(5) = 0
         a(6) = 0
         a(7) = 0
         goto 1003
      endif

         
c ... fit the data
      key = 3
      np = ihi
      n = 5
      do 6 i=1,10
         e(i) = 1.0e-3
 6    continue
      ef = 0.
      escale = 0.07 / e(1)
      iprint = 1
      icon = 1
      maxit = 70
      lunbtm = 18
      
      u(1) = 10
      u(2) = 0.5
      u(3) = 10
      u(4) = 10
      u(5) = 1
      if(iidd.gt.1 .and. a(3).ne.0) then
         u(1) = a(3)
         u(2) = a(4)
         u(3) = a(5)
         u(4) = a(6)
         u(5) = a(7)
      endif
      if(iidd.eq.1 .and. iitt.gt.1 .and. olda(3).ne.0) then
         u(1) = olda(3)
         u(2) = olda(4)
         u(3) = olda(5)
         u(4) = olda(6)
         u(5) = olda(7)
      endif
      
      call botm ( u,e,n,ef,escale,iprint,icon,maxit,w,lunbtm )
      
      a(3) = u(1)
      a(4) = u(2)
      a(5) = u(3)
      a(6) = u(4)
      a(7) = u(5)
c     
c ... now, fit high energy  data     

 1003 ihi = 0
      do 9 i=1,nme
         if(eamu(i).ge.ehigh) then
            ihi = ihi + 1
            x(ihi) = eamu(i)
            y(ihi) = dedx(i)
         endif
 9    continue
      if(ihi.lt.1) then
         a(8) = 0
         a(9) = 0
         a(10)= 0
         goto 1004
      endif

c ... fit the data
      key = 2
      np = ihi
      n = 3
      do 10 i=1,10
         e(i) = 1.e-3
 10   continue
      ef = 0.
      escale = 0.07 / e(1)
      iprint = 1
      icon = 1
      maxit = 70
      lunbtm = 18
      
      u(1) = 10
      u(2) = 1
      u(3) = 1
      if(iidd.gt.1 .and. a(8).ne.0) then
         u(1) = a(8)
         u(2) = a(9)
         u(3) = a(10)
      endif
      if(iidd.eq.1 .and. iitt.gt.1 .and.olda(8).ne.0) then
         u(1) = olda(8)
         u(2) = olda(9)
         u(3) = olda(10)
      endif
      
      call botm ( u,e,n,ef,escale,iprint,icon,maxit,w,lunbtm )
      
      a(8) = u(1)
      a(9) = u(2)
      a(10)= u(3)
      
      if(iidd.eq.1) then
         do 11 i=1,10
            olda(i) = a(i)
 11      continue
      endif
      
 1004 return
      end



c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      subroutine calcfx ( n,u,f )
      implicit double precision (a-h,o-z)
      parameter (mxep=100)

      common / fitss / np, x(mxep), y(mxep)
      dimension u(10)
 
      sum = 0.
      do 100 i=1,np
 
         fit = fitfun ( x(i), u )
c        err = fit / y(i) - 1.
         err = fit - y(i)
         sum = sum + err*err
 
 100  continue
 
      f = sum

      return
      end

c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      function fitfun ( x, u )
      implicit double precision (a-h,o-z)
      dimension u(10)
      common /hlkey/key
c
c ... key=1: fitting low energy data
c     key=2: fitting high energy data      

      if(key.eq.1) then
        fitfun = u(1)*x**u(2)                
      endif
      if(key.eq.2) then
         fitfun = u(1)/x * log(1.0+u(2)/x + u(3)*x)
      endif
      if(key.eq.3) then
         sl = u(1)*x**u(2)
         factor = 1.0+u(4)/x + u(5)*x
         if(factor.gt.0) sh = u(3)/x * log(factor)
         if(factor.le.0) sh = 1.0e30
         if(sl+sh .ne. 0) then
            fitfun = sl*sh/(sl+sh)
         else
            fitfun = 0
         endif
      endif
      return
      end
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 
      SUBROUTINE BOTM ( X,E,N,EF,ESCALE,IPRINT,ICON,MAXIT,W,LUNBTM )
      implicit real*8 (a-h,o-z)

      DIMENSION W(100),X(10),E(10)

      DDMAG = 0.1*ESCALE
      SCER = 0.05/ESCALE
      JJ = N*(N+1)
      JJJ=JJ+N
      K=N+1
      NFCC=1
      IND=1
      INN=1
      DO 4 I=1,N
        W(I)=ESCALE
      DO 4 J=1,N
        W(K)=0.
        IF (I-J) 4,3,4
    3   W(K)=ABS(E(I))
    4   K=K+1
      ITERC=1
      ISGRAD=2
      CALL CALCFX ( N,X,F )
      FKEEP=2.*ABS(F)
    5 ITONE=1
      FP=F
      SUM=0.
      IXP=JJ
      DO 6 I=1,N
        IXP=IXP+1
    6   W(IXP)=X(I)
      IDIRN=N+1
      ILINE=1
    7 DMAX=W(ILINE)
      DACC=DMAX*SCER
      DMAG=MIN(DDMAG,0.1*DMAX)
      DMAG=MAX(DMAG,20.*DACC)
      DDMAX=10.*DMAG
      GO TO (70,70,71),ITONE
   70 DL=0.
      D=DMAG
      FPREV=F
      IS=5
      FA=FPREV
      DA=DL
    8 DD=D-DL
      DL=D
   58 K=IDIRN
      DO 9 I=1,N
        X(I)=X(I)+DD*W(K)
    9   K=K+1
      CALL CALCFX ( N,X,F )
      NFCC=NFCC+1
      GO TO (10,11,12,13,14,96),IS
   14 IF (F-FA) 15,16,24
   16 IF (ABS(D)-DMAX) 17,17,18
   17 D=D+D
      GO TO 8
 18   continue
c  18 WRITE (LUNBTM,19)
   19 FORMAT (5X,'MAXIMUM CHANGE DOES NOT ALTER FUNCTION')
      GO TO 20
   15 FB=F
      DB=D
      GO TO 21
   24 FB=FA
      DB=DA
      FA=F
      DA=D
   21 GO TO (83,23),ISGRAD
   23 D=DB+DB-DA
      IS=1
      GO TO 8
   83 D=0.5*(DA+DB-(FA-FB)/(DA-DB))
      IS=4
      IF ((DA-D)*(D-DB)) 25,8,8
   25 IS=1
      IF (ABS(D-DB)-DDMAX) 8,8,26
   26 D=DB+SIGN(DDMAX,DB-DA)
      IS=1
      DDMAX=DDMAX+DDMAX      
      DDMAG=DDMAG+DDMAG
      ddmag = min(ddmag,1.0e10)
      IF (DDMAX-DMAX) 8,8,27
   27 DDMAX=DMAX
      GO TO 8
   13 IF (F-FA) 28,23,23
   28 FC=FB
      DC=DB
   29 FB=F
      DB=D
      GO TO 30
   12 IF (F-FB) 28,28,31
   31 FA=F
      DA=D
      GO TO 30
   11 IF (F-FB) 32,10,10
   32 FA=FB
      DA=DB
      GO TO 29
   71 DL=1.
      DDMAX=5.
      FA=FP
      DA=-1.
      FB=FHOLD
      DB=0.
      D=1.
   10 FC=F
      DC=D
   30 A=(DB-DC)*(FA-FC)
      B=(DC-DA)*(FB-FC)
      IF ((A+B)*(DA-DC)) 33,33,34
   33 FA=FB
      DA=DB
      FB=FC
      DB=DC
      GO TO 26
   34 D=0.5*(A*(DB+DC)+B*(DA+DC))/(A+B)
      DI=DB
      FI=FB
      IF (FB-FC) 44,44,43
   43 DI=DC
      FI=FC
   44 GO TO (86,86,85),ITONE
   85 ITONE=2
      GO TO 45
   86 IF (ABS(D-DI)-DACC) 41,41,93
   93 IF (ABS(D-DI)-0.03*ABS(D)) 41,41,45
   45 IF ((DA-DC)*(DC-D)) 47,46,46
   46 FA=FB
      DA=DB
      FB=FC
      DB=DC
      GO TO 25
   47 IS=2
      IF ((DB-D)*(D-DC)) 48,8,8
   48 IS=3
      GO TO 8
   41 F=FI
      D=DI-DL
      DD=SQRT((DC-DB)*(DC-DA)*(DA-DB)/(A+B))
      DO 49 I=1,N
        X(I)=X(I)+D*W(IDIRN)
        W(IDIRN)=DD*W(IDIRN)
   49   IDIRN=IDIRN+1
      W(ILINE)=W(ILINE)/DD
      ILINE=ILINE+1
      IF (IPRINT-1) 51,50,51
 50   continue
c  50 WRITE (LUNBTM,52) ITERC,NFCC,F,(X(I),I=1,N)
   52 FORMAT (/' ITERATION',I5,I15,'FUNCTION VALUES',10X,'F =',D15.8/
     &       (2X,8D16.8))
      GO TO (51,53),IPRINT
   51 GO TO (55,38),ITONE
   55 IF (FPREV-F-SUM) 94,95,95
   95 SUM=FPREV-F
      JIL=ILINE
   94 IF (IDIRN-JJ) 7,7,84
   84 GO TO (92,72),IND
   92 FHOLD=F
      IS=6
      IXP=JJ
      DO 59 I=1,N
        IXP=IXP+1
   59   W(IXP)=X(I)-W(IXP)
      DD=1.
      GO TO 58
   96 GO TO (112,87),IND
  112 IF (FP-F) 37,37,91
   91 D=2.*(FP+F-2.*FHOLD)/(FP-F)**2
      IF (D*(FP-FHOLD-SUM)**2-SUM) 87,37,37
   87 J=JIL*N+1
      IF (J-JJ) 60,60,61
   60 DO 62 I=J,JJ
        K=I-N
   62   W(K)=W(I)
      DO 97 I=JIL,N
   97   W(I-1)=W(I)
   61 IDIRN=IDIRN-N
      ITONE=3
      K=IDIRN
      IXP=JJ
      AAA=0.
      DO 67 I=1,N
        IXP=IXP+1
        W(K)=W(IXP)
        IF (AAA-ABS(W(K)/E(I))) 66,67,67
   66   AAA=ABS(W(K)/E(I))
   67   K=K+1
      DDMAG=1.
      W(N)=ESCALE/AAA
      ILINE=N
      GO TO 7
   37 IXP=JJ
      AAA=0.
      F=FHOLD
      DO 99 I=1,N
        IXP=IXP+1
        X(I)=X(I)-W(IXP)
        IF (AAA*ABS(E(I))-ABS(W(IXP))) 98,99,99
   98   AAA=ABS(W(IXP)/E(I))
   99   CONTINUE
      GO TO 72
   38 AAA=AAA*(1.+DI)
      GO TO (72,106),IND
   72 IF (IPRINT-2) 53,50,50
   53 GO TO (109,88),IND
  109 IF (AAA-0.1) 89,89,76
   89 GO TO (20,116),ICON
  116 IND=2
      GO TO (100,101),INN
  100 INN=2
      K=JJJ
      DO 102 I=1,N
        K=K+1
        W(K)=X(I)
  102   X(I)=X(I)+10.*E(I)
      FKEEP=F
      CALL CALCFX ( N,X,F )
      NFCC=NFCC+1
      DDMAG=0.
      GO TO 108
   76 IF (F-FP) 35,78,78
 78   continue
c  78 WRITE (LUNBTM,80)
   80 FORMAT (5X,'ACCURACY LIMITED BY ERRORS IN F')
      GO TO 20
   88 IND=1
   35 DDMAG=0.4*SQRT(FP-F)
      ISGRAD=1
  108 ITERC=ITERC+1
      IF (ITERC-MAXIT) 5,5,81
 81   continue
c  81 WRITE (LUNBTM,82) MAXIT
   82 FORMAT (I5,' ITERATIONS COMPLETED BY BOTM')
      IF (F-FKEEP) 20,20,110
  110 F=FKEEP
      DO 111 I=1,N
        JJJ=JJJ+1
  111   X(I)=W(JJJ)
      GO TO 20
  101 JIL=1
      FP=FKEEP
      IF (F-FKEEP) 105,78,104
  104 JIL=2
      FP=F
      F=FKEEP
  105 IXP=JJ
      DO 113 I=1,N
        IXP=IXP+1
        K=IXP+N
        GO TO (114,115),JIL
  114   W(IXP)=W(K)
        GO TO 113
  115   W(IXP)=X(I)
        X(I)=W(K)
  113   CONTINUE
      JIL=2
      GO TO 92
  106 IF (AAA-0.1) 20,20,107
   20 EF=F
      RETURN
  107 INN=1
      GO TO 33
      END



      subroutine getrange(kzz,it,id, mep, ee, elog, fdedx)
      implicit real*8 (a-h,o-z)
      parameter (npx=1500)
      
      dimension ee(npx), elog(npx), fdedx(npx)
      dimension eq(100), wq(100)
      common / hlbdy / elow, ehig, epeak      
      common /rrang/range(50,50,200)

      do 100 iee=1,mep
         emax = ee(iee)

         rag = 0.0
         je = 0
         
         eh = emax*1.0e6
         el = 0.1 * eh
         np = 30

 15      call nodewt(el,eh,eq,wq,np)
         
         add = 0
         do 150 ie=np,1,-1
            
c ----------------------- dE/dX
            eion = eq(ie)/1.0e6
c            if(eion.le.elow) then
c               dedx = a(1)*eion**a(2)
c            endif
c            if(eion.gt.elow .and.eion.lt.ehig) then
c               sl = a(3)*eion**a(4)
c               sh = (a(5)/eion)*log(1.0+a(6)/eion+a(7)*eion)
c               dedx = sl*sh/(sl+sh)
c            endif
c            if(eion.ge.ehig) then
c               dedx = (a(8)/eion)*log(1.0+a(9)/eion+a(10)*eion)
c     endif
            mx = -1
            dedx = cubint(dlog(eion), mx, elog, fdedx, mep)
            dedx = 1.0e-15 * dedx
c -----------------------
            if(dedx.le.0.) then
               see = 1
            endif
            add = add +  wq(ie)/dedx
 150     continue

         rag = rag + add

         ratio = abs(add/rag)
         if(ratio .gt.1.0e-3 .and. el.gt.101) then
            eh = el
            el = 0.1*eh
            if(el.lt.100) el=100
            goto 15
         endif
c
c ... for low energy region : using analytical expression
c -------------------------------------------------------

         c = dedx/sqrt(eion)
         rag = rag + 2*sqrt(eion)/c
         
c
c .... for gold and Aluminum, convert to mg/cm2
c     ---------------------------------
c         aucoef = 6.023e23/amass * 1d-3
c         rag = rag/aucoef         
         range(it,id,iee) = rag

 100  continue


      return
      end









