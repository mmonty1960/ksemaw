**********    kSEMAW IbridOne sources   ***************
*
*   hereinafter subroutines used in IbridOne method
C
C
C   SPECTRO-ELLIPSOMETRIC MEAUREMENT ANALYSIS WORKBENCH
C
C  Author: Marco Montecchi
C  Department of Energy Technologies
C  ENEA C.R. Casaccia
C  Roma - Italy
C
C  Copyright (C) 2020  Marco Montecchi
C
C   This program is free software: you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation, either version 3 of the License, or
C   (at your option) any later version.
C
C   This program is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program.  If not, see <https://www.gnu.org/licenses/>.
C
      
      
      
      SUBROUTINE FITSELMQ(rc)
      PARAMETER(nparmax=17) ! Max numero di parametri del fit
      INTEGER IXW(20),s1p2,dato(14),PPM(nparmax),mwl,md,mda
      INTEGER ieon,jie,jiemax,nmisure
      INTEGER, allocatable :: IWA(:)
      REAL, allocatable :: P(:)
      REAL, allocatable :: FVEC(:)
      REAL, allocatable :: WA(:)
      REAL, allocatable :: fjac(:,:) 
      REAL SOL(999,6),PF(7,21),PM(200,5),MIS(17,201,2),ELI(8,201,2),
     *     PAR(60,5),CNK(16,3),RXY(30,4),VOT(5,2),ARSE(500,2)
!      REAL ENORM,SPMPAR,TOL,FNORM
      REAL cov(nparmax,nparmax),cinv(nparmax,nparmax),
     +     PS(200),wwl(999),nn(999),kk(999),kklog(999),rn(201),rp(201),
     +     r1(201),Pot(200),ncent(999),ndev(999),kcent(999),kdev(999),
     +     MISS(5,201),pst(nparmax,2),e1(999),e2(999)
      DOUBLE PRECISION vs(99),k,dT
      CHARACTER R*2,NANK(16)*256,SPM(200)*17,SFMT(200)*3,ST3*3,
     *          SFU(9)*24,st63*63,rc*2
      COMMON /VALORI/ MIS,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      EXTERNAL FSQ,FRCK,DELTAT
      SAVE Pot,ifirstcall,chi2
      DATA ifirstcall,chi2,icol,chi2min,jobtot,ifit 
     +    / 0, -1.0, 0, 1.e+36, 0, 0/

     
      length=len_trim(nank(16))
**** inizializzazione stringhe
      sfu(1)='Lorentz'
      sfu(2)='Quant-omo'
      sfu(3)='Quant-inomo' 
      sfu(4)='Flat'
      sfu(5)='Drude'
      sfu(6)='Indirect Gap'
      sfu(7)='Direct Gap'
      sfu(8)='M0M3-Cody'
      sfu(9)='M0M3-Tauc'
**** impostazione pannello di fit e puntatore PM
      npp=nint(par(34,5))
      do i=1,npp
       ppm(i)=nint(par(35+i,5))
      end do
**** definizione pannello interfaccia fit      
      call PANSTRING(SPM,SFMT)
**** controllo allocazione
      iAlloIWA=0
      iAlloP=0
      iAlloFVEC=0
      iAlloWA=0
      iAlloFJAC=0
**** conteggio dati
      mwl=201 !relativo a MIS ed ELI
      mwla=0  !wl misure abilitate
      md=nint(sol(1,1)) !relativo a sol-nk
      mda=0 !dati nk abilitati
      ial=0
      par(38,1)=1.e+20
      par(38,2)=0.
      do i=2,md+1
       if(nint(sol(i,6)).gt.0) then
        if(sol(i,1).ge.par(4,1).and.sol(i,1).le.par(4,2)) then
          mda=mda+1
          par(38,1)=MIN(par(38,1),sol(i,1))
          par(38,2)=MAX(par(38,2),sol(i,1))
        else
         sol(i,6)=0.
        end if
       end if
       if(sol(i,4).le.0.) then
        if(ial.eq.0) then
         write(*,*)'ERRn = 0 @',sol(i,1),' => ERRn = Dn =',par(21,3)
         ial=1
        end if
        sol(i,4)=abs(par(21,3))
       end if
      end do
      par(38,3)=1.e+20
      par(38,4)=0.
      do i=1,201
       if(nint(mis(17,i,1)).eq.1) then
        mwla=mwla+1
        par(38,3)=MIN(par(38,3),mis(7,i,1))
        par(38,4)=MAX(par(38,4),mis(7,i,1))
       end if
      end do
**** inizializzazione parametri di fit
      n=nint(par(35,5))
      if(iAlloP.eq.0) then
       ALLOCATE(p(n),stat=ierror)
       if (ierror.ne.0) then
        print*,'error: could not allocate memory for p(',n,')'
        iAlloP=0
       else
        iAlloP=1
       endif
      end if
      do i=1,n
       ipm=nint(pm(i,3))
       p(i)=pm(ipm,1)
       write(*,*)i,' ipm = ',ipm,' p = ',p(i),pm(ipm,2)
      end do
      if(ifirstcall.eq.0) then !settaggio alla prima chiamata
       do i=1,200
        Pot(i)=pm(i,1)
       end do
       TOL = SQRT(SPMPAR(1))
       par(25,1)=tol
      end if
**** caricamento parametri utili da PAR
      do i=1,14
       dato(i)=nint(par(i,4))
      end do
      ioptf=nint(pm(100,1))
c      jobtot=nint(par(8,1))
      jobview=nint(par(8,2))
      write(*,*) 'jobview=',jobview
      LEV=nint(RXY(25,3))
      chi2fin=par(55,2)
      fredeg=par(56,2)
**** inizializzazione minimale
      cnk(1,1)=0. !nk incogniti
      te=0.  !angolo
      imR=3  !Rn
      iwr=13 !Rn
      s1p2=nint(par(27,2)) !polarizzazione
      tol=par(25,1)
      nmisure=0
      do i=1,5
       if(dato(i).eq.2) nmisure=nmisure+1 
      end do
      jie=0
      jiemax=2**nmisure
c      icol=0
      ieon=0 !calcolo errore su OFF
c      chi2min=1.e+36
!       write(*,*)
!       write(*,*)'******************************************************'
!       write(*,*)'                 Subroutine FITSELMQ'
!       write(*,*)'     1)  fit n - k con SUM f_j'
!       write(*,*)'     2)  Ibridone: k da T & BF di R con n da SUM f_j'
!       write(*,*)
!       write(*,*)'           f_j(eV)=epr_j(eV)+i*epi_j(eV):'
!       write(*,*)'            - classica'
!       write(*,*)'            - quantistica omogenea'
!       write(*,*)'            - quantistica inomogenea'
!       write(*,*)'            - flat (1+i*1)'
!       write(*,*)
!       write(*,*)'    n = sqrt(sqrt(epr*epr+epi*epi)+epr)/sqrt(2.)'
!       write(*,*)
!       write(*,*)'    k = sqrt(sqrt(epr*epr+epi*epi)-epr)/sqrt(2.)'
!       write(*,*)'      if K_j >= 0 epi=SUM K_j*epi_j'
!       write(*,*)
!       write(*,*)'          precisione del calcolo = ',spmpar(1)
!       write(*,*)
!       write(*,*)'    N. wl abilitate  = ',mwla,' su ',mwl,' disponibili'
!       write(*,*)'    N. nk abilitati  = ',mda,' su ',md,' disponibili'
!       write(*,*)'     parametri di fit = ',n
!       write(*,*)'******************************************************'
      
***** loop azione
      write(*,*)'|IbridOne> rc=',rc 
      r='? '
      DO WHILE(R.NE.'X ')
       st63=
     + 'N.  Fit  parametro             valore      ERR(99%)  cor. glob.'
*** controllo ciclo ibridone con errore
       if(ieon.eq.1) then
        jie=jie+1
        if(jie.gt.jiemax) then !stop
         ieon=2
         do ii=1,201
          do i=1,5
           mis(i,ii,1)=miss(i,ii) !ripristino valori misura 
          end do
         end do
         do i=1,n
          ip=nint(pm(i,3))
          p(i)=pst(i,1)          !valore centrale (jie=0)
          pm(ip,1)=p(i)
          if(pst(i,2).gt..0) then
           pm(ip,4)=sqrt(pst(i,2)) !errore = rms deviazioni
          else
           pm(ip,4)=0.
          end if
         end do
         do i=1,201
*** memorizza le soluzioni centrali
          nn(i)=ncent(i)
          kk(i)=kcent(i)
          mis(16,i,1)=nn(i)
          mis(16,i,2)=kk(i)
          if(ndev(i).gt.0.) then
           ndev(i)=sqrt(ndev(i))
          else
           ndev(i)=0.
          end if
          if(kdev(i).gt.0.) then
           kdev(i)=sqrt(kdev(i))
          else
           kdev(i)=0.
          end if
         end do
         st63=
     + 'N.  Fit  parametro             valore      ERR(rms)  cor. glob.'
        end if
       end if
                      
       if(ieon.eq.0) then
        st3(1:2)=rc
       else if(ieon.eq.1) then
        st3='i  '
       else if(ieon.eq.2) then
        st3='M  '
       end if
       r=st3(1:2)
             
!  2     if(r.eq.'h ') then
!         WRITE(*,*)'"0 " azzera le soluzioni memorizzate'
!         write(*,*)'"e " edit abilitazione_fit n'
!         write(*,*)'"em" edit abilitazione_fit R'
!         write(*,*)'"eR" abilita tutte le wl per fit R'
!         write(*,*)'"eN" abilita tutti gli NK per fit n'
!         write(*,*)'"g " grafica curva con FT# di lavoro'
!         write(*,*)'"gi" grafica curva con n-k ibridone'
!         write(*,*)'"G " refresh grafici mantenendo i range'
!         write(*,*)'"M " ibridone con errori e mantieni nk'
c!         write(*,*)'"R " riporta i parametri alla situazione migliore'
!         write(*,*)'"s " salva/carica file.nk'
c!         write(*,*)'"- " recedi a Job# precedente'
c!         write(*,*)'"+ " torna a Job# successivo'
!        end if

       
       if(st3(1:2).eq.'G ') then
*** refresh grafici
        do i=3,5
         if(dato(i).eq.2) then
          imR=i
          iwr=10+i
          CALL PLOTMIS(imR,iwr,3)
         end if
        end do
        CALL PLOTNK(1,3)
        CALL PLOTNK(2,3)
        if(nint(par(10,1)).eq.1) then
         CALL PLOTNK(3,3) !plot epsi1
         CALL PLOTNK(4,3) !plot epsi2
        end if
       end if
       
       
       if(r.eq.'M '.and.ieon.eq.0) then
*** impostazione ibridone con calcolo errore
        ieon=1
        jie=0
        if(nmisure.ge.1) then
         r='i ' !ibridone
         do ii=1,201
          do i=1,5
           miss(i,ii)=mis(i,ii,1) !storage valori
          end do
         end do
        else
         r='  '    !annulla ibridone con calcolo errore
         st3='   '
        end if
       end if

       if(r.eq.'i '.and.ieon.eq.1.and.jie.gt.0) then
*** alterazione misure sperimentali
        imisura=0
        write(*,*)'..... Ibridone con calcolo errore: passo',jie,
     +  ' di',jiemax
        do i=1,5
         if(dato(i).eq.2) then
          iespo=int(jie/2.**imisura)
          write(*,*)'  misura =',i,' (-1)^espo =',(-1.)**iespo
          do ii=1,201
           mis(i,ii,1)=miss(i,ii)+mis(i,ii,2)*(-1.)**iespo
          end do
          imisura=imisura+1
         end if
        end do 
       end if

       if(r.eq.'f '.or.r.eq.'i ') then 
*** lancio fit di n o di ibridone
        s1p2=nint(par(27,2))
        if(r.eq.'i ') then !impostazione calcolo IBRIDO
         ibridok=0
         if(dato(1).eq.2) then
          write(*,*)'-> lancio IBRIDONE: calcolo k da Tn'
          ibridok=1          
         else if(dato(2).eq.2) then
          write(*,*)'-> lancio IBRIDONE: calcolo k da Tp teta =',
     +     par(6,1),' deg  S1P2 =',s1p2
          ibridok=1
         else if(dato(1).ne.2.and.dato(2).ne.2.and.
     +           nint(par(32,2)).eq.2) then
          write(*,*)'-> lancio IBRIDONE: k imposto = selmq_k'
          ibridok=1
         end if
         if(dato(3).eq.2) then
          write(*,*)'                     fit di Rn'
          ibridok=ibridok+1
         end if
         if(dato(4).eq.2) then
          write(*,*)'                     fit di Rp  teta =',
     +     par(6,1),' deg  S1P2 =',s1p2
          ibridok=ibridok+1
         end if
         if(dato(5).eq.2) then
          write(*,*)'                     fit di R1'
          ibridok=ibridok+1
         end if
         if(ibridok.lt.2) then
          write(*,*)'Attenzione IBRIDONE richiede R&T o R&k_imposto!'
          r='  '!annulla calcolo IBRIDONE
         else
          m=(ibridok-1)*mwla
          write(*,*)'     ... per un totale di',m,' punti e ',n,
     +              ' parametri'
         end if
        else if(r.eq.'f ') then ! elimina p_fit ><  SELMQ 
         if(nint(par(32,5)).eq.0) then
          m=mda
         else
          m=2*mda
         end if
         npfit=n
         do ipf=npfit,1,-1
          ipm=nint(pm(ipf,3))
          write(*,*)'ipf = ',ipf,' ipm = ',ipm
          if(ipm.lt.100) then
           pm(ipm,2)=.0
           do iii=ipf,n
            pm(iii,3)=pm(iii+1,3) !sposto di un posto in basso
            pm(iii,4)=pm(iii+1,4)
            pm(iii,5)=pm(iii+1,5)
           end do
           n=n-1
           write(*,*)'n = ',n
          end if
          par(35,5)=n
         end do
        end if
        
****Fitta, calcola matrici jac & cov, err. e corr.
        ifit=1
        OPEN(1,STATUS='UNKNOWN',FILE='temp/semaw.log')
        if(iAlloIWA.eq.0) then
         ALLOCATE(IWA(n),stat=ierror)
         if (ierror.ne.0) then
          print*,'error: could not allocate memory for IWA(',n,')'
          iAlloIWA=0
         else
          iAlloIWA=1
         endif
        end if
        if(iAlloFVEC.eq.0) then
         ALLOCATE(FVEC(m),stat=ierror)
         if (ierror.ne.0) then
          print*,'error: could not allocate memory for FVEC(',m,')'
          iAlloFVEC=0
         else
          iAlloFVEC=1
         endif
        end if
        if(iAlloWA.eq.0) then
         lwa=M*N+5*N+M
         ALLOCATE(WA(lwa),stat=ierror)
         if (ierror.ne.0) then
          print*,'error: could not allocate memory for WA(',lwa,')'
          iAlloWA=0
         else
          iAlloWA=1
         endif
        end if
        if(iAlloFJAC.eq.0) then
         ALLOCATE(FJAC(m,n),stat=ierror)
         if (ierror.ne.0) then
          print*,'error: could not allo mem. for FJAC(',m,',',n,')'
          iAlloFJAC=0
         else
          iAlloFJAC=1
         endif
        end if
        write(*,*) 'n= ',n,' m= ',m,' lwa= ',lwa
        if(r.eq.'f '.and.n.gt.0) then
         CALL LMDIF1(FSQ,M,N,P,FVEC,TOL,INFO,IWA,WA,LWA)
        else if(r.eq.'i '.and.n.gt.0) then
         CALL LMDIF1(FRCK,m,N,P,FVEC,TOL,INFO,IWA,WA,LWA)
        end if
        CLOSE(1)
        
**** analisi chi2
        if((r.eq.'f '.or.(r.eq.'i '.and.ieon.eq.0)).and.n.gt.0) then
         OPEN(1,STATUS='OLD',FILE='temp/semaw.log',IOSTAT=IOS,ERR=19)
         read(1,*) chi2ini
         chi2fin=chi2ini
         i=1
         imin=1
         do while(ios.eq.0)
          read(1,*,IOSTAT=IOS,ERR=18) chi2
          i=i+1
          if(chi2.lt.chi2fin) then
           chi2fin=chi2
           imin=i
          end if
 18       continue
         end do
         REWIND(1)
         do i=1,imin
          read(1,*,IOSTAT=IOS,ERR=19) chi2fin,(p(j),j=1,n)
         end do
         CLOSE(1)
         par(55,2)=chi2fin
 19      continue
        end if
        
**** imponi valori di BF-MIN
        FNORM = ENORM(M,FVEC)
        do i=1,n
         ip=nint(pm(i,3))
         pm(ip,1)=p(i)
         pm(ip,4)=.0
        end do
        
        if(m.gt.0) then 
         
*** calcolo matrice jabobiana fjac
         do i=1,200 !salvataggio situazione di best fit
          ps(i)=pm(i,1)
         end do
         iflag=2
         epsfcn=0.0e0
         if(r.eq.'f '.and.n.gt.0) then
          call fdjac2(fsq,m,n,P,fvec,fjac,m,iflag,epsfcn,wa)
         else if(r.eq.'i '.and.n.gt.0) then
          call fdjac2(frck,m,n,P,fvec,fjac,m,iflag,epsfcn,wa)
         end if
         do i=1,200 !ripristino situazione di best fit
          pm(i,1)=ps(i)
         end do
         do i=1,n
          p(i)=pm(nint(pm(i,3)),1)
         end do 
         
*** calcolo matrice covarianza inversa (fjac_trasposta)*(fjac)
         do i=1,n
          do j=1,n
           cinv(i,j)=.0
           do i2=1,m
            cinv(i,j)=cinv(i,j)+fjac(i2,i)*fjac(i2,j)
           end do
          end do
         end do
         
*** calcolo matrice varianza-covarianza
         call matinv(n,nparmax,cinv,cov)
         
         
*** calcolo errori sui singoli parametri al 99% di confidenza
         do i=1,n 
          if(cov(i,i).lt.1.e3) then
           pm(nint(pm(i,3)),4)=sqrt(abs(6.635*cov(i,i)))
          end if
         end do
         
*** calcolo global correlation coefficients
         if(n.gt.1) then
          do i=1,n
           if(abs(cov(i,i)*cinv(i,i)).gt.0.) then
            pm(nint(pm(i,3)),5)=sqrt(1.-1./(cov(i,i)*cinv(i,i)))
           else
            pm(nint(pm(i,3)),5)=1.
           end if
          end do
c         else
c          pm(nint(pm(i,3)),5)=0.0000
         end if
        end if !(m.gt.0)
        
       end if !(r.eq.'f '.or.r.eq.'i ')
             
       if(r.eq.'f '.or.r.eq.'i '.or.r.eq.'g '.or.r.eq.'gi') then
*** aggiorna grafici
        s1p2=nint(par(27,2))
        icol=icol+1
        if(icol.gt.15) icol=1
        ymin=1.e6
        ymax=-1.e6
        chi2=.0
        fredeg=real(m-n)
        par(56,2)=fredeg
        do i=1,mwl
         wl=mis(7,i,1)
         par(7,1)=wl
         ev=12400./wl
         if(LEV.EQ.1) then
          wwl(i)=wl
         else
          wwl(i)=12400./wl
         end if
         if(r.eq.'f '.or.r.eq.'i '.or.r.eq.'g ') then
          call FDISP(ioptf,ev,PF,PM,sqn,sqk)
          cnk(1,2)=sqn
          cnk(1,3)=sqk
          if(r.eq.'i '.or.r.eq.'g ') then! ibridone, calcola k con i parametri di BF
           par(24,2)=real(i)
           mis(16,i,1)=sqn
           if(nint(par(32,2)).ne.2) then !k non imposto
*** calcolo di k da T
            k=DBLE(mis(16,i,2))
            cnk(1,3)=mis(16,i,2)
            vs(1)=DBLE(sqn)
            if(nint(par(50+nint(par(53,2)),3)).eq.1) then
             if(dato(1).gt.0) Trasm=mis(1,i,1)
             if(dato(2).gt.0) Trasm=mis(2,i,1)
             dinco=pm(nint(par(53,2)),1)
             if(iw.eq.3) write(*,*) 'Trasm = ',Trasm
             vs(21)=DBLE(-0.1*wl/4./3.14/dinco*log(Trasm))
             vs(22)=1.e-10   !dlim
            else
             vs(21)=0.001   !DELTA_k per film
             vs(22)=1.e-5   !dlim
            end if
            vs(23)=DBLE(par(25,1)) !tol
            iw=0 !verbosit
            call cercazero(DELTAT,iw,k,dT,vs)
            mis(16,i,2)=SNGL(k)
            cnk(1,3)=SNGL(k)
            kk(i)=SNGL(k)
            if(kk(i).ge.rxy(17,1)) then
             kklog(i)=log10(kk(i))
            else
             kklog(i)=rxy(17,1)
            end if
           else
            mis(16,i,2)=cnk(1,3) 
           end if
          end if
         else if(r.eq.'gi') then
          cnk(1,2)=mis(16,i,1)
          cnk(1,3)=mis(16,i,2)  
         end if
         nn(i)=cnk(1,2)
         kk(i)=cnk(1,3)
         e1(i)=cnk(1,2)*cnk(1,2)-cnk(1,3)*cnk(1,3)
         e2(i)=2.*cnk(1,2)*cnk(1,3)
         if(kk(i).ge.rxy(17,1)) then
          kklog(i)=log10(kk(i))
         else
          kklog(i)=rxy(17,1)
         end if
         ymax=MAX(kk(i),ymax)
         ymin=MIN(kk(i),ymin)
*** calcolo di R
         if(dato(3).eq.2.or.dato(5).eq.2) then
          te=0.
          CALL ASSEMBLER(i,wl,1,te,par,pm,vot)
          rn(i)=vot(2,s1p2)*100.
          r1(i)=vot(3,s1p2)*100.
          if(dato(3).eq.2) chi2=chi2+
     +       ((mis(3,i,1)-rn(i)/100.)/mis(3,i,2))**2./fredeg
          if(dato(5).eq.2) chi2=chi2+
     +       ((mis(5,i,1)-r1(i)/100.)/mis(5,i,2))**2./fredeg
         end if
         if(dato(4).eq.2) then
          te=par(6,1)/180.*3.141592654
          CALL ASSEMBLER(i,wl,1,te,par,pm,vot)
          rp(i)=vot(2,s1p2)*100.
          if(dato(4).eq.2) chi2=chi2+
     +       ((mis(4,i,1)-rp(i)/100.)/mis(4,i,2))**2./fredeg
         end if
        end do
        if(r.eq.'g '.or.r.eq.'gi '.or.(r.eq.'i '.and.n.eq.0)) then
         chi2fin=chi2
         par(55,2)=chi2
         par(27,1)=chi2
        end if
***grafico n
        CALL CXY(3,RXY,1,IXW,0)
        call pgsls(1)
        call pgsci(icol)
        call pgline(mwl,wwl,nn)
***grafico k
        iklog=0
        if(par(31,5).gt.0.5) then
         iklog=1
        end if
        CALL CXY(4,RXY,1,IXW,iklog)
        call pgsls(1)
        call pgsci(icol)
        if(iklog.eq.0) then
         call pgline(mwl,wwl,kk)
        else
         call pgline(mwl,wwl,kklog)
        end if
***grafico epsi1 e epsi2
        if(nint(par(10,1)).eq.1) then
         CALL CXY(8,RXY,1,IXW,0)
         call pgsls(1)
         call pgsci(icol)
         call pgline(mwl,wwl,e1)
         CALL CXY(9,RXY,1,IXW,0)
         call pgsls(1)
         call pgsci(icol)
         call pgline(mwl,wwl,e2)
        end if
***grafico delle R
        if(dato(3).eq.2) then
         CALL CXY(13,RXY,1,IXW,0)
         call pgsls(3)
         call pgsci(icol)
         call pgline(201,wwl,rn)
        end if
        if(dato(4).eq.2) then
         CALL CXY(14,RXY,1,IXW,0)
         call pgsls(3)
         call pgsci(icol)
         call pgline(201,wwl,rp)
        end if
        if(dato(5).eq.2) then
         CALL CXY(15,RXY,1,IXW,0)
         call pgsls(3)
         call pgsci(icol)
         call pgline(201,wwl,r1)
        end if
       end if
       
       if(ieon.eq.1) then
*** salvataggio valore centrale e deviazioni quadratiche
        if(jie.eq.0) then
         do i=1,201
          ncent(i)=nn(i)
          ndev(i)=0.
          kcent(i)=kk(i)
          kdev(i)=0.
         end do
         do i=1,n
          pst(i,1)=p(i) !valore centrale
          pst(i,2)=0.   !azzeramento sommatore
         end do
        else
         do i=1,201
          ndev(i)=(ncent(i)-nn(i))**2./jiemax
          kdev(i)=(kcent(i)-kk(i))**2./jiemax
         end do
         do i=1,n
          pst(i,2)=(pst(i,1)-p(i))**2./jiemax
         end do
        end if
       end if
       
              
       
       if(r.eq.'M '.and.ieon.eq.2) then
*** salvataggio interno di nk_ibridone
        ISV=nint(par(18,2))
        isol=nint(SOL(1,1))
        irima=998-isol
        write(*,*)'<<<<<<<<<<<< Salvataggio nk Ibridone >>>>>>>>>>>>>>'
        if(isv.eq.0) then
         sol(1,1)=mwl
         write(*,*)'Salvataggio soluzioni in temp:'
         WRITE(*,*)' N. soluzioni presenti = ',isol
         write(*,*)' N. nuove soluzioni =', mwl
         write(*,*)' spazio disponibile=', irima
        end if
        do i=1,mwl
         if(isv.eq.0) then
          if(LEV.EQ.1) then
           sol(i+1,1)=wwl(i)
          else
           sol(i+1,1)=12400./wwl(i)
          end if
          sol(i+1,2)=ncent(i)
          sol(i+1,3)=kcent(i)
          sol(i+1,4)=ndev(i)
          sol(i+1,5)=kdev(i)
          sol(i+1,6)=1.
         else 
          isv=isv+7
          mis(isv,i,1)=ncent(i)
          mis(isv,i,2)=kcent(i)
          nank(isv-7)='nk ibridone'
         end if
        end do
c        ieon=0
       end if
       
              
*** salva Job# in /temp
       if(ieon.eq.2.or.(ieon.eq.0.and.(r.eq.'f '.or.r.eq.'i '))) then
        jobtot=jobtot+1
        jobview=jobtot
        par(8,1)=jobtot
        par(8,2)=jobtot
        call PROGETTO(3)
       end if
       
       
       
       if(r.eq.'e ') then
*** edit abilitazione dati n
        call pgslct(ixw(3))
        call EDINK(SOL,LEV)
        mda=0
        par(38,1)=1.e+20
        par(38,2)=0.
        do i=2,md+1
         if(nint(sol(i,6)).gt.0) then
          mda=mda+1
          par(38,1)=MIN(par(38,1),sol(i,1))
          par(38,2)=MAX(par(38,2),sol(i,1))
         end if
        end do
        par(21,2)=mda
       end if
       
*** edit abilitazione dati R
       if(r.eq.'em') then
        call EDIMIS(lev,IXW,dato,MIS,mwla,par(38,3),par(38,4))
        par(21,1)=mwla
       end if
       
       if(r.eq.'eR') then
*** abilita tutta la R
        mwla=mwl
        do i=1,201
         mis(17,i,1)=1.
        end do
        par(21,1)=mwla
        if(dato(3).eq.2) then
         imR=3
        else if(dato(4).eq.2) then
         imR=4
        else if(dato(5).eq.2) then
         imR=5
        end if
        iwr=10+imR
        CALL PGSLS(1)
        call PLOTMIS(imR,iwr,3)
        par(38,3)=par(4,1)
        par(38,4)=par(4,2)
       end if
       
       if(r.eq.'eN') then
*** abilita tutti NK
        ial=0
        mda=0
        call pgslct(ixw(3))
        CALL pgsci(3)
        par(38,1)=1.e+20
        par(38,2)=0.
        do i=2,md+1
         if(sol(i,1).ge.par(4,1).and.sol(i,1).le.par(4,2)) then
          par(38,1)=MIN(par(38,1),sol(i,1))
          par(38,2)=MAX(par(38,2),sol(i,1))
          mda=mda+1
          sol(i,6)=1
          if(sol(i,4).le.0.) then
           if(ial.eq.0) then
            write(*,*)'ERRn = 0 @',sol(i,1),' => ERRn = Dn =',par(21,3)
            ial=1
           end if
           sol(i,4)=abs(par(21,3))
          end if
          XX=SOL(I,1)
          IF(LEV.EQ.2) XX=12400./XX
          CALL pgmove(XX,SOL(I,2)-SOL(I,4))
          IF(nint(SOL(I,6)).EQ.0) CALL pgsci(2)
          IF(nint(SOL(I,6)).EQ.1) CALL pgsci(4)
          CALL pgdraw(XX,SOL(I,2)+SOL(I,4))
         else
          sol(i,6)=0.
         end if
        end do
        par(21,2)=mda
       end if       

       
       if(r.eq.'R ') then 
*** riporta i parametri alla situazione migliore
        do i=1,200
         pm(i,1)=Pot(i)
        end do
       end if
       
       
       if(r.eq.'0 ') then
*** svuota il file SOL
        sol(1,1)=.0
        md=0
        mda=0
       end if    
          
*** navigazione Job#
c       if(r.eq.'+ '.or.r.eq.'- ') then
c        write(*,*) 'inizio: jobview=',jobview
c        if(r.eq.'+ ') jobview=jobview+1
c        if(r.eq.'- ') jobview=jobview-1
c        if(jobview.lt.1) jobview=1
c        if(jobview.gt.jobtot) jobview=jobtot
c        write(*,*) 'jobtot=',jobtot
c        write(*,*) 'jobview=',jobview
c        par(8,1)=jobtot
c        par(8,2)=jobview
c        call PROGETTO(4)
c**** re-impostazione pannello di fit e puntatore PM
c        npp=nint(par(34,5))
c        do i=1,npp
c         ppm(i)=nint(par(35+i,5))
c        end do
c**** ri-definizione pannello interfaccia fit      
c        call PANSTRING(SPM,SFMT)
c**** re-inizializzazione parametri di fit
c        n=nint(par(35,5))
c        do i=1,n
c         ipm=nint(pm(i,3))
c         p(i)=pm(ipm,1)
c        end do
c        do i=n+1,nparmax
c         p(i)=0.
c        end do
c**** ri-caricamento parametri utili da PAR
c        do i=1,14
c         dato(i)=nint(par(i,4))
c        end do
c        LEV=nint(RXY(25,3))
c        chi2fin=par(55,2)
c        fredeg=par(56,2)
c       end if
       
       
       write(*,*)
       if(ifit.eq.1) then ! Scrivi informazioni esito fit
        write(*,*)'********************INFO FIT***********************'
        write(*,*)'Punti abilitati  = ',m
        write(*,*)'Numero parametri = ',n,' => ',m-n,' gradi liberta"'
        if(n.gt.0) then
         write(*,*)'Chi^2: ',chi2ini,' --->',chi2fin
         par(27,1)=chi2ini
        else
         write(*,*)'Chi^2            = ',chi2fin
         par(27,1)=chi2fin
        end if
        write(*,*)'Final L2 norm of the residual =',fnorm
        if(info.lt.0) then
         write(*,*)'excetution terminated by user'
        elseif(info.eq.0) then
         write(*,*)'info = 0:  improper input parameters'
        elseif(info.eq.1) then
         write(*,*)'info = 1:  algorithm estimates that the relative err
     +or'
         write(*,*)'in the sum of squares is at most tol'
        elseif(info.eq.2) then
         write(*,*)'info = 2:  algorithm estimates that the relative err
     +or'
         write(*,*)'between x and the solution is at most tol'
        elseif(info.eq.3) then
         write(*,*)'info = 3:  conditions for info = 1 and info = 2 both
     + hold'
        elseif(info.eq.4) then
         write(*,*)'info = 4:  fvec is orthogonal to the columns of the
     +jacobian to machine precision'
        elseif(info.eq.5) then
         write(*,*)'info = 5:  number of calls to fcn has reached or exc
     +eeded 200*(n+1)'
        elseif(info.eq.6) then
         write(*,*)'info = 6: tol is too small. no further reduction in
     +the sum of squares is possible'
        elseif(info.eq.7) then
         write(*,*)'info = 7:  tol is too small. no further improvement
     +in the approximate solution x is possible'
        end if
        write(*,*)'*************************************************'
        write(*,*)
        ifit=0
        if(chi2fin.lt.chi2min) then !salva PM
         do i=1,200
          Pot(i)=PM(i,1)
         end do
         chi2min=chi2fin
        end if
       end if
       write(*,*)
     +  '--------- SEMAW - IBRIDONE - '//nank(16)(1:length)
     +//'.Spj ---------'
       write(*,*)'  N. wl abilitate  = ',mwla,' / ',mwl
       write(*,*)'  N. nk abilitati  = ',mda,' / ',md      
       write(*,*)'  Job# = ',jobview,'/',jobtot,'  Chi^2 = ',chi2fin
       write(*,*) st63
       par(21,1)=mwla
       par(21,2)=mda
       par(28,2)=mwl
       par(29,2)=md
       do j=1,npp
        ivp=nint(pm(ppm(j),2))
        if(ivp.eq.0) then
         if(SFMT(ppm(j)).eq.'110') then
          write(*,FMT=110)
     +    ' "',j,'"     ',spm(ppm(j)),pm(ppm(j),1)
         else if(SFMT(ppm(j)).eq.'112') then
          write(*,FMT=112)
     +    ' "',j,'"     ',spm(ppm(j)),pm(ppm(j),1)
         else if(SFMT(ppm(j)).eq.'114') then
          write(*,FMT=114)
     +    ' "',j,'"     ',spm(ppm(j)),pm(ppm(j),1)
         end if
        else
         if(SFMT(ppm(j)).eq.'110') then
          write(*,FMT=110)
     +    ' "',j,'"  X  ',spm(ppm(j)),pm(ppm(j),1),
     +                  ' +- ',pm(ppm(j),4),pm(ppm(j),5)
         else if(SFMT(ppm(j)).eq.'112') then
          write(*,FMT=112)
     +    ' "',j,'"  X  ',spm(ppm(j)),pm(ppm(j),1),
     +                  ' +- ',pm(ppm(j),4),pm(ppm(j),5)
         else if(SFMT(ppm(j)).eq.'114') then
          write(*,FMT=114)
     +    ' "',j,'"  X  ',spm(ppm(j)),pm(ppm(j),1),
     +                  ' +- ',pm(ppm(j),4),pm(ppm(j),5)
         end if
        end if
       end do
       nosc=nint(pf(ioptf,1))
       write(*,'(1x,a27,i1,1x,a12)')
     +  '"OF" opzione di lavoro FIT#',ioptf,'composta da:'
       do i=2,nosc+1
        iosc=nint(pf(ioptf,i))
        iotype=nint(pm(101+(iosc-1)*5,1))
        write(*,'(4x,a2,i2.2,a2,a24)')'"o',iosc,'" ',sfu(iotype)
       end do
          
c       if(ieon.eq.0.or.(ieon.eq.2.and.r.eq.'M ')) r='X '
       if(ieon.eq.0) then
        r='X '
       else if(ieon.eq.2.and.r.eq.'M ') then
        rc='gi'
        ieon=0
       end if
       
       par(8,1)=jobtot
       par(8,2)=jobview
       
      END DO
      if(iAlloIWA.eq.1) then
       DEALLOCATE(IWA,stat=ierror)
       if (ierror.ne.0) then
        print*,'error in deallocating IWA'
       end if
      endif
      if(iAlloP.eq.1) then
       DEALLOCATE(p,stat=ierror)
       if (ierror.ne.0) then
        print*,'error in deallocating p'
       end if
      end if
      if(iAlloFVEC.eq.1) then
       DEALLOCATE(FVEC,stat=ierror)
       if (ierror.ne.0) then
        print*,'error in deallocating FVEC'
       end if
      end if
      if(iAlloWA.eq.1) then
       DEALLOCATE(WA,stat=ierror)
       if (ierror.ne.0) then
        print*,'error in deallocating WA'
       end if
      end if
      if(iAlloFJAC.eq.1) then
       DEALLOCATE(FJAC,stat=ierror)
       if (ierror.ne.0) then
        print*,'error in deallocating FJAC'
       end if
      end if

 110  FORMAT(a2,i2,a6,a17,3x,f8.4,a4,f8.4,3x,f7.4,:)
 112  FORMAT(a2,i2,a6,a17,2x,e9.3,a4,f8.4,3x,f7.4,:)
 114  FORMAT(a2,i2,a6,a17,4x,f7.1,a4,1x,f7.1,3x,f7.4,:)
c 116  FORMAT(a2,i2,a6,a17,i2)
      
      RETURN
      END
      
      
      
      SUBROUTINE FSQ(M,N,P,FVEC,IFLAG)
c      PARAMETER(npmax=1001) ! Max numero dati sperimentali
c      PARAMETER(nparmax=17)! Max numero parametri di fit
c      PARAMETER(nv=16091) ! Dimensione dei vettori (lwa=m*n+5*n+m)
c      PARAMETER(nm=1001)  ! Dimensione delle matrici
      INTEGER M,N,IFLAG,IXW(20),iiflag
      REAL SOL(999,6),PF(7,21),PM(200,5),MIS(17,201,2),ELI(8,201,2),
     *     PAR(60,5),CNK(16,3),RXY(30,4),ARSE(500,2)
      REAL P(n),FVEC(m)
      CHARACTER NANK(16)*256
      LOGICAL logic
      COMMON /VALORI/ MIS,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      iiflag=iflag!per evitare warning unused
*** subroutine per la costruzione di fvec fittando n
      ioptf=nint(pm(100,1)) !FIT# di lavoro
      
*** indirizzamento parametri di fit su PM
      do i=1,n
       ip=nint(pm(i,3))
       pm(ip,1)=p(i)
      end do

*** calcolo fvec
      md=nINT(sol(1,1))
      j=0
      chi2=0.
      fredeg=real(m-n)
      erynMIN=(par(16,2)-par(16,1))/1.e5
      erykMIN=(par(17,2)-par(17,1))/1.e5
      do i=1,md
       if(nint(sol(i,6)).eq.1) then !il punto è abilitato
        wl=sol(i,1)
        ev=12400./wl
        call FDISP(ioptf,ev,PF,PM,sqn,sqk)
        if(nint(par(32,5)).le.1) then
         j=j+1
         y=sol(i,2)
         ery=sol(i,4)
         if(ery.lt.erynMIN) ery=erynMIN
         fvec(j)=(y-sqn)/ery
         chi2=chi2+(fvec(j)**2.)/fredeg
        end if
        if(nint(par(32,5)).eq.1) then
         j=j+1
         y=sol(i,3)
         ery=sol(i,5)
         if(ery.lt.erykMIN) ery=erykMIN
         fvec(j)=(y-sqk)/ery
         chi2=chi2+(fvec(j)**2.)/fredeg
        end if
        if(nint(par(32,5)).eq.2) then
         j=j+1
         fvec(j)=((sol(i,2)*sol(i,2)-sol(i,3)*sol(i,3))-
     +    (sqn*sqn-sqk*sqk))/0.01
         chi2=chi2+(fvec(j)**2.)/fredeg
         j=j+1
         fvec(j)=(2.*sol(i,2)*sol(i,3)-2.*sqn*sqk)/0.01
         chi2=chi2+(fvec(j)**2.)/fredeg
        end if
       end if
      end do
      INQUIRE(1,OPENED=logic)
      if(logic.eqv..TRUE.) write(1,*) chi2,(p(i),i=1,n)
      RETURN
      END
      


      
      
      SUBROUTINE FRCK(M,N,P,FVEC,IFLAG)
c      PARAMETER(npmax=1001) ! Max numero dati sperimentali
c      PARAMETER(nparmax=17)! Max numero parametri di fit
c      PARAMETER(nv=16091) ! Dimensione dei vettori (lwa=m*n+5*n+m)
c      PARAMETER(nm=1001)  ! Dimensione delle matrici
      INTEGER M,N,IFLAG,s1p2,dato(14),ixw(20)
      REAL SOL(999,6),PF(7,21),PM(200,5),MIS(17,201,2),ELI(8,201,2),
     *     PAR(60,5),CNK(16,3),RXY(30,4),VOT(5,2),ARSE(500,2)
      REAL P(n),FVEC(m)
      DOUBLE PRECISION vs(99),k,dT
      CHARACTER NANK(16)*256
      LOGICAL logic
      COMMON /VALORI/ MIS,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      EXTERNAL DELTAT
***** Subroutine per IBRIDONE: calcolo di k da T dato n
****                           poi fit di R      
***
**
      iflag=iflag !to avoid warning
      ioptf=nint(pm(100,1)) !FIT# di lavoro
      iw=0 !verbosit
      mwl=201
      s1p2=nint(par(27,2))
      do i=1,14
       dato(i)=nint(par(i,4))
      end do
      if(dato(4).eq.2) te=par(6,1)/180.*3.141592654
*** indirizzamento parametri di fit su PM
      do i=1,n
       ip=nint(pm(i,3))
       pm(ip,1)=p(i)
      end do
      if(iw.gt.0) write(*,*)'|FRCK> ',(p(i),i=1,n)
      chi2=.0
      fredeg=real(m-n)
*** calcolo fvec
      ivec=1
      do i=1,mwl
       if(nint(mis(17,i,1)).eq.1) then
        wl=mis(7,i,1)
        ev=12400./wl
        par(7,1)=wl
        par(24,2)=real(i)
        call FDISP(ioptf,ev,PF,PM,sqn,sqk)
        mis(16,i,1)=sqn
        cnk(1,2)=sqn
        if(nint(par(32,2)).eq.2) then !impongo k = FDISP
         mis(16,i,2)=sqk 
         cnk(1,3)=sqk
         k=DBLE(sqk)
        else !valore k iniziale per Ibridone
         cnk(1,3)=mis(16,i,2) 
         k=DBLE(mis(16,i,2))
        end if
        if(nint(par(32,2)).ne.2) then !k libero
*** calcolo di k da T
         vs(1)=DBLE(sqn)
**   DELTA_k
         if(nint(par(50+nint(par(53,2)),3)).eq.1) then !layer incoerente
          if(dato(1).gt.0) Trasm=mis(1,i,1)
          if(dato(2).gt.0) Trasm=mis(2,i,1)
          dinco=pm(nint(par(53,2)),1)
          if(iw.eq.3) write(*,*) 'Trasm = ',Trasm
          vs(21)=DBLE(-0.1*wl/4./3.14/dinco*log(Trasm))
          vs(22)=1.e-10   !dlim
         else
          vs(21)=0.001
          vs(22)=1.e-5   !dlim
         end if
         if(iw.eq.3) write(*,*)'DELTA_k = ',vs(21)
c        vs(23)=DBLE(par(25,1)) !tol
         vs(23)=vs(21)/1000. !tol
         call cercazero(DELTAT,iw,k,dT,vs)
         mis(16,i,2)=SNGL(k)
        end if
        cnk(1,3)=SNGL(k)   
*** calcolo di R
        if(dato(3).eq.2.or.dato(5).eq.2) then
         CALL ASSEMBLER(i,wl,1,0.,par,pm,vot)
         if(dato(3).eq.2) then
          fvec(ivec)=(mis(3,i,1)-vot(2,s1p2))/mis(3,i,2)
          chi2=chi2+(fvec(ivec)**2.)/fredeg
          ivec=ivec+1
         end if
         if(dato(5).eq.2) then
          fvec(ivec)=(mis(5,i,1)-vot(3,s1p2))/mis(5,i,2)
          chi2=chi2+(fvec(ivec)**2.)/fredeg
          ivec=ivec+1
         end if
        end if
        if(dato(4).eq.2) then
         CALL ASSEMBLER(i,wl,1,te,par,pm,vot)
         fvec(ivec)=(mis(4,i,1)-vot(2,s1p2))/mis(4,i,2)
         chi2=chi2+(fvec(ivec)**2.)/fredeg
         ivec=ivec+1
        end if
       end if
      end do
*** monitor valori
      write(*,*) chi2,(p(i),i=1,n)
      INQUIRE(1,OPENED=logic)
      if(logic.eqv..TRUE.) write(1,*) chi2,(p(i),i=1,n)
      RETURN
      END
      
      
      
      SUBROUTINE DELTAT(iw,k,vs,dT)
      INTEGER ixw(20),dato(14),s1p2
      REAL SOL(999,6),PF(7,21),PM(200,5),MIS(17,201,2),ELI(8,201,2),
     *     PAR(60,5),CNK(16,3),RXY(30,4),VOT(5,2),ARSE(500,2)
      DOUBLE PRECISION vs(99),k,dT
      CHARACTER NANK(16)*256
      COMMON /VALORI/ MIS,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
c subroutine per calcolare la differenza tra Texp e Tcalc

      do i=1,14
       dato(i)=nint(par(i,4))
      end do
      wl=par(7,1)
      iwl=nint(par(24,2))
      if(dato(1).eq.2) then
       texp=mis(1,iwl,1)
       te=0.
       s1p2=1
      else if(dato(2).eq.2) then
       texp=mis(2,iwl,1)
       te=par(6,1)/180.*3.141592654
       s1p2=nint(par(27,2))
      end if
      cnk(1,2)=SNGL(vs(1))
      cnk(1,3)=SNGL(k)
      CALL ASSEMBLER(iwl,wl,1,te,par,pm,vot)
      tcalc=vot(1,s1p2)
      dT=DBLE(tcalc-texp)
      if(iw.eq.3) then
       write(*,*)'|DELTAT> wl = ',wl,' iwl =',iwl,' Texp = ',texp
       write(*,*)'          n = ',cnk(1,2),' k = ',k
       write(*,*)'      Tcalc = ',tcalc,' dT = ',dT
      end if
      RETURN
      END
      


      SUBROUTINE EDINK(SOL,LEV)
      DIMENSION SOL(999,6)
      CHARACTER R*1
      R='A'
      DO WHILE(R.EQ.'A')
       write(*,*)'DELIMITAZIONE AREA DA DISABILITARE:'
       write(*,*)'right-button to ESC'
       WRITE(*,*)'1) posiziona il cursore su n_max wl_min ed INVIO'
       cursore=pgcurs(xmin,ymax,r)
       if(r.ne.'A') goto 1
       IF(LEV.EQ.2) xmin=12400./xmin
       WRITE(*,*)'wl_min = ',xmin,' n_max = ',ymax
       WRITE(*,*)'2) posiziona il cursore su n_min wl_max ed INVIO'
       cursore=pgcurs(xmax,ymin,r)
       IF(LEV.EQ.2) xmax=12400./xmax
       WRITE(*,*)'wl_max = ',xmax,' n_min = ',ymin
       DO I=2,nINT(SOL(1,1))+1
        IF(SOL(I,1).GE.xmin.AND.SOL(I,1).LE.xmax.and.
     +      SOL(I,2).GE.ymin.AND.SOL(I,2).LE.ymax) THEN
          SOL(I,6)=0.
          XX=SOL(I,1)
          IF(LEV.EQ.2) XX=12400./XX
          CALL pgmove(XX,SOL(I,2)-SOL(I,4))
          IF(nint(SOL(I,6)).EQ.0) CALL pgsci(2)
          IF(nint(SOL(I,6)).EQ.1) CALL pgsci(4)
          CALL pgdraw(XX,SOL(I,2)+SOL(I,4))
        END IF
       END DO
1     END DO
      RETURN
      END


      
!       SUBROUTINE CFREMA
!       REAL PF(7,21),M(17,201,2),E(8,201,2),PAR(60,5),CNK(16,3),
!      *     RXY(30,4),VNK(16,2),na,ka,nb,kb,ne,ke,v(201,3),SOL(999,6),
!      +     PM(200,5),ARSE(500,2)
!       INTEGER IXW(20)
!       CHARACTER R*1,NANK(16)*256,SOUT*19
!       COMMON /VALORI/ M,E,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
!       COMMON /STRINGHE/ NANK
!       fa=.5
!       icol=1
!       R='?'
!       DO WHILE(R.NE.'x')
!        write(*,*)'Calcolo n,k EMA con i materiali A e B:'
!        CALL GECNK(3,'s','"a":  ',sout)
!        CALL GECNK(4,'s','"b":  ',sout)
!        WRITE(*,*) '"f": f_a = ',1.-fb
!        write(*,*) '     f_b = ',fb
!        write(*,*) '"n": calcolo EMA e plot di n'
!        write(*,*) '"k": calcolo EMA e plot di k'
!        write(*,*) '"s": salva n,k_EMA sul vettore SOLnk'
!        write(*,*) '"x": torna a |FITnk>'
!        WRITE(*,*)'|CFRema> Che fare ?'
!        READ(*,'(a1)') R
! 
!        IF(R.EQ.'a') CALL GECNK(3,'y','a mate',sout)
! 
!        IF(R.EQ.'b') CALL GECNK(4,'y','b mate',sout)
! 
!        IF(R.eq.'f') THEN
! 1       write(*,*)'Nuovo valore di f_b?'
!         read(*,*,ERR=1) fb
!        END IF
!        
!        IF(R.eq.'n') THEN
!         R='c'
!         iw=3
!        END IF
! 
! 
!        IF(R.eq.'k') THEN
!         R='c'
!         iw=4
!        END IF
! 
!        IF(R.EQ.'c') THEN
!         ICOL=ICOL+1
!         IF(ICOL.GE.15) ICOL=1
!         CALL pgsci(ICOL)
!         CALL pgsls(1)
!         CALL CXY(iw,RXY,0,IXW)
!         CALL PGQWIN(XL,XH,YL,YH)
!         LEV=nINT(RXY(25,3))
!         j=0
!         do i=1,201
!          call cosvnk(1,vnk,i)
!          na=vnk(3,1)
!          ka=vnk(3,2)
!          nb=vnk(4,1)
!          kb=vnk(4,2)
!          call ema(na,ka,nb,kb,fb,ne,ke)
!          v(i,1)=m(7,i,1)
!          v(i,2)=ne
!          v(i,3)=ke
!          X=v(i,1)
!          IF(nINT(RXY(25,3)).EQ.2) X=12400./X      
!          if(nINT(RXY(25,2)).eq.16) then
!           y=v(i,2)
!          else
!           y=v(i,3)
!          end if
!          IF(X.GE.XL.AND.X.LE.XH.AND.Y.GE.YL.AND.Y.LE.YH) THEN
!           IF(J.EQ.0) CALL pgmove(X,Y)
!           IF(J.EQ.1) CALL pgdraw(X,Y)
!           J=1
!          ELSE
!           J=0
!          END IF
!         end do
!        END IF
! 
!        IF(R.EQ.'s') THEN
!         sol(1,1)=201.
!         do i=1,201
!          do jj=1,3
!           sol(1+i,jj)=v(i,jj)
!          end do
!          sol(1+i,4)=.0
!          sol(1+i,5)=.0
!          sol(1+i,6)=1.
!         end do
!        END IF
! 
!       END DO
!       RETURN
!       END



      SUBROUTINE CERCAZERO(FT,iw,t,d1,vs)
      IMPLICIT DOUBLE PRECISION (a - h,o - z)
      DOUBLE PRECISION vs(99)
      CHARACTER r*1
**********
c     iw: comando modalit di esecuzione del calcolo
c         iw=0 -> calcolo silente'
c         iw=1 -> calcolo verboso
c         iw=2 -> calcolo passo-passo
c         iw=3 -> calcolo passo anche nella subroutine FT
c
c     vs(99): vettore passaggio valori
c         1-20: valori di INPUT OUTPUT della subroutine FT
c        21-30: parametri ricerca ZERO
c           21=dt    incremento numerico per la delimitazione di ZERO
c           22=dlim  livello di soglia
c           23=tol   precisione di calcolo
c            ....
c           26=parametro di controllo fit: 
c                                      1 -> prosegui
c                                     -1 -> stop: iflag negativo in fcn
**********
      dt=vs(21)
      dlim=vs(22)
      tol=vs(23)
      iww=0
      if(iw.eq.3) iww=3
      call FT(iww,t,vs,d1)
      t0=t
      d0=d1
      if(iw.ge.2) then
       write(*,*)'C0: Valori iniziali:'
       write(*,*)'C0: t, d1 = ',t,d1
      end if
      istop=0
      irif=0
      inul=0
      imi=0
      ilim=1000
      if(abs(d1).le.dlim) istop=1
      do while(istop.eq.0)
       t=t+dt
       if(iw.eq.3) iww=3
       call FT(iww,t,vs,d1)
       if(iw.ge.2) then
        write(*,*)'C0: t0,d0 = ',t0,d0
        write(*,*)'C0: t,d1 = ',t,d1
       end if
       if(abs(d1).le.dlim) then
        istop=1
       else if(d0*d1.lt..0.and.abs(d1).ge.(10.*abs(d0))) then  
        t=t0
        d1=d0
        dt=dt/10.
        if(iw.ge.2) write(*,*)'|d1| > 10*|d0| -> dt decimato = ',dt
        if(abs(dt).le.dlim) then
         istop=1
         if(iw.ge.2) then
          write(*,*)'   dt = ',dt,' dlim = ',dlim
          write(*,*)'   STOP per raggiunto limite precisione'
         end if
        end if
       else if(d0*d1.lt..0.and.abs(d1).lt.(10.*abs(d0))) then ! Brent
        if(iw.ge.2) write(*,*)'Inizio ricerca con metodo di Brent'
        a=t
        b=t0
        c=b
        fa=d1
        fb=d0
        fc=d0
        itmax=1000
        eps=3.e-9
        iter=0
        do while(istop.eq.0)
         iter=iter+1
         if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a !Rename a, b, c and adjust bounding interval d.
          fc=fa
          d=b-a
          e=d
         endif
         if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
         endif
         if(iw.ge.2) then
          write(*,*)'Brent: N_iterazione= ',iter
          write(*,*)'   a,b,c=',a,b,c
          write(*,*)'   fa,fb,fc=',fa,fb,fc
          write(*,*)'Continuo?'
          read(*,'(a1)') r
         end if
         tol1=2.*EPS*abs(b)+0.5*tol ! Convergence check.
         xm=.5*(c-b)
         if(abs(xm).le.tol1.or.fb.eq.0..or.iter.ge.itmax)then
          if(iw.ge.2) then
           write(*,*)'eps = ',eps,' tol = ',tol,' tol1 = ',tol1
           write(*,*)'xm = ',xm
          end if
          t=b
          d1=fb
          istop=1              
          if(iw.ge.2) then
           if(iter.lt.itmax) then
            write(*,*)'STOP per raggiunto limite di precisione'
           else
            write(*,*)'STOP per raggiunto N.max iterazioni'
           end if
          end if  
         endif
         if(istop.eq.0) then
          if(abs(e).ge.tol1.and.abs(fa).gt.abs(fb)) then
           ss=fb/fa !Attempt inverse quadratic interpolation.
           if(a.eq.c) then
            pp=2.*xm*ss
            q=1.-ss
           else
            q=fa/fc
            rr=fb/fc
            pp=ss*(2.*xm*q*(q-rr)-(b-a)*(rr-1.))
            q=(q-1.)*(rr-1.)*(ss-1.)
           endif
           if(pp.gt.0.) q=-q !Check whether in bounds.
           pp=abs(pp)
           if(2.*pp.lt.min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d !Accept interpolation.
            d=pp/q
            if(iw.ge.2) write(*,*)'...interpolazione quadratica...'
           else
            if(iw.ge.2) write(*,*)'...bisezione..'
            d=xm !Interpolation failed, use bisection.
            e=d
           endif
          else !Bounds decreasing too slowly, use bisection.
           if(iw.ge.2) write(*,*)'...bisezione per velocizzare..'
           d=xm
           e=d
          endif
          a=b !Move last best guess to a.
          fa=fb
          if(abs(d).gt.tol1) then !Evaluate new trial root.
           b=b+d
          else
           b=b+sign(tol1,xm)
          endif
          t=sngl(b)
          if(iw.eq.3) iww=3
          call FT(iww,t,vs,fb)
          if(iw.ge.2) write(*,*)'b,fb = ',b,fb  
         end if           
        end do
       else if(abs(d1).lt.abs(d0)) then
        if(iw.ge.2) write(*,*)'miglioramento'   
        t0=t
        d0=d1
        irif=0
        imi=imi+1
       else if(abs(d1).gt.abs(d0).and.d0*d1.gt.0.) then
        if(iw.ge.2) write(*,*)'peggioramento'
        t=t0
        dt=-dt
        irif=irif+1
        imi=0
        if(irif.gt.2) then 
         dt=dt/2.
         if(iw.ge.2) write(*,*)'Incremento dimezzato!'
        end if
        d1=d0
       else       
        inul=inul+1
       end if
       if(irif.gt.ilim.or.imi.gt.ilim.or.inul.gt.ilim) then
        write(*,*)'C0: irif o imi o inul > 1000! -> stop'
c        write(*,*)'Prosegui (ret) verbosamente (w) o interrompi (i)?'
c        read(*,'(a1)')r
c        if(r.eq.'w') then 
c         iw=3
c        else if(r.eq.'i') then
         istop=1
         vs(26)=-1
c        else if(r.eq.' ') then
c         irif=0
c         imi=0
c         inul=0
c        end if
       end if
       if(iw.eq.3.and.istop.eq.0) then
        write(*,*)'C0: t, d1 = ',t, d1
        write(*,*)'Continuo?'
        read(*,'(a1)') r
       end if
      end do !istop.eq.0
      RETURN
      END
      
      
      
      SUBROUTINE PANSTRING(SPM,SFMT)
c subroutine per la costruzione del pannello interfaccia fit
      CHARACTER SPM(200)*17,SFMT(200)*3,sn*2
**** definizione di SPM(200)
      spm(1) ='d_1              '
      spm(2) ='d_2              '
      spm(3) ='d_3              '
      spm(4) ='d_4              '
      spm(5) ='d_5              '
      spm(6) ='d_6              '
      spm(7) ='d_7              '
      spm(8) ='d_8              '
      spm(9) ='d_9              '
      spm(10)='/                '
      spm(11)='Dn/<n>_1         '
      spm(12)='Dn/<n>_2         '
      spm(13)='Dn/<n>_3         '
      spm(14)='Dn/<n>_4         '
      spm(15)='Dn/<n>_5         '
      spm(16)='Dn/<n>_6         '
      spm(17)='Dn/<n>_7         '
      spm(18)='Dn/<n>_8         '
      spm(19)='Dn/<n>_9         '
      spm(20)='/                '
      spm(21)='(nav/<n>-1)_1    '
      spm(22)='(nav/<n>-1)_2    '
      spm(23)='(nav/<n>-1)_3    '
      spm(24)='(nav/<n>-1)_4    '
      spm(25)='(nav/<n>-1)_5    '
      spm(26)='(nav/<n>-1)_6    '
      spm(27)='(nav/<n>-1)_7    '
      spm(28)='(nav/<n>-1)_8    '
      spm(29)='(nav/<n>-1)_9    '
      spm(30)='/                '
      spm(31)='Dk/<k>_1         '
      spm(32)='Dk/<k>_2         '
      spm(33)='Dk/<k>_3         '
      spm(34)='Dk/<k>_4         '
      spm(35)='Dk/<k>_5         '
      spm(36)='Dk/<k>_6         '
      spm(37)='Dk/<k>_7         '
      spm(38)='Dk/<k>_8         '
      spm(39)='Dk/<k>_9         '
      spm(40)='/                '
      spm(41)='(kav/<k>-1)_1    '
      spm(42)='(kav/<k>-1)_2    '
      spm(43)='(kav/<k>-1)_3    '
      spm(44)='(kav/<k>-1)_4    '
      spm(45)='(kav/<k>-1)_5    '
      spm(46)='(kav/<k>-1)_6    '
      spm(47)='(kav/<k>-1)_7    '
      spm(48)='(kav/<k>-1)_8    '
      spm(49)='(kav/<k>-1)_9    '
      spm(50)='/                '
      spm(51)='sigma_roughness_1'
      spm(52)='sigma_roughness_2'
      spm(53)='sigma_roughness_3'
      spm(54)='sigma_roughness_4'
      spm(55)='sigma_roughness_5'
      spm(56)='sigma_roughness_6'
      spm(57)='sigma_roughness_7'
      spm(58)='sigma_roughness_8'
      spm(59)='sigma_roughness_9'
      spm(60)='                 '
      spm(61)='                 '
      spm(61)='slope_Dn/<n>_1   '
      spm(62)='slope_Dn/<n>_2   '
      spm(63)='slope_Dn/<n>_3   '
      spm(64)='slope_Dn/<n>_4   '
      spm(65)='slope_Dn/<n>_5   '
      spm(66)='slope_Dn/<n>_6   '
      spm(67)='slope_Dn/<n>_7   '
      spm(68)='slope_Dn/<n>_8   '
      spm(69)='slope_Dn/<n>_9   '
      spm(70)='                 '
      spm(71)='                 '
      spm(72)='                 '
      spm(73)='                 '
      spm(74)='                 '
      spm(75)='                 '
      spm(76)='                 '
      spm(77)='                 '
      spm(78)='                 '
      spm(79)='                 '
      spm(80)='                 '
      spm(81)='                 '
      spm(82)='                 '
      spm(83)='                 '
      spm(84)='                 '
      spm(85)='                 '
      spm(86)='                 '
      spm(87)='                 '
      spm(88)='                 '
      spm(89)='                 '
      spm(90)='                 '
      spm(91)='                 '
      spm(92)='                 '
      spm(93)='                 '
      spm(94)='                 '
      spm(95)='                 '
      spm(96)='                 '
      spm(97)='                 '
      spm(98)='                 '
      spm(99)='                 '
      spm(100)='                 '
      do i=1,20
       write(sn,'(i2.2)') i
       spm(101+(i-1)*5)='                 '     
       spm(102+(i-1)*5)='C'//sn//'              '
       spm(103+(i-1)*5)='E'//sn//' (eV)         '
       spm(104+(i-1)*5)='D'//sn//' (eV)         '
       spm(105+(i-1)*5)='K'//sn//'              '
      end do

      
**** impostazione indice del formato
      do i=1,10
       SFMT(i)='114'
      end do
      do i=11,50
       SFMT(i)='110'
      end do
      do i=51,60
       SFMT(i)='114'
      end do
      do i=1,20
       SFMT(102+(i-1)*5)='110'
       SFMT(103+(i-1)*5)='110'
       SFMT(104+(i-1)*5)='110'
       SFMT(105+(i-1)*5)='112'
      end do
      
      RETURN
      END
      
      
      
!       SUBROUTINE PANFIT(npp,SPM,PPM,PAR,PM)
!       PARAMETER(nparmax=17) ! Max numero di parametri del fit
!       INTEGER PPM(nparmax),IOS
!       REAL PAR(60,5),PM(200,5)
!       CHARACTER SPM(200)*17,r*1,st4*4
!       
! 
!       r='?'
!       do while(r.ne.' ')
!        write(*,*)'**************************************************'
!        write(*,*)'N.   parametro'
!        do i=1,npp
!         write(*,'(a1,i2,a1,4x,a17)') '"',i,'"',spm(ppm(i))
!        end do
!        write(*,*)'**************************************************'
!        write(*,*)'Che fare ("k" kill,"a" add, "m" modify, "s" swap)'
!        write(*,*)'         ("k#"    " aX#"    "m#"                )?'
!        read(*,'(a4)') st4
!        r=st4(1:1)
!        
!        if(r.eq.'k'.or.r.eq.'m') then
!         read(st4(2:4),*,IOSTAT=IOS) ip
!         if(ios.ne.0) then
!  4       write(*,*)'N. parametro?'
!          read(*,*,ERR=4) ip
!         end if
!         if(r.eq.'k') then
!          do i=ip,npp-1
!           ppm(i)=ppm(i+1)
!          end do
!          ppm(npp)=0
!          npp=npp-1
!         else if(r.eq.'m') then
!          write(*,*)'Codice -> Parametro                 valore attuale'
!          do i=1,200
!           if(spm(i).ne.'                ') 
!      +     write(*,'(4x,i3,a4,a17,3x,e10.4)')i,' -> ',spm(i),pm(i,1)
!          end do
!  1       write(*,*) 'Codice nuovo parametro?'
!          read(*,*,ERR=1) icod
!          ppm(ip)=icod
!         end if
!         
!        else if(r.eq.'a'.and.npp.lt.nparmax) then
!         read(st4(3:4),'(i2)',IOSTAT=IOS,ERR=5) iosc
!  5      if(ios.eq.0.and.iosc.ge.1) then
!          if(st4(2:2).eq.'C'.and.iosc.le.20) then
!           icod=101+(iosc-1)*5+1
!          else if(st4(2:2).eq.'E'.and.iosc.le.20) then
!           icod=101+(iosc-1)*5+2
!          else if(st4(2:2).eq.'D'.and.iosc.le.20) then
!           icod=101+(iosc-1)*5+3
!          else if(st4(2:2).eq.'K'.and.iosc.le.20) then
!           icod=101+(iosc-1)*5+4
!          else if(st4(2:2).eq.'d'.and.iosc.le.9) then
!           icod=iosc
!          else if(st4(2:2).eq.'g'.and.iosc.le.9) then
!           icod=10+iosc
!          else if(st4(2:2).eq.'u'.and.iosc.le.9) then
!           icod=20+iosc
!          else if(st4(2:2).eq.'G'.and.iosc.le.9) then
!           icod=30+iosc
!          else if(st4(2:2).eq.'U'.and.iosc.le.9) then
!           icod=40+iosc
!          else if(st4(2:2).eq.'r'.and.iosc.le.9) then
!           icod=50+iosc
!          end if
!         else
!          write(*,*)'Codice -> Parametro'
!          do i=1,200
!           if(spm(i).ne.'                ') 
!      +     write(*,'(4x,i3,a4,a17,3x,e10.4)')i,' -> ',spm(i),pm(i,1)
!          end do
!  2       write(*,*) 'Codice nuovo parametro?'
!          read(*,*,ERR=2) icod
!         end if
!          if(icod.ge.1.and.icod.le.200) then
!           npp=npp+1
!           ppm(npp)=icod
!         end if
!         
!        else if(r.eq.'a'.and.npp.eq.nparmax) then
!         write(*,*)'ATTENZIONE: raggiunto limite MAX N. parametri = ',
!      +            nparmax
!         
!        else if(r.eq.'s') then
!  3      write(*,*)'Indici parametri da scambiare (i1,i2)?'
!         read(*,*,ERR=3) i1,i2
!         if(i1.lt.1.or.i1.gt.npp.or.i2.lt.1.or.i2.gt.npp) goto 3
!         itmp=ppm(i1)
!         ppm(i1)=ppm(i2)
!         ppm(i2)=itmp
!        end if   
!       end do
!       
! *** salvataggio su PAR
!       par(34,5)=npp
!       do i=1,nparmax
!        par(35+i,5)=ppm(i)
!       end do
! 
!       RETURN
!       END
!       


      
      
      SUBROUTINE EDIMIS(lev,IXW,dato,MIS,mwla,WLmin,WLmax)
      CHARACTER R*1
      INTEGER LEV,IXW(20),dato(14)
      REAL MIS(17,201,2),WLmin,WLmax
*** subroutine per abilitazione/disabilitazione misura per fit ibridone
      if(dato(3).eq.2) then
       imR=3
      else if(dato(4).eq.2) then
       imR=4
      else if(dato(5).eq.2) then
       imR=5
      end if
      iwr=10+imR
      call pgslct(IXW(iwr))
      write(*,*)'Finestra attiva = ',IXW(iwr)
      write(*,*)'Selezione intervallo da abilitare/disabilitare'
      WRITE(*,*)'POSIZIONA IL CURS SU ESTREMO SINISTRO E invio !'
      cursore=pgcurs(xval,yval,r)
      RGLMIN=XVAL
      IF(LEV.EQ.2) THEN
       RGLMIN=12400./RGLMIN
       WRITE(*,*)'LAMBDA MIN = ',RGLMIN
      END IF
      WRITE(*,*)'POSIZIONA IL CURS SU ESTREMO DESTRO E invio !'
      cursore=pgcurs(xval,yval,r)
      RGLMAX=XVAL
      IF(LEV.EQ.2) THEN
       RGLMAX=12400./RGLMAX
       WRITE(*,*)'LAMBDA MAX = ',RGLMAX
      END IF
      mwla=0
      WLmin=1.e+20
      WLmax=0.
      do i=1,201
       if(mis(7,i,1).ge.rglmin.and.mis(7,i,1).le.rglmax) 
     +    mis(17,i,1)=1.-mis(17,i,1)
       if(nint(mis(17,i,1)).eq.1) then
        mwla=mwla+1
        WLmin=MIN(WLmin,mis(7,i,1))
        WLmax=MAX(WLmax,mis(7,i,1))
       end if
      end do
      CALL PGSLS(1)
      call PLOTMIS(imR,iwr,3)
      RETURN
      END

