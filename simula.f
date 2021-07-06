*********    kSEMAW SIMULA sources   ***************
*
*   hereinafter subroutines to simulate the experimental measurements
*   on the basis of nk-files and oscillator analytical functions
*
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


      SUBROUTINE CNFR(r)
c      DOUBLE PRECISION mrtp(99,99),mblade(99,99),ri,re,alpha,areat,
c     +  areac,alphai,fraz
      REAL M(17,201,2),E(8,201,2),PAR(60,5),PM(200,5),CNK(16,3),
     +     SOL(999,6),PF(7,21),RXY(30,4),MC(14,201),wwl(201),y(201),
     +     vot(5,2),vt(3,2),av(6),avs(6),tet(201),tau(201),rho(201),
     +     rho1(201),wp(2501,2),val(6,2),ARSE(500,2)
      INTEGER DATI(14),S1P2,IXW(20),S1P2U3
      CHARACTER NANK(16)*256,SMIS(19)*24,SOUT*19,R*2,
     +     st76*76,tab*1,S2*2
      PARAMETER(d2r=0.017453293)
c      PARAMETER(pig=3.14159265359)
      COMMON /VALORI/ M,E,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      
*** inizializzazione stringhe
      SMIS(1)='Tn'
      SMIS(2)='Tp'
      SMIS(3)='Rn'
      SMIS(4)='Rp'
      SMIS(5)='R1'
      SMIS(6)='Apds'
      SMIS(7)='DELTA_1'
      SMIS(8)='PSI_1'
      SMIS(9)='DELTA_2'
      SMIS(10)='PSI_2'
      SMIS(11)='DELTA_3'
      SMIS(12)='PSI_3'
      SMIS(13)='DELTA_4'
      SMIS(14)='PSI_4'
      SMIS(15)='           '
      SMIS(16)='           '
      SMIS(17)='           '
      SMIS(18)='           '
      SMIS(19)='           '
      
c       typm(1)='Notch                               '
c       ity(1)=5
c       st12(1)='            '
c       
c       typm(2)='VIS 10 wl                           '
c       ity(2)=9
c       st12(2)='ufo_/media.1'
c       
c       typm(3)='VIS 30 wl                           '
c       ity(3)=9
c       st12(3)='ufo_/media.2'
c 
c       typm(4)='VIS D65                             '
c       ity(4)=7
c       st12(4)='ufo_/media.3'
c       
c       typm(5)='solar 10 wl                         '
c       ity(5)=11
c       st12(5)='ufo_/media.4'
c       
c       typm(6)='solar global radiation              '
c       ity(6)=22
c       st12(6)='ufo_/media.5'
c       
c       typm(7)='solar direct normal radiation 100 wl'
c       ity(7)=36
c       st12(7)='ufo_/media.6'
c       
c       typm(8)='PWO scintillation                   '
c       ity(8)=17
c       st12(8)='ufo_/pwtdr.1'
c       
c       typm(9)='PWO_scintillation * APD_iqe         '
c       ity(9)=27
c       st12(9)='ufo_/pwiqe.1'
c 
c       typm(10)='direct-circumsolar radiation 441 wl'
c       ity(10)=36
c       st12(10)='ufo_/media.7'
c       
c       typm(11)='direct-circumsolar rad ASTM G173-03'
c       ity(11)=36
c       st12(11)='ufo_/media.9'

      tab=char(9)
      
      length=len_trim(nank(16))
      length2=len_trim(nank(11))

**** caricamento parametri utili da PAR
      do i=1,14
       dati(i)=nint(par(i,4))
      end do
      LEV=nint(RXY(25,3))
      s1p2=nint(par(27,2))   !polarizzazione Tp e/o Rp
      imed=nint(par(35,1))   !tipo di media
      s1p2u3=nint(par(35,2)) !polarizzazione per la media
      wlnotch=par(35,3)      !WL nella media notch 
      cnk(1,1)=par(35,4)     !opzione nk incogniti
      if(wlnotch.lt.m(7,1,1).or.wlnotch.gt.m(7,201,1)) then
       wlnotch=m(7,1,1)
       par(35,3)=wlnotch
      end if
      
***azzeramento mc
      do i=1,201
       do j=1,14
        mc(j,i)=0.
       end do
      end do

***** loop azione
      icol=nint(par(7,2))
       write(*,*)
     +'------------- SEMAW - SIMULA - '//nank(16)(1:length)
     +//'.Spj --------------'
       WRITE(*,*)'"o " opzione nk-incogniti:'
       CALL GECNK(1,'s','      ',sout)
       write(*,*)'"me" media = ',imed,' S1P2U3 = ',s1p2u3
       write(*,*) nank(10)
c       if(imed.eq.1) 
c     +  write(*,*)'  wl_notch = ',wlnotch,'  S1P2U3 = ',s1p2u3
       write(*,*)
     +  '--------------------------------------------------------------'
       

!         write(*,*)'"g " grafica curve simulate'
!         write(*,*)'"gt" grafica VS teta la media dei valori SF simulati'
!         write(*,*)'"G " refresh grafici mantenendo i range'
!         write(*,*)'"pa" plot Absorptance at each layer
!         write(*,*)'"ms" matrice rho_near_specular thetaXphi'
       write(*,*)'comando (r)= ',r
       

**** refresh grafici misure EXP
       if(r.eq.'G ') then
!        CALL PGSLS(1)
        do i=1,14
         if(dati(i).ne.0) then
          imR=i
          iwr=10+i
          CALL PLOTMIS(imR,iwr,3) !plot mis
         end if
        end do
        CALL PLOTNK(1,3)  !plot n
        CALL PLOTNK(2,3)  !plot k
        if(nint(par(10,1)).eq.1) then
         CALL PLOTNK(3,3) !plot epsi1
         CALL PLOTNK(4,3) !plot epsi2
        end if
*** salvataggio su HD misure sperimentali splined
        open(10,status='unknown',file='expo/MisSFexp.dat')
        open(11,status='unknown',file='expo/MisSFexpErr.dat')
        open(12,status='unknown',file='expo/MisELIexp.dat')
        open(13,status='unknown',file='expo/MisELIexpErr.dat')
        write(10,*)'Tn Tp Rn Rp R1 APDS wl(A) n1 n2 n3 n4 n5 n6 n7 n8 n_
     +ibri'
        write(11,*)'eTn eTp eRn eRp eR1 eAPDS eLETT k1 k2 k3 k4 k5 k6 k7
     + k8 k_ibri'
        write(12,*)'DEL_1 PSI_1 DEL_2 PSI_2 DEL_3 PSI_3 DEL_4 PSI_4'
        write(13,*)'eDEL_1 ePSI_1 eDEL_2 ePSI_2 eDEL_3 ePSI_3 eDEL_4 ePS
     +I_4'
        do i=1,201
         write(10,FMT=100) (M(j,i,1),tab,j=1,16)
         write(11,FMT=101) (M(j,i,2),tab,j=1,16)
         write(12,FMT=102) (E(j,i,1),tab,j=1,8)
         write(13,FMT=102) (E(j,i,2),tab,j=1,8)
        end do
        close(10)
        close(11)
        close(12)
        close(13)        
       end if

       
**** calcolo e grafico curve calcolate
       if(r.eq.'g ') then
        irfh=1
        icol=icol+1
        if(icol.gt.15) icol=1
        par(7,2)=icol
        call CALMIS(MC)
        do i=1,201
         if(LEV.EQ.1) then
          wwl(i)=m(7,i,1)
         else
          wwl(i)=12400./m(7,i,1)
         end if
        end do
        if(nint(par(19,1)).eq.1) then !misure_calcolate -> misure_sperimentali
         i1=1
         i2=14
         do ic=i1,i2,1
          if(dati(ic).ne.0) then
           dati(ic)=1
           par(ic,4)=1.
           do i=1,201
            if(ic.le.5) then
             m(ic,i,1)=mc(ic,i)
            else if(ic.eq.6) then
             m(ic,i,1)=mc(ic,i)
            else if(ic.ge.7.and.ic.le.14) then
             e(ic-6,i,1)=mc(ic,1)
            end if
           end do
          end if
         end do
        end if
        write(*,*)'<<<<<<<<<<<<<<<< confronto EXP-SIM >>>>>>>>>>>>>>>>>'
        write(*,'(3x,a6,20x,a14)')'Misura','Scarto/Err RMS'
        do ic=1,14
         if(dati(ic).ne.0) then
          write(*,*)'ic=',ic,' dati(ic)=',dati(ic)
          if(ic.ge.1.and.ic.le.6) iw=ic+10
          if(ic.eq.7.or.ic.eq.9.or.ic.eq.11.or.ic.eq.13) iw=17
          if(ic.eq.8.or.ic.eq.10.or.ic.eq.12.or.ic.eq.14) iw=18
          if(ic.ge.7.and.ic.le.14) rxy(25,2)=real(ic)
          CALL CXY(iw,RXY,irfh,IXW,0)
          call pgsls(3)
          call pgsci(icol)
          FM=.0
          do i=1,201
           y(i)=mc(ic,i)
           if(ic.le.5) then
            y(i)=y(i)*100.
           end if
           if(ic.le.6) then
            FM=FM+((mc(ic,i)-m(ic,i,1))/m(ic,i,2))**2./201.
           else
            FM=FM+((mc(ic,i)-e(ic-6,i,1))/e(ic-6,i,2))**2./201.
           end if
          end do
          write(*,'(3x,a24,3x,F10.3)') smis(ic),FM
          par(35+ic,1)=FM
          call pgline(201,wwl,y)          
         end if
        end do        
        write(*,*)'<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>'
*** salvataggio su HD curve simulate      
        OPEN(10,STATUS='UNKNOWN',FILE='expo/MisSim.dat',ERR=5)
        write(10,*) 'wl',tab,'Tn',tab,'Tp',tab,'Rn',tab,
     +       'Rp',tab,'R1',tab,'Apds ',tab,'DELTA',tab,'PSI  '
        do i=1,201
         write(10,FMT=14) m(7,i,1),(tab,mc(ii,i),ii=1,14)
        end do
        CLOSE(10)
5     end if

**** calcolo valor medio SF sperimentali e simulate
       if(r.eq.'g '.or.r.eq.'G ') then
c   setta la media
        open(10,status='old',file=nank(10),iostat=ios,err=13)
        read(10,'(a76)') st76
        write(*,'(a76)') st76
        read(10,'(i5)') ipr
        write(*,*) 'N dati = ',ipr
        do i=1,ipr
         read(10,*) wp(i,1),wp(i,2)
        end do
        close(10)
 13     if(ios.ne.0) write(*,*)'ATTENZIONE file_media non trovato!!!'
c   peso media
        speso=0.
        jjmin=1
        do while(wp(jjmin,1).lt.m(7,1,1).and.jjmin.lt.ipr)
         jjmin=jjmin+1
        end do
        jj=jjmin
        jjmax=jj
        do while(wp(jj,1).ge.m(7,1,1).and.wp(jj,1).le.m(7,201,1).and.
     +          jj.le.ipr)
         speso=speso+wp(jj,2)
         jjmax=jj
         jj=jj+1
        end do
        write(*,*)'jjmin=',jjmin,'wl=',wp(jjmin,1)
        write(*,*)'jjmax=',jjmax,'wl=',wp(jjmax,1)
        write(*,*)'SUMpeso=',speso
c   azzero i sommatori
        do i=1,6
         av(i)=.0
         avs(i)=.0
        end do
c   calcolo media a teta
        do ii=jjmin,jjmax
         ij=1
         do while(wp(ii,1).gt.m(7,ij,1).and.ij.le.200)
          ij=ij+1
         end do
         if(ij.gt.1) ij=ij-1
         x=(wp(ii,1)-m(7,ij,1))/(m(7,ij+1,1)-m(7,ij,1))
         do j=1,2
          ik=ij+j-1
          if(j.eq.1) then
            do iv1=1,6
             val(iv1,1)=m(iv1,ik,1)*(1.-x)
             val(iv1,2)=mc(iv1,ik)*(1.-x)
            end do
          else
            do iv1=1,6
             val(iv1,1)=val(iv1,1)+m(iv1,ik,1)*x
             val(iv1,2)=val(iv1,2)+mc(iv1,ik)*x
            end do
          end if
         end do !j
         do j=1,6
          av(j)=av(j)+val(j,1)*wp(ii,2)/speso
          avs(j)=avs(j)+val(j,2)*wp(ii,2)/speso
         end do
        end do !ii
        write(*,*)'<<<<<<<<<<<<<<<< Valori medi (%) >>>>>>>>>>>>>>>>>>>'
        write(*,'(23x,a12,3x,a10)')'Experimental','Simulation'
        do i=1,6
         if(dati(i).ne.0) then
          write(*,'(17x,a4,2x,f9.2,8x,f9.2)') smis(i)(1:4),av(i)*100.,
     +                                avs(i)*100.
          par(35+i,2)=av(i)*100.
          par(35+i,3)=avs(i)*100.
         end if
        end do
        write(*,*)'<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>'
       end if

**** calcolo media misure SF VS teta
       if(r.eq.'gt') then
c   setta la media
        open(10,status='old',file=nank(10),iostat=ios,err=11)
        read(10,'(a76)') st76
        write(*,'(a76)') st76
        read(10,'(i5)') ipr
        write(*,*) 'N dati = ',ipr
        do i=1,ipr
         read(10,*) wp(i,1),wp(i,2)
        end do
        close(10)
 11     if(ios.ne.0) write(*,*)'ATTENZIONE file_media non trovato!!!'
        
c   peso media
        speso=0.
        do jj=1,ipr
         speso=speso+wp(jj,2)
        end do
c   loop su teta
        NT=91
        DO I=1,NT
         teta=rxy(21,1)+(rxy(21,2)-rxy(21,1))*(i-1)/real(nt-1)
         if(teta.gt.90.) teta=90.
         tetar=teta*d2r
         av(1)=.0
         av(2)=.0
         av(3)=.0
         ymin=1.e+36
         ymax=-1.e+36
c   calcolo media a teta
         do ii=1,ipr
          if(wp(ii,1).le.m(7,1,1)) then
           ij=1
           x=0.
           if(wp(ii,1).lt.m(7,1,1).and.i.eq.1)
     *      write(*,*)'Attenzione wl_media = ',wp(ii,1),'< WLmin !!!'
          elseif(wp(ii,1).gt.m(7,1,1).and.wp(ii,1).lt.m(7,201,1))
     *    then
           ij=1
           do while(wp(ii,1).gt.m(7,ij,1))
            ij=ij+1
           end do
           if(ij.gt.1) ij=ij-1
           x=(wp(ii,1)-m(7,ij,1))/(m(7,ij+1,1)-m(7,ij,1))
          elseif(wp(ii,1).ge.m(7,201,1)) then
           ij=200
           x=1.
           if(wp(ii,1).gt.m(7,201,1).and.i.eq.1)
     *      write(*,*)'Attenzione wl_media = ',wp(ii,1),'> WLmax !!!'
          end if
          do j=1,2
           ik=ij+j-1
           par(7,1)=M(7,ik,1)
           if(teta.lt.90.) then
            call ASSEMBLER(ik,M(7,ik,1),1,tetar,par,pm,vot)
           else
            vot(1,1)=.0
            vot(1,2)=.0
            vot(2,1)=1.
            vot(2,2)=1.
            vot(3,1)=1.
            vot(3,2)=1.
           end if
           if(j.eq.1) then
            do iv2=1,2
             do iv1=1,3
              vt(iv1,iv2)=vot(iv1,iv2)*(1.-x)
             end do
            end do
           else
            do iv2=1,2
             do iv1=1,3
              vt(iv1,iv2)=vt(iv1,iv2)+vot(iv1,iv2)*x
             end do
            end do
           end if
          end do !j
          do j=1,3
           if(s1p2u3.le.2) then
            av(j)=av(j)+vt(j,s1p2u3)*wp(ii,2)/speso
           else
            av(j)=av(j)+(vt(j,1)+vt(j,2))/2.*wp(ii,2)/speso
           end if
          end do
         end do !ii
         do j=1,3
          ymin=MIN(ymin,av(j))
          ymax=MAX(ymax,av(j))
         end do
         tet(i)=teta
         tau(i)=av(1)*100.
         rho(i)=av(2)*100.
         rho1(i)=av(3)*100.
        end do !i
        rxy(23,3)=ymin*100.
        rxy(23,4)=ymax*100.
        CALL PGSLS(1)
        call cxy(20,rxy,3,IXW,0)
        CALL pgsci(icol)
        CALL PGLINE(NT,tet,tau)
        icol=icol+1
        IF(ICOL.GE.15) ICOL=2
        CALL pgsci(icol)
        CALL PGLINE(NT,tet,rho)
        icol=icol+1
        IF(ICOL.GE.15) ICOL=2
        CALL pgsci(icol)
        CALL PGLINE(NT,tet,rho1)
        icol=icol+1
        IF(ICOL.GE.15) ICOL=2
        par(7,2)=icol
*** salvataggio su file curve medie VS teta
        OPEN(10,STATUS='UNKNOWN',FILE='expo/ThetaRhoVStheta.dat',
     +    IOSTAT=IOS,ERR=12)
        write(10,*)'% Teta(deg)  Tau  Rho  Rho1'
        write(10,*) '%',NT
        do i=1,NT
         write(10,*) tet(i),tab,tau(i),tab,rho(i),tab,rho1(i)
        end do
        CLOSE(10)
        write(*,*)'data vs theta saved in: expo/ThetaRhoVStheta.dat'
 12    end if
              
       if(r.eq.'pA') then
        write(*,*)'plot Abs@0deg at each layer of ',nint(par(51,2))
        icol=1
        CALL PGSLS(1)
        call cxy(19,rxy,3,IXW,0)
        iLayer0=nint(par(51,2))
        do iLayer=1,iLayer0
         par(51,2)=iLayer
         if(iLayer.lt.iLayer0) then
          cnk(16,1)=cnk(nint(par(50+iLayer+1,1)),1)
          cnk(16,2)=cnk(nint(par(50+iLayer+1,1)),2)
          cnk(16,3)=cnk(nint(par(50+iLayer+1,1)),3)
         else
          cnk(16,1)=cnk(9,1) ! mezzo_OUT per SF
          cnk(16,2)=cnk(9,2)
          cnk(16,3)=cnk(9,3)
         end if
** salvataggio Abs su file distinti
         write(s2,'(I2.2)') iLayer
         OPEN(10,STATUS='UNKNOWN',FILE='expo/Abs#'//S2//'.dat',
     +    IOSTAT=IOS,ERR=23)
         write(10,*)'Lambda Absorptance Layer#'//S2
         write(10,*) 201
         do i=1,201
          if(LEV.EQ.1) then
           wwl(i)=m(7,i,1)
          else
           wwl(i)=12400./m(7,i,1)
          end if
          par(7,1)=m(7,i,1)
          CALL ASSEMBLER(i,m(7,i,1),1,0.,par,pm,vot)
          y(i)=1.-vot(1,1)-vot(2,1)
          write(10,*) m(7,i,1),y(i)
         end do
         close(10)
         CALL pgsci(icol)
         CALL PGLINE(201,wwl,y)
         icol=icol+1
         if(icol.gt.15) icol=1
        end do
 23    end if

c*** matrice rho_near_specular thetaXphi
c       if(r.eq.'ms') then     
c  setta la media
c        if(imed.eq.1) then
c         ipr=1
c         wp(1,1)=wlnotch
c         wp(1,2)=1.
c         ios=0
c        else
c         open(10,status='old',file=st12(imed),iostat=ios,err=17)
c         read(10,'(a76)') st76
c         write(*,'(a76)') st76
c         read(10,'(i5)') ipr
c         write(*,*) 'N dati = ',ipr
c         do i=1,ipr
c          read(10,*) wp(i,1),wp(i,2)
c         end do
c         lose(10)
c 17      if(ios.ne.0) write(*,*)'ATTENZIONE file_media non trovato!!!'
c        end if
c  peso media
c        speso=0.
c        do jj=1,ipr
c         speso=speso+wp(jj,2)
c        end do
c        write(*,*)
c        write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c        write(*,*)'Calcolo matrice rho_near_specular teta X phi:'
c        write(*,*)'R~s(phi)=Rh*exp(-(4*pi*sigma*cos(teta))^2/wl^2)'
c        write(*,*)'sigma(phi)=A1*exp(-(phi/s1)^2)+A2*exp(-(phi/s2)^2)'
c        write(*,*)'"1" A1 (nm)  = ',par(10,1)
c        write(*,*)'"2" s1 (mrad)= ',par(10,2)
c        write(*,*)'"3" A2 (nm)  = ',par(11,1)
c        write(*,*)'"4" s2 (mrad)= ',par(11,2)
c        write(*,*)'"5" phiLimit(mrad) = ',par(9,2)
c        write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c        phiL=nint(par(9,2))
c        open(10,status='UNKNOWN',file=nank(16)(1:length)//'.mrtp')
c        write(10,*)'Matrice rho_near_specular; theta(deg) X phi(mrad)'
c        write(10,*)'sample= ',nank(11)(1:length2)
c  loop su phi
c        do ji=0,phiL
c         phi=ji
c         if(ji.eq.0) then
c          sigma=0.
c         else
c          sigma=(par(10,1)*exp(-(phi/par(10,2))**2.)+
c     +           par(11,1)*exp(-(phi/par(11,2))**2.))*10.
c         end if
c         write(*,*)'ji=',ji,'phi=',phi,'sigma=',sigma
c  loop su teta
c         NT=91
c         DO I=1,NT
c          teta=rxy(21,1)+(rxy(21,2)-rxy(21,1))*(i-1)/real(nt-1)
c          if(teta.gt.90.) teta=90.
c          tetar=teta*d2r
c          av(1)=.0
c          av(2)=.0
c          av(3)=.0
c          ymin=1.e+36
c          ymax=-1.e+36
c  calcolo media a teta
c          do ii=1,ipr
c           if(wp(ii,1).le.m(7,1,1)) then
c            ij=1
c            x=0.
c            if(wp(ii,1).lt.m(7,1,1).and.i.eq.1)
c     *    write(*,*)'Attenzione wl_media = ',wp(ii,1),'< WLmin !!!'
c           elseif(wp(ii,1).gt.m(7,1,1).and.wp(ii,1).lt.m(7,201,1))
c     *     then
c            ij=1
c            do while(wp(ii,1).gt.m(7,ij,1))
c             ij=ij+1
c            end do
c            if(ij.gt.1) ij=ij-1
c            x=(wp(ii,1)-m(7,ij,1))/(m(7,ij+1,1)-m(7,ij,1))
c           elseif(wp(ii,1).ge.m(7,201,1)) then
c            ij=200
c            x=1.
c            if(wp(ii,1).gt.m(7,201,1).and.i.eq.1)
c     *       write(*,*)'Attenzione wl_media = ',wp(ii,1),'> WLmax !!!'
c           end if
c           do j=1,2
c            ik=ij+j-1
c            par(7,1)=M(7,ik,1)
c            if(teta.lt.90.) then
c             all ASSEMBLER(ik,1,tetar,par,pm,vot)
c             vot(2,1)=vot(2,1)*
c     +exp(-(4.*pig*sigma*cos(tetar))**2./par(7,1)**2.)
c             vot(2,2)=vot(2,2)*
c     +exp(-(4.*pig*sigma*cos(tetar))**2./par(7,1)**2.)
c            else
c             vot(1,1)=.0
c             vot(1,2)=.0
c             vot(2,1)=1.
c             vot(2,2)=1.
c             vot(3,1)=1.
c             vot(3,2)=1.
c            end if
c            if(j.eq.1) then
c             do iv2=1,2
c              do iv1=1,3
c               vt(iv1,iv2)=vot(iv1,iv2)*(1.-x)
c              end do
c             end do
c            else
c             do iv2=1,2
c              do iv1=1,3
c               vt(iv1,iv2)=vt(iv1,iv2)+vot(iv1,iv2)*x
c              end do
c             end do
c            end if
c           end do !j
c           do j=1,3
c            if(s1p2u3.le.2) then
c             av(j)=av(j)+vt(j,s1p2u3)*wp(ii,2)/speso
c            else
c             av(j)=av(j)+(vt(j,1)+vt(j,2))/2.*wp(ii,2)/speso
c            end if
c           end do
c          end do !ii
c          do j=1,3
c           ymin=MIN(ymin,av(j))
c           ymax=MAX(ymax,av(j))
c          end do
c          tet(i)=teta
c          tau(i)=av(1)*100.
c          rho(i)=av(2)*100.
c          rho1(i)=av(3)*100.
c          mrtp(ji+1,i)=DBLE(rho(i)/100.)
c         end do !i
c         if(ji.eq.0) 
c     +  write(10,'(a3,91(1a,f9.4))')"phi",(tab,tet(i),i=1,91)
c       write(10,'(92(f9.4,1a))')phi,tab,(mrtp(ji+1,i),tab,i=1,91)
c         rxy(23,3)=ymin*100.
c         rxy(23,4)=ymax*100.
c         CALL PGSLS(1)
c         all cxy(20,rxy,1,IXW,0)
c         CALL pgsci(3)
c         CALL PGLINE(NT,tet,rho)
c        end do  !ji
c        lose(10)
c*** calcolo matrice mblade
c        do i=1,NT
c         mblade(1,i)=DBLE(0.5)*mrtp(phiL+1,i) 
c         do ji=phiL+1,3,-1
c          mrtp(ji,i)=mrtp(ji,i)-mrtp(ji-1,i) !converto in incrementi
c         end do
c         if(i.eq.1) then
c          do ji=1,phiL+1
c           write(*,*) ji,mrtp(ji,i)
c          end do
c         end if
c         do ji=1,phiL-1
c          mblade(ji+1,i)=DBLE(0.)
c          rp=REAL(ji)
c          do jj=ji,phiL-1
c           ri=REAL(jj)
c           re=REAL(jj+1)
c           alpha=acos(rp/re)
c           areat=re**2.*(alpha-sin(alpha)*cos(alpha))
c           areac=REAL(0.)
c           if(ri.gt.rp) then
c            alphai=acos(rp/ri)
c            areac=ri**2.*(alphai-sin(alphai)*cos(alphai))
c           end if
c           fraz=(areat-areac)/3.14159265359/(re**2.-ri**2.)
c           mblade(ji+1,i)=mblade(ji+1,i)+mrtp(jj+2,i)*fraz
c          end do
c         end do
c        end do
c        open(10,status='UNKNOWN',file=nank(16)(1:length)//'.mblade')
c        write(10,*)'Matrice blade; theta(deg) X phi(mrad)'
c        write(10,*)'sample= ',nank(11)(1:length2)
c        write(10,'(a3,91(1a,f9.4))')"phi",(tab,tet(i),i=1,91)
c        do ji=1,phiL
c         write(10,'(I3,1a,91(E12.5,1a))')ji-1,tab,
c     +  (mblade(ji,i),tab,i=1,91)
c        end do
c        lose(10)
c       end if
       
  14  FORMAT (f9.1,6(1a,f9.4),8(1a,f9.2))
 100  FORMAT(6(f9.4,1a),f9.1,1a,9(f9.4,1a),:)
 101  FORMAT(16(f9.4,1a),:)
 102  FORMAT(8(f9.4,1a),:)

      RETURN
      END
      
      
      
      SUBROUTINE CALMIS(MC)
      REAL M(17,201,2),E(8,201,2),PAR(60,5),PM(200,5),CNK(16,3),
     + SOL(999,6),PF(7,21),RXY(30,4),MC(14,201),teta(5),vot(5,2),
     + ARSE(500,2),wwl(999),nn(201),kk(201),kklog(201),VNK(16,2),
     + e1(201),e2(201)
      INTEGER DATI(14),S1P2,IXW(20),IKIND(6)
      CHARACTER NANK(16)*256
      PARAMETER(d2r=0.017453293)
      COMMON /VALORI/ M,E,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      DATA IKIND,icol / 1, 2, 3, 4, 5, 6, 1 /
*** subroutine per generazione della matrice con le misure simulate
      
**** caricamento parametri utili da PAR
      do i=1,14
       dati(i)=nint(par(i,4))
      end do
      s1p2=nint(par(27,2)) !polarizzazione
      teta(1)=par(6,1)*d2r
      teta(2)=par(14,1)*d2r
      teta(3)=par(15,1)*d2r
      teta(4)=par(16,1)*d2r
      teta(5)=par(17,1)*d2r
      LEV=nint(RXY(25,3))
      icol=icol+1
      if(icol.gt.15) icol=1
      isvnk=nint(par(9,2))
      if(isvnk.eq.1) sol(1,1)=201
      do i=1,201
       wl=m(7,i,1)
       ev=12400./wl
       if(LEV.EQ.1) then
        wwl(i)=wl
       else
        wwl(i)=ev
       end if
       call COSVNK(1,VNK,i)
       nn(i)=vnk(1,1)
       kk(i)=vnk(1,2)
       e1(i)=vnk(1,1)*vnk(1,1)-vnk(1,2)*vnk(1,2)
       e2(i)=2.*vnk(1,1)*vnk(1,2)
       if(isvnk.eq.1) then
        sol(i+1,1)=wl
        sol(i+1,2)=vnk(1,1)
        sol(i+1,3)=vnk(1,2)
        sol(i+1,4)=vnk(1,1)*0.001
        sol(i+1,5)=vnk(1,2)*0.001
       end if
       if(kk(i).ge.rxy(17,1)) then
        kklog(i)=log10(kk(i))
       else
        kklog(i)=rxy(17,1)
       end if
       par(7,1)=wl
       CALL ASSEMBLER(i,wl,ikind(1),0.,par,pm,vot)
       mc(1,i)=vot(1,1)
       mc(3,i)=vot(2,1)
       mc(5,i)=vot(3,1)
       if(dati(2).ne.0.or.dati(4).ne.0) then
        CALL ASSEMBLER(i,wl,ikind(1),teta(1),par,pm,vot)
        mc(2,i)=vot(1,s1p2)
        mc(4,i)=vot(2,s1p2)
       else
        mc(2,i)=.0
        mc(4,i)=.0
       end if
       if(dati(6).ne.0) then
        CALL ASSEMBLER(i,wl,ikind(2),0.,par,pm,vot)
        mc(6,i)=vot(4,1)
       else
        mc(6,i)=.0
       end if
       do j=1,4
        if(dati(7+2*(j-1)).ne.0.or.dati(8+2*(j-1)).ne.0) then
         CALL ASSEMBLER(i,wl,ikind(2+j),teta(j+1),par,pm,vot)
         mc(7+2*(j-1),i)=vot(5,2) !Delta
         mc(8+2*(j-1),i)=vot(5,1) !Psi
        else
         mc(7+2*(j-1),i)=.0
         mc(8+2*(j-1),i)=.0
        end if
       end do
      end do
***grafico n
      CALL CXY(3,RXY,1,IXW,0)
      call pgsls(1)
      call pgsci(icol)
      call pgline(201,wwl,nn)
***grafico k
      iklog=0
      if(par(31,5).gt.0.5) then
       iklog=1
      end if
      CALL CXY(4,RXY,1,IXW,iklog)
      call pgsls(1)
      call pgsci(icol)
      if(iklog.eq.0) then
       call pgline(201,wwl,kk)
      else
       call pgline(201,wwl,kklog)
      end if
***grafico epsi1 epsi2
      if(nint(par(10,1)).eq.1) then
       call cxy(8,rxy,1,ixw,0)
       call pgsls(1)
       call pgsci(icol)
       call pgline(201,wwl,e1)
       call cxy(9,rxy,1,ixw,0)
       call pgsls(1)
       call pgsci(icol)
       call pgline(201,wwl,e2)
      end if
      RETURN
      END
