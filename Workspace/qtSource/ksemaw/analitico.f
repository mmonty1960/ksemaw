*********    kSEMAW ANALITICO sources   ***************
*
*   hereinafter subroutines used in solution-tracking method
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



      SUBROUTINE AUTOS(q)
      REAL LAM,PAR(60,5),NK(999,9),AMIN(4),AMAX(4),NKF(4),PM(200,5),
     *NKM,LAM1,LAM2,SOL(999,6),NKCL(4),INSTR(20),ELI(8,201,2),DNKF(4),
     *NK0,T(5),EM(14),MIS(17,201,2),NKNEW(999,5),RXY(30,4),
     *TE(5),CNK(16,3),VL(5),BMIN(4),BMAX(4),ARSE(500,2),
     *PF(7,21),minmax(3,2),ps(5,4),ys(5),pss(4),yso(5)
      INTEGER OK,LEV,SG,SIG(512,9),ES,H,IXW(20),icol
      CHARACTER R*1,Q*1,NANK(16)*256,Q1*1,st4*4
      COMMON /VALORI/ MIS,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      EXTERNAL funk
      DATA nk(1,1),icol /0.0, 0/
      do i=1,20
       instr(i)=par(i,3)
      end do
      DO I=1,512
       SIG(I,9)=(-1)**int(I/256.-.01)
       SIG(I,8)=(-1)**int(I/128.-.01)
       SIG(I,7)=(-1)**int(I/64.-.01)
       SIG(I,6)=(-1)**int(I/32.-.01)
       SIG(I,5)=(-1)**int(I/16.-.01)
       SIG(I,4)=(-1)**int(I/8.-.01)
       SIG(I,3)=(-1)**int(I/4.-.01)
       SIG(I,2)=(-1)**int(I/2.-.01)
       SIG(I,1)=-((-1)**I)
      END DO
      DO I=1,14
       EM(I)=0.
      END DO
      VL(1)=MIS(7,1,1)
      VL(2)=MIS(7,201,1)
      VL(3)=100.
      VL(4)=.0
      VL(5)=.0
      T(1)=PAR(6,1)/180.*3.141592654
      do i=1,4
       TE(i+1)=PAR(13+i,1)/180.*3.141592654
       t(i+1)=te(i+1)
      end do 
      OK=0
      iw=nint(par(18,1)) !non verboso 0=n0, 1=yes
      il1=nint(par(26,1))
      il2=nint(par(26,2))
      lam1=mis(7,il1,1)
      lam2=mis(7,il2,1)
      cnk(1,1)=0. ! nk-incogniti
      LEV=NINT(RXY(25,3))
      write(*,*)'|autos> con q=',q
**** quadro comandi
!        WRITE(*,*)'"i" salva soluzione nk da spazio n-k
!        WRITE(*,*)'"a" Esegui ricerca sequenziale'
!        WRITE(*,*)'"0" Azzera le soluzioni memorizzate'
!        WRITE(*,*)'"G" refresh grafici mantenendo i range'
!        WRITE(*,*)'"M" Mantieni le soluzioni '
!        WRITE(*,*)'"x" Ricerca soluzioni nel piano n-k'
!         if(iw.eq.0) then
!          WRITE(*,*)'"w" imposta calcolo verboso'
!         else if(iw.eq.1) then
!          WRITE(*,*)'"w" imposta calcolo NON verboso'
!         end if
       

       IF(Q.EQ.'x') CALL MSR(IIL)

       
       IF(q.eq.'G') THEN
        call pgsci(1)
        do ip=1,2
         ikplot=0
         iww=ip+2
         if(iww.eq.4.and.(par(31,5).gt.0.5)) then
          ikplot=1
         end if
!        if(ip.le.2) then !plot soluzioni memorizzate
          call PLOTNK(ip,3)
!         else !erase merit function graph
!          CALL CXY(iww,RXY,3,IXW,ikplot)
!         end if
!          if(ip.eq.1) then
!           nrfl=nint(par(17,4))
!           IF(nrfl.gt.0) then
!            call PLOTMIS(nrfl,iww,10+nrfl)
!            CALL CXY(3,RXY,1,IXW)!set wl-n window
!           end if
!          end if
        END do
        lam1=par(20,1)
        lam2=par(20,2)
        IL1=1
        IL2=201
        DO I=1,200
         IF(MIS(7,I,1).LE.LAM1.AND.MIS(7,I+1,1).GE.LAM1) IL1=I
         IF(MIS(7,I,1).LE.LAM2.AND.MIS(7,I+1,1).GE.LAM2) IL2=I+1
        END DO
        par(26,1)=il1
        par(26,2)=il2 
        PAR(1,1)=RXY(16,1)
        PAR(1,2)=RXY(16,2)
        PAR(2,1)=RXY(17,1)
        PAR(2,2)=RXY(17,2)
        ICOL=0
       END IF
       
       
*** controllo presenza misure selezionate prima di autocalcola
       if(q.eq.'a'.or.q.eq.'A') then
        nmsel=0
        do ij=1,14
         IF(nint(par(ij,4)).eq.2) nmsel=nmsel+1
        end do
        write(*,*)'nmsel=',nmsel
        if(nmsel.eq.0) then
         write(*,*)'ATTENTION!!!!'
         write(*,*)'  Please select some experimetal measurement'
         write(*,*)
         return
        end if
c       end if
       
       
c       If(q.eq.'a'.or.q.eq.'A') then
        PFM=PAR(5,1)
        NK(1,1)=0.
        NMIS=nmsel
c        NMIS=nint(par(15,4))
        NANGOLI=nint(par(16,4))
        write(*,*)'nmis=',nmis
        do ip=1,3
         minmax(ip,1)=1.E36
         minmax(ip,2)=-1.
        end do
        do i=1,14
         em(i)=.0
        end do
        IF(par(25,3).GT..00001) THEN
         ES=1
        ELSE
         ES=0
        END IF
        IF(par(33,1).LT.0) THEN
         ILIM=2
         JMAX=2
         NKCL(3)=PAR(1,1)
         NKCL(4)=0.
        ELSE
         ILIM=4
         JMAX=4
        END IF
        DO J=1,JMAX
         IF(nint(par(30+J,2)).GE.2) ILIM=ILIM-1
        END DO
        IF(ILIM.LT.NMIS) WRITE(*,*)'N.incognite < N. mis : FM non atten
     *dibile !'
        I2I=2
        if(q.eq.'a') then
         IL=IL1
         ISTEP=1
        else
         IL=IL2
         ISTEP=-1
        end if
        SG=0
        FAV=0.
        NAV=0
        ICOL=ICOL+1
        IF(ICOL.GE.8) ICOL=1
        
*** loop IL
       DO WHILE(IL.GE.IL1.and.il.le.il2)
        LAM=MIS(7,IL,1)
        par(7,1)=MIS(7,IL,1)
        if(iw.eq.1) write(*,*)'***** IL = ',IL,' Lambda = ',LAM,'*****'
        DO J=1,JMAX
         AMIN(J)=10000.
         AMAX(J)=-10000.
         IF(nint(par(30+J,1)).EQ.0) THEN
          IF(IL.EQ.IL1) NKCL(J)=par(30+J,3)
          IF(IL.EQ.IL2.AND.ISTEP.LT.0) NKCL(J)=par(30+J,4)
          IF(nint(par(30+J,2)).GT.0) THEN
           IF(ISTEP.GT.0) THEN
            NKCL(J)=par(30+J,3)
           ELSE
            NKCL(J)=par(30+J,4)
           END IF
          END IF
         END IF
         IF(nint(par(30+J,1)).GE.1.AND.nint(par(30+J,1)).LE.7)THEN
          IF(SG.EQ.0.OR.nint(par(30+J,2)).GT.0) THEN
           IF(J.EQ.1.OR.J.EQ.3) THEN
            ev=12400./lam
            CALL FDISP(nint(par(30+J,1)),ev,PF,PM,CX,CY)
            NKCL(J)=CX
           ELSE
            ev=12400./lam
            CALL FDISP(nint(par(30+J,1)),ev,PF,PM,CX,CY)
            NKCL(J)=CY
           END IF
          END IF
         END IF
         IF(nint(par(30+J,1)).GE.8.AND.
     +      nint(par(30+J,1)).LE.15) THEN
          IF(SG.EQ.0.OR.int(par(30+J,2)+.01).GT.0) THEN
           IF(J.EQ.1.OR.J.EQ.3) THEN
            NKCL(J)=MIS(nint(par(30+J,1)),IL,1)
           ELSE
            NKCL(J)=MIS(nint(par(30+J,1)),IL,2)
           END IF
          END IF
         END IF
         IF(J.EQ.3.AND.int(par(33,2)+.01).EQ.3) par(33,3)=NKCL(3)
         if(nint(par(30+J,2)).eq.2) then
          BMIN(J)=NKCL(J)
          BMAX(J)=NKCL(J)
         end if
        END DO
        FMINIMO=1.E36
        FMASSI=-.1
        if(iw.gt.0) write(*,*)'n_ini = ',nkcl(1),'  k_ini = ',nkcl(2)
        DO J=1,14
         IF(J.LE.6.and.nint(par(j,4)).eq.2) THEN
          EM(j)=MIS(J,IL,1)
         ELSE IF(J.ge.7.and.nint(par(j,4)).eq.2) THEN
          EM(J)=ELI(J-6,IL,1)
          IF(EM(J).LT.-180.) EM(J)=EM(J)+360.
          IF(EM(J).GT.180.) EM(J)=EM(J)-360.
         END IF
        END DO
        DO JJ=1,(2**(NANGOLI+ES))
         NAS=1
         DO J=1,4
          IF(nint(par(2*J+5,4)).EQ.2.OR.nint(par(2*J+6,4)).EQ.2) THEN
           T(J+1)=TE(J+1)+SIG(JJ,NAS)*5.24E-4
           NAS=NAS+1
          END IF
         END DO
         cNK(2,2)=cNK(2,2)+SIG(JJ,NAS)*par(25,3)
         IF(NKCL(1).GE.PAR(1,1).AND.NKCL(1).LE.PAR(1,2).AND.
     *      NKCL(3).GE.PAR(1,1).AND.NKCL(3).LE.PAR(1,2)) THEN
          IFLAG=0
         ELSE
          BMIN(1)=NKCL(1)
          BMAX(1)=NKCL(1)
          BMIN(3)=NKCL(3)
          BMAX(3)=NKCL(3)
          IFLAG=ILIM
         END IF
         IF(ILIM.EQ.0) THEN
          cNK(1,2)=NKCL(1)
          cNK(1,3)=NKCL(2)
          BMIN(1)=NKCL(1)
          BMAX(1)=NKCL(1)
          BMIN(2)=NKCL(2)
          BMAX(2)=NKCL(2)
          IF(nint(par(33,1)).GT.0) THEN
           IF(nint(par(33,2)).EQ.3) THEN
            IF(NKCL(1).GT.par(33,3)) THEN
             NKCL(3)=par(33,3)
            ELSE
             NKCL(3)=NKCL(1)
            END IF
           END IF
           cNK(I2I,2)=NKCL(3)
           cNK(I2I,3)=NKCL(4)
           BMIN(3)=NKCL(3)
           BMAX(3)=NKCL(3)
           BMIN(4)=NKCL(4)
           BMAX(4)=NKCL(4)
          END IF
          CALL FMER(TE,EM,FMIN,IL)
          if(iw.gt.0) write(*,*)'ILIM = 0 FM = ',FMIN
         END IF
         ispx=1
         if(nint(par(25,2)).eq.1) then ! Pre-minimizzazione con Simplex 
*** inizializzazione del simplex      
 6969     if(iw.gt.0)
     +     write(*,*)'Valori iniziali Simplex alla iterazione',ispx,':'
          do jx=1,jmax+1
           do jy=1,jmax
            ps(jx,jy)=nkcl(jy)
            if(jy.eq.jx-1) ps(jx,jy)=nkcl(jy)*par(25,2)
            pss(jy)=ps(jx,jy)
           end do
           ys(jx)=funk(pss,PAR,T,EM,CNK,IL)
           if(ispx.eq.1) yso(jx)=ys(jx)
           if(iw.gt.0) write(*,*) (pss(jy),jy=1,jmax),ys(jx)
          end do
*** lancio simplex
          if(iw.gt.0) write(*,*)'lancio simplex ...'
          ftol=par(26,3)
          call SIMPLEXFMER(ps,ys,jmax,ftol,funk,iter,PAR,T,EM,cNK,IL)
          if(iw.gt.0) then
           write(*,*)'....minimizzazione terminata:'
           do jx=1,jmax+1
            write(*,*)(ps(jx,jy),jy=1,jmax),ys(jx)
           end do
          end if
          do jy=1,jmax
           nkcl(jy)=ps(1,jy)
          end do
*** controllo miglioramento
          dfm=0.
          fmedio=.0
          do jx=1,jmax+1
           dfm=dfm+(yso(jx)-ys(jx))
           fmedio=fmedio+ys(jx)/Real(jmax+1)
          end do
          Den=(ys(2)-ys(1))/(nkcl(1)*par(25,2))
          Dek=(ys(3)-ys(1))/(nkcl(2)*par(25,2))
          if(iw.gt.0) then 
           write(*,*)'dfm=',dfm,' fmedio=',fmedio
           write(*,*)'dFM/dn = ',Den,' dFM/dk = ',Dek
           write(*,*)'dFM/dn*Dn = ',Den*par(21,3),'  dFM/dk*Dk = ',
     +                              Dek*par(22,3)
          end if
          if(fmedio.ge.1..and.dfm.gt.0.1.and.ispx.le.50 ) then
           ispx=ispx+1
           do jx=1,jmax+1
            yso(jx)=ys(jx)
           end do
           goto 6969
          end if
         end if
*** inizio ricerca incrementale
         if(jmax.gt.1) then
          j=1 !inizia a ottimizzare k
         else
          j=0 !altrimenti inizia dall'unica variabile
         end if
         iloop=0
         DO WHILE(IFLAG.LT.ILIM)
          J=J+1
          iloop=iloop+1
          IF(J.GT.JMAX) J=1
          IF(nint(par(30+J,2)).LT.2) THEN
           cNK(1,2)=NKCL(1)
           cNK(1,3)=NKCL(2)
           IF(nint(par(33,1)).GE.0) THEN
            IF(nint(par(33,2)).EQ.3) THEN
             IF(NKCL(1).GT.par(33,3)) THEN
              NKCL(3)=par(33,3)
             ELSE
              NKCL(3)=NKCL(1)
             END IF
            END IF
            cNK(I2I,2)=NKCL(3)
            cNK(I2I,3)=NKCL(4)
           END IF
           CALL FMER(T,EM,FO,IL)
           if(iw.gt.0) then
            write(*,*)'      inizio ricerca incrementale ....'
            write(*,*)'J = ',J,' ILIM = ',ILIM,'  FO = ',FO,
     +       ' con ',NKCL(J)
           end if
           MIN=0
           MAX=0
           NK0=NKCL(J)
           NKM=NKCL(J)
           OK=1
           MIGO=-1
           FINIZ=FO
           FMIN=FO
           icont=0
           icontmax=10000
           q1='y'
           write(st4,'(f4.0)') FO
           if(st4.eq.' nan') then
            write(*,*)'FM = NAN -> STOP a IL = ',IL,'  Lambda = ',LAM
            OK=0
            BMIN(J)=NKCL(J)
            BMAX(J)=NKCL(J)
           end if
           DO WHILE(OK.EQ.1)
            icont=icont+1
            NKCL(J)=NKCL(J)+par(20+J,3)
            if(iw.gt.0) write(*,*)'nkcl -> ',NKCL(J)
            IF(NKCL(1).GE.PAR(1,1).AND.NKCL(1).LE.PAR(1,2).AND.
     *          NKCL(3).GE.PAR(1,1).AND.NKCL(3).LE.PAR(1,2)) THEN
             cNK(1,2)=NKCL(1)
             cNK(1,3)=NKCL(2)
             IF(nint(par(33,1)).GE.0) THEN
              IF(nint(par(33,2)).EQ.3) THEN
               IF(NKCL(1).GT.par(33,3)) THEN
                NKCL(3)=par(33,3)
               ELSE
                NKCL(3)=NKCL(1)
               END IF
              END IF
              cNK(I2I,2)=NKCL(3)
              cNK(I2I,3)=NKCL(4)
             END IF
             CALL FMER(T,EM,FM,IL)
             write(st4,'(f4.0)') FM
             if(iw.gt.0) write(*,*) 'FM = ',FM,' con ',NKCL(J)
             
**** Dimensionamento incremento
             idim=0
             if(FO.gt.1.) then
              den=FO
             else
              den=1.
             end if
             do while((abs(FM-FO)/den.gt.0.5.or.
     +                 abs(FM-FO)/den.lt.0.005).and.
     +                 st4.ne.' nan'.and.idim.lt.1000)
              if(iw.gt.0) write(*,*)'DFM/FO = ',abs(FM-FO)/den
              NKCL(J)=NKCL(J)-par(20+J,3)
              if(abs(FM-FO)/den.gt.0.5.and.st4.ne.' nan') then
               par(20+J,3)=par(20+J,3)/2.
               if(iw.gt.0) write(*,*)'DNK(',J,')/2 = ',par(20+J,3)
              else if(abs(FM-FO)/den.lt.0.005.and.st4.ne.' nan') then
!              else if(abs(FM-FO)/den.lt.0.005.and.st4.ne.' nan'.
!     +               and.abs(par(20+J,3)).lt.0.01) then
               par(20+J,3)=2.*par(20+J,3)
               if(iw.gt.0) write(*,*)'DNK(',J,')*2 = ',par(20+J,3)
              end if
              NKCL(J)=NKCL(J)+par(20+J,3)
              cNK(1,2)=NKCL(1)
              cNK(1,3)=NKCL(2)
              IF(nint(par(33,1)).GE.0) THEN
               IF(nint(par(33,2)).EQ.3) THEN
                IF(NKCL(1).GT.par(33,3)) THEN
                 NKCL(3)=par(33,3)
                ELSE
                 NKCL(3)=NKCL(1)
                END IF
               END IF
               cNK(I2I,2)=NKCL(3)
               cNK(I2I,3)=NKCL(4)
              END IF
              CALL FMER(T,EM,FM,IL)
              write(st4,'(f4.0)') FM
              if(iw.gt.0) write(*,*) 'FM_dim = ',FM,' con ',NKCL(J)
              idim=idim+1
              if(idim.eq.1000) 
     +         write(*,*)'idim=1000 -> STOP dimensionamento'
             end do
****         
             if(st4.eq.' nan') then
              write(*,*)'FM = NAN a IL = ',IL,'  Lambda = ',LAM
              nkcl(j)=nkcl(j)-par(20+j,3)
              if(par(20+j,3).gt..0) then
               max=1
               if(iw.gt.0) write(*,*)'assumo NK_precedente = MAX'
              else
               mim=1
               if(iw.gt.0) write(*,*)'assumo NK_precedente = MIN'
              end if
             else
              IF((FM-1.)*(FO-1.).LE.0.) THEN
               IF(par(20+j,3)*(FO-FM).GT.0.) THEN
                BMIN(J)=NKCL(J)-par(20+j,3)+(1.-FO)/(FM-FO)*par(20+j,3)
                MIN=1
                if(iw.gt.0) write(*,*) 'trovato MIN!'
               ELSE
                BMAX(J)=NKCL(J)-par(20+j,3)+(1.-FO)/(FM-FO)*par(20+j,3)
                MAX=1
                if(iw.gt.0) write(*,*) 'trovato MAX!'
               END IF
              END IF
              IF(FM.LT.FO) THEN
               MIG=1
              ELSE
               MIG=-1
              END IF
              IF(FM.LT.FMIN) THEN
               FMIN=FM
               NKM=NKCL(J)
              END IF
             end if
             IF((MIG.EQ.1.OR.FM.LT.1.).and.st4.ne.' nan') THEN
              OK=1
              FO=FM
              MIGO=MIG
             ELSE
              IF((MIN+MAX).EQ.2) then
               OK=0
               if(iw.gt.0) write(*,*) 'trovati MIN e MAX -> STOP'
              END IF
              IF((MIN+MAX).EQ.1) THEN
               OK=1
               NKCL(J)=NK0
               par(20+j,3)=-par(20+j,3)
               MIGO=-1
               FO=FINIZ
               if(iw.gt.0) write(*,*) 'trovato MIN(MAX) -> inizio'
              END IF
              IF((MIN+MAX).EQ.0) THEN
               IF(MIGO.EQ.1.AND.MIG.EQ.-1) THEN
                OK=0
                BMIN(J)=NKCL(J)-par(20+j,3)
                BMAX(J)=BMIN(J)
                if(iw.gt.0) write(*,*) 'soluzione puntuale -> STOP'
               ELSE
                OK=1
                par(20+j,3)=-par(20+j,3)
                NKCL(J)=NKCL(J)+par(20+j,3)
                MIGO=1
                if(iw.gt.0) write(*,*) 'inverto la direzione di ricerca'
               END IF
              END IF
             END IF
            ELSE
             OK=0
             BMIN(J)=NKCL(J)
             BMAX(J)=NKCL(J)
             if(iw.gt.0) write(*,*) 'soluzione fuori range -> STOP'
            END IF
            if(icont.gt.icontmax) then
             write(*,*)'Stop: N. iterazioni > 10000 !'
             OK=0
            end if
           END DO !do while(OK.eq.1)
           NKCL(J)=(BMIN(J)+BMAX(J))/2
           if(iw.gt.0) write(*,*) 'Min = ',BMIN(J),' Max = ',BMAX(J) 
           IF(NKCL(1).GE.PAR(1,1).AND.NKCL(1).LE.PAR(1,2).AND.
     *         NKCL(3).GE.PAR(1,1).AND.NKCL(3).LE.PAR(1,2)) THEN
c            TOL=(BMAX(J)-BMIN(J))*par(26,3)
c            IF(TOL.LT.2.*ABS(par(20+J,3))) TOL=2.*ABS(par(20+J,3))
             tol=par(26,3)*ABS(par(20+J,3))
            if(iw.gt.0) write(*,*)'|NKnew-NKold| = ',
     + ABS((BMIN(J)+BMAX(J))/2.-NK0),' tol = ',tol
            IF(ABS((BMIN(J)+BMAX(J))/2.-NK0).LE.TOL) THEN
             IFLAG=IFLAG+1
             if(iw.gt.0) write(*,*)'|NKnew-NKold| < tol -> iflag+1'
            ELSE
             IFLAG=0
             if(iw.gt.0) write(*,*)'|NKnew-NKold| > tol -> iflag=0'
            END IF
           ELSE
            IFLAG=ILIM
            if(iw.gt.0) write(*,*)'NK fuori range -> iflag=ilim'
           END IF
           if(q1.ne.'y') IFLAG=ILIM
           if(iw.gt.0) write(*,*)'iflag = ',iflag
          ELSE
           BMIN(J)=NKCL(J)
           BMAX(J)=NKCL(J)
          END IF
          if(iloop.gt.1E+03) then
           write(*,*)'Stop: iflag=0 dopo ',iloop,' loops!'
           iflag=ilim
          end if
         END DO
         DO J=1,JMAX
          IF(BMIN(J).LT.AMIN(J)) AMIN(J)=BMIN(J)
          IF(BMAX(J).GT.AMAX(J)) AMAX(J)=BMAX(J)
         END DO
         IF(FMIN.LT.FMINIMO) FMINIMO=FMIN
         IF(FMIN.GT.FMASSI) FMASSI=FMIN
        END DO
        FMEDIO=PFM*FMINIMO+(1.-PFM)*FMASSI
        IF(SQRT(FMEDIO).GT.minmax(3,2)) minmax(3,2)=SQRT(FMEDIO)
        IF(SQRT(FMEDIO).LT.minmax(3,1)) minmax(3,1)=SQRT(FMEDIO)
        FAV=FAV+FMEDIO
        NAV=NAV+1
        DO J=1,JMAX
         NKF(J)=(AMIN(J)+AMAX(J))/2.
         DNKF(J)=(AMAX(J)-AMIN(J))/2.
         if(iw.gt.0)
     +    write(*,*) 'J = ',j,' Min = ',amin(j),' Max = ',amax(j)
         NKCL(J)=NKF(J)
         IF((NKF(J)+DNKF(J)).LT..0.AND.
     *       NKF(1).GE.PAR(1,1).AND.NKF(1).LE.PAR(1,2).AND.
     *       NKF(3).GE.PAR(1,1).AND.NKF(3).LE.PAR(1,2)) THEN
          XX=LAM
          IF(nINT(RXY(25,3)).EQ.2) XX=12400./XX
         END IF
        END DO
        IF(nint(par(33,1)).LT.0) NKF(3)=PAR(1,1)
        IF(NKF(1).GE.PAR(1,1).AND.NKF(1).LE.PAR(1,2).AND.
     *      NKF(3).GE.PAR(1,1).AND.NKF(3).LE.PAR(1,2)) THEN
         NK(1,1)=NK(1,1)+1
         ISOL=nINT(NK(1,1))
         NK(ISOL+1,1)=LAM
         NK(ISOL+1,2)=NKF(1)
         NK(ISOL+1,3)=NKF(2)
         NK(ISOL+1,4)=DNKF(1)
         NK(ISOL+1,5)=DNKF(2)
         NK(ISOL+1,6)=NKF(3)
         NK(ISOL+1,7)=NKF(4)
         NK(ISOL+1,8)=DNKF(3)
         NK(ISOL+1,9)=DNKF(4)
         do ip=1,2
          if(nkf(ip).lt.minmax(ip,1)) minmax(ip,1)=nkf(ip)
          if(nkf(ip).gt.minmax(ip,2)) minmax(ip,2)=nkf(ip)
         end do
         if(iw.gt.0)write(*,*)'n = ',nkf(1),' k = ',nkf(2)
*** grafico soluzioni calcolate
         XX=LAM
         IF(nINT(RXY(25,3)).EQ.2) XX=12400./XX
         do ip=1,2!3
          ikplot=0
          iww=ip+2
          if(iww.eq.4.and.par(31,5).gt.0.5) then
            ikplot=1
           end if
          if(ixw(iww).le.0) then
           CALL CXY(iww,RXY,0,IXW,ikplot)
          else
           call pgslct(ixw(iww))
          end if
          CALL PGQWIN(XL,XH,YL,YH)
          CALL pgsci(ICOL)
          IF(IP.EQ.1.OR.IP.EQ.2) THEN
           Y11=NKF(ip)-DNKF(ip)
           Y21=NKF(ip+2)-DNKF(ip+2)
           if(ikplot.eq.1) then
            Y11=log10(Y11)
            Y21=log10(Y21)
           end if
           IF(Y11.LT.YL) Y11=YL
           IF(Y21.LT.YL) Y21=YL
           Y12=NKF(ip)+DNKF(ip)
           Y22=NKF(ip+2)+DNKF(ip+2)
           if(ikplot.eq.1) then
            Y22=log10(Y22)
            Y12=log10(Y12)
           end if
           IF(Y12.GT.YH) Y12=YH
           IF(Y22.GT.YH) Y22=YH
           CALL pgmove(XX,Y11)
           IF(Y11.LE.YH.AND.Y12.GE.YL) CALL pgdraw(XX,Y12)
           IF(nint(par(33,1)).GE.0) THEN
            CALL pgmove(XX,Y21)
            IF(Y21.LE.YH.AND.Y22.GE.YL) CALL pgdraw(XX,Y22)
           END IF
          ELSE if(ip.eq.3) then
           FMS=SQRT(FMEDIO)
           IF(IL.EQ.IL1) THEN
            XXOLD=XX
            FMSOLD=FMS
           END IF
           IF(FMS.GE.YL.AND.FMS.LE.YH) THEN
            CALL pgmove(XXOLD,FMSOLD)
            CALL pgdraw(XX,FMS)
            XXOLD=XX
            FMSOLD=FMS
           END IF
          END IF
         end do
         IL=IL+ISTEP
         SG=1
        ELSE
         SG=0
         DO J=1,3,2
          NDOMANDA=0
          IF(NKF(J).LT.PAR(1,1).OR.NKF(J).GT.PAR(1,2).AND.
     *                                            NDOMANDA.EQ.0) THEN
           IF(nint(par(30+J,1)).EQ.0.AND.nint(par(30+J,2)).EQ.0) THEN
            WRITE(*,*)'n(WL) is out of range! WL=',lam
            WRITE(*,*)'Right-button = Sospendi la ricerca '
            WRITE(*,*)'Left-button = Riprendi dal valore del puntatore'
            call pgslct(ixw(3))
            cursore=pgcurs(xval,yval,r)
            write(*,*)'r=',r
            if(r.eq.'A') then
             NKCL(1)=yval
             SG=1
             IL=IL+ISTEP
            else
             IL=-100
            end if
           ELSE
            IL=IL+ISTEP
           END IF
          END IF
         END DO
        END IF
       END DO
       RXY(16,3)=minmax(1,1)
       RXY(16,4)=minmax(1,2)
       RXY(17,3)=minmax(2,1)
       RXY(17,4)=minmax(2,2)
       RXY(15,3)=minmax(3,1)
       RXY(15,4)=minmax(3,2)
       WRITE(*,*)'<FM> = ',SQRT(FAV/(NAV-0.999))
       end if !(q.eq.'a'
        
        
        
       IF(Q.EQ.'i') THEN
        iiww=1
        iiok=0
        do while(iiok.ne.1.and.iiww.le.20)
         if(ixw(iiww).eq.1) then
          iiok=1
         else
          iiww=iiww+1
         end if
        end do
        if(iiok.eq.0) goto 9
        call pgslct(iiww)
        write(*,*) 'select n-k_cross with 2 OBLIQUE points Nw=',iiww
        cursore=pgcurs(xval1,yval1,r)
        cursore=pgcurs(xval2,yval2,r)
        SOL(1,1)=SOL(1,1)+1
        LSOL=nINT(SOL(1,1))
        if(LSOL.le.999) then
         SOL(LSOL+1,1)=par(7,1)
c         write(*,*) (xval1+xval2)/2.
         SOL(LSOL+1,2)=(xval1+xval2)/2.
c         write(*,*) (yval1+yval2)/2.
         SOL(LSOL+1,3)=(yval1+yval2)/2.
c         write(*,*) abs(xval1-xval2)/2.
         SOL(LSOL+1,4)=abs(xval1-xval2)/2.
c         write(*,*) abs(yval1-yval2)/2.
         SOL(LSOL+1,5)=abs(yval1-yval2)/2.
         do i=1,5
          par(30,i)=SOL(LSOL+1,i)
         end do
        else
         write(*,*)'nk-file is full! Please make space before save!'
        end if
 9     END IF
        
        
        
       IF(Q.EQ.'M'.and.nint(NK(1,1)).gt.0) THEN
        ISV=nint(par(18,2))
        isol=nint(SOL(1,1))
        isolnew=nint(NK(1,1))
        irima=998-isol
        IF(ISV.EQ.0) THEN
         write(*,*)'Salvataggio soluzioni in temp:'
         WRITE(*,*)' N. soluzioni presenti = ',isol
         write(*,*)' N. nuove soluzioni =', isolnew
         write(*,*)' spazio disponibile=', irima
         rss=isol
         if(isolnew.gt.irima) then
          isolnew=irima
          write(*,*)'ATTENZIONE: salvataggio parziale!'
         end if
**** ordinamento wl crescente delle soluzioni         
         DO I=2,isolnew+1
          WWM=1.E30
          DO H=2,isolnew+1
           IF(NK(H,1).GE..0.AND.NK(H,1).LT.WWM) THEN
            WWM=NK(H,1)
            NORD=H
           END IF
          END DO
          NKNEW(I,1)=NK(NORD,1)
          NKNEW(I,2)=NK(NORD,2)
          NKNEW(I,3)=NK(NORD,3)
          NKNEW(I,4)=NK(NORD,4)
          NKNEW(I,5)=NK(NORD,5)
          NK(NORD,1)=-1.
         END DO
         DO I=2,isolnew+1
          NK(I,1)=NKNEW(I,1)
          NK(I,2)=NKNEW(I,2)
          NK(I,3)=NKNEW(I,3)
          NK(I,4)=NKNEW(I,4)
          NK(I,5)=NKNEW(I,5)
         END DO
**** fine ordinamento
         DO I=1,isolnew
          SOL(1,1)=SOL(1,1)+1
          LSOL=nINT(SOL(1,1))
          SOL(LSOL+1,1)=NK(i+1,1)
          SOL(LSOL+1,2)=NK(i+1,2)
          SOL(LSOL+1,3)=NK(i+1,3)
          SOL(LSOL+1,4)=NK(i+1,4)
          SOL(LSOL+1,5)=NK(i+1,5)
          SOL(LSOL+1,6)=1. !dato abilitato nel fit    
         END DO
         par(38,1)=sol(2,1)
         par(38,2)=sol(LSOL+1,1)
         par(21,2)=sol(1,1)
         par(29,2)=sol(1,1)
         WRITE(*,*)' => Lmin =',sol(2,1),' Lmax =',sol(LSOL+1,1),
     +     ' Nsol = ',nint(sol(1,1))
         
        else 
         isv=isv+7
         DO I=IL1,IL2
          IF(ABS(MIS(7,I,1)-NK(I-IL1+2,1)).LT..001) THEN
           MIS(isv,I,1)=NK(I-IL1+2,2)
           MIS(isv,I,2)=NK(I-IL1+2,3)
          else
           WRITE(*,*)'ATTENZIONE : MISMATCHING TRA GLI INDICI'
          END IF
         END DO
        END IF
       END IF
       
       if(q.eq.'0') then
        sol(1,1)=0.
        nk(1,1)=0.
       end if
       
      RETURN
      END



      SUBROUTINE MSR(N)
      REAL MIS(17,201,2),ELI(8,201,2),PAR(60,5),ARSE(500,2),
     *     NK(999,6),CNK(16,3),PF(7,21),RXY(30,4),PM(200,5)
      INTEGER DATO(14),DF(3),LEV,IXW(20)
      CHARACTER  NANK(16)*256,swl*10,SMIS(6)*2
      COMMON /VALORI/ MIS,ELI,PAR,PM,CNK,NK,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      DATA SMIS /'Tn','Tp','Rn','Rp','R1','Ap'/
      LEV=nINT(RXY(25,3))
      df(1)=nint(par(22,2))
      df(2)=nint(par(23,1))
      df(3)=nint(par(23,2))
      do i=1,14
       if(nint(par(i,4)).eq.2) then
        dato(i)=1
       else
        dato(i)=0
       end if
      end do
      CALL pgsls(1)
      DBSB=.0022
      DRCSRC=.005
      n=1
      DO WHILE((MIS(7,N,1).LT.PAR(7,1)).and.(n.lt.201))
       N=N+1
      END DO
      PAR(7,1)=MIS(7,N,1)
      WRITE(*,*)'DATO DISPONIBILE A LAMBDA = ',PAR(7,1)
      do i=1,6
       if(dato(i).eq.1) then
        write(*,*) smis(i),' =',mis(i,n,1)*100.,
     +             ' +-',mis(i,n,2)*100.,' %'
       end if
      end do
      IF(DF(3).EQ.1) THEN
       WRITE(*,*)'     LAMBDA = ',PAR(7,1)
       DO I=7,14
        IPD=I-INT(I/2.)*2
        IF(DATO(I).EQ.1) THEN
         IF(IPD.EQ.1) THEN
          WRITE(*,*)'DELTA = ',ELI(I-6,N,1),' +- ',ELI(I-6,N,2)
         ELSE
          WRITE(*,*)'PSI = ',ELI(I-6,N,1),' +- ',ELI(I-6,N,2)
         END IF
        END IF
       END DO
      END IF
      Y=mis(1,n,1)+mis(3,n,1)
      IF(DATO(1).ge.1.AND.DATO(3).ge.1) WRITE(*,*)'Tn + Rn = ',Y
      Y=mis(1,n,1)+mis(5,n,1)
      IF(DATO(1).ge.1.AND.DATO(5).ge.1) WRITE(*,*)'Tn + R1 = ',Y
      Y=mis(2,n,1)+mis(4,n,1)
      IF(DATO(2).ge.1.AND.DATO(4).ge.1) WRITE(*,*)'Tp + Rp = ',Y
      CALL pgsls(1)
      call CXY(1,RXY,3,IXW,0)
      PAR(1,1)=RXY(16,1)
      PAR(1,2)=RXY(16,2)
      PAR(2,1)=RXY(17,1)
      PAR(2,2)=RXY(17,2)
      PAR(20,1)=RXY(22,1)
      PAR(20,2)=RXY(22,2)
*** scrittura label con wl
      write(swl,'(f7.1)') par(7,1)
      xmin=(par(1,1)+par(1,2))/2.
      ymin=(par(2,1)+par(2,2))/2.
      call PGPTXT (xmin,ymin,0.,0.5,'wl = '//swl)
      CALL CURVA(dato,N)
      RETURN
      END
      
      
      
      SUBROUTINE CURVA(dato,N)
      CHARACTER NANK(16)*256
      REAL PM(200,5),ELI(8,201,2),PAR(60,5),SOL(999,6),M(17,201,2),
     *     CNK(16,3),PF(7,21),RXY(30,4),ARSE(500,2)
      INTEGER dato(14),IXW(20)
      COMMON /VALORI/ M,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      do imis=1,6
       if(dato(imis).ge.1) then
        CALL pgsci(imis)
        DT=M(imis,N,2)
        M(imis,N,1)=M(imis,N,1)-DT
        CALL SOLVE(IMIS,N)
        M(imis,N,1)=M(imis,N,1)+2.*DT
        CALL SOLVE(IMIS,N)
        M(imis,N,1)=M(imis,N,1)-DT
       end if
      end do
      DO imis=7,14
       IF(dato(imis).EQ.1) THEN
        i=imis-6
        IF(I.EQ.1.OR.I.EQ.3.OR.I.EQ.5.OR.I.EQ.7) CALL pgsci(7)
        IF(I.EQ.2.OR.I.EQ.4.OR.I.EQ.6.OR.I.EQ.8) CALL pgsci(8)
        dt=eli(i,n,2)
        ELI(I,N,1)=ELI(I,N,1)-dt
        IF(ELI(I,N,1).LT.-180.) ELI(I,N,1)=ELI(I,N,1)+360.
        IF(ELI(I,N,1).GT.180.) ELI(I,N,1)=ELI(I,N,1)-360.
        CALL SOLVE(IMIS,N)
        ELI(I,N,1)=ELI(I,N,1)+2.*dt
        IF(ELI(I,N,1).LT.-180.) ELI(I,N,1)=ELI(I,N,1)+360.
        IF(ELI(I,N,1).GT.180.) ELI(I,N,1)=ELI(I,N,1)-360.
        CALL SOLVE(IMIS,N)
        ELI(I,N,1)=ELI(I,N,1)-dt
        IF(ELI(I,N,1).LT.-180.) ELI(I,N,1)=ELI(I,N,1)+360.
        IF(ELI(I,N,1).GT.180.) ELI(I,N,1)=ELI(I,N,1)-360.
       END IF
      END DO
      RETURN
      END










      SUBROUTINE SOLVE(IMIS,NN)
      REAL PAR(60,5),LAM,CNK(16,3),M(17,201,2),NCL,KCL,PF(7,21),
     +     RXY(30,4),ELI(8,201,2),SOL(999,6),PM(200,5),VOT(5,2),
     +     ARSE(500,2)
      INTEGER S1P2,IP(2),IXW(20)
      CHARACTER  NANK(16)*256
      DOUBLE PRECISION F(3)
      COMMON /VALORI/ M,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
C   NN indice wl

*** tipologia misura per selezione mezzo di ingresso
c         (1,2,3,4,5,6,7,8,9,10,11,12,13,14)
      goto(1,2,3,4,5,6,7,7,9, 9,11,11,13,13) imis
 1     ikind=1
       tetar=.0
       ivot=1
      goto 14
 2     ikind=1
       tetar=par(6,1)
       ivot=1
      goto 14
 3     ikind=1
       tetar=.0
       ivot=2
      goto 14
 4     ikind=1
       tetar=par(6,1)
       ivot=2
      goto 14
 5     ikind=1
       tetar=.0
       ivot=3
      goto 14
 6     ikind=2
       tetar=.0
       ivot=4
      goto 14
 7     ikind=3
       tetar=par(14,1)
       ivot=5
       IPD=IMIS-INT(IMIS/2.)*2
      goto 14
 9     ikind=4
       tetar=par(15,1)
       ivot=5
       IPD=IMIS-INT(IMIS/2.)*2
      goto 14
 11    ikind=5
       tetar=par(16,1)
       ivot=5
       IPD=IMIS-INT(IMIS/2.)*2
      goto 14
 13    ikind=6
       tetar=par(17,1)
       ivot=5
       IPD=IMIS-INT(IMIS/2.)*2
 14   S1P2=nint(par(27,2))
c      LAM=PAR(7,1)
      LAM=M(7,NN,1)
      TETAR=TETAR/180.*3.141592654
      IP(1)=1
      IP(2)=2
      DO IZ=2,1,-1
       ZLIM=(PAR(IP(IZ),2)-PAR(IP(IZ),1))/500.
       DW=(PAR(IP(3-IZ),2)-PAR(IP(3-IZ),1))/50.
       J=0
       DO I=1,51
        W=PAR(IP(3-IZ),1)+(I-1)*DW
        DZ=PAR(IP(IZ),2)-PAR(IP(IZ),1)
        DO L=1,2
         Z=PAR(IP(IZ),1)+(L-1)*DZ
         IF(IZ.EQ.2) THEN
          NCL=W
          KCL=Z
         ELSE
          NCL=Z
          KCL=W
         END IF
         cNK(1,2)=NCL
         cNK(1,3)=KCL
         call ASSEMBLER(nn,LAM,ikind,tetar,par,pm,vot)
         F(L)=0.
         IF(IMIS.LE.6) THEN
          F(L)=F(L)+DBLE(m(imis,nn,1)-vot(ivot,s1p2))
         ELSE
          F(L)=F(L)+SIND1((ELI(IMIS-6,NN,1)-vot(ivot,IPD+1))/2.)
         END IF
        END DO
        IF(F(1)*F(2).LT.0) THEN
         Z=PAR(IP(IZ),1)
         DZ=DZ/2.
         DO WHILE(ABS(DZ).GT.ZLIM)
          Z=Z+DZ
          IF(IZ.EQ.2) THEN
           NCL=W
           KCL=Z
          ELSE
           NCL=Z
           KCL=W
          END IF
          cNK(1,2)=NCL
          cNK(1,3)=KCL
          call ASSEMBLER(nn,LAM,ikind,tetar,par,pm,vot)
          F(3)=0.
          IF(IMIS.LE.6) THEN
           F(3)=F(3)+DBLE(m(imis,nn,1)-vot(ivot,s1p2))
          ELSE
           F(3)=F(3)+SIND1((ELI(IMIS-6,NN,1)-vot(ivot,IPD+1))/2.)
          END IF
          IF(F(3)*F(1).LT.0) THEN
           DZ=-DZ/2.
          ELSE
           DZ=DZ/2
          END IF
          F(1)=F(3)
         END DO
         YCL=KCL
         IF(J.EQ.0) CALL pgmove(NCL,YCL)
         IF(J.EQ.1) CALL pgdraw(NCL,YCL)
         J=1
        ELSE
         J=0
        END IF
       END DO
      END DO
      RETURN
      END




      SUBROUTINE FMER(T,EXME,FM,IL)
      REAL wl,PM(200,5),ELI(8,201,2),EXME(14),CNK(16,3),MIS(17,201,2),
     * PAR(60,5),T(5),PF(7,21),RXY(30,4),SOL(999,6),vot(5,2),
     * ARSE(500,2)
      INTEGER S1P2,IXW(20)
      CHARACTER NANK(16)*256
      COMMON /VALORI/ MIS,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      s1p2=nint(par(27,2))
      wl=MIS(7,IL,1)
      FM=.0
      NMIS=0
      IF(nint(par(1,4)).eq.2.OR.nint(par(3,4)).eq.2.OR.
     +   nint(par(5,4)).eq.2) THEN
       CALL ASSEMBLER(il,wl,1,0.,par,pm,vot)
       do i=1,3
        j=2*i-1
        if(nint(par(j,4)).eq.2) then
         fm=fm+((vot(i,s1p2)-exme(j))/mis(j,il,2))**2.
         nmis=nmis+1
        end if
       end do 
      END IF
      IF(nint(par(2,4)).eq.2.OR.nint(par(4,4)).eq.2) THEN
       CALL ASSEMBLER(il,wl,1,t(1),par,pm,vot)
       do i=1,2
        j=2*i
        if(nint(par(j,4)).eq.2) then
         fm=fm+((vot(i,s1p2)-exme(j))/mis(j,il,2))**2.
         nmis=nmis+1
        end if
       end do 
      END IF
      IF(nint(par(6,4)).eq.2) THEN
       CALL ASSEMBLER(il,wl,2,0.,par,pm,vot)
       FM=FM+(vot(4,1)-EXME(6))**2.
       NMIS=NMIS+1
      END IF
      DO I=1,4
       IF(nint(par(5+2*i,4)).EQ.2.OR.nint(par(6+2*i,4)).EQ.2) THEN
        CALL ASSEMBLER(il,wl,i+2,t(1+i),par,pm,vot)
        DE=EXME(5+2*I)
        DC=vot(5,2)
        DDE=ELI(5+2*I-6,IL,2)*3.1416/180.
        PE=EXME(6+2*I)
        PC=vot(5,1)
        DP=ELI(6+2*I-6,IL,2)*3.1416/180.
        DFMD=(SIND1(DC)-SIND1(DE))*(SIND1(DC)-SIND1(DE))
        DFMD=DFMD+(COSD1(DC)-COSD1(DE))*(COSD1(DC)-COSD1(DE))
        DFMD=DFMD/DDE/DDE
        DFMP=(SIND1(PC)-SIND1(PE))*(SIND1(PC)-SIND1(PE))
        DFMP=DFMP+(COSD1(PC)-COSD1(PE))*(COSD1(PC)-COSD1(PE))
        DFMP=DFMP/DP/DP
        if(nint(par(5+2*i,4)).eq.2) then
         FM=FM+DFMD
         nmis=nmis+1
        end if
        if(nint(par(6+2*i,4)).eq.2) then
         FM=FM+DFMP
         nmis=nmis+1
        end if
       END IF
      END DO
      FM=FM/(NMIS+.0001)
      RETURN
      END



      SUBROUTINE SIMPLEXFMER(p,y,ndim,ftol,funk,iter,PAR,T,EM,CNK,IL)
      INTEGER ITMAX
      REAL TINY
      PARAMETER (ITMAX=500,TINY=1.e-10) 
c Maximum allowed dimensions and function evaluations, and a small 
c number.
      INTEGER iter,ndim,IL
      REAL ftol,p(5,4),y(5),funk,T(5),EM(14),CNK(16,3),PAR(60,5)
      EXTERNAL funk
C     USES amotry,funk
c  Multidimensional minimization of the function funk(x) where 
c  x(1:ndim) is a vector in ndim dimensions, by the downhill simplex 
c  method of Nelder and Mead. The matrix p(1:ndim+1,1:ndim) is input.
c  Its ndim+1 rows are ndim-dimensional vectors which are the vertices
c  of the starting simplex. 
c  Also input is the vector y(1:ndim+1), whose components must be 
c  pre-initialized to the values of funk evaluated at the ndim+1 
c  vertices (rows) of p; and ftol the fractional convergence tolerance
c  to be achieved in the function value (n.b!). 
c  On output, p and y will have been reset to ndim+1 new points all
c  within ftol of a minimum function value, and iter gives the number 
c  of function evaluations taken.
      INTEGER i,ihi,ilo,inhi,j,m,n
      REAL rtol,sum,swap,ysave,ytry,psum(4),amotry
      iw=0
      iter=0
 1    do n=1,ndim !Enter here when starting or have just overall 
c                  contracted.
       sum=0. !Recompute psum.
       do m=1,ndim+1
        sum=sum+p(m,n)
       end do
       psum(n)=sum
      end do
      
 2    ilo=1 !Enter here when have just changed a single point.
      if (y(1).gt.y(2)) then !Determine which point is the highest
c                             (worst), next-highest,
       ihi=1 !and lowest (best),
       inhi=2
      else
       ihi=2
       inhi=1
      end if
      do i=1,ndim+1 !by looping over the points in the simplex.
       if(y(i).le.y(ilo)) ilo=i
       if(y(i).gt.y(ihi)) then
        inhi=ihi
        ihi=i
       else if(y(i).gt.y(inhi)) then
        if(i.ne.ihi) inhi=i
       end if
      end do
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
c Compute the fractional range from highest to lowest and return if 
c satisfactory.
      if (rtol.lt.ftol) then !If returning, put best point and value 
c                             in slot 1.
       swap=y(1)
       y(1)=y(ilo)
       y(ilo)=swap
       do n=1,ndim
        swap=p(1,n)
        p(1,n)=p(ilo,n)
        p(ilo,n)=swap
       end do
       return
      endif
c      if(iter.ge.ITMAX) pause 'ITMAX exceeded in amoeba'
      if(iter.ge.ITMAX) return
      iter=iter+2
c  Begin a new iteration. First extrapolate by a factor -1 through the 
c  face of the simplex across from the high point, i.e., reflect the 
c  simplex from the high point.
      ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0,PAR,T,EM,cNK,IL)
      if(iw.gt.0) write(*,*)'    reflecting from HP ytry = ',rtol
      if (ytry.le.y(ilo)) then
c Gives a result better than the best point, so try an additional 
c extrapolation by a factor 2.
       ytry=amotry(p,y,psum,ndim,funk,ihi,2.0,PAR,T,EM,cNK,IL)
       if(iw.gt.0) write(*,*)'     extrapolation*2 ytry = ',rtol
      else if (ytry.ge.y(inhi)) then
c  The reflected point is worse than the second-highest, so look for 
c  an intermediate lower point,i.e., do a one-dimensional contraction.
       ysave=y(ihi)
       ytry=amotry(p,y,psum,ndim,funk,ihi,0.5,PAR,T,EM,cNK,IL)
       if(iw.gt.0) write(*,*)'     contraction ytry = ',rtol
       if (ytry.ge.ysave) then !Can not seem to get rid of that high 
c                               point. Better contract
        do i=1,ndim+1 !around the lowest (best) point.
         if(i.ne.ilo)then
          do j=1,ndim
           psum(j)=0.5*(p(i,j)+p(ilo,j))
           p(i,j)=psum(j)
          end do
          y(i)=funk(psum,PAR,T,EM,VNK,IL)
          if(iw.gt.0) write(*,*)'   ',(p(i,j),j=1,ndim),y(i)
         end if
        end do
        iter=iter+ndim !Keep track of function evaluations.
        goto 1 !Go back for the test of doneness and the next iteration.
       end if
      else
       iter=iter-1 !Correct the evaluation count.
      endif
      goto 2
      END



      FUNCTION amotry(p,y,psum,ndim,funk,ihi,fac,PAR,T,EM,VNK,IL)
      INTEGER ihi,mp,ndim,np,IL
      REAL amotry,fac,p(5,4),psum(4),y(5),funk,T(5),EM(14),
     +     VNK(16,2),PAR(60,5)
      EXTERNAL funk
C     USES funk
c  Extrapolates by a factor fac through the face of the simplex 
c  across from the high point, tries it, and replaces the high point 
c  if the new point is better.
      INTEGER j
      REAL fac1,fac2,ytry,ptry(4)
      mp=ndim+1
      np=ndim
      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do j=1,ndim
       ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
      end do
      ytry=funk(ptry,PAR,T,EM,VNK,IL) !Evaluate the function at the 
c                                         trial point.
c      write(*,*)' |Amotry> ptry:',(ptry(j),j=1,ndim),ytry
      if (ytry.lt.y(ihi)) then !If it is better than the highest, 
c                               then replace the highest.
       y(ihi)=ytry
       do j=1,ndim
        psum(j)=psum(j)-p(ihi,j)+ptry(j)
        p(ihi,j)=ptry(j)
       end do
      endif
      amotry=ytry
      return
      END
      
      
      
      FUNCTION funk(p,PAR,T,EM,CNK,IL)
      INTEGER IL
      REAL p(4),T(5),EM(14),PAR(60,5),NKCL(4),CNK(16,3)
**
*** funzione di merito utilizzata da simplex in AUTOS
**
      IF(par(33,1).LT.0) THEN
       jmax=2
      else
       jmax=4
      end if
      do i=1,jmax
       nkcl(i)=p(i)
      end do
*** scrittura variabili su VNK 
      cNK(1,2)=NKCL(1)
      cNK(1,3)=NKCL(2)
      IF(nint(par(33,1)).GE.0) THEN
       IF(nint(par(33,2)).EQ.3) THEN
        IF(NKCL(1).GT.par(33,3)) THEN
         NKCL(3)=par(33,3)
        ELSE
         NKCL(3)=NKCL(1)
        END IF
       END IF
       i2i=2
       if(abs(nint(par(19,1))).eq.4) i2i=3
       cNK(I2I,2)=NKCL(3)
       cNK(I2I,3)=NKCL(4)
      END IF
      CALL FMER(T,EM,FM,IL)
      funk=FM
      RETURN
      END
