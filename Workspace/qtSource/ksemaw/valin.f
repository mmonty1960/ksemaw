*********    kSEMAW VALIN sources   ***************
*
*   hereinafter MAIN and subroutines for
*     - loading of measurements and files
*     - global resampling over 201 points
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
C
C   PAR(60,5) matrice di trasporto parametri vari
C
C      (1,1),(1,2)........NMIN, NMAX (range di ricerca soluzioni)
C      (2,1),(2,2)........KMIN, KMAX (  "   "     "        "    )
C      (3,1),(3,2)........nSUB,kSUB (quando sono noti)
C      (4,1),(4,2)........LAMBDAmin,LAMBDAmax files dati
C      (5,1),(5,2)........soglia accettazione soluzione (AUTOS), incr. Df
C      (6,1),(6,2)........TETA, incr. TETA
C      (7,1),(7,2)........LAM, icolor
C      (8,1),(8,2)........jobtot , jobview
C      (9,1),(9,2)........if 1 k>=0, if 1 nk_sim->temp from scratch
C     (10,1),(10,2).......if 1 plot eps1 eps2 , /
C     (11,1),(11,2)......./ , /
C     (12,1),(12,2).......ifinestra , x(Angstrom) !leggi cursore
C     (13,1),(13,2).......x(eV) , y               !leggi cursore
C     (14,1),(14,2).......TETAELIS1 , incr.
C     (15,1),(15,2).......TETAELIS2 , incr.
C     (16,1),(16,2).......TETAELIS3 , incr.
C     (17,1),(17,2).......TETAELIS4 , incr.
C     (18,1),(18,2).......verboso (0=no, 1=si), salvataggio nk in 0 -> sol
c                                                               1 -> mis( 8,i,1=n 2=k)
c                                                               .............
c                                                               8 -> mis(15,i,1=n 2=k)
C     (19,1),(19,2).......0 => nulla 1 => Sim -> Exp , dPSI
C     (20,1),(20,2).......dMIN, dMAX (range di ricerca soluzioni)
C     (21,1),(21,2).......mwla(N. wl abilitate), mda(N. nk abilitati)
C     (22,1),(22,2).......mr(0=>R assolute, 1=> R*Rrif) , DF(1)
C     (23,1),(23,2).......DF(2), DF(3)
C     (24,1),(24,2).......IUVIR,i-lambda
C     (25,1),(25,2).......tol,controllo SIMPLEX (= 0 no, >1 Si)
C     (26,1),(26,2).......IL1,IL2 indici wl MIN MAX in autos
C     (27,1),(27,2).......Delta_chi2 ,S1P2
C         0->exit   1->VALIN   2->SIMULA   3->IBRIDONE
C     (28,1) .............N per il calcolo di <f(x)>
C            (28,2).......mwl (N. tot wl)
C     (29,1),.............N sublayer discretizzazione strato inomogeneo
C            (29,2).......md (N. tot nk)
C     (30,1-5) ...........wl,n,k,ERRn,ERRk
C     (31,1-    4)...Procedura_j, Rigidit_j, VALini_j, VALfin_j
C     (32,1-    4)...   "           "          "         "
C     (33,1-    4)...   "           "          "         "
C     (34,1-    4)...   "           "          "         "
C     (35,1-    4)... IMED   S1P2U3  WL_nothc  CNK(1,1) in simula
C
C     (51-60,1)      puntatore su VNK -> nk layer #j=1,9
C         (51,2)     N. complessivo degli strati film & bulk (max 9)
C         (52,2)     multistrato simmetrico su 2a faccia ultimo strato:
C                     0 = no, 1 = si
C         (53,2)     puntatore strato incognito
C         (54,2)     tipo misure SF: speculare 0, emisferiche 1
C         (55,2)     chi2fin (last)
C         (56,2)     fredeg (degree of freedom)
C       (51-60,3)    tipologia strato: 1 = bulk (incoerente)
C                                      2 = film omogeneo
C                                      3 = film inomogeneo
C       ( 1-20,3)    INSTR(20)
C       (21-26,3)    DNK(6)
C        ( 1-14,4)   DATO(14)
C        (21-26,4)   Nmis
C            (31,5)  scala log k-plot (0->lineare 1->log)
C            (32,5)  0->fit_n 1->fit_nk 2->fit epsi1 and epsi2
C            (33,5)  /
C            (34,5)  N parametri del pannello di fit
C            (35,5)  N parametri abilitati al fit
C
C      (36,1-3) Tn:    Delta/error(rms) <Experimental>  <Simulation>
C      (37,1-3) Tp:    Delta/error(rms) <Experimental>  <Simulation>
C      (38,1-    4)... WLmin_solnk, WLmax_solnk,WLmis_mis,WLmax_mis enabled 
C        .....
C      (49,1-3) PSI-4: Delta/error(rms) <Experimental>  <Simulation>
C
C         (36-53,5)  PPM(17) 
C
C
C        DF(1) : n,k del sub = 0 (k-bord) , 1 (file)
C        DF(2) : misure spettrofoto.= -1 (no), 0 (k-bord) , 1 (file)
C        DF(3) : misure ellissome.  = -1 (no), 0 (k-bord) , 1 (file)
C
C        IUVIR : range spettrale misura = 1 UV-VIS-NIR
C                                         2 IR
C       
C
C   INSTR(20) vettore parametri strumentali trasportato da PAR:
C                         INSTR(i)=PAR(i,3)
C   INSTR(1) : a(deff)= efficienza di raccolta seconda faccia al 1ord
C   INSTR(2) : a(2*deff)= efficienza di raccolta seconda faccia al 2ord
C   INSTR(3) : DBK/BK= fluttuazione back correction SF
C   INSTR(4) : DRcal/Rcal=errore su R specchio di calibrazione
C   INSTR(5) : Errore di Lettura
C   INSTR(6) :   "    "    "     (in)dipendente da wl (0)1    
C   INSTR(7) : cte per calcolo contributo SUB al segnale PDS
C   INSTR(8) : /
C   INSTR(9) : /
C   INSTR(10): 0<-> V-W;
C              1<-> RIF05
C              2<-> RIF06 prima del 5/dic/94
C              3<-> RIF06 dopo  il  5/dic/94
C              4<-> RIF05 dopo  il 18/mar/96
C              5<-> RIF08
C              6<-> RIF08 dopo 17/gen/2018
C   INSTR(11): /
C   INSTR(12): /
C   INSTR(13): /
C   INSTR(14): /
C   INSTR(15): /
C   INSTR(16): /
C   INSTR(17): /
C   INSTR(18): /
C   INSTR(19): /
C   INSTR(20): /
C
C
C
C   DNK(6) incrementi ricerca soluzioni di autos
C          DNK(J)=PAR(20+J,3) J=1,2,..6
C
C   DATO(14) misure sperimetali caricate, trasportato da PAR:
C                         DATO(i)=PAR(i,4)
C       i=1: Tn
C       i=2: Tp                            
C       i=3: Rn
C       i=4: Rp
C       i=5: R1
C       i=6: Apds
C       i=7: DELTA_1
C       i=8: PSI_1
C      ....
C       i=14: PSI_4
C       #i=15: N. misure su cui operare
C       #i=16: N. angoli diversi su cui operare
C       #i=17: N. misura da graficare su plot di n VS wl
C
C       dato(i) = 0 -> misura assente
C                 1 -> misura caricata
C                 2 -> misura selezionata su cui operare
C
C
C
C
C
c **** PM(200,5) Multistrato generalizzato & Sellmeier quantistica ****
c        (i,1) = valore parametro
c        (i,2) = gestione fit: 0 -> costante
c                              j -> P(j) dove P parametro di fit
c        (j,3) = indirizzamento P: P(j) = PM(int(PM(j,3),1)
c        (i,4) = ERROR
c        (i,5) = global correlation
c
c   1- 9 spessori
c  11-19 Dn/<n>
c  21-29 curvatura di n
c  31-39 Dk/<k>
c  41-49 curvatura di k
c  51-59 roughness  
c  61-69 slope_Dn/<n> in 1/eV  
c  71-79 /
c  81-89 /
c  91-99 /
c 100  FitOp# considerata in Ibridone (1,2 ...,7)
c 101  fdisp#1
c 102  C1
c 103  E1
c 104  D1
c 105  K1
c etc etc fino a 20 oscillatori
C   
C
C
C   PPM(nparmax=17) puntatore PM trasportato da PAR:
C                       PPM(i)=par(35+i,5) i=1,17
C                       n_parametri_fit=par(35,5) 
C
C
C
C   MIS(17,201,2) matrice misure spettrofotometriche e nk materiali noti
C      (i,j,k) 
C       i=1: Tn     j=1,201 <-> Lambda   k=1 <-> valore
C       i=2: Tp                            2 <-> errore
C       i=3: Rn
C       i=4: Rp
C       i=5: R1
C       i=6: Apds
C       i=7: Lambda
C       i=8-15:  materiale #1-8   k=1 <-> n     
C                                   2 <-> k
C       i=16: n(k=1),k(k=2) ibridone
C       i=17: abilitazione dato per fit in Ibridone
C
C   ELI(8,201,2) matrice misure ellissometriche
C      (i,j,k) 
C       i=1: DELTA_1   j=1,201 <-> Lambda   k=1 <-> valore
C       i=2: PSI_1                             2 <-> errore
C      .. ..
C       i=7: DELTA_4
C       i=8: PSI_4
C
C
C   NANK(1,..,8) : nome "mate/aa999.9" nk-noti
C   NANK(9)      : nome "mate/aa999.9" nk-soluzioni
C   NANK(10)     : nome "mate/aa999  " media custom
C   NANK(11)     : nome "mate/aa999.9" misure SF
C   NANK(12,.,15): nome "mate/aa999.9" misure EL
C   NANK(16)     : nome "mate/aa999.9" progetto
C
C   CNK(J,K): matrice di gestione di VNK(J,H)
C       J=1    <-> VNK(1,1 o 2)    = n,k DI LAVORO
C       J=2..9 <-> VNK(3..8,1 o 2) = n,k NOTI
C       J=10   <-> VNK(10,1 o 2)   = n,k MEZZO_ingresso Tn,..,R1
C       J=11   <-> VNK(11,1 o 2)   = n,k MEZZO_ingresso Apds
C       J=12   <-> VNK(12,1 o 2)   = n,k MEZZO_ingresso PSI,DELTA #1
C       J=13   <-> VNK(13,1 o 2)   = n,k MEZZO_ingresso PSI,DELTA #2
C       J=14   <-> VNK(14,1 o 2)   = n,k MEZZO_ingresso PSI,DELTA #3
C       J=15   <-> VNK(15,1 o 2)   = n,k MEZZO_ingresso PSI,DELTA #
C       J=16   <-> VNK(9,1 o 2)    = n,k MEZZO_uscita   Tn,..,R1
C    CNK(J,1)=17      =>  n,k  INCOGNITI
C    CNK(J,1)=0       =>  n,k cte  n=CNK(J,2)  k=CNK(J,3)
C    CNK(J,1)=1..5    =>  n,k da  #.ft
C    CNK(J,1)=8..15   =>n,k da file.nk caricati in MIS(K,L,NK)
C
C
C    RXY(30,4): matrice di gestione della window (coordinate reali)
C       (J ,K)
C 
C       J= 1 <-> Tn
C       J= 2 <-> Tp
C       J= 3 <-> Rn
C       J= 4 <-> Rp
C       J= 5 <-> R1
C       J= 6 <-> Apds
C       J= 7 <-> Delta_1
C       J= 8 <-> PSI_1
C       J= 9 <-> Delta_2
C       J=10 <-> PSI_2
C       J=11 <-> Delta_3
C       J=12 <-> PSI_3
C       J=13 <-> Delta_4
C       J=14 <-> PSI_4
C       J=15 <-> FM, funzione di merito
C       J=16 <-> n
C       J=17 <-> k
C       J=18 <-> Aaria=1-Tn-Rn
C       J=19 <-> Derivata di n 
C       J=20 <-> Lambda (Angstrom)
C       J=21 <-> Teta (deg)
C       J=22 <-> Thickness (Angstrom)
C       J=23 <-> Tau,Rho,Rho1
C       J=24 <-> Size of graph (width & aspect)
C       J=25 <-> informazioni sul grafico
C       J=26 <-> eps1
C       J=27 <-> eps2
C
C          per 1=< J <=23 && 26,27 informazioni sulle quantita' fisiche
C              RXY(J,1)  = Xmin, estremo INF asse grafico
C              RXY(J,2)  = Xmax, estremo SUP asse grafico
C              RXY(J,3)  = Vmin, valore MIN quantita' J-esima
C              RXY(J,4)  = Vmax, valore MAX quantita' J-esima
C          per J=24 settaggi finestre grafiche
C          per J=25  informazioni sul grafico
C         int(RXY(25,1)) = indice quantita' graficata su l'asse X
C         int(RXY(25,2)) = indice quantita' graficata su l'asse Y   
C         int(RXY(25,3)) = 1 => fotone in WL
C         int(RXY(25,3)) = 2 => fotone in eV
C         int(RXY(25,4)) = 1 => step spline WL
C         int(RXY(25,4)) = 2 => step spline eV
C
C    ARSE(500,2) array di servizio 

      PROGRAM SEMAW
      INTEGER ixw(20)
      REAL MIS(17,201,2),NK(999,6),PM(200,5),PAR(60,5),PF(7,21),
     *     ELI(8,201,2),CNK(16,3),RXY(30,4),ARSE(500,2)
      CHARACTER NANK(16)*256,scmd*5,fctrl*255,r2*2,r*1,st76*76,user*255
      COMMON /VALORI/ MIS,ELI,PAR,PM,CNK,NK,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      WRITE(*,*)
      WRITE(*,*)'****************************************************'
      WRITE(*,*)'              Program kSEMAW' 
      write(*,*)
      write(*,*)'Spectro-Ellipsometric Measurement Analysis Workbench'
      WRITE(*,*)'  (spectrophotometric, ellipsometric and PDS)'
      WRITE(*,*)
      write(*,*)'       version 0.9.6 5th May 2021'
      write(*,*)'       Author: Marco Montecchi, ENEA (Italy) '
      write(*,*)'       email: marco.montecchi@enea.it '
      WRITE(*,*)'****************************************************'

***** inizializzazione matrici e vettori
      NK(1,1)=0. !nessuna soluzione n-k caricata
      DO I=1,17
       DO J=1,201
        MIS(I,J,1)=0.
        if(i.le.15) then
         MIS(I,J,2)=1.
        else
         MIS(I,J,2)=0.
        end if
       END DO
      END DO
      DO I=1,8
       DO J=1,201
        ELI(I,J,1)=0.
        ELI(I,J,2)=1.
       END DO
      END DO
      DO I=1,7
       do j=1,21
        PF(I,j)=0.
       end do
      END DO

**** filename
      DO I=1,15
       NANK(I)='mate/aa999.9'
      END DO
      NANK(16)='temp/defau.1'

*** inizializzazione grafica
      do i=1,20
       ixw(i)=0 !nessuna finestra grafica aperta
      end do

**** menu principale di SEMAW
      call getenv("USER", user)
      length=len_trim(user)
      write(*,*) 'FORTRAN: user= ',user(1:length)
c      fctrl='/home/'//user(1:length)//
c     +'/Workspace/Qt4-prove/ksemaw/ksemaw.ctrl'
      fctrl='/home/'//user(1:length)//
     +'/Workspace/qtSource/ksemaw/ksemaw.ctrl'
      write(*,*) 'FORTRAN: fctrl= ',fctrl
      DO WHILE(scmd.NE.'exit ')
       call sleep(1)
       scmd='pause'
       do while(scmd.eq.'pause'.or.scmd.eq.'done!') !in attesa di un comando
 1      call sleep(1)
        open(1,status='old',file=fctrl,ERR=1)
        read(1,*) scmd
        if(scmd.ne.'pause'.and.scmd.ne.'done!') then
         if(scmd.eq.'simul'.or.scmd.eq.'ibrid') read(1,*) r2
         if(scmd.eq.'autos') read(1,*) r
         if(scmd.eq.'iofnk') then
          read(1,*) is
          if(is.eq.1) read(1,'(a76)') st76
         end if
         rewind(1)
         write(1,*)'proce'
        end if
        close(1)
       end do
       write(*,*)
       write(*,*)'++++++++++++++++++++++ ksemaw +++++++++++++++++++++++'
       write(*,*)'New job:',scmd
       
       if(scmd.eq.'valin') then
        call valin !carica file misure, nk e parametri
        call progetto(1)
       end if
       
       if(scmd.eq.'simul') then !simula
        call progetto(2)
        call cnfr(r2)
        call progetto(1)
       end if
       
       if(scmd.eq.'autos') then !calcolo nk ad inseguimento
        call progetto(2)
        call autos(r)
        call progetto(1)
       end if
       
       if(scmd.eq.'curso') then !leggi posizione cursore
        call progetto(2)
        ifinestra=nint(par(12,1))
        call pgslct(ifinestra)
        cursore=pgcurs(xval,yval,r)
        if(nint(rxy(25,3)).eq.1) then
         par(12,2)=xval
         par(13,1)=12400./xval
        else
         par(12,2)=12400./xval
         par(13,1)=xval
        end if
        par(13,2)=yval
        call progetto(1)
       end if
       
       if(scmd.eq.'iofnk') then !carica/salva file nk
        call progetto(2)
        call IONK(is,NK,NANK,RXY,st76)
        call progetto(1)
       end if
       
       if(scmd.eq.'ibrid') then !IbridOne
        call progetto(2)
        call FITSELMQ(r2)
        par(32,5)=0. !reset to fit only n
        call progetto(1)
       end if
      
       write(*,*)'scrivo done!'
       open(1,status='old',file=fctrl) !segnala fine task
       write(1,'(a5)') 'done!'
       close(1)
      END DO
      CALL pgend
      STOP
      END
      
      
      
      SUBROUTINE IONK(ic,NK,NANK,RXY,st76)
******Salva/carica i file nk (file.nk)
c      ic =
c           1 -> salva
c           2 -> carica
      REAL NK(999,6),NKNEW(999,5),RXY(30,4)
      INTEGER H
      CHARACTER NANK(16)*256,R*1,ST76*76,tab*1
      tab=char(9)
      length=len_trim(nank(9))
      if(ic.eq.1) then
       write(*,*)'Save file.nk= '//nank(9)(1:length)
       ios=0
       OPEN(30,STATUS='UNKNOWN',FILE=nank(9)(1:length)//'.nk',
     +     IOSTAT=IOS,ERR=2)
2      IF(IOS.NE.0) 
     +   WRITE(*,*) nank(9)(1:length), 'NON PUO^ ESSERE APERTO'
       IF(IOS.EQ.0) THEN
c        write(*,'(a76)') st76
        WRITE(30,'(A76)') ST76
        N=nINT(NK(1,1))
        WRITE(30,'(I5)') N
        DO L=2,N+1
         WRITE(30,*) NK(L,1),tab,NK(L,2),tab,NK(L,3),tab,NK(L,4),tab,
     +               NK(L,5)
        END DO
        CLOSE(30)
       END IF
      ELSE
       write(*,*)'Load file.nk= '//nank(9)(1:length)
       IOS=0
       OPEN(30,FILE=nank(9)(1:length)//'.nk',STATUS='OLD',
     +  IOSTAT=IOS,ERR=3)
3      IF(IOS.EQ.0) THEN
        READ(30,'(A76)') ST76
c        WRITE(*,'(A76)') ST76
        READ(30,'(I5)') N
        WRITE(*,*)'N. (L,n,k) = ',N
        WRITE(*,*)'SPAZIO DISPONIBILE = ',nINT(998-NK(1,1))
        IF(R.EQ.'y'.AND.N.LE.(998-nINT(NK(1,1)))) THEN
         NOF=nint(NK(1,1))
        ELSE
         NOF=0
         RXY(16,3)=1000.
         RXY(16,4)=-1000.
         RXY(17,3)=1000.
         RXY(17,4)=-1000.
        END IF
        NK(1,1)=NOF+N
        DO L=NOF+2,N+NOF+1
         READ(30,*,IOSTAT=IOS2) NK(L,1),NK(L,2),NK(L,3),NK(L,4),NK(L,5)
         if(ios2.ne.0) then
          nk(L,4)=.0
          nk(L,5)=.0
         end if
         RXY(20,3)=MIN(RXY(20,3),NK(L,1))
         RXY(20,4)=MAX(RXY(20,4),NK(L,1))
         RXY(16,3)=MIN(RXY(16,3),NK(L,2)-NK(L,4))
         RXY(16,4)=MAX(RXY(16,4),NK(L,2)+NK(L,4))
         RXY(17,3)=MIN(RXY(17,3),NK(L,3)-NK(L,5))
         RXY(17,4)=MAX(RXY(17,4),NK(L,3)+NK(L,5))
        END DO
        CLOSE(30)
        DO I=2,nint(NK(1,1))+1
         WWM=1.E30
         DO H=2,nint(NK(1,1))+1
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
        JN=-1
        JK=-1
        ERRN=.0
        ERRK=.0
        i=2
        inew=2
        DO while(I.le.(nint(NK(1,1))+1))
         IF(nknew(I,4).LE.0.) nknew(I,4)=0.001*nknew(I,2)
         IF(nk(I,5).LE.0) THEN
          if(nknew(I,3).gt.0.) then
           nknew(I,5)=0.001*nknew(I,3)
          else
           nknew(I,5)=0.001*(rxy(17,4)-rxy(17,3))
          end if
         END IF
         NK(inew,1)=NKNEW(I,1)
         NK(inew,2)=NKNEW(I,2)
         NK(inew,3)=NKNEW(I,3)
         NK(inew,4)=NKNEW(I,4)
         NK(inew,5)=NKNEW(I,5)
         nk(inew,6)=1 !punto abilitato per il fit
         inew=inew+1
         i=i+1
        END DO
        nk(1,1)=inew-2
        N=nint(NK(1,1))
        WRITE(*,*)'RANGE : Lmin = ',NK(2,1),'  Lmax = ',NK(N+1,1)
        rxy(20,3)=MIN(rxy(20,3),nk(2,1))
        rxy(20,4)=MAX(rxy(20,4),nk(n+1,1))        
       ELSE
        WRITE(*,*)
     +  'IL FILE = ',nank(9)(1:length)//'.nk',' NON ESISTE !'
       END IF
      END IF
      RETURN
      END
      
      
      
      SUBROUTINE PROGETTO(ia)
      INTEGER ixw(20),length,irxy
      REAL MIS(17,201,2),NK(999,6),PM(200,5),PAR(60,5),PF(7,21),
     *     ELI(8,201,2),CNK(16,3),RXY(30,4),ARSE(500,2)
      CHARACTER NANK(16)*256,sp*260,st76*76,s4*4,st7*7
      LOGICAL answer
      COMMON /VALORI/ MIS,ELI,PAR,PM,CNK,NK,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
*** subroutine per il salvataggio/carica progetto
c     ia=0 salvataggio interattivo
c        1 salva su nank(16).Spj
c        2 carica 
c        3 salva progetto su temp/semaw.#jobtot.Spj
c        4 recupera progetto temp/semaw.#jobview.Spj
      write(*,*)'Progetto with ia=',ia
      length=len_trim(nank(16))
      sp=nank(16)(1:length)//'.Spj'
      irxy=30
      if(ia.le.1) then
       write(*,*)'save '//sp(1:length+4)
      else if(ia.eq.2) then
       write(*,*)'load '//sp(1:length+4)
      end if
      if(ia.eq.0) then
       OPEN(1,STATUS='UNKNOWN',FILE=sp(1:length+4),IOSTAT=IOS,ERR=3) 
 3     if(ios.ne.0) write(*,*)sp,'NON puo" essere aperto!!!'
      else if(ia.eq.1) then
       ios=0
      else if(ia.eq.3) then
       ios=0
       jobtot=nint(par(8,1))
       write(s4,'(I4.4)') jobtot
      end if
      if((ia.eq.0.or.ia.eq.1.or.ia.eq.3).and.ios.eq.0) then
       INQUIRE(1,OPENED=answer)
       if(answer.eqv..FALSE.) then
        if(ia.eq.0.or.ia.eq.1) then
         OPEN(1,STATUS='UNKNOWN',FILE=sp(1:length+4))
        else
         OPEN(1,STATUS='UNKNOWN',FILE='temp/semaw'//s4//'.Spj')
         st76='Job# = '//s4//' di  '//sp(1:length+4)
        end if
       end if
       if(ia.eq.1) then
        read(1,'(a7)') st7
        read(1,'(a76)') st76
       else
        write(1,'(a7)') 'iVspj=1'
        write(1,'(a76)') st76
       end if
       do i=1,15
        write(1,'(a256)') nank(i)
       end do
       do i=1,30
        write(1,*) (rxy(i,j), j=1,4)
       end do
       do i=1,15
        write(1,*) (cnk(i,j), j=1,3)
       end do
       do i=1,60
        write(1,*) (par(i,j), j=1,5)
       end do
       do i=1,200
        write(1,*) (pm(i,j), j=1,5)
       end do
       do i=1,21
        write(1,*) (pf(j,i), j=1,7)
       end do
       close(1)
       if(ia.eq.1) write(*,*)'..... salvato!'
      else if(ia.eq.2.or.ia.eq.4) then
       nriga=0
       if(ia.eq.2) then
        OPEN(1,STATUS='OLD',FILE=sp(1:length+4),IOSTAT=IOS,ERR=4)
       else if(ia.eq.4) then
        jobview=nint(par(8,2))
        write(s4,'(I4.4)') jobview
        write(*,*)'load temp/semaw'//s4//'.Spj'
        OPEN(1,STATUS='OLD',FILE='temp/semaw'//s4//'.Spj',
     +       IOSTAT=IOS,ERR=4)
       end if
 4     if(ios.eq.0) then
        read(1,'(a76)',err=5) st76
        if(st76(1:7).eq.'iVspj=1') then
         write(*,*) 'iVspj=1'
         irxy=30
         read(1,'(a76)',err=5) st76
        else
         write(*,*) 'old Spj Project'
         irxy=25
        end if
        nriga=nriga+1
        write(*,*) st76
        do i=1,15
         read(1,'(a256)',err=5) nank(i)
         length=len_trim(nank(i))
         nriga=nriga+1
c        write(*,*) 'nriga= ',nriga,' nank',i,' ',nank(i)(1:length)
        end do
        do i=1,irxy
         read(1,*,err=5) (rxy(i,j), j=1,4)
         nriga=nriga+1
c         write(*,*) nriga,'rxy',i,' ',(rxy(i,j), j=1,4)
        end do
        do i=1,15
         read(1,*,err=5) (cnk(i,j), j=1,3)
         nriga=nriga+1
         if(nint(cnk(i,1)).eq.16) then
          write(*,*)'*************************************************'
          write(*,*)'ATTENZIONE: opzione "16" nk <-> ibridone'
          write(*,*)'            eseguire subito ibridone!!!!'
          write(*,*)'*************************************************'
         end if
c        write(*,*) nriga,'cnk',i,' ',(cnk(i,j), j=1,3)
        end do
        cnk(16,1)=cnk(9,1) ! mezzo_OUT per SF
        cnk(16,2)=cnk(9,2)
        cnk(16,3)=cnk(9,3)
c        write(*,*)'cnk(16,1)= ',cnk(16,1)
        do i=1,60
         read(1,*,err=5) (par(i,j), j=1,5)
         nriga=nriga+1
c        write(*,*) nriga,'par',i,' ',(par(i,j), j=1,5)
        end do
c        ifu=nint(par(33,5))
c        if(ifu.lt.1.or.ifu.gt.3) par(33,5)=2.
        do i=1,200
         read(1,*,ERR=6) (pm(i,j), j=1,5)
c         pm(i,4)=.0
c         pm(i,5)=.0
 6       nriga=nriga+1
c         write(*,*) nriga,'pm',i,' ',(pm(i,j), j=1,5)
        end do
        do i=1,21
         read(1,*,err=5,end=5) (pf(j,i), j=1,7)
         nriga=nriga+1
c        write(*,*) nriga,'pf',i,' ',(pf(i,j), j=1,7)
        end do
        close(1)
        if(ia.eq.2) write(*,*)'..... letto!'
  5     continue
       else
        write(*,*)'   ATTENZIONE progetto inesistente!'
       end if
      end if
      if(ia.eq.4) par(8,1)=jobtot
      RETURN
      END



!       SUBROUTINE CODICI
!       WRITE(*,*)'LISTA DEI SUBSTRATI STANDART                          '
!       WRITE(*,*)'                                                      '
!       WRITE(*,*)'"v1":Vetrino Balzer      11/11/91   D=25mm    d=1mm   '
!       WRITE(*,*)'"v2":Vetro SIV (Rn=R1)   28/2/91    frammento d=5.9mm '
!       WRITE(*,*)'                                                      '
!       WRITE(*,*)'"q1":Quarzo Tetrasil                D=25mm    d=1mm   '
!       WRITE(*,*)'"q2":Quarzo Tetrasil     12/7/90    D=25.5mm  d=1mm   '
!       WRITE(*,*)'                                                      '
!       WRITE(*,*)'"s1":Silicio 1.02-2.5um  14/1/91    frammento d=1.082 '
!       WRITE(*,*)'"s2":Silicio 2.5-10.0um    "            "       "     '
!       WRITE(*,*)'                                                      '
!       WRITE(*,*)'                                                      '
!       WRITE(*,*)'"g1":Germanio 1.75-45um  22/7/91    ORIEL     d=3.03mm'
!       WRITE(*,*)'                                                      '
!       RETURN
!       END
! 
! 
! 
! 
!       SUBROUTINE SIGLE
!       WRITE(*,*)'    MATERIALI <---->  mate                            '
!       WRITE(*,*)'                                                      '
!       WRITE(*,*)'                                                      '
!       WRITE(*,*)' a-C     = a-c_               SiO2    = sio2          '
!       WRITE(*,*)' a-C:H   = a-ch               Ta2O5   = tao5          '
!       WRITE(*,*)' Ag      = ag__               TiO2    = tio2          '
!       WRITE(*,*)' BaF2    = baf2               Vetro   = vetr          '
!       WRITE(*,*)' CeO2    = ceo2               ZnO     = zno_          '
!       WRITE(*,*)' Cr      = cr__               ZnS     = zns_          '
!       WRITE(*,*)' Cu      = cu__               ZnSe    = znse          '
!       WRITE(*,*)' GaAs    = gaas               ZrO2    = zro2          '
!       WRITE(*,*)' Ge      = ge__               KCl     = kcl_          '
!       WRITE(*,*)' HfO2    = hfo2               Y2O35   = y2o3          '
!       WRITE(*,*)' MgF2    = mgf2                                       '
!       WRITE(*,*)' Ni      = ni__                                       '
!       WRITE(*,*)' Quarzo  = quar                                       '
!       WRITE(*,*)' Silicio = sili                                       '
!       WRITE(*,*)'                                                      '
!       RETURN
!       END
! 



      SUBROUTINE VALIN
      REAL M(17,201,2),PM(200,5),ELI(8,201,2),PAR(60,5),
     *     CNK(16,3),PF(7,21),RXY(30,4),SOL(999,6),ARSE(500,2)
      INTEGER DATO(14),ixw(20)
      CHARACTER NANK(16)*256,DIS(6)*3,ST2*2,ST76*76,EST(7)*6,
     *          ST1*1
      COMMON /VALORI/ M,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
 
***** inizializzazione matrice misure e soluzioni temp
      SOL(1,1)=0. !nessuna soluzione n-k caricata
      DO I=1,17
       DO J=1,201
        M(I,J,1)=0.
        if(i.le.15) then
         M(I,J,2)=1.
        else
         M(I,J,2)=0.
        end if
       END DO
      END DO
      DO I=1,8
       DO J=1,201
        ELI(I,J,1)=0.
        ELI(I,J,2)=1.
       END DO
      END DO
      call PROGETTO(2)
      iuvir=nint(par(24,1))
      IF(IUVIR.EQ.1) ST2='.v'
      IF(IUVIR.EQ.2) ST2='.i'
      DIS(1)='.tn'
      DIS(2)='.tp'
      DIS(3)='.rn'
      DIS(4)='.rp'
      DIS(5)='.r1'
      DIS(6)='.an'
      do i=1,14
       dato(i)=nint(par(i,4))
       if(dato(i).ge.1.and.i.le.6) then
        WRITE(ST1,'(I1)') nint(par(20+i,4))
        EST(I)=ST2//ST1//DIS(I)
       end if
      end do
      if(NANK(9)(1:12).ne.'mate/aa999.9') then
       length=len_trim(nank(9))
       write(*,*)'Caricamento soluzioni nk dal file ',nank(9)(1:length)
       call IONK(2,SOL,NANK,RXY,st76)
      end if 
      CALL SPADA(EST,IUVIR)
      par(38,3)=par(4,1)
      par(38,4)=par(4,2)
      par(21,1)=201
      par(28,2)=201
      RETURN
      END



      
      SUBROUTINE SPADA(EST,IUVIR)
      INTEGER DATO(14),STEP,DF(3),H,H0,HRID,OK,IXW(20),ios,MR
      REAL M(17,201,2),E(8,201,2),PAR(60,5),PM(200,5),
     *     CNK(16,3),SOL(999,6),PF(7,21),RXY(30,4),
     +     ARSE(500,2)
      REAL X(10000),Y(10000),Z(10000),LMIN,LMAX,LAM,
     *  EPR(2,201,2),INSTR(20)
      CHARACTER ST6*6,ST8*8,ST76*76,R*1,ST10*10,NANK(16)*256,ST40*40,
     *   EST(7)*6,ST4*4,Q*1,SAL*1,SPA*1,st15*15,
     +   s76*76,st9*9,st12*12,ER*1,st23*23,st256*256,user*255
      COMMON /VALORI/ M,E,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      df(1)=nint(par(22,2))
      df(2)=nint(par(23,1))
      df(3)=nint(par(23,2))
      do i=1,20
       instr(i)=par(i,3)
      end do
      if(instr(6).lt..5) then
       ER='y'
      else
       ER=' '
      end if
      do i=1,14
       dato(i)=nint(par(i,4))
      end do
      DMIN=.0
      DMAX=1.e6
      SPA='L'
      if(nint(rxy(25,4)).eq.2) SPA='E'
      SAL='n'
      MR=nint(par(22,1))
      nrif=nint(par(10,3))
      iexp=0 
      
*** controllo file nk
      DO J=1,8
       length=len_trim(nank(j))
       IF(NANK(J)(1:length).NE.'mate/aa999.9') THEN
        write(*,*)'load_nk '//NANK(J)(1:length)//'.nk'
        OPEN(10,FILE=NANK(J)(1:length)//'.nk',STATUS='OLD',ERR=1,
     +       IOSTAT=IOS)
        READ(10,'(A76)') ST76
        READ(10,'(I5)') NDATI
        if(ndati.eq.0) goto 1
        DO I=1,NDATI
         READ(10,*) X(I),Y(I),Z(I)
        END DO
        CLOSE(10)
        DMIN=MAX(X(1),DMIN)
        DMAX=MIN(X(NDATI),DMAX)
1       if(ios.ne.0) 
     +   write(*,*)'ERRORE aprendo il file '//nank(j)(1:length)
       END IF
      END DO
      
*** controllo file SF
      IF(DF(2).EQ.1) THEN
       DO I=1,6
        IF(DATO(I).ge.1) THEN
         length=len_trim(nank(11))
         write(*,*)'load_SF '//nank(11)(1:length)//EST(I)
         OPEN(10,FILE=nank(11)(1:length)//EST(I),STATUS='OLD',
     +        ERR=204,IOSTAT=IOS)
         READ(10,'(A76)') ST76
         READ(10,'(A12)') ST12
         if(st76(1:5).eq.'PE UV') then
c trattasi di file Perkin Elmer L900 - L950
          do ir=1,11
           read(10,*)
          end do
          READ(10,'(A23)') st23
          if(st23.eq.'PerkinElmer UV WinLab 5') then
           ilinrim=65
          else
           ilinrim=69
          end if
          do ir=1,ilinrim
           read(10,*)
          end do
          read(10,*) wmax
          read(10,*) ridelta
          idelta=nint(ridelta)
          read(10,*) ndati
          wmin=wmax+idelta*(ndati-1)
          wmax=wmax*10.
          wmin=wmin*10.
         else if(st12.eq.'#####SCALED '.or.st12.eq.'#####SCALEDA'.or.
     *      st12.eq.'#####SCALED%') then
c         trattasi di file scalato
          read(10,*) wini
          read(10,*) wfin
          if(wini.lt.wfin) then
           wmin=wini
           wmax=wfin
          else
           wmin=wfin
           wmax=wini
          end if
          wmax=wmax*10.
          wmin=wmin*10.         
         elseif(ST12.EQ.'##XYUNITS= W') THEN
c         trattasi di file l/9 o ftir vecchio tipo
          READ(10,'(a10,i6,a11)') ST10,IMAX,ST76
          READ(10,'(a10,i6,a11)') ST10,IMIN,ST76
          LMIN=real(IMIN)
          LMAX=real(IMAX)
          IF(IUVIR.EQ.1) THEN
           WMIN=LMIN*10.
           WMAX=LMAX*10.
          ELSE
           WMIN=1.E8/LMAX
           WMAX=1.E8/LMIN
          END IF
          
         ELSE
          READ(10,'(A15)') ST15
          
          IF(ST15.EQ.'##DATA TYPE= UL') THEN
c          trattasi file l/9 nuovo tipo
           DO L=1,8
            READ(10,'(A76)') ST76
           END DO
           READ(10,*) ST9,LMIN
           READ(10,*) ST8,LMAX
           WMIN=LMIN*10.
           WMAX=LMAX*10.
           
          ELSE IF(ST15.EQ.'##DATA TYPE= UV') THEN
c          trattasi file l/19
           DO L=1,11
            READ(10,'(A76)') ST76
           END DO
           READ(10,*) ST9,LMAX
           READ(10,*) ST8,LMIN
           WMIN=LMIN*10.
           WMAX=LMAX*10.
           
          ELSE IF(ST15.EQ.'##DATA TYPE= IN') THEN
c          trattasi file ftir nuovo tipo
           DO L=1,14
            READ(10,'(A76)') ST76
           END DO
           READ(10,*) ST9,LMIN
           READ(10,*) ST8,LMAX
           WMIN=1.E8/LMIN
           WMAX=1.E8/LMAX
           
          END IF
         END IF
         DMIN=MAX(DMIN,WMIN)
         DMAX=MIN(DMAX,WMAX)
         CLOSE(10)
204      if(ios.ne.0) write(*,*) 'ERRORE aprendo il file '
     +        //ST10//EST(I)
        END IF
       END DO
      END IF
      
*** controllo file ELI
      IF(DF(3).EQ.1) THEN
       DO I=1,4
        IF(NANK(11+I)(1:12).NE.'mate/aa999.9') THEN
         write(*,*)'load_EL '//NANK(11+I)//'.el'
         OPEN(10,FILE=NANK(11+I)//'.el',STATUS='OLD')
         READ(10,'(A40)') ST40
         READ(10,'(I3)') NANG
         DO J=1,NANG
          READ(10,'(F6.2,1X,I4,1X,F7.1,1X,F7.1)') TE,NDAT,LMIN,LMAX
          IF(ABS(TE-PAR(13+J,1)).LT..001) THEN
           WMIN=LMIN*10.
           WMAX=LMAX*10.
          END IF
         END DO
         DMIN=MAX(DMIN,WMIN)
         DMAX=MIN(DMAX,WMAX)
         CLOSE(10)
        END IF
       END DO
      END IF
      Q='?'
c      rxy(20,1)=dmin
c      rxy(20,2)=dmax
      dmin=par(4,1)
      dmax=par(4,2)
      write(*,*)'******************** SPADA *************************'
      WRITE(*,*)'r) WLmin = ',dmin,'  WLmax = ',dmax
      IF(DATO(3).ge.1.OR.DATO(4).ge.1.OR.DATO(5).ge.1) THEN
       IF(MR.EQ.0) then
        WRITE(*,*)'M) Rn,Rp,R1 assolute'
       else
        call getenv("USER", user)
        length=len_trim(user)
        OPEN(11,FILE='/home/'//user(1:length)//
     +'/Workspace/qtSource/ksemaw/referenceMirrors.txt',STATUS='OLD')
c     +'/Workspace/Qt4-prove/ksemaw/referenceMirrors.txt',STATUS='OLD')
        do i=1,nrif
         read(11,'(a76)') st76
         read(11,'(a256)') st256
         length=len_trim(st256)
        end do
        close(11)
        WRITE(*,*)'M) Rn,Rp,R1 da * Rrif con'
        write(*,*)'Rif = ',st76
        write(*,*) st256(1:length)
c         if(nrif.eq.1) write(*,*)'Rrif = RIF05'
c         if(nrif.eq.2) write(*,*)'Rrif = Rif06 before 5 December 1994'
c         if(nrif.eq.3) write(*,*)'Rrif = Rif06 after 5 Decembre 1994'
c         if(nrif.eq.4) write(*,*)'Rrif = Rif05 after 18 March 1996'
c         if(nrif.eq.5) write(*,*)'Rrif = Rif08 since 13 December 2011'
c         if(nrif.eq.6) write(*,*)'Rrif = Rif08 after 17 January 2018'
c         if(nrif.eq.7) write(*,*)'Rrif = RifMir PV ENEA'
       end if
      END IF
      WRITE(*,*)'e) Parametri per il calcolo dell"errore SF:'
      WRITE(*,*)'        DBK/BK     = ',INSTR(3)
      WRITE(*,*)'        DRrif/Rrif = ',INSTR(4)
      IF(ER.NE.'y') THEN
       WRITE(*,*)'        DLettura   = ',INSTR(5),' dipendente da WL'
      ELSE
       WRITE(*,*)'        DLettura   = ',INSTR(5),' indipendente da WL'
      END IF
      IF(SPA.EQ.'E')THEN
       WRITE(*,*)'s) Equispaziatura 201 punti in ENERGIA'
      ELSE
       WRITE(*,*)'s) Equispaziatura 201 punti in WL'
      END IF
      write(*,*)'****************************************************'

      STEP=1
*** costruzione WL e fluttuazione linea di base
      DO L=1,201
       IF(SPA.EQ.'L') THEN
        M(7,L,1)=DMIN+(L-1)/200.*(DMAX-DMIN)
       ELSE
        M(7,L,1)=1./(1./DMIN+(L-1)/200.*(1./DMAX-1./DMIN))
       END IF
       M(6,L,1)=1.!riflettanza di correzione
       M(7,L,2)=INSTR(5)
       M(17,L,1)=1. !dato abilitato per Ibridone
       IF(ER.NE.'y') THEN
        IF(IUVIR.EQ.1)THEN
         IF(M(7,L,1).GT.1849..AND.M(7,L,1).LT.8608.1) THEN
          M(7,L,2)=.0005+.838951/(1+((M(7,L,1)-1835.)/7.)**2)+
     *                .0026/(1+((M(7,L,1)-8700.)/250.)**4)
         END IF
         IF(M(7,L,1).GT.31000..AND.M(7,L,1).LT.32001.) THEN
          M(7,L,2)=-4.66E-04+.00476/(1+((M(7,L,1)-35000)/2000)**2)
         END IF
        END IF
        IF(M(7,L,2).LT.INSTR(5)) M(7,L,2)=INSTR(5)
       END IF
      END DO
*** caricamento file nk
      DO I=1,8
       length=len_trim(nank(i))
       IF(NANK(I)(1:length).NE.'mate/aa999.9') THEN
        OPEN(10,FILE=NANK(I)(1:length)//'.nk',STATUS='OLD',ERR=10,
     +       IOSTAT=ios)
        READ(10,'(A76)') ST76
        READ(10,'(I5)') NDATI
        DO J=1,NDATI
         READ(10,*) X(J),Y(J),Z(J)
        END DO
        CLOSE(10)
        N=1
        CALL CONVER(X,Y,NDATI,N,STEP,7+I,1)
        CALL CONVER(X,Z,NDATI,N,STEP,7+I,2)
10      if(ios.ne.0) write(*,*)'ERRORE aprendo il file '//
     +    nank(i)(1:length)//'.nk'
       END IF
      END DO
*** caricamento file SF
      IF(DF(2).EQ.1) THEN
       DO I=0,6
        IF(I.EQ.0) then
         if(MR.ne.0) THEN
           length=len_trim(st256)
           write(*,*) st256(1:length)
           OPEN(10,FILE=st256(1:length),STATUS='OLD',IOSTAT=ios,ERR=200)
c           if(nrif.ge.1.and.nrif.le.4) then
c            OPEN(10,FILE='al__/ri006.v1.rn',STATUS='OLD',IOSTAT=ios,
c      +         ERR=200)
c            write(*,*)'Caricamento file al__/ri006.v1.rn'
c           else if(nrif.eq.5) then
c            OPEN(10,FILE='al__/ri008.v2.rn',STATUS='OLD',IOSTAT=ios,
c      +        ERR=200)
c            write(*,*)'Caricamento file al__/ri008.v2.rn'
c           else if(nrif.eq.6) then
c            OPEN(10,FILE='al__/ri008.v3.rn',STATUS='OLD',IOSTAT=ios,
c      +        ERR=200)
c            write(*,*)'Caricamento file al__/ri008.v3.rn'
c           else if(nrif.eq.7) then
c            OPEN(10,FILE='al__/rifPVenea.v1.rn',STATUS='OLD',IOSTAT=ios,
c      +        ERR=200)
c            write(*,*)'Caricamento file al__/rifPVenea.v1.rn'
c           end if
200       if(ios.ne.0) write(*,*)'ERRORE aprendo il file '//
     +                 'R_reference'
         else
          goto 202
         end if
        ELSE if(dato(i).ge.1) then
         length=len_trim(NANK(11))
         write(*,*)'Caricamento file ',nank(11)(1:length)//EST(I)
         OPEN(10,FILE=nank(11)(1:length)//EST(I),STATUS='OLD',
     +    IOSTAT=ios,ERR=201)
201      if(ios.ne.0) write(*,*)'ERRORE aprendo il file '//
     +                 nank(11)(1:length)//EST(I)
        ELSE
         goto 202
        END IF
         if(ios.ne.0) goto 202
         READ (10,'(76a)') st76
         READ(10,'(A12)') ST12
         write(*,*) 'ST12=',ST12
         if(st76(1:5).eq.'PE UV') then
*** file Perkin Elmer L900 - L950
          div=100.
          do ir=1,11
           read(10,*)
          end do
          READ(10,'(A23)') st23
          if(st23.eq.'PerkinElmer UV WinLab 5') then
           ilinrim=66
          else
           ilinrim=70
          end if
          do ir=1,ilinrim
           read(10,*)
          end do
          read(10,*) rstep
          step=nint(rstep)
          read(10,*) ndati
          n=ndati
          do ir=1,4
           read(10,*)
          end do
          do j=1,ndati
           read(10,*) x(j),y(j)
          end do 
         else if(st12.eq.'#####SCALED '.or.st12.eq.'#####SCALEDA'.or.
     *      st12.eq.'#####SCALED%') then
c   trattasi di file scalato
          read(10,*) xin
          read(10,*) xfi
          read(10,*) ndati
c          write(*,*) xin,xfi,ndati
          if(xin.lt.xfi) then
           step=1
           n=1
          else
           step=-1
           n=ndati
          end if
          if(st12.eq.'#####SCALED ') then
3          write(*,*)'N. di colonna da considerare (2,..,6)?'
           read(*,*,ERR=3) ncolo
           do j=1,ndati
            if(ncolo.eq.2) read(10,*,ERR=3) x(j),y(j)
            if(ncolo.eq.3) read(10,*,ERR=3) x(j),zz,y(j)
            if(ncolo.eq.4) read(10,*,ERR=3) x(j),zz,zz,y(j)
            if(ncolo.eq.5) read(10,*,ERR=3) x(j),zz,zz,zz,y(j)
            if(ncolo.eq.6) read(10,*,ERR=3) x(j),zz,zz,zz,zz,y(j)
           end do
           write(*,*)'Unita" assolute (y) o % (ret)'
           read(*,'(a)')r
           if(r.eq.'y') then
            div=1
           else
            div=100.
           end if
           
          elseif(st12.eq.'#####SCALEDA') then
           do j=1,ndati
            read(10,*) x(j),y(j)
c            write(*,*) x(j),y(j)
           end do
           div=1.
           
          elseif(st12.eq.'#####SCALED%') then
           do j=1,ndati
            read(10,*) x(j),y(j)
           end do
           div=100.
          end if
          
          
         elseIF(ST12.EQ.'##XYUNITS= W') THEN
c trattasi di file l/9 o ftir vecchio tipo
          READ(10,'(A10,A6,4X,A8)') ST76,ST6,ST8
          READ(10,'(A10,A6,4X,A8)') ST76,ST6,ST8
          READ(10,*) ST12,dlam,div
          DIV=div*100.
          READ(10,*) ST10,ndati
          READ(10,'(A20)') ST76
          NRIGHE=INT(NDATI/14.)
          DE=14.*((NDATI/14.)-INT(NDATI/14.))
          IF(DE.GT..5) NRIGHE=NRIGHE+1
          DO J=1,NRIGHE
           K=(J-1)*14
       READ(10,'(F6.1,14F5.0)') X1,Y(K+1),Y(K+2),Y(K+3),Y(K+4),Y(K+5),
     *                           Y(K+6),Y(K+7),Y(K+8),Y(K+9),Y(K+10),
     *                           Y(K+11),Y(K+12),Y(K+13),Y(K+14)
           DO L=1,14
            X(K+L)=X1+(L-1)*DLAM
           END DO
          END DO
          IF(IUVIR.EQ.1) THEN
           N=NDATI
           STEP=-1
          ELSE
           N=1
           STEP=1
          END IF
          
         ELSE
          READ(10,'(A15)') ST15
          
          IF(ST15.EQ.'##DATA TYPE= UL') THEN
c trattasi di file l/9 nuovo tipo
           DO J=1,6
            READ(10,'(A76)') ST76
           END DO
           READ(10,*) ST10,xf
           READ(10,*) ST10,yf
           READ(10,*) ST9,wmin
           READ(10,*) ST8,wmax
           READ(10,*) ST10,ndati
           READ(10,'(A76)') ST76
           READ(10,'(A76)') ST76
           DLAM=(WMAX-WMIN)/(NDATI-1)
           DIV=1./YF
           NRIGHE=INT(NDATI/10.)
           DE=10.*((NDATI/10.)-INT(NDATI/10.))
           NP=nint(DE)
           DO J=1,NRIGHE
            K=(J-1)*10
            READ(10,*) X1,Y(K+1),Y(K+2),Y(K+3),Y(K+4),
     *          Y(K+5),Y(K+6),Y(K+7),Y(K+8),Y(K+9),Y(K+10)
            DO L=1,10
             X(K+L)=X1*XF+(L-1)*DLAM
            END DO
           END DO
           K=(J-1)*10
           IF(NP.EQ.1) READ(10,*) X1,Y(K+1)
           IF(NP.EQ.2) READ(10,*) X1,Y(K+1),Y(K+2)
           IF(NP.EQ.3) READ(10,*) X1,Y(K+1),Y(K+2),Y(K+3)
           IF(NP.EQ.4) READ(10,*) X1,Y(K+1),Y(K+2),Y(K+3),Y(K+4)
           IF(NP.EQ.5) READ(10,*) X1,Y(K+1),Y(K+2),Y(K+3),Y(K+4),Y(K+5)
           IF(NP.EQ.6) READ(10,*) X1,Y(K+1),Y(K+2),Y(K+3),Y(K+4),Y(K+5),
     *                              Y(K+6)
           IF(NP.EQ.7) READ(10,*) X1,Y(K+1),Y(K+2),Y(K+3),Y(K+4),Y(K+5),
     *                              Y(K+6),Y(K+7)
           IF(NP.EQ.8) READ(10,*) X1,Y(K+1),Y(K+2),Y(K+3),Y(K+4),Y(K+5),
     *                              Y(K+6),Y(K+7),Y(K+8)
           IF(NP.EQ.9) READ(10,*) X1,Y(K+1),Y(K+2),Y(K+3),Y(K+4),Y(K+5),
     *                              Y(K+6),Y(K+7),Y(K+8),Y(K+9)
           DO L=1,NP
            X(K+L)=X1*XF+(L-1)*DLAM
           END DO
           N=1
           STEP=1
           
          elseIF(ST15.EQ.'##DATA TYPE= UV') THEN
c trattasi di file l/19
           DO J=1,9
            READ(10,'(A76)') ST76
           END DO
           READ(10,'(A10,F5.5)') ST10,xf
           READ(10,'(A10,F23.20)') ST10,yf
           READ(10,*) ST9,wmax
           READ(10,*) ST8,wmin
           READ(10,*) ST10,ndati
           READ(10,'(A76)') ST76
           READ(10,'(A76)') ST76
           READ(10,'(A76)') ST76
           READ(10,'(A76)') ST76
           DLAM=(WMAX-WMIN)/(NDATI-1)
           DIV=100./YF
           NRIGHE=INT(NDATI/5.)
           DE=5.*((NDATI/5.)-INT(NDATI/5.))
           NP=nint(DE)
           DO J=1,NRIGHE
            K=(J-1)*5
            READ(10,*) X1,Y(K+1),Y(K+2),Y(K+3),Y(K+4),Y(K+5)
            DO L=1,5
             X(K+L)=X1*XF-(L-1)*DLAM
            END DO
           END DO
           K=(J-1)*5
           IF(NP.EQ.1) READ(10,*) X1,Y(K+1)
           IF(NP.EQ.2) READ(10,*) X1,Y(K+1),Y(K+2)
           IF(NP.EQ.3) READ(10,*) X1,Y(K+1),Y(K+2),Y(K+3)
           IF(NP.EQ.4) READ(10,*) X1,Y(K+1),Y(K+2),Y(K+3),Y(K+4)
           DO L=1,NP
            X(K+L)=X1*XF-(L-1)*DLAM
           END DO           
           N=NDATI
           STEP=-1
           
          ELSE IF(ST15.EQ.'##DATA TYPE= IN') THEN
c trattasi di file IR nuovo tipo
           DO J=1,9
            READ(10,'(A76)') ST76
           END DO
           READ(10,*) ST9,DLAM
           READ(10,*)
           READ(10,*)
           READ(10,*) ST10,XF
           READ(10,*) ST10,yf
           READ(10,*) ST9,WMIN
           READ(10,*) ST8,WMAX
           READ(10,*)ST10,ndati
           DO J=1,4
            READ(10,'(A76)') ST76
           END DO
           DIV=1./YF
           K=0
           DO WHILE(K.LT.NDATI)
            K=K+1
            read(10,'(a76)') st76
            if(st76(1:1).ne.' ') then
c FTIR1760            
              H=1
              DO WHILE(st76(H:H).NE.' ')
               H=H+1
              END DO
              read(s76(1:h),*) x1
              X(K)=X1*XF
              DO WHILE(st76(H:H).EQ.' ')
               H=H+1
              END DO
              H0=H
              OK=1
              DO WHILE(OK.EQ.1)
               H=2
               DO WHILE((H+H0).LT.76.AND.st76(H+H0:H+H0).NE.'+'.AND.
     *            st76(H+H0:H+H0).NE.'-'.AND.st76(H+H0:H+H0).NE.' ')
                H=H+1
               END DO
               IF(st76(H+H0:H+H0).EQ.'+'.OR.st76(H+H0:H+H0).EQ.'-'.OR.
     *           st76(H+H0:H+H0).EQ.' ') THEN
                H=H-1
                HRID=1
               ELSE
                HRID=0
               END IF
               read(s76(h0:h0+h),*) y(k)
               IF(HRID.EQ.1) H=H+1
               IF(st76(H0+H:H0+H).EQ.'+'.OR.st76(H0+H:H0+H).EQ.'-') THEN
                H0=H0+H
                K=K+1
                X(K)=X(K-1)+DLAM
                OK=1
               ELSE
                OK=0
               END IF
              END DO
            else
c FTIR SPECTRUM GS
             if(k.eq.1) div=div*100.
             read(st76,*) x1,y(k),y(k+1),y(k+2),y(k+3),y(k+4)
             do iiii=1,5
              x(k+iiii-1)=x1*xf+(iiii-1)*dlam
             end do
             k=k+4
            end if
           END DO
           N=1
           STEP=1
           
          END IF
         END IF
         CLOSE(10)
         DO L=1,NDATI
          IF(IUVIR.EQ.2) X(L)=1./X(L)*1.E7
          X(L)=X(L)*10.
          Y(L)=Y(L)/DIV
         END DO
         IF(I.EQ.0) THEN
          CALL CONVER(X,Y,NDATI,N,STEP,6,IUVIR)
          rxy(6,3)=0.
          rxy(6,4)=0.
         ELSE
          CALL CONVER(X,Y,NDATI,N,STEP,I,IUVIR)
         END IF
202      continue
       END DO
      END IF
      IF(DF(3).EQ.1) THEN
       DO J=1,4
        IF(NANK(11+J)(1:12).NE.'mate/aa999.9') THEN
         length=len_trim(NANK(11+J))
         OPEN(30,FILE=NANK(11+J)(1:length)//'.el',STATUS='OLD')
         write(*,*)'Caricamento file ',NANK(11+J)(1:length)//'.el'
         READ(30,'(A40)') ST40
         write(*,*) ST40
c        READ(30,'(I3,1X,I4,3x,a4)') NANG,NDATI,st4
         read(30,*)  NANG,NDATI,st4
         write(*,*)'Nang= ',NANG,' NdatiEL= ',NDATI,' st4= ',st4
         DO I=1,NANG
          READ(30,'(A76)') ST76
         END DO
C        OFFSET=.0
C        WRITE(*,*)'Mis. EL ',NANK(11+J),' : offset in WL ?'
C        READ(*,'(A)') R
C        IF(R.EQ.'y') THEN
C         WRITE(*,*)'Offset = ?'
C         READ(*,*) OFFSET
C        END IF
         NDEL=0
         iformat=0
         if(st4.eq.'tldp') then
          iformat=1
         elseif(st4.eq.'etpd') then
          iformat=2
         elseif(st4.eq.'elia') then
          iformat=3
         end if
         if(iformat.eq.0) then
          write(*,*)'ATTENZIONE. imposta formato file ELI!!!'
        end if
         DO I=1,NDATI
          if(iformat.eq.1) then
c           READ(30,'(F6.2,1X,F7.1,1X,F7.3,1X,F6.3,1X,F5.3,1X,F5.3)')
            read(30,*)
     *      TE,LAM,DE,PS,DDE,DPS
          elseif(iformat.eq.2) then
           read(30,*) lam,te,ps,de,dps,dde
           lam=1240./lam
          elseif(iformat.eq.3) then
           read(30,*) lam,te,ps,de,dps,dde
          end if
c          write(*,*) TE,LAM,DE,PS,DDE,DPS
          IF(ABS(TE-PAR(13+J,1)).LT..001) THEN
           NDEL=NDEL+1
           X(NDEL)=LAM*10.
           EPR(1,NDEL,1)=DE
           EPR(1,NDEL,2)=DDE
           EPR(2,NDEL,1)=PS
           EPR(2,NDEL,2)=DPS
          END IF
         END DO
         CLOSE(30)
         write(*,*)'Theta= ',PAR(13+J,1),' NdatEl= ',NDEL
         DO I=1,2
          DO II=1,2
           DO L=1,NDEL
            Y(L)=EPR(I,L,II)
           END DO
           CALL CONVER(X,Y,NDEL,1,1,20+2*(J-1)+I,II)
          END DO
         END DO
        END IF
       END DO
      END IF

      RETURN
      END




      SUBROUTINE CONVER(X,Y,NDATI,N1,STEP,I,IUVIR)
      INTEGER STEP,H,ALARM,NS(1000),IXW(20)
      REAL M(17,201,2),E(8,201,2),PAR(60,5),PM(200,5),
     *     CNK(16,3),SOL(999,6),PF(7,21),RXY(30,4),
     +     ARSE(500,2)
      REAL X(10000),Y(10000),XS(500),YS(500),INSTR(20)
!      REAL P(nv),FVEC(nv),WA(nv),fjac(nm,nm) ! sono passati a SMINPACK
      INTEGER, allocatable :: IWA(:)
      REAL, allocatable :: P(:)
      REAL, allocatable :: FVEC(:)
      REAL, allocatable :: WA(:)
      REAL, allocatable :: fjac(:,:) 
      REAL TOL
      DOUBLE PRECISION XD1,XD2,XD3,YD1,YD2,YD3,DETD,AD,BD,CD
      CHARACTER R*1,NANK(16)*256
      EXTERNAL FPAR
      COMMON /VALORI/ M,E,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      iw=0 !verbosit (0 silente, 1 verboso)
      ALARM=0
      ymin=1.e+36
      ymax=1.e-36
      do j=1,20
       instr(j)=par(j,3)
      end do
      
      DO H=1,201
       WL=M(7,H,1)
       
*** fattore correttivo misure SF
C   INSTR(10): 0<-> V-W;
C              1<-> RIF05
C              2<-> RIF06 prima del 5/dic/94
C              3<-> RIF06 dopo  il  5/dic/94
C              4<-> RIF05 dopo  il 18/mar/96
C              5<-> RIF08
C              6<-> RIF08 dopo 17/gen/2018
       FACO=1.
       IF((I.EQ.3.OR.I.EQ.4.OR.I.EQ.5).AND.IUVIR.EQ.1) THEN
c        IF(nint(INSTR(10)).EQ.1) THEN
c         IF(WL.GT.2400..AND.WL.LT.4250.) FACO=-6.33E-6*WL+1.0245
c         IF(WL.GE.4250..AND.WL.LT.8000.) FACO=.9976
c         IF(WL.GE.8000..AND.WL.LT.8700.) FACO=.987+
c     *                                  .011/(1+(WL/1000.-7.6)**10.)
c         IF(WL.GE.8700.) FACO=.987
c        END IF
c        IF(nint(INSTR(10)).EQ.2) THEN
c         IF(WL.GE.1900..AND.WL.LE.8250)FACO=.9994+.04191/(1.+((WL-1110.)
c     *                                                      /1670.)**4.)
c         IF(WL.GT.8250.) FACO=.97197+.02796/(1.+((WL-8250.)/350.)**10.)
c        END IF
c        IF(nint(INSTR(10)).EQ.3) THEN
c         IF(WL.GT.8100.) FACO=.975+.025/(1.+((WL-8100.)/550.)**10.)
c        END IF
c        IF(nint(INSTR(10)).EQ.4) THEN
c         IF(WL.GT.8100.) FACO=.979+.021/(1.+((WL-8100.)/550.)**10.)
c        END IF
        if(nint(instr(10)).gt.0) then
         FACO=FACO*M(6,H,1)
        end if
       END IF
       
*** individuazione punti da mediare
       if(h.gt.1) then
        WLa=WL-(M(7,H,1)-M(7,H-1,1))/2.
       else
        WLa=WL-(M(7,H+1,1)-M(7,H,1))/2.
       end if
       if(h.lt.201) then
        WLb=WL+(M(7,H+1,1)-M(7,H,1))/2.
       else
        WLb=WL+(M(7,H,1)-M(7,H-1,1))/2.
       end if
       N=N1
       DO WHILE(X(N).LT.WLa.AND.(N+STEP).LE.NDATI.AND.(N+STEP).GE.1)
        N=N+STEP
       END DO
       if(iw.eq.1) then
        write(*,*)
        write(*,*)'WLa=',WLa,' WLb=',WLb
        write(*,*)'N=',n,' X(N)=',x(n)
        write(*,*)'Dati compresi tra WLa e WLb:'
       end if
       NPM=0
       J=1
       DO WHILE(X(n).ge.WLa.and.X(n).le.WLb.AND.n.LE.NDATI.and.n.ge.1)
        NPM=NPM+1
        ns(j)=n
        xs(j)=x(n)
        ys(j)=y(n)
        if(iw.eq.1)
     +   write(*,*)xs(j),ys(j)
        j=j+1
        n=n+step
        if(n.gt.ndati.or.n.lt.1) goto 200
       END DO
200    continue       
       
*** analisi situazione
       if(npm.lt.3) then
        if(iw.eq.1)
     +   write(*,*)'... npm =',npm,' < 3 -> impongo npm=3'
        if(npm.ge.1) then
         n=ns(1)
        else if(npm.eq.0) then
         npm=1
         ns(1)=n
        end if
        if(npm.eq.1.and.n.gt.1.and.n.lt.ndati.and.
     +     x(n+step).gt.WLb) n=n-step
        do while((n+step*2).lt.1.or.(n+step*2).gt.ndati)
         n=n-step
        end do
        npm=3
        do j=1,npm
         ns(j)=n
         xs(j)=x(n)
         ys(j)=y(n)
         if(iw.eq.1)
     +    write(*,*)xs(j),ys(j)
         n=n+step
        end do
       else if(npm.gt.3.and.i.ge.21) then !misure ellissometriche
        if(iw.eq.1)
     +   write(*,*)'...misure ELLI -> impongo npm=3'
        do while(x(n+step).lt.wl)
         n=n+step
        end do
        npm=3
        do j=1,npm
         ns(j)=n
         xs(j)=x(n)
         ys(j)=y(n)
         if(iw.eq.1)
     +    write(*,*)xs(j),ys(j)
         n=n+step
        end do
       end if
       
*** raccordo salto misure ellissometriche
       IF(I.GE.21) THEN
        OS=0.
        DO J=1,NPM-1
         IF(ABS(YS(J+1)-YS(J)).GT.250.) THEN
          IF((YS(J+1)-YS(J)).LT.0.) THEN
           OS=360.
          ELSE
           OS=-360.
          END IF
          YS(J+1)=YS(J+1)+OS
         END IF
        END DO
       END IF
        
*** fit/interpolazione     
       IF(NPM.GT.3) THEN !fit
        mFit=npm
        nFit=3
        lwa=mFit*nFit+5*nFit+mFit
        ALLOCATE(IWA(nFit),stat=ierror)
        ALLOCATE(p(nFit),stat=ierror)
        ALLOCATE(FVEC(mFit),stat=ierror)
        ALLOCATE(WA(lwa),stat=ierror)
        ALLOCATE(FJAC(mFit,nFit),stat=ierror)
        P(3)=0.
        do j=1,npm
         arse(j,1)=xs(j)
         arse(j,2)=ys(j)
         P(3)=P(3)+ys(j)
        end do
        P(1)=0.
        P(2)=0.
        P(3)=P(3)/real(npm)
        tol=par(25,1)
        CALL LMDIF1(FPAR,mFit,nFit,P,FVEC,TOL,INFO,IWA,WA,LWA)
        c=P(3)
        b=P(2)
        a=P(1)
        if((iw.eq.1)) then
         write(*,*)'Coefficienti di BF parabolico:'
         write(*,*)'A = ',A
         write(*,*)'B = ',B
         write(*,*)'C = ',C
         do iii=1,npm
          write(*,*) xs(iii),' ',ys(iii)
         end do
        end if
        DEALLOCATE(IWA,stat=ierror)
        DEALLOCATE(p,stat=ierror)
        DEALLOCATE(FVEC,stat=ierror)
        DEALLOCATE(WA,stat=ierror)
        DEALLOCATE(FJAC,stat=ierror)
       ELSE !interpolazione su 3 punti
        IF(ABS(XS(1)-XS(2)).LT..01.OR.ABS(XS(2)-XS(3)).LT..01) THEN
         WRITE(*,*)
     *    'ATTENZIONE IL FILE ',I,' CONTIENE 2 DATI DI EGUALE ASCISSA'
        END IF
        XD1=DBLE(XS(1))
        XD2=DBLE(XS(2))
        XD3=DBLE(XS(3))
        YD1=DBLE(YS(1))
        YD2=DBLE(YS(2))
        YD3=DBLE(YS(3))
        DETD=XD1*XD1*(XD2-XD3)-XD1*(XD2*XD2-XD3*XD3)
        DETD=DETD+XD2*XD2*XD3-XD3*XD3*XD2
        AD=YD1*(XD2-XD3)-XD1*(YD2-YD3)
        AD=(AD+YD2*XD3-YD3*XD2)/DETD
        BD=XD1*XD1*(YD2-YD3)-YD1*(XD2*XD2-XD3*XD3)
        BD=(BD+XD2*XD2*YD3-XD3*XD3*YD2)/DETD
        CD=XD1*XD1*(XD2*YD3-XD3*YD2)
        CD=CD-XD1*(XD2*XD2*YD3-XD3*XD3*YD2)
        CD=(CD+YD1*(XD2*XD2*XD3-XD3*XD3*XD2))/DETD
        A=SNGL(AD)
        B=SNGL(BD)
        C=SNGL(CD)
        if((iw.eq.1)) then
         write(*,*)'Coefficienti parabola per 3 punti:'
         write(*,*)'X: ',XD1,' ',XD2,' ',XD3
         write(*,*)'Y: ',YD1,' ',YD2,' ',YD3
         write(*,*)'A = ',A
         write(*,*)'B = ',B
         write(*,*)'C = ',C
        end if
       END IF
       YFIN=(A*M(7,H,1)*M(7,H,1)+B*M(7,H,1)+C)*FACO
       IF(I.GE.7.AND.I.LE.15) THEN
        IF(YFIN.LT..0.AND.ALARM.EQ.0) THEN
         length=len_trim(nank(i-7))
         IF(IUVIR.EQ.1) 
     +    WRITE(*,*)'n-file ',nank(i-7)(1:length),' < 0 !!!'
         IF(IUVIR.EQ.2) 
     +    WRITE(*,*)'k-file ',nank(i-7)(1:length),' < 0 !!!'
         IF(nint(par(9,1)).EQ.1)
     +    write(*,*)'         ... but k will be set >=0 !!!'
         ALARM=1
        END IF
        IF(nint(par(9,1)).EQ.1.AND.YFIN.LT..0)  YFIN=.0
        M(I,H,IUVIR)=YFIN
       END IF
       IF(I.LE.6) THEN
        M(17,H,1)=1. !misura abilitata
        M(I,H,1)=YFIN
        DL=M(7,H,2)
        IF(I.EQ.1.OR.I.EQ.2)           M(I,H,2)=M(I,H,1)*INSTR(3)+DL
        IF(I.EQ.3.OR.I.EQ.4.OR.I.EQ.5) M(I,H,2)=M(I,H,1)*
     *                                        (INSTR(3)+INSTR(4))+DL
        if(iw.eq.1)write(*,*)'i=',i,' h=',h,' M=',M(i,h,1)
       END IF
       IF(I.GE.21.AND.I.LE.28) THEN
        IF(YFIN.LT.-180.) YFIN=YFIN+360.
        IF(YFIN.GT.180.) YFIN=YFIN-360.
        E(I-20,H,IUVIR)=YFIN
       END IF
       ymin=MIN(ymin,yfin)
       ymax=MAX(ymax,yfin)
      END DO
      if(i.ge.1.and.i.le.6) then
       rxy(i,3)=ymin*100.
       rxy(i,4)=ymax*100.
      else if(i.ge.21.and.i.le.28.and.iuvir.eq.1) then
       rxy(i-14,3)=ymin
       rxy(i-14,4)=ymax
      end if
      if(iw.eq.1) then
       write(*,*)'Proseguo ?'
       read(*,'(a1)') r
      end if
      RETURN
      END



      SUBROUTINE CXY(IW,RXY,IC,IXW,iklog)
      INTEGER ic,pgopen
      INTEGER q(2),ixw(20)
      REAL RXY(30,4)
      CHARACTER r*1,LABXY(2)*8
      
c  IW  : n. di finestra grafica
c  ixw : registro finestre grafiche     N.   spazio     
c                                        1    n-k 
c                                        2    /
c                                        3    wl-n
c                                        4    wl-k
c                                        5    wl-FM
c                                        6    wl-derivata
c                                        7    d-n,k,T,R
c                                        8    wl-eps1
c                                        9    wl-eps2
c                                       ..    ....
c                                       11    wl-Tn 
c                                       12    wl-Tp
c                                       13    wl-Rn
c                                       14    wl-Rp
c                                       15    wl-R1   
c                                       16    wl-Apds
c                                       17    wl-DELTA
c                                       18    wl-PSI
c                                       19    wl-Aaria
c                                       20    teta-T,R,tau,rho  
c  IC : controllo richieste: 
c       
c       1 -> mantieni grafico e imponi viewport
c       
c       3 -> imponi viewport, cancella e traccia senza chiedere consenso
c 10+imis -> dato iw mantieni grafico,imponi valori estremi Y
c            relativi imis e grafica
!      call pgqci(icol)
!      call pgqls(ils)
      
      if(ic.le.3) then !impongo rxy(25,1) e rxy(25,2) dato iw
c           1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20    
       goto(1,2,3,4,5,6,7,8,9,19,11,11,11,11,11,11,12,12,13,14) iw   
 1      rxy(25,1)=16.
        rxy(25,2)=17.
       goto 19
 2      rxy(25,1)=16.
        rxy(25,2)=22.
       goto 19
 3      rxy(25,1)=20.
        rxy(25,2)=16.
       goto 19
 4      rxy(25,1)=20.
        rxy(25,2)=17. 
       goto 19
 5      rxy(25,1)=20.
        rxy(25,2)=15.
       goto 19
 6      rxy(25,1)=20.
        rxy(25,2)=19.     
       goto 19
 7      rxy(25,1)=22
c       rxy(25,2) imposto esternamente     
       goto 19
 8      rxy(25,1)=20.
        rxy(25,2)=26.
       goto 19
 9      rxy(25,1)=20.
        rxy(25,2)=27. 
       goto 19
 11     rxy(25,1)=20.
        rxy(25,2)=real(iw-10)
       goto 19
 12     rxy(25,1)=20. 
        rxy(25,2)=real(iw-10)
c       rxy(25,2) imposto esternamente 
       goto 19
 13     rxy(25,1)=20.
        rxy(25,2)=18.
       goto 19
 14     rxy(25,1)=21.
        rxy(25,2)=23
 19    continue
      else !considero "ic-10" come la misura da graficare su iw
       if(ic.ge.11.and.ic.le.24) then
        rxy(25,1)=20.
        rxy(25,2)=real(ic-10)
       end if
      end if
      write(*,*)'ixw(',iw,') = ',ixw(iw)
      r=' '
      if(ixw(iw).le.0) then
       ixw(iw) = pgopen( '/XWINDOW' )
       if(ixw(iw).le.0) then
        write(*,*)'Re-launch the project and use less windows!'
        return
       end if
       call PGSCR (0,1.0,1.0,1.0) !background color = white
       call PGSCR (1,0.0,0.0,0.0) !foreground color = black
       call PGPAP (rxy(24,1), rxy(24,2)) !size of view surface
       call PGSCH (2.0)!character size
       CALL PGSLW(nint(rxy(24,3)))
       r='t'
       write(*,*)'CXY new window: ixw(',iw,') = ',ixw(iw)
       if(ixw(iw).le.0) stop
       inew=1
       call pgsvp(0.08,0.96,0.17,0.96)!(.055,.7,.62,.96) size of viewport
      else
       call pgslct(ixw(iw))
       inew=0
      end if
      q(1)=nint(rxy(25,1))
      q(2)=nint(rxy(25,2))
!       if(ic.eq.3) then
!        rxy(q(2),1)=rxy(q(2),3)
!        rxy(q(2),2)=rxy(q(2),4)
!       end if
      if(ic.eq.3) r='t'
      if(q(1).eq.20.or.q(2).eq.20) then
       LEV=nint(rxy(25,3))
      else
       LEV=1
      end if
      do I=1,2
c            1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
       goto(21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,49,40,
c           21 22 23,24,25,26,27
     +      41,42,43,49,49,46,47) q(i)
 21     LABXY(I)='Tn      '
       goto 49
 22     LABXY(I)='Tp      ' 
       goto 49
 23     LABXY(I)='Rn      '
       goto 49
 24     LABXY(I)='Rp      '
       goto 49
 25     LABXY(I)='R1      '
       goto 49
 26     LABXY(I)='Aps     '
       goto 49
 27     LABXY(I)='Delta_1 '
       goto 49
 28     LABXY(I)='Psi_1   '
       goto 49
 29     LABXY(I)='Delta_2 '
       goto 49
 30     LABXY(I)='Psi_2   '
       goto 49
 31     LABXY(I)='Delta_3 '
       goto 49
 32     LABXY(I)='Psi_3   '
       goto 49
 33     LABXY(I)='Delta_4 '
       goto 49
 34     LABXY(I)='Psi_4   '
       goto 49
 35     LABXY(I)='Merit F.'
       goto 49
 36     LABXY(I)='n       '
       goto 49
 37     LABXY(I)='k       '
       goto 49
 38     LABXY(I)='Asf     '
       goto 49
 40     if(LEV.eq.1) then
         LABXY(I)='Angstrom'
        else
         LABXY(I)='eV      '
        end if
       goto 49
 41     LABXY(I)='Teta (d)'
       goto 49
 42     LABXY(I)='d (Angs)'
       goto 49
 43     LABXY(I)='TaRhoRh1'
       goto 49
 46     LABXY(I)='epsi1'
       goto 49
 47     LABXY(I)='epsi2'
 49    CONTINUE
      end do
      call pgsci(1)
      call pgsls(1)
      if(iklog.eq.1) then
       if(rxy(q(2),1).lt.0.) then
        rxy(q(2),1)=rxy(q(2),4)/1000.
       end if
      end if
c      do while(r.ne.' '.and.r.ne.'t')
       write(*,*)'**** PARAMETRI GRAFICO '//LABXY(2)//' Vs '//LABXY(1)
       write(*,*)'Asse  Quantita"      Min       Max'
       write(*,*)'X:    ',LABXY(1),'    ',rxy(q(1),1),'  ',rxy(q(1),2)
       write(*,*)'  valori estremi  ',rxy(q(1),3),'  ',rxy(q(1),4)
       write(*,*)'Y:    ',LABXY(2),'    ',rxy(q(2),1),'  ',rxy(q(2),2)
       write(*,*)'  valori estremi  ',rxy(q(2),3),'  ',rxy(q(2),4)
       write(*,*)
c       if(LEV.ne.0) then
c        write(*,*)'"L" -> Angstrom (eV) -> eV (Angstrom)'
c       end if
       if(rxy(q(1),1).ge.rxy(q(1),2)) then
        if(rxy(q(1),3).ge.rxy(q(1),4)) then
         rxy(q(1),4)=rxy(q(1),3)+1.;
        end if
        rxy(q(1),1)=rxy(q(1),3)
        rxy(q(1),2)=rxy(q(1),4)
       end if
       if(rxy(q(2),1).ge.rxy(q(2),2)) then
        if(rxy(q(2),3).ge.rxy(q(2),4)) then
         rxy(q(2),4)=rxy(q(2),3)+1.;
        end if
        rxy(q(2),1)=rxy(q(2),3)
        rxy(q(2),2)=rxy(q(2),4)
       end if
c      end do
      xmin=rxy(q(1),1)
      xmax=rxy(q(1),2)
      if(q(1).eq.20.and.LEV.eq.2) then
       xmin=12400./xmin
       xmax=12400./xmax
      end if
      ymin=rxy(q(2),1)
      ymax=rxy(q(2),2)
      if(iklog.eq.1) then
       ymin=log10(ymin)
       ymax=log10(ymax)
      end if
      if(q(2).eq.20.and.LEV.eq.2) then
       ymin=12400./ymin
       ymax=12400./ymax
      end if 
      if(r.eq.' ') then
       call pgswin(xmin,xmax,ymin,ymax)
      else
       call pgeras
       call pgswin(xmin,xmax,ymin,ymax)
       if(iklog.eq.0) then
        call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
       else
        call pgbox('BCNST',0.0,0,'BCLNST',0.0,0)
       end if
       call pglab(labxy(1),labxy(2),'')
       IF(YMIN.LT.0..AND.YMAX.GT.0..AND.iklog.ne.1) THEN
        CALL pgmove(xmin,0.)
        CALL pgdraw(xmax,0.)
       END IF
       IF(YMIN.LT.180..AND.YMAX.GT.180.) THEN
        CALL pgmove(XMIN,180.)
        CALL pgdraw(XMAX,180.)
       END IF
      end if
      if(q(1).eq.20.or.q(2).eq.20) rxy(25,3)=LEV
!      call pgsci(icol)
!      call pgsls(ils)
      return
      end



      SUBROUTINE PLTDAT
      REAL M(17,201,2),SOL(999,6),XV(201),YV(201),ELI(8,201,2),PAR(60,5)
     *    ,CNK(16,3),PF(7,21),RXY(30,4),PM(200,5),ARSE(500,2)
      INTEGER LEV,IXW(20)
      CHARACTER NANK(16)*256
      COMMON /VALORI/ M,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
*** subroutine per il grafico delle misure sperimentali selezionate in
***            autos
      iw1=0
      icur=3 !imponi vieport, cancella e traccia
      LEV=nINT(RXY(25,3))
      DO I=1,14
       IF(nint(par(I,4)).EQ.2) THEN
        if(i.ge.1.and.i.le.6) iw=i+10
        if(i.eq.7.or.i.eq.9.or.i.eq.11.or.i.eq.13) iw=17
        if(i.eq.8.or.i.eq.10.or.i.eq.12.or.i.eq.14) iw=18
        rxy(25,2)=i
        if(iw1.eq.0) iw1=iw
        call CXY(iw1,rxy,icur,IXW,0)
        icur=4 !mantieni il grafico, imponi viewport, ma non rxy(25,2)
        LEV=nint(rxy(25,3))
        CALL pgsls(1)
        IF(I.LE.6) CALL pgsci(I)
        IF(I.EQ.7.OR.I.EQ.9.OR.I.EQ.11.OR.I.EQ.13) CALL pgsci(7)
        IF(I.EQ.8.OR.I.EQ.10.OR.I.EQ.12.OR.I.EQ.14) CALL pgsci(8)
        DO L=1,201
         XV(L)=M(7,L,1)
         IF(LEV.EQ.2) XV(L)=12400./XV(L)
         IF(I.LE.6) THEN
          YV(L)=M(I,L,1)*100.
         ELSE
          YV(L)=ELI(I-6,L,1)
         END IF
        END DO
        CALL pgline(201,XV,YV)
       END IF
      END DO
      RETURN
      END
      
      

      SUBROUTINE PLOTNK(ink,ir)
      INTEGER IXW(20)
      REAL SOL(999,6),PF(7,21),PM(200,5),M(17,201,2)
     *    ,E(8,201,2),PAR(60,5),CNK(16,3),RXY(30,4),ARSE(500,2)
      COMMON /VALORI/ M,E,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      CHARACTER NANK(16)*256
*** subroutine per il grafico di nk
c     ir = parametro di controllo richieste di CXY
      
      iklog=0
      if(ink.eq.1) then
       nk=1
       iw=3
      else if(ink.eq.2) then
       nk=2
       iw=4
       if(par(31,5).gt.0.5) then
        iklog=1
       end if
      else if(ink.eq.3) then
       iw=8
       rxy(26,3)=1000.
       rxy(26,4)=-1000.
      else if(ink.eq.4) then
       iw=9
       rxy(27,3)=1000.
       rxy(27,4)=-1000.
      end if
      CALL CXY(iw,RXY,ir,IXW,iklog)
      CALL PGQWIN(XL,XH,YL,YH)
      call pgsls(1)
      LEV=nINT(RXY(25,3))
      DO I=2,nINT(SOL(1,1))+1
       XX=SOL(I,1)
       IF(LEV.EQ.2) XX=12400./XX
       if(ink.eq.1.or.ink.eq.2) then
        YD=SOL(I,NK+1)-SOL(I,NK+3)
        YU=SOL(I,NK+1)+SOL(I,NK+3)
        if(nk.eq.2.and.(par(31,5).gt.0.5)) then
         yd=log10(yd)
         yu=log10(yu)
        end if
       else if(ink.eq.3) then
        v1=(SOL(I,1+1)-SOL(I,1+3))**2.-(SOL(I,2+1)-SOL(I,2+3))**2.
        v2=(SOL(I,1+1)-SOL(I,1+3))**2.-(SOL(I,2+1)+SOL(I,2+3))**2.
        v3=(SOL(I,1+1)+SOL(I,1+3))**2.-(SOL(I,2+1)-SOL(I,2+3))**2.
        v4=(SOL(I,1+1)+SOL(I,1+3))**2.-(SOL(I,2+1)+SOL(I,2+3))**2.
        yd=MIN(v1,v2)
        yd=MIN(yd,v3)
        yd=MIN(yd,v4)
        yu=MAX(v1,v2)
        yu=MAX(yd,v3)
        yu=MAX(yd,v4)
        rxy(26,3)=MIN(rxy(26,3),yd)
        rxy(26,4)=MAX(rxy(26,4),yu)
       else if(ink.eq.4) then
        v1=2.*(SOL(I,1+1)-SOL(I,1+3))*(SOL(I,2+1)-SOL(I,2+3))
        v2=2.*(SOL(I,1+1)-SOL(I,1+3))*(SOL(I,2+1)+SOL(I,2+3))
        v3=2.*(SOL(I,1+1)+SOL(I,1+3))*(SOL(I,2+1)-SOL(I,2+3))
        v4=2.*(SOL(I,1+1)+SOL(I,1+3))*(SOL(I,2+1)+SOL(I,2+3))
        yd=MIN(v1,v2)
        yd=MIN(yd,v3)
        yd=MIN(yd,v4)
        yu=MAX(v1,v2)
        yu=MAX(yd,v3)
        yu=MAX(yd,v4)
        rxy(27,3)=MIN(rxy(27,3),yd)
        rxy(27,4)=MAX(rxy(27,4),yu)
       end if
       IF(((XX.GE.XL.AND.XX.LE.XH).or.(XX.LE.XL.AND.XX.GE.XH)).AND.
     +     YD.GE.YL.AND.YU.LE.YH) THEN
        CALL pgmove(XX,YD)
        IF(nint(SOL(I,6)).EQ.0) CALL pgsci(2)
        IF(nint(SOL(I,6)).EQ.1) CALL pgsci(3)
        CALL pgdraw(XX,YU)
       END IF
      END DO
      if(ink.eq.3) then
       write(*,*) 'epsi1:',rxy(26,3),'<->',rxy(26,4)
      else if(ink.eq.4) then
       write(*,*) 'epsi2:',rxy(27,3),'<->',rxy(27,4)
      end if
      RETURN
      END
        
        
      
      SUBROUTINE PLOTMIS(ic,iw,ir)
      INTEGER IXW(20)
      REAL M(17,201,2),E(8,201,2),PAR(60,5),PM(200,5),CNK(16,3),
     +     SOL(999,6),PF(7,21),RXY(30,4),ARSE(500,2)
      COMMON /VALORI/ M,E,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      CHARACTER NANK(16)*256
*** subroutine per il grafico delle misure sperimentali
c     ic = indice della misura
c     iw = N. della finestra grafica su cui operare
c     ir = parametro di controllo richieste di CXY

      if(ic.lt.1.or.ic.gt.15) goto 1
!      call pgqci(icol)
!      call pgqls(ilstyle)
      call cxy(iw,rxy,ir,IXW,0)
      call pgqwin(XL,XH,YL,YH)
      J=0
      OFF=0.
      IF(IC.GE.1.AND.IC.LE.6) icolnew=IC
      IF(IC.EQ.7.OR.IC.EQ.9.OR.IC.EQ.11.OR.IC.EQ.13) icolnew=7
      IF(IC.EQ.8.OR.IC.EQ.10.OR.IC.EQ.12.OR.IC.EQ.14) icolnew=8
      if(ic.eq.15) icolnew=1
      IF(IC.GE.7.AND.IC.LE.14) YOLD=E(IC-6,1,1)
      call pgsls(1)
      DO L=1,201
       XX=M(7,L,1)
       IF(nint(rxy(25,3)).eq.2) XX=12400./XX
       IF((XX.GE.XL.AND.XX.LE.XH).or.(XX.GE.XH.AND.XX.LE.XL)) THEN
        IF(IC.LE.6) then
         Y=M(IC,L,1)*100
         if(nint(m(17,L,1)).eq.1) then
          call pgsci(icolnew)
         else
          call pgsci(icolnew+1)
         end if
        end if
!        IF(IC.EQ.6) Y=M(IC,L,1)
        IF(IC.GE.7.AND.IC.LE.14) THEN
         Y=E(IC-6,L,1)+OFF
         DO WHILE(ABS(Y-YOLD).GT.300.)
          IF((Y-YOLD).GT.300.) OFF=OFF-360.
          IF((Y-YOLD).LT.-300.) OFF=OFF+360.
          Y=E(IC-6,L,1)+OFF
         END DO
         YOLD=Y
        END IF
        IF(IC.EQ.15) Y=1-M(1,L,1)-M(3,L,1)
        IF(Y.GE.YL.AND.Y.LE.YH) THEN
         IF(J.EQ.0) CALL pgmove(XX,Y)
!         IF(J.EQ.1) CALL pgdraw(XX,Y)
         IF(J.EQ.1) CALL pgpt(1,XX,Y,17)
         J=1
        ELSE
         J=0
        END IF
       END IF
      END DO
!      call pgsci(icol)
!      call pgsls(ilstyle)
 1    RETURN
      END
      
      
      SUBROUTINE FPAR(M,N,P,FVEC,IFLAG)
      INTEGER M,N,IFLAG,IXW(20),iiflag
      REAL SOL(999,6),PF(7,21),PM(200,5),MIS(17,201,2),ELI(8,201,2),
     *     PAR(60,5),CNK(16,3),RXY(30,4),ARSE(500,2)
      REAL P(n),FVEC(m)
      CHARACTER NANK(16)*256
      COMMON /VALORI/ MIS,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
      iiflag=iflag!per evitare warning unused
*** subroutine per la costruzione di fvec
      do i=1,m
        fvec(i)=ARSE(i,2)-P(1)*ARSE(i,1)*ARSE(i,1)-P(2)*ARSE(i,1)-P(3)
      end do
      RETURN
      END
