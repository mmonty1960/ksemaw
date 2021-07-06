*********    kSEMAW KERNEL sources   ***************
*
*   hereinafter subroutines to calculate
*   T R R1 Apds PSI DELTA for s and p polarization
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
C     VNK(16,2) : vettore contenente n,k relativi a
C     VNK(1,1),VNK(1,2) : n,k incogniti
C     VNK(2,1),VNK(2,2) : n,k NOTI
C     VNK(3,1),VNK(3,2) : n,k NOTI
C     VNK(4,1),VNK(4,2) : n,k NOTI
C     VNK(5,1),VNK(5,2) : n,k NOTI
C     VNK(6,1),VNK(6,2) : n,k NOTI
C     VNK(7,1),VNK(7,2) : n,k NOTI
C     VNK(8,1),VNK(8,2) : n,k NOTI
C     VNK(9,1),VNK(9,2) : n,k NOTI
C     VNK(10,1),VNK(10,2) : n,k MEZZO_ingresso Tn,..,R1
C     VNK(11,1),VNK(11,2) : n,k MEZZO_ingresso Apds
C     VNK(12,1),VNK(12,2) : n,k MEZZO_ingresso PSI,DELTA #1
C     VNK(13,1),VNK(13,2) : n,k MEZZO_ingresso PSI,DELTA #2
C     VNK(14,1),VNK(14,2) : n,k MEZZO_ingresso PSI,DELTA #3
C     VNK(15,1),VNK(15,2) : n,k MEZZO_ingresso PSI,DELTA #4
C     VNK(16,1),VNK(10,2) : n,k MEZZO_uscita   Tn,..,R1


      SUBROUTINE ASSEMBLER(iwl,wl,ikind,teta,par,pm,vot)
      REAL wl,VNK(16,2),PAR(60,5),PM(200,5),VOT(5,2),VOSI(5,2)
      COMPLEX IRUP,IRDW,IRUP2,IRDW2,MUPS,MUPP,MDWS,MDWP,PQ,RR
      DATA PIG /3.141592654/
*** Subroutine di assemblaggio di tutte le interfacce del
**  multistrato generalizzato composto da strati coerenti e non
C     iwl: indice di wl
C
C     IKIND tipo di misura
C         1 SF
C         2 PDS
C         3 ELI_1
C         .......
C         6 ELI_4
C
C     VOT(5,2): vettore misure ottiche calcolate
C              1     2
C         1    Ts    Tp
C         2    Rs    Rp
C         3    R1s   R1p
C         4    Apds  /
C         5    PSI   DELTA
      
*** inizializzazione VOT e ALFA
      vot(1,1)=1.
      vot(1,2)=1.
      do j=2,5
       vot(j,1)=0.
       vot(j,2)=0.
      end do
      alfa=1.
      
*** indici di rifrazione e inizializzazione n_ingresso
      call COSVNK(1,VNK,iwl)
      ivnkup=9+ikind !mezzo di ingresso (dipende dal tipo di misura)
      irup=CMPLX(vnk(ivnkup,1),-vnk(ivnkup,2)) 
      irup2=irup*irup
      pq=(irup*sin(teta))**2. ! (n_ingresso*sin(teta))**2.
      
*** parametri multistrato
      imax=nint(par(51,2))
      
*** loop
      i=1
      do while(i.le.imax+1)
         
*** calcolo R T alla i-esima interfaccia
       if(i.le.imax.and.
     +    (nint(par(50+i,3)).gt.1.or.pm(50+i,1).gt..0)) then
          
***** film          
        if(nint(par(50+i,3)).gt.1) then
         ifst=i
         ncoe=1
         do while(ifst+ncoe.lt.imax.and.nint(par(50+ifst+ncoe,3)).gt.1)
          i=i+1
          ncoe=ncoe+1
         end do
         call BUILDER(iwl,wl,ikind,ifst,ncoe,pq,par,pm,vosi)
         vot(4,1)=vosi(4,1) !salvo Apds
         vot(5,1)=vosi(5,1) !salva PSI
         vot(5,2)=vosi(5,2) !salva DELTA
         i=i+1 !film(s) accorpato(i) allo strato bulk successivo
         
***** bulk rugoso
        else if(pm(50+i,1).gt..0) then !bulk rugoso
         ifst=i
         ncoe=0
         call BUILDER(iwl,wl,ikind,ifst,ncoe,pq,par,pm,vosi)!pq era teta
        end if
        
***** settaggio irdw
        if(i.le.imax) then
         ivnkdw=nint(par(50+i,1))
        else !il mezzo successivo è quello di uscita SF
         ivnkdw=16
        end if
        irdw=CMPLX(vnk(ivnkdw,1),-vnk(ivnkdw,2))
c        if(i.eq.imax+1.and.iwl.eq.201) then
c         write(*,*)'@1 vnk(',ivnkdw,'q1)=',vnk(ivnkdw,1),vnk(ivnkdw,2)
c        end if
        irdw2=irdw*irdw
        mdws=SQRT(irdw2-pq)
        
***** bulk liscio
       else
        if(i.le.imax.or.(i.eq.(imax+1).and.nint(par(52,2)).eq.0)) then
         if(i.gt.1) then
          ivnkup=nint(par(50+i-1,1))
          irup=CMPLX(vnk(ivnkup,1),-vnk(ivnkup,2))
          irup2=irup*irup
         end if
         if(i.le.imax) then
          ivnkdw=nint(par(50+i,1))
         else !il mezzo successivo è quello di uscita SF
          ivnkdw=16
         end if
c         if(i.eq.imax+1.and.iwl.eq.201) then
c          write(*,*)'@2 vnk(',ivnkdw,'q2)=',vnk(ivnkdw,1),vnk(ivnkdw,2)
c         end if
         irdw=CMPLX(vnk(ivnkdw,1),-vnk(ivnkdw,2))
         irdw2=irdw*irdw
         mups=SQRT(irup2-pq)
         mupp=irup2/mups
         mdws=SQRT(irdw2-pq)
         mdwp=irdw2/mdws
         vosi(1,1)=4.*REAL(mups)*REAL(mdws)/(ABS(mups+mdws))**2.
         vosi(1,2)=4.*REAL(mupp)*REAL(mdwp)/(ABS(mupp+mdwp))**2.
         vosi(2,1)=(ABS((mups-mdws)/(mups+mdws)))**2.
         vosi(2,2)=(ABS((mupp-mdwp)/(mupp+mdwp)))**2.
         vosi(3,1)=vosi(2,1)
         vosi(3,2)=vosi(2,2)
         if(i.eq.1.and.imax.eq.1) then !BARE SUBSTRATE
*** Apds
          vot(4,1)=4.*REAL(irup*(CONJG(mdws)-irdw))/(ABS(irup+mdws))**2.
*** PSI e DELTA
          rr=((mupp-mdwp)/(mupp+mdwp))/((mups-mdws)/(mups+mdws))
          vot(5,1)=ATAND1(ABS(RR))
          vot(5,2)=ATAN2D1(-AIMAG(RR),-REAL(RR))
         end if
        else !multistrato simmetrico su 2a faccia ultimo strato
         do ii=1,2
          vosi(1,ii)=vot(1,ii)
          vosi(2,ii)=vot(3,ii)
          vosi(3,ii)=vot(2,ii)
         end do
        end if
       end if
       
*** combinazione con R T dell'interfaccia precedente
       uguales=alfa/(1.-vot(3,1)*vosi(2,1)*alfa*alfa)
       ugualep=alfa/(1.-vot(3,2)*vosi(2,2)*alfa*alfa)
       ts=vot(1,1)*vosi(1,1)*uguales
       tp=vot(1,2)*vosi(1,2)*ugualep
       rs=vot(2,1)+vot(1,1)**2.*vosi(2,1)*alfa*uguales
       rp=vot(2,2)+vot(1,2)**2.*vosi(2,2)*alfa*ugualep
       r1s=vosi(3,1)+vosi(1,1)**2.*vot(3,1)*alfa*uguales
       r1p=vosi(3,2)+vosi(1,2)**2.*vot(3,2)*alfa*ugualep
       vot(1,1)=ts
       vot(1,2)=tp
       vot(2,1)=rs
       vot(2,2)=rp
       vot(3,1)=r1s
       vot(3,2)=r1p
       
*** calcolo attenuazione lungo lo spessore dello strato bulk
       alfa=EXP(-4.*PIG*pm(i,1)/wl*AIMAG(-mdws))       
       
       i=i+1 
      end do
      
      RETURN
      END
      
      
      
      SUBROUTINE BUILDER(iwl,wl,ikind,ifst,ncoe,pq,par,pm,vosi)
      INTEGER NFA,iv(10,2),irougfa(10),irougms(10)
      REAL VNK(16,2),PAR(60,5),PM(200,5),VOSI(5,2),ds(10),gn(10),cn(10),
     +     gk(10),ck(10),ru(10),nmedio,kmedio,w(5,6),dx(5,6),
     +     al(10),d(999),wl
      COMPLEX IR(999),PQ,NQ,MUS,MUP,rhoS,rhopS,tauS,rhoP,rhopP,tauP,
     +        nq1,mus1,mup1,out(8,2)
c      PARAMETER(PIG=3.141592654,CI=(0.,1.))
c w=1/sqrt(2.*pig)*exp(-0.5*(N*dx)**2)*3./No N=0,1,2..,No
c all data are renormalized to satisfy SUM weigth=1.0
      DATA w(1,1),w(1,2),w(1,3),w(1,4),w(1,5),w(1,6)
     +     /9.78264917E-1,1.08675416E-2,0.0,0.0,0.0,0.0/
      DATA w(2,1),w(2,2),w(2,3),w(2,4),w(2,5),w(2,6)
     +     /5.9826E-1,1.9423E-1,6.6460E-3,0.0,0.0,0.0/
      DATA w(3,1),w(3,2),w(3,3),w(3,4),w(3,5),w(3,6)
     +     /0.39905,0.242036,5.40056E-2,4.43305E-3,0.0,0.0/
      DATA w(4,1),w(4,2),w(4,3),w(4,4),w(4,5),w(4,6)
     +     /2.99373E-1,2.25978E-1,9.7192E-2,2.38179E-2,3.32573E-3,0.0/
      DATA w(5,1),w(5,2),w(5,3),w(5,4),w(5,5),w(5,6)
     +     /2.39559E-1,2.00097E-1,1.16606E-1,4.74085E-2,1.34476E-2,
     +      2.66126E-3/
      DATA dx(1,1),dx(1,2),dx(1,3),dx(1,4),dx(1,5),dx(1,6)
     +     /0.0,3.0,0.0,0.0,0.0,0.0/
      DATA dx(2,1),dx(2,2),dx(2,3),dx(2,4),dx(2,5),dx(2,6)
     +     /0.0,1.5,3.0,0.0,0.0,0.0/
      DATA dx(3,1),dx(3,2),dx(3,3),dx(3,4),dx(3,5),dx(3,6)
     +     /0.0,1.0,2.0,3.0,0.0,0.0/
      DATA dx(4,1),dx(4,2),dx(4,3),dx(4,4),dx(4,5),dx(4,6)
     +     / 0.0,0.75,1.5,2.25,3.0,0.0/
      DATA dx(5,1),dx(5,2),dx(5,3),dx(5,4),dx(5,5),dx(5,6)
     +     /0.0,0.6,1.2,1.8,2.4,3.0/
*** Subroutine di costruzione della porzione coerente del multistrato
**             e gestione della rugosità
c     
c   iv(10,3) vettore abilitazione e indicizzazione
c     (k,1) : 0 (interfaccia piana), 1 (intrefaccia rugosa)
c     (k,2) : ik indice del termine della sommatoria
      
*** parametri vari
      nino=nint(par(29,1)) !discretizzazione inomogeneit
      N=nint(par(28,1))    !integrale numerico rigosit
      
*** parametri identificativi strati      
      do i=1,10
       ds(i)=pm(i,1)      !spessore strato
       al(i)=par(50+i,3)  !tipo di strato
       gn(i)=pm(10+i,1)+12400./wl*pm(60+i,1)   !Dn/<n>
       cn(i)=pm(20+i,1)   !curv_n
       gk(i)=pm(30+i,1)   !Dk/<k>
       ck(i)=pm(40+i,1)   !curv_k
       ru(i)=pm(50+i,1)   !roughness
       iv(i,1)=0
       if(ru(i).gt..0) iv(i,1)=1
      end do      
       
*** parametri multistrato
      imax=nint(par(51,2))
      
*** indici di rifrazione e inizializzazione n_ingresso
      call COSVNK(1,VNK,iwl)
      ivnk1=9+ikind
      ir(1)=CMPLX(vnk(ivnk1,1),-vnk(ivnk1,2))
      if(ifst.gt.1) then
       ivnk1=nint(par(50+ifst-1,1))
       ir(1)=CMPLX(vnk(ivnk1,1),-vnk(ivnk1,2))
      end if
      nfa=1
*     
**       
******* costruzione vettore indici di rifrazione e spessori
**    
*     
      ims=ifst
      nroug=0
      do while(ims.le.(ifst+ncoe-1).or.ru(ims).gt..0)
         
** rugosità
       if(ru(ims).gt.0.) then !aggiungi bistrato per rugosità
        nroug=nroug+1
        nfa=nfa+1
        irougfa(nroug)=nfa
        irougms(nroug)=ims
        ir(nfa)=ir(nfa-1)
        d(nfa)=3.*ru(ims)
        nfa=nfa+1
        ivnk=nint(par(50+ims,1))
        ir(nfa)=CMPLX(vnk(ivnk,1),-vnk(ivnk,2))
        d(nfa)=d(nfa-1)
        ds(ims)=ds(ims)-3.*ru(ims) !spessore al netto della rugosit
        if(ds(ims).le..0) goto 1 !solo bistrato rugoso
         
       end if
       
** inomogeneità
       if(nint(al(ims)).eq.5) then !strato inomogeneo exp _/
        ivnk=nint(par(50+ims,1))
        nmedio=vnk(ivnk,1)
        kmedio=vnk(ivnk,2)
        dz=1./REAL(nino)
        BpN=gn(ims)*nmedio
        do i=1,nino
         nfa=nfa+1
         zz=1.+dz/2.-real(i)*dz
         ir(nfa)=CMPLX(nmedio+bpn*exp(-3.*(1.-zz)),-kmedio)
         d(nfa)=dz*ds(ims)
        end do
        if(ru(ims).gt.0.) ir(nfa-nino)=ir(nfa-nino+1) !nk per roughness
       else if(nint(al(ims)).eq.4) then !strato inomogeneo exp \_
        ivnk=nint(par(50+ims,1))
        nmedio=vnk(ivnk,1)
        kmedio=vnk(ivnk,2)
        dz=1./REAL(nino)
        BpN=gn(ims)*nmedio
        do i=1,nino
         nfa=nfa+1
         zz=1.+dz/2.-real(i)*dz
         ir(nfa)=CMPLX(nmedio+bpn*exp(-3.*zz),-kmedio)
         d(nfa)=dz*ds(ims)
        end do
        if(ru(ims).gt.0.) ir(nfa-nino)=ir(nfa-nino+1) !nk per roughness
       else if(nint(al(ims)).eq.3) then !strato inomogeneo
        ivnk=nint(par(50+ims,1))
        nmedio=vnk(ivnk,1)
        kmedio=vnk(ivnk,2)
        dz=1./REAL(nino)
        ApN=6.*cn(ims)*nmedio
        ApK=6.*ck(ims)*kmedio
        BpN=gn(ims)*nmedio-ApN
        BpK=gk(ims)*kmedio-ApK
        CpN=nmedio-ApN/3.-BpN/2.
        CpK=kmedio-ApK/3.-BpK/2.
        do i=1,nino
         nfa=nfa+1
         zz=1.+dz/2.-real(i)*dz
         ir(nfa)=CMPLX(apn*zz*zz+bpn*zz+cpn,-(apk*zz*zz+bpk*zz+cpk))
         d(nfa)=dz*ds(ims)
        end do
        if(ru(ims).gt.0.) ir(nfa-nino)=ir(nfa-nino+1) !nk per roughness
       else if(nint(al(ims)).eq.2) then  !strato omogeneo
        nfa=nfa+1
        ivnk=nint(par(50+ims,1))
        ir(nfa)=CMPLX(vnk(ivnk,1),-vnk(ivnk,2))
        d(nfa)=ds(ims)        
       end if
       
 1     ims=ims+1
      end do

** substrato
      nfa=nfa+1
      if(ru(ims-1).gt..0.and.nint(al(ims-1)).eq.1) ims=ims-1 !bulk rugoso
      if(ims.le.imax) then
       ivnk=nint(par(50+ims,1))
      else !senza sub
       ivnk=16
      end if
      ir(nfa)=CMPLX(vnk(ivnk,1),-vnk(ivnk,2))
      
**** gestione roughness e calcolo coefficienti di Fresnel
      
      rhoS=CMPLX(.0,.0)
      rhopS=CMPLX(.0,.0)
      tauS=CMPLX(.0,.0)
      rhoP=CMPLX(.0,.0)
      rhopP=CMPLX(.0,.0)
      tauP=CMPLX(.0,.0)
      Rs=0.
      R1s=0.
      Rp=0.
      R1p=0.
      Ts=0.
      Tp=0.
      As=0.
      Ap=0.
      
      do i1=1,(2*N*iv(1,1)+1)
       do i2=1,(2*N*iv(2,1)+1)
        do i3=1,(2*N*iv(3,1)+1)
         do i4=1,(2*N*iv(4,1)+1)
          do i5=1,(2*N*iv(5,1)+1)
           do i6=1,(2*N*iv(6,1)+1)
            do i7=1,(2*N*iv(7,1)+1)
             do i8=1,(2*N*iv(8,1)+1)
              do i9=1,(2*N*iv(9,1)+1)
               do i10=1,(2*N*iv(10,1)+1)
                iv(1,2)=i1
                iv(2,2)=i2
                iv(3,2)=i3
                iv(4,2)=i4
                iv(5,2)=i5
                iv(6,2)=i6
                iv(7,2)=i7
                iv(8,2)=i8
                iv(9,2)=i9
                iv(10,2)=i10
                ilyru=1
                wtot=1.
                do i=1,nroug
                 do while(iv(ilyru,1).ne.1.and.ilyru.lt.10)
                  ilyru=ilyru+1
                 end do
                 ra=(-1)**iv(ilyru,2)*dx(N,int(iv(ilyru,2)/2.+1))
                 wtot=wtot*w(N,int(iv(ilyru,2)/2.+1))**iv(ilyru,1)
                 d(irougfa(i))=(3.+ra)*ru(irougms(i))
                 d(irougfa(i)+1)=(3.-ra)*ru(irougms(i))
                end do
                call CALFRE(nfa,wl,pq,ir,d,out)
                tauS=tauS+out(1,1)*wtot
                rhoS=rhoS+out(2,1)*wtot
                rhopS=rhopS+out(3,1)*wtot
                Ts=Ts+REAL(out(4,1))*wtot
                Rs=Rs+REAL(out(5,1))*wtot
                R1s=R1s+REAL(out(6,1))*wtot
                As=As+REAL(out(7,1))*wtot
                tauP=tauP+out(1,2)*wtot
                rhoP=rhoP+out(2,2)*wtot
                rhopP=rhopP+out(3,2)*wtot
                Tp=Tp+REAL(out(4,2))*wtot
                Rp=Rp+REAL(out(5,2))*wtot
                R1p=R1p+REAL(out(6,2))*wtot
                Ap=Ap+REAL(out(7,2))*wtot
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
       end do
      end do   
      
*** salvataggio su VOSI

      if(nint(par(54,2)).eq.0) then !misure SF speculari
       NQ=IR(NFA)*IR(NFA)
       MUS=SQRT(NQ-PQ)
       MUP=NQ/MUS
       nq1=ir(1)*ir(1)
       mus1=SQRT(nq1-pq)
       mup1=nq1/mus1
       vosi(1,1)=real(mus)/real(mus1)*abs(tauS)**2.
       vosi(1,2)=real(mup)/real(mup1)*abs(tauP)**2.
       vosi(2,1)=ABS(rhoS)**2.
       vosi(2,2)=ABS(rhoP)**2.
       vosi(3,1)=ABS(rhopS)**2.
       vosi(3,2)=ABS(rhopP)**2.
      else !misure SF emisferiche
       vosi(1,1)=Ts
       vosi(1,2)=Tp
       vosi(2,1)=Rs
       vosi(2,2)=Rp
       vosi(3,1)=R1s
       vosi(3,2)=R1p
      end if
      vosi(4,1)=As
      vosi(4,2)=Ap
      vosi(5,1)=ATAND1(ABS(rhoP/rhoS))
      vosi(5,2)=ATAN2D1(-AIMAG(rhoP/rhoS),-REAL(rhoP/rhoS))

      RETURN
      END
      
      
      SUBROUTINE CALFRE(nfa,wl,pq,ir,d,out)
      REAL d(999),pig,wl
      COMPLEX IR(999),SC(2,2),PC(2,2),PQ,MUS,MUP,NQ,CI,S1,S2,P1,P2,
     +        DELTA,C,S,B,BP,CP,out(8,2),deno,nq1,mus1,mup1,
     +        rhoS,rhopS,tauS,rhoP,rhopP,tauP
      DATA CI,PIG/(0.,1.),3.141592654/
*** Subroutine per il calcolo della matrice carateristica e dei
**             coefficienti di Fresnel
      
*** inizializzazione SC PC
      SC(1,1)=CMPLX(1.,0.)
      SC(2,2)=CMPLX(1.,0.)
      PC(1,1)=CMPLX(1.,0.)
      PC(2,2)=CMPLX(1.,0.)
      SC(1,2)=CMPLX(0.,0.)
      SC(2,1)=CMPLX(0.,0.)
      PC(1,2)=CMPLX(0.,0.)
      PC(2,1)=CMPLX(0.,0.)          
*** calcolo matrice caratteristica
      DO I=2,NFA-1
       NQ=IR(I)*IR(I)
       MUS=SQRT(NQ-PQ)
       MUP=NQ/MUS
       DELTA=2.*PIG*D(I)*MUS/wl
       C=COS(DELTA)
       S=SIN(DELTA)*CI
       S1=SC(1,1)
       S2=SC(2,2)
       P1=PC(1,1)
       P2=PC(2,2)
       SC(1,1)=C*S1+S*MUS*SC(1,2)
       SC(2,2)=C*S2+S/MUS*SC(2,1)
       PC(1,1)=C*P1+S*MUP*PC(1,2)
       PC(2,2)=C*P2+S/MUP*PC(2,1)
       SC(1,2)=S/MUS*S1+C*SC(1,2)
       SC(2,1)=S*MUS*S2+C*SC(2,1)
       PC(1,2)=S/MUP*P1+C*PC(1,2)
       PC(2,1)=S*MUP*P2+C*PC(2,1)
      END DO
      NQ=IR(NFA)*IR(NFA)
      MUS=SQRT(NQ-PQ)
      MUP=NQ/MUS
      nq1=ir(1)*ir(1)
      mus1=SQRT(nq1-pq)
      mup1=nq1/mus1
C ***
C *** CALCOLO COEFFICIENTI DI FRESNEL e R R1 T A
C ***
*** polarizzazione s
      B=SC(1,1)+SC(1,2)*MUS
      C=SC(2,1)+SC(2,2)*MUS
      BP=-SC(1,1)+SC(1,2)*MUS
      CP=-SC(2,1)+SC(2,2)*MUS
      deno=mus1*B+C
      rhoS=(mus1*B-C)/deno
      rhopS=(mus1*BP+CP)/deno
      tauS=2.*mus1/deno
      Rs=ABS(rhoS)**2.
      R1s=ABS(rhopS)**2.
      Ts=4.*REAL(mus1)*REAL(mus)/ABS(deno)**2.
      As=4.*REAL(mus1)*REAL(B*CONJG(C)-mus)/ABS(deno)**2.
       
*** polarizzazione p
      B=PC(1,1)+PC(1,2)*MUP
      C=PC(2,1)+PC(2,2)*MUP
      BP=-PC(1,1)+PC(1,2)*MUP
      CP=-PC(2,1)+PC(2,2)*MUP
      deno=mup1*B+C
      rhoP=(mup1*B-C)/deno
      rhopP=(mup1*BP+CP)/deno
      tauP=2.*mup1/deno
      Rp=ABS(rhoP)**2.
      R1p=ABS(rhopP)**2.
      Tp=4.*REAL(mup1)*REAL(mup)/ABS(deno)**2.
      Ap=4.*REAL(mup1)*REAL(B*CONJG(C)-mup)/ABS(deno)**2.
       
*** PSI e DELTA
      PSI=ATAND1(ABS(rhoP/rhoS))
      DEL=ATAN2D1(-AIMAG(rhoP/rhoS),-REAL(rhoP/rhoS))
      
*** salvataggio su OUT
      out(1,1)=tauS
      out(2,1)=rhoS
      out(3,1)=rhopS
      out(4,1)=Ts
      out(5,1)=Rs
      out(6,1)=R1s
      out(7,1)=As
      out(8,1)=PSI
      out(1,2)=tauP
      out(2,2)=rhoP
      out(3,2)=rhopP
      out(4,2)=Tp
      out(5,2)=Rp
      out(6,2)=R1p
      out(7,2)=Ap
      out(8,2)=DEL
      
      RETURN
      END
      
      
!       SUBROUTINE MODELWIZARD
!       INTEGER IXW(20)
!       REAL PM(200,5),PAR(60,5),CNK(16,3),PF(7,21),M(17,201,2),
!      *     E(8,201,2),RXY(25,4),SOL(999,6),ARSE(500,2)
!       CHARACTER NANK(16)*256,SOUT*19,st1*1,st6*6,r*1
!       COMMON /VALORI/ M,E,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
!       COMMON /STRINGHE/ NANK
! *** settaggio guidato del modello
!       
! **    assegnazione file.nk a CNK
!       do i=8,15
!        if(nank(i-7).ne.'mate/aa999.9') then
!         cnk(i-6,1)=i
!        end if
!       end do
!       
!       write(*,*)'********************* MODEL WIZARD *******************'
!       Write(*,*)'Inizio procedura guidata impostazione modello ...'
!  1    write(*,*)'Numero strati (1-9)?'
!       read(*,*,ERR=1) par(51,2)
!       if(nint(par(51,2)).lt.1) par(51,2)=1.
!       if(nint(par(51,2)).gt.9) par(51,2)=9.
!       do i=1,nint(par(51,2))
!        write(*,*)'   strato N.',i,':'
!  2     write(*,*)'tipo (1->bulk, 2->film_omo, 3->film_inomo) e spessore
!      +(bulk mm; film Angstrom)?'
!        read(*,*,ERR=2) itype,pm(i,1)
!        if(itype.lt.1) itype=1
!        if(itype.gt.5) itype=5
!        par(50+i,3)=itype
!        if(itype.eq.1) pm(i,1)=pm(i,1)*1.e+7
!        write(*,*)'<<<<<<<<< Selezione indice del materiale >>>>>>>>>'
!        write(*,*)'1)  x '
!        do jj=2,9
!         write(st1,'(i1.1)') jj
!         st6=st1//')noto'
!         call gecnk(jj,'s',st6,sout)
!        end do
!  3     write(*,*)'Quale indice scegli (dopo con M# puoi modificare il ma 
!      +teriale)?'
!        read(*,*,ERR=3) par(50+i,1)
!        if(par(50+i,1).lt.1..or.par(50+i,1).gt.9.) goto 3
!        if(nint(par(50+i,1)).eq.1) par(53,2)=i !strato incognito
!       end do
!       write(*,*)'multistrato simmetrico su faccia posteriore (y o ret)?'
!       read(*,'(a1)') r
!       par(52,2)=.0
!       if(r.eq.'y') par(52,2)=1.
!       write(*,*)'... procedura completata!!!'
!       write(*,*)'"m" per accedere al pannello di controllo completo' 
!       write(*,*)'*****************************************************'
!       RETURN
!       END
!       
!       
!       SUBROUTINE MODEL(ic)
!       INTEGER IXW(20)
!       REAL PM(200,5),PAR(60,5),CNK(16,3),PF(7,21),M(17,201,2),
!      *     E(8,201,2),RXY(25,4),SOL(999,6),ARSE(500,2)
!       CHARACTER R*1,R1*1,R2*2,NANK(16)*256,st1*1,st6*6,
!      +          st2*2,TY(5)*5,SOUT*19,SINP*6
!       COMMON /VALORI/ M,E,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
!       COMMON /STRINGHE/ NANK
!       DATA TY(1)/'bulk '/, TY(2)/'f-omo'/, TY(3)/'f-ino'/,
!      +TY(4)/'fexpL'/, TY(5)/'fexpR'/
! c
! c   ic = -1 -> nessuna richiesta
! c         0 -> richieste complete
! c         # -> vai direttamente a M# e poi lascia senza conferma
! c
!       
! *** parametri multistrato
!       nlayer=nint(par(51,2))
!       
! *** controllo N per il calcolo della rugosità
!       if(par(28,1).lt.1.) par(28,1)=1.
!       if(par(28,1).gt.5.) par(28,1)=5.
! 
!       
!       R='?'
!       do while(R.ne.' ')
!  17    write(*,*)
!      +'************************* MULTISTRATO GENERALIZZATO *************
!      +************'
!        write(*,'(a60)')
!      +'layer type  d(mm-A) rough(A)   Gn       Cn       Gk       Ck    
!      +nk'
!        write(*,'(a68)')
!      +' "#"  "#t"    "#d"    "#r"    "#g"     "#u"     "#G"     "#U"  "#  
!      +m"'
!        do i=1,nlayer
!         call gecnk(nint(par(50+i,1)),'?','      ',sout)
!         if(nint(par(50+i,3)).eq.1) then !strato incoerente
!          write(*,FMT=100) ' "',i,'" ',TY(nint(par(50+i,3))),
!      +            PM(i,1)*1.E-7,PM(50+i,1),nint(par(50+i,1)),sout(1:15)
!  100     FORMAT(a2,i1,a2,a5,f9.2,2x,f6.1,38x,i1,1x,a15)
!         else if(nint(par(50+i,3)).eq.2) then!strato coerente omogeneo
!          write(*,FMT=101) ' "',i,'" ',TY(nint(par(50+i,3))),PM(i,1),
!      +         PM(50+i,1),nint(par(50+i,1)),sout(1:15)
!  101     FORMAT(a2,i1,a2,a5,f9.2,2x,f6.1,38x,i1,1x,a15)
!         else if(nint(par(50+i,3)).eq.3) then!strato coerente inomogeneo
!          write(*,FMT=102) ' "',i,'" ',TY(nint(par(50+i,3))),PM(i,1),
!      +         PM(50+i,1),PM(10+i,1),PM(20+i,1),PM(30+i,1),PM(40+i,1),
!      +         nint(par(50+i,1)),sout(1:15)
!  102     FORMAT(a2,i1,a2,a5,f9.2,2x,f6.1,4(2x,f7.4),2x,i1,1x,a15)
!         else if(nint(par(50+i,3)).eq.4) then!strato coerente inomogeneo exp\_
!          write(*,FMT=103) ' "',i,'" ',TY(nint(par(50+i,3))),PM(i,1),
!      +         PM(50+i,1),PM(10+i,1),PM(20+i,1),PM(30+i,1),PM(40+i,1),
!      +         nint(par(50+i,1)),sout(1:15)
!  103     FORMAT(a2,i1,a2,a5,f9.2,2x,f6.1,4(2x,f7.4),2x,i1,1x,a15)
!          else if(nint(par(50+i,3)).eq.4) then!strato coerente inomogeneo exp\_
!          write(*,FMT=104) ' "',i,'" ',TY(nint(par(50+i,3))),PM(i,1),
!      +         PM(50+i,1),PM(10+i,1),PM(20+i,1),PM(30+i,1),PM(40+i,1),
!      +         nint(par(50+i,1)),sout(1:15)
!  104     FORMAT(a2,i1,a2,a5,f9.2,2x,f6.1,4(2x,f7.4),2x,i1,1x,a15)       
!         end if
!        end do
!        if(nint(par(54,2)).eq.0) then
!         write(*,*)'"S" R&T speculari'
!        else
!         write(*,*)'"S" R&T emisferiche'
!        end if
!        if(nint(par(52,2)).eq.0) then
!         write(*,*)
!      +        '"s" faccia posteriore interfacciata al mezzo di uscita'
!        else
!         write(*,*)'"s" multistrato simmetrico sulla faccia posteriore!'
!        end if  
!        write(*,*)
!      +'*****************************************************************
!      +************'     
!        write(*,*)'Comando - riga da modificare (h = lista comandi)?'
!        read(*,'(a2)') st2
!        r=st2(1:1)
!        
!        if(r.eq.'h') then
!         write(*,*)'"n"  -> numero layers'
!         write(*,*)'"N"  -> calcolo di <f(x)>, N =',nint(par(28,1))
!         write(*,*)'"z"  -> lista mezzi ingresso e uscita'
!         write(*,*)'"i"  -> profilo parabolico, Nfilm  =',nint(par(29,1))
!         write(*,*)'"#+" -> inserisci uno strato al posto #'
!         write(*,*)'"#-" -> elimina lo strato #'
!         write(*,*)'"M#" -> modifica il materiale #'
!         
!        else if(r.eq.'S') then
!         par(54,2)=1.-par(54,2)
!            
!        else if(r.eq.'s') then
!         par(52,2)=1.-par(52,2)
!         
!        else if(r.eq.'M') then
!         read(st2(2:2),'(i1)',ERR=5) j
!         if(j.ge.2.and.j.le.9) then
!          sinp(1:5)='mat. '
!          sinp(6:6)=st2(2:2)
!          call gecnk(j,'y',SINP,SOUT)
!         end if
!  5      continue          
!        else if(r.eq.'n') then
!  14     write(*,*)'Nuovo numero strati(max 9)?'
!         read(*,*,ERR=14) nlayer
!         if(nlayer.lt.0.or.nlayer.gt.9) goto 14
!         par(51,2)=nlayer
!                      
!        else if(r.eq.'N') then
!  21     write(*,*)'Nuovo valore di N (1<= N <=5)?'
!         read(*,*,ERR=21) par(28,1)
!         if(par(28,1).lt.1.) par(28,1)=1.
!         if(par(28,1).gt.5.) par(28,1)=5.
!            
!        else if(r.eq.'i') then
!  20     write(*,*)'Nuovo valore di N (1<= N <=99)?'
!         read(*,*,ERR=20) par(29,1)
!         if(par(29,1).lt.1.) par(29,1)=1.
!         if(par(29,1).gt.99.) par(29,1)=99.
!         
!        else if(r.eq.'z') then
!         do while(r2.ne.'  ')
!          CALL GECNK(10,'s','i)inSF',sout)
!          CALL GECNK(11,'s','a)Apds',sout)
!          CALL GECNK(12,'s','1)EL1 ',sout)
!          CALL GECNK(13,'s','2)EL2 ',sout)
!          CALL GECNK(14,'s','3)EL3 ',sout)
!          CALL GECNK(15,'s','4)EL4 ',sout)
!          CALL GECNK(16,'s','o)ouSF',sout)
!          read(*,'(a2)')r2
!          if(r2.eq.'i ') CALL GECNK(10,'t','inp-SF',sout)
!          if(r2.eq.'a ') CALL GECNK(11,'t','i-Apds',sout)
!          if(r2.eq.'1 ') CALL GECNK(12,'t','in-EL1',sout)
!          if(r2.eq.'2 ') CALL GECNK(13,'t','in-EL2',sout)
!          if(r2.eq.'3 ') CALL GECNK(14,'t','in-EL3',sout)
!          if(r2.eq.'4 ') CALL GECNK(15,'t','in-EL4',sout)
!          if(r2.eq.'o ') CALL GECNK(16,'t','out-SF',sout)
!         end do
!         
!        else if(r.ne.' ') then
!         read(st2(1:1),*,ERR=17) j
!         r1=st2(2:2)
!         
!         if(r1.eq.'t') then
!  4       write(*,*)'Nuovo tipo (1->bulk, 2->film_omo, 3->film_inomo, 4->
!      *film_expL, 5->film_expR)?'
!          read(*,*,ERR=4) itype
!          if(itype.lt.1) itype=1
!          if(itype.gt.5) itype=5
!          par(50+j,3)=itype
!          
!         else if(r1.eq.'d') then
!  15      write(*,*)'Nuovo valore spessore(bulk mm, film A)?'
!          read(*,*,ERR=15) PM(j,1)
!          if(nint(par(50+j,3)).eq.1) pm(j,1)=pm(j,1)*1.e+7
!          
!         else if(r1.eq.'r') then
!  19      write(*,*)'Nuovo valore sigma_roughness (A)?'
!          read(*,*,ERR=19) PM(50+j,1)
!          if(6.*pm(50+j,1).gt.pm(j,1)) then
!           write(*,*)'ATTENZIONE sigma-rough <= spessore !!!!'
!           goto 19
!          end if
!          
!         else if(r1.eq.'g') then
!  18      write(*,*)'Nuovo valore Dn/<n> (gradiente di n)?'
!          read(*,*,ERR=18) PM(10+j,1)
!          
!         else if(r1.eq.'u') then
!  1       write(*,*)'Nuovo valore nav/<n>-1 (curvatura di n)?'
!          read(*,*,ERR=1) PM(20+j,1)
!          
!         else if(r1.eq.'G') then
!  2       write(*,*)'Nuovo valore Dk/<k> (gradiente di k)?'
!          read(*,*,ERR=2) PM(30+j,1)
!          
!         else if(r1.eq.'U') then
!  3       write(*,*)'Nuovo valore kav/<k>-1 (cruvatura di k)?'
!          read(*,*,ERR=3) PM(40+j,1)
!          
!         elseif(r1.eq.'m') then
!          if(nint(par(53,2)).eq.j) par(53,2)=0. !azzera il puntatore
!          write(*,*)'<<<<<<<<< Selezione indice del materiale >>>>>>>>>'
!          write(*,*)'1)  x '
!          do jj=2,9
!           write(st1,'(i1.1)') jj
!           st6=st1//')noto'
!           call gecnk(jj,'s',st6,sout)
!          end do
!  16      write(*,*)'Quale indice scegli (dopo con M# puoi modificare il 
!      +materiale)?'
!          read(*,*,ERR=16) par(50+j,1)
!          if(par(50+j,1).lt.1..or.par(50+j,1).gt.9.) goto 16
!          if(nint(par(50+j,1)).eq.1) par(53,2)=j !strato incognito
!          
!         elseif(r1.eq.'+') then
!          if(j.ge.1.and.j.le.nlayer) then
!           do jj=nlayer,j,-1
!            do ij=1,6
!             do jjj=1,3
!              pm(10*(ij-1)+jj+1,jjj)=pm(10*(ij-1)+jj,jjj)
!             end do
!            end do
!            par(50+jj+1,1)=par(50+jj,1) !puntatore VNK
!            par(50+jj+1,3)=par(50+jj,3) !Tipologia
!           end do
!           nlayer=nlayer+1
!          end if
!          
!         elseif(r1.eq.'-') then
!          if(j.ge.1.and.j.le.nlayer) then
!           do jj=j,nlayer,j
!            do ij=1,6
!             do jjj=1,3
!              pm(10*(ij-1)+jj,jjj)=pm(10*(ij-1)+jj+1,jjj)
!             end do
!            end do
!            par(50+jj,1)=par(50+jj+1,1) !puntatore VNK
!            par(50+jj,3)=par(50+jj+1,3) !Tipologia
!           end do
!           nlayer=nlayer-1
!          end if
!          
!         end if
!            
!        end if
!        
! *** controllo selezione materiali
!        do j=1,nlayer
!         i=nint(par(50+j,1))
!         if(i.lt.1.or.i.gt.9) then
!          write(*,*)'Scelta indice del materiale del layer',j
!          if(nint(par(53,2)).eq.j) par(53,2)=0. !azzera il puntatore
!          write(*,*)'1)  x '
!          do jj=2,9
!           write(st1,'(i1.1)') jj
!           st6=st1//')noto'
!           call gecnk(jj,'s',st6,sout)
!          end do
!   6      write(*,*)'Quale indice scegli(dopo con M# puoi modificare il 
!      +materiale)?'
!          read(*,*,ERR=6) par(50+j,1)
!          if(par(50+j,1).lt.1..or.par(50+j,1).gt.9.) goto 6
!          if(nint(par(50+j,1)).eq.1) par(53,2)=j !strato incognito
!         end if
!        end do
!        
!       end do
!       par(51,2)=nlayer
! 
!       RETURN
!       END
!       
      
      SUBROUTINE GECNK(I,R,CH,SOUT)
      REAL CNK(16,3),PF(7,21),MIS(17,201,2),ELI(8,201,2),PAR(60,5),
     *     RXY(30,4),SOL(999,6),PM(200,5),ARSE(500,2)
      INTEGER OK(18),IXW(20)
      CHARACTER NANK(16)*256,R*1,CH*6,st76*76,SOUT*19,st1,sn*4,sk*4,
     + sf*11
      COMMON /VALORI/ MIS,ELI,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
C     R='y' : stampa e modifica l'opzione prescelta
C     R='s' : stampa l'opzione prescelta
C     R='t' : modifica l'opzione prescelta
C     R='?' : scrivi solo su sout l'opzione prescelta
C     CH    : character*6 simbolo del soggetto di cui si modificano n,k
      
*** controllo CH
      jch=6
      do while(ch(jch:jch).eq.' ')
       jch=jch-1
       if(jch.eq.0) goto 200
      end do
200   continue      
      sout='                   '

      OK(1)=1
      DO J=2,17
       OK(J)=0
      END DO
      OK(18)=1
      if(r.eq.'?') then
       io=nint(cnk(i,1))
       if(io.eq.0) then
        if(i.eq.1) then
         sout=ch(1:jch)//'nk incogniti'
        else
         write(sn,'(f4.2)') cnk(i,2)
         write(sk,'(f4.2)') cnk(i,3)
         sout=ch(1:jch)//'n='//sn//' k='//sk
        end if
       else if(io.ge.1.and.io.le.7) then
        write(st1,'(I1)')io
        sout=ch(1:jch)//'FIT#'//st1
       else if(io.ge.8.and.io.le.15) then
        length=len_trim(NANK(IO-7))
        sout=ch(1:jch)//NANK(IO-7)(1:length)
       else if(io.eq.16) then
        sout=ch(1:jch)//'ibridone'
       else if(io.eq.17) then
        sout=ch(1:jch)//'nk incogniti'
       else
        i1=nint(cnk(i,1)/1000.)
        i2=nint((cnk(i,1)-i1*1000.)/10.)
        f2=cnk(i,1)-i1*1000.-i2*10.
        write(sf,'(i2.2,1x,i2.2,f6.3)') i1,i2,f2
        sout=ch(1:jch)//'EMA '//sf
       end if
      end if
      IF(r.eq.'y'.or.r.eq.'s') THEN
       IO=nint(CNK(I,1))
       IF(IO.EQ.0) then
        WRITE(*,*) CH,' n,k  <-> cte: n=',CNK(I,2),' k=',CNK(I,3)
       else IF(IO.GE.1.AND.IO.LE.7) then
        WRITE(*,*) CH,' n,k  <-> FIT#',IO,' Osc:',
     +   (nint(pf(io,j)),j=2,nint(pf(io,1)+1))
       else IF(IO.GE.8.AND.IO.LE.15) then
        length=len_trim(NANK(IO-7))
        WRITE(*,*) CH,' n,k  <-> file.nk   = ',NANK(IO-7)(1:length)
       else if(io.eq.16) then
        write(*,*) CH,' n,k  <-> calcolati con Ibridone'  
       else
        i1=nint(cnk(i,1)/1000.)
        i2=nint((cnk(i,1)-i1*1000.)/10.)
        f2=cnk(i,1)-i1*1000.-i2*10.
        write(*,'(1x,a6,a6,i2.2,a5,i2.2,a7,f6.3)') 
     +       CH,' EMA "',i1,'" & "',i2,'" f2 = ',f2 
        IF(i1.EQ.0) then
         WRITE(*,'(8x,a1,i2.2,a21,f8.4,a3,e9.3)') '"',i1,
     +      '" =  n,k  <-> cte: n=',CNK(I,2),' k=',CNK(I,3)
        else IF(i1.GE.1.AND.I1.LE.7) then
         WRITE(*,'(8x,a1,i2.2,a17,i1,a5,20(1x,i2),:)')'"',i1,
     +'" = n,k  <-> FIT#',I1,' Osc:',
     +   (nint(pf(i1,j)),j=2,nint(pf(i1,1))+1)
        else IF(I1.GE.8.AND.I1.LE.15) then
        length=len_trim(NANK(I1-7))
         WRITE(*,'(8x,a1,i2.2,a25,a16)') '"',i1,
     +      '" = n,k  <-> file.nk   = ',NANK(I1-7)(1:length)
        else if(i1.eq.17) then
         write(*,'(8x,a1,i2.2,a17)')'"',i1,'" = n,k incogniti'
        end if
        IF(i2.EQ.0) then
         WRITE(*,'(8x,a1,i2.2,a21,f8.4,a3,e9.3)') '"',i2,
     +       '" = n,k  <-> cte: n=',CNK(I,2),' k=',CNK(I,3)
        else IF(i2.GE.1.AND.I2.LE.7) then
         WRITE(*,'(8x,a1,i2.2,a17,i1,a5,20(1x,i2),:)')'"',i2,
     +'" = n,k  <-> FIT#',I2,' Osc:',
     +   (nint(pf(i2,j)),j=2,nint(pf(i2,1))+1)
        else IF(I2.GE.8.AND.I2.LE.15) then
        length=len_trim(NANK(I2-7))
         WRITE(*,'(8x,a1,i2.2,a25,a16)')'"', i2,
     +       '" = n,k  <-> file.nk   = ',NANK(I2-7)(1:length)
        else if(i2.eq.17) then
         write(*,'(8x,a1,i2.2,a17)')'"',i2,'" = n,k incogniti'
        end if
       end if
      END IF
      IF(r.eq.'y'.or.r.eq.'t') THEN
       WRITE(*,*)'************** OPZIONI DISPONIBILI *****************'
       WRITE(*,'(1x,a17)')'"00" n,k <-> cte '
       DO M=1,7
        OK(M+1)=1
        WRITE(*,'(1x,a1,i2.2,a14,i1,a5,20(1x,i2),:)')
     +   '"',M,'" n,k <-> FIT#',M,
     +   ' Osc:',(nint(pf(m,j)),j=2,nint(pf(m,1))+1)
       END DO
       DO M=8,15
        IF(NANK(M-7).NE.'mate/aa999.9') THEN
         OK(M+1)=1
        length=len_trim(NANK(M-7))
         WRITE(*,'(1x,a1,i2.2,a22,a16)')'"',M,'" n,k <-> file.nk   = ',
     +        NANK(M-7)(1:length)
        END IF
       END DO
       if(mis(16,1,1).gt.0.) then
        ok(16+1)=1
        write(*,*)'"16" n,k <-> calcolati con Ibridone'
       end if
       write(*,*)'"17" n,k <-> incognite'
       WRITE(*,*)'***************************************************'
 1     WRITE(*,*)'Opzione o opz#1 opz#2 f2 per EMA (o1o2f.fff)'
       READ (*,*,ERR=1) val
       if(val.ge.0.and.val.le.16.) then 
        IOPT=nint(val)
        IF(OK(IOPT+1).EQ.1) THEN
         CNK(I,1)=IOPT
         IF(IOPT.EQ.0) THEN
          WRITE(*,*)'n = ',CNK(I,2),' k = ',CNK(I,3)
 2        WRITE(*,*)'Nuovi valori o ret!'
          READ(*,'(A76)') st76
          IF(st76(1:1).ne.' ') then
           READ(st76,*,ERR=2,IOSTAT=ios) CNK(I,2),CNK(I,3)
           if(ios.ne.0) goto 2
          end if
         else if(iopt.ge.1.and.iopt.le.7) then
          call OPFITCTRL(iopt,PF,PM)
         END IF
        END IF
       else
        i1=nint(val/1000.)
        i2=nint((val-i1*1000.)/10.)
        f2=val-i1*1000.-i2*10.
        if(ok(i1+1).eq.1.and.ok(i2+1).eq.1.and.(f2.le.1.and.f2.ge..0))
     +   cnk(i,1)=val
       end if
      END IF
      RETURN
      END
      
      
      SUBROUTINE COSVNK(Imin,VNK,L)
      INTEGER IXW(20)
      REAL CNK(16,3),VNK(16,2),PF(7,21),M(17,201,2),E(8,201,2),
     *     par(60,5),RXY(30,4),SOL(999,6),PM(200,5),ARSE(500,2)
      REAL n1,k1,n2,k2
      CHARACTER NANK(16)*256
      COMMON /VALORI/ M,E,PAR,PM,CNK,SOL,PF,RXY,ARSE,IXW
      COMMON /STRINGHE/ NANK
C     Imin=1 quando n,k di lavoro sono noti (ad esempio in simula)
C          2 quando n,k di lavoro sono incogniti
C
C     L = indice wl

      WL=M(7,L,1)
      DO J=Imin,16
       io=nint(CNK(J,1))
       if(io.ge.0.and.io.le.17) then
        call SETVNK(io,J,CNK,VNK,M,PF,PM,L)
       else if(io.ge.0) then
        i1=nint(cnk(j,1)/1000.)
        i2=nint((cnk(j,1)-i1*1000.)/10.)
        f2=cnk(j,1)-i1*1000.-i2*10.
        call SETVNK(i1,J,CNK,VNK,M,PF,PM,L)
        n1=VNK(J,1)
        k1=VNK(J,2)
        call SETVNK(i2,J,CNK,VNK,M,PF,PM,L)
        n2=VNK(J,1)
        k2=VNK(J,2)
        CALL EMA(n1,k1,n2,k2,f2,vnk(j,1),vnk(j,2))
       END IF
      END DO
      RETURN
      END
      
      
      SUBROUTINE SETVNK(io,J,CNK,VNK,M,PF,PM,L)
      REAL CNK(16,3),VNK(16,2),M(17,201,2),PF(7,21),PM(200,5)
*** set valori VNK(J,1) e VNK(J,2) dato CNK(J,k)
c         (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) io    
      GOTO(1,2,2,2,2,2,2,2,3,3, 3, 3, 3, 3, 3, 3, 3, 4) io+1
 1     VNK(J,1)=CNK(J,2)
       VNK(J,2)=CNK(J,3)
      goto 9
 2     WL=M(7,L,1)
       ev=12400./wl
       CALL FDISP(io,ev,PF,PM,VNK(J,1),VNK(J,2)) 
      goto 9
 3     VNK(J,1)=M(io,L,1)
       VNK(J,2)=M(io,L,2)
      goto 9
 4     VNK(J,1)=vnk(1,1)
       VNK(J,2)=vnk(1,2)
 9    RETURN
      END
      
      
      SUBROUTINE EMA(N,K,NA,KA,FA,NE,KE)
      REAL N,K,NA,KA,NE,KE
      COMPLEX EA,EB,EE,E
      IF(FA.LT..0001) THEN
       NE=N
       KE=K
      ELSE
       IF(FA.GE..0001.AND.FA.LE..9999) THEN
        EA=CMPLX(NA*NA-KA*KA,-2*NA*KA)
        EB=CMPLX(N*N-K*K,-2*N*K)
        SEGNO=(1.-FA)*AIMAG(EB)+FA*AIMAG(EA)
        EE=(1.-FA)*(EA-2.*EB)+FA*(EB-2.*EA)
        E=(-EE+CSQRT(EE*EE+8.*EA*EB))/4.
        ER=REAL(E)
        EI=AIMAG(E)
        IF(EI*SEGNO.LT.0.) THEN
         E=(-EE-CSQRT(EE*EE+8.*EA*EB))/4.
         ER=REAL(E)
         EI=AIMAG(E)
        END IF
        KE=SQRT((-ER+SQRT(ER*ER+EI*EI))*.5)*SIGN(1.,-EI)
        NE=SQRT(ER+KE*KE)
       ELSE
        NE=NA
        KE=KA
       END IF
      END IF
      RETURN
      END      
      
      
      SUBROUTINE OPFITCTRL(iopt,PF,PM)
      REAL PF(7,21),PM(200,5)
      CHARACTER r*1,sfu(6)*12,r3*3
      sfu(1)='classico'
      sfu(2)='quant. omo'
      sfu(3)='quant. inomo' 
      sfu(4)='flat'
      sfu(5)='drude'
      sfu(6)='free carrier'
      r='?'
      do while(r.ne.' ')
       noscFT=nint(pf(iopt,1))
       write(*,*)'****** Opzione fit = ',iopt,'***********************'
       write(*,*) ' N.  Osc.    Tipo       C         E         D      K'
       do io=1,noscFT
        i=nint(pf(iopt,io+1))
        ifu=nint(pm(101+(I-1)*5,1))
        if(ifu.lt.1.or.ifu.gt.6) then
         ifu=1
         pm(101+(I-1)*5,1)=ifu
        end if
        write(*,FMT=100) io,i,sfu(ifu),(pm(101+(I-1)*5+j,1),j=1,4)
       end do
       write(*,*)'****************************************************'
       write(*,*)'"a " -> aggiungi un termine'
       write(*,*)'"k#" -> elimina il termine #'
       write(*,*)'"o#" -> modifica il tipo del termine #'
       write(*,*)'"X#" -> modifica il parametro "X" del termine #'
       write(*,*)'Che fare ?'
       read(*,'(a3)') r3
       r=r3(1:1)

       if(r.eq.'a') then
        if(noscFT.le.19) then
         noscFT=noscFT+1
         pf(iopt,1)=noscFT
         write(*,*)'N.osc   stato         tipo: '
         do ii=1,20
          ifu=nint(pm(101+(ii-1)*5,1))
          if(ifu.ge.1.and.ifu.le.6) then
           write(*,'(1x,i2,2x,a9,2x,a24)') ii,'esistente',sfu(ifu)
          else
           write(*,'(1x,i2,2x,a6)') ii,'libero'
          end if
         end do
 1       write(*,*)'Numero oscillatore da associare al nuovo termine?'
         read(*,'(i2)',ERR=1) i
         if(i.ge.1.and.i.le.20) then
          pf(iopt,noscFT+1)=i
          ifu=nint(pm(101+(i-1)*5,1))
          if(ifu.lt.1.or.ifu.gt.6) then
 3         write(*,*)'Tipo (1=Class, 2=Qomo, 3=Qinomo, 4=Flat, 5=Drude, 
     + 6=FreeCarriers) ?'
           read(*,'(i2)',ERR=3) ifu
           if(ifu.ge.1.and.ifu.le.6) then
            pm(101+(I-1)*5,1)=ifu
           else
            goto 3
           end if
          end if
         else
          goto 1
         end if
        end if
       end if
       
       if(r.eq.'k') then
        read(r3(2:3),'(i2)',ERR=10) io
        if(io.le.noscFT) then
         noscFT=noscFT-1         
         pf(iopt,1)=noscFT
         pf(iopt,io+1)=0
        end if
10      continue
       end if
       
       if(r.eq.'o') then
        read(r3(2:3),'(i2)',ERR=102) io
        i=nint(pf(iopt,io+1))
        pm(101+(I-1)*5,1)=pm(101+(I-1)*5,1)+1.
        if(pm(101+(I-1)*5,1).ge.5.) pm(101+(I-1)*5,1)=1.
102     continue
       end if
       
       if(r.eq.'C'.or.r.eq.'E'.or.r.eq.'D'.or.r.eq.'K') then
        read(r3(2:3),'(i2)',ERR=101) io
        i=nint(pf(iopt,io+1))
        if(r.eq.'C') then
         j=1
        else if(r.eq.'E') then
         j=2
        else if(r3(1:1).eq.'D') then
         j=3
        else if(r3(1:1).eq.'K') then
         j=4
        else
         j=10
        end if
        if(j.lt.10) then
 2       write(*,*)'Nuovo valore?'
         read(*,*,ERR=2) pm(101+(I-1)*5+j,1)
        end if
101     continue
       end if
       
      end do
 100  FORMAT(1x,i2,2x,i2.2,2x,a12,2x,f8.4,2x,f8.4,2x,f8.4,2x,e9.3)
      RETURN
      END
      
      
      SUBROUTINE FDISP(iopt,eV,PF,PM,sqn,sqk)
c      COMPLEX b,dc,xmin,xmax,aim
      REAL PM(200,5),PF(7,21),K,cln2,sqrln2
      REAL I2CODYRE,ICCODYRE,I2CODYIM,ICCODYIM
      REAL I2TAUCRE,ICTAUCRE,I2TAUCIM,ICTAUCIM
      DATA cln2,sqrln2 /0.693147181,0.832554611/
      
*** epsilon e nk da sommatoria oscillatori
      noscFT=nint(pf(iopt,1))
      E=eV
      epr=0.
      epi=0.
      epik=0.
      do io=2,noscFT+1
       i=nint(pf(iopt,io))
       ifu=nint(pm(101+(I-1)*5,1))
       C=abs(pm(102+(I-1)*5,1))
       E0=pm(103+(I-1)*5,1)
       D=pm(104+(I-1)*5,1)
       K=pm(105+(I-1)*5,1)
       goto(1,2,3,4,5,6,7,8,9) ifu
 1      den=((E0**2.-E**2.)**2.+(E*D)**2.) !oscillatore classico
        eprj=C*(E0**2.-E**2.)/den
        epij= K*C*E*D/den
       goto 88
 2      x=(E0-E)/D !oscillatore quantistico omogeneo
        eprj=C*x/(1+x**2.)
        epij= K*C/(1+x**2.)
       goto 88
 3      x=(E0-E)/D !oscillatore quantistico inomogeneo
        eprj=C*DAWS(sqrln2*x)
        epij= K*C*Exp(-cln2*x**2.)
       goto 88
 4      eprj=C*C !solo in questo caso C=n_cte ! n=costante
        epij=0.
       goto 88
 5      eprj=-E0*E0/(D*D+E*E) !drude model
        epij= K*E0*E0*D/(E**3.+E*D*D)
       goto 88
 6      De=K !indirectGapCody
        eprj=C*(2.*ICCODYRE(E,E0+3.*De/4.,D)-2.*ICCODYRE(E,E0+De/4.,D)+
     +   16./De**2.*(I2CODYRE(E,E0,E0+De/4.,D)-I2CODYRE(E,E0,E0,D)-
     +   I2CODYRE(E,E0+De/2.,E0+3.*De/4.,D)+
     +   I2CODYRE(E,E0+De/2.,E0+De/4.,D)+I2CODYRE(E,E0+De,E0+De,D)-
     +   I2CODYRE(E,E0+De,E0+3.*De/4.,D)))
        
        epij=C*(2.*ICCODYIM(E,E0+3.*De/4.,D)-2.*ICCODYIM(E,E0+De/4.,D)+
     +   16./De**2.*(I2CODYIM(E,E0,E0+De/4.,D)-I2CODYIM(E,E0,E0,D)-
     +   I2CODYIM(E,E0+De/2.,E0+3.*De/4.,D)+
     +   I2CODYIM(E,E0+De/2.,E0+De/4.,D)+I2CODYIM(E,E0+De,E0+De,D)-
     +   I2CODYIM(E,E0+De,E0+3.*De/4.,D)))
       goto 88
 7      De=K !indirectGapTauc
        eprj=C*(2.*ICTAUCRE(E,E0+3.*De/4.,D)-2.*ICTAUCRE(E,E0+De/4.,D)+
     +   16./De**2.*(I2TAUCRE(E,E0,E0+De/4.,D)-I2TAUCRE(E,E0,E0,D)-
     +   I2TAUCRE(E,E0+De/2.,E0+3.*De/4.,D)+
     +   I2TAUCRE(E,E0+De/2.,E0+De/4.,D)+I2TAUCRE(E,E0+De,E0+De,D)-
     +   I2TAUCRE(E,E0+De,E0+3.*De/4.,D)))
        
        epij=C*(2.*ICTAUCIM(E,E0+3.*De/4.,D)-2.*ICTAUCIM(E,E0+De/4.,D)+
     +   16./De**2.*(I2TAUCIM(E,E0,E0+De/4.,D)-I2TAUCIM(E,E0,E0,D)-
     +   I2TAUCIM(E,E0+De/2.,E0+3.*De/4.,D)+
     +   I2TAUCIM(E,E0+De/2.,E0+De/4.,D)+I2TAUCIM(E,E0+De,E0+De,D)-
     +   I2TAUCIM(E,E0+De,E0+3.*De/4.,D)))
     
       goto 88
 8      E3=K+E0 !directGapCody
        eprj=C*ReM0M3(E,E0,E3,D)
        epij=C*aImM0M3(E,E0,E3,D)
       goto 88 
 9      E3=K+E0 !directGapTauc
        eprj=C*ReTaucL(E,E0,D,E3)
        epij=C*aImTaucL(E,E0,D,E3)
 88    epr=epr+eprj
       epi=epi+epij
      end do
      sqn=sqrt(sqrt(epr*epr+epi*epi)+epr)/sqrt(2.)
      sqk=sqrt(sqrt(epr*epr+epi*epi)-epr)/sqrt(2.)
      
      RETURN
      END
      
      FUNCTION I2CODYRE(E,E0,Ep,D)
      REAL I2CODYRE
      DATA PIG /3.141592654/
      I2CODYRE=0.5/PIG*((Ep-E)*(3.*E-4.*E0+Ep)+4.*D*(E-E0)*
     + atan2(E-Ep,D)+((E-E0)**2.-D*D)*log(D*D+(E-Ep)**2.))
      RETURN
      END
      
      FUNCTION ICCODYRE(E,Ep,D)
      REAL ICCODYRE
      DATA PIG /3.141592654/
      ICCODYRE=0.5/PIG*(log(D*D+(E-Ep)**2.))
      RETURN
      END
      
      
      FUNCTION I2CODYIM(E,E0,Ep,D)
      REAL I2CODYIM
      DATA PIG /3.141592654/
      I2CODYIM=1./PIG*((D*D-(E-E0)**2.)*
     + atan2(E-Ep,D)+D*(Ep-E+(E-E0)*log(D*D+(E-Ep)**2.)))
      RETURN
      END
      
      FUNCTION ICCODYIM(E,Ep,D)
      REAL ICCODYIM
      DATA PIG /3.141592654/
      ICCODYIM=-1./PIG*atan2(E-Ep,D)
      RETURN
      END
      
      FUNCTION ICTAUCRE(E,Ep,D)
      REAL ICTAUCRE
      DATA PIG /3.141592654/
      ICTAUCRE=(-4.*D*E*atan2(E-Ep,D)+2.*E*(D*D+E*E)/Ep+(E*E-D*D)*
     + log((D*D+(E-Ep)**2.)/(Ep*Ep)))/(2.*PIG*(D*D+E*E)**2.)
      RETURN
      END
      
      FUNCTION I2TAUCRE(E,E0,Ep,D)
      REAL I2TAUCRE
      DATA PIG /3.141592654/
      I2TAUCRE=(4.*D*E0*(D*D+E*E-E*E0)*atan2(E-Ep,D)+(D**4.+
     + E*E*(E-E0)**2.+D*D*(2.*E*E-2.*E*E0-E0*E0))*log(D*D+(E-Ep)**2.)+
     + 2.*E0*(E*(D*D+E*E)*E0/Ep+(E*E*(2.*E-E0)+D*D*(2.*E+E0))*log(Ep)))/
     + (2.*PIG*(D*D+E*E)**2.)
      RETURN
      END
      
      FUNCTION ICTAUCIM(E,Ep,D)
      REAL ICTAUCIM
      DATA PIG /3.141592654/
      ICTAUCIM=-((E*E-D*D)*atan2(E-Ep,D)+D*((D*D+E*E)/Ep+
     + E*log((D*D+(E-Ep)**2.)/(Ep*Ep))))/(PIG*(D*D+E*E)**2.)
      RETURN
      END
      
      FUNCTION I2TAUCIM(E,E0,Ep,D)
      REAL I2TAUCIM
      DATA PIG /3.141592654/
      I2TAUCIM=-((D**4.+E*E*(E-E0)**2.+D*D*(2.*E*E-2.*E*E0-E0*E0))*
     + atan2(E-Ep,D)+D*E0*(E0*(D*D+E*E)/Ep+(D*D+E*E-E*E0)*
     + log(Ep*Ep/(D*D+(E-Ep)**2.))))/(PIG*(D*D+E*E)**2.)
      RETURN
      END
      
      FUNCTION DAWS(x)
      INTEGER NMAX
      REAL daws,x,H,A1,A2,A3
      PARAMETER (NMAX=6,H=0.4,A1=2./3.,A2=0.4,A3=2./7.)
      INTEGER i,init,n0
      REAL d1,d2,e1,e2,sum,x2,xp,xx,c(NMAX)
      SAVE init,c
      DATA init/0/  ! Flag is 0 if we need to initialize, else 1
      if(init.eq.0)then
       init=1
       do i=1,NMAX
        c(i)=exp(-((2.*float(i)-1.)*H)**2.)
       enddo
      endif
      if(abs(x).lt.0.2)then ! Use series expansion
       x2=x**2.
       DAWS=x*(1.-A1*x2*(1.-A2*x2*(1.-A3*x2)))
      else ! Use sampling theorem representation
       xx=abs(x)
       n0=2*nint(0.5*xx/H)
       xp=xx-float(n0)*H
       e1=exp(2.*xp*H)
       e2=e1**2.
       d1=float(n0+1)
       d2=d1-2.
       sum=0.
       do i=1,NMAX
        sum=sum+c(i)*(e1/d1+1./(d2*e1))
        d1=d1+2.
        d2=d2-2.
        e1=e2*e1
       end do                                    
       DAWS=0.5641895835*sign(exp(-xp**2.),x)*sum ! Constant is 1/sqrt(pig)
      end if
      RETURN
      END
      
      
      FUNCTION ReFdirGap(E,E0,D,x)
      COMPLEX b,dc,xc,aim
      b=CMPLX(E-E0,0.)
      dc=CMPLX(D,0.)
      xc=CMPLX(x,0.)
      aim=CMPLX(0.,1.)
      ReFdirGap=REAL(
     +   2.*sqrt(dc*xc+b)-sqrt(aim*dc-b)*atan(sqrt((dc*xc+b)/
     +  (aim*dc-b)))-sqrt(b+aim*dc)*atanh(sqrt((dc*xc+b)/(aim*dc+b))))
      RETURN
      END
      
      
      FUNCTION aImFdirGap(E,E0,D,x)
      COMPLEX b,dc,xc,aim
      b=CMPLX(E-E0,0.)
      dc=CMPLX(D,0.)
      xc=CMPLX(x,0.)
      aim=CMPLX(0.,1.)
      aImFdirGap=REAL(
     +   aim*(-sqrt(aim*dc-b)*atan(sqrt((dc*xc+b)/
     +  (aim*dc-b)))+sqrt(b+aim*dc)*atanh(sqrt((dc*xc+b)/(aim*dc+b)))))
      RETURN
      END
      
      
      FUNCTION ReM0M3(E,E0,E3,reD)
      COMPLEX aim,a,b,d,pi,ctmp
      a=CMPLX((E0-E)/reD,0.)
      b=CMPLX((E3-E)/reD,0.)
      d=CMPLX(reD,0.)
      aim=CMPLX(0.,1.)
      pi=CMPLX(acos(-1.),0.)
      ctmp=0.5*d*pi*(a+b-2.*REAL(sqrt(a+aim)*sqrt(b+aim))) 
      ReM0M3=REAL(ctmp)
      RETURN
      END
      
      
      FUNCTION aImM0M3(E,E0,E3,reD)
      COMPLEX aim,a,b,d,ctmp
      a=CMPLX((E0-E)/reD,0.)
      b=CMPLX((E3-E)/reD,0.)
      d=CMPLX(reD,0.)
      aim=CMPLX(0.,1.)
      pi=CMPLX(acos(-1.),0.)
      ctmp=-pi*d*(1.-IMAG(sqrt(a+aim)*sqrt(b+aim)))
      aImM0M3=REAL(ctmp)
      RETURN
      END
      
      FUNCTION ReTaucL(E,E0,D,E3)
      COMPLEX aim,cE,cE0,cD,cE3,cTmp
      aim=CMPLX(0.,1.)
      cE =CMPLX(E ,0.)
      cE0=CMPLX(E0,0.)
      cD =CMPLX(D ,0.)
      cE3=CMPLX(E3,0.)
      cTmp=
     + 0.5/(cD-aim*cE)*(2.*aim*cD+2.*cE-2.*sqrt(cE0*cE3)-(aim-1.)*
     + sqrt(2.)*sqrt(cD-aim*(cE-cE0))*sqrt(-aim*cD-cE+cE3) )
      cTmp=cTmp+
     + 0.5/(cD+aim*cE)*(-2.*aim*cD+2.*cE-2.*sqrt(cE0*cE3)+(aim+1.)*
     + sqrt(2.)*sqrt(cD+aim*(cE-cE0))*sqrt(aim*cD-cE+cE3) )
      cTmp=cTmp+
     + cE*(cD-aim*cE)**2.*sqrt(cD+aim*(cE-cE0))*sqrt(cD+aim*(cE-cE3))/
     + (cD*cD+cE*cE)**2.
      cTmp=cTmp+
     + cE*(cD+aim*cE)**2.*sqrt(cD-aim*(cE-cE0))*sqrt(cD-aim*(cE-cE3))/
     + (cD*cD+cE*cE)**2.
      cTmp=cTmp-
     + cE*cD*(cD*cD*(cE0+cE3)+cE*(-4.*cE0*cE3+cE*(cE0+cE3)))/
     + ((cD*cD+cE*cE)**2.*sqrt(cE0*cE3))
      cTmp=0.5/cD*cTmp
      ReTaucL=REAL(cTmp)
      RETURN
      END
      
      FUNCTION aImTaucL(E,E0,D,E3)
      COMPLEX aim,cE,cE0,cD,cE3,cTmp
      aim=CMPLX(0.,1.)
      cE =CMPLX(E ,0.)
      cE0=CMPLX(E0,0.)
      cD =CMPLX(D ,0.)
      cE3=CMPLX(E3,0.)
      cTmp=-(cD-aim*cE)**2.*sqrt(cD+aim*(cE-cE0))*sqrt(cD+aim*(cE-cE3))/
     + (cD*cD+cE*cE)**2.
      cTmp=cTmp-
     + (cD+aim*cE)**2.*sqrt(cD-aim*(cE-cE3))*sqrt(cD-aim*(cE-cE0))/
     + (cD*cD+cE*cE)**2.
      cTmp=cTmp+
     + cD*(cD*cD*(cE0+cE3)+cE*(-4.*cE0*cE3+cE*(cE0+cE3)))/
     + (cD*cD+cE*cE)**2./sqrt(cE0*cE3)
      aImTaucL=REAL(0.5*cTmp)
      RETURN
      END
            

      FUNCTION SIND1(a)
      real a
      SIND1=SIN(a/180.*acos(-1.))
      END


      FUNCTION COSD1(a)
      real a
      COSD1=COS(a/180.*acos(-1.))
      END




      FUNCTION ATAND1(a)
      real a
      ATAND1=180.*atan(a)/acos(-1.)
      END
      
      
      FUNCTION ATAN2D1(a,b)
      real a,b
      ATAN2D1=180.*atan2(a,b)/acos(-1.)
      END
