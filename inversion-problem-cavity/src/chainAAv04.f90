program chain
    implicit none
    character*80 outname,outname1,outname2,outname3,outname4,outname5
    character*80 outname6,outname7
    complex*16, allocatable :: H(:,:), SigmaS(:,:),SigmaD(:,:)
    complex*16, allocatable :: GgammaS(:,:), GgammaD(:,:), GF(:,:)

    integer Lx,Nph,NEpoints,Ndisorder,seed
    real*8 t, gam, omega, beta, gAA
    common /block1/ Lx,Nph,NEpoints,Ndisorder,t,gam,omega,beta,gAA,seed

    real*8 tcS,tcD,tlS,tlD,muS,muD
    common /block2/ tcS,tcD,tlS,tlD,muS,muD

    integer ii, NN, nei
    real*8, allocatable :: bTt(:,:),logbTt(:,:),lcbTt(:,:),LDOS(:,:,:),G2(:,:)
    real*8, allocatable :: vIPR(:)
    real*8 Energ, avgTt, errTt, varTt, desvTt

    real*8 IPR

    read(5,*) outname
    read(5,*) Lx, Nph, NEpoints, Ndisorder
    read(5,*) seed
    read(5,*) gAA, beta
    read(5,*) t, gam, omega
    read(5,*) tcS,tcD,tlS,tlD,muS,muD

    outname1='dados'//outname
    outname2='TT'//outname
    outname3='logTT'//outname
! outname4='avgGIPR'//outname
! outname5='bareGIPR'//outname
! outname6='LDOSNph0'//outname
! outname7='LDOSNph1'//outname
    open(unit = 1011, file =adjustl(adjustr(outname1)//'.dat'))
    open(unit = 1012, file =adjustl(adjustr(outname2)//'.dat'))
    open(unit = 1013, file =adjustl(adjustr(outname3)//'.dat'))
! open(unit = 1014, file =adjustl(adjustr(outname4)//'.dat'))
! open(unit = 1015, file =adjustl(adjustr(outname5)//'.dat'))
! open(unit = 1016, file =adjustl(adjustr(outname6)//'.dat'))
! open(unit = 1017, file =adjustl(adjustr(outname7)//'.dat'))

    call writeInput()

    NN=Lx*(Nph+1)
    allocate(H(NN,NN))
    allocate(SigmaS(NN,NN))
    allocate(SigmaD(NN,NN))
    allocate(GgammaS(NN,NN))
    allocate(GgammaD(NN,NN))
    allocate(GF(NN,NN))
    allocate(bTt(2*NEpoints-1,Ndisorder))
    allocate(logbTt(2*NEpoints-1,Ndisorder))
    allocate(LDOS(NN,2*NEpoints-1,Ndisorder))
    allocate(G2(2*NEpoints-1,Ndisorder))
    allocate(vIPR(Ndisorder))


    bTt(:,:)=0.d0
    logbTt(:,:)=0.d0
    LDOS(:,:,:)=0.d0
    do ii=1, Ndisorder
        write(*,*) "Disorder config.:", ii
        call Hamiltonian(H,NN,IPR)
        vIPR(ii) = IPR
        call Transmissivity(H,GgammaS,GgammaD,SigmaS,SigmaD,GF,NN,bTt,logbTt,LDOS,ii)
    end do

! call GIPR(NN,LDOS,G2)

! do ii=1,Lx
! write(1016,*) LDOS(ii,:, 1)
! write(1017,*) LDOS(Lx + ii,:, 1)
! end do

    do ii=1,2*NEpoints-1
        Energ= -2.d0 + 2.d0*dfloat(ii)/dfloat(NEpoints)
        call values1(bTt(ii,:),avgTt,varTt,desvTt,errTt,Ndisorder)
        !	write(1012,*) Energ, avgTt,varTt,desvTt,errTt
        ! write(1012,*) Energ, avgTt,,errTt   Verificar se é isso mesmo
        write(1012,*) Energ, avgTt, errTt
        call values1(logbTt(ii,:),avgTt,varTt,desvTt,errTt,Ndisorder)
!	 write(1013,*) Energ, avgTt,varTt,desvTt,errTt
        write(1013,*) Energ, avgTt,errTt
!	 call values1(G2(ii,:),avgTt,varTt,desvTt,errTt,Ndisorder)
!	 write(1014,*) Energ, avgTt, errTt
    end do

!	 call values1(vIPR(:),avgTt,varTt,desvTt,errTt,Ndisorder)
!	 write(1011,*) "IPR(0) =", avgTt, errTt

! do ii=1,2*NEpoints-1
! Energ= -2.d0 + 2.d0*dfloat(ii)/dfloat(NEpoints)
!   do nei=1,Ndisorder
!	 write(1015,*) Energ, G2(ii,nei)
!   end do
! end do

end program

!-------------------!
subroutine writeInput()
!-------------------!
    implicit none

    integer Lx,Nph,NEpoints,Ndisorder,seed
    real*8 t, gam, omega, beta, gAA
    common /block1/ Lx,Nph,NEpoints,Ndisorder,t,gam,omega,beta,gAA,seed

    real*8 tcS,tcD,tlS,tlD,muS,muD
    common /block2/ tcS,tcD,tlS,tlD,muS,muD

    write(1011,*) "Input data"
    write(1011,*) "Lx=",Lx, "Nph=", Nph
    write(1011,*) "gAA=",gAA, "beta=", beta, "seed =", seed
    write(1011,*) "Energy grid=", NEpoints, "Number of disorder conf.=", Ndisorder
    write(1011,*) "t=",t, "gam=", gam, "Omega=", omega
    write(1011,*) "tcS=",tcS,"tcD=",tcD,"tlS=",tlS,"tlD=",tlD
    write(1011,*) "muD=",muS,"muD=",muD

    return
end subroutine writeInput
!---------------------!

!-------------------!
subroutine values1(vm,aval,varval,desval,erval,bins)
!-------------------!
    implicit none
    integer bins
    real*8 aval,varval,desval,erval
    real*8 vm(0:bins-1), vr(bins)

    aval= sum(vm(:))/dfloat(bins)

    vr(:) = (vm(:) - aval)**2

    varval = sum(vr)/dfloat(bins - 1)
    desval = dsqrt( sum(vr)/dfloat(bins - 1) )
    erval = desval/dsqrt(dfloat(bins))

    return
end subroutine values1
!---------------------!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!-------------------!
subroutine GIPR(NN,LDOS,G2)
!-------------------!
    implicit none

    integer Lx,Nph,NEpoints,Ndisorder,seed
    real*8 t, gam, omega, beta, gAA
    common /block1/ Lx,Nph,NEpoints,Ndisorder,t,gam,omega,beta,gAA,seed

    real*8 LDOS(NN,2*NEpoints-1,Ndisorder)
    real*8 G2(2*NEpoints-1,Ndisorder), A0, A1, B
    integer i,j,ne,ii,NN


    G2(:,:)=0.d0
    do ii=1,Ndisorder
        do ne=1,2*NEPoints-1
            A1= 0.d0
            B = 0.d0

            do i= 1, Lx
                A0 = 0.d0
                do j=0,Nph
!	write(*,*) ii,ne,i,j,LDOS(j*Lx + i, ne, ii)
                    A0 = A0 + LDOS(j*Lx + i, ne, ii)
                    B  = B  + LDOS(j*Lx + i, ne, ii)
                end do!j
                A1 = A1 + A0**2
            end do! i
            G2(ne,ii)=A1/(B**2)
!	write(*,*) A1, (B**2)
        end do !ne
    end do!ii

    return
end subroutine GIPR
!---------------------!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Hamiltonian(H,NN,IPR)
    implicit none
! lista de variáveis usadas pelo Lapack para diagonalização de matrizes. Ver comentários no final do prog.
    character JOBZ, UPLO
    integer NN, LDA,LWORK, info
    complex*16 H(NN,NN), WORK(2*NN-1)
    real*8 lambda(NN),RWORK(3*NN-2)

    integer Lx,Nph,NEpoints,Ndisorder,seed
    real*8 t, gam, omega, beta, gAA
    common /block1/ Lx,Nph,NEpoints,Ndisorder,t,gam,omega,beta,gAA,seed

    real(8), allocatable :: P(:,:), delta(:,:)

    integer i,j,j1,j2,fats,fatj,s,N,M
    real*8 U(Lx),g
    complex*16 CI, hNM
    real*8 IPR

    complex*16 Hd(NN,NN)

    g=gam/t
    CI=(0.d0,1.d0)

    LDA=NN
    LWORK=2*NN - 1
! INFO=0
    JOBZ='V'
    UPLO='U'

    H(:,:) = 0.d0

    allocate(P(0:Nph,0:Nph))
    allocate(delta(0:Nph,0:Nph))

    call randomdistribution(U,Lx,gAA,beta,seed)
    call funcoesP(Nph,P,delta)

!%%%%%%%%%%%%% diagonal %%%%%%%%%%%%
! write(*,*) "escrevendo a matrix"

    do j=0,Nph
        do i= 1, Lx
            H(j*Lx + i, j*Lx +i) = dfloat(j)*omega + U(i)
        end do
    end do

!%%%%%%%%%%%%% off-diagonal %%%%%%%%%%%%
! write(*,*) "escrevendo a matrix"

    do j1=0,Nph !j1=N
        do j2=j1,Nph !j2=M

            hNM=0.d0

            fats=1
            do s=0,j1

                fatj=1
                do j=0,j2
                    hNM = hNM + dexp(-0.5d0*g**2)*((CI*g)**s)*((CI*g)**j)* &
                    & delta(j1-s,j2-j)*P(s,j1)*P(j,j2)/dfloat(fats*fatj)
                    fatj=fatj*(j+1)
                end do

                fats=fats*(s+1)
            end do

            do i=1,Lx-1
                H(j2*Lx + i, j1*Lx + i+1) = hNM*t
                H(j2*Lx + i+1, j1*Lx + i) = conjg(hNM)*t

                H(j1*Lx + i+1, j2*Lx + i) = conjg(hNM)*t
                H(j1*Lx + i, j2*Lx + i+1) = hNM*t
            end do

        end do
    end do

! write(*,*) H(1+4,1+4), H(1+4,2+4), H(1+4,3+4), H(1+4,4+4)
! write(*,*) H(2+4,1+4), H(2+4,2+4), H(2+4,3+4), H(2+4,4+4)
! write(*,*) H(3+4,1+4), H(3+4,2+4), H(3+4,3+4), H(3+4,4+4)
! write(*,*) H(4+4,1+4), H(4+4,2+4), H(4+4,3+4), H(4+4,4+4)


!%%%%%%%%%%%%% off-diagonal %%%%%%%%%%%%
! write(*,*) "Start: Diagonalizando a matrix"

    Hd(:,:)=H(:,:)

    call ZHEEV(JOBZ,UPLO,NN,Hd,LDA,lambda,WORK,LWORK,RWORK,INFO)



! write(*,*) "Finished: Diagonalizando a matrix"

    return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine funcoesP(Nph,P,delta)
    implicit none
    integer i,M,j,Nph
    real(8) P(0:Nph,0:Nph), delta(0:Nph,0:Nph)

    P(:,:)=0.d0
    P(0,:)=1.d0

    delta(:,:)=0.d0
    delta(0,0)=1.d0

    do M=1,Nph
        delta(M,M)=1.d0
        do j=1,M
            P(j,M) = dsqrt(dfloat(M-(j-1)))*P(j-1,M)
        end do
    end do



    return
end
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Transmissivity(H,GgammaS,GgammaD,SigmaS,SigmaD,GF,NN,bTt,logbTt,LDOS,ii)
    implicit none
! lista de variáveis usadas pelo Lapack para diagonalização de matrizes. Ver comentários no final do prog.
    integer NN, LDA, IPIV(NN),LWORK, info
    complex*16 WORK(NN)

    integer Lx,Nph,NEpoints,Ndisorder,seed
    real*8 t, gam, omega, beta, gAA
    common /block1/ Lx,Nph,NEpoints,Ndisorder,t,gam,omega,beta,gAA,seed

    real*8 tcS,tcD,tlS,tlD,muS,muD
    common /block2/ tcS,tcD,tlS,tlD,muS,muD

    complex*16 GgammaS(NN,NN), GgammaD(NN,NN), SigmaS(NN,NN), SigmaD(NN,NN)
    complex*16 GF(NN,NN), H(NN,NN)
    complex*16 Id(NN,NN)
    real*8 E, xD,xS,Tt, bTt(2*NEpoints-1,Ndisorder), logbTt(2*NEpoints-1,Ndisorder)
    real*8 LDOS(NN,2*NEpoints-1,Ndisorder)
    real*8 PI
    integer ne, i, j,ii,jj
    complex*16 CI

    LDA=NN
    LWORK=NN
    INFO=0

    PI=dacos(-1.d0)
    CI=(0.d0,1.d0)

    Id(:,:)=0.d0
    do j=0,Nph
        do i= 1, Lx
            Id(j*Lx + i, j*Lx +i) = 1.d0
        end do
    end do


    do ne=1,2*NEPoints-1
        E= -2.d0 + 2.d0*dfloat(ne)/dfloat(NEpoints)
        xS=(E-muS)/tlS
        xD=(E-muD)/tlD

        SigmaS(:,:)=0.d0
        SigmaD(:,:)=0.d0
        j=0
        i= 1
        SigmaS(j*Lx + i, j*Lx +i)=((tcS**2)/tlS)*0.5d0*(xS - CI*dsqrt(4.d0 - xS**2))

        i= Lx
        SigmaD(j*Lx + i, j*Lx +i)=((tcD**2)/tlD)*0.5d0*(xD - CI*dsqrt(4.d0 - xD**2))

        GF(:,:) =0.d0
        GF(:,:) = E*Id(:,:) - H(:,:) - SigmaS(:,:) - SigmaD(:,:)

        call ZGETRF(NN,NN,GF,NN,IPIV,INFO)
        if (INFO .ne. 0) then
            WRITE(*,*) INFO,"Wrong inversion"
            stop
        end if

        call zgetri(NN,GF,LDA,IPIV,WORK,LWORK,INFO)
        if (INFO .ne. 0) then
            WRITE(*,*) INFO,"Wrong inversion"
            stop
        end if

        GgammaS(:,:) = 0.d0
        GgammaS(:,:) = CI*( SigmaS(:,:) - CONJG(SigmaS(:,:)) )
        GgammaD(:,:) = 0.d0
        GgammaD(:,:) = CI*( SigmaD(:,:) - CONJG(SigmaD(:,:)) )


        do jj=1,NN
            LDOS(jj,ne,ii)=-aimag(GF(jj,jj))/PI
        end do


        call transm(Tt,GgammaS,GgammaD,GF,NN,Lx,Nph)

!  write(*,*) E,Tt, ((2.d0*dsin(dacos(E/2.d0)))**2)*(cdabs(GF(1,80)))**2
        BTt(ne,ii)=Tt
        logBTt(ne,ii)=dlog(Tt)
    end do

    return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine transm(Tt,GgammaS,GgammaD,GF,NN,Lx,Nph)
    implicit none
    integer NN,i,j,Lx,Nph
    complex*16 GgammaS(NN,NN), GgammaD(NN,NN),GF(NN,NN)
    complex*16 A(NN,NN),B(NN,NN),C(NN,NN),cGF(NN,NN)
    real*8 Tt

    A(:,:)=0.d0
    B(:,:)=0.d0
    C(:,:)=0.d0

    call zgemm('n','n',NN,NN,NN,1.0d0,GgammaS,NN,GF,NN,0.0d0,A,NN)

! cGF=CONJG(GF)
! cGF=TRANSPOSE(cGF)

    do i=1,NN
        do j=1,NN
            cGF(i,j) = CONJG(GF(j,i))
        end do
    end do

    call zgemm('n','n',NN,NN,NN,1.0d0,GgammaD,NN,cGF,NN,0.0d0,B,NN)
    call zgemm('n','n',NN,NN,NN,1.0d0,A,NN,B,NN,0.0d0,C,NN)


    Tt=0.d0
    do j=0,Nph
        do i= 1, Lx
            Tt = Tt + real(C(j*Lx + i, j*Lx +i))
        end do
    end do

    return
end
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine randomdistribution(U,Lx,gAA,beta,seed)
    implicit none
    real*8 U(Lx),PI,gAA,beta,phi,ran2
    integer Lx, i,seed

    PI=dacos(-1.d0)

    phi=2.d0*PI*ran2(seed)

!!  write(*,*) ">>>> Phi =", phi

    do i=1,Lx
        U(i) = gAA*dcos(2.d0*PI*beta*dfloat(i) + phi)
    end do

    return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DOUBLE PRECISION FUNCTION RAN2(IDUM)
    implicit none
    double precision rm
    integer ia,ic,idum,iff,ir,iy,j,m
    save
    PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112d-6)
    DIMENSION IR(97)
    DATA IFF /0/
    IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IC-IDUM,M)
        DO 11 J=1,97
            IDUM=MOD(IA*IDUM+IC,M)
            IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
    ENDIF
    J=1+(97*IY)/M
!c      IF(J.GT.97.OR.J.LT.1)PAUSE
    IY=IR(J)
    RAN2=IY*RM
    IDUM=MOD(IA*IDUM+IC,M)
    IR(J)=IDUM
    RETURN
END




