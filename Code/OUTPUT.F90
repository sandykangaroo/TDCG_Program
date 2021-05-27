!======================================================================
    module ModOutput
        use ModTypDef
        implicit none
        real(R8),ALLOCATABLE :: Nodes(:,:)      ! Nodes coordinate
        real(R8),ALLOCATABLE :: cVariables(:,:) ! Cell variables value
        integer ,ALLOCATABLE :: cNodes(:,:)     ! Cell nodes number
        integer :: nNodes

    endmodule ModOutput
!======================================================================

!======================================================================
    subroutine OutputFlowFieldASCII(TimeStepStr)
    use ModInpGlobal
    use ModOutput
    use ModInpMesh
    use ModMesh,    only : nCells
    implicit none
    character(*),INTENT(IN) :: TimeStepStr
    character(80)           :: FileName
    integer                 :: iost, i, j
    integer,PARAMETER       :: nVars = 8 ! Number of cVariables for output
                                ! node 1-8, U, V, W, Rou, T, P, Ma, Cross
    integer                 :: ios ! If most of cells is NAN, output mesh only

    ALLOCATE(Nodes(nCells*8,3))  ! Reserve enough space for Nodes-array.
    ALLOCATE(cNodes(nCells,8))
    print*, 'Outputting ASCII data......'
    FileName=trim(OutputName)//'-'//TimeStepStr//'.'//OutputFormat
    print*, 'Save to file: ',FileName
    call InitNodeInfo
    print"(A7,I10,/,A7,I10)", "Nodes:", nNodes, "Cells:", nCells

    open(21, file=FileName, iostat=iost, status="replace", action="write")
        if ( iost /= 0 ) stop ("Error====> Error opening file ")
        write(21,*) 'TITLE="TDCG-program Results"'

        if (MeshOnly) then
            ALLOCATE(cVariables(nCells,1))
            call initTmpStorageVar(MeshOnly)  ! Temporary Storage cVariables
            write(21,*) 'Variables="X","Y","Z","Cross"'
            write(21,*) 'ZONE N=',nNodes,'E=',nCells
            write(21,*) 'DATAPACKING=BLOCK  ','ZONETYPE=FEbrick'
            write(21,"(1X,A)") 'VARLOCATION=([1-3]=NODAL,[4]=CELLCENTERED)'
            write(21,"(F20.10)") ((real(Nodes(i,j),R8),i=1,nNodes),j=1,3)
            write(21,"(F20.10)") (cVariables(i,1),i=1,nCells)
        else
            ALLOCATE(cVariables(nCells,nVars))
            call initTmpStorageVar(MeshOnly)  ! Temporary Storage cVariables
            write(21,*) 'Variables="X","Y","Z","U","V","W","Rou"',      &
                        ',"T","P","Ma","Cross"'
            write(21,*) 'ZONE N=',nNodes,'E=',nCells
            write(21,*) 'DATAPACKING=BLOCK  ','ZONETYPE=FEbrick'
            write(21,"(1X,A28,I2,A15)")                    &
                'VARLOCATION=([1-3]=NODAL,[4-',nVars+3,']=CELLCENTERED)'
            write(21,"(F20.10)") ((real(Nodes(i,j),R8),i=1,nNodes),j=1,3)
            write(21,"(F20.10)") ((cVariables(i,j),i=1,nCells),j=1,nVars)
        endif
        write(21,"(8(I9,1X))") ((cNodes(i,j),j=1,8),i=1,nCells)
    close(21)
    DEALLOCATE(Nodes)
    DEALLOCATE(cVariables)
    DEALLOCATE(cNodes)
    print*, 'Done'
    end subroutine OutputFlowFieldASCII
!======================================================================
    subroutine OutputFlowFieldBinary(TimeStepStr,FileFormat)
    use ModInpGlobal
    use ModOutput
    use ModInpMesh
    use ModMesh,    only : nCells
    use ModSolve,   only : step, TimeStep
    implicit none
    include '../lib/tecio.f90'
    character(*),INTENT(IN) :: TimeStepStr
    character(80)           :: FileName
    integer                 :: ios, i
    integer,PARAMETER       :: nVars = 8 ! Number of cVariables for output
    character*1             :: NULLCHR
    real(R8)                :: SolTime
    Integer*4               :: FileFormat
    Integer*4               :: VarLocation(11)
    INTEGER*4,TARGET        :: NULL(11)
    integer*4,pointer       :: NullPtr(:)
    
    NullPtr => Null
    NullPtr =  0
    VarLocation=(/1,1,1,0,0,0,0,0,0,0,0/)
    NULLCHR = CHAR(0)
    SolTime = real(step,R8)*TimeStep
    
    ALLOCATE(Nodes(nCells*8,3))  ! Reserve enough space for Nodes-array.
    ALLOCATE(cNodes(nCells,8))
    print*, 'Outputting binary data......'
    FileName=trim(OutputName)//'-'//TimeStepStr//'.'//OutputFormat
    print*, 'Save to file: ',FileName
    call InitNodeInfo
    print"(A7,I10,/,A7,I10)", "Nodes:", nNodes, "Cells:", nCells
    
    if (MeshOnly) then
        ! Open the file and write the tecplot datafile header information.
        ios = TecIni142('TDCG-program Results'//NULLCHR, &
                        'X Y Z Cross'//NULLCHR, &
                        FileName//NULLCHR, &
                        NULLCHR, &
                        FileFormat, &        ! FileFormat
                        0, &        ! FileType
                        Debug, &
                        0)          ! VIsDouble     = 0 Single
                                    !               = 1 Double
            if (ios/=0) stop "Error value returned in TecIni142"
        ! Write the zone header information.
        ios = TecZne142('Zone'//NULLCHR, &
                        5, &        ! ZoneType
                        nNodes, &   ! NumPts
                        nCells, &   ! NumElements
                        0, &        ! Not used for FEbrick Zone type.
                        0, 0, 0, &  ! For future use.
                        SolTime, &  ! SolutionTime
                        0, &        ! StrandID
                        0, &        ! ParentZn
                        1, &        ! IsBlock
                        0, &        ! NumFaceConnections
                        3, &        ! FaceNeighborMode
                        0, 0, 0, &  ! Not used for FEbrick Zone type.
                        Null, &     ! PassiveVarList
                        [1,1,1,0], & ! VarLocation
                        Null, &     ! ShareVarFromZone
                        0)          ! ShareConnectivityFromZone
            if (ios/=0) stop "Error value returned in TecZne142"
            
        do i=1,3
        ios = TecDat142(nNodes,real(Nodes(1:nNodes,i),R4),0)
            if (ios/=0) stop "Error value returned in TecDat142"
        enddo
        DEALLOCATE(Nodes)

        ALLOCATE(cVariables(nCells,1))
        call initTmpStorageVar(MeshOnly)  ! Temporary Storage cVariables
        ios = TecDat142(nCells,real(cVariables(1:nCells,1),R4),0)
        DEALLOCATE(cVariables)
            if (ios/=0) stop "Error value returned in TecDat142"
    else
        ! Open the file and write the tecplot datafile header information.
        ios = TecIni142('TDCG-program Results'//NULLCHR, &
                        'X Y Z U V W Rou T P Ma Cross'//NULLCHR, &
                        FileName//NULLCHR, &
                        NULLCHR, &
                        FileFormat, &        ! FileFormat
                        0, &        ! FileType
                        Debug, &
                        0)          ! VIsDouble     = 0 Single
                                    !               = 1 Double
            if (ios/=0) stop "Error value returned in TecIni142"
        ! Write the zone header information.
        ios = TecZne142('Zone'//NULLCHR, &
                        5, &        ! ZoneType
                        nNodes, &   ! NumPts
                        nCells, &   ! NumElements
                        0, &        ! Not used for FEbrick Zone type.
                        0, 0, 0, &  ! For future use.
                        SolTime, &  ! SolutionTime
                        0, &        ! StrandID
                        0, &        ! ParentZn
                        1, &        ! IsBlock
                        0, &        ! NumFaceConnections
                        3, &        ! FaceNeighborMode
                        0, 0, 0, &  ! Not used for FEbrick Zone type.
                        Null, &     ! PassiveVarList
                        [1,1,1,0,0,0,0,0,0,0,0], & ! VarLocation
                        Null, &     ! ShareVarFromZone
                        0)          ! ShareConnectivityFromZone
            if (ios/=0) stop "Error value returned in TecZne142"

        do i=1,3
        ios = TecDat142(nNodes,real(Nodes(1:nNodes,i),R4),0)
            if (ios/=0) stop "Error value returned in TecDat142"
        enddo
        DEALLOCATE(Nodes)

        ALLOCATE(cVariables(nCells,nVars))
        call initTmpStorageVar(MeshOnly)  ! Temporary Storage cVariables
        do i=1,nVars
        ios = TecDat142(nCells,real(cVariables(1:nCells,i),R4),0)
            if (ios/=0) stop "Error value returned in TecDat142"
        enddo
        DEALLOCATE(cVariables)

    endif
    ios = TecNod142(transpose(cNodes))
    DEALLOCATE(cNodes)
        if (ios/=0) stop "Error value returned in TecNod142"
    ios = TecEnd142()
        if (ios/=0) stop "Error value returned in TecEnd142"
    print*, 'Done'
    endsubroutine OutputFlowFieldBinary
!======================================================================
    subroutine InitNodeInfo    ! nCellst=nCells, as a parameter form.
    use ModPrecision
    use ModInpMesh
    use ModMesh
    use ModOutput
    use ModTypDef
    implicit none
    type(octCell),pointer :: t
    integer :: i, j, k
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        t=>Cell(i, j, k)
        call NodeInfo(t)
    enddo
    enddo
    enddo
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "InitNodeInfo time: ", tEnd-tStart

    contains
!----------------------------------------------------------------------
    recursive subroutine NodeInfo(c)
    use ModTypDef
    use ModMesh
    use ModOutput
    implicit none
    type(octCell),pointer :: c
    real(R8):: dx, dy, dz, tN(6) ! Temp-Nodes
    integer,SAVE:: n=0

    if(ASSOCIATED(c%son8))then
        call NodeInfo(c%son1)
        call NodeInfo(c%son2)
        call NodeInfo(c%son3)
        call NodeInfo(c%son4)
        call NodeInfo(c%son5)
        call NodeInfo(c%son6)
        call NodeInfo(c%son7)
        call NodeInfo(c%son8)
        return
    elseif(ASSOCIATED(c%son4))then
        call NodeInfo(c%son1)
        call NodeInfo(c%son2)
        call NodeInfo(c%son3)
        call NodeInfo(c%son4)
        return
    elseif(ASSOCIATED(c%son2))then
        call NodeInfo(c%son1)
        call NodeInfo(c%son2)
        return
    endif

    n=n+1
    dx=BGCellSize(1)/2**(c%lvl(1)+1)
    dy=BGCellSize(2)/2**(c%lvl(2)+1)
    dz=BGCellSize(3)/2**(c%lvl(3)+1)
    ! Initial node number.
    tN(1)=c%Center(1)-dx   ! x -
    tN(2)=c%Center(2)-dy   ! y -
    tN(3)=c%Center(3)-dz   ! z -
    tN(4)=c%Center(1)+dx   ! x +
    tN(5)=c%Center(2)+dy   ! y +
    tN(6)=c%Center(3)+dz   ! z +

    nNodes=nNodes+1
    ! c%Node(1)=nNodes
    Nodes(nNodes,1)=tN(1)
    Nodes(nNodes,2)=tN(2)
    Nodes(nNodes,3)=tN(3)
    cNodes(n,1)=nNodes

    nNodes=nNodes+1
    ! c%Node(2)=nNodes
    Nodes(nNodes,1)=tN(4)
    Nodes(nNodes,2)=tN(2)
    Nodes(nNodes,3)=tN(3)
    cNodes(n,2)=nNodes

    nNodes=nNodes+1
    ! c%Node(3)=nNodes
    Nodes(nNodes,1)=tN(4)
    Nodes(nNodes,2)=tN(5)
    Nodes(nNodes,3)=tN(3)
    cNodes(n,3)=nNodes

    nNodes=nNodes+1
    ! c%Node(4)=nNodes
    Nodes(nNodes,1)=tN(1)
    Nodes(nNodes,2)=tN(5)
    Nodes(nNodes,3)=tN(3)
    cNodes(n,4)=nNodes

    nNodes=nNodes+1
    ! c%Node(5)=nNodes
    Nodes(nNodes,1)=tN(1)
    Nodes(nNodes,2)=tN(2)
    Nodes(nNodes,3)=tN(6)
    cNodes(n,5)=nNodes

    nNodes=nNodes+1
    ! c%Node(6)=nNodes
    Nodes(nNodes,1)=tN(4)
    Nodes(nNodes,2)=tN(2)
    Nodes(nNodes,3)=tN(6)
    cNodes(n,6)=nNodes

    nNodes=nNodes+1
    ! c%Node(7)=nNodes
    Nodes(nNodes,1)=tN(4)
    Nodes(nNodes,2)=tN(5)
    Nodes(nNodes,3)=tN(6)
    cNodes(n,7)=nNodes

    nNodes=nNodes+1
    ! c%Node(8)=nNodes
    Nodes(nNodes,1)=tN(1)
    Nodes(nNodes,2)=tN(5)
    Nodes(nNodes,3)=tN(6)
    cNodes(n,8)=nNodes

    endsubroutine NodeInfo
!----------------------------------------------------------------------
    endsubroutine InitNodeInfo
!======================================================================
!======================================================================
    subroutine initTmpStorageVar(meshonly)
    ! ios = 0, mesh only; = 1, all vars
    use ModPrecision
    use ModInpMesh
    use ModMesh
    use ModOutput
    implicit none
    type(octCell),pointer :: t
    integer :: i, j, k
    logical,INTENT(IN) :: meshonly ! If most of cells is NAN, output mesh only
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        t=>Cell(i, j, k)
        call TmpStorageVar(t,meshonly)
    enddo
    enddo
    enddo
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "initTmpStorageVar time: ", tEnd-tStart
        contains
!----------------------------------------------------------------------
        recursive subroutine TmpStorageVar(c,ios)
        use ModInpInflow,only : Rgas, Gama00
        implicit none
        type(octCell),pointer :: c
        logical,INTENT(IN)    :: ios
        ! If most of cells is NAN, output mesh only
        integer,save :: n=0
        REAL(R8):: u, v, w,p

        if(ASSOCIATED(c%son8))then
            call TmpStorageVar(c%son1,ios)
            call TmpStorageVar(c%son2,ios)
            call TmpStorageVar(c%son3,ios)
            call TmpStorageVar(c%son4,ios)
            call TmpStorageVar(c%son5,ios)
            call TmpStorageVar(c%son6,ios)
            call TmpStorageVar(c%son7,ios)
            call TmpStorageVar(c%son8,ios)
            return
        elseif(ASSOCIATED(c%son4))then
            call TmpStorageVar(c%son1,ios)
            call TmpStorageVar(c%son2,ios)
            call TmpStorageVar(c%son3,ios)
            call TmpStorageVar(c%son4,ios)
            return
        elseif(ASSOCIATED(c%son2))then
            call TmpStorageVar(c%son1,ios)
            call TmpStorageVar(c%son2,ios)
            return
        endif

        ! n=c%nCell
        ! do jj=1,8
        !     cNodes(n,jj)=c%Node(jj)
        ! enddo
        n=n+1
        if (ios)then
            cVariables(n,1)= c%cross    ! Cross
        else
            u=c%U(1)/c%U(4); v=c%U(2)/c%U(4); w=c%U(3)/c%U(4)
            p=Rgas*c%U(4)*c%U(5)
            cVariables(n,1)= u          ! u
            cVariables(n,2)= v          ! v
            cVariables(n,3)= w          ! w
            cVariables(n,4)= c%U(4)     ! rou
            cVariables(n,5)= c%U(5)     ! T
            cVariables(n,6)= p          ! P
            cVariables(n,7)= sqrt((u*u+v*v+w*w)/abs(Gama00*p/c%U(4)))!Ma
            cVariables(n,8)= c%cross    ! Cross
        endif
        endsubroutine TmpStorageVar
!----------------------------------------------------------------------
    endsubroutine initTmpStorageVar
!======================================================================
!======================================================================
!======================================================================
!----------------------------------------------------------------------
!----------------------------------------------------------------------
