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
    integer                 :: ios, i, j, k
    integer,PARAMETER       :: nVars = 8 ! Number of cVariables for output
                                ! node 1-8, U, V, W, Rou, T, P, Ma, Cross

    ALLOCATE(Nodes(nCells*8,3))  ! Reserve enough space for Nodes-array.
    ALLOCATE(cVariables(nCells,nVars))
    ALLOCATE(cNodes(nCells,8))
    print*, 'Outputting ASCII data......'
    FileName=trim(OutputNameStr)//'-'//TimeStepStr//OutputFormat
    print*, 'Save to file: ',FileName
    call InitNodeInfo
    print"(A7,I10,/,A7,I10)", "Nodes:", nNodes, "Cells:", nCells
    call initTmpStorageVar  ! Temporary Storage cVariables

    open(21, file=FileName, iostat=ios, status="replace", action="write")
        if ( ios /= 0 ) stop ("Error====> Error opening file "//FileName)
        write(21,*) 'TITLE="TDCG-program Results"'
        !write(21,*) 'cVariables="X","Y","Z","U","V"'
        write(21,*) 'Variables="X","Y","Z","U","V","W","Rou"',      &
                    ',"T","P","Ma","Cross"'
        write(21,*) 'ZONE N=',nNodes,'E=',nCells
        write(21,*) 'DATAPACKING=BLOCK  ','ZONETYPE=FEbrick'
        write(21,"(1X,A28,I2,A15)")                    &
            'VARLOCATION=([1-3]=NODAL,[4-',nVars+3,']=CELLCENTERED)'

        write(21,"(F20.10)") ((real(Nodes(i,j),R8),i=1,nNodes),j=1,3)
        write(21,"(F20.10)") ((cVariables(i,j),i=1,nCells),j=1,nVars)
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
    ALLOCATE(cVariables(nCells,nVars))
    ALLOCATE(cNodes(nCells,8))
    print*, 'Outputting binary data......'
    FileName=trim(OutputNameStr)//'-'//TimeStepStr//OutputFormat
    print*, 'Save to file: ',FileName
    call InitNodeInfo
    print"(A7,I10,/,A7,I10)", "Nodes:", nNodes, "Cells:", nCells
    call initTmpStorageVar  ! Temporary Storage cVariables

    ! Open the file and write the tecplot datafile header information.
    ios = TecIni142('TDCG-program Results'//NULLCHR, &
                    'X Y Z U V W Rou T P Ma Cross'//NULLCHR, &
                    FileName//NULLCHR, &
                    NULLCHR, &
                    FileFormat, &        ! FileFormat
                    0, &        ! FileType
                    Debug, &
                    1)          ! VIsDouble     = 0 Single
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
                    VarLocation, &
                    Null, &     ! ShareVarFromZone
                    0)          ! ShareConnectivityFromZone
        if (ios/=0) stop "Error value returned in TecZne142"
    do i=1,3
    ios = TecDat142(nNodes,real(Nodes(1:nNodes,i),R4),0)
        if (ios/=0) stop "Error value returned in TecDat142"
    enddo
    do i=1,nVars
    ios = TecDat142(nCells,real(cVariables(1:nCells,i),R4),0)
        if (ios/=0) stop "Error value returned in TecDat142"
    enddo
    ios = TecNod142(transpose(cNodes))
        if (ios/=0) stop "Error value returned in TecNod142"
    ios = TecEnd142()
        if (ios/=0) stop "Error value returned in TecEnd142"
    DEALLOCATE(Nodes)
    DEALLOCATE(cVariables)
    DEALLOCATE(cNodes)
    print*, 'Done'
    endsubroutine OutputFlowFieldBinary
!======================================================================
    subroutine InitNodeInfo    ! nCellst=nCells, as a parameter form.
    use ModMesh
    use ModOutput
    use ModTypDef
    implicit none
    type(typCell),pointer :: ct
    integer :: i
    
    !nNodes=1
    !Nodes(1,1:3)=0.0

    do i=1,nBGCells
        ct=>Cell(i)
        call NodeInfo(ct)
    enddo

    contains
!----------------------------------------------------------------------
    recursive subroutine NodeInfo(c)
    use ModTypDef
    use ModMesh
    use ModOutput
    implicit none
    type(typCell),pointer :: c
    real(R8):: step(3), tN(6) ! Temp-Nodes
    logical :: Mark(8)
    integer :: ii

    if(ASSOCIATED(c%son8))then
        call NodeInfo(c%son1)
        call NodeInfo(c%son2)
        call NodeInfo(c%son3)
        call NodeInfo(c%son4)
        call NodeInfo(c%son5)
        call NodeInfo(c%son6)
        call NodeInfo(c%son7)
        call NodeInfo(c%son8)
    elseif(ASSOCIATED(c%son4))then
        call NodeInfo(c%son1)
        call NodeInfo(c%son2)
        call NodeInfo(c%son3)
        call NodeInfo(c%son4)
    elseif(ASSOCIATED(c%son2))then
        call NodeInfo(c%son1)
        call NodeInfo(c%son2)
    else
        Mark=.true.
        do ii=1,3
            step(ii)=BGStep(ii)/(2**(c%lvl(ii)+1))
        enddo
        ! Initial node number.
        tN(1)=c%Center(1)-step(1)   ! x -
        tN(2)=c%Center(2)-step(2)   ! y -
        tN(3)=c%Center(3)-step(3)   ! z -
        tN(4)=c%Center(1)+step(1)   ! x +
        tN(5)=c%Center(2)+step(2)   ! y +
        tN(6)=c%Center(3)+step(3)   ! z +
        ! Nodes array in xyz:
        ! 1 --- 2 +-- 3 ++- 4 -+- 5 --+ 6 +-+ 7 +++ 8 -++
        ! do ii=1,nNodes
        !     if (Nodes(ii,1)==tN(1)) then     ! x -
        !         if (Nodes(ii,2)==tN(2)) then     ! y -
        !             if (Nodes(ii,3)==tN(3)) then     ! z -
        !                 c%Node(1)=ii; mark(1)=.false.          ! 1 ---
        !             elseif (Nodes(ii,3)==tN(6)) then ! z +
        !                 c%Node(5)=ii; mark(5)=.false.          ! 5 --+
        !             endif
        !         elseif (Nodes(ii,2)==tN(5)) then ! y +
        !             if (Nodes(ii,3)==tN(3)) then     ! z -
        !                 c%Node(4)=ii; mark(4)=.false.          ! 4 -+-
        !             elseif (Nodes(ii,3)==tN(6)) then ! z +
        !                 c%Node(8)=ii; mark(8)=.false.          ! 8 -++
        !             endif
        !         endif
        !     elseif (Nodes(ii,1)==tN(4)) then ! x +
        !         if (Nodes(ii,2)==tN(2)) then     ! y -
        !             if (Nodes(ii,3)==tN(3)) then     ! z -
        !                 c%Node(2)=ii; mark(2)=.false.          ! 2 +--
        !             elseif (Nodes(ii,3)==tN(6)) then ! z +
        !                 c%Node(6)=ii; mark(6)=.false.          ! 6 +-+
        !             endif
        !         elseif (Nodes(ii,2)==tN(5)) then ! y +
        !             if (Nodes(ii,3)==tN(3)) then     ! z -
        !                 c%Node(3)=ii; mark(3)=.false.          ! 3 ++-
        !             elseif (Nodes(ii,3)==tN(6)) then ! z +
        !                 c%Node(7)=ii; mark(7)=.false.          ! 7 +++
        !             endif
        !         endif
        !     endif
        ! enddo

        if (mark(1)) then   ! Node 1 ---
            nNodes=nNodes+1
            c%Node(1)=nNodes
            Nodes(nNodes,1)=tN(1)
            Nodes(nNodes,2)=tN(2)
            Nodes(nNodes,3)=tN(3)
        endif
        if (mark(2)) then   ! Node 2 +--
            nNodes=nNodes+1
            c%Node(2)=nNodes
            Nodes(nNodes,1)=tN(4)
            Nodes(nNodes,2)=tN(2)
            Nodes(nNodes,3)=tN(3)
        endif
        if (mark(3)) then   ! Node 3 ++-
            nNodes=nNodes+1
            c%Node(3)=nNodes
            Nodes(nNodes,1)=tN(4)
            Nodes(nNodes,2)=tN(5)
            Nodes(nNodes,3)=tN(3)
        endif
        if (mark(4)) then   ! Node 4 -+-
            nNodes=nNodes+1
            c%Node(4)=nNodes
            Nodes(nNodes,1)=tN(1)
            Nodes(nNodes,2)=tN(5)
            Nodes(nNodes,3)=tN(3)
        endif
        if (mark(5)) then   ! Node 5 --+
            nNodes=nNodes+1
            c%Node(5)=nNodes
            Nodes(nNodes,1)=tN(1)
            Nodes(nNodes,2)=tN(2)
            Nodes(nNodes,3)=tN(6)
        endif
        if (mark(6)) then   ! Node 6 +-+
            nNodes=nNodes+1
            c%Node(6)=nNodes
            Nodes(nNodes,1)=tN(4)
            Nodes(nNodes,2)=tN(2)
            Nodes(nNodes,3)=tN(6)
        endif
        if (mark(7)) then   ! Node 7 +++
            nNodes=nNodes+1
            c%Node(7)=nNodes
            Nodes(nNodes,1)=tN(4)
            Nodes(nNodes,2)=tN(5)
            Nodes(nNodes,3)=tN(6)
        endif
        if (mark(8)) then   ! Node 8 -++
            nNodes=nNodes+1
            c%Node(8)=nNodes
            Nodes(nNodes,1)=tN(1)
            Nodes(nNodes,2)=tN(5)
            Nodes(nNodes,3)=tN(6)
        endif
    endif
    endsubroutine NodeInfo
!----------------------------------------------------------------------
    endsubroutine InitNodeInfo
!======================================================================
!======================================================================
    subroutine initTmpStorageVar
    use ModMesh
    use ModOutput
    implicit none
    type(typCell),pointer :: ct
    integer :: i
    do i = 1, nBGCells
        ct=>Cell(i)
        call TmpStorageVar(ct)
    enddo
        contains
!----------------------------------------------------------------------
        recursive subroutine TmpStorageVar(c)
        use ModInpInflow,only : Rgas, Gama00
        implicit none
        type(typCell),pointer :: c
        integer :: n, j
        REAL(R8):: u, v, w,p

        if(ASSOCIATED(c%son8))then
            call TmpStorageVar(c%son1)
            call TmpStorageVar(c%son2)
            call TmpStorageVar(c%son3)
            call TmpStorageVar(c%son4)
            call TmpStorageVar(c%son5)
            call TmpStorageVar(c%son6)
            call TmpStorageVar(c%son7)
            call TmpStorageVar(c%son8)
        elseif(ASSOCIATED(c%son4))then
            call TmpStorageVar(c%son1)
            call TmpStorageVar(c%son2)
            call TmpStorageVar(c%son3)
            call TmpStorageVar(c%son4)
        elseif(ASSOCIATED(c%son2))then
            call TmpStorageVar(c%son1)
            call TmpStorageVar(c%son2)
        else
            n=c%nCell
            do j=1,8; cNodes(n,j)=c%Node(j); enddo
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
            if (c%U(4)==0 .or. c%U(5)==0) cVariables(n,:)=-9999999
        endif
        endsubroutine TmpStorageVar
!----------------------------------------------------------------------
    endsubroutine initTmpStorageVar
!======================================================================
!======================================================================
!======================================================================
!----------------------------------------------------------------------
!----------------------------------------------------------------------
