!======================================================================
    module ModOutput
        use ModTypDef
        implicit none
        real(R8),ALLOCATABLE :: Nodes(:,:)
        integer :: nNodes

    endmodule ModOutput
!======================================================================
    subroutine OutputFlowField(TimeStepStr)
    use ModInpGlobal
    use ModOutput
    use ModInpMesh
    use ModMesh,    only : nCells, Cell
    implicit none
    character(*),INTENT(IN):: TimeStepStr
    character(80):: FileName
    integer :: ios, tp  ! Temp Precision
    integer :: i, j
    !real(R8):: minSize

    ALLOCATE(Nodes(nCells*8,3))  ! Reserve enough space for Nodes-array.
    print*, 'Outputting data......'
    FileName=trim(NameStr)//'-'//TimeStepStr//'.dat'
    print*, 'Saving to file: ',FileName
    call InitNodeInfo
    print"(A7,I10,/,A7,I10)", "Nodes:", nNodes, "Cells:", nCells

    open(21, file=FileName, iostat=ios, status="replace", action="write")
        if ( ios /= 0 ) stop ("Error====> Error opening file "//FileName)
        write(21,*) 'TITLE="3D Results"'
        write(21,*) 'VARIABLES="X","Y","Z"'
        write(21,*) 'ZONE N=',nNodes,'E=',nCells,'ZONETYPE=FEbrick'
        write(21,*) 'DATAPACKING=BLOCK'
        write(21,*) 'VARLOCATION=([1-3]=NODAL)'
        !write(21,*) 'VARIABLES="X","Y","Z","U","V","W","Rou"',      &
        !            ',"T","P","Ma","Cross"'
        ! write(21,*) 'ZONE N=',nNodes,'E=',nCells,'ZONETYPE=FEbrick'
        ! write(21,*) 'DATAPACKING=BLOCK'
        ! write(21,*) 'VARLOCATION=([1-3]=NODAL,[4-11]=CELLCENTERED)'

        ! Single precision -- croase mesh
        ! Double precision -- fine mesh
        if (InitRefineLVL <= 5) then
            write(21,"(F12.6)") ((real(Nodes(i,j),R4),i=1,nNodes),j=1,3)
        else
            write(21,"(F16.10)") ((real(Nodes(i,j),R8),i=1,nNodes),j=1,3)
        endif
        call initTSData  ! TS -- Temporary Storage
        write(21,"(8I10)") ((Cell(j)%Node(i),i=1,8),j=1,nCells)
    close(21)
    DEALLOCATE(Nodes)
    print*, 'Done'
    end subroutine OutputFlowField
!======================================================================
    subroutine InitNodeInfo    ! nCellst=nCells, as a parameter form.
    use ModMesh
    use ModOutput
    use ModTypDef
    implicit none
    type(typCell),pointer :: ct
    integer :: i
    
    nNodes=1
    Nodes(1,1:3)=0.0

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
        do ii=1,3; step(ii)=BGStep(ii)/(2**c%lvl(ii)+1); enddo
        ! Initial node number.
        tN(1)=c%Center(1)-step(1)
        tN(2)=c%Center(2)-step(2)
        tN(3)=c%Center(3)-step(3)
        tN(4)=c%Center(1)+step(1)
        tN(5)=c%Center(2)+step(2)
        tN(6)=c%Center(3)+step(3)
! Nodes array in xyz:  1 --- 2 +-- 3 ++- 4 -+- 5 --+ 6 +-+ 7 +++ 8 -++
        do ii=1,nNodes
            if (Nodes(ii,1)==tN(1)) then     ! x -
                if (Nodes(ii,2)==tN(2)) then     ! y -
                    if (Nodes(ii,3)==tN(3)) then     ! z -
                        c%Node(1)=ii; mark(1)=.false.          ! 1 ---
                    elseif (Nodes(ii,3)==tN(6)) then ! z +
                        c%Node(5)=ii; mark(5)=.false.          ! 5 --+
                    endif
                elseif (Nodes(ii,2)==tN(5)) then ! y +
                    if (Nodes(ii,3)==tN(3)) then     ! z -
                        c%Node(4)=ii; mark(4)=.false.          ! 4 -+-
                    elseif (Nodes(ii,3)==tN(6)) then ! z +
                        c%Node(8)=ii; mark(8)=.false.          ! 8 -++
                    endif
                endif
            elseif (Nodes(ii,1)==tN(4)) then ! x +
                if (Nodes(ii,2)==tN(2)) then     ! y -
                    if (Nodes(ii,3)==tN(3)) then     ! z -
                        c%Node(2)=ii; mark(2)=.false.          ! 2 +--
                    elseif (Nodes(ii,3)==tN(6)) then ! z +
                        c%Node(6)=ii; mark(6)=.false.          ! 6 +-+
                    endif
                elseif (Nodes(ii,2)==tN(5)) then ! y +
                    if (Nodes(ii,3)==tN(3)) then     ! z -
                        c%Node(3)=ii; mark(3)=.false.          ! 3 ++-
                    elseif (Nodes(ii,3)==tN(6)) then ! z +
                        c%Node(7)=ii; mark(7)=.false.          ! 7 +++
                    endif
                endif
            endif
        enddo

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
    subroutine initTSData
    use ModMesh
    implicit none
    type(typCell),pointer :: ct
    integer :: i
    do i = 1, nBGCells
        ct=>Cell(i)
        call TSData(ct)
    enddo
        contains
!----------------------------------------------------------------------
        subroutine TSData(c)
        implicit none
        type(typCell),pointer :: c

        !write(21,*), (c%nCell(i), i=1,8)
        endsubroutine TSData
!----------------------------------------------------------------------
    endsubroutine initTSData
!======================================================================
!======================================================================
!======================================================================
!----------------------------------------------------------------------
!----------------------------------------------------------------------
