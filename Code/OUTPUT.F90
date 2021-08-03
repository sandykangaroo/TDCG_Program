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
    integer,PARAMETER       :: nVars = 9 ! Number of cVariables for output
                                ! U, V, W, Rou, T, P, Ma, Cross,
                                !walldistance,node 1-8, 
    integer                 :: ios ! If most of cells is NAN, output mesh only

    ALLOCATE(Node(nCells*8,3))  ! Reserve enough space for Nodes-array.
    ALLOCATE(cNodes(nCells,8))
    print*, 'Outputting ASCII data......'
    FileName=trim(OutputName)//'-'//TimeStepStr//'.'//OutputFormat
    print*, 'Save to file: ',FileName
    call InitNodeInfo
    print"(A7,I10,/,A7,I10)", "Nodes:", nNode, "Cells:", nCells

    open(21, file=FileName, iostat=iost, status="replace", action="write")
        if ( iost /= 0 ) stop ("Error====> Error opening file ")
        write(21,*) 'TITLE="TDCG-program Results"'

        if (MeshOnly) then
            ALLOCATE(cVariables(nCells,2))
            call initTmpStorageVar(MeshOnly)  ! Temporary Storage cVariables
            write(21,*) 'Variables="X","Y","Z","Cross","WallDistance"'
            write(21,*) 'ZONE N=',nNode,'E=',nCells
            write(21,*) 'DATAPACKING=BLOCK  ','ZONETYPE=FEbrick'
            write(21,"(1X,A)") 'VARLOCATION=([1-3]=NODAL,[4-5]=CELLCENTERED)'
            write(21,"(F20.10)") ((real(Node(i,j),R8),i=1,nNode),j=1,3)
            write(21,"(F20.10)") ((cVariables(i,j),i=1,nCells),j=1,2)
        else
            ALLOCATE(cVariables(nCells,nVars))
            call initTmpStorageVar(MeshOnly)  ! Temporary Storage cVariables
            write(21,*) 'Variables="X","Y","Z","U","V","W","Rou"',      &
                        ',"T","P","Ma","Cross","WallDistance"'
            write(21,*) 'ZONE N=',nNode,'E=',nCells
            write(21,*) 'DATAPACKING=BLOCK  ','ZONETYPE=FEbrick'
            write(21,"(1X,A28,I2,A15)")                    &
                'VARLOCATION=([1-3]=NODAL,[4-',nVars+3,']=CELLCENTERED)'
            write(21,"(F20.10)") ((real(Node(i,j),R8),i=1,nNode),j=1,3)
            write(21,"(F20.10)") ((cVariables(i,j),i=1,nCells),j=1,nVars)
        endif
        write(21,"(8(I9,1X))") ((cNodes(i,j),j=1,8),i=1,nCells)
    close(21)
    DEALLOCATE(Node)
    DEALLOCATE(cVariables)
    !DEALLOCATE(cNodes)
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
    integer,PARAMETER       :: nVars = 9 ! Number of cVariables for output
    character*1             :: NULLCHR
    real(R8)                :: SolTime
    Integer*4               :: FileFormat
    Integer*4               :: VarLocation(12)
    INTEGER*4,TARGET        :: NULL(11)
    integer*4,pointer       :: NullPtr(:)

    NullPtr => Null
    NullPtr =  0
    VarLocation=(/1,1,1,0,0,0,0,0,0,0,0,0/)
    NULLCHR = CHAR(0)
    SolTime = real(step,R8)*TimeStep

    ALLOCATE(Node(nCells*8,3))  ! Reserve enough space for Nodes-array.
    ALLOCATE(cNodes(nCells,8))
    print*, 'Outputting binary data......'
    FileName=trim(OutputName)//'-'//TimeStepStr//'.'//OutputFormat
    print*, 'Save to file: ',FileName
    call InitNodeInfo
    print"(A7,I10,/,A7,I10)", "Nodes:", nNode, "Cells:", nCells

    if (MeshOnly) then
        ! Open the file and write the tecplot datafile header information.
        ios = TecIni142('TDCG-program Results'//NULLCHR, &
                        'X Y Z Cross WallDistance'//NULLCHR, &
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
                        nNode, &   ! NumPts
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
                        [1,1,1,0,0], & ! VarLocation
                        Null, &     ! ShareVarFromZone
                        0)          ! ShareConnectivityFromZone
            if (ios/=0) stop "Error value returned in TecZne142"

        do i=1,3
        ios = TecDat142(nNode,real(Node(1:nNode,i),R4),0)
            if (ios/=0) stop "Error value returned in TecDat142"
        enddo
        DEALLOCATE(Node)

        ALLOCATE(cVariables(nCells,2))
        call initTmpStorageVar(MeshOnly)  ! Temporary Storage cVariables
        do i=1,2
        ios = TecDat142(nCells,real(cVariables(1:nCells,i),R4),0)
            if (ios/=0) stop "Error value returned in TecDat142"
        enddo
        DEALLOCATE(cVariables)
    else
        ! Open the file and write the tecplot datafile header information.
        ios = TecIni142('TDCG-program Results'//NULLCHR, &
                        'X Y Z U V W Rou T P Ma Cross WallDistance'//NULLCHR, &
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
                        nNode, &   ! NumPts
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
                        [1,1,1,0,0,0,0,0,0,0,0,0], & ! VarLocation
                        Null, &     ! ShareVarFromZone
                        0)          ! ShareConnectivityFromZone
            if (ios/=0) stop "Error value returned in TecZne142"

        do i=1,3
        ios = TecDat142(nNode,real(Node(1:nNode,i),R4),0)
            if (ios/=0) stop "Error value returned in TecDat142"
        enddo
        DEALLOCATE(Node)

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
    use ModNeighbor
    use ModOutput
    use ModTypDef
    implicit none
    type(FTTCell),pointer :: t,cc
    integer :: i, j, k,a
    real(R8):: dx,dy,dz,tN(6)
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    
    nNode=8
    
    cc=>OctCell(1,1,1)
    cc%cNode(1)=1
    cc%cNode(2)=2
    cc%cNode(3)=3
    cc%cNode(4)=4
    cc%cNode(5)=5
    cc%cNode(6)=6
    cc%cNode(7)=7
    cc%cNode(8)=8
    
    dx=BGCellSize(1)/2**(cc%LVL(1)+1)
    dy=BGCellSize(2)/2**(cc%LVL(2)+1)
    dz=BGCellSize(3)/2**(cc%LVL(3)+1)
    tN(1)=cc%Center(1)-dx   ! x -
    tN(2)=cc%Center(2)-dy   ! y -
    tN(3)=cc%Center(3)-dz   ! z -
    tN(4)=cc%Center(1)+dx   ! x +
    tN(5)=cc%Center(2)+dy   ! y +
    tN(6)=cc%Center(3)+dz   ! z +
    
    Node(1,1)=tN(1)
    Node(1,2)=tN(2)
    Node(1,3)=tN(3)
    
    Node(2,1)=tN(4)
    Node(2,2)=tN(2)
    Node(2,3)=tN(3)
    
    Node(3,1)=tN(4)
    Node(3,2)=tN(5)
    Node(3,3)=tN(3)
    
    Node(4,1)=tN(1)
    Node(4,2)=tN(5)
    Node(4,3)=tN(3)
    
    Node(5,1)=tN(1)
    Node(5,2)=tN(2)
    Node(5,3)=tN(6)
    
    Node(6,1)=tN(4)
    Node(6,2)=tN(2)
    Node(6,3)=tN(6)
    
    Node(7,1)=tN(4)
    Node(7,2)=tN(5)
    Node(7,3)=tN(6)
    
    Node(8,1)=tN(1)
    Node(8,2)=tN(5)
    Node(8,3)=tN(6)
    
    
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        if(i==1.and.j==1.and.k==1)then
            goto 66
        else
        t=>OctCell(i, j, k)
        call InitialValueNodeInfo(t)
66        a=1
      endif
    enddo
    enddo
    enddo
    
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        if(i==1.and.j==1.and.k==1)then
            goto 68
        else
        t=>OctCell(i, j, k)
        call NodeInfo(t)
68        a=1
      endif
    enddo
    enddo
    enddo
    
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "InitNodeInfo time: ", tEnd-tStart

    contains
!----------------------------------------------------------------------
     recursive subroutine InitialValueNodeInfo(c)
        implicit none
        type(FTTCell),pointer :: c,cs
        integer:: is
        
        if(ASSOCIATED(c%Octson))then
            do is=1,c%Octson%nSon
                cs=>c%Octson%son(is)
                call InitialValueNodeInfo(cs)
            enddo
        else
         c%cnode(1)=-1
         c%cnode(2)=-1
         c%cnode(3)=-1
         c%cnode(4)=-1
         c%cnode(5)=-1
         c%cnode(6)=-1
         c%cnode(7)=-1
         c%cnode(8)=-1
        endif
       return 
     endsubroutine InitialValueNodeInfo
!----------------------------------------------------------------------    
        recursive subroutine NodeInfo(c)
        implicit none
        type(FTTCell),pointer :: c,cX1,cX2,cY1,cZ1
        real(R8):: dx, dy, dz, tN(6) ! Temp-Nodes
        type(FTTCell),pointer :: cs
        integer:: is

        if(ASSOCIATED(c%Octson))then
            do is=1,c%Octson%nSon
                cs=>c%Octson%son(is)
                call NodeInfo(cs)
            enddo
            return
        endif

        
        if(c%Location==0)then
            if(c%nBGCell(1)==1 .and. c%nBGCell(2)==1)then
                cZ1=>NeighborZ1(c)
                c%cnode(1)=cZ1%cnode(5)   
                c%cnode(2)=cZ1%cnode(6)
                c%cnode(3)=cZ1%cnode(7)
                c%cnode(4)=cZ1%cnode(8)
                
                nNode=nNode+1
                c%cnode(5)=nNode
                Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
                
                nNode=nNode+1
                c%cnode(6)=nNode
                Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
                
                nNode=nNode+1
                c%cnode(7)=nNode
                Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
                
                nNode=nNode+1
                c%cnode(8)=nNode
                Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            elseif(c%nBGCell(2)==1)then
                cX1=>NeighborX1(c)
                if(ASSOCIATED(cX1%Octson))then
                    call corner(cX1)
                endif
                
                if(c%nBGCell(2)==1 .and. c%nBGCell(3)==1)then
                    c%cnode(1)=cX1%cnode(2)
                    c%cnode(4)=cX1%cnode(3)
                    c%cnode(5)=cX1%cnode(6)
                    c%cnode(8)=cX1%cnode(7)
                    
                    nNode=nNode+1
                    c%cnode(2)=nNode
                    Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
                    
                    nNode=nNode+1
                    c%cnode(3)=nNode
                    Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
                    
                    nNode=nNode+1
                    c%cnode(6)=nNode
                    Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
                    
                    nNode=nNode+1
                    c%cnode(7)=nNode
                    Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
                else
                    cZ1=>NeighborZ1(c)
                    
                    c%cnode(1)=cZ1%cnode(5)
                    c%cnode(2)=cZ1%cnode(6)
                    c%cnode(3)=cZ1%cnode(7)
                    c%cnode(4)=cZ1%cnode(8)
    
                    c%cnode(5)=cX1%cnode(6)
                    c%cnode(8)=cX1%cnode(7)
                    
                    nNode=nNode+1
                    c%cnode(6)=nNode
                    Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
                    
                    nNode=nNode+1
                    c%cnode(7)=nNode
                    Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
                endif
            elseif(c%nBGCell(1)==1)then
                cY1=>NeighborY1(c)
                if(ASSOCIATED(cY1%Octson))then
                    call corner(cY1)
                endif
                
                if(c%nBGCell(1)==1 .and. c%nBGCell(3)==1)then
                    c%cnode(1)=cY1%cnode(4)
                    c%cnode(2)=cY1%cnode(3)
                    c%cnode(5)=cY1%cnode(8)
                    c%cnode(6)=cY1%cnode(7)
                    
                    nNode=nNode+1
                    c%cnode(3)=nNode
                    Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
                    
                    nNode=nNode+1
                    c%cnode(4)=nNode
                    Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
                    
                    nNode=nNode+1
                    c%cnode(7)=nNode
                    Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
                    
                    nNode=nNode+1
                    c%cnode(8)=nNode
                    Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
                else
                    cZ1=>NeighborZ1(c)
                    
                    c%cnode(1)=cZ1%cnode(5)
                    c%cnode(2)=cZ1%cnode(6)
                    c%cnode(3)=cZ1%cnode(7)
                    c%cnode(4)=cZ1%cnode(8)
                    
                    c%cnode(5)=cY1%cnode(8)
                    c%cnode(6)=cY1%cnode(7)
                    
                    nNode=nNode+1
                    c%cnode(7)=nNode
                    Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
                    
                    nNode=nNode+1
                    c%cnode(8)=nNode
                    Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
                endif
            elseif(c%nBGCell(3)==1)then
                cY1=>NeighborY1(c)
                cX1=>NeighborX1(c)
                
                if(ASSOCIATED(cY1%Octson))then
                    call corner(cY1)
                endif
                if(ASSOCIATED(cX1%Octson))then
                    call corner(cX1)
                endif
                
                c%cnode(1)=cY1%cnode(4)
                c%cnode(2)=cY1%cnode(3)
                c%cnode(5)=cY1%cnode(8)
                c%cnode(6)=cY1%cnode(7)
                
                c%cnode(4)=cX1%cnode(3)
                c%cnode(8)=cX1%cnode(7)
                
                nNode=nNode+1
                c%cnode(3)=nNode
                Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
                
                nNode=nNode+1
                c%cnode(7)=nNode
                Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            else
                cZ1=>NeighborZ1(c)
                cY1=>NeighborY1(c)
                cX1=>NeighborX1(c)
                
                if(ASSOCIATED(cZ1%Octson))then
                    call corner(cZ1)
                endif
                if(ASSOCIATED(cY1%Octson))then
                    call corner(cY1)
                endif
                if(ASSOCIATED(cX1%Octson))then
                    call corner(cX1)
                endif
                
                c%cnode(1)=cZ1%cnode(5)
                c%cnode(2)=cZ1%cnode(6)
                c%cnode(3)=cZ1%cnode(7)
                c%cnode(4)=cZ1%cnode(8)
                c%cnode(5)=cY1%cnode(8)
                c%cnode(6)=cY1%cnode(7)
                c%cnode(8)=cX1%cnode(7)
                
                nNode=nNode+1
                c%cnode(7)=nNode
                Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif
        elseif(c%Location==1)then                           !son1
            cZ1=>NeighborZ1(c)
            cY1=>NeighborY1(c)
            cX1=>NeighborX1(c)
            
            if(ASSOCIATED(cZ1))then
            if(ASSOCIATED(cZ1%Octson))then
                call corner(cZ1)
            endif
            endif
            if(ASSOCIATED(cY1))then
            if(ASSOCIATED(cY1%Octson))then
                call corner(cY1)
            endif
            endif
            if(ASSOCIATED(cX1))then
            if(ASSOCIATED(cX1%Octson))then
                call corner(cX1)
            endif
            endif
            
            if(ASSOCIATED(cZ1) .and.cZ1%cnode(5)/=-1.and.cZ1%LVL(1)+1>=c%LVL(1)&  ! node1
               &.and. cZ1%LVL(2)+1>=c%LVL(2).and. cZ1%LVL(3)+1>=c%LVL(3))then
                c%cnode(1)=cZ1%cnode(5)                    
            elseif(ASSOCIATED(cY1) .and.cY1%cnode(4)/=-1.and.cY1%LVL(1)+1>=c%LVL(1)& 
               &.and. cY1%LVL(2)+1>=c%LVL(2).and. cY1%LVL(3)+1>=c%LVL(3))then
                c%cnode(1)=cY1%cnode(4)
            elseif(ASSOCIATED(cX1) .and.cX1%cnode(2)/=-1.and.cX1%LVL(1)+1>=c%LVL(1)& 
               &.and. cX1%LVL(2)+1>=c%LVL(2).and. cX1%LVL(3)+1>=c%LVL(3))then
                c%cnode(1)=cX1%cnode(2) 
            else
                nNode=nNode+1
                c%cnode(1)=nNode
                Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
                
           if(ASSOCIATED(cZ1).and. cZ1%cnode(6)/=-1.and.cZ1%LVL(1)==c%LVL(1)&  ! node2
               &.and. cZ1%LVL(2)==c%LVL(2).and. cZ1%LVL(3)==c%LVL(3))then
               c%cnode(2)=cZ1%cnode(6)
           elseif(ASSOCIATED(cY1).and. cY1%cnode(3)/=-1.and.cY1%LVL(1)==c%LVL(1)&
               &.and. cY1%LVL(2)==c%LVL(2).and. cY1%LVL(3)==c%LVL(3))then
               c%cnode(2)=cY1%cnode(3)
           else
                nNode=nNode+1
                c%cnode(2)=nNode
                Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
           endif
           
           if(ASSOCIATED(cZ1).and. cZ1%cnode(7)/=-1.and.cZ1%LVL(1)==c%LVL(1)&  ! node3
               &.and. cZ1%LVL(2)==c%LVL(2).and. cZ1%LVL(3)==c%LVL(3))then
               c%cnode(3)=cZ1%cnode(7)
           else
                nNode=nNode+1
                c%cnode(3)=nNode
                Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1) 
           endif
           
            
           if(ASSOCIATED(cZ1).and.cZ1%cnode(8)/=-1.and.cZ1%LVL(1)==c%LVL(1)& ! node4
               &.and. cZ1%LVL(2)==c%LVL(2).and. cZ1%LVL(3)==c%LVL(3))then
                    c%cnode(4)=cZ1%cnode(8)   
           elseif(ASSOCIATED(cX1).and.cX1%cnode(3)/=-1.and.cX1%LVL(1)==c%LVL(1)& 
               &.and. cX1%LVL(2)==c%LVL(2).and. cX1%LVL(3)==c%LVL(3))then
                    c%cnode(4)=cX1%cnode(3) 
           else
                    nNode=nNode+1
                    c%cnode(4)=nNode
                    Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
                    Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
                    Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)  
           endif
           if(ASSOCIATED(cY1).and. cY1%cnode(8)/=-1.and.cY1%LVL(1)==c%LVL(1)&  ! node5
               &.and. cY1%LVL(2)==c%LVL(2).and. cY1%LVL(3)==c%LVL(3))then
               c%cnode(5)=cY1%cnode(8) 
           elseif(ASSOCIATED(cX1).and. cX1%cnode(6)/=-1.and.cX1%LVL(1)==c%LVL(1)& 
               &.and. cX1%LVL(2)==c%LVL(2).and. cX1%LVL(3)==c%LVL(3))then
               c%cnode(5)=cX1%cnode(6) 
           else
                nNode=nNode+1
                c%cnode(5)=nNode
                Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1) 
           endif
           
           if(ASSOCIATED(cY1).and. cY1%cnode(7)/=-1.and.cY1%LVL(1)==c%LVL(1)&  ! node6
               &.and. cY1%LVL(2)==c%LVL(2).and. cY1%LVL(3)==c%LVL(3))then   
               c%cnode(6)=cY1%cnode(7) 
           else
                nNode=nNode+1
                c%cnode(6)=nNode
                Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
                Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
                Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
           endif
           
           nNode=nNode+1                                                   ! node7
           c%cnode(7)=nNode
           Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
           Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
           Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
           
           if(ASSOCIATED(cX1).and. cX1%cnode(7)/=-1.and.cX1%LVL(1)==c%LVL(1)&  ! node8
               &.and. cX1%LVL(2)==c%LVL(2).and. cX1%LVL(3)==c%LVL(3))then
               c%cnode(8)=cX1%cnode(7)
           else
               nNode=nNode+1
               c%cnode(8)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
           endif
        elseif(c%Location==2)then                       !son2
           cZ1=> NeighborZ1(c)
           cY1=> NeighborY1(c)
           cX1=> NeighborX1(c)
           cX2=> NeighborX2(c)
           
            if(ASSOCIATED(cZ1))then
            if(ASSOCIATED(cZ1%Octson))then
                call corner(cZ1)
            endif
            endif
            if(ASSOCIATED(cY1))then
            if(ASSOCIATED(cY1%Octson))then
                call corner(cY1)
            endif
            endif
            if(ASSOCIATED(cX1))then
            if(ASSOCIATED(cX1%Octson))then
                call corner(cX1)
            endif
            endif
            if(ASSOCIATED(cX2))then
            if(ASSOCIATED(cX2%Octson))then
                call corner(cX2)
            endif
            endif
            
            if(ASSOCIATED(cX1).and. cX1%cnode(2)/=-1.and.cX1%LVL(1)==c%LVL(1)& 
               &.and. cX1%LVL(2)==c%LVL(2).and. cX1%LVL(3)==c%LVL(3))then              ! node1
                c%cnode(1)=cX1%cnode(2)
            else
               nNode=nNode+1
               c%cnode(1)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(6)/=-1.and.cZ1%LVL(1)+1>=c%LVL(1)& 
               &.and. cZ1%LVL(2)+1>=c%LVL(2).and. cZ1%LVL(3)+1>=c%LVL(3))then              ! node2
                c%cnode(2)=cZ1%cnode(6)
            elseif(ASSOCIATED(cY1).and. cY1%cnode(3)/=-1.and.cY1%LVL(1)+1>=c%LVL(1)& 
               &.and. cY1%LVL(2)+1>=c%LVL(2).and. cY1%LVL(3)+1>=c%LVL(3))then 
                c%cnode(2)=cY1%cnode(3)
            else
               nNode=nNode+1
               c%cnode(2)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
               
            if(ASSOCIATED(cZ1).and. cZ1%cnode(7)/=-1.and.cZ1%LVL(1)==c%LVL(1)&  ! node3
               &.and. cZ1%LVL(2)==c%LVL(2).and. cZ1%LVL(3)==c%LVL(3))then
               c%cnode(3)=cZ1%cnode(7)
            elseif(ASSOCIATED(cX2).and. cX2%cnode(4)/=-1.and.cX2%LVL(1)==c%LVL(1)& 
               &.and. cX2%LVL(2)==c%LVL(2).and. cX2%LVL(3)==c%LVL(3))then
               c%cnode(3)=cX2%cnode(4)  
            else
               nNode=nNode+1
               c%cnode(3)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cX1).and. cX1%cnode(3)/=-1.and.cX1%LVL(1)==c%LVL(1)& 
               &.and. cX1%LVL(2)==c%LVL(2).and. cX1%LVL(3)==c%LVL(3))then                 !node4
                c%cnode(4)=cX1%cnode(3)  
            else
               nNode=nNode+1
               c%cnode(4)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cX1).and. cX1%cnode(6)/=-1.and.cX1%LVL(1)==c%LVL(1)& 
               &.and. cX1%LVL(2)==c%LVL(2).and. cX1%LVL(3)==c%LVL(3))then                 !node5
                c%cnode(5)=cX1%cnode(6)  
            else
               nNode=nNode+1
               c%cnode(5)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(7)/=-1.and.cY1%LVL(1)==c%LVL(1)&  ! node6
               &.and. cY1%LVL(2)==c%LVL(2).and. cY1%LVL(3)==c%LVL(3))then
               c%cnode(6)=cY1%cnode(7)
            elseif(ASSOCIATED(cX2).and. cX2%cnode(5)/=-1.and.cX2%LVL(1)==c%LVL(1)& 
               &.and. cX2%LVL(2)==c%LVL(2).and. cX2%LVL(3)==c%LVL(3))then
               c%cnode(6)=cX2%cnode(5)  
            else
               nNode=nNode+1
               c%cnode(6)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cX2).and. cX2%cnode(8)/=-1.and.cX2%LVL(1)==c%LVL(1)& ! node7
               &.and. cX2%LVL(2)==c%LVL(2).and. cX2%LVL(3)==c%LVL(3))then
               c%cnode(7)=cX2%cnode(8)  
            else
               nNode=nNode+1
               c%cnode(7)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif  
            
            if(ASSOCIATED(cX1).and. cX1%cnode(7)/=-1.and.cX1%LVL(1)==c%LVL(1)& ! node8
               &.and. cX1%LVL(2)==c%LVL(2).and. cX1%LVL(3)==c%LVL(3))then
               c%cnode(8)=cX1%cnode(7)  
            else
               nNode=nNode+1
               c%cnode(8)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif 

    elseif(c%Location==3)then                                      !son3
            cZ1=>NeighborZ1(c)
            cY1=>NeighborY1(c)
            cX2=>NeighborX2(c)
            
            if(ASSOCIATED(cZ1))then
            if(ASSOCIATED(cZ1%Octson))then
                call corner(cZ1)
            endif
            endif
            if(ASSOCIATED(cY1))then
            if(ASSOCIATED(cY1%Octson))then
                call corner(cY1)
            endif
            endif
            if(ASSOCIATED(cX2))then
            if(ASSOCIATED(cX2%Octson))then
                call corner(cX2)
            endif
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(4)/=-1)then         !node1
                c%cnode(1)=cY1%cnode(4)  
            else
               nNode=nNode+1
               c%cnode(1)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(3)/=-1)then         !node2
                c%cnode(2)=cY1%cnode(3)  
            else
               nNode=nNode+1
               c%cnode(2)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
               
            if(ASSOCIATED(cZ1).and. cZ1%cnode(7)/=-1.and.cZ1%LVL(1)+1>=c%LVL(1)& 
               &.and. cZ1%LVL(2)+1>=c%LVL(2).and. cZ1%LVL(3)+1>=c%LVL(3))then         !node3
                c%cnode(3)=cZ1%cnode(7)  
            elseif(ASSOCIATED(cX2).and. cX2%cnode(4)/=-1.and.cX2%LVL(1)+1>=c%LVL(1)& 
               &.and. cX2%LVL(2)+1>=c%LVL(2).and. cX2%LVL(3)+1>=c%LVL(3))then 
                c%cnode(3)=cX2%cnode(4)  
            else
               nNode=nNode+1
               c%cnode(3)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif  
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(8)/=-1.and.cZ1%LVL(1)==c%LVL(1)&  ! node4
               &.and. cZ1%LVL(2)==c%LVL(2).and. cZ1%LVL(3)==c%LVL(3))then
               c%cnode(4)=cZ1%cnode(8)
            else
               nNode=nNode+1
               c%cnode(4)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)  
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(8)/=-1)then         !node5
                c%cnode(5)=cY1%cnode(8)  
            else
               nNode=nNode+1
               c%cnode(5)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif  
             
            if(ASSOCIATED(cY1).and. cY1%cnode(7)/=-1)then         !node6
                c%cnode(6)=cY1%cnode(7)  
            else
               nNode=nNode+1
               c%cnode(6)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif 
               
            if(ASSOCIATED(cX2).and. cX2%cnode(8)/=-1.and.cX2%LVL(1)==c%LVL(1)&  ! node7
               &.and. cX2%LVL(2)==c%LVL(2).and. cX2%LVL(3)==c%LVL(3))then
               c%cnode(7)=cX2%cnode(8)
            else
               nNode=nNode+1
               c%cnode(7)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)  
            endif
            
            nNode=nNode+1                                         ! node8
            c%cnode(8)=nNode
            Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
            Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
            Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)    
        elseif(c%Location==4)then                             !son4
            cZ1=>NeighborZ1(c)
            cY1=>NeighborY1(c)
            cX1=>NeighborX1(c)
            cX2=>NeighborX2(c)
            
            if(ASSOCIATED(cZ1))then
            if(ASSOCIATED(cZ1%Octson))then
                call corner(cZ1)
            endif
            endif
            if(ASSOCIATED(cY1))then
            if(ASSOCIATED(cY1%Octson))then
                call corner(cY1)
            endif
            endif
            if(ASSOCIATED(cX1))then
            if(ASSOCIATED(cX1%Octson))then
                call corner(cX1)
            endif
            endif
            if(ASSOCIATED(cX2))then
            if(ASSOCIATED(cX2%Octson))then
                call corner(cX2)
            endif
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(4)/=-1)then         !node1
                c%cnode(1)=cY1%cnode(4)  
            else
               nNode=nNode+1
               c%cnode(1)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif 
            
            if(ASSOCIATED(cY1).and. cY1%cnode(3)/=-1)then         !node2
                c%cnode(2)=cY1%cnode(3)  
            else
               nNode=nNode+1
               c%cnode(2)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cX2).and. cX2%cnode(4)/=-1)then         !node3
                c%cnode(3)=cX2%cnode(4)  
            else
               nNode=nNode+1
               c%cnode(3)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(8)/=-1.and.cZ1%LVL(1)+1>=c%LVL(1)& 
               &.and. cZ1%LVL(2)+1>=c%LVL(2).and. cZ1%LVL(3)+1>=c%LVL(3))then         !node4
                c%cnode(4)=cZ1%cnode(8)  
            elseif(ASSOCIATED(cX1).and. cX1%cnode(3)/=-1.and.cX1%LVL(1)+1>=c%LVL(1)& 
               &.and. cX1%LVL(2)+1>=c%LVL(2).and. cX1%LVL(3)+1>=c%LVL(3))then
                c%cnode(4)=cX1%cnode(3) 
            else
               nNode=nNode+1
               c%cnode(4)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(8)/=-1.and.cY1%LVL(1)==c%LVL(1)& 
               &.and. cY1%LVL(2)==c%LVL(2).and. cY1%LVL(3)==c%LVL(3))then         !node5
                c%cnode(5)=cY1%cnode(8)  
            else
               nNode=nNode+1
               c%cnode(5)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(7)/=-1.and.cY1%LVL(1)==c%LVL(1)&  ! node6
               &.and. cY1%LVL(2)==c%LVL(2).and. cY1%LVL(3)==c%LVL(3))then
               c%cnode(6)=cY1%cnode(7)
            else
               nNode=nNode+1
               c%cnode(6)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)  
            endif
            
            if(ASSOCIATED(cX2).and. cX2%cnode(8)/=-1.and.cX2%LVL(1)==c%LVL(1)&  ! node7
               &.and. cX2%LVL(2)==c%LVL(2).and. cX2%LVL(3)==c%LVL(3))then
               c%cnode(7)=cX2%cnode(8)
            else
               nNode=nNode+1
               c%cnode(7)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)  
            endif
            
            if(ASSOCIATED(cX1).and. cX1%cnode(7)/=-1.and.cX1%LVL(1)==c%LVL(1)&  ! node8
               &.and. cX1%LVL(2)==c%LVL(2).and. cX1%LVL(3)==c%LVL(3))then
               c%cnode(8)=cX1%cnode(7)
            else
               nNode=nNode+1
               c%cnode(8)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)  
            endif
        elseif(c%Location==5)then                                         !son5
            cZ1=>NeighborZ1(c)
            cY1=>NeighborY1(c)
            cX1=>NeighborX1(c)
            cX2=>NeighborX2(c)
            
            if(ASSOCIATED(cZ1))then
            if(ASSOCIATED(cZ1%Octson))then
                call corner(cZ1)
            endif
            endif
            if(ASSOCIATED(cY1))then
            if(ASSOCIATED(cY1%Octson))then
                call corner(cY1)
            endif
            endif
            if(ASSOCIATED(cX1))then
            if(ASSOCIATED(cX1%Octson))then
                call corner(cX1)
            endif
            endif
            if(ASSOCIATED(cX2))then
            if(ASSOCIATED(cX2%Octson))then
                call corner(cX2)
            endif
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(5)/=-1)then         !node1
                c%cnode(1)=cZ1%cnode(5)  
            else
               nNode=nNode+1
               c%cnode(1)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(6)/=-1)then         !node2
                c%cnode(2)=cZ1%cnode(6)  
            else
               nNode=nNode+1
               c%cnode(2)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(7)/=-1)then         !node3
                c%cnode(3)=cZ1%cnode(7)  
            else
               nNode=nNode+1
               c%cnode(3)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(8)/=-1)then         !node4
                c%cnode(4)=cZ1%cnode(8)  
            else
               nNode=nNode+1
               c%cnode(4)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(8)/=-1.and.cY1%LVL(1)+1>=c%LVL(1)& 
               &.and. cY1%LVL(2)+1>=c%LVL(2).and. cY1%LVL(3)+1>=c%LVL(3))then         !node5
                c%cnode(5)=cY1%cnode(8)  
            elseif(ASSOCIATED(cX1).and. cX1%cnode(6)/=-1.and.cX1%LVL(1)+1>=c%LVL(1)& 
               &.and. cX1%LVL(2)+1>=c%LVL(2).and. cX1%LVL(3)+1>=c%LVL(3))then      
                c%cnode(5)=cX1%cnode(6)      
            else
               nNode=nNode+1
               c%cnode(5)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(7)/=-1.and.cY1%LVL(1)==c%LVL(1)&  ! node6
               &.and. cY1%LVL(2)==c%LVL(2).and. cY1%LVL(3)==c%LVL(3))then
               c%cnode(6)=cY1%cnode(7)
            else
               nNode=nNode+1
               c%cnode(6)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)  
            endif
            
            nNode=nNode+1                                                 ! node7
            c%cnode(7)=nNode
            Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
            Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
            Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            
            if(ASSOCIATED(cX1).and. cX1%cnode(7)/=-1.and.cX1%LVL(1)==c%LVL(1)&  ! node8
               &.and. cX1%LVL(2)==c%LVL(2).and. cX1%LVL(3)==c%LVL(3))then
               c%cnode(8)=cX1%cnode(7)
            else
               nNode=nNode+1
               c%cnode(8)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)  
            endif
        elseif(c%Location==6)then                                       !son1
            cZ1=>NeighborZ1(c)
            cY1=>NeighborY1(c)
            cX1=>NeighborX1(c)
            cX2=>NeighborX2(c)
            
            if(ASSOCIATED(cZ1))then
            if(ASSOCIATED(cZ1%Octson))then
                call corner(cZ1)
            endif
            endif
            if(ASSOCIATED(cY1))then
            if(ASSOCIATED(cY1%Octson))then
                call corner(cY1)
            endif
            endif
            if(ASSOCIATED(cX1))then
            if(ASSOCIATED(cX1%Octson))then
                call corner(cX1)
            endif
            endif
            if(ASSOCIATED(cX2))then
            if(ASSOCIATED(cX2%Octson))then
                call corner(cX2)
            endif
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(5)/=-1.and.cZ1%LVL(1)==c%LVL(1)& 
               &.and. cZ1%LVL(2)==c%LVL(2).and. cZ1%LVL(3)==c%LVL(3))then         !node1
                c%cnode(1)=cZ1%cnode(5)  
            else
               nNode=nNode+1
               c%cnode(1)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(6)/=-1.and.cZ1%LVL(1)==c%LVL(1)& 
               &.and. cZ1%LVL(2)==c%LVL(2).and. cZ1%LVL(3)==c%LVL(3))then         !node2
                c%cnode(2)=cZ1%cnode(6)  
            else
               nNode=nNode+1
               c%cnode(2)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(7)/=-1.and.cZ1%LVL(1)==c%LVL(1)& 
               &.and. cZ1%LVL(2)==c%LVL(2).and. cZ1%LVL(3)==c%LVL(3))then         !node3
                c%cnode(3)=cZ1%cnode(7)  
            else
               nNode=nNode+1
               c%cnode(3)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(8)/=-1.and.cZ1%LVL(1)==c%LVL(1)& 
               &.and. cZ1%LVL(2)==c%LVL(2).and. cZ1%LVL(3)==c%LVL(3))then         !node4
                c%cnode(4)=cZ1%cnode(8)  
            else
               nNode=nNode+1
               c%cnode(4)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cX1).and. cX1%cnode(6)/=-1.and.cX1%LVL(1)==c%LVL(1)& 
               &.and. cX1%LVL(2)==c%LVL(2).and. cX1%LVL(3)==c%LVL(3))then         !node5
                c%cnode(5)=cX1%cnode(6)  
            else
               nNode=nNode+1
               c%cnode(5)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(7)/=-1.and.cY1%LVL(1)+1>=c%LVL(1)& 
               &.and. cY1%LVL(2)+1>=c%LVL(2).and. cY1%LVL(3)+1>=c%LVL(3))then         !node6
                c%cnode(6)=cY1%cnode(7)  
            elseif(ASSOCIATED(cX2).and. cX2%cnode(5)/=-1.and.cX2%LVL(1)+1>=c%LVL(1)& 
               &.and. cX2%LVL(2)+1>=c%LVL(2).and. cX2%LVL(3)+1>=c%LVL(3))then
                c%cnode(6)=cX2%cnode(5)  
            else
               nNode=nNode+1
               c%cnode(6)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cX2).and. cX2%cnode(8)/=-1.and.cX2%LVL(1)==c%LVL(1)&  ! node7
               &.and. cX2%LVL(2)==c%LVL(2).and. cX2%LVL(3)==c%LVL(3))then
               c%cnode(7)=cX2%cnode(8)
            else
               nNode=nNode+1
               c%cnode(7)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)  
            endif

            if(ASSOCIATED(cX1).and. cX1%cnode(7)/=-1.and.cX1%LVL(1)==c%LVL(1)&  ! node8
               &.and. cX1%LVL(2)==c%LVL(2).and. cX1%LVL(3)==c%LVL(3))then
               c%cnode(8)=cX1%cnode(7)
            else
               nNode=nNode+1
               c%cnode(8)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)  
            endif
        elseif(c%Location==7)then
            cZ1=>NeighborZ1(c)
            cY1=>NeighborY1(c)
            cX1=>NeighborX1(c)
            cX2=>NeighborX2(c)

            if(ASSOCIATED(cZ1))then
            if(ASSOCIATED(cZ1%Octson))then
                call corner(cZ1)
            endif
            endif
            if(ASSOCIATED(cY1))then
            if(ASSOCIATED(cY1%Octson))then
                call corner(cY1)
            endif
            endif
            if(ASSOCIATED(cX1))then
            if(ASSOCIATED(cX1%Octson))then
                call corner(cX1)
            endif
            endif
            if(ASSOCIATED(cX2))then
            if(ASSOCIATED(cX2%Octson))then
                call corner(cX2)
            endif
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(5)/=-1)then         !node1
                c%cnode(1)=cZ1%cnode(5)  
            else
               nNode=nNode+1
               c%cnode(1)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(6)/=-1)then         !node2
                c%cnode(2)=cZ1%cnode(6)  
            else
               nNode=nNode+1
               c%cnode(2)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(7)/=-1)then         !node3
                c%cnode(3)=cZ1%cnode(7)  
            else
               nNode=nNode+1
               c%cnode(3)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(8)/=-1)then         !node4
                c%cnode(4)=cZ1%cnode(8)  
            else
               nNode=nNode+1
               c%cnode(4)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(8)/=-1)then         !node5
                c%cnode(5)=cY1%cnode(8)  
            else
               nNode=nNode+1
               c%cnode(5)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cY1).and. cY1%cnode(7)/=-1)then         !node6
                c%cnode(6)=cY1%cnode(7)  
            else
               nNode=nNode+1
               c%cnode(6)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            if(ASSOCIATED(cX2).and. cX2%cnode(8)/=-1.and.cX2%LVL(1)+1>=c%LVL(1)& 
               &.and. cX2%LVL(2)+1>=c%LVL(2).and. cX2%LVL(3)+1>=c%LVL(3))then         !node7
                c%cnode(7)=cX2%cnode(8)  
            else
               nNode=nNode+1
               c%cnode(7)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif
            
            nNode=nNode+1                                     !node8
            c%cnode(8)=nNode
            Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
            Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
            Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
        elseif(c%Location==8)then
            cZ1=>NeighborZ1(c)
            cY1=>NeighborY1(c)
            cX1=>NeighborX1(c)
            cX2=>NeighborX2(c)
            
            if(ASSOCIATED(cZ1))then
            if(ASSOCIATED(cZ1%Octson))then
                call corner(cZ1)
            endif
            endif
            if(ASSOCIATED(cY1))then
            if(ASSOCIATED(cY1%Octson))then
                call corner(cY1)
            endif
            endif
            if(ASSOCIATED(cX1))then
            if(ASSOCIATED(cX1%Octson))then
                call corner(cX1)
            endif
            endif
            if(ASSOCIATED(cX2))then
            if(ASSOCIATED(cX2%Octson))then
                call corner(cX2)
            endif
            endif
            
            if(ASSOCIATED(cZ1).and. cZ1%cnode(5)/=-1)then          !node1
                c%cnode(1)=cZ1%cnode(5)  
            else
               nNode=nNode+1
               c%cnode(1)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif

            if(ASSOCIATED(cZ1).and. cZ1%cnode(6)/=-1)then          !node2
                c%cnode(2)=cZ1%cnode(6)  
            else
               nNode=nNode+1
               c%cnode(2)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif

            if(ASSOCIATED(cZ1).and. cZ1%cnode(7)/=-1)then          !node3
                c%cnode(3)=cZ1%cnode(7)  
            else
               nNode=nNode+1
               c%cnode(3)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif

            if(ASSOCIATED(cZ1).and. cZ1%cnode(8)/=-1)then          !node4
                c%cnode(4)=cZ1%cnode(8)  
            else
               nNode=nNode+1
               c%cnode(4)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
            endif

            if(ASSOCIATED(cY1).and. cY1%cnode(8)/=-1)then          !node5
                c%cnode(5)=cY1%cnode(8)  
            else
               nNode=nNode+1
               c%cnode(5)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif

            if(ASSOCIATED(cY1).and. cY1%cnode(7)/=-1)then          !node6
                c%cnode(6)=cY1%cnode(7)  
            else
               nNode=nNode+1
               c%cnode(6)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)-BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif

            if(ASSOCIATED(cX2).and. cX2%cnode(8)/=-1)then          !node7
                c%cnode(7)=cX2%cnode(8)  
            else
               nNode=nNode+1
               c%cnode(7)=nNode
               Node(nNode,1)=c%Center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif

            if(ASSOCIATED(cX1).and. cX1%cnode(7)/=-1.and.cX1%LVL(1)+1>=c%LVL(1)& 
               &.and. cX1%LVL(2)+1>=c%LVL(2).and. cX1%LVL(3)+1>=c%LVL(3))then          !node8
                c%cnode(8)=cX1%cnode(7)  
            else
               nNode=nNode+1
               c%cnode(8)=nNode
               Node(nNode,1)=c%Center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
               Node(nNode,2)=c%Center(2)+BGCellSize(2)/2**(c%LVL(2)+1) 
               Node(nNode,3)=c%Center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
            endif
        endif
        endsubroutine NodeInfo
!----------------------------------------------------------------------
        recursive subroutine corner(c)
        implicit none
        type(FTTCell),pointer :: c,c1,c2,c3,c4,c5,c6,c7,c8
        
        if(associated(c%Octson%son(1)%Octson))then
            c1=>c%Octson%son(1)
            call corner(c1)
        endif
        
        if(associated(c%Octson%son(2)%Octson))then
            c2=>c%Octson%son(2)
            call corner(c2)
        endif
        
        if(associated(c%Octson%son(3)%Octson))then
            c3=>c%Octson%son(3)
            call corner(c3)
        endif
        
        if(associated(c%Octson%son(4)%Octson))then
            c4=>c%Octson%son(4)
            call corner(c4)
        endif
        
        if(associated(c%Octson%son(5)%Octson))then
            c5=>c%Octson%son(5)
            call corner(c5)
        endif
        
        if(associated(c%Octson%son(6)%Octson))then
            c6=>c%Octson%son(6)
            call corner(c6)
        endif
        
        if(associated(c%Octson%son(7)%Octson))then
            c7=>c%Octson%son(7)
            call corner(c7)
        endif
        
        if(associated(c%Octson%son(8)%Octson))then
            c8=>c%Octson%son(8)
            call corner(c8)
        endif
        
        c%cnode(1)=c%Octson%son(1)%cnode(1)
        c%cnode(2)=c%Octson%son(2)%cnode(2)
        c%cnode(3)=c%Octson%son(3)%cnode(3)
        c%cnode(4)=c%Octson%son(4)%cnode(4)
        c%cnode(5)=c%Octson%son(5)%cnode(5)
        c%cnode(6)=c%Octson%son(6)%cnode(6)
        c%cnode(7)=c%Octson%son(7)%cnode(7)
        c%cnode(8)=c%Octson%son(8)%cnode(8)
        return
        end subroutine corner
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
    type(FTTCell),pointer :: t
    integer :: i, j, k, l
    logical,INTENT(IN) :: meshonly ! If most of cells is NAN, output mesh only
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    l=0
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        t=>OctCell(i, j, k)
        call TmpStorageVar(t,meshonly,l)
    enddo
    enddo
    enddo
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "initTmpStorageVar time: ", tEnd-tStart
        contains
!----------------------------------------------------------------------
        recursive subroutine TmpStorageVar(c,ios,n)
        use ModInpInflow,only : Rgas, Gama00
        implicit none
        type(FTTCell),pointer :: c
        logical,INTENT(IN)    :: ios
        ! If most of cells is NAN, output mesh only
        integer,INTENT(INOUT) :: n
        REAL(R8):: u, v, w,p
        type(FTTCell),pointer :: cs
        integer:: is

        if(ASSOCIATED(c%Octson))then
            do is=1,c%Octson%nSon
                cs=>c%Octson%son(is)
                call TmpStorageVar(cs,ios,n)
            enddo
            return
        endif

        n=n+1
        if (ios)then
            cVariables(n,1)= c%Cross    ! Cross
            cVariables(n,2)= c%WallDistance
            cNodes(n,1)= c%cnode(1)
            cNodes(n,2)= c%cnode(2)
            cNodes(n,3)= c%cnode(3)
            cNodes(n,4)= c%cnode(4)
            cNodes(n,5)= c%cnode(5)
            cNodes(n,6)= c%cnode(6)
            cNodes(n,7)= c%cnode(7)
            cNodes(n,8)= c%cnode(8)
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
            cVariables(n,8)= c%Cross    ! Cross
            cVariables(n,9)= c%WallDistance
            cNodes(n,1)= c%cnode(1)
            cNodes(n,2)= c%cnode(2)
            cNodes(n,3)= c%cnode(3)
            cNodes(n,4)= c%cnode(4)
            cNodes(n,5)= c%cnode(5)
            cNodes(n,6)= c%cnode(6)
            cNodes(n,7)= c%cnode(7)
            cNodes(n,8)= c%cnode(8)
        endif
        endsubroutine TmpStorageVar
!----------------------------------------------------------------------
    endsubroutine initTmpStorageVar
!======================================================================
!======================================================================
!======================================================================
!----------------------------------------------------------------------
!----------------------------------------------------------------------
