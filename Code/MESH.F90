!======================================================================
    module ModMesh
        use ModTypDef
        implicit none
        ! Mesh
        integer :: nBGCells     ! Number of the background Cells.
        integer :: nCells=0     ! Number of total Cells.
        real(R8):: BGStep(3)    ! Step size for background cell.
        type(typCell),pointer :: Cell(:)
        ! Geometry
        integer :: nPoints
        real(R8),ALLOCATABLE :: Geometry(:,:)
    endmodule ModMesh
!======================================================================
    module ModCellInsect
    use ModMesh
    use ModTypDef
    implicit none
    contains
        integer function initCellIntersect(c)
        use ModInpMesh
        type(typCell),pointer :: c

        select case (cIntersectMethod)
        case(1)
            initCellIntersect=RayCast(c)
        end select
        endfunction initCellIntersect
!----------------------------------------------------------------------
        integer function RayCast(c)
        type(typCell),pointer :: c
        endfunction RayCast
!----------------------------------------------------------------------
    end module ModCellInsect
!======================================================================
    subroutine ReadGeometry
    use ModInpGlobal
    use ModMesh
    implicit none
    integer :: ios, i, j 
    character(10)::FileForm

    print*,'Read geometry file: ', GeometryName
    open(unit=31, file=GeometryName, iostat=ios, status="old", action="read")
    if ( ios /= 0 ) stop ("Error opening file: "//GeometryName)
    read(31, fmt="(G10.1)", iostat=ios) FileForm
    if ( ios /= 0 ) stop ("Error reading file: "//GeometryName)
    if (.NOT.FileForm=="FACET FILE")                                &
        stop ("Error file header: "//GeometryName)
    read(31,"(///I9)") nPoints
    print*,'Geometry points count: ', nPoints
    ALLOCATE(Geometry(nPoints,3))
    read(31,"(3(E23.15,1X))") ((Geometry(i,j),j=1,3),i=1,nPoints)
    close(31)
    endsubroutine ReadGeometry
    subroutine GenerateBGMesh   ! BG -- back-ground
    use ModPrecision
    use ModTypDef
    use ModMesh
    use ModInpMesh
    use ModCellInsect
    implicit none
    integer :: i
    type(typCell),pointer :: t

    BGStep(1)=Domain(1)/nCell(1)
    BGStep(2)=Domain(2)/nCell(2)
    BGStep(3)=Domain(3)/nCell(3);
    nBGCells=nCell(1)*nCell(2)*nCell(3)
    ALLOCATE(Cell(nBGCells))
    do i = 1, nBGCells
        t=>Cell(i)
            t%nBGCell     = i
            t%nCell       = i
            t%lvl         = 0
            t%cross       = initCellIntersect(t)
            t%fSplitType  = 0
            t%Location    = 0
            t%Node        = 0
            t%Center      = GBGMFindCellCenter(i)
            t%U           = 0
            NULLIFY(t%Father,                                       &
                    t%son1, t%son2, t%son3, t%son4,                 &
                    t%son5, t%son6, t%son7, t%son8,                 &
                    t%NeighborX1, t%NeighborX2, t%NeighborY1,       &
                    t%NeighborY2, t%NeighborZ1, t%NeighborZ2)
        nCells=nCells+1
    enddo
    contains
!----------------------------------------------------------------------
    function GBGMFindCellCenter(num)
    use ModPrecision
    implicit none
    real(R8),DIMENSION(3)   :: GBGMFindCellCenter
    integer,INTENT(IN)      :: num
    integer                 :: xx, yy, zz

    xx=mod(num,nCell(1))
    if (xx==0) xx=nCell(1)
    yy=mod(int((num-xx)/nCell(1)+1),nCell(2))
    if (yy==0) yy=nCell(2)
    zz=ceiling(real(num)/real(nCell(1)*nCell(2)))
    GBGMFindCellCenter=(/(xx-0.5)*BGStep(1),                    &
                         (yy-0.5)*BGStep(2),                    &
                         (zz-0.5)*BGStep(3)/)
    endfunction GBGMFindCellCenter
    endsubroutine GenerateBGMesh
!======================================================================
!======================================================================
    subroutine SurfaceAdapt
    use ModMesh
    implicit none

    endsubroutine SurfaceAdapt
!======================================================================
!======================================================================
!======================================================================
!----------------------------------------------------------------------
