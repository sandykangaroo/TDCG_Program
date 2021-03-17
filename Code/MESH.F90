!======================================================================
    module ModMesh
        use ModTypDef
        implicit none
        integer :: nBGCells     ! Number of the background Cells.
        integer :: nCells=0     ! Number of total Cells.
        real(R8):: BGStep(3)    ! Step size for background cell.
        type(typCell),pointer :: Cell(:)
    endmodule ModMesh
!======================================================================
    subroutine GenerateBGMesh   ! BG -- back-ground
    use ModPrecision
    use ModTypDef
    use ModMesh
    use ModInpMesh
    implicit none
    integer :: i

    BGStep(1)=Domain(1)/nCell(1)
    BGStep(2)=Domain(2)/nCell(2)
    BGStep(3)=Domain(3)/nCell(3);
    nBGCells=nCell(1)*nCell(2)*nCell(3)
    ALLOCATE(Cell(nBGCells))
    do i = 1, nBGCells
        ASSOCIATE(t=>Cell(i))
            t%nBGCell     = i
            t%nCell       = i
            t%lvl         = 0
            t%cross       = 0 ! CellInsect(Cell(i))
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
        end ASSOCIATE
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
    function CellInsect(c)
    use ModPrecision
    use ModTypDef
    implicit none
    integer :: CellInsect
    type(typCell),pointer :: c
    end function CellInsect
!======================================================================
!======================================================================
!======================================================================
!======================================================================
!======================================================================
!----------------------------------------------------------------------
