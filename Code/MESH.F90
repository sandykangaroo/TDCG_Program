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
        integer :: nGeoPoints   ! Number of the geometry points.
        integer :: nGeoFaces    ! Number of the geometry faces.
        real(R8),ALLOCATABLE :: Geometry(:,:)
        integer ,ALLOCATABLE :: GeoFace(:,:)
        type(typPoint)       :: GeoBBOX(2) ! Bounding box of geometry.
    endmodule ModMesh
!======================================================================
    module ModMeshTools
    use ModMesh
    use ModTypDef
    implicit none
    contains
        integer function initCellIntersect(c)
        use ModInpMesh
        implicit none
        type(typCell),pointer :: c

        select case (cIntersectMethod)
        case(1)
            initCellIntersect=CellCast(c)
        end select
        endfunction initCellIntersect
!----------------------------------------------------------------------
        logical function BBOX(p,box)
        implicit none
        real(R8),INTENT(IN)::p(3)
        real(R8),INTENT(IN)::box(6)

        if (p(1) < box(1) .or.    &
            p(2) < box(2) .or.    &
            p(3) < box(3) .or.    &
            p(1) > box(4) .or.    &
            p(2) > box(5) .or.    &
            p(3) > box(6)) then
            BBOX = .False.
        else
            BBOX = .True.
        endif
        endfunction BBOX
!----------------------------------------------------------------------
        integer function CellCast(c)
        use ModKDTree
        use ModInpGlobal
        implicit none
        type(typCell),pointer :: c
        type(typPoint)        :: p(8)
        real(R8):: x, y, z, dx, dy, dz
        real(R8):: box(6)
        integer :: i, ii, Pintersect

        x=c%Center(1); y=c%Center(2); z=c%Center(3)
        dx=BGStep(1)/2**(c%lvl(1)+1)
        dy=BGStep(2)/2**(c%lvl(2)+1)
        dz=BGStep(3)/2**(c%lvl(3)+1)
        p(1)%P=(/x-dx,y-dy,z-dz/)
        p(2)%P=(/x+dx,y-dy,z-dz/)
        p(3)%P=(/x+dx,y+dy,z-dz/)
        p(4)%P=(/x-dx,y+dy,z-dz/)
        p(5)%P=(/x-dx,y-dy,z+dz/)
        p(6)%P=(/x+dx,y-dy,z+dz/)
        p(7)%P=(/x+dx,y+dy,z+dz/)
        p(8)%P=(/x-dx,y+dy,z+dz/)

        Pintersect=0
        do i = 1,8
            do ii = 1,nGeometry
            ! Quick BBOX identify
            ! if (BBOX(p(i)%P,KDTree(ii)%root%box)) then ! inside
                Pintersect=Pintersect+RayCast(p(i)%P,KDTree(ii)%root)
            ! endif
            enddo
        enddo

        if (Pintersect==0) then
            CellCast = 0    ! outside
        elseif (Pintersect==8) then
            CellCast = 3    ! inside
        else
            Pintersect=0
            do ii = 1,nGeometry
                Pintersect=Pintersect+RayCast(c%Center,KDTree(ii)%root)
            enddo
            if (Pintersect==0) then
                CellCast=1  ! intersect but outside
            else
                CellCast=2  ! intersect but inside
            endif
        endif
        endfunction CellCast
!----------------------------------------------------------------------
        integer function RayCast(point,tree)
        use ModGeometry
        use ModKDTree
        use ModInpGlobal, only: nGeometry
        implicit none
        real(R8),INTENT(IN):: point(3)   ! point(3) = x, y, z.
        type(KDT_node),pointer:: tree
        integer,PARAMETER  :: nRays = 3   ! Number of Rays.
        type(typPoint)     :: k(nRays)  ! Three rays
        integer :: i, j, jj, ng
        integer :: nIntersect ! Save the intersections number of any ray.
        integer :: nRayCast   ! Once a ray intersections, nRay=nRay+1.
        type(typPoint):: triFace(3)

        k(1)%P = [1., 0., 0.]
        k(2)%P = [0., 1., 0.]
        k(3)%P = [0., 0., 1.]

        nRayCast   = 0
        do ng= 1,nGeometry
        do i = 1,nRays
            nIntersect = 0
            do j = 1,body(ng)%nse
                do jj = 1,3
                    triFace(jj)%P = body(ng)%se3d(j)%P(jj)%P
                enddo
                ! Quick BBOX identify
                ! if (.not. BBOX(p(i)%P,body(ng)%box)) cycle ! outside
                nIntersect=nIntersect+MollerTrumbore(Point,k(i)%P,triFace)
            enddo
            if (nIntersect == 0) then
                RayCast=0; return
            elseif (mod(nIntersect,2)==0) then
                nRayCast = nRayCast + 0
            else
                nRayCast = nRayCast + 1
            endif
        enddo
        enddo
        
        if (nRayCast > nRays/2)then ! nRays/2=INT(real.nRays/2.0)
            RayCast=1; return
        else
            RayCast=0; return
        endif
        endfunction RayCast
!----------------------------------------------------------------------
        integer function MollerTrumbore(Point,D,tri)
        use ModTools
        implicit none
        real(R8)      ,INTENT(IN):: point(3)   ! point(3) = x, y, z.
        real(R8)      ,INTENT(IN):: D(3)
        type(typPoint),INTENT(IN):: tri(3)
        real(R8):: V0(3), V1(3), V2(3)
        real(R8):: P(3), Q(3)
        real(R8):: t, u, v
        real(R8):: A

        V0 = point-tri(1)%P
        V1 = tri(2)%P-tri(1)%P
        V2 = tri(3)%P-tri(1)%P

        ! Parallel identify
        P = CROSS_PRODUCT_3(V1,V2)
        t = DOT_PRODUCT(D,P)
        if (t==0) then  ! Parallel
            MollerTrumbore=0; return
        endif

        P = CROSS_PRODUCT_3(D,V2)
        Q = CROSS_PRODUCT_3(V0,V1)
        A = DOT_PRODUCT(P,V1)   ! A=0 is impossible.
        t = DOT_PRODUCT(Q,V2)
        u = DOT_PRODUCT(P,V0)
        v = DOT_PRODUCT(Q,D)

        if (A>0 .and. u>=0 .and. v>=0 .and. u+v<=A .and. t>0) then
            MollerTrumbore=1; return
        elseif (A<0 .and. u<=0 .and. v<=0 .and. u+v>=A .and. t<0) then
            MollerTrumbore=1; return
        else
        MollerTrumbore=0; return
        endif
        endfunction MollerTrumbore
!----------------------------------------------------------------------
        integer function RayCast2D(point,triFace)
        ! Not used. 2021.3.22
        ! Barycentric Coordinates Method (in the Theory Guide)
        real(R8)      ,INTENT(IN):: point(3)   ! point(3) = x, y, z.
        type(typPoint),INTENT(IN):: triFace(3)
        !real(R8):: A(2), B(2), C(2) ! Vertices's coordinate of triangle
        !real(R8):: P(2) ! Point got by dimensionality reduction
        real(R8):: V0(3), V1(3), V2(3)    ! Vectors AP, AB, AC
        real(R8):: c1, c2   ! Coefficient u, v in the Theory Guide.
        real(R8):: dot01, dot02, dot11, dot12, dot22 ! Dot product
        real(R8):: t        ! Temp value

        ! Quick BBOX detremine
        if     (triFace(1)%P(1) > point(1) .and.    &
                triFace(2)%P(1) > point(1) .and.    &
                triFace(3)%P(1) > point(1)) then
                RayCast2D=0; return
        elseif (triFace(1)%P(1) < point(1) .and.    &
                triFace(2)%P(1) < point(1) .and.    &
                triFace(3)%P(1) < point(1)) then
                RayCast2D=0; return
        elseif (triFace(1)%P(2) > point(2) .and.    &
                triFace(2)%P(2) > point(2) .and.    &
                triFace(3)%P(2) > point(2)) then
                RayCast2D=0; return
        elseif (triFace(1)%P(2) < point(2) .and.    &
                triFace(2)%P(2) < point(2) .and.    &
                triFace(3)%P(2) < point(2)) then
                RayCast2D=0; return
        elseif (triFace(1)%P(3) > point(3) .and.    &
                triFace(2)%P(3) > point(3) .and.    &
                triFace(3)%P(3) > point(3)) then
                RayCast2D=0; return
        elseif (triFace(1)%P(3) < point(3) .and.    &
                triFace(2)%P(3) < point(3) .and.    &
                triFace(3)%P(3) < point(3)) then
                RayCast2D=0; return
        endif
        ! Dimensionality reduction
        ! if (triFace(1)%P(3) == triFace(2)%P(3).and.   &
        !     triFace(2)%P(3) == triFace(3)%P(3))then
        !     A=[triFace(1)%P(1),triFace(1)%P(2)]
        !     B=[triFace(2)%P(1),triFace(2)%P(2)]
        !     C=[triFace(3)%P(1),triFace(3)%P(2)]
        !     P=[point(1),point(2)]
        ! else
        !     A=[triFace(1)%P(3),triFace(1)%P(2)]
        !     B=[triFace(2)%P(3),triFace(2)%P(2)]
        !     C=[triFace(3)%P(3),triFace(3)%P(2)]
        !     P=[point(3),point(2)]
        ! endif

        ! Barycentric Coordinates Method (in the Theory Guide)
        V0 = point - triFace(1)%P
        V1 = triFace(2)%P - triFace(1)%P
        V2 = triFace(3)%P - triFace(1)%P
        ! dot01 = DPROD(V0,V1)
        ! dot02 = DPROD(V0,V2)
        ! dot11 = DPROD(V1,V1)
        ! dot12 = DPROD(V1,V2)
        ! dot22 = DPROD(V2,V2)
        dot01 = DOT_PRODUCT(V0,V1)
        dot02 = DOT_PRODUCT(V0,V2)
        dot11 = DOT_PRODUCT(V1,V1)
        dot12 = DOT_PRODUCT(V1,V2)
        dot22 = DOT_PRODUCT(V2,V2)
        t=1/(dot22*dot11-dot12*dot12)

        c1=(dot11*dot02-dot12*dot01)*t
        if (c1 < 0 .or. c1 >1) then
            RayCast2D=0; return
        endif

        c2=(dot22*dot01-dot12*dot02)*t
        if (c2 < 0 .or. c2 >1) then
            RayCast2D=0; return
        endif

        if (c1+c2 <= 1) then
            RayCast2D=0; return
        else
            RayCast2D=1; return
        endif
        endfunction RayCast2D
!----------------------------------------------------------------------
        !function IntersectPoint(point,normal,triFace)
        !!Not used 2021.3.22
        !! Method (in the Theory Guide)
        !implicit none
        !real(8),DIMENSION(3):: IntersectPoint
        !real(8)       ,INTENT(IN):: point
        !real(8)       ,INTENT(IN):: normal
        !type(typPoint),INTENT(IN):: triFace(3)
        !real(R8):: n(3)       ! Normal vector: n=ABxAC
        !real(R8):: V1(3), V2(3)
        !real(R8):: d    ! Temp value
        !! Find the triangle surface Ax+By+Cz+D=0
        !!   turned to formula n(1)*x+n(2)*y+n(3)*z=c
        !V1 = triFace(2)%P - triFace(1)%P
        !V2 = triFace(3)%P - triFace(1)%P
        !d = n(1)*triFace(1)%P(1)+ &
        !    n(2)*triFace(1)%P(2)+ &
        !    n(3)*triFace(1)%P(3)
        !
        !endfunction IntersectPoint
!----------------------------------------------------------------------
        subroutine NewCell(c,split)
        ! Called by: many
        ! Calls    : initNeighbor
        implicit none
        integer(I4),INTENT(IN):: split
        type(typCell),pointer :: c
        real(R8)              :: x, y, z, dx, dy, dz

        select case (split)
        case (0)
            ALLOCATE( c%son1, c%son2, c%son3, c%son4,   &
                      c%son5, c%son6, c%son7, c%son8 )
            c%son1%nBGCell=c%nBGCell
            c%son2%nBGCell=c%nBGCell
            c%son3%nBGCell=c%nBGCell
            c%son4%nBGCell=c%nBGCell
            c%son5%nBGCell=c%nBGCell
            c%son6%nBGCell=c%nBGCell
            c%son7%nBGCell=c%nBGCell
            c%son8%nBGCell=c%nBGCell
            c%son1%nCell=c%nCell
            c%son2%nCell=nCells+1
            c%son3%nCell=nCells+2
            c%son4%nCell=nCells+3
            c%son5%nCell=nCells+4
            c%son6%nCell=nCells+5
            c%son7%nCell=nCells+6
            c%son8%nCell=nCells+7
            nCells=nCells+7
            c%son1%lvl=c%lvl+1
            c%son2%lvl=c%lvl+1
            c%son3%lvl=c%lvl+1
            c%son4%lvl=c%lvl+1
            c%son5%lvl=c%lvl+1
            c%son6%lvl=c%lvl+1
            c%son7%lvl=c%lvl+1
            c%son8%lvl=c%lvl+1
            c%son1%fSplitType=0
            c%son2%fSplitType=0
            c%son3%fSplitType=0
            c%son4%fSplitType=0
            c%son5%fSplitType=0
            c%son6%fSplitType=0
            c%son7%fSplitType=0
            c%son8%fSplitType=0
            c%son1%Location=1
            c%son2%Location=2
            c%son3%Location=3
            c%son4%Location=4
            c%son5%Location=5
            c%son6%Location=6
            c%son7%Location=7
            c%son8%Location=8
            c%son1%Node=0
            c%son2%Node=0
            c%son3%Node=0
            c%son4%Node=0
            c%son5%Node=0
            c%son6%Node=0
            c%son7%Node=0
            c%son8%Node=0
            x=c%Center(1); y=c%Center(2); z=c%Center(3)
            dx=BGStep(1)/2**(c%lvl(1)+2)
            dy=BGStep(2)/2**(c%lvl(2)+2)
            dz=BGStep(3)/2**(c%lvl(3)+2)
            c%son1%Center=(/x-dx,y-dy,z-dz/)
            c%son2%Center=(/x+dx,y-dy,z-dz/)
            c%son3%Center=(/x+dx,y+dy,z-dz/)
            c%son4%Center=(/x-dx,y+dy,z-dz/)
            c%son5%Center=(/x-dx,y-dy,z+dz/)
            c%son6%Center=(/x+dx,y-dy,z+dz/)
            c%son7%Center=(/x+dx,y+dy,z+dz/)
            c%son8%Center=(/x-dx,y+dy,z+dz/)
            c%son1%U=c%U
            c%son2%U=c%U
            c%son3%U=c%U
            c%son4%U=c%U
            c%son5%U=c%U
            c%son6%U=c%U
            c%son7%U=c%U
            c%son8%U=c%U
            if (c%cross == 1 .or. c%cross == 2) then
                c%son1%cross=initCellIntersect(c%son1)
                c%son2%cross=initCellIntersect(c%son2)
                c%son3%cross=initCellIntersect(c%son3)
                c%son4%cross=initCellIntersect(c%son4)
                c%son5%cross=initCellIntersect(c%son5)
                c%son6%cross=initCellIntersect(c%son6)
                c%son7%cross=initCellIntersect(c%son7)
                c%son8%cross=initCellIntersect(c%son8)
            else
                c%son1%cross=c%cross
                c%son2%cross=c%cross
                c%son3%cross=c%cross
                c%son4%cross=c%cross
                c%son5%cross=c%cross
                c%son6%cross=c%cross
                c%son7%cross=c%cross
                c%son8%cross=c%cross
            endif
            c%son1%Father=>c
            c%son2%Father=>c
            c%son3%Father=>c
            c%son4%Father=>c
            c%son5%Father=>c
            c%son6%Father=>c
            c%son7%Father=>c
            c%son8%Father=>c
            !TODO: Neighbor
            call NullifyCell(c%son1)
            call NullifyCell(c%son2)
            call NullifyCell(c%son3)
            call NullifyCell(c%son4)
            call NullifyCell(c%son5)
            call NullifyCell(c%son6)
            call NullifyCell(c%son7)
            call NullifyCell(c%son8)
        case (1:3)
            ALLOCATE(c%son1,c%son2)
            c%son1%nBGCell=c%nBGCell
            c%son2%nBGCell=c%nBGCell
            c%son1%nCell=c%nCell
            c%son2%nCell=nCells+1
            nCells=nCells+1
            c%son1%fSplitType=split
            c%son2%fSplitType=split
            c%son1%Location=1
            c%son2%Location=2
            c%son1%Node=0
            c%son2%Node=0
            x=c%Center(1); y=c%Center(2); z=c%Center(3)
            dx=BGStep(1)/2**(c%lvl(1)+2)
            dy=BGStep(2)/2**(c%lvl(2)+2)
            dz=BGStep(3)/2**(c%lvl(3)+2)
            select case (split)
            case(1)
                c%son1%lvl=c%lvl+(/1,0,0/)
                c%son2%lvl=c%lvl+(/1,0,0/)
                c%son1%Center=(/x-dx,y,z/)
                c%son2%Center=(/x+dx,y,z/)
            case(2)
                c%son1%lvl=c%lvl+(/0,1,0/)
                c%son2%lvl=c%lvl+(/0,1,0/)
                c%son1%Center=(/x,y-dy,z/)
                c%son2%Center=(/x,y+dy,z/)
            case(3)
                c%son1%lvl=c%lvl+(/0,0,1/)
                c%son2%lvl=c%lvl+(/0,0,1/)
                c%son1%Center=(/x,y,z-dz/)
                c%son2%Center=(/x,y,z+dz/)
            end select
            c%son1%U=c%U
            c%son2%U=c%U
            if (c%cross == 1 .or. c%cross == 2) then
                c%son1%cross=initCellIntersect(c%son1)
                c%son2%cross=initCellIntersect(c%son2)
            else
                c%son1%cross=c%cross
                c%son2%cross=c%cross
            endif
            c%son1%Father=>c
            c%son2%Father=>c
            !TODO: Neighbor
            call NullifyCell(c%son1)
            call NullifyCell(c%son2)
        case (4:6)
            ALLOCATE(c%son1,c%son2,c%son3,c%son4)
            c%son1%nBGCell=c%nBGCell
            c%son2%nBGCell=c%nBGCell
            c%son3%nBGCell=c%nBGCell
            c%son4%nBGCell=c%nBGCell
            c%son1%nCell=c%nCell
            c%son2%nCell=nCells+1
            c%son3%nCell=nCells+2
            c%son4%nCell=nCells+3
            nCells=nCells+3
            c%son1%fSplitType=split
            c%son2%fSplitType=split
            c%son3%fSplitType=split
            c%son4%fSplitType=split
            c%son1%Location=1
            c%son2%Location=2
            c%son3%Location=3
            c%son4%Location=4
            c%son1%Node=0
            c%son2%Node=0
            c%son3%Node=0
            c%son4%Node=0
            x=c%Center(1); y=c%Center(2); z=c%Center(3)
            dx=BGStep(1)/2**(c%lvl(1)+2)
            dy=BGStep(2)/2**(c%lvl(2)+2)
            dz=BGStep(3)/2**(c%lvl(3)+2)
            select case (split)
            case (4)
                c%son1%lvl=c%lvl+(/1,1,0/)
                c%son2%lvl=c%lvl+(/1,1,0/)
                c%son3%lvl=c%lvl+(/1,1,0/)
                c%son4%lvl=c%lvl+(/1,1,0/)
                c%son1%Center=(/x-dx,y-dy,z/)
                c%son2%Center=(/x+dx,y-dy,z/)
                c%son3%Center=(/x+dx,y+dy,z/)
                c%son4%Center=(/x-dx,y+dy,z/)
            case (5)
                c%son1%lvl=c%lvl+(/1,0,1/)
                c%son2%lvl=c%lvl+(/1,0,1/)
                c%son3%lvl=c%lvl+(/1,0,1/)
                c%son4%lvl=c%lvl+(/1,0,1/)
                c%son1%Center=(/x-dx,y,z-dz/)
                c%son2%Center=(/x+dx,y,z-dz/)
                c%son3%Center=(/x+dx,y,z-dz/)
                c%son4%Center=(/x-dx,y,z-dz/)
            case (6)
                c%son1%lvl=c%lvl+(/0,1,1/)
                c%son2%lvl=c%lvl+(/0,1,1/)
                c%son3%lvl=c%lvl+(/0,1,1/)
                c%son4%lvl=c%lvl+(/0,1,1/)
                c%son1%Center=(/x,y-dy,z-dz/)
                c%son2%Center=(/x,y-dy,z-dz/)
                c%son3%Center=(/x,y+dy,z-dz/)
                c%son4%Center=(/x,y+dy,z-dz/)
            end select
            c%son1%U=c%U
            c%son2%U=c%U
            c%son3%U=c%U
            c%son4%U=c%U
            if (c%cross == 1 .or. c%cross == 2) then
                c%son1%cross=initCellIntersect(c%son1)
                c%son2%cross=initCellIntersect(c%son2)
                c%son3%cross=initCellIntersect(c%son3)
                c%son4%cross=initCellIntersect(c%son4)
            else
                c%son1%cross=c%cross
                c%son2%cross=c%cross
                c%son3%cross=c%cross
                c%son4%cross=c%cross
            endif
            c%son1%Father=>c
            c%son2%Father=>c
            c%son3%Father=>c
            c%son4%Father=>c
            !TODO: Neighbor
            call NullifyCell(c%son1)
            call NullifyCell(c%son2)
            call NullifyCell(c%son3)
            call NullifyCell(c%son4)
        end select
        endsubroutine NewCell
!----------------------------------------------------------------------
        subroutine DeletCell(c)
        implicit none
        type(typCell),pointer :: c
        endsubroutine DeletCell
!----------------------------------------------------------------------
        subroutine NullifyCell(c)
        implicit none
        type(typCell),pointer :: c
        NULLIFY(c%Father,                                       &
                c%son1, c%son2, c%son3, c%son4,                 &
                c%son5, c%son6, c%son7, c%son8,                 &
                c%NeighborX1, c%NeighborX2, c%NeighborY1,       &
                c%NeighborY2, c%NeighborZ1, c%NeighborZ2)
        endsubroutine NullifyCell
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    end module ModMeshTools
!======================================================================
!======================================================================
    subroutine GenerateBGMesh   ! BG -- back-ground
    use ModPrecision
    use ModTypDef
    use ModMesh
    use ModInpMesh
    use ModMeshTools
    implicit none
    integer :: i
    type(typCell),pointer :: t

    BGStep=abs(Domain1-Domain2)/nCell
    ! BGStep(2)=Domain(2)/nCell(2)
    ! BGStep(3)=Domain(3)/nCell(3);
    nBGCells=nCell(1)*nCell(2)*nCell(3)
    nCells=nBGCells
    ALLOCATE(Cell(nBGCells))
    do i = 1, nBGCells
        t=>Cell(i)
            t%nBGCell     = i
            t%nCell       = i
            t%lvl         = 0
            t%fSplitType  = 0
            t%Location    = 0
            t%Node        = 0
            t%Center      = GBGMFindCellCenter(i)
            t%U           = 0
            t%cross       = initCellIntersect(t)
            NULLIFY(t%Father,                                       &
                    t%son1, t%son2, t%son3, t%son4,                 &
                    t%son5, t%son6, t%son7, t%son8,                 &
                    t%NeighborX1, t%NeighborX2, t%NeighborY1,       &
                    t%NeighborY2, t%NeighborZ1, t%NeighborZ2)
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
    GBGMFindCellCenter=(/Domain1(1)+(xx-0.5)*BGStep(1),             &
                         Domain1(2)+(yy-0.5)*BGStep(2),             &
                         Domain1(3)+(zz-0.5)*BGStep(3)/)
    endfunction GBGMFindCellCenter
    endsubroutine GenerateBGMesh
!======================================================================
!======================================================================
    subroutine initSurfaceAdapt
    use ModMesh
    use ModMeshTools
    use ModInpMesh
    implicit none
    type(typCell),pointer :: t
    integer :: i, ii

    do i=1,nBGCells
        t=>Cell(i)
        call SurfaceAdapt(t)
    enddo
        contains
!----------------------------------------------------------------------
        recursive subroutine SurfaceAdapt(c)
    implicit none
        type(typCell),pointer :: c

        do ii=1,3
            if (c%lvl(ii) >= InitRefineLVL) return
        enddo
        if (c%cross == 1 .or. c%cross ==2)then
            call NewCell(c,0)
            call SurfaceAdapt(c%son1)
            call SurfaceAdapt(c%son2)
            call SurfaceAdapt(c%son3)
            call SurfaceAdapt(c%son4)
            call SurfaceAdapt(c%son5)
            call SurfaceAdapt(c%son6)
            call SurfaceAdapt(c%son7)
            call SurfaceAdapt(c%son8)
        endif
    endsubroutine SurfaceAdapt
    endsubroutine initSurfaceAdapt
!======================================================================
!======================================================================
!======================================================================
!----------------------------------------------------------------------
