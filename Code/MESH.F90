!======================================================================
    module ModMeshTools
    use ModMesh
    use ModTypDef
    implicit none
    contains
!----------------------------------------------------------------------
        subroutine initCellCross(c)
        use ModInpMesh
        use ModKDTree
        use ModInpGlobal, only: nGeometry
        implicit none
        type(FTTCell),pointer :: c

        select case (cIntersectMethod)
        case (1:2)
            call CellCast(c)
        case (3:4)
            call AABB(c)
        ! case (5)
        !     call AABBTraverse(c)
        ! case (6)
        !     call CellCast(c)
        case (7) ! After surface refine, change cIntersectMethod=1,2 to =7
            call CellCast(c)
            if (c%Cross == -4) call CellInout(c)
        case (8) ! After surface refine, change cIntersectMethod=3,4 to =8
            call AABB(c)
            if (c%Cross == -4) call CellInout(c)
        end select
        endsubroutine initCellCross
!----------------------------------------------------------------------
        subroutine CellCast(c)
        use ModKDTree
        use ModInpGlobal
        use ModInpMesh,only: cIntersectMethod
        implicit none
        type(FTTCell),pointer :: c
        type(typPoint)        :: p(8)
        real(R8):: x, y, z, dx, dy, dz
        integer :: i, ii, iii
        integer :: Pintersect ! Number of the intersect point
        integer :: nIntersect ! Intersect times for one point
        integer :: Rintersect ! Once ray intersects, nRay++
        real(R8):: k(3)

        x=c%Center(1); y=c%Center(2); z=c%Center(3)
        dx=BGCellSize(1)/2**(c%LVL(1)+1)
        dy=BGCellSize(2)/2**(c%LVL(2)+1)
        dz=BGCellSize(3)/2**(c%LVL(3)+1)
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
            loop2:do ii = 1,nGeometry
                Rintersect=0
                loop3:do iii = 1, 3  !nRays
                    nIntersect=0
                    ! Rays direction:
                    ! k(1)%P = [1., 0., 0.]
                    ! k(2)%P = [0., 1., 0.]
                    ! k(3)%P = [0., 0., 1.]
                    k = CSHIFT([0.,0.,1.],4-iii)
                    call RayCast(CSHIFT(p(i)%P,iii-1), &
                                 KDTree(ii)%root,k,nIntersect,iii)
                    ! If no intersect, must be a outside point, so return.
                    if (nIntersect == 0) then
                        Pintersect=Pintersect+0
                        exit loop2
                    elseif (mod(nIntersect,2)==0) then
                        Rintersect = Rintersect + 0
                    else
                        Rintersect = Rintersect + 1
                    endif
                    ! If twice back a same result, not do the third.
                    if ( iii==2 .and. Rintersect/=1 ) exit loop3
                enddo loop3
                if ( Rintersect > 1 ) Pintersect=Pintersect+1 ! inside
            enddo loop2
            if (Pintersect/=0.and.Pintersect<i) then ! Intersect
                c%Cross = -3
                return
            endif
        enddo
        ! Not intersect
        c%Cross = -4
        endsubroutine CellCast
!----------------------------------------------------------------------
        subroutine CellInout(c)
        use ModKDTree
        use ModInpGlobal
        use ModInpMesh,only: cIntersectMethod
        use ModTools, only: BBOX
        implicit none
        type(FTTCell),pointer :: c
        integer :: nIntersect ! Intersect times for one point
        integer :: Rintersect ! Once ray intersects, Rintersect++
        integer :: iii ! Ray
        real(R8):: k(3)

        if (c%Cross == -4) then
            c%Cross = 0
        elseif (c%Cross == -3) then
            c%Cross = 1
        else
            stop 'Subroutine CellInout error1'
        endif
        ! Quick Bounding-BOX identify
        if (BBOX(c%Center,KDTree(1)%root%box)>=0) then
            Rintersect=0
            do iii = 1, 3  !nRays
                nIntersect=0
                ! Rays direction:
                ! k(1)%P = [1., 0., 0.]
                ! k(2)%P = [0., 1., 0.]
                ! k(3)%P = [0., 0., 1.]
                k = CSHIFT([0.,0.,1.],4-iii)
                call RayCast(CSHIFT(c%Center,iii-1), &
                             KDTree(1)%root,k,nIntersect,iii)
                ! If no intersect, must be a outside point, so return.
                if (nIntersect == 0) then
                    return
                elseif (mod(nIntersect,2)==0) then
                    Rintersect = Rintersect + 0
                else
                    Rintersect = Rintersect + 1
                endif
                ! If twice back a same result, not do the third.
                if ( iii==2 .and. Rintersect/=1 ) exit
            enddo
            if ( Rintersect > 1 ) then ! inside
                if (c%Cross == 0) then
                    c%Cross = 3
                elseif (c%Cross == 1) then
                    c%Cross = 2
                endif
            endif
        endif
        endsubroutine CellInout
!----------------------------------------------------------------------
        recursive subroutine RayCast(point,tree,k,nIntersect,i)
        ! Back nIntersect -> times of intersect
        use ModGeometry
        use ModKDTree
        use ModTools, only: MollerTrumbore
        implicit none
        real(R8),INTENT(IN)    :: point(3)     ! point(3) = x, y, z.
        type(KDT_node),pointer :: tree
        real(R8),INTENT(IN)    :: k(3)
        integer, INTENT(INOUT) :: nIntersect ! Save the intersect number
        integer, INTENT(IN)    :: i ! Which Ray
        real(R8)               :: box(6)

        if (.not.ASSOCIATED(tree)) return

        ! Quick Bounding-BOX identify
        select case (i) ! faster than cshift.
        case (1)
            box(1) = tree%box(1)
            box(2) = tree%box(2)
            box(4) = tree%box(4)
            box(5) = tree%box(5)
            box(6) = tree%box(6)
        case (2)
            box(1) = tree%box(2)
            box(2) = tree%box(3)
            box(4) = tree%box(5)
            box(5) = tree%box(6)
            box(6) = tree%box(4)
        case (3)
            box(1) = tree%box(3)
            box(2) = tree%box(1)
            box(4) = tree%box(6)
            box(5) = tree%box(4)
            box(6) = tree%box(5)
        end select
        if (point(1)>box(1).and.point(2)>box(2).and. &
            point(1)<box(4).and.point(2)<box(5).and. & ! Box max
            point(3)<box(6)) then ! Positive direction
            if (MollerTrumbore(CSHIFT(point,1-i),k, &
                               tree%the_data,.false.)) &
                nIntersect = nIntersect + 1
                ! Into the next tree
            call RayCast(point,tree%right,k,nIntersect,i)
            call RayCast(point,tree%left,k,nIntersect,i)
        endif
        ! box = CSHIFT(tree%box(1:3),i-1)
        ! if (point(1)>box(1).and.point(2)>box(2)) then ! Box min
        !     box = CSHIFT(tree%box(4:6),i-1)
        !     if(point(1)<box(1).and.point(2)<box(2).and. & ! Box max
        !        point(3)<box(3)) then ! Positive direction
        !         if (MollerTrumbore(CSHIFT(point,1-i),&
                        !k,tree%the_data,.false.)) then
        !             nIntersect = nIntersect + 1
        !         endif
        !         ! Into the next tree
        !         call RayCast(point,tree%right,k,nIntersect,i)
        !         call RayCast(point,tree%left,k,nIntersect,i)
        !     endif
        ! endif

        endsubroutine RayCast
!----------------------------------------------------------------------
        ! recursive subroutine RayCastTraverse(point,tree,k,nIntersect,i)
        ! use ModGeometry
        ! use ModKDTree
        ! implicit none
        ! real(R8),INTENT(IN)    :: point(3)     ! point(3) = x, y, z.
        ! type(KDT_node),pointer,INTENT(IN) :: tree
        ! real(R8),INTENT(IN)    :: k(3)
        ! integer, INTENT(INOUT) :: nIntersect ! Save the intersect number
        ! integer, INTENT(IN)    :: i ! Which Ray
        ! real(R8)               :: box(6)

        ! if (.not.ASSOCIATED(tree)) return

        ! if (MollerTrumbore(CSHIFT(point,1-i),k,tree%the_data,.false.)) then
        !     nIntersect = nIntersect + 1
        ! endif
        ! ! Into the next tree
        ! call RayCastTraverse(point,tree%right,k,nIntersect,i)
        ! call RayCastTraverse(point,tree%left,k,nIntersect,i)
        ! endsubroutine RayCastTraverse
!----------------------------------------------------------------------
        subroutine AABB(c)
        use ModKDTree
        use ModInpGlobal
        use ModGeometry
        implicit none
        type(FTTCell),pointer :: c
        type(triangle)           :: tri
        real(R8)                 :: boxCell(6)
        integer                  :: aaa

        aaa=0
        boxCell(1)=c%center(1)-BGCellSize(1)/2.**(c%LVL(1)+1)
        boxCell(2)=c%center(2)-BGCellSize(2)/2.**(c%LVL(2)+1)
        boxCell(3)=c%center(3)-BGCellSize(3)/2.**(c%LVL(3)+1)
        boxCell(4)=c%center(1)+BGCellSize(1)/2.**(c%LVL(1)+1)
        boxCell(5)=c%center(2)+BGCellSize(2)/2.**(c%LVL(2)+1)
        boxCell(6)=c%center(3)+BGCellSize(3)/2.**(c%LVL(3)+1)

        call AABBFindTri(c,boxCell,c%Father%CrossTri,aaa)
        ! c%nCrossTri = aaa
        if (aaa==0) then
            c%Cross = -4
        else
            c%Cross = -3
        endif
        return
        end subroutine AABB
!----------------------------------------------------------------------
        ! subroutine AABBTraverse(c)!1=intersect,0=no intersect
        ! use ModKDTree
        ! use ModInpGlobal
        ! use ModGeometry
        ! implicit none
        ! type(FTTCell),pointer :: c
        ! type(triangle)        :: tri
        ! integer                     :: ng, i
        ! real(R8)              ::boxCell(6)

        ! c%Cross=-4
        ! Loop1:do ng=1,nGeometry
        !     do i=1, body(ng)%nse
        !         tri= body(ng)%se3d(i)
        !         if(TriBoxOverlap (c,tri))then
        !             boxCell(1)=c%center(1)-BGCellSize(1)/2**(c%LVL(1)+1)
        !             boxCell(2)=c%center(2)-BGCellSize(2)/2**(c%LVL(2)+1)
        !             boxCell(3)=c%center(3)-BGCellSize(3)/2**(c%LVL(3)+1)
        !             boxCell(4)=c%center(1)+BGCellSize(1)/2**(c%LVL(1)+1)
        !             boxCell(5)=c%center(2)+BGCellSize(2)/2**(c%LVL(2)+1)
        !             boxCell(6)=c%center(3)+BGCellSize(3)/2**(c%LVL(3)+1)
        !             exit Loop1
        !         endif
        !     enddo
        ! enddo Loop1
        ! return
        ! end subroutine AABBTraverse
!----------------------------------------------------------------------
        recursive subroutine AABBFindTri(c,boxCell,node,aaa)
        use ModKDTree
        implicit none
        type(FTTCell),pointer    :: c
        real(R8),INTENT(IN)      :: boxCell(6)
        type(tCrossTri), pointer :: node
        integer, INTENT(INOUT)   :: aaa
        integer                  :: i
        type(tCrossTri), pointer :: CrossTri

        if (.not. ASSOCIATED(node)) return

        if (boxCell(1)<=node%tri%box(4) .and. &
            boxCell(4)>=node%tri%box(1) .and. &
            boxCell(2)<=node%tri%box(5) .and. &
            boxCell(5)>=node%tri%box(2) .and. &
            boxCell(3)<=node%tri%box(6) .and. &
            boxCell(6)>=node%tri%box(3) ) then
            if (TriBoxOverlap(c,node%tri%the_data)) then
                aaa = aaa + 1
                CrossTri  => c%CrossTri
                if (aaa==1) then
                    ALLOCATE(c%CrossTri)
                    c%CrossTri%tri => node%tri
                    ! CrossTri%prev=> null()
                    c%CrossTri%next=> null()
                else
                    do i = 1, aaa-2
                        CrossTri => CrossTri%next
                    enddo
                    ALLOCATE(CrossTri%next)
                    CrossTri%next%tri => node%tri
                    ! CrossTri%next%prev=> CrossTri
                    CrossTri%next%next=> null()
                endif
            endif
        endif

        call AABBFindTri(c,boxCell,node%next,aaa)

        end subroutine AABBFindTri
!----------------------------------------------------------------------
    ! use separating axis theorem to test overlap between triangle and box
    ! need to test for overlap in these directions:
    ! 1) the {x,y,z}-directions (actually, since we use the AABB of
    !    the triangle we do not even need to test these)
    ! 2) normal of the triangle
    ! 3) crossproduct(edge from tri, {x,y,z}-directin)
    ! this gives 3x3=9 more tests
        Logical function TriBoxOverlap (c,tri)
            use ModTools
            use ModKDTree
            use ModInpGlobal
            use ModGeometry
            implicit none

            type(FTTCell),pointer :: c
            type(triangle)        :: tri
            real(R8)              :: v0(3), v1(3), v2(3)
            real(R8)              :: normal(3)
            real(R8)              :: e0(3), e1(3), e2(3)
            real(R8)              :: boxhalfsize(3)
            real(R8)              :: minvlue, maxvlue
            real(R8)              :: rad, d, dd, fex, fey, fez, a
            real(R8)              :: X01P0, X01P2
            real(R8)              :: Y02P0, Y02P2
            real(R8)              :: Z12P2,Z12P1
            real(R8)              :: Z0P0, Z0P1
            real(R8)              :: X2P0, X2P1
            real(R8)              :: Y1P0, Y1P1

            boxhalfsize(1)= BGCellSize(1)/2**(c%LVL(1)+1)
            boxhalfsize(2)= BGCellSize(2)/2**(c%LVL(2)+1)
            boxhalfsize(3)= BGCellSize(3)/2**(c%LVL(3)+1)

            v0(:)=tri%p(1)%P-c%center
            v1(:)=tri%p(2)%P-c%center
            v2(:)=tri%p(3)%P-c%center

            !compute triangle edges
            e0(:)=v1-v0
            e1(:)=v2-v1
            e2(:)=v0-v2

            !/* Bullet 3: */
            !/* test the 9 tests first (this was faster) */
            fex=abs(e0(1))
            fey=abs(e0(2))
            fez=abs(e0(3))
            !AxisTestX01
            X01P0=e0(3)*v0(2)-e0(2)*v0(3)
            X01P2=e0(3)*v2(2)-e0(2)*v2(3)
            if(X01P0<X01P2)then
                minvlue=X01P0
                maxvlue=X01P2
            else
                minvlue=X01P2
                maxvlue=X01P0
            endif
            rad=fez*boxhalfsize(2)+fey*boxhalfsize(3)
            if(minvlue>rad.or.maxvlue<-rad)then
                TriBoxOverlap=.false.
                return
            endif
            !AxisTestY02
            Y02P0=-e0(3)*v0(1)+e0(1)*v0(3)
            Y02P2=-e0(3)*v2(1)+e0(1)*v2(3)
            if(Y02P0<Y02P2)then
                minvlue=Y02P0
                maxvlue=Y02P2
            else
                minvlue=Y02P2
                maxvlue=Y02P0
            endif
            rad=fez*boxhalfsize(1)+fex*boxhalfsize(3)
            if(minvlue>rad.or.maxvlue<-rad)then
                TriBoxOverlap=.false.
                return
            endif
            !AxisTestZ12
            Z12P1=e0(2)*v1(1)-e0(1)*v1(2)
            Z12P2=e0(2)*v2(1)-e0(1)*v2(2)
            if(Z12P2<Z12P1)then
                minvlue=Z12P2
                maxvlue=Z12P1
            else
                minvlue=Z12P1
                maxvlue=Z12P2
            endif
            rad=fey*boxhalfsize(1)+fex*boxhalfsize(2)
            if(minvlue>rad.or.maxvlue<-rad)then
                TriBoxOverlap=.false.
                return
            endif

            fex=abs(e1(1))
            fey=abs(e1(2))
            fez=abs(e1(3))
            !AxisTestX01
            X01P0=e1(3)*v0(2)-e1(2)*v0(3)
            X01P2=e1(3)*v2(2)-e1(2)*v2(3)
            if(X01P0<X01P2)then
                minvlue=X01P0
                maxvlue=X01P2
            else
                minvlue=X01P2
                maxvlue=X01P0
            endif
            rad=fez*boxhalfsize(2)+fey*boxhalfsize(3)
            if(minvlue>rad.or.maxvlue<-rad)then
                TriBoxOverlap=.false.
                return
            endif
            !AxisTestY02
            Y02P0=-e1(3)*v0(1)+e1(1)*v0(3)
            Y02P2=-e1(3)*v2(1)+e1(1)*v2(3)
            if(Y02P0<Y02P2)then
                minvlue=Y02P0
                maxvlue=Y02P2
            else
                minvlue=Y02P2
                maxvlue=Y02P0
            endif
            rad=fez*boxhalfsize(1)+fex*boxhalfsize(3)
            if(minvlue>rad.or.maxvlue<-rad)then
                TriBoxOverlap=.false.
                return
            endif
            !AxisTestZ0
            Z0P0=e1(2)*v0(1)-e1(1)*v0(2)
            Z0P1=e1(2)*v1(1)-e1(1)*v1(2)
            if(Z0P0<Z0P1)then
                minvlue=Z0P0
                maxvlue=Z0P1
            else
                minvlue=Z0P1
                maxvlue=Z0P0
            endif
            rad=fey*boxhalfsize(1)+fex*boxhalfsize(2)
            if(minvlue>rad.or.maxvlue<-rad)then
                TriBoxOverlap=.false.
                return
            endif


            fex=abs(e2(1))
            fey=abs(e2(2))
            fez=abs(e2(3))
            !AxisTestX2
            X2P0=e2(3)*v0(2)-e2(2)*v0(3)
            X2P1=e2(3)*v1(2)-e2(2)*v1(3)
            if(X2P0<X2P1)then
                minvlue=X2P0
                maxvlue=X2P1
            else
                minvlue=X2P1
                maxvlue=X2P0
            endif
            rad=fez*boxhalfsize(2)+fey*boxhalfsize(3)
            if(minvlue>rad.or.maxvlue<-rad)then
                TriBoxOverlap=.false.
                return
            endif

            !AxisTestY1
            Y1P0=-e2(3)*v0(1)+e2(1)*v0(3)
            Y1P1=-e2(3)*v1(1)+e2(1)*v1(3)
            if(Y1P0<Y1P1)then
                minvlue=Y1P0
                maxvlue=Y1P1
            else
                minvlue=Y1P1
                maxvlue=Y1P0
            endif
            rad=fez*boxhalfsize(1)+fex*boxhalfsize(3)
            if(minvlue>rad.or.maxvlue<-rad)then
                TriBoxOverlap=.false.
                return
            endif
            !AxisTestZ12
            Z12P1=e2(2)*v1(1)-e2(1)*v1(2)
            Z12P2=e2(2)*v2(1)-e2(1)*v2(2)
            if(Z12P2<Z12P1)then
                minvlue=Z12P2
                maxvlue=Z12P1
            else
                minvlue=Z12P1
                maxvlue=Z12P2
            endif
            rad=fey*boxhalfsize(1)+fex*boxhalfsize(2)
            if(minvlue>rad.or.maxvlue<-rad)then
                TriBoxOverlap=.false.
                return
            endif
            !/* Bullet 1: */
            !/* first test overlap in the {x,y,z}-directions */
            !/* find min, max of the triangle each direction, and test for */
            !/* overlap in that direction -- this is equivalent to testing */
            !/* a minimal AABB around the triangle against the AABB */
            !/* test in X-direction */
            minvlue=Min(v0(1),v1(1),v2(1))
            maxvlue=Max(v0(1),v1(1),v2(1))
            if(minvlue>boxhalfsize(1).or.maxvlue<-boxhalfsize(1))then
                TriBoxOverlap=.false.
                return
            endif

            !/* test in Y-direction */
            minvlue=Min(v0(2),v1(2),v2(2))
            maxvlue=Max(v0(2),v1(2),v2(2))
            if(minvlue>boxhalfsize(2).or.maxvlue<-boxhalfsize(2))then
                TriBoxOverlap=.false.
                return
            endif

            !/* test in Z-direction */
            minvlue=Min(v0(3),v1(3),v2(3))
            maxvlue=Max(v0(3),v1(3),v2(3))
            if(minvlue>boxhalfsize(3).or.maxvlue<-boxhalfsize(3)) then
                TriBoxOverlap=.false.
                return
            endif

            !/* Bullet 2: */
            !/* test if the box intersects the plane of the triangle */
            !/* compute plane equation of triangle: normal*x+d=0 */
            normal = CROSS_PRODUCT_3(e0,e1)
            if(.not.planeBoxOverlap(normal,v0,boxhalfsize)) then
                TriBoxOverlap=.false.
                return
            endif

            TriBoxOverlap=.true.
            return
        end function TriBoxOverlap
!----------------------------------------------------------------------
        Logical function planeBoxOverlap(normal,vert,maxbox)
            use ModPrecision
            implicit none

            integer  :: i,ii,iii
            real(R8) :: v,vmin(3),vmax(3),normal(3),vert(3),maxbox(3)

            do i=1,3
                v=vert(i)
                if(normal(i)>0)then
                    vmin(i)=-maxbox(i)-v
                    vmax(i)=maxbox(i)-v
                else
                    vmin(i)=maxbox(i)-v
                    vmax(i)=-maxbox(i)-v
                endif
            enddo

            ii = DOT_PRODUCT(normal, vmin)
            if(ii>0)then
                planeBoxOverlap =.false.
                return
            endif
            iii= DOT_PRODUCT(normal, vmax)
            if(iii>=0)then
                planeBoxOverlap =.true.
                return
            endif

            planeBoxOverlap =.false.
            return
        end function
!----------------------------------------------------------------------
        subroutine NewCell(c,split)
        ! Input parameter: split: fSplitTyp in module ModTypDef
        use ModInpGlobal
        use ModNeighbor
        implicit none
        type(FTTCell),pointer:: c
        integer(I4),INTENT(IN)  :: split
        integer(I4)             :: spl
        real(R8)                :: x, y, z, dx, dy, dz
        type(FTTCell),pointer   :: cs
        integer                 :: is

        ALLOCATE(c%Octson)
        if (Aniso) then
            spl=split
        else
            spl=0
        endif
        x=c%Center(1)
        y=c%Center(2)
        z=c%Center(3)
        dx=BGCellSize(1)/2**(c%LVL(1)+2)
        dy=BGCellSize(2)/2**(c%LVL(2)+2)
        dz=BGCellSize(3)/2**(c%LVL(3)+2)
        select case (spl)
        case (0)
            c%Octson%nSon = 8
            c%Octson%son(1)%Center(1:3)=(/x-dx,y-dy,z-dz/)
            c%Octson%son(2)%Center(1:3)=(/x+dx,y-dy,z-dz/)
            c%Octson%son(3)%Center(1:3)=(/x+dx,y+dy,z-dz/)
            c%Octson%son(4)%Center(1:3)=(/x-dx,y+dy,z-dz/)
            c%Octson%son(5)%Center(1:3)=(/x-dx,y-dy,z+dz/)
            c%Octson%son(6)%Center(1:3)=(/x+dx,y-dy,z+dz/)
            c%Octson%son(7)%Center(1:3)=(/x+dx,y+dy,z+dz/)
            c%Octson%son(8)%Center(1:3)=(/x-dx,y+dy,z+dz/)
        case(1)
            c%Octson%nSon = 2
            c%Octson%son(1)%Center(1:3)=(/x-dx,y,z/)
            c%Octson%son(2)%Center(1:3)=(/x+dx,y,z/)
        case(2)
            c%Octson%nSon = 2
            c%Octson%son(1)%Center(1:3)=(/x,y-dy,z/)
            c%Octson%son(2)%Center(1:3)=(/x,y+dy,z/)
        case(3)
            c%Octson%nSon = 2
            c%Octson%son(1)%Center(1:3)=(/x,y,z-dz/)
            c%Octson%son(2)%Center(1:3)=(/x,y,z+dz/)
        case (4)
            c%Octson%nSon = 4
            c%Octson%son(1)%Center(1:3)=(/x-dx,y-dy,z/)
            c%Octson%son(2)%Center(1:3)=(/x+dx,y-dy,z/)
            c%Octson%son(3)%Center(1:3)=(/x+dx,y+dy,z/)
            c%Octson%son(4)%Center(1:3)=(/x-dx,y+dy,z/)
        case (5)
            c%Octson%nSon = 4
            c%Octson%son(1)%Center(1:3)=(/x-dx,y,z-dz/)
            c%Octson%son(2)%Center(1:3)=(/x+dx,y,z-dz/)
            c%Octson%son(3)%Center(1:3)=(/x+dx,y,z-dz/)
            c%Octson%son(4)%Center(1:3)=(/x-dx,y,z-dz/)
        case (6)
            c%Octson%nSon = 4
            c%Octson%son(1)%Center(1:3)=(/x,y-dy,z-dz/)
            c%Octson%son(2)%Center(1:3)=(/x,y-dy,z-dz/)
            c%Octson%son(3)%Center(1:3)=(/x,y+dy,z-dz/)
            c%Octson%son(4)%Center(1:3)=(/x,y+dy,z-dz/)
        end select

            c%Octson%son(1)%nCell  = c%nCell
        do is = 2, c%Octson%nSon
            nCells          = nCells+1
            c%Octson%son(is)%nCell = nCells
        enddo

        do is = 1, c%Octson%nSon
            cs => c%Octson%son(is)
            cs%Father       =>c
            cs%nBGCell(1:3) = c%nBGCell(1:3)
            cs%Location     = is
            cs%LVL          = c%LVL+1
            cs%fSplitTyp    = spl
            cs%Mark(1:3)    = .false.
            cs%Cross        = c%Cross
            NULLIFY(cs%Octson)
            NULLIFY(cs%CrossTri)
        enddo
        if (c%Cross==1 .or. c%Cross==2 .or. c%Cross==-3) then
            do is = 1, c%Octson%nSon
                cs=>c%Octson%son(is)
                call initCellCross(cs)
            enddo
        endif
        return
        endsubroutine NewCell
!----------------------------------------------------------------------
        ! recursive subroutine FindNeighbor(c)
        ! use ModNeighbor
        ! implicit none
        ! type(FTTCell),pointer :: c
        ! type(FTTCell),pointer :: cs
        ! integer:: is

        ! c%Octson%NeighborX1=>NeighborX1(c)
        ! c%Octson%NeighborX2=>NeighborX2(c)
        ! c%Octson%NeighborY1=>NeighborY1(c)
        ! c%Octson%NeighborY2=>NeighborY2(c)
        ! c%Octson%NeighborZ1=>NeighborZ1(c)
        ! c%Octson%NeighborZ2=>NeighborZ2(c)

        ! if(ASSOCIATED(c%Octson))then
        !     do is=1,c%Octson%nSon
        !         cs=>c%Octson%son(is)
        !         call FindNeighbor(cs)
        !     enddo
        ! endif
        ! endsubroutine FindNeighbor
!----------------------------------------------------------------------
        recursive subroutine UpdateNeighbor(c,dirct)
        USE ModNeighbor, only: Neighbor
        implicit none
        type(FTTCell),pointer :: c
        integer,INTENT(IN)    ::dirct
        type(FTTCell),pointer :: cs
        integer:: is

        if (.not.ASSOCIATED(c%Octson)) return

        select case (dirct)
        case (0)
            c%Octson%Neighbor1=>Neighbor(c,1)
            c%Octson%Neighbor2=>Neighbor(c,2)
            c%Octson%Neighbor3=>Neighbor(c,3)
            c%Octson%Neighbor4=>Neighbor(c,4)
            c%Octson%Neighbor5=>Neighbor(c,5)
            c%Octson%Neighbor6=>Neighbor(c,6)
        case (1)
            c%Octson%Neighbor1=>Neighbor(c,dirct)
        case (2)
            c%Octson%Neighbor2=>Neighbor(c,dirct)
        case (3)
            c%Octson%Neighbor3=>Neighbor(c,dirct)
        case (4)
            c%Octson%Neighbor4=>Neighbor(c,dirct)
        case (5)
            c%Octson%Neighbor5=>Neighbor(c,dirct)
        case (6)
            c%Octson%Neighbor6=>Neighbor(c,dirct)
        end select

        do is=1,c%Octson%nSon
            cs=>c%Octson%son(is)
            call UpdateNeighbor(cs,dirct)
        enddo

        endsubroutine UpdateNeighbor
!----------------------------------------------------------------------
        subroutine DeletCell(c,split)
        implicit none
        type(FTTCell),pointer :: c
        integer(I4),INTENT(IN):: split
        endsubroutine DeletCell
!----------------------------------------------------------------------
        subroutine BGCellAABB(c)
        use ModKDTree
        use ModInpGlobal
        use ModGeometry
        implicit none
        type(FTTCell),pointer    :: c
        type(triangle)              :: tri
        real(R8)                    :: boxCell(6)
        integer                     :: aaa
        integer                     :: ng
        type(tCrossTri), pointer  :: CrossTri
        INTEGER                     :: ii
        type(FTTCell), pointer   :: cc

        aaa      =  0
        CrossTri => null()

        boxCell(1)=c%center(1)-BGCellSize(1)/2.
        boxCell(2)=c%center(2)-BGCellSize(2)/2.
        boxCell(3)=c%center(3)-BGCellSize(3)/2.
        boxCell(4)=c%center(1)+BGCellSize(1)/2.
        boxCell(5)=c%center(2)+BGCellSize(2)/2.
        boxCell(6)=c%center(3)+BGCellSize(3)/2.

        do ii = 1, nGeometry
            call BGCellFindTri(c,boxCell,kdtree(ii)%root,aaa)
        enddo
        ! c%nCrossTri = aaa
        if (aaa==0) then
            c%Cross = -4
        else
            c%Cross = -3
        endif

        end subroutine BGCellAABB
!----------------------------------------------------------------------
        recursive subroutine BGCellFindTri(c,boxCell,node,aaa)
        use ModKDTree
        implicit none
        type(FTTCell),pointer    :: c
        real(R8),INTENT(IN)         :: boxCell(6)
        type(KDT_node), pointer     :: node
        integer, INTENT(INOUT)      :: aaa
        integer                     :: ii
        type(tCrossTri), pointer  :: CrossTri

        if (.not. ASSOCIATED(node)) return

        if (boxCell(1)>node%box(4).or.boxCell(4)<node%box(1).or.&
            boxCell(2)>node%box(5).or.boxCell(5)<node%box(2).or.&
            boxCell(3)>node%box(6).or.boxCell(6)<node%box(3)) &
            return

        if (TriBoxOverlap(c,node%the_data)) then
            aaa = aaa + 1
            CrossTri  => c%CrossTri
            if (aaa==1) then
                ALLOCATE(c%CrossTri)
                c%CrossTri%tri => node
                ! c%CrossTri%prev=> null()
                c%CrossTri%next=> null()
            else
                do ii = 1, aaa-2
                    CrossTri => CrossTri%next
                enddo
                ALLOCATE(CrossTri%next)
                CrossTri%next%tri => node
                CrossTri%next%next=> null()
            endif
        endif

        call BGCellFindTri(c,boxCell,node%right,aaa)
        call BGCellFindTri(c,boxCell,node%left ,aaa)

        end subroutine BGCellFindTri
!----------------------------------------------------------------------
        recursive subroutine initCrossCellInout(c)
        implicit none
        type(FTTCell),pointer :: c
        type(FTTCell),pointer :: cs
        integer:: is

        if(ASSOCIATED(c%Octson))then
            do is=1,c%Octson%nSon
                cs=>c%Octson%son(is)
                call initCrossCellInout(cs)
            enddo
            return
        endif

        if (c%Cross==-3) call CrossCellInout(c)

        endsubroutine initCrossCellInout
!----------------------------------------------------------------------
        subroutine CrossCellInout(c)
        use ModTools, only: BBOX, MollerTrumbore
        use ModTools, only: DistPointToTri
        implicit none
        type(FTTCell),pointer    :: c
        type(tCrossTri),pointer  :: tri, minTri
        real(R8)                 :: Dist, minDist
        type(typPoint)           :: tar

        minDist = SUM(BGCellSize)
        tar%p(1:3) = c%Center(1:3)
        tri => c%CrossTri
        do while (ASSOCIATED(tri))
            dist = DistPointToTri(tar,tri%tri%the_data)
            if (dist < minDist) then
                minDist = Dist
                minTri  =>tri
            endif
            tri => tri%next
        enddo
        if (DotProductInout(c%Center,minTri%tri%the_data)) then
            c%Cross = 1
        else
            c%Cross = 2
        endif
        endsubroutine CrossCellInout
!----------------------------------------------------------------------
        logical function DotProductInout(p,tri)
        ! Dot Poduct Method to discriminate in/out cell.
        ! T: outside; F: inside.
        use ModPrecision
        use ModKDTree
        use ModTypDef
        use ModTools, only: CROSS_PRODUCT_3
        implicit none
        real(R8),INTENT(IN)       :: p(3)
        type(triangle),INTENT(IN) :: tri
        real(R8)                  :: n(3), v0(3), v1(3), v2(3)

        v1 = tri%P(2)%P-tri%P(1)%P
        v2 = tri%P(3)%P-tri%P(2)%P
        v0 = p-tri%P(4)%P
        n  = CROSS_PRODUCT_3(v1,v2)

        if (DOT_PRODUCT(n,v0)>=0) then
            DotProductInout=.true.
        else
            DotProductInout=.false.
        endif
        endfunction DotProductInout
    end module ModMeshTools
!======================================================================
    subroutine GenerateBGMesh   ! BG -- back-ground
    use ModTypDef
    use ModMesh
    use ModInpMesh
    use ModMeshTools
    implicit none
    integer :: i, j, k
    type(FTTCell),pointer :: t

    BGCellSize=abs(DomainMin-DomainMax)/nCell
    nBGCells=nCell(1)*nCell(2)*nCell(3)
    nCells=0
    ALLOCATE(OctCell(nCell(1),nCell(2),nCell(3)))
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        t=>OctCell(i, j, k)
        nCells=nCells+1
            t%nBGCell     = [i, j, k]
            t%nCell       = nCells
            t%Location    = 0
            t%LVL(1:3)    = 0
            t%fSplitTyp  = 0
            t%cross       = -5
            t%Mark(1:3)   = .false.
            t%Center      = (/DomainMin(1)+(i-0.5)*BGCellSize(1),   &
                              DomainMin(2)+(j-0.5)*BGCellSize(2),   &
                              DomainMin(3)+(k-0.5)*BGCellSize(3)/)
            NULLIFY (t%Father, t%Octson, t%CrossTri)
    enddo
    enddo
    enddo
    endsubroutine GenerateBGMesh
!======================================================================
    subroutine initFindNeighbor
    use ModPrecision
    use ModTypDef
    use ModMesh
    use ModInpMesh
    use ModMeshTools
    use ModNeighbor
    implicit none
    integer :: i, j, k
    type(FTTCell),pointer :: t
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        t=>OctCell(i,j,k)
        call UpdateNeighbor(t,0)
    enddo
    enddo
    enddo
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "FindNeighbor time: ", tEnd-tStart
    endsubroutine initFindNeighbor
!======================================================================
    subroutine initCellInout
    use ModPrecision
    use ModTypDef
    use ModMesh
    use ModInpMesh
    use ModMeshTools
    use ModNeighbor
    implicit none
    integer :: i, j, k
    type(FTTCell),pointer :: t
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time
    real(R8):: p        ! Used to print precentage
    integer :: step=0   ! Counter, print progress precentage

    call CPU_TIME(tStart)
    write(6,'(1X,A,12X,A)',advance='no') 'OctCell inout progress:', ''
    flush(6)
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        t=>OctCell(i,j,k)
        call initCellInout2(t)
            step=step+1
            p=step/real(nBGCells,R8)*100
            write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
            flush(6)
    enddo
    enddo
    enddo
    print*,''
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "CellInout time: ", tEnd-tStart
    contains
!----------------------------------------------------------------------
        recursive subroutine initCellInout2(c)
        implicit none
        type(FTTCell),pointer :: c
        type(FTTCell),pointer :: cs
        integer:: is

        if(ASSOCIATED(c%Octson))then
            do is=1,c%Octson%nSon
                cs=>c%Octson%son(is)
                call initCellInout2(cs)
            enddo
            return
        endif

        if (c%Cross==-4 .or. c%Cross==-3) call CellInout(c)

        endsubroutine initCellInout2
    endsubroutine initCellInout
!======================================================================
    subroutine initSurfaceAdapt
    use ModMesh
    use ModMeshTools
    use ModInpMesh
    use ModNeighbor
    implicit none
    type(FTTCell),pointer :: t
    integer :: i, j, k
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time
    real(R8):: p        ! Used to print precentage
    integer :: step=0   ! Counter, print progress precentage

    call CPU_TIME(tStart)
    write(6,'(1X,A,12X,A)',advance='no') 'Cell cross progress:', ''
    flush(6)
    select case (cIntersectMethod)
    case (1)    ! Ray-cast with Painting Algorithm
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t       =>OctCell(i,j,k)
            call CellCast(t)
            if (t%cross == -3) call SurfaceAdapt(t)
            step=step+1
            p=step/real(nBGCells,R8)*100
            write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
            flush(6)
        enddo
        enddo
        enddo
        write(*,*) '' ! Stop write with advance='no'
        call CPU_TIME(tEnd)
        write(*,'(1X,A,F10.2)') "CrossCell time: ", tEnd-tStart
        call initFindNeighbor
        call initPaintingAlgorithm ! Painting Algorithm
        call initCellInout
        cIntersectMethod = 7 ! Close the Painting Algorithm
        call initSmoothMesh

    case (2)    ! Ray-cast only
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t       =>OctCell(i,j,k)
            call CellCast(t)
            if (t%cross == -3) call SurfaceAdapt(t)
            step=step+1
            p=step/real(nBGCells,R8)*100
            write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
            flush(6)
        enddo
        enddo
        enddo
        write(*,*) '' ! Stop write with advance='no'
        call CPU_TIME(tEnd)
        write(*,'(1X,A,F10.2)') "CrossCell time: ", tEnd-tStart
        call initFindNeighbor
        call initCellInout
        cIntersectMethod = 7 ! Close the Painting Algorithm
        call initSmoothMesh

    case (3)    ! AABB with Painting Algorithm
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t       =>OctCell(i,j,k)
            call BGCellAABB(t)
            if (t%cross == -3) call SurfaceAdapt(t)
            call initCrossCellInout(t)
            step=step+1
            p=step/real(nBGCells,R8)*100
            write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
            flush(6)
        enddo
        enddo
        enddo
        write(*,*) '' ! Stop write with advance='no'
        call CPU_TIME(tEnd)
        write(*,'(1X,A,F10.2)') "CrossCell time: ", tEnd-tStart
        call initFindNeighbor
        call initPaintingAlgorithm ! Painting Algorithm
        cIntersectMethod = 8 ! Close the Painting Algorithm
        call initSmoothMesh

    case (4)    ! AABB only
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t       =>OctCell(i,j,k)
            call BGCellAABB(t)
            if (t%cross == -3) call SurfaceAdapt(t)
            call initCrossCellInout(t)
            step=step+1
            p=step/real(nBGCells,R8)*100
            ! if (p>50) return
            write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
            flush(6)
        enddo
        enddo
        enddo
        !call TDCGOutput('OK')
        write(*,*) '' ! Stop write with advance='no'
        call CPU_TIME(tEnd)
        write(*,'(1X,A,F10.2)') "CrossCell time: ", tEnd-tStart
        call initFindNeighbor
        call initCellInout
        cIntersectMethod = 8 ! Close the Painting Algorithm
        call initSmoothMesh

    ! case (5)    ! AABBTraverse
    !     do k = 1, nCell(3)
    !     do j = 1, nCell(2)
    !     do i = 1, nCell(1)
    !         t       =>OctCell(i,j,k)
    !         call AABBTraverse(t)
    !         if (t%cross == -3) call SurfaceAdapt(t)
    !         call initCrossCellInout(t)
    !         step=step+1
    !         p=step/real(nBGCells,R8)*100
    !         write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
    !         flush(6)
    !     enddo
    !     enddo
    !     enddo
    !     write(*,*) '' ! Stop write with advance='no'
    !     call CPU_TIME(tEnd)
    !     write(*,'(1X,A,F10.2)') "CrossCell time: ", tEnd-tStart
    !     call initFindNeighbor
    !     call initCellInout
    !     call initSmoothMesh

    ! case (6)    ! Ray-cast only
    !     do k = 1, nCell(3)
    !     do j = 1, nCell(2)
    !     do i = 1, nCell(1)
    !         t       =>OctCell(i,j,k)
    !         call CellCast(t)
    !         if (t%cross == -3) call SurfaceAdapt(t)
    !         call initCrossCellInout(t)
    !         step=step+1
    !         p=step/real(nBGCells,R8)*100
    !         write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
    !         flush(6)
    !     enddo
    !     enddo
    !     enddo
    !     write(*,*) '' ! Stop write with advance='no'
    !     call CPU_TIME(tEnd)
    !     write(*,'(1X,A,F10.2)') "CrossCell time: ", tEnd-tStart
    !     call initFindNeighbor
    !     call initCellInout
    !     call initSmoothMesh
    end select

    contains
!----------------------------------------------------------------------
        recursive subroutine SurfaceAdapt(c)
        implicit none
        type(FTTCell),pointer :: c
        type(FTTCell),pointer :: cs
        integer:: is
        integer :: ii

        do ii=1,3
            if (c%LVL(ii) >= InitRefineLVL) return
        enddo
        if (c%Cross == 1 .or. c%Cross == 2 .or. c%Cross == -3)then
            call NewCell(c,0)
            do is=1,c%Octson%nSon
                cs=>c%Octson%son(is)
                call SurfaceAdapt(cs)
            enddo
        endif
        endsubroutine SurfaceAdapt
    endsubroutine initSurfaceAdapt
!======================================================================
    subroutine initPaintingAlgorithm
    use ModMesh
    use ModMeshTools
    use ModInpMesh
    use ModNeighbor
    use ModPrecision
    implicit none
    type(FTTCell),pointer :: t
    integer :: ii, jj, kk
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time
    real(R8):: p        ! Used to print precentage
    integer :: step=0   ! Counter, print progress precentage

    call CPU_TIME(tStart)

    loop:  do kk = 1, nCell(3)
    do jj = 1, nCell(2)
    do ii = 1, nCell(1)
        t       =>OctCell(ii,jj,kk)
        if (PAFindSeed(t)) exit loop
    enddo
    enddo
    enddo loop

    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.8)') "Painting Algorithm time: ", tEnd-tStart
    contains
!----------------------------------------------------------------------
        recursive function PAFindSeed(c) result(PA)
        ! PA: Have run/not run the subroutine PaintingAlgorithm
        implicit none
        logical               :: PA
        type(FTTCell),pointer :: c
        type(FTTCell),pointer :: cn
        type(FTTCell),pointer :: cs
        integer               :: is
        integer               :: i
        logical,SAVE          :: PAused=.false.
        ! If used PaintingAlgorithm, PAused=.T.

        if (PAused) return

        if(ASSOCIATED(c%Octson))then
            do is=1,c%Octson%nSon
                cs=>c%Octson%son(is)
                PA = PAFindSeed(cs)
            enddo
            return
        endif

        PA = .false.
        call CellInout(c)
        if (c%Cross==0) then
            do i = 1, 6
                cn => Neighbor(c,i)
                if (ASSOCIATED(cn)) call PaintingAlgorithm(cn,i)
            enddo
            PA = .true.
            PAused = .true.
        endif
        endfunction PAFindSeed
!----------------------------------------------------------------------
        recursive subroutine PaintingAlgorithm(c,dirct)
        implicit none
        type(FTTCell),pointer :: c, cc
        type(FTTCell),pointer :: cn, cnn
        type(FTTCell),pointer :: cs
        integer,INTENT(IN)    :: dirct
        integer:: is
        integer:: i, j, k

        if(ASSOCIATED(c%Octson))then
            do is=1,c%Octson%nSon
                cs=>c%Octson%son(is)
                call PaintingAlgorithm(cs,dirct)
            enddo
            return
        endif

        if (c%Cross/=-4) return ! OctCell has been paintted

        loop1: do i = 1, 6
            cn=>Neighbor(c,i)
            if (.not.ASSOCIATED(cn)) cycle loop1
            if (cn%cross==0) then
                exit loop1
            elseif (cn%cross/=-4) then
                if (ASSOCIATED(cn%Octson)) then
                    if (mod(i,2)==1) then
                        k = i + 1
                    else
                        k = i - 1
                    endif
                    do j = 1, cn%Octson%nSon
                        if (cn%Octson%son(j)%cross==0) then
                            cnn => cn%Octson%son(j)
                            cnn => Neighbor(cnn,k)
                            if (cnn%nCell==c%nCell) exit loop1
                        endif
                    enddo
                endif
            endif
            if (i == 6) return
        enddo loop1

        c%Cross=0
        do i = 1, 6
            if (mod(i,2)==1) then
                k = i + 1
            else
                k = i - 1
            endif
            if (dirct /= k) then
                cc => Neighbor(c,i)
                if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,i)
            endif
        enddo
        endsubroutine PaintingAlgorithm
!----------------------------------------------------------------------
        recursive subroutine PAMarkCrossCell(c)
        implicit none
        type(FTTCell),pointer :: c
        type(FTTCell),pointer :: cs
        integer:: is

        if (c%Cross/=-3) return

        if(ASSOCIATED(c%Octson))then
            do is=1,c%Octson%nSon
                cs=>c%Octson%son(is)
                call PAMarkCrossCell(cs)
            enddo
            return
        endif

        call CellInout(c)

        endsubroutine PAMarkCrossCell
    endsubroutine initPaintingAlgorithm
!======================================================================
    subroutine initSmoothMesh
    use ModPrecision
    use ModMesh
    use ModTypDef
    use ModMeshTools
    use ModInpMesh
    implicit none
    type(FTTCell),POINTER :: t
    integer               :: ii, jj, kk, iii
    logical               :: iost
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    do iii=1, InitRefineLVL-1

        do kk = 1, nCell(3)
        do jj = 1, nCell(2)
        do ii = 1, nCell(1)
            t=>OctCell(ii, jj, kk)
            call PreSmoothMesh(t) ! Mark neighbor OctCell if need refine.
        enddo
        enddo
        enddo

        iost = .false.
        do kk = 1, nCell(3)
        do jj = 1, nCell(2)
        do ii = 1, nCell(1)
            t=>OctCell(ii, jj, kk)
            call SmoothMesh(t,iost) ! Do smooth.
        enddo
        enddo
        enddo

        if (.not.iost) exit
    enddo
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "Smooth Mesh time: ", tEnd-tStart
    contains
!----------------------------------------------------------------------
        recursive subroutine PreSmoothMesh(c)
        use ModNeighbor, only: Neighbor
        implicit none
        type(FTTCell),POINTER :: c
        type(FTTCell),POINTER :: cn
        type(FTTCell),pointer :: cs
        integer:: is
        integer:: i
        ! mark(3)  1 x refine; 2 y refine; 3 z refine;

        if(ASSOCIATED(c%Octson))then
            do is=1,c%Octson%nSon
                cs=>c%Octson%son(is)
                call PreSmoothMesh(cs)
            enddo
            return
        endif

        do i = 1, 6
            cn => Neighbor(c,i)
            if (ASSOCIATED(cn)) then
                if (cn%LVL(1)+1<c%LVL(1)) cn%mark(1)=.true.
            endif
        enddo

        ! Hole OctCell
        ! if (.not.mark(1)) then
        !     if (ASSOCIATED(c%Octson%NeighborX1).and.&
                  !ASSOCIATED(c%Octson%NeighborX2))then
        !         if (c%Octson%NeighborX1%LVL(1)>c%LVL(1) .and. &
        !             c%Octson%NeighborX2%LVL(1)>c%LVL(1)) then
        !             mark(1)=.true.
        !         elseif (c%Octson%NeighborX1%LVL(1)<c%LVL(1) .and. &
        !             c%Octson%NeighborX2%LVL(1)<c%LVL(1)) then
        !             mark(4)=.true.
        !         endif
        !     endif
        ! endif
        endsubroutine PreSmoothMesh
!----------------------------------------------------------------------
        recursive subroutine SmoothMesh(c,iosout)
        implicit none
        type(FTTCell),POINTER :: c
        logical,INTENT(INOUT) :: iosout
        logical               :: ios
        type(FTTCell),pointer :: cs
        integer:: is
        integer:: k
        ! mark(3)  1 x refine; 2 y refine; 3 z refine;

        if (maxval(c%LVL)>=InitRefineLVL) return


        if(ASSOCIATED(c%Octson))then
            do is=1,c%Octson%nSon
                cs=>c%Octson%son(is)
                call SmoothMesh(cs,ios)
            enddo
            return
        endif

        if (c%Cross/=0) return
        ! Refine
        ios = .false.
        if     ( c%mark(1) .and. c%mark(2) .and. c%mark(3)) then
            call NewCell(c,0)
            ios = .true.
        elseif ( c%mark(1) .and. .not.c%mark(2) .and. .not.c%mark(3) ) then
            call NewCell(c,1)
            ios = .true.
        elseif ( .not.c%mark(1) .and. c%mark(2) .and. .not.c%mark(3) ) then
            call NewCell(c,2)
            ios = .true.
        elseif ( .not.c%mark(1) .and. .not.c%mark(2) .and. c%mark(3) ) then
            call NewCell(c,3)
            ios = .true.
        elseif ( c%mark(1) .and. c%mark(2) .and. .not.c%mark(3) ) then
            call NewCell(c,4)
            ios = .true.
        elseif ( c%mark(1) .and. .not.c%mark(2) .and. c%mark(3) ) then
            call NewCell(c,5)
            ios = .true.
        elseif ( .not.c%mark(1) .and. c%mark(2) .and. c%mark(3) ) then
            call NewCell(c,6)
            ios = .true.
        endif

        if (ios) then
            c%mark(1:3) = .false.
            iosout      = .true.
            call UpdateNeighbor(c,0)
            cs=>c%Octson%Neighbor1
            if (ASSOCIATED(cs)) call UpdateNeighbor(cs,2)
            cs=>c%Octson%Neighbor2
            if (ASSOCIATED(cs)) call UpdateNeighbor(cs,1)
            cs=>c%Octson%Neighbor3
            if (ASSOCIATED(cs)) call UpdateNeighbor(cs,4)
            cs=>c%Octson%Neighbor4
            if (ASSOCIATED(cs)) call UpdateNeighbor(cs,3)
            cs=>c%Octson%Neighbor5
            if (ASSOCIATED(cs)) call UpdateNeighbor(cs,6)
            cs=>c%Octson%Neighbor6
            if (ASSOCIATED(cs)) call UpdateNeighbor(cs,5)
        endif
        endsubroutine SmoothMesh
    endsubroutine initSmoothMesh
!======================================================================
!======================================================================
!----------------------------------------------------------------------
