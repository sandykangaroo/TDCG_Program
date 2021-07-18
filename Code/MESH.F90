!======================================================================
    module ModMesh
        use ModTypDef
        implicit none
        integer :: nBGCells     ! Number of the background Cells.
        integer :: nCells     ! Number of total Cells.
        real(R8):: BGCellSize(3)    ! Step size for background OctCell.
        type(typOctCell),pointer :: OctCell(:,:,:)
        ! type(typStructCell),pointer:: StruCell(:,:,:,:)
    endmodule ModMesh
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
        type(typOctCell),pointer :: c

        select case (cIntersectMethod)
        case (1:2)
            call CellCast(c)
        case (3:4)
            call AABB(c)
        case (5)
            call AABBTraverse(c)
        case (6)
            call CellCast(c)
        case (7) ! After surface refine, change cIntersectMethod=1,2 to =7
            call CellCast(c)
            if (c%cross == -4) call CellInout(c)
        case (8) ! After surface refine, change cIntersectMethod=3,4 to =8
            call AABB(c)
            if (c%cross == -4) call CellInout(c)
        end select
        endsubroutine initCellCross
!----------------------------------------------------------------------
        subroutine CellCast(c)
        use ModKDTree
        use ModInpGlobal
        use ModInpMesh,only: cIntersectMethod
        implicit none
        type(typOctCell),pointer :: c
        type(typPoint)        :: p(8)
        real(R8):: x, y, z, dx, dy, dz
        integer :: i, ii, iii
        integer :: Pintersect ! Number of the intersect point
        integer :: nIntersect ! Intersect times for one point
        integer :: Rintersect ! Once ray intersects, nRay++
        real(R8):: k(3)

        x=c%Center(1); y=c%Center(2); z=c%Center(3)
        dx=BGCellSize(1)/2**(c%lvl(1)+1)
        dy=BGCellSize(2)/2**(c%lvl(2)+1)
        dz=BGCellSize(3)/2**(c%lvl(3)+1)
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
                    if (cIntersectMethod/=6) then
                        call RayCast(CSHIFT(p(i)%P,iii-1), &
                                    KDTree(ii)%root,k,nIntersect,iii)
                    else
                        call RayCastTraverse(CSHIFT(p(i)%P,iii-1), &
                                    KDTree(ii)%root,k,nIntersect,iii)
                    endif
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
                call CellInout(c)
                return
            endif
        enddo
        ! Not intersect
        c%cross = -4
        endsubroutine CellCast
!----------------------------------------------------------------------
        subroutine CellInout(c)
        use ModKDTree
        use ModInpGlobal
        use ModInpMesh,only: cIntersectMethod
        use ModTools, only: BBOX
        implicit none
        type(typOctCell),pointer :: c
        integer :: nIntersect ! Intersect times for one point
        integer :: Rintersect ! Once ray intersects, Rintersect++
        integer :: iii ! Ray
        real(R8):: k(3)

        ! if (c%cross == -4) then
        !     CellInout = 0
        ! else
        !     CellInout = 1
        ! endif
        ! ! Quick Bounding-BOX identify
        ! if (BBOX(c%Center,KDTree(1)%root%box)) then
        !     Rintersect=0
        !     do iii = 1, 3  !nRays
        !         nIntersect=0
        !         k = CSHIFT([0.,0.,1.],4-iii)
        !             call RayCast(CSHIFT(c%Center,iii-1), &
        !                         KDTree(1)%root,k,nIntersect,iii)
        !         if (mod(nIntersect,2)==0) then
        !             Rintersect = Rintersect + 0
        !         else
        !             Rintersect = Rintersect + 1
        !         endif
        !     enddo
        !     do iii = 1, 3  !nRays
        !         nIntersect=0
        !         k = CSHIFT([0.,0.,-1.],4-iii)
        !             call RayCastTraverse(CSHIFT(c%Center,iii-1), &
        !                         KDTree(1)%root,k,nIntersect,iii)
        !         if (mod(nIntersect,2)==0) then
        !             Rintersect = Rintersect + 0
        !         else
        !             Rintersect = Rintersect + 1
        !         endif
        !     enddo
        !         nIntersect=0
        !         k = [1.,1.,1.]
        !             call RayCastTraverse(CSHIFT(c%Center,0), &
        !                         KDTree(1)%root,k,nIntersect,1)
        !         if (mod(nIntersect,2)==0) then
        !             Rintersect = Rintersect + 0
        !         else
        !             Rintersect = Rintersect + 1
        !         endif
        !     if ( Rintersect > 3 ) then ! inside
        !         if (c%cross == -4) then
        !             CellInout = 3
        !         else
        !             CellInout = 2
        !         endif
        !     endif
        ! endif

        if (c%cross == -4) then
            c%cross = 0
        elseif (c%cross == -3) then
            c%cross = 1
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
                if (cIntersectMethod/=6) then
                    call RayCast(CSHIFT(c%Center,iii-1), &
                                KDTree(1)%root,k,nIntersect,iii)
                else
                    call RayCastTraverse(CSHIFT(c%Center,iii-1), &
                                KDTree(1)%root,k,nIntersect,iii)
                endif
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
                if (c%cross == 0) then
                    c%cross = 3
                elseif (c%cross == 1) then
                    c%cross = 2
                else
                    stop 'Subroutine CellInout error2'
                endif
            endif
        endif
        endsubroutine CellInout
!----------------------------------------------------------------------
        recursive subroutine RayCast(point,tree,k,nIntersect,i)
        ! Back nIntersect -> times of intersect
        use ModGeometry
        use ModKDTree
        implicit none
        real(R8),INTENT(IN)    :: point(3)     ! point(3) = x, y, z.
        type(KDT_node),pointer,INTENT(IN) :: tree
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
                if (MollerTrumbore(CSHIFT(point,1-i),k,tree%the_data,.false.)) then
                    nIntersect = nIntersect + 1
                endif
                ! Into the next tree
                call RayCast(point,tree%right,k,nIntersect,i)
                call RayCast(point,tree%left,k,nIntersect,i)
        endif
        ! box = CSHIFT(tree%box(1:3),i-1)
        ! if (point(1)>box(1).and.point(2)>box(2)) then ! Box min
        !     box = CSHIFT(tree%box(4:6),i-1)
        !     if(point(1)<box(1).and.point(2)<box(2).and. & ! Box max
        !        point(3)<box(3)) then ! Positive direction
        !         if (MollerTrumbore(CSHIFT(point,1-i),k,tree%the_data,.false.)) then
        !             nIntersect = nIntersect + 1
        !         endif
        !         ! Into the next tree
        !         call RayCast(point,tree%right,k,nIntersect,i)
        !         call RayCast(point,tree%left,k,nIntersect,i)
        !     endif
        ! endif

        endsubroutine RayCast
!----------------------------------------------------------------------
        recursive subroutine RayCastTraverse(point,tree,k,nIntersect,i)
        use ModGeometry
        use ModKDTree
        implicit none
        real(R8),INTENT(IN)    :: point(3)     ! point(3) = x, y, z.
        type(KDT_node),pointer,INTENT(IN) :: tree
        real(R8),INTENT(IN)    :: k(3)
        integer, INTENT(INOUT) :: nIntersect ! Save the intersect number
        integer, INTENT(IN)    :: i ! Which Ray
        real(R8)               :: box(6)

        if (.not.ASSOCIATED(tree)) return

        if (MollerTrumbore(CSHIFT(point,1-i),k,tree%the_data,.false.)) then
            nIntersect = nIntersect + 1
        endif
        ! Into the next tree
        call RayCastTraverse(point,tree%right,k,nIntersect,i)
        call RayCastTraverse(point,tree%left,k,nIntersect,i)
        endsubroutine RayCastTraverse
!----------------------------------------------------------------------
        logical function MollerTrumbore(Point,D,tri,ios)
        use ModTools
        implicit none
        real(R8)      ,INTENT(IN):: point(3)   ! point(3) = x, y, z.
        real(R8)      ,INTENT(IN):: D(3)
        type(triangle),INTENT(IN):: tri
        logical       ,INTENT(IN):: ios ! T: Called by CCIspecial.
                                        ! F: Called by others.
        real(R8):: V0(3), V1(3), V2(3)
        real(R8):: P(3), Q(3)
        real(R8):: t, u, v
        real(R8):: A

        V0 = point      - tri%P(1)%P
        V1 = tri%P(2)%P - tri%P(1)%P
        V2 = tri%P(3)%P - tri%P(1)%P
        P = CROSS_PRODUCT_3(D,V2)
        Q = CROSS_PRODUCT_3(V0,V1)
        A = DOT_PRODUCT(P,V1)   ! A=0 is impossible.
        u = DOT_PRODUCT(P,V0)
        v = DOT_PRODUCT(Q,D)
        t = DOT_PRODUCT(Q,V2)
        ! Add if() judgment before each operation is no necessary, for 
        ! the reason that BBOX judgment leads to a high probability of 
        ! return value = .True.
        if (ios) then
            if (A>0 .and. u>=0 .and. v>=0 .and. u+v<=A .and. &
                t>=0 .and. t<=A) then
                MollerTrumbore=.true.; return
            elseif (A<0 .and. u<=0 .and. v<=0 .and. u+v>=A .and. &
                t<=0 .and. t>=A) then
                MollerTrumbore=.true.; return
            else
                MollerTrumbore=.false.; return
            endif
        endif
        if (A>0 .and. u>=0 .and. v>=0 .and. u+v<=A .and. t>=0) then
            MollerTrumbore=.true.; return
        elseif (A<0 .and. u<=0 .and. v<=0 .and. u+v>=A .and. t<=0) then
            MollerTrumbore=.true.; return
        else
            MollerTrumbore=.false.; return
        endif
        endfunction MollerTrumbore
!----------------------------------------------------------------------
        subroutine AABB(c)
        use ModKDTree
        use ModInpGlobal
        use ModGeometry
        implicit none
        type(typOctCell),pointer :: c
        type(triangle)           :: tri
        real(R8)                 :: boxCell(6)
        integer                  :: aaa

        aaa=0
        boxCell(1)=c%center(1)-BGCellSize(1)/2.**(c%lvl(1)+1)
        boxCell(2)=c%center(2)-BGCellSize(2)/2.**(c%lvl(2)+1)
        boxCell(3)=c%center(3)-BGCellSize(3)/2.**(c%lvl(3)+1)
        boxCell(4)=c%center(1)+BGCellSize(1)/2.**(c%lvl(1)+1)
        boxCell(5)=c%center(2)+BGCellSize(2)/2.**(c%lvl(2)+1)
        boxCell(6)=c%center(3)+BGCellSize(3)/2.**(c%lvl(3)+1)

        call AABBFindTri(c,boxCell,c%Father%CrossTri,aaa)
        ! c%nCrossTri = aaa
        if (aaa==0) then
            c%cross = -4
        else
            c%cross = -3
        endif
        return
        end subroutine AABB
!----------------------------------------------------------------------
        subroutine AABBTraverse(c)!1=intersect,0=no intersect
        use ModKDTree
        use ModInpGlobal
        use ModGeometry
        implicit none
        type(typOctCell),pointer :: c
        type(triangle)        :: tri
        integer                     :: ng, i
        real(R8)              ::boxCell(6)

        c%cross=-4
        Loop1:do ng=1,nGeometry
            do i=1, body(ng)%nse
                tri= body(ng)%se3d(i)
                if(TriBoxOverlap (c,tri))then
                    boxCell(1)=c%center(1)-BGCellSize(1)/2**(c%lvl(1)+1)
                    boxCell(2)=c%center(2)-BGCellSize(2)/2**(c%lvl(2)+1)
                    boxCell(3)=c%center(3)-BGCellSize(3)/2**(c%lvl(3)+1)
                    boxCell(4)=c%center(1)+BGCellSize(1)/2**(c%lvl(1)+1)
                    boxCell(5)=c%center(2)+BGCellSize(2)/2**(c%lvl(2)+1)
                    boxCell(6)=c%center(3)+BGCellSize(3)/2**(c%lvl(3)+1)
                    exit Loop1
                endif
            enddo
        enddo Loop1
        return         
        end subroutine AABBTraverse
!----------------------------------------------------------------------
        recursive subroutine AABBFindTri(c,boxCell,node,aaa)
        use ModKDTree
        implicit none
        type(typOctCell),pointer    :: c
        real(R8),INTENT(IN)         :: boxCell(6)
        type(typCrossTri), pointer  :: node
        integer, INTENT(INOUT)      :: aaa
        integer                     :: i
        type(typCrossTri), pointer  :: CrossTri

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
            
            type(typOctCell),pointer :: c
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
        
            boxhalfsize(1)= BGCellSize(1)/2**(c%lvl(1)+1)
            boxhalfsize(2)= BGCellSize(2)/2**(c%lvl(2)+1)
            boxhalfsize(3)= BGCellSize(3)/2**(c%lvl(3)+1)
            
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
    !/* find min, max of the triangle each direction, and test for overlap in */
    !/* that direction -- this is equivalent to testing a minimal AABB around */
    !/*  the triangle against the AABB */
            
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
        use ModInpGlobal
        use ModNeighbor
        implicit none
        integer(I4),INTENT(IN):: split
        integer(I4)           :: spl
        type(typOctCell),pointer :: c
        real(R8)              :: x, y, z, dx, dy, dz

        spl=split
        if (.not.Aniso) spl=0
        select case (spl)
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
            ! c%son1%Node=0
            ! c%son2%Node=0
            ! c%son3%Node=0
            ! c%son4%Node=0
            ! c%son5%Node=0
            ! c%son6%Node=0
            ! c%son7%Node=0
            ! c%son8%Node=0
            c%son1%Mark        = .false.
            c%son2%Mark        = .false.
            c%son3%Mark        = .false.
            c%son4%Mark        = .false.
            c%son5%Mark        = .false.
            c%son6%Mark        = .false.
            c%son7%Mark        = .false.
            c%son8%Mark        = .false.
            x=c%Center(1); y=c%Center(2); z=c%Center(3)
            dx=BGCellSize(1)/2**(c%lvl(1)+2)
            dy=BGCellSize(2)/2**(c%lvl(2)+2)
            dz=BGCellSize(3)/2**(c%lvl(3)+2)
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
            c%son1%Father=>c
            c%son2%Father=>c
            c%son3%Father=>c
            c%son4%Father=>c
            c%son5%Father=>c
            c%son6%Father=>c
            c%son7%Father=>c
            c%son8%Father=>c
            ! c%son1%nCrossTri=0
            ! c%son2%nCrossTri=0
            ! c%son3%nCrossTri=0
            ! c%son4%nCrossTri=0
            ! c%son5%nCrossTri=0
            ! c%son6%nCrossTri=0
            ! c%son7%nCrossTri=0
            ! c%son8%nCrossTri=0
            if (c%cross==1 .or. c%cross==2 .or. c%cross==-3) then
                call initCellCross(c%son1)
                call initCellCross(c%son2)
                call initCellCross(c%son3)
                call initCellCross(c%son4)
                call initCellCross(c%son5)
                call initCellCross(c%son6)
                call initCellCross(c%son7)
                call initCellCross(c%son8)
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
            call NullifyCell(c%son1)
            call NullifyCell(c%son2)
            call NullifyCell(c%son3)
            call NullifyCell(c%son4)
            call NullifyCell(c%son5)
            call NullifyCell(c%son6)
            call NullifyCell(c%son7)
            call NullifyCell(c%son8)
            return
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
            ! c%son1%Node=0
            ! c%son2%Node=0
            c%son1%Mark        = .false.
            c%son2%Mark        = .false.
            x=c%Center(1); y=c%Center(2); z=c%Center(3)
            dx=BGCellSize(1)/2**(c%lvl(1)+2)
            dy=BGCellSize(2)/2**(c%lvl(2)+2)
            dz=BGCellSize(3)/2**(c%lvl(3)+2)
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
            c%son1%Father=>c
            c%son2%Father=>c
            ! c%son1%nCrossTri=0
            ! c%son2%nCrossTri=0
            if (c%cross==1 .or. c%cross==2 .or. c%cross==-3) then
                call initCellCross(c%son1)
                call initCellCross(c%son2)
            else
                c%son1%cross=c%cross
                c%son2%cross=c%cross
            endif
            call NullifyCell(c%son1)
            call NullifyCell(c%son2)
            return
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
            ! c%son1%Node=0
            ! c%son2%Node=0
            ! c%son3%Node=0
            ! c%son4%Node=0
            c%son1%Mark        = .false.
            c%son2%Mark        = .false.
            c%son3%Mark        = .false.
            c%son4%Mark        = .false.
            x=c%Center(1); y=c%Center(2); z=c%Center(3)
            dx=BGCellSize(1)/2**(c%lvl(1)+2)
            dy=BGCellSize(2)/2**(c%lvl(2)+2)
            dz=BGCellSize(3)/2**(c%lvl(3)+2)
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
            c%son1%Father=>c
            c%son2%Father=>c
            c%son3%Father=>c
            c%son4%Father=>c
            ! c%son1%nCrossTri=0
            ! c%son2%nCrossTri=0
            ! c%son3%nCrossTri=0
            ! c%son4%nCrossTri=0
            if (c%cross==1 .or. c%cross==2 .or. c%cross==-3) then
                call initCellCross(c%son1)
                call initCellCross(c%son2)
                call initCellCross(c%son3)
                call initCellCross(c%son4)
            else
                c%son1%cross=c%cross
                c%son2%cross=c%cross
                c%son3%cross=c%cross
                c%son4%cross=c%cross
            endif
            call NullifyCell(c%son1)
            call NullifyCell(c%son2)
            call NullifyCell(c%son3)
            call NullifyCell(c%son4)
            return
        end select
        endsubroutine NewCell
!----------------------------------------------------------------------
        recursive subroutine FindNeighbor(c)
        use ModNeighbor
        implicit none
        type(typOctCell),pointer :: c

        if (.not.ASSOCIATED(c)) return

        c%NeighborX1=>NeighborX1(c)
        c%NeighborX2=>NeighborX2(c)
        c%NeighborY1=>NeighborY1(c)
        c%NeighborY2=>NeighborY2(c)
        c%NeighborZ1=>NeighborZ1(c)
        c%NeighborZ2=>NeighborZ2(c)

        call FindNeighbor(c%son1)
        call FindNeighbor(c%son2)
        call FindNeighbor(c%son3)
        call FindNeighbor(c%son4)
        call FindNeighbor(c%son5)
        call FindNeighbor(c%son6)
        call FindNeighbor(c%son7)
        call FindNeighbor(c%son8)
        endsubroutine FindNeighbor
!----------------------------------------------------------------------
        recursive subroutine UpdateNeighbor(c,dirct)
        USE ModNeighbor
        implicit none
        type(typOctCell),pointer :: c
        integer,INTENT(IN)    ::dirct

        if (.not.ASSOCIATED(c)) return

        select case (dirct)
        case (0)
            c%NeighborX1=>NeighborX1(c)
            c%NeighborX2=>NeighborX2(c)
            c%NeighborY1=>NeighborY1(c)
            c%NeighborY2=>NeighborY2(c)
            c%NeighborZ1=>NeighborZ1(c)
            c%NeighborZ2=>NeighborZ2(c)
        case (1)
            c%NeighborX1=>NeighborX1(c)
        case (2)
            c%NeighborX2=>NeighborX2(c)
        case (3)
            c%NeighborY1=>NeighborY1(c)
        case (4)
            c%NeighborY2=>NeighborY2(c)
        case (5)
            c%NeighborZ1=>NeighborZ1(c)
        case (6)
            c%NeighborZ2=>NeighborZ2(c)
        end select

        call UpdateNeighbor(c%son1,dirct)
        call UpdateNeighbor(c%son2,dirct)
        call UpdateNeighbor(c%son3,dirct)
        call UpdateNeighbor(c%son4,dirct)
        call UpdateNeighbor(c%son5,dirct)
        call UpdateNeighbor(c%son6,dirct)
        call UpdateNeighbor(c%son7,dirct)
        call UpdateNeighbor(c%son8,dirct)

        endsubroutine UpdateNeighbor
!----------------------------------------------------------------------
        subroutine DeletCell(c,split)
        implicit none
        type(typOctCell),pointer :: c
        integer(I4),INTENT(IN):: split
        endsubroutine DeletCell
!----------------------------------------------------------------------
        subroutine NullifyCell(c)
        implicit none
        type(typOctCell),pointer :: c
        NULLIFY(c%son1, c%son2, c%son3, c%son4,                 &
                c%son5, c%son6, c%son7, c%son8,                 &
                c%NeighborX1, c%NeighborX2, c%NeighborY1,       &
                c%NeighborY2, c%NeighborZ1, c%NeighborZ2)
        endsubroutine NullifyCell
!----------------------------------------------------------------------
        ! recursive subroutine InitialCellMark(c)
        ! implicit none
        ! type(typOctCell),pointer :: c
        ! if(ASSOCIATED(c%son8))then
        !     call InitialCellMark(c%son1)
        !     call InitialCellMark(c%son2)
        !     call InitialCellMark(c%son3)
        !     call InitialCellMark(c%son4)
        !     call InitialCellMark(c%son5)
        !     call InitialCellMark(c%son6)
        !     call InitialCellMark(c%son7)
        !     call InitialCellMark(c%son8)
        !     return
        ! elseif(ASSOCIATED(c%son4))then
        !     call InitialCellMark(c%son1)
        !     call InitialCellMark(c%son2)
        !     call InitialCellMark(c%son3)
        !     call InitialCellMark(c%son4)
        !     return
        ! elseif(ASSOCIATED(c%son2))then
        !     call InitialCellMark(c%son1)
        !     call InitialCellMark(c%son2)
        !     return
        ! endif
        ! c%mark=.false.
        ! endsubroutine InitialCellMark
!----------------------------------------------------------------------
        subroutine BGCellAABB(c)
        use ModKDTree
        use ModInpGlobal
        use ModGeometry
        implicit none
        type(typOctCell),pointer    :: c
        type(triangle)              :: tri
        real(R8)                    :: boxCell(6)
        integer                     :: aaa
        integer                     :: ng
        type(typCrossTri), pointer  :: CrossTri
        INTEGER                     :: ii
        type(typOctCell), pointer   :: cc

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
            c%cross = -4
        else
            c%cross = -3
        endif

        end subroutine BGCellAABB
!----------------------------------------------------------------------
        recursive subroutine BGCellFindTri(c,boxCell,node,aaa)
        use ModKDTree
        implicit none
        type(typOctCell),pointer    :: c
        real(R8),INTENT(IN)         :: boxCell(6)
        type(KDT_node), pointer     :: node
        integer, INTENT(INOUT)      :: aaa
        integer                     :: ii
        type(typCrossTri), pointer  :: CrossTri

        if (.not. ASSOCIATED(node)) return

        if (boxCell(1)>node%box(4).or.boxCell(4)<node%box(1).or.&
            boxCell(2)>node%box(5).or.boxCell(5)<node%box(2).or.&
            boxCell(3)>node%box(6).or.boxCell(6)<node%box(3)) then
            return
        else
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
                    ! CrossTri%next%prev=> CrossTri
                    CrossTri%next%next=> null()
                endif
            endif

            call BGCellFindTri(c,boxCell,node%right,aaa)
            call BGCellFindTri(c,boxCell,node%left ,aaa)
        endif

        end subroutine BGCellFindTri
!----------------------------------------------------------------------
        recursive subroutine initCrossCellInout(t)
        implicit none
        type(typOctCell),pointer :: t

        if(ASSOCIATED(t%son8))then
            call initCrossCellInout(t%son1)
            call initCrossCellInout(t%son2)
            call initCrossCellInout(t%son3)
            call initCrossCellInout(t%son4)
            call initCrossCellInout(t%son5)
            call initCrossCellInout(t%son6)
            call initCrossCellInout(t%son7)
            call initCrossCellInout(t%son8)
            return
        elseif(ASSOCIATED(t%son4))then
            call initCrossCellInout(t%son1)
            call initCrossCellInout(t%son2)
            call initCrossCellInout(t%son3)
            call initCrossCellInout(t%son4)
            return
        elseif(ASSOCIATED(t%son2))then
            call initCrossCellInout(t%son1)
            call initCrossCellInout(t%son2)
            return
        endif

        if (t%cross==-3) call CrossCellInout(t)

        endsubroutine initCrossCellInout
!----------------------------------------------------------------------
        subroutine CrossCellInout(c)
        use ModGlobalConstants, only : epsR8
        use ModTools, only: BBOX
        implicit none
        type(typOctCell),pointer :: c
        INTEGER                  :: nTri ! Num. of CrossTri
        INTEGER                  :: nDP  ! Num. of DotProductInout = F.
        type(typCrossTri),pointer:: p
        logical                  :: ios
        type(typCrossTri),pointer:: pp
        real(R8)                 :: n(3)
        integer                  :: b ! counter for loop22
        real(R8)                 :: boxCell(6)
        integer                  :: i
        integer                  :: pointOrder
        real(R8)                 :: point(3)
        type(triangle), pointer  :: tri
        real(R8)                 :: derta
        integer:: aia,aib
        aia=0
        aib=0

        nTri = 1
        p    => c%CrossTri
        boxCell(1)=c%center(1)-BGCellSize(1)/2.**(c%lvl(1)+1)
        boxCell(2)=c%center(2)-BGCellSize(2)/2.**(c%lvl(2)+1)
        boxCell(3)=c%center(3)-BGCellSize(3)/2.**(c%lvl(3)+1)
        boxCell(4)=c%center(1)+BGCellSize(1)/2.**(c%lvl(1)+1)
        boxCell(5)=c%center(2)+BGCellSize(2)/2.**(c%lvl(2)+1)
        boxCell(6)=c%center(3)+BGCellSize(3)/2.**(c%lvl(3)+1)

        loop11: do while (ASSOCIATED(p))
            tri => p%tri%the_data
            ! Find a tri-vertex inside cell
            pointOrder = 0
            loop3: do i = 1, 3
                if (BBOX(tri%P(i)%P, boxCell(:))==-1) then
                    cycle loop3
                elseif(BBOX(tri%P(i)%P, boxCell(:))==0) then
                    if (pointOrder == 0)then
                        pointOrder = -i
                    else
                        pointOrder = pointOrder - i - 5
                    endif
                else
                    pointOrder = i
                    exit loop3
                endif
            enddo loop3

            if (pointOrder == 0) then ! No point inside cell.
                derta = BGCellSize(3)/2.**c%lvl(3)
                point = FindPointInsideCell(boxCell,tri,derta)
            elseif (pointOrder < 0) then ! point on the cell face.
                if (pointOrder >= -3) then ! only one point on the cell.
                    if (pointOrder == -1) then
                        point = tri%P(1)%P + &
                                (tri%P(2)%P - tri%P(1)%P) * epsR8 + &
                                (tri%P(3)%P - tri%P(1)%P) * epsR8
                    elseif (pointOrder == -2) then
                        point = tri%P(2)%P + &
                                (tri%P(1)%P - tri%P(2)%P) * epsR8 + &
                                (tri%P(3)%P - tri%P(2)%P) * epsR8
                    else ! pointOrder == -3
                        point = tri%P(3)%P + &
                                (tri%P(1)%P - tri%P(3)%P) * epsR8 + &
                                (tri%P(2)%P - tri%P(3)%P) * epsR8
                    endif
                elseif (pointOrder >= -15) then ! two point on the cell.
                    if (pointOrder == -8) then ! point 1/2
                        point = (tri%P(1)%P+tri%P(2)%P)/2
                        point = (tri%P(3)%P-point)*epsR8+point
                    elseif (pointOrder == -9) then ! point 1/3
                        point = (tri%P(1)%P+tri%P(3)%P)/2
                        point = (tri%P(2)%P-point)*epsR8+point
                    else ! point 2/3
                        point = (tri%P(2)%P+tri%P(3)%P)/2
                        point = (tri%P(1)%P-point)*epsR8+point
                    endif
                else ! three point on the cell.
                    point = tri%P(4)%P
                endif
            else ! Have a point inside cell.
                if (pointOrder == 1) then
                    point = tri%P(1)%P + &
                            (tri%P(2)%P - tri%P(1)%P) * epsR8 + &
                            (tri%P(3)%P - tri%P(1)%P) * epsR8
                elseif (pointOrder == 2) then
                    point = tri%P(2)%P + &
                            (tri%P(1)%P - tri%P(2)%P) * epsR8 + &
                            (tri%P(3)%P - tri%P(2)%P) * epsR8
                else ! pointOrder == 3
                    point = tri%P(3)%P + &
                            (tri%P(1)%P - tri%P(3)%P) * epsR8 + &
                            (tri%P(2)%P - tri%P(3)%P) * epsR8
                endif
            endif

            n = point - c%Center
            pp=> c%CrossTri
            b = 1
            loop22: do while (ASSOCIATED(pp))
                ! if line_n cross with other tri, return.
                if (b /= nTri) then
                    tri => pp%tri%the_data
                    if (MollerTrumbore(c%Center,n,tri,.true.)) then
                        if (.not.ASSOCIATED(p%next)) print*,"error"
                        p => p%next
                        nTri = nTri + 1
                        cycle loop11
                    endif
                endif
                pp => pp%next
                b = b + 1
            enddo loop22
            if (DotProductInout(c%Center,tri)) then
                c%cross = 1
            else
                c%cross = 2
            endif
            return
        enddo loop11
        ! if (nTri<=1) stop 'error subroutine CrossCellInout'! nTri must >=1.
        ! if (nDP == nTri) then
        !     c%cross = 2 ! inside
        ! else
        !     c%cross = 1 ! outside
        ! endif
        ! return
            ! If the normal vector not cross with any other tri,
            ! this is the normal we are looking for.
            ! c%cross = 2 ! inside
            ! aia=aia+1
            ! return




                ! nDP = nDP + 1
                ! pp=> c%CrossTri
                ! b = 1
                ! n = p%tri%the_data%P(4)%P - c%Center
                ! ! n = n * 1.0000001 ! might have floot error
                ! do while (ASSOCIATED(pp))
                !     ! if normal cross with other tri, return.
                !     if (b /= nTri) then
                !         if (MollerTrumbore( c%Center, n, pp%tri%the_data, &
                !                             .true.)) then
                !             p => p%next
                !             nTri = nTri + 1
                !             cycle loop11
                !         endif
                !     endif
                !     pp => pp%next
                !     b = b + 1
                ! enddo
                ! ! If the normal vector not cross with any other tri,
                ! ! this is the normal we are looking for.
                ! c%cross = 2 ! inside
                ! aib=aib+1
                ! return
        !     p => p%next
        !     nTri = nTri + 1
        ! enddo loop11
        ! if (aia>=1.and.aib>=1)then
        ! print*,'1252346'
        ! endif
        ! if (nTri<=1) stop 'error subroutine CrossCellInout'! nTri must >=1.
        ! if (nDP == nTri) then
        !     c%cross = 2 ! inside
        ! else
        !     c%cross = 1 ! outside
        ! endif
        ! return
        endsubroutine CrossCellInout
!----------------------------------------------------------------------
        function FindPointInsideCell(box,tri,derta)
        use ModPrecision
        use ModTools, only: CrossPoint
        implicit none
        real(R8)                :: FindPointInsideCell(3)
        real(R8),INTENT(IN)     :: box(6)
        type(triangle),INTENT(IN):: tri
        real(R8),INTENT(IN)     :: derta
        real(R8)                :: V1(3), V2(3)
        real(R8)                :: D(3)
        real(R8)                :: point1(3,3) ! point cell cross with tri.
        INTEGER                 :: i
        real(R8)                :: p(4,3)
        
        i = 1
        D = (/1,0,0/)
        p(:,i) = CrossPoint(box(1:3),D*derta,tri)
        if (p(4,i)==1) then
            i = i + 1
        endif
        D = (/0,1,0/)
        p(:,i) = CrossPoint(box(1:3),D*derta,tri)
        if (p(4,i)==1) then
            i = i + 1
        endif
        D = (/0,0,1/)
        p(:,i) = CrossPoint(box(1:3),D*derta,tri)
        if (p(4,i)==1) then
            i = i + 1
            if (i == 4) goto 10
        endif
        D = (/-1,0,0/)
        p(:,i) = CrossPoint(box(4:6),D*derta,tri)
        if (p(4,i)==1) then
            i = i + 1
            if (i == 4) goto 10
        endif
        D = (/0,-1,0/)
        p(:,i) = CrossPoint(box(4:6),D*derta,tri)
        if (p(4,i)==1) then
            i = i + 1
            if (i == 4) goto 10
        endif
        D = (/0,0,-1/)
        p(:,i) = CrossPoint(box(4:6),D*derta,tri)
        if (p(4,i)==1) then
            i = i + 1
        endif

10      if (i /= 4) print*, 'error subroutine FindPointInsideCell i=', i
        FindPointInsideCell(1) = (p(1,1) + p(1,2) + p(1,3)) / 3
        FindPointInsideCell(2) = (p(2,1) + p(2,2) + p(2,3)) / 3
        FindPointInsideCell(3) = (p(3,1) + p(3,2) + p(3,3)) / 3

        endfunction FindPointInsideCell
!----------------------------------------------------------------------
        ! subroutine CCIspecial(c,a,tri,ios)
        ! use ModTools, only: PtF_normal
        ! implicit none
        ! type(typOctCell),pointer :: c
        ! integer,INTENT(IN)       :: a ! counter for loop11
        ! type(typCrossTri),pointer:: tri
        ! logical,INTENT(OUT)      :: ios ! statu, T: has changed c%cross
        ! type(typCrossTri),pointer:: ptri
        ! real(R8)                 :: n(3)
        ! integer                  :: b ! counter for loop22

        ! ios = .false.
        ! ptri=> c%CrossTri
        ! b = 1
        ! n = tri%tri%the_data%P(4)%P - c%Center
        ! n = n * 1.0000001 ! might have floot error
        ! ! n = PtF_normal(c%center,tri%tri%the_data)
        ! ! If foot point not on this tri, return.
        ! ! if (.not.MollerTrumbore(c%Center,n,tri%tri%the_data,.true.)) then
        ! !     print*,'!'
        ! !     return
        ! ! endif
        ! loop22: do while (ASSOCIATED(ptri))
        !     ! if normal cross with other tri, return.
        !     if (a /= b) then
        !     if (MollerTrumbore(c%Center,n,ptri%tri%the_data,.true.)) return
        !     endif
        !     ptri => ptri%next
        !     b = b + 1
        ! enddo loop22

        ! ! If the normal vector not cross with any other tri, this is the
        ! ! normal we are looking for.
        ! c%cross = 2 ! inside
        ! ios = .true.
        ! return

        ! ! stop 'subroutine CCIspecial special error'
        ! endsubroutine CCIspecial
! !----------------------------------------------------------------------
!         subroutine CCIspecial(c)
!         use ModTools, only: PtF_normal
!         implicit none
!         type(typOctCell),pointer :: c
!         type(typCrossTri),pointer:: p, pp
!         real(R8)                 :: n(3)
!         real(R8)                 :: dist
!         integer                  :: a, b ! counter for loop11 and loop22
!         real(R8)                 :: maxdist
!         type(triangle),pointer   :: maxtri


!         p => c%CrossTri
!         a = 1
!         maxdist = 0
!         loop11: do while (ASSOCIATED(p))
!             pp=> c%CrossTri
!             b = 1
!             n = PtF_normal(c%center,p%tri%the_data)
!             if (.not.MollerTrumbore(c%Center,n,p%tri%the_data,.true.)) return
!             loop22: do while (.true.)
!                 if (a /= b) then
!                 if (MollerTrumbore(c%Center,n,pp%tri%the_data,.true.)) &
!                     exit loop22
!                 endif
!                 pp => pp%next
!                 b = b + 1
!                 if (.not.ASSOCIATED(pp)) then !'That is my good boy!'
!                     if (MollerTrumbore(c%Center,n,p%tri%the_data,.false.))then
!                         ! Foot point on the tri.
!                         if (DotProductInout(c%Center,p%tri%the_data))then
!                             c%cross = 1
!                         else
!                             c%cross = 2
!                         endif
!                         return
!                     else
!                         dist = sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
!                         if (maxdist<dist) then
!                             maxdist = dist
!                             maxtri  =>p%tri%the_data
!                         endif
!                         exit loop22
!                     endif
!                 endif
!             enddo loop22
!             p => p%next
!             a = a + 1
!         enddo loop11

!         if (a == 1) stop 'subroutine CCIspecial no c%CrossTri'

!         if (DotProductInout(c%Center,maxtri))then
!             c%cross = 1
!         else
!             c%cross = 2
!         endif
!         return

!         ! stop 'subroutine CCIspecial special error'
!         endsubroutine CCIspecial
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
    type(typOctCell),pointer :: t

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
            t%lvl         = 0
            t%fSplitType  = 0
            t%Location    = 0
            ! t%Node        = 0
            t%Mark        = .false.
            t%Center      = (/DomainMin(1)+(i-0.5)*BGCellSize(1),   &
                              DomainMin(2)+(j-0.5)*BGCellSize(2),   &
                              DomainMin(3)+(k-0.5)*BGCellSize(3)/)
            t%U           = 0
            t%cross       = -5
            ! t%nCrossTri   = 0
            NULLIFY(t%Father,                                       &
                    t%son1, t%son2, t%son3, t%son4,                 &
                    t%son5, t%son6, t%son7, t%son8,                 &
                    t%NeighborX1, t%NeighborX2, t%NeighborY1,       &
                    t%NeighborY2, t%NeighborZ1, t%NeighborZ2)
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
    type(typOctCell),pointer :: t
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        t=>OctCell(i,j,k)
        call FindNeighbor(t)
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
    type(typOctCell),pointer :: t
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
        type(typOctCell),pointer :: c

        if(ASSOCIATED(c%son8))then
            call initCellInout2(c%son1)
            call initCellInout2(c%son2)
            call initCellInout2(c%son3)
            call initCellInout2(c%son4)
            call initCellInout2(c%son5)
            call initCellInout2(c%son6)
            call initCellInout2(c%son7)
            call initCellInout2(c%son8)
            return
        elseif(ASSOCIATED(c%son4))then
            call initCellInout2(c%son1)
            call initCellInout2(c%son2)
            call initCellInout2(c%son3)
            call initCellInout2(c%son4)
            return
        elseif(ASSOCIATED(c%son2))then
            call initCellInout2(c%son1)
            call initCellInout2(c%son2)
            return
        endif

        if (c%cross==-4 .or. c%cross==-3) call CellInout(c)

        endsubroutine initCellInout2
    endsubroutine initCellInout
!======================================================================
    subroutine initSurfaceAdapt
    use ModMesh
    use ModMeshTools
    use ModInpMesh
    use ModNeighbor
    implicit none
    type(typOctCell),pointer :: t
    integer :: i, j, k
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time
    real(R8):: p        ! Used to print precentage
    integer :: step=0   ! Counter, print progress precentage

    call CPU_TIME(tStart)
    write(6,'(1X,A,12X,A)',advance='no') 'Cell cross progress:', ''
    flush(6)
    select case (cIntersectMethod)
    case (1)    ! Ray-cast with Painting Algorithm Method
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
        call initPaintingAlgorithm2 ! Painting Algorithm Method
        call initCellInout
        cIntersectMethod = 7 ! Close the Painting Algorithm Method
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
        cIntersectMethod = 7 ! Close the Painting Algorithm Method
        call initSmoothMesh

    case (3)    ! AABB with Painting Algorithm Method
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
        call initPaintingAlgorithm2 ! Painting Algorithm Method
        cIntersectMethod = 8 ! Close the Painting Algorithm Method
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
        write(*,*) '' ! Stop write with advance='no'
        call CPU_TIME(tEnd)
        write(*,'(1X,A,F10.2)') "CrossCell time: ", tEnd-tStart
        call initFindNeighbor
        call initCellInout
        cIntersectMethod = 8 ! Close the Painting Algorithm Method
        call initSmoothMesh

    case (5)    ! AABBTraverse
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t       =>OctCell(i,j,k)
            call AABBTraverse(t)
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
        call initCellInout
        call initSmoothMesh

    case (6)    ! Ray-cast only
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t       =>OctCell(i,j,k)
            call CellCast(t)
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
        call initCellInout
        call initSmoothMesh
    end select

    contains
!----------------------------------------------------------------------
        recursive subroutine SurfaceAdapt(c)
        implicit none
        type(typOctCell),pointer :: c
        integer :: ii

        do ii=1,3
            if (c%lvl(ii) >= InitRefineLVL) return
        enddo
        if (c%cross == 1 .or. c%cross == 2 .or. c%cross == -3)then
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
    subroutine initPaintingAlgorithm2
    use ModMesh
    use ModMeshTools
    use ModInpMesh
    use ModNeighbor
    use ModPrecision
    implicit none
    type(typOctCell),pointer :: t
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
        recursive function PAFindSeed(c1) result(PA)
        ! PA: Have run/not run the subroutine PaintingAlgorithm
        implicit none
        logical               :: PA
        type(typOctCell),pointer :: c1
        logical,SAVE          :: PAused=.false.
        ! If used PaintingAlgorithm, PAused=.T.

        if (PAused) return
        if(ASSOCIATED(c1%son8))then
            PA = PAFindSeed(c1%son1)
            PA = PAFindSeed(c1%son2)
            PA = PAFindSeed(c1%son3)
            PA = PAFindSeed(c1%son4)
            PA = PAFindSeed(c1%son5)
            PA = PAFindSeed(c1%son6)
            PA = PAFindSeed(c1%son7)
            PA = PAFindSeed(c1%son8)
            return
        elseif(ASSOCIATED(c1%son4))then
            PA = PAFindSeed(c1%son1)
            PA = PAFindSeed(c1%son2)
            PA = PAFindSeed(c1%son3)
            PA = PAFindSeed(c1%son4)
            return
        elseif(ASSOCIATED(c1%son2))then
            PA = PAFindSeed(c1%son1)
            PA = PAFindSeed(c1%son2)
            return
        endif

        PA = .false.
        call CellInout(c1)
        if (c1%cross==0) then
            if (ASSOCIATED(c1%NeighborX1)) &
                call PaintingAlgorithm(c1%NeighborX1,1)
            if (ASSOCIATED(c1%NeighborX2)) &
                call PaintingAlgorithm(c1%NeighborX2,4)
            if (ASSOCIATED(c1%NeighborY1)) &
                call PaintingAlgorithm(c1%NeighborY1,2)
            if (ASSOCIATED(c1%NeighborY2)) &
                call PaintingAlgorithm(c1%NeighborY2,5)
            if (ASSOCIATED(c1%NeighborZ1)) &
                call PaintingAlgorithm(c1%NeighborZ1,3)
            if (ASSOCIATED(c1%NeighborZ2)) &
                call PaintingAlgorithm(c1%NeighborZ2,6)
            PA = .true.
            PAused = .true.
        endif
        endfunction PAFindSeed
!----------------------------------------------------------------------
        recursive subroutine PaintingAlgorithm(c,dirct)
        implicit none
        type(typOctCell),pointer :: c, cc
        integer,INTENT(IN)    :: dirct

        if(ASSOCIATED(c%son8))then
            call PaintingAlgorithm(c%son1,dirct)
            call PaintingAlgorithm(c%son2,dirct)
            call PaintingAlgorithm(c%son3,dirct)
            call PaintingAlgorithm(c%son4,dirct)
            call PaintingAlgorithm(c%son5,dirct)
            call PaintingAlgorithm(c%son6,dirct)
            call PaintingAlgorithm(c%son7,dirct)
            call PaintingAlgorithm(c%son8,dirct)
            return
        elseif(ASSOCIATED(c%son4))then
            call PaintingAlgorithm(c%son1,dirct)
            call PaintingAlgorithm(c%son2,dirct)
            call PaintingAlgorithm(c%son3,dirct)
            call PaintingAlgorithm(c%son4,dirct)
            return
        elseif(ASSOCIATED(c%son2))then
            call PaintingAlgorithm(c%son1,dirct)
            call PaintingAlgorithm(c%son2,dirct)
            return
        endif

        if (c%cross/=-4) return ! OctCell has been paintted

        if (ASSOCIATED(c%NeighborX1)) then
        if (c%NeighborX1%cross==0) then
            goto 10
        elseif (c%NeighborX1%cross/=-4) then
            if (ASSOCIATED(c%NeighborX1%son1)) then
                if (c%NeighborX1%son1%cross==0) then
                if (c%NeighborX1%son1%NeighborX2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborX1%son2%cross==0) then
                if (c%NeighborX1%son2%NeighborX2%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborX1%son3)) then
                if (c%NeighborX1%son3%cross==0) then
                if (c%NeighborX1%son3%NeighborX2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborX1%son4%cross==0) then
                if (c%NeighborX1%son4%NeighborX2%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborX1%son5)) then
                if (c%NeighborX1%son5%cross==0) then
                if (c%NeighborX1%son5%NeighborX2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborX1%son6%cross==0) then
                if (c%NeighborX1%son6%NeighborX2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborX1%son7%cross==0) then
                if (c%NeighborX1%son7%NeighborX2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborX1%son8%cross==0) then
                if (c%NeighborX1%son8%NeighborX2%nCell==c%nCell) goto 10
                endif
            endif
        endif
        endif

        if (ASSOCIATED(c%NeighborX2)) then
        if (c%NeighborX2%cross==0) then
            goto 10
        elseif (c%NeighborX2%cross/=-4) then
            if (ASSOCIATED(c%NeighborX2%son1)) then
                if (c%NeighborX2%son1%cross==0) then
                if (c%NeighborX2%son1%NeighborX1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborX2%son2%cross==0) then
                if (c%NeighborX2%son2%NeighborX1%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborX2%son3)) then
                if (c%NeighborX2%son3%cross==0) then
                if (c%NeighborX2%son3%NeighborX1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborX2%son4%cross==0) then
                if (c%NeighborX2%son4%NeighborX1%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborX2%son5)) then
                if (c%NeighborX2%son5%cross==0) then
                if (c%NeighborX2%son5%NeighborX1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborX2%son6%cross==0) then
                if (c%NeighborX2%son6%NeighborX1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborX2%son7%cross==0) then
                if (c%NeighborX2%son7%NeighborX1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborX2%son8%cross==0) then
                if (c%NeighborX2%son8%NeighborX1%nCell==c%nCell) goto 10
                endif
            endif
        endif
        endif

        if (ASSOCIATED(c%NeighborY1)) then
        if (c%NeighborY1%cross==0) then
            goto 10
        elseif (c%NeighborY1%cross/=-4) then
            if (ASSOCIATED(c%NeighborY1%son1)) then
                if (c%NeighborY1%son1%cross==0) then
                if (c%NeighborY1%son1%NeighborY2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborY1%son2%cross==0) then
                if (c%NeighborY1%son2%NeighborY2%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborY1%son3)) then
                if (c%NeighborY1%son3%cross==0) then
                if (c%NeighborY1%son3%NeighborY2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborY1%son4%cross==0) then
                if (c%NeighborY1%son4%NeighborY2%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborY1%son5)) then
                if (c%NeighborY1%son5%cross==0) then
                if (c%NeighborY1%son5%NeighborY2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborY1%son6%cross==0) then
                if (c%NeighborY1%son6%NeighborY2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborY1%son7%cross==0) then
                if (c%NeighborY1%son7%NeighborY2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborY1%son8%cross==0) then
                if (c%NeighborY1%son8%NeighborY2%nCell==c%nCell) goto 10
                endif
            endif
        endif
        endif

        if (ASSOCIATED(c%NeighborY2)) then
        if (c%NeighborY2%cross==0) then
            goto 10
        elseif (c%NeighborY2%cross/=-4) then
            if (ASSOCIATED(c%NeighborY2%son1)) then
                if (c%NeighborY2%son1%cross==0) then
                if (c%NeighborY2%son1%NeighborY1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborY2%son2%cross==0) then
                if (c%NeighborY2%son2%NeighborY1%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborY2%son3)) then
                if (c%NeighborY2%son3%cross==0) then
                if (c%NeighborY2%son3%NeighborY1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborY2%son4%cross==0) then
                if (c%NeighborY2%son4%NeighborY1%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborY2%son5)) then
                if (c%NeighborY2%son5%cross==0) then
                if (c%NeighborY2%son5%NeighborY1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborY2%son6%cross==0) then
                if (c%NeighborY2%son6%NeighborY1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborY2%son7%cross==0) then
                if (c%NeighborY2%son7%NeighborY1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborY2%son8%cross==0) then
                if (c%NeighborY2%son8%NeighborY1%nCell==c%nCell) goto 10
                endif
            endif
        endif
        endif

        if (ASSOCIATED(c%NeighborZ1)) then
        if (c%NeighborZ1%cross==0) then
            goto 10
        elseif (c%NeighborZ1%cross/=-4) then
            if (ASSOCIATED(c%NeighborZ1%son1)) then
                if (c%NeighborZ1%son1%cross==0) then
                if (c%NeighborZ1%son1%NeighborZ2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborZ1%son2%cross==0) then
                if (c%NeighborZ1%son2%NeighborZ2%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborZ1%son3)) then
                if (c%NeighborZ1%son3%cross==0) then
                if (c%NeighborZ1%son3%NeighborZ2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborZ1%son4%cross==0) then
                if (c%NeighborZ1%son4%NeighborZ2%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborZ1%son5)) then
                if (c%NeighborZ1%son5%cross==0) then
                if (c%NeighborZ1%son5%NeighborZ2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborZ1%son6%cross==0) then
                if (c%NeighborZ1%son6%NeighborZ2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborZ1%son7%cross==0) then
                if (c%NeighborZ1%son7%NeighborZ2%nCell==c%nCell) goto 10
                endif
                if (c%NeighborZ1%son8%cross==0) then
                if (c%NeighborZ1%son8%NeighborZ2%nCell==c%nCell) goto 10
                endif
            endif
        endif
        endif

        if (ASSOCIATED(c%NeighborZ2)) then
        if (c%NeighborZ2%cross==0) then
            goto 10
        elseif (c%NeighborZ2%cross/=-4) then
            if (ASSOCIATED(c%NeighborZ2%son1)) then
                if (c%NeighborZ2%son1%cross==0) then
                if (c%NeighborZ2%son1%NeighborZ1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborZ2%son2%cross==0) then
                if (c%NeighborZ2%son2%NeighborZ1%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborZ2%son3)) then
                if (c%NeighborZ2%son3%cross==0) then
                if (c%NeighborZ2%son3%NeighborZ1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborZ2%son4%cross==0) then
                if (c%NeighborZ2%son4%NeighborZ1%nCell==c%nCell) goto 10
                endif
            endif
            if (ASSOCIATED(c%NeighborZ2%son5)) then
                if (c%NeighborZ2%son5%cross==0) then
                if (c%NeighborZ2%son5%NeighborZ1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborZ2%son6%cross==0) then
                if (c%NeighborZ2%son6%NeighborZ1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborZ2%son7%cross==0) then
                if (c%NeighborZ2%son7%NeighborZ1%nCell==c%nCell) goto 10
                endif
                if (c%NeighborZ2%son8%cross==0) then
                if (c%NeighborZ2%son8%NeighborZ1%nCell==c%nCell) goto 10
                endif
            endif
        endif
        endif
        return

10      c%cross=0
        if (dirct /= 4) then
            cc => c%NeighborX1
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,1)
        endif
        if (dirct /= 1) then
            cc => c%NeighborX2
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,4)
        endif
        if (dirct /= 5) then
            cc => c%NeighborY1
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,2)
        endif
        if (dirct /= 2) then
            cc => c%NeighborY2
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,5)
        endif
        if (dirct /= 6) then
            cc => c%NeighborZ1
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,3)
        endif
        if (dirct /= 3) then
            cc => c%NeighborZ2
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,6)
        endif
        endsubroutine PaintingAlgorithm
!----------------------------------------------------------------------
        recursive subroutine PAMarkCrossCell(c)
        implicit none
        type(typOctCell),pointer :: c

        if (c%cross/=-3) return
        if(ASSOCIATED(c%son8))then
            call PAMarkCrossCell(c%son1)
            call PAMarkCrossCell(c%son2)
            call PAMarkCrossCell(c%son3)
            call PAMarkCrossCell(c%son4)
            call PAMarkCrossCell(c%son5)
            call PAMarkCrossCell(c%son6)
            call PAMarkCrossCell(c%son7)
            call PAMarkCrossCell(c%son8)
            return
        elseif(ASSOCIATED(c%son4))then
            call PAMarkCrossCell(c%son1)
            call PAMarkCrossCell(c%son2)
            call PAMarkCrossCell(c%son3)
            call PAMarkCrossCell(c%son4)
            return
        elseif(ASSOCIATED(c%son2))then
            call PAMarkCrossCell(c%son1)
            call PAMarkCrossCell(c%son2)
            return
        endif

        call CellInout(c)

        endsubroutine PAMarkCrossCell
    endsubroutine initPaintingAlgorithm2
!======================================================================
    subroutine initPaintingAlgorithm
    use ModMesh
    use ModMeshTools
    use ModInpMesh
    use ModNeighbor
    use ModPrecision
    implicit none
    type(typOctCell),pointer :: t
    integer :: ii, jj, kk
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    write(6,'(1X,A,12X,A)',advance='no') 'PaintingAlgorithm progress:', ''
    flush(6)

    loop:  do kk = 1, nCell(3)
    do jj = 1, nCell(2)
    do ii = 1, nCell(1)
        t       =>OctCell(ii,jj,kk)
        if (PAFindSeed(t)) exit loop
    enddo
    enddo
    enddo loop

    write(*,*) '' ! Stop write with advance='no'
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "Painting Algorithm time: ", tEnd-tStart

    contains
!----------------------------------------------------------------------
        recursive function PAFindSeed(c1) result(PA)
        ! PA: Have run/not run the subroutine PaintingAlgorithm
        implicit none
        logical               :: PA
        type(typOctCell),pointer :: c1
        logical,SAVE          :: PAused=.false.
        ! If used PaintingAlgorithm, PAused=.T.

        if (PAused) return
        if(ASSOCIATED(c1%son8))then
            PA = PAFindSeed(c1%son1)
            PA = PAFindSeed(c1%son2)
            PA = PAFindSeed(c1%son3)
            PA = PAFindSeed(c1%son4)
            PA = PAFindSeed(c1%son5)
            PA = PAFindSeed(c1%son6)
            PA = PAFindSeed(c1%son7)
            PA = PAFindSeed(c1%son8)
            return
        elseif(ASSOCIATED(c1%son4))then
            PA = PAFindSeed(c1%son1)
            PA = PAFindSeed(c1%son2)
            PA = PAFindSeed(c1%son3)
            PA = PAFindSeed(c1%son4)
            return
        elseif(ASSOCIATED(c1%son2))then
            PA = PAFindSeed(c1%son1)
            PA = PAFindSeed(c1%son2)
            return
        endif

        PA = .false.
        call CellInout(c1)
        if (c1%cross==0) then
            if (ASSOCIATED(c1%NeighborX1)) &
                call PaintingAlgorithm(c1%NeighborX1,1,.true.)
            if (ASSOCIATED(c1%NeighborX2)) &
                call PaintingAlgorithm(c1%NeighborX2,4,.true.)
            if (ASSOCIATED(c1%NeighborY1)) &
                call PaintingAlgorithm(c1%NeighborY1,2,.true.)
            if (ASSOCIATED(c1%NeighborY2)) &
                call PaintingAlgorithm(c1%NeighborY2,5,.true.)
            if (ASSOCIATED(c1%NeighborZ1)) &
                call PaintingAlgorithm(c1%NeighborZ1,3,.true.)
            if (ASSOCIATED(c1%NeighborZ2)) &
                call PaintingAlgorithm(c1%NeighborZ2,6,.true.)
            PA = .true.
            PAused = .true.
        endif
        endfunction PAFindSeed
!----------------------------------------------------------------------
        recursive subroutine PaintingAlgorithm(c,dirct,iosIn)
        implicit none
        type(typOctCell),pointer :: c, cc
        integer               :: dirct
        LOGICAL,INTENT(IN)    :: iosIn  ! = F , come from cross=-3 OctCell
                                        ! = T , come from cross= 0 OctCell
        LOGICAL               :: iosOut
        integer,SAVE          :: progress=0 ! counter, print precentage

        iosOut=iosIn
        if (c%cross==-3) iosOut = .false.
        if(ASSOCIATED(c%son8))then
            call PaintingAlgorithm(c%son1,dirct,iosOut)
            call PaintingAlgorithm(c%son2,dirct,iosOut)
            call PaintingAlgorithm(c%son3,dirct,iosOut)
            call PaintingAlgorithm(c%son4,dirct,iosOut)
            call PaintingAlgorithm(c%son5,dirct,iosOut)
            call PaintingAlgorithm(c%son6,dirct,iosOut)
            call PaintingAlgorithm(c%son7,dirct,iosOut)
            call PaintingAlgorithm(c%son8,dirct,iosOut)
            return
        elseif(ASSOCIATED(c%son4))then
            call PaintingAlgorithm(c%son1,dirct,iosOut)
            call PaintingAlgorithm(c%son2,dirct,iosOut)
            call PaintingAlgorithm(c%son3,dirct,iosOut)
            call PaintingAlgorithm(c%son4,dirct,iosOut)
            return
        elseif(ASSOCIATED(c%son2))then
            call PaintingAlgorithm(c%son1,dirct,iosOut)
            call PaintingAlgorithm(c%son2,dirct,iosOut)
            return
        endif

        if (c%cross/=-3 .and. c%cross/=-4) return ! Cell has been paintted

        if (iosIn) then
            if (c%cross==-4)then ! This OctCell is a outside OctCell too.
                c%cross=0
                iosOut = .True.
            elseif (c%cross==-3) then
                progress = progress + 1 ! Once call cellinout, progress++
                call CellInout(c)
                iosOut = .false.
                !if (c%cross==2) return
            endif
        else
            progress = progress + 1
            call CellInout(c)
            if (c%cross==3) then
                iosOut = .false.
                return ! Inside OctCell, return.
            elseif (c%cross==2) then
                iosOut = .false.
                !return
            elseif (c%cross==1) then
                iosOut = .false.
            elseif (c%cross==0) then
                iosOut = .true.
            endif
        endif

        ! counter, print progress precentage
        if (mod(progress,100)==0) then
            write(6,'(A,F5.1,A)',advance='no') &
                    '\b\b\b\b\b\b', progress/real(nCells,R8)*100, '%'
            flush(6)
        endif

        if (dirct /= 4) then
            cc => c%NeighborX1
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,1,iosOut)
        endif
        if (dirct /= 1) then
            cc => c%NeighborX2
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,4,iosOut)
        endif
        if (dirct /= 5) then
            cc => c%NeighborY1
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,2,iosOut)
        endif
        if (dirct /= 2) then
            cc => c%NeighborY2
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,5,iosOut)
        endif
        if (dirct /= 6) then
            cc => c%NeighborZ1
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,3,iosOut)
        endif
        if (dirct /= 3) then
            cc => c%NeighborZ2
            if (ASSOCIATED(cc)) call PaintingAlgorithm(cc,6,iosOut)
        endif
        endsubroutine PaintingAlgorithm
    endsubroutine initPaintingAlgorithm
!======================================================================
    subroutine initSmoothMesh
    use ModPrecision
    use ModMesh
    use ModTypDef
    use ModMeshTools
    use ModInpMesh
    implicit none
    type(typOctCell),POINTER :: t
    integer               :: i, j, k, ii
    logical               :: iost
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    do ii=1, InitRefineLVL-1

        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t=>OctCell(i, j, k)
            call PreSmoothMesh(t) ! Mark neighbor OctCell if need refine.
        enddo
        enddo
        enddo

        iost = .false.
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t=>OctCell(i, j, k)
            call SmoothMesh(t,iost) ! Do smooth.
        enddo
        enddo
        enddo

        if (.not.iost) exit

        ! call initFindNeighbor

    enddo
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "Smooth Mesh time: ", tEnd-tStart
    contains
!----------------------------------------------------------------------
        recursive subroutine PreSmoothMesh(c)
        implicit none
        type(typOctCell),POINTER :: c,cx1,cx2,cy1,cy2,cz1,cz2
        ! mark(6)  1 x refine; 2 y refine; 3 z refine; 
        !          4 x coarse; 5 y coarse; 6 z coarse; 
        if(ASSOCIATED(c%son8))then
            call PreSmoothMesh(c%son1)
            call PreSmoothMesh(c%son2)
            call PreSmoothMesh(c%son3)
            call PreSmoothMesh(c%son4)
            call PreSmoothMesh(c%son5)
            call PreSmoothMesh(c%son6)
            call PreSmoothMesh(c%son7)
            call PreSmoothMesh(c%son8)
            return
        elseif(ASSOCIATED(c%son4))then
            call PreSmoothMesh(c%son1)
            call PreSmoothMesh(c%son2)
            call PreSmoothMesh(c%son3)
            call PreSmoothMesh(c%son4)
            return
        elseif(ASSOCIATED(c%son2))then
            call PreSmoothMesh(c%son1)
            call PreSmoothMesh(c%son2)
            return
        endif
        
        if (ASSOCIATED(c%NeighborX1)) then
            if (c%NeighborX1%lvl(1)+1<c%lvl(1)) c%NeighborX1%mark(1)=.true.
        endif
        if (ASSOCIATED(c%NeighborX2)) then
            if (c%NeighborX2%lvl(1)+1<c%lvl(1)) c%NeighborX2%mark(1)=.true.
        endif
        if (ASSOCIATED(c%NeighborY1)) then
            if (c%NeighborY1%lvl(2)+1<c%lvl(2)) c%NeighborY1%mark(2)=.true.
        endif
        if (ASSOCIATED(c%NeighborY2)) then
            if (c%NeighborY2%lvl(2)+1<c%lvl(2)) c%NeighborY2%mark(2)=.true.
        endif
        if (ASSOCIATED(c%NeighborZ1)) then
            if (c%NeighborZ1%lvl(3)+1<c%lvl(3)) c%NeighborZ1%mark(3)=.true.
        endif
        if (ASSOCIATED(c%NeighborZ2)) then
            if (c%NeighborZ2%lvl(3)+1<c%lvl(3)) c%NeighborZ2%mark(3)=.true.
        endif
        ! Hole OctCell
        ! if (.not.mark(1)) then
        !     if (ASSOCIATED(c%NeighborX1).and.ASSOCIATED(c%NeighborX2))then
        !         if (c%NeighborX1%lvl(1)>c%lvl(1) .and. &
        !             c%NeighborX2%lvl(1)>c%lvl(1)) then
        !             mark(1)=.true.
        !         elseif (c%NeighborX1%lvl(1)<c%lvl(1) .and. &
        !             c%NeighborX2%lvl(1)<c%lvl(1)) then
        !             mark(4)=.true.
        !         endif
        !     endif
        ! endif
        
        endsubroutine PreSmoothMesh
!----------------------------------------------------------------------
        recursive subroutine SmoothMesh(c,iosout)
        implicit none
        type(typOctCell),POINTER :: c
        logical,INTENT(INOUT) :: iosout
        logical               :: ios
        ! mark(6)  1 x refine; 2 y refine; 3 z refine;
        !          4 x coarse; 5 y coarse; 6 z coarse;

        if (maxval(c%lvl)>=InitRefineLVL) return

        if(ASSOCIATED(c%son8))then
            call SmoothMesh(c%son1,ios)
            call SmoothMesh(c%son2,ios)
            call SmoothMesh(c%son3,ios)
            call SmoothMesh(c%son4,ios)
            call SmoothMesh(c%son5,ios)
            call SmoothMesh(c%son6,ios)
            call SmoothMesh(c%son7,ios)
            call SmoothMesh(c%son8,ios)
            return
        elseif(ASSOCIATED(c%son4))then
            call SmoothMesh(c%son1,ios)
            call SmoothMesh(c%son2,ios)
            call SmoothMesh(c%son3,ios)
            call SmoothMesh(c%son4,ios)
            return
        elseif(ASSOCIATED(c%son2))then
            call SmoothMesh(c%son1,ios)
            call SmoothMesh(c%son2,ios)
            return
        endif

        if (c%cross/=0) return
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
            c%mark=.false.
            iosout = .true.
            call UpdateNeighbor(c%NeighborX1,2)
            call UpdateNeighbor(c%NeighborX2,1)
            call UpdateNeighbor(c%NeighborY1,4)
            call UpdateNeighbor(c%NeighborY2,3)
            call UpdateNeighbor(c%NeighborZ1,6)
            call UpdateNeighbor(c%NeighborZ2,5)
            call FindNeighbor(c%son1)
            call FindNeighbor(c%son2)
            call FindNeighbor(c%son3)
            call FindNeighbor(c%son4)
            call FindNeighbor(c%son5)
            call FindNeighbor(c%son6)
            call FindNeighbor(c%son7)
            call FindNeighbor(c%son8)
        endif
        ! Initial OctCell mark.

        ! Coarse
        ! TODO
        ! elseif ( c%mark(4) .and. c%mark(5) .and. c%mark(6) ) then
        !     call DeletCell(c,0)
        ! elseif ( c%mark(4) .and. .not.c%mark(5) .and. .not.c%mark(6) ) then
        !     call DeletCell(c,1)
        ! elseif ( .not.c%mark(4) .and. c%mark(5) .and. .not.c%mark(6) ) then
        !     call DeletCell(c,2)
        ! elseif ( .not.c%mark(4) .and. .not.c%mark(5) .and. c%mark(6) ) then
        !     call DeletCell(c,3)
        ! elseif ( c%mark(4) .and. c%mark(5) .and. .not.c%mark(6) ) then
        !     call DeletCell(c,4)
        ! elseif ( c%mark(4) .and. .not.c%mark(5) .and. c%mark(6) ) then
        !     call DeletCell(c,5)
        ! elseif ( .not.c%mark(4) .and. c%mark(5) .and. c%mark(6) ) then
        !     call DeletCell(c,6)
        endsubroutine SmoothMesh
!----------------------------------------------------------------------
    endsubroutine initSmoothMesh
!======================================================================
!======================================================================
!----------------------------------------------------------------------
