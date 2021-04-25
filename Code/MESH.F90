!======================================================================
    module ModMesh
        use ModTypDef
        implicit none
        ! Mesh
        integer :: nBGCells     ! Number of the background Cells.
        integer :: nCells     ! Number of total Cells.
        real(R8):: BGCellSize(3)    ! Step size for background cell.
        type(octCell),pointer :: Cell(:,:,:)
        ! Geometry
        ! integer :: nGeoPoints   ! Number of the geometry points.
        ! integer :: nGeoFaces    ! Number of the geometry faces.
        ! real(R8),ALLOCATABLE :: Geometry(:,:)
        ! integer ,ALLOCATABLE :: GeoFace(:,:)
        ! type(typPoint)       :: GeoBBOX(2) ! Bounding box of geometry.
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
        type(octCell),pointer :: c

        select case (cIntersectMethod)
        case (1)
            if (CellCast(c)) then
                c%cross = -3
            else
                c%cross = -4
            endif
        case (2)
            if (CellCast(c)) then
                c%cross = -3
                c%cross = CellInout(c)
            else
                c%cross = -4
                c%cross = CellInout(c)
            endif
        case (3)
            if (AABB(c)) then
                c%cross = -3
            else
                c%cross = -4
            endif
        case (4)
            if (AABB(c)) then
                c%cross = -3
                c%cross = CellInout(c)
            else
                c%cross = -4
                c%cross = CellInout(c)
            endif
        case (5)
            if (AABBTraverse(c)) then
                c%cross = -3
                c%cross = CellInout(c)
            else
                c%cross = -4
                c%cross = CellInout(c)
            endif
        end select
        endsubroutine initCellCross
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
        logical function CellCast(c)    ! return -3, 0, 3
        use ModKDTree
        use ModInpGlobal
        implicit none
        type(octCell),pointer :: c
        type(typPoint)        :: p(8)
        real(R8):: x, y, z, dx, dy, dz
        integer :: i, ii, Pintersect

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
            do ii = 1,nGeometry
            ! Quick Bounding-BOX identify
            if (BBOX(p(i)%P,KDTree(ii)%root%box)) then ! inside
                Pintersect=Pintersect+RayCast(p(i)%P,KDTree(ii)%root)
            endif
            enddo
            if (Pintersect/=0.and.Pintersect<i) then
                CellCast = .true.  ! Intersect
                return
            endif
        enddo
                CellCast = .false.  ! Not intersect
        endfunction CellCast
!----------------------------------------------------------------------
        integer function CellInout(c)
        use ModKDTree
        use ModInpGlobal
        implicit none
        type(octCell),pointer :: c
        type(typPoint)        :: point
        integer :: ii, Pintersect

        Pintersect=0
        point%P=(/c%Center(1),c%Center(2),c%Center(3)/)
        do ii = 1,nGeometry
            ! Quick Bounding-BOX identify
            if (BBOX(point%P,KDTree(ii)%root%box)) then
            Pintersect=Pintersect+RayCast(point%P,KDTree(ii)%root)
            endif
        enddo
        if (Pintersect == 0) then ! inside
            if (c%cross == -4) then
                CellInout = 0 ; return
            else
                CellInout = 1 ; return
            endif
        else
            if (c%cross == -4) then
                CellInout = 3 ; return
            else
                CellInout = 2 ; return
            endif
        endif
        endfunction CellInout
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
                ! Quick BBOX identify
                do jj = 1,3
                    triFace(jj)%P = body(ng)%se3d(j)%P(jj)%P
                enddo
                select case (i)
                case (1)
                if ((point(2)<min(triFace(1)%P(2), &
                                  triFace(2)%P(2), &
                                  triFace(3)%P(2))  .or.  &
                     point(2)>max(triFace(1)%P(2), &
                                  triFace(2)%P(2), &
                                  triFace(3)%P(2))) .and. &
                    (point(3)<min(triFace(1)%P(3), &
                                  triFace(2)%P(3), &
                                  triFace(3)%P(3))  .or.  &
                     point(3)>min(triFace(1)%P(3), &
                                  triFace(2)%P(3), &
                                  triFace(3)%P(3)))) cycle ! outside
                case (2)
                if ((point(1)<min(triFace(1)%P(1), &
                                  triFace(2)%P(1), &
                                  triFace(3)%P(1))  .or.  &
                     point(1)>max(triFace(1)%P(1), &
                                  triFace(2)%P(1), &
                                  triFace(3)%P(1))) .and. &
                    (point(3)<min(triFace(1)%P(3), &
                                  triFace(2)%P(3), &
                                  triFace(3)%P(3))  .or.  &
                     point(3)>min(triFace(1)%P(3), &
                                  triFace(2)%P(3), &
                                  triFace(3)%P(3)))) cycle ! outside
                case (3)
                if ((point(1)<min(triFace(1)%P(1), &
                                  triFace(2)%P(1), &
                                  triFace(3)%P(1))  .or.  &
                     point(1)>max(triFace(1)%P(1), &
                                  triFace(2)%P(1), &
                                  triFace(3)%P(1))) .and. &
                    (point(2)<min(triFace(1)%P(2), &
                                  triFace(2)%P(2), &
                                  triFace(3)%P(2))  .or.  &
                     point(2)>min(triFace(1)%P(2), &
                                  triFace(2)%P(2), &
                                  triFace(3)%P(2)))) cycle ! outside
                end select
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
            RayCast=1; return   ! inside
        else
            RayCast=0; return   ! outside
        endif
        endfunction RayCast
!----------------------------------------------------------------------
        Logical function AABB(c)!1=intersect,0=no intersect
            use ModKDTree
            use ModInpGlobal
            use ModGeometry
            implicit none
            type(octCell),pointer :: c
            type(triangle)        :: tri
            real(R8)              :: boxCell(6)
            type(KDT_node),pointer:: node
            logical               :: aaa
            type(typKDTtree), pointer   :: tp => null()
            integer                     :: ng, i
            
            aaa=.false.
            boxCell(1)=c%center(1)-BGCellSize(1)/2**(c%lvl(1)+1)
            boxCell(2)=c%center(2)-BGCellSize(2)/2**(c%lvl(2)+1)
            boxCell(3)=c%center(3)-BGCellSize(3)/2**(c%lvl(3)+1)
            boxCell(4)=c%center(1)+BGCellSize(1)/2**(c%lvl(1)+1)
            boxCell(5)=c%center(2)+BGCellSize(2)/2**(c%lvl(2)+1)
            boxCell(6)=c%center(3)+BGCellSize(3)/2**(c%lvl(3)+1)
            
            tp => kdtree(1)
            node=>tp%root
            if( boxCell(1)>node%box(4).or.boxCell(2)>node%box(5).or. &
                boxCell(3)>node%box(6).or.boxCell(4)<node%box(1).or. &
                boxCell(5)<node%box(2).or.boxCell(6)<node%box(3))then
                 AABB=.false.
            else
                call KdFindTri(c,boxCell,node,aaa)
                if(aaa)then
                    AABB=.true.
                else
                    AABB=.false.
                endif        
            endif
            return
        end function AABB
!----------------------------------------------------------------------
        Logical function AABBTraverse(c)!1=intersect,0=no intersect
            use ModKDTree
            use ModInpGlobal
            use ModGeometry
            implicit none
            type(octCell),pointer :: c
            type(triangle)        :: tri
            integer                     :: ng, i
          
            Loop1:do ng=1,nGeometry
                do i=1, body(ng)%nse
                    tri= body(ng)%se3d(i)
                    AABBTraverse=TriBoxOverlap (c,tri)
                    if(AABBTraverse) exit Loop1
                enddo
            enddo Loop1
            return         
        end function AABBTraverse
!----------------------------------------------------------------------
        recursive subroutine KdFindTri(c,boxCell,node,aaa)
            use ModKDTree
            implicit none
            
            real(R8),INTENT(IN)                     :: boxCell(6)
            integer                                 :: i,split,splitmax,a
            type(KDT_node), pointer,INTENT(IN)      :: node
            type(KDT_node), pointer                 :: leftside,rightside
            type(triangle)                          :: triright,trileft,tri
            logical, INTENT(OUT)                    :: aaa
            type(octCell),pointer                   :: c
            
            tri = node%the_data
            aaa = TriBoxOverlap(c,tri)
            !write(*,*) node%level
            if(aaa)then
                return
            endif
            
            if(associated(node%right))then
            if( boxCell(1)<=node%right%box(4).and.boxCell(4)>=node%right%box(1).and.&
                boxCell(2)<=node%right%box(5).and.boxCell(5)>=node%right%box(2).and.&
                boxCell(3)<=node%right%box(6).and.boxCell(6)>=node%right%box(3))then
                    call KdFindTri(c,boxCell,node%right,aaa)
                    if(aaa)then
                      return
                    endif
            endif
            endif
            
            if(associated(node%left))then
            if( boxCell(1)<=node%left%box(4).and.boxCell(4)>=node%left%box(1).and.&
                boxCell(2)<=node%left%box(5).and.boxCell(5)>=node%left%box(2).and.&
                boxCell(3)<=node%left%box(6).and.boxCell(6)>=node%left%box(3))then
                    call KdFindTri(c,boxCell,node%left,aaa)
                    if(aaa)then
                      return
                    endif
            endif
            endif
        end subroutine KdFindTri    
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
            
            type(octCell),pointer :: c
            type(triangle)        :: tri
            real(R8)              :: v0(3), v1(3), v2(3)
            real(R8)              :: boxcenter(3), normal(3)
            real(R8)              :: e0(3), e1(3), e2(3)
            real(R8)              :: boxhalfsize(3)
            real(R8)              :: min, max
            real(R8)              :: rad, d, dd, fex, fey, fez, a
            real(R8)              :: X01P0, X01P2
            real(R8)              :: Y02P0, Y02P2
            real(R8)              :: Z12P2,Z12P1
            real(R8)              :: Z0P0, Z0P1
            real(R8)              :: X2P0, X2P1
            real(R8)              :: Y1P0, Y1P1
        
            boxcenter(1)=c%center(1)
            boxcenter(2)=c%center(2)
            boxcenter(3)=c%center(3)
            boxhalfsize(1)= BGCellSize(1)/2**(c%lvl(1)+1)
            boxhalfsize(2)= BGCellSize(2)/2**(c%lvl(2)+1)
            boxhalfsize(3)= BGCellSize(3)/2**(c%lvl(3)+1)
            
            
            v0(:)=tri%p(1)%P-boxcenter
            v1(:)=tri%p(2)%P-boxcenter
            v2(:)=tri%p(3)%P-boxcenter
            
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
                min=X01P0
                max=X01P2
            else
                min=X01P2
                max=X01P0
            endif
            rad=fez*boxhalfsize(2)+fey*boxhalfsize(3)
            if(min>rad.or.max<-rad)then
                TriBoxOverlap=.false.
                return 
            endif
            !AxisTestY02
            Y02P0=-e0(3)*v0(1)+e0(1)*v0(3)
            Y02P2=-e0(3)*v2(1)+e0(1)*v2(3)
            if(Y02P0<Y02P2)then
                min=Y02P0
                max=Y02P2
            else
                min=Y02P2
                max=Y02P0
            endif
            rad=fez*boxhalfsize(1)+fex*boxhalfsize(3)
            if(min>rad.or.max<-rad)then
                TriBoxOverlap=.false.
                return 
            endif
            !AxisTestZ12
            Z12P1=e0(2)*v1(1)-e0(1)*v1(2)
            Z12P2=e0(2)*v2(1)-e0(1)*v2(2)
            if(Z12P2<Z12P1)then
                min=Z12P2
                max=Z12P1
            else
                min=Z12P1
                max=Z12P2
            endif
            rad=fey*boxhalfsize(1)+fex*boxhalfsize(2)
            if(min>rad.or.max<-rad)then
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
                min=X01P0
                max=X01P2
            else
                min=X01P2
                max=X01P0
            endif
            rad=fez*boxhalfsize(2)+fey*boxhalfsize(3)
            if(min>rad.or.max<-rad)then
                TriBoxOverlap=.false.
                return 
            endif
            !AxisTestY02
            Y02P0=-e1(3)*v0(1)+e1(1)*v0(3)
            Y02P2=-e1(3)*v2(1)+e1(1)*v2(3)
            if(Y02P0<Y02P2)then
                min=Y02P0
                max=Y02P2
            else
                min=Y02P2
                max=Y02P0
            endif
            rad=fez*boxhalfsize(1)+fex*boxhalfsize(3)
            if(min>rad.or.max<-rad)then
                TriBoxOverlap=.false.
                return 
            endif
            !AxisTestZ0
            Z0P0=e1(2)*v0(1)-e1(1)*v0(2)
            Z0P1=e1(2)*v1(1)-e1(1)*v1(2)
            if(Z0P0<Z0P1)then
                min=Z0P0
                max=Z0P1
            else
                min=Z0P1
                max=Z0P0
            endif
            rad=fey*boxhalfsize(1)+fex*boxhalfsize(2)
            if(min>rad.or.max<-rad)then
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
                min=X2P0
                max=X2P1
            else
                min=X2P1
                max=X2P0
            endif
            rad=fez*boxhalfsize(2)+fey*boxhalfsize(3)
            if(min>rad.or.max<-rad)then
                TriBoxOverlap=.false.
                return
            endif
            
            !AxisTestY1
            Y1P0=-e2(3)*v0(1)+e2(1)*v0(3)
            Y1P1=-e2(3)*v1(1)+e2(1)*v1(3)
            if(Y1P0<Y1P1)then
                min=Y1P0
                max=Y1P1
            else
                min=Y1P1
                max=Y1P0
            endif
            rad=fez*boxhalfsize(1)+fex*boxhalfsize(3)
            if(min>rad.or.max<-rad)then
                TriBoxOverlap=.false.
                return
            endif
            !AxisTestZ12
            Z12P1=e2(2)*v1(1)-e2(1)*v1(2)
            Z12P2=e2(2)*v2(1)-e2(1)*v2(2)
            if(Z12P2<Z12P1)then
                min=Z12P2
                max=Z12P1
            else
                min=Z12P1
                max=Z12P2
            endif
            rad=fey*boxhalfsize(1)+fex*boxhalfsize(2)
            if(min>rad.or.max<-rad)then
                TriBoxOverlap=.false.
                return
            endif
            
            
!/* Bullet 1: */
!/* first test overlap in the {x,y,z}-directions */
!/* find min, max of the triangle each direction, and test for overlap in */
!/* that direction -- this is equivalent to testing a minimal AABB around */
!/*  the triangle against the AABB */
            
            !/* test in X-direction */ 
            min=FindMin(v0(1),v1(1),v2(1))
            max=FindMax(v0(1),v1(1),v2(1))
            if(min>boxhalfsize(1).or.max<-boxhalfsize(1))then
                TriBoxOverlap=.false.
                return
            endif
                        
            !/* test in Y-direction */ 
            min=FindMin(v0(2),v1(2),v2(2))
            max=FindMax(v0(2),v1(2),v2(2))
            if(min>boxhalfsize(2).or.max<-boxhalfsize(2))then
                TriBoxOverlap=.false.
                return
            endif
            
            !/* test in Z-direction */ 
            min=FindMin(v0(3),v1(3),v2(3))
            max=FindMax(v0(3),v1(3),v2(3))
            if(min>boxhalfsize(3).or.max<-boxhalfsize(3)) then
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
          function Sub(a,b)
            use ModPrecision
            implicit none
            real(R8)           :: SUB(3)
            real(R8),INTENT(IN):: a(3), b(3)
            SUB =              [a(1)-b(1), &
                                a(2)-b(2), &
                                a(3)-b(3)]
        endfunction Sub
!----------------------------------------------------------------------        
        Real(R8)function FindMin(x0,x1,x2)
            use ModPrecision
            implicit none
            real(R8),INTENT(IN)  :: x0,x1,x2
                       
            FindMin = x0
            
            if (x1 < FindMin) then
                FindMin = x1
            endif 
            if (x2 < FindMin) then
                FindMin = x2
            endif
        end function FindMin
       Real(R8) function FindMax(x0,x1,x2)
            use ModPrecision
            implicit none
            real(R8),INTENT(IN)  ::x0,x1,x2
            
            FindMax = x0
            if (x1 > FindMax) then
                FindMax = x1
            endif 
            if (x2 > FindMax) then
                FindMax = x2
            endif
        end function FindMax
        
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
        subroutine NewCell(c,split)
        use ModInpGlobal
        use ModNeighbor
        ! Called by: many
        ! Calls    : initNeighbor
        implicit none
        integer(I4),INTENT(IN):: split
        integer(I4)           :: spl
        type(octCell),pointer :: c
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
            c%son1%Node=0
            c%son2%Node=0
            c%son3%Node=0
            c%son4%Node=0
            c%son5%Node=0
            c%son6%Node=0
            c%son7%Node=0
            c%son8%Node=0
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
            c%son1%Father=>c
            c%son2%Father=>c
            c%son3%Father=>c
            c%son4%Father=>c
            c%son5%Father=>c
            c%son6%Father=>c
            c%son7%Father=>c
            c%son8%Father=>c
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
            c%son1%Node=0
            c%son2%Node=0
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
            if (c%cross==1 .or. c%cross==2 .or. c%cross==-3) then
                call initCellCross(c%son1)
                call initCellCross(c%son2)
            else
                c%son1%cross=c%cross
                c%son2%cross=c%cross
            endif
            c%son1%Father=>c
            c%son2%Father=>c
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
            c%son1%Node=0
            c%son2%Node=0
            c%son3%Node=0
            c%son4%Node=0
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
            c%son1%Father=>c
            c%son2%Father=>c
            c%son3%Father=>c
            c%son4%Father=>c
            call NullifyCell(c%son1)
            call NullifyCell(c%son2)
            call NullifyCell(c%son3)
            call NullifyCell(c%son4)
            return
        end select
        endsubroutine NewCell
!----------------------------------------------------------------------
        subroutine DeletCell(c,split)
        implicit none
        type(octCell),pointer :: c
        integer(I4),INTENT(IN):: split
        endsubroutine DeletCell
!----------------------------------------------------------------------
        subroutine NullifyCell(c)
        implicit none
        type(octCell),pointer :: c
        NULLIFY(c%son1, c%son2, c%son3, c%son4,                 &
                c%son5, c%son6, c%son7, c%son8,                 &
                c%NeighborX1, c%NeighborX2, c%NeighborY1,       &
                c%NeighborY2, c%NeighborZ1, c%NeighborZ2)
        endsubroutine NullifyCell
!----------------------------------------------------------------------
        recursive subroutine InitialCellMark(c)
        implicit none
        type(octCell),pointer :: c
        if(ASSOCIATED(c%son8))then
            call InitialCellMark(c%son1)
            call InitialCellMark(c%son2)
            call InitialCellMark(c%son3)
            call InitialCellMark(c%son4)
            call InitialCellMark(c%son5)
            call InitialCellMark(c%son6)
            call InitialCellMark(c%son7)
            call InitialCellMark(c%son8)
            return
        elseif(ASSOCIATED(c%son4))then
            call InitialCellMark(c%son1)
            call InitialCellMark(c%son2)
            call InitialCellMark(c%son3)
            call InitialCellMark(c%son4)
            return
        elseif(ASSOCIATED(c%son2))then
            call InitialCellMark(c%son1)
            call InitialCellMark(c%son2)
            return
        endif
        c%mark=.false.
        endsubroutine InitialCellMark
!----------------------------------------------------------------------
    end module ModMeshTools
!======================================================================
    subroutine GenerateBGMesh   ! BG -- back-ground
    use ModTypDef
    use ModMesh
    use ModInpMesh
    use ModMeshTools
    implicit none
    integer :: i, j, k
    type(octCell),pointer :: t

    BGCellSize=abs(DomainMin-DomainMax)/nCell
    nBGCells=nCell(1)*nCell(2)*nCell(3)
    nCells=0
    ALLOCATE(Cell(nCell(1),nCell(2),nCell(3)))
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        t=>Cell(i, j, k)
        nCells=nCells+1
            t%nBGCell     = [i, j, k]
            t%nCell       = nCells
            t%lvl         = 0
            t%fSplitType  = 0
            t%Location    = 0
            t%Node        = 0
            t%Mark        = .false.
            t%Center      = (/DomainMin(1)+(i-0.5)*BGCellSize(1),   &
                              DomainMin(2)+(j-0.5)*BGCellSize(2),   &
                              DomainMin(3)+(k-0.5)*BGCellSize(3)/)
            t%U           = 0
            t%cross       = -5
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
    use ModTypDef
    use ModMesh
    use ModInpMesh
    use ModMeshTools
    use ModNeighbor
    implicit none
    integer :: i, j, k
    type(octCell),pointer :: t
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    do k = 1, nCell(3)
    do j = 1, nCell(2)
    do i = 1, nCell(1)
        t=>Cell(i,j,k)
        call FindNeighbor(t)
    enddo
    enddo
    enddo
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "FindNeighbor time: ", tEnd-tStart
    contains
!----------------------------------------------------------------------
        recursive subroutine FindNeighbor(c)
        implicit none
        type(octCell),pointer :: c

        if(ASSOCIATED(c%son8))then
            call FindNeighbor(c%son1)
            call FindNeighbor(c%son2)
            call FindNeighbor(c%son3)
            call FindNeighbor(c%son4)
            call FindNeighbor(c%son5)
            call FindNeighbor(c%son6)
            call FindNeighbor(c%son7)
            call FindNeighbor(c%son8)
            return
        elseif(ASSOCIATED(c%son4))then
            call FindNeighbor(c%son1)
            call FindNeighbor(c%son2)
            call FindNeighbor(c%son3)
            call FindNeighbor(c%son4)
            return
        elseif(ASSOCIATED(c%son2))then
            call FindNeighbor(c%son1)
            call FindNeighbor(c%son2)
            return
        endif
        c%NeighborX1=>NeighborX1(c)
        c%NeighborX2=>NeighborX2(c)
        c%NeighborY1=>NeighborY1(c)
        c%NeighborY2=>NeighborY2(c)
        c%NeighborZ1=>NeighborZ1(c)
        c%NeighborZ2=>NeighborZ2(c)
        endsubroutine FindNeighbor
    endsubroutine initFindNeighbor
!======================================================================
    subroutine initSurfaceAdapt
    use ModMesh
    use ModMeshTools
    use ModInpMesh
    use ModNeighbor
    implicit none
    type(octCell),pointer :: t
    integer :: i, j, k
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time
    real(R8):: p        ! Used to print precentage
    integer :: step=0   ! Counter, used to print precentage

    call CPU_TIME(tStart)
    write(6,'(1X,A,12X,A)',advance='no') 'SurfaceAdapt progress:', ''
    flush(6)
    select case (cIntersectMethod)
    case (1)    ! Ray-cast with Painting Algorithm Method
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t       =>Cell(i,j,k)
            if (CellCast(t)) then
                t%cross = -3
                call SurfaceAdapt(t)
            else
                t%cross = -4
            endif
            step=step+1
            p=step/real(nBGCells,R8)*100
            write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
            flush(6)
        enddo
        enddo
        enddo
        write(*,*) '' ! Stop write with advance='no'
        call initFindNeighbor
        call initSmoothMesh
        call initFindNeighbor
        call initPaintingAlgorithm ! Painting Algorithm Method
        cIntersectMethod = 2 ! Close the Painting Algorithm Method

    case (2)    ! Ray-cast only
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t       =>Cell(i,j,k)
            if (CellCast(t)) then
                t%cross = -3
                t%cross = CellInout(t)
                call SurfaceAdapt(t)
            else
                t%cross = -4
                t%cross = CellInout(t)
            endif
            step=step+1
            p=step/real(nBGCells,R8)*100
            write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
            flush(6)
        enddo
        enddo
        enddo
        write(*,*) '' ! Stop write with advance='no'
        call initFindNeighbor
        call initSmoothMesh
        call initFindNeighbor

    case (3)    ! AABB with Painting Algorithm Method
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t       =>Cell(i,j,k)
            if (AABB(t)) then
                t%cross = -3
                call SurfaceAdapt(t)
            else
                t%cross = -4
            endif
            step=step+1
            p=step/real(nBGCells,R8)*100
            write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
            flush(6)
        enddo
        enddo
        enddo
        write(*,*) '' ! Stop write with advance='no'
        call initFindNeighbor
        call initSmoothMesh
        call initFindNeighbor
        call initPaintingAlgorithm ! Painting Algorithm Method
        cIntersectMethod = 4 ! Close the Painting Algorithm Method

    case (4)    ! AABB only
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t       =>Cell(i,j,k)
            if (AABB(t)) then
                t%cross = -3
                t%cross = CellInout(t)
                call SurfaceAdapt(t)
            else
                t%cross = -4
                t%cross = CellInout(t)
            endif
            step=step+1
            p=step/real(nBGCells,R8)*100
            write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
            flush(6)
        enddo
        enddo
        enddo
        write(*,*) '' ! Stop write with advance='no'
        call initFindNeighbor
        call initSmoothMesh
        call initFindNeighbor

    case (5)    ! AABBTraverse
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t       =>Cell(i,j,k)
            if (AABBTraverse(t)) then
                t%cross = -3
                t%cross = CellInout(t)
                call SurfaceAdapt(t)
            else
                t%cross = -4
                t%cross = CellInout(t)
            endif
            step=step+1
            p=step/real(nBGCells,R8)*100
            write(6,'(A,F5.1,A)',advance='no') '\b\b\b\b\b\b', p, '%'
            flush(6)
        enddo
        enddo
        enddo
        write(*,*) '' ! Stop write with advance='no'
        call initFindNeighbor
        call initSmoothMesh
        call initFindNeighbor
    end select


    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "SurfaceAdapt time: ", tEnd-tStart
    contains
!----------------------------------------------------------------------
        recursive subroutine SurfaceAdapt(c)
        implicit none
        type(octCell),pointer :: c
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
    subroutine initPaintingAlgorithm
    use ModMesh
    use ModMeshTools
    use ModInpMesh
    use ModNeighbor
    use ModPrecision
    implicit none
    type(octCell),pointer :: t
    integer :: ii, jj, kk
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    loop:  do kk = 1, nCell(3)
    do jj = 1, nCell(2)
    do ii = 1, nCell(1)
        t       =>Cell(ii,jj,kk)
        if (initPaintingAlgorithm2(t)) exit loop
    enddo
    enddo
    enddo loop
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "Painting Algorithm time: ", tEnd-tStart
    contains
!----------------------------------------------------------------------
        recursive function initPaintingAlgorithm2(c1) result(PA)
        ! PA: Have run/not run the subroutine PaintingAlgorithm
        implicit none
        logical               :: PA
        type(octCell),pointer :: c1
        logical,SAVE          :: PAused=.false.
        ! If used PaintingAlgorithm, PAused=.T.

        if (PAused) return
        if(ASSOCIATED(c1%son8))then
            PA = initPaintingAlgorithm2(c1%son1)
            PA = initPaintingAlgorithm2(c1%son2)
            PA = initPaintingAlgorithm2(c1%son3)
            PA = initPaintingAlgorithm2(c1%son4)
            PA = initPaintingAlgorithm2(c1%son5)
            PA = initPaintingAlgorithm2(c1%son6)
            PA = initPaintingAlgorithm2(c1%son7)
            PA = initPaintingAlgorithm2(c1%son8)
            return
        elseif(ASSOCIATED(c1%son4))then
            PA = initPaintingAlgorithm2(c1%son1)
            PA = initPaintingAlgorithm2(c1%son2)
            PA = initPaintingAlgorithm2(c1%son3)
            PA = initPaintingAlgorithm2(c1%son4)
            return
        elseif(ASSOCIATED(c1%son2))then
            PA = initPaintingAlgorithm2(c1%son1)
            PA = initPaintingAlgorithm2(c1%son2)
            return
        endif

        PA = .false.
        c1%cross = CellInout(c1)
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
        endfunction initPaintingAlgorithm2
!----------------------------------------------------------------------
        recursive subroutine PaintingAlgorithm(c,dirct)
        implicit none
        type(octCell),pointer :: c, cc
        integer               :: dirct

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

        if (c%cross/=-3 .and. c%cross/=-4) return ! Cell has been paintted.
        c%cross = CellInout(c)
        if (c%cross==3) return ! Inside cell, return.

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
    endsubroutine initPaintingAlgorithm
!======================================================================
    subroutine initSmoothMesh
    use ModPrecision
    use ModMesh
    use ModTypDef
    use ModMeshTools
    use ModInpMesh
    implicit none
    type(octCell),POINTER :: t
    integer               :: i, j, k, ii
    real(R8):: tStart   ! Start time
    real(R8):: tEnd     ! End time

    call CPU_TIME(tStart)
    do ii=1, InitRefineLVL+1

        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t=>Cell(i, j, k)
            call InitialCellMark(t)
        enddo
        enddo
        enddo
        
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t=>Cell(i, j, k)
            call PreSmoothMesh(t)
        enddo
        enddo
        enddo

        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t=>Cell(i, j, k)
            call SmoothMesh(t)
        enddo
        enddo
        enddo

        call initFindNeighbor

    enddo
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "Smooth Mesh time: ", tEnd-tStart
    contains
!----------------------------------------------------------------------
        recursive subroutine PreSmoothMesh(c)
        implicit none
        type(octCell),POINTER :: c,cx1,cx2,cy1,cy2,cz1,cz2
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
        ! Hole Cell
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
        recursive subroutine SmoothMesh(c)
        implicit none
        type(octCell),POINTER :: c
        ! mark(6)  1 x refine; 2 y refine; 3 z refine; 
        !          4 x coarse; 5 y coarse; 6 z coarse; 
        if (maxval(c%lvl)>=InitRefineLVL) return
        if(ASSOCIATED(c%son8))then
            call SmoothMesh(c%son1)
            call SmoothMesh(c%son2)
            call SmoothMesh(c%son3)
            call SmoothMesh(c%son4)
            call SmoothMesh(c%son5)
            call SmoothMesh(c%son6)
            call SmoothMesh(c%son7)
            call SmoothMesh(c%son8)
            return
        elseif(ASSOCIATED(c%son4))then
            call SmoothMesh(c%son1)
            call SmoothMesh(c%son2)
            call SmoothMesh(c%son3)
            call SmoothMesh(c%son4)
            return
        elseif(ASSOCIATED(c%son2))then
            call SmoothMesh(c%son1)
            call SmoothMesh(c%son2)
            return
        endif
        

        ! Refine
        if     ( c%mark(1) .and. c%mark(2) .and. c%mark(3)) then
            call NewCell(c,0)
        elseif ( c%mark(1) .and. .not.c%mark(2) .and. .not.c%mark(3) ) then
            call NewCell(c,1)
        elseif ( .not.c%mark(1) .and. c%mark(2) .and. .not.c%mark(3) ) then
            call NewCell(c,2)
        elseif ( .not.c%mark(1) .and. .not.c%mark(2) .and. c%mark(3) ) then
            call NewCell(c,3)
        elseif ( c%mark(1) .and. c%mark(2) .and. .not.c%mark(3) ) then
            call NewCell(c,4)
        elseif ( c%mark(1) .and. .not.c%mark(2) .and. c%mark(3) ) then
            call NewCell(c,5)
        elseif ( .not.c%mark(1) .and. c%mark(2) .and. c%mark(3) ) then
            call NewCell(c,6)
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
        endif
        endsubroutine SmoothMesh
!----------------------------------------------------------------------
    endsubroutine initSmoothMesh
!======================================================================
!======================================================================
!----------------------------------------------------------------------
