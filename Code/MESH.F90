!======================================================================
    module ModMesh
        use ModTypDef
        implicit none
        ! Mesh
        integer :: nBGCells     ! Number of the background Cells.
        integer :: nCells     ! Number of total Cells.
        real(R8):: BGStep(3)    ! Step size for background cell.
        type(typCell),pointer :: Cell(:,:,:)
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
        use ModKDTree
        use ModInpGlobal, only: nGeometry
        implicit none
        type(typCell),pointer :: c
        integer ::i

        select case (cIntersectMethod)
        case(1)
            initCellIntersect=CellCast(c)
        case(2)
            if (AABB(c)) then
                do i = 1,nGeometry
                    if (RayCast(c%Center,KDTree(i)%root)==1) initCellIntersect=2
                    if (RayCast(c%Center,KDTree(i)%root)==0) initCellIntersect=1
                enddo
            else
                do i = 1,nGeometry
                    if (RayCast(c%Center,KDTree(i)%root)==1) initCellIntersect=3
                    if (RayCast(c%Center,KDTree(i)%root)==0) initCellIntersect=0
                enddo
            endif
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
            ! Quick Bounding-BOX identify
            if (BBOX(p(i)%P,KDTree(ii)%root%box)) then ! inside
                Pintersect=Pintersect+RayCast(p(i)%P,KDTree(ii)%root)
            endif
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
            RayCast=1; return
        else
            RayCast=0; return
        endif
        endfunction RayCast
!----------------------------------------------------------------------
        Logical function AABB(c)!1=intersect,0=no intersect
            use ModKDTree
            use ModInpGlobal
            use ModGeometry
            implicit none
            type(typCell),pointer :: c
            type(triangle)        :: tri
            integer               :: i,ii,ng
            real(R8)              :: boxCell(6)
            type(KDT_node),pointer:: node
            logical               :: aaa
            type(typKDTtree), pointer   :: tp => null()
            
            aaa=.false.
            boxCell(1)=c%center(1)-BGStep(1)/2**(c%lvl(1)+1)
            boxCell(2)=c%center(2)-BGStep(2)/2**(c%lvl(2)+1)
            boxCell(3)=c%center(3)-BGStep(3)/2**(c%lvl(3)+1)
            boxCell(4)=c%center(1)+BGStep(1)/2**(c%lvl(1)+1)
            boxCell(5)=c%center(2)+BGStep(2)/2**(c%lvl(2)+1)
            boxCell(6)=c%center(3)+BGStep(3)/2**(c%lvl(3)+1)
            
            tp => kdtree(1)
            node=>tp%root
            if(boxCell(1)>node%box(4).or.boxCell(2)>node%box(5).or.boxCell(3)>node%box(6).or.&
                &boxCell(4)<node%box(1).or.boxCell(5)<node%box(2).or.boxCell(6)<node%box(3))then
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

!Traverse            
            !Loop1:do ng=1,nGeometry
            !    do i=1, body(ng)%nse
            !        tri= body(ng)%se3d(i)
            !        AABB=TriBoxOverlap (c,tri)
            !        if(AABB) exit Loop1
            !    enddo
            !enddo Loop1
            !return         
        end function AABB
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
            type(typCell),pointer                   :: c
            
            tri = node%the_data
            aaa = TriBoxOverlap(c,tri)
         
            if(aaa)then
                return
            endif
            
            split=node%splitaxis 
            splitmax=node%splitaxis+3
            if(boxCell(split)>node%the_data%p(4)%p(split))then
                if(associated(node%right))then
                    call KdFindTri(c,boxCell,node%right,aaa)
                    if(aaa)then
                      return
                    endif
                endif
            elseif(boxCell(splitmax)<node%the_data%p(4)%p(split))then
                if(associated(node%left))then
                    call KdFindTri(c,boxCell,node%left,aaa)
                    if(aaa)then
                        return
                    endif
                endif
            else
                if(associated(node%right))then
                    call KdFindTri(c,boxCell,node%right,aaa)
                    if(aaa)then
                        return
                    elseif(associated(node%left))then
                        call KdFindTri(c,boxCell,node%left,aaa)
                        if(aaa)then
                            return
                        endif
                    endif
                endif
            endif
        end subroutine KdFindTri    
!----------------------------------------------------------------------
!/* use separating axis theorem to test overlap between triangle and box */
!/* need to test for overlap in these directions: */
!/* 1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
!/* we do not even need to test these) */
!/* 2) normal of the triangle */
!/* 3) crossproduct(edge from tri, {x,y,z}-directin) */
!/* this gives 3x3=9 more tests */
        Logical function TriBoxOverlap (c,tri)
            use ModKDTree
            use ModInpGlobal
            use ModGeometry
            use ModTools
            implicit none
            
            type(typCell),pointer :: c
            type(triangle)        :: tri
            real(R8)              :: v0(3), v1(3), v2(3), boxcenter(3), normal(3), e0(3), e1(3), e2(3),boxhalfsize(3)
            real(R8)              :: min, max, rad, d, dd, fex, fey, fez,a
            real(R8)              :: X01P0, X01P2,Y02P0,Y02P2,Z12P2,Z12P1,Z0P0,Z0P1,X2P0,X2P1,Y1P0,Y1P1
        
            boxcenter(1)=c%center(1)
            boxcenter(2)=c%center(2)
            boxcenter(3)=c%center(3)
            boxhalfsize(1)= BGStep(1)/2**(c%lvl(1)+1)
            boxhalfsize(2)= BGStep(2)/2**(c%lvl(2)+1)
            boxhalfsize(3)= BGStep(3)/2**(c%lvl(3)+1)
            
            
            v0(:)=Sub(tri%p(1)%P(:),boxcenter(:))
            v1(:)=Sub(tri%p(2)%P(:),boxcenter(:))
            v2(:)=Sub(tri%p(3)%P(:),boxcenter(:))
            
            !compute triangle edges
            e0(:)=Sub(v1(:),v0(:))
            e1(:)=Sub(v2(:),v1(:))
            e2(:)=Sub(v0(:),v2(:))
            
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
    end module ModMeshTools
!======================================================================
    subroutine GenerateBGMesh   ! BG -- back-ground
    use ModPrecision
    use ModTypDef
    use ModMesh
    use ModInpMesh
    use ModMeshTools
    implicit none
    integer :: i, j, k
    type(typCell),pointer :: t

    BGStep=abs(Domain1-Domain2)/nCell
    nBGCells=nCell(1)*nCell(2)*nCell(3)
    nCells=0
    ALLOCATE(Cell(nCell(1),nCell(2),nCell(3)))
    do i = 1, nCell(1)
    do j = 1, nCell(2)
    do k = 1, nCell(3)
        t=>Cell(i, j, k)
        nCells=nCells+1
            t%nBGCell     = [i, j, k]
            t%nCell       = nCells
            t%lvl         = 0
            t%fSplitType  = 0
            t%Location    = 0
            t%Node        = 0
            t%Center      = (/Domain1(1)+(i-0.5)*BGStep(1),         &
                              Domain1(2)+(j-0.5)*BGStep(2),         &
                              Domain1(3)+(k-0.5)*BGStep(3)/)
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
    call BGMeshCross
    contains
!----------------------------------------------------------------------
        subroutine BGMeshCross
        implicit none

        if (PaintingAlgorithmMethod) then
            loop1: do i = 1, nCell(1)
                    do j = 1, nCell(2)
                    do k = 1, nCell(3)
                        t=>Cell(i,j,k)
                        t%cross=initCellIntersect(t)
                        if (t%cross==0) then
                            call PaintingAlgorithm(t,0)
                            exit loop1
                        endif
                    enddo
                    enddo
            enddo loop1
        else
            do i = 1, nCell(1)
            do j = 1, nCell(2)
            do k = 1, nCell(3)
                t=>Cell(i,j,k)
                t%cross=initCellIntersect(t)
            enddo
            enddo
            enddo
        endif
        endsubroutine BGMeshCross
!----------------------------------------------------------------------
        recursive subroutine PaintingAlgorithm(c,dirct)
        implicit none
        type(typCell),pointer :: c, cc
        integer               :: dirct

        if (dirct /= 0) then
            if (c%cross /=-5) return    ! Cell has been paintted, return.
            c%cross = initCellIntersect(c)
            if (c%cross == 3) return    ! Inside cell, return.
        endif

        if (dirct /= 4) then
            if (c%nBGCell(1)<nCell(1)) then
                cc=>Cell(c%nBGCell(1)+1,c%nBGCell(2),c%nBGCell(3))
                call PaintingAlgorithm(cc,1)
            endif
        endif
        if (dirct /= 1) then
            if (c%nBGCell(1)>1) then
                cc=>Cell(c%nBGCell(1)-1,c%nBGCell(2),c%nBGCell(3))
                call PaintingAlgorithm(cc,4)
            endif
        endif
        if (dirct /= 5) then
            if (c%nBGCell(2)<nCell(2)) then
                cc=>Cell(c%nBGCell(1),c%nBGCell(2)+1,c%nBGCell(3))
                call PaintingAlgorithm(cc,2)
            endif
        endif
        if (dirct /= 2) then
            if (c%nBGCell(2)>1) then
                cc=>Cell(c%nBGCell(1),c%nBGCell(2)-1,c%nBGCell(3))
                call PaintingAlgorithm(cc,5)
            endif
        endif
        if (dirct /= 6) then
            if (c%nBGCell(3)<nCell(3)) then
                cc=>Cell(c%nBGCell(1),c%nBGCell(2),c%nBGCell(3)+1)
                call PaintingAlgorithm(cc,3)
            endif
        endif
        if (dirct /= 3) then
            if (c%nBGCell(3)>1) then
                cc=>Cell(c%nBGCell(1),c%nBGCell(2),c%nBGCell(3)-1)
                call PaintingAlgorithm(cc,6)
            endif
        endif
        endsubroutine PaintingAlgorithm
    endsubroutine GenerateBGMesh
!======================================================================
    subroutine initSurfaceAdapt
    use ModMesh
    use ModMeshTools
    use ModInpMesh
    implicit none
    type(typCell),pointer :: t
    integer :: i, j, k

    do i = 1, nCell(1)
    do j = 1, nCell(2)
    do k = 1, nCell(3)
        t=>Cell(i, j, k)
        call SurfaceAdapt(t)
    enddo
    enddo
    enddo
        contains
!----------------------------------------------------------------------
        recursive subroutine SurfaceAdapt(c)
    implicit none
        type(typCell),pointer :: c
        integer :: ii

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
