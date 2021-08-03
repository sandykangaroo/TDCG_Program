!======================================================================
 module ModWalldistance
    use ModPrecision
    use geometry_mod2
    implicit none
    integer :: step
    real(R8):: TimeStep
    
    contains   
!---------------------------------------------------------------------- 
    recursive subroutine GetWallDistance
        use ModPrecision
        use ModMesh
        use ModTypDef
        use ModMeshTools
        use ModInpMesh
        implicit none
        type(FttCell),POINTER :: t
        integer               :: i, j, k
    
        do k = 1, nCell(3)
        do j = 1, nCell(2)
        do i = 1, nCell(1)
            t=>OctCell(i, j, k)
            call GetMinDistanceKDT(t)
        enddo
        enddo
        enddo
    endsubroutine GetWallDistance  
!----------------------------------------------------------------------     
    recursive subroutine GetMinDistanceKDT(c)
        use ModPrecision
        use ModTypDef
        use ModMesh
        use ModKDTree 
        use ModInpMesh
        use ModMeshTools
        use ModGeometry
        use geometry_mod2
        implicit none
    
        integer ::  l, h, is
        type(FttCell),pointer           :: c,cs
        type(typKDTtree),pointer        :: tp => null()
        real(R8)                        :: mindis
        type(KDT_node), pointer         :: nearest
        type(typpoint)                  :: tem
        type (triangle)                 :: tri
        
        
            if(ASSOCIATED(c%Octson))then
                do is=1,c%Octson%nSon
                    cs=>c%Octson%son(is)
                    call GetMinDistanceKDT(cs)
                enddo
                return
            endif
            
            tp=>kdtree(1)
            tem%P(:)=c%center(:)
            nearest => tp%root
            h=0
            call nearest_search(tem, tp%root, nearest, h)
            tri%P(:) = nearest%the_data%p(:)
            mindis=DisBetweenPointAndTri(tem,tri) 
            !write(102,*)mindis
            c%WallDistance=mindis
            
  end subroutine GetMinDistanceKDT 
    
    
    
    
!     subroutine GetMinDistance
!     use ModPrecision
!     use ModTypDef
!     use ModMesh
!     use ModKDTree
!     use ModInpMesh
!     use ModMeshTools
!     use ModGeometry
!     use geometry_mod2
!     implicit none

!     integer :: i, j, k
!     type(FTTCell),pointer           :: t
!     real(R8)                        :: tstart,tend

!     call CPU_TIME(tstart)
!     do i=1,nCell(1)
!     do j=1,nCell(2)
!     do k=1,nCell(3)
!         t=>OctCell(i,j,k)
!         call GetMinDistanceKDT(t)
!     enddo
!     enddo
!     enddo
!     call CPU_TIME (tend)
!     write(*,'(1X,A,F10.2)') "KDT find min distance time", tend-tstart
!     contains
! !----------------------------------------------------------------------
!         recursive subroutine GetMinDistanceKDT(c)
!         implicit none

!         integer ::  l, h
!         type(FTTCell),pointer           :: c
!         type(typKDTtree),pointer        :: tp => null()
!         real(R8)                        :: tempoint(3),mindis
!         type(KDT_node), pointer         :: nearest
!         type(typpoint)                  :: tem

!         if(ASSOCIATED(c%son8))then
!             call GetMinDistanceKDT(c%son1)
!             call GetMinDistanceKDT(c%son2)
!             call GetMinDistanceKDT(c%son3)
!             call GetMinDistanceKDT(c%son4)
!             call GetMinDistanceKDT(c%son5)
!             call GetMinDistanceKDT(c%son6)
!             call GetMinDistanceKDT(c%son7)
!             call GetMinDistanceKDT(c%son8)
!             return
!         elseif(ASSOCIATED(c%son4))then
!             call GetMinDistanceKDT(c%son1)
!             call GetMinDistanceKDT(c%son2)
!             call GetMinDistanceKDT(c%son3)
!             call GetMinDistanceKDT(c%son4)
!             return
!         elseif(ASSOCIATED(c%son2))then
!             call GetMinDistanceKDT(c%son1)
!             call GetMinDistanceKDT(c%son2)
!             return
!         endif

!         tp=>kdtree(1)
!         tem%P(:)=c%center(:)
!         nearest => tp%root
!         h=0
!         call nearest_search(tem, tp%root, nearest, h)
!         mindis=distance(tem,nearest%the_data%p(4))
!         write(102,*)h
!         end subroutine GetMinDistanceKDT
!     end subroutine GetMinDistance
!======================================================================
endmodule ModWalldistance
