!======================================================================
    module ModSolve
    use ModPrecision
    implicit none
    integer :: step
    real(R8):: TimeStep
    endmodule ModSolve
!======================================================================
    subroutine GetMinDistance
    use ModPrecision
    use ModTypDef
    use ModMesh
    use ModKDTree 
    use ModInpMesh
    use ModMeshTools
    use ModGeometry
    use geometry_mod2
    implicit none
    
    integer :: i, j, k, l, h, m
    type(octCell),pointer           :: t
    type(typKDTtree),pointer        :: tp => null()
    real(R8)                        :: tempoint(3),mindis
    real(R8)                        :: tstart,tend
    type(KDT_node), pointer         :: nearest
    
    
    
    !call CPU_TIME(tstart)
    !do i=1,nCell(1)
    !    do j=1,nCell(2)
    !        do k=1,nCell(3)
    !            t=>Cell(i,j,k)
    !            tempoint(:)=t%center(:)
    !            tem%P(:)= tempoint(:)
    !            do l=1,body(1)%nse
    !                mindis=min(mindis,distance(tem,body(1)%se3d(l)%p(4)))
    !                write(5,*) mindis
    !                close(5)
    !            enddo
    !        enddo
    !    enddo
    !enddo
    !call CPU_TIME (tend)
    !write(*,*) "traverse find min distance time", tend-tstart
    !
    !
    !call CPU_TIME(tstart)
    !do i=1,nCell(1)
    !    do j=1,nCell(2)
    !        do k=1,nCell(3)
    !            t=>Cell(i,j,k)
    !            tempoint(:)=t%center(:)
    !            tem%P(:)= tempoint(:)
    !            tp=>kdtree(1)
    !            nearest => tp%root
    !            h=0
    !            call nearest_search(tem, tp%root, nearest, h)
    !            mindis=distance(tem,nearest%the_data%p(4))
    !
    !                write(6,*) mindis
    !                close(6)
    !        enddo
    !    enddo
    !enddo
    !call CPU_TIME (tend)
    !write(*,*) "KDT find min distance time", tend-tstart
        
    
    
    
    end subroutine GetMinDistance
!======================================================================
!======================================================================
    