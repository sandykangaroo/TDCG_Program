!======================================================================
    module geometry_mod2
    use ModPrecision
    use ModTypDef
    use ModGeometry
    use ModKDTree
    implicit none

    contains
!----------------------------------------------------------------------
    subroutine create_KDT_tree_for_body_i(Tree, i)
    ! .. create the actual tree structure, given an input array of data.
    ! .. Return Value ..
        type(typKDTtree)       :: Tree
    ! .. Input Arguments ..
        Integer, Intent (In)   :: i      ! the ith body
    ! .. Local Arguments ..
        integer                :: n      ! number of triangles of body i
        integer                :: j      ! loop parameter
        integer                :: depth  ! level of node in the tree
        type(triangle), pointer:: res(:) ! store triangles of body i
        type(KDT_node), pointer:: p
        real(R8)               :: box(6)

        n =  body(i)%nse
        depth = 1
        res => body(i)%se3d
        p => null()
        box(:) = body(i)%box(:)
        Tree%root => build_tree(res, p, box, depth)
        nullify(res)

    end subroutine create_KDT_tree_for_body_i
!----------------------------------------------------------------------

    recursive function build_tree(res, p, box, depth) result (td)
        use ModInpMesh,only : useKDT
    ! .. Return Value ..
        type(KDT_node), pointer              :: td
    ! .. Input Arguments ..
        type(KDT_node), pointer              :: p, pt
        type(triangle), pointer              :: res(:), resB(:)
        integer                              :: depth, depth2
        real(R8)                             :: box(6), b(6), bt
    ! .. Local Arguments ..
        integer                              :: n, j, mid

        n = size(res)
        allocate(td)
        allocate(td%the_data)
        mid = (1+n)/2
        td%level = depth
        td%parent => p
        td%box(:) = box(:)

        b(:) = box(:)
        depth2 = depth + 1
        pt => td
        if ( n == 1) then
            td%the_data = res(mid)
            td%splitaxis = 1
            td%left => null()
            td%right => null()
        else if ( n == 2 ) then
            td%the_data = res(mid)
            td%splitaxis = 1
            resB => res(2:2)
            if ( res(2)%p(4)%P(1) < res(1)%p(4)%P(1) ) then
                b(4) = max(resB(1)%p(1)%P(1), resB(1)%p(2)%P(1), resB(1)%p(3)%P(1))
                td%left => build_tree(resB, pt, b, depth2)
                td%right => null()
            else ! if ( res(2)%p(4)%P(1) >= res(2)%p(4)%P(1) ) then
                b(1) = min(resB(1)%p(1)%P(1), resB(1)%p(2)%P(1), resB(1)%p(3)%P(1))
                td%left => null()
                td%right => build_tree(resB, pt, b, depth2)
            end if
        else
            if (useKDT) then
                !find the most_spread_direction and define it as split direction
                ! td%SplitAxis = mod(depth,3)+1
                call find_split_direction(res, td%splitaxis)
                call sort_under_split_direction(res, td%splitaxis)
            else
            ! define the split direction alternative
                if ( mod(td%level,3) == 1 ) then
                td%splitaxis = 1
                else if ( mod(td%level,3) == 2 ) then
                td%splitaxis = 2
                else if ( mod(td%level,3) == 0 ) then
                    td%splitaxis = 3
                end if
            endif
            td%the_data = res(mid)
            !       build left child
            resB => res(1:mid - 1)
            if ( td%splitaxis == 1 ) then
                bt = b(4)
                b(4) = max(resB(1)%p(1)%P(1), resB(1)%p(2)%P(1), resB(1)%p(3)%P(1))
                do j = 2, mid - 1
                    b(4) = max(b(4), resB(j)%p(1)%P(1), resB(j)%p(2)%P(1), resB(j)%p(3)%P(1))
                end do
            elseif ( td%splitaxis == 2 ) then
                bt = b(5)
                b(5) = max(resB(1)%p(1)%P(2), resB(1)%p(2)%P(2), resB(1)%p(3)%P(2))
                do j = 2, mid - 1
                    b(5) = max(b(5), resB(j)%p(1)%P(2), resB(j)%p(2)%P(2), resB(j)%p(3)%P(2))
                end do
            elseif ( td%splitaxis == 3 ) then
                bt = b(6)
                b(6) = max(resB(1)%p(1)%P(3), resB(1)%p(2)%P(3), resB(1)%p(3)%P(3))
                do j = 2, mid - 1
                    b(6) = max(b(6), resB(j)%p(1)%P(3), resB(j)%p(2)%P(3), resB(j)%p(3)%P(3))
                end do
            end if
            td%left => build_tree(resB, pt, b, depth2)
!       build right child
           resB => res(mid + 1 : n)
            if ( td%splitaxis == 1 ) then
                b(4) = bt
                b(1) = min(resB(1)%p(1)%P(1), resB(1)%p(2)%P(1), resB(1)%p(3)%P(1))
                do j = 2, n - mid
                    b(1) = min(b(1), resB(j)%p(1)%P(1), resB(j)%p(2)%P(1), resB(j)%p(3)%P(1))
                end do
            elseif ( td%splitaxis == 2 ) then
                b(5) = bt
                b(2) = min(resB(1)%p(1)%P(2), resB(1)%p(2)%P(2), resB(1)%p(3)%P(2))
                do j = 2, n - mid
                    b(2) = min(b(2), resB(j)%p(1)%P(2), resB(j)%p(2)%P(2), resB(j)%p(3)%P(2))
                end do
            elseif ( td%splitaxis == 3 ) then
                b(6) = bt
                b(3) = min(resB(1)%p(1)%P(3), resB(1)%p(2)%P(3), resB(1)%p(3)%P(3))
                do j = 2, n - mid
                    b(3) = min(b(3), resB(j)%p(1)%P(3), resB(j)%p(2)%P(3), resB(j)%p(3)%P(3))
                end do
            end if
            td%right => build_tree(resB, pt, b, depth2)
        end if

!        write(*,*) depth
        return
    end function build_tree
!----------------------------------------------------------------------

    subroutine  find_split_direction(res, split)
    ! find the most_spread_direction and define it as split direction
    ! .. Input Arguments ..
        type(triangle), pointer :: res(:)  !
    ! .. Output Arguments ..
        integer                 :: split ! split direction, 1 x; 2 y; 3 z.
    ! .. Local Arguments ..
        integer                 :: n, j
        real(R8)                :: summ(3), average(3), variance(3)
        real(R8)                :: var

        n = size(res)
        summ = 0.d0
        average = 0.d0
        variance = 0.d0
        do j = 1, n
            summ(:) = summ(:) + res(j)%p(4)%P(:)
        end do
        average = summ / n
        do j = 1, n
            variance(:) = variance(:) + (res(j)%p(4)%P(:) - average(:))*(res(j)%p(4)%P(:)- average(:))
        end do
        var = -1.d0
        do j = 1, 3
            if (variance(j) > var) then
                var = variance(j)
                split = j
            end if
        end do
!        write(*,*) variance, split

    end subroutine  find_split_direction
!----------------------------------------------------------------------

    subroutine sort_under_split_direction(res, split)
    ! sort the data (res) under the direction (split), and store the
    ! results into (sort)
    ! .. Input Arguments ..
        integer                 :: split    ! split direction
    ! .. Input/Output Arguments ..
        type(triangle), pointer :: res(:)
    ! .. Local Arguments ..
        real(R8), allocatable   :: A(:)
        integer                 :: i, n

        n = size(res)
        allocate(A(n))
        do i =1, n
            A(i) = res(i)%p(4)%P(split)
        end do

        call quicksort(res, A)

    end subroutine sort_under_split_direction
!----------------------------------------------------------------------

    recursive subroutine quicksort(res, A)
        type(triangle), pointer                   :: res(:), B(:)
        integer                                   :: iq, n, i
        real(R8), intent(in out)                :: A(:)

        B => res

        if(size(A) > 1) then
            call Partition(B, A, iq)
            !write(*,*) size(res),size(a),size(b),iq
            B => res(:iq-1)
             !write(*,*) 'hi'
            call quicksort(B, A(:iq-1))
             !write(*,*) 'hello'
            B => res(iq:)
            call quicksort(B, A(iq:))
        endif

    end subroutine quicksort
!----------------------------------------------------------------------

    subroutine Partition(res, A, marker)
        type(triangle), pointer :: res(:)
        real(R8), intent(inout) :: A(:)
        integer, intent(out)    :: marker
        integer                 :: i, j, k, l
        type(triangle), pointer :: temp, Ai, Aj
        real(R8)                :: tempx, x ! pivot point

        allocate(Ai, Aj)

        x = A(1)
        i= 0
        j= size(A) + 1
        do
            j = j-1
            do
                if (A(j) <= x) exit
                j = j-1
            end do
            i = i+1
            do
                if (A(i) >= x) exit
                i = i+1
            end do
            if (i < j) then
            ! exchange A(i) and A(j)
                do l =1, 4
                    do k =1,3
                        Ai%p(l)%P(k) = res(i)%p(l)%P(k)
                    end do
                    Ai%p(l)%label = res(i)%p(l)%label
                end do
                do l =1, 4
                    do k =1,3
                        Aj%p(l)%P(k) = res(j)%p(l)%P(k)
                    end do
                    Aj%p(l)%label = res(j)%p(l)%label
                end do
 !               write(*,*) A(i)%p(4)%label, A(j)%p(4)%label, Ai%p(4)%label, Aj%p(4)%label
                temp => Ai
  !              write(*,*) A(i)%p(4)%label, A(j)%p(4)%label, Ai%p(4)%label, Aj%p(4)%label,  temp%p(4)%label
                Ai => Aj
                Aj => temp
   !             write(*,*) A(i)%p(4)%label, A(j)%p(4)%label, Ai%p(4)%label, Aj%p(4)%label,  temp%p(4)%label
                do l =1, 4
                    do k =1,3
                        res(i)%p(l)%P(k) = Ai%p(l)%P(k)
                    end do
                    res(i)%p(l)%label = Ai%p(l)%label
                end do
                do l =1, 4
                    do k =1,3
                        res(j)%p(l)%P(k) = Aj%p(l)%P(k)
                    end do
                    res(j)%p(l)%label = Aj%p(l)%label
                end do
    !            write(*,*) A(i)%p(4)%label, A(j)%p(4)%label, Ai%p(4)%label, Aj%p(4)%label,  temp%p(4)%label
                tempx = A(i)
                A(i) = A(j)
                A(j) = tempx
            elseif (i == j) then
                marker = i+1
                return
            else
                marker = i
                return
            endif
        end do
        deallocate(Ai, Aj)

    end subroutine Partition
!----------------------------------------------------------------------
    recursive subroutine KDTree_out(node)
        type(KDT_node), pointer                          :: node, p
        real(R8)                                       :: xmid, x(6), aa
        integer                                          :: i, mid
        logical                                          :: opd

        inquire(file='./Data/split2.dat', opened=opd)
        if ( .not. opd) then
            open(1, file='./Data/split2.dat', status='unknown', form='formatted')
        end if

        if( .not. associated(node%left) .and. &
            .not. associated(node%right) ) return

        p => node
        mid = p%splitaxis
        xmid = p%the_data%p(4)%P(mid)
!        write(1,*) ' ZONE T = "1" '
        write(1,*) ' VARIABLES = "X", "Y", "Z" '
        write(1,*) ' ZONE NODES=4, ELEMENTS=1, DATAPACKING=POINT, ZONETYPE=FEBRICK '
        if ( mid == 1) then
            write(1,*) xmid, p%box(2), p%box(3)
            write(1,*) xmid, p%box(2), p%box(6)
            write(1,*) xmid, p%box(5), p%box(6)
            write(1,*) xmid, p%box(5), p%box(3)
            write(1,*) '1 2 3 4 1 2 3 4'
        else if ( mid == 2) then
            write(1,*) p%box(1), xmid, p%box(3)
            write(1,*) p%box(1), xmid, p%box(6)
            write(1,*) p%box(4), xmid, p%box(6)
            write(1,*) p%box(4), xmid, p%box(3)
            write(1,*) '1 2 3 4 1 2 3 4'
        else if ( mid == 3) then
            write(1,*) p%box(1), p%box(2), xmid
            write(1,*) p%box(1), p%box(5), xmid
            write(1,*) p%box(4), p%box(5), xmid
            write(1,*) p%box(4), p%box(2), xmid
            write(1,*) '1 2 3 4 1 2 3 4'
        end if

        if( associated(node%left) .OR. associated(node%right) ) then
            if ( associated(node%left) .and. associated(node%right) ) then
                p => node%left
                call KDTree_out(p)
                p => node%right
                call KDTree_out(p)
            end if
            if ( associated(node%left) .and. .not. associated(node%right) ) then
                p => node%left
                call KDTree_out(p)
            end if
            if ( .not. associated(node%left) .and. associated(node%right) ) then
                p => node%right
                call KDTree_out(p)
            end if
        end if

        if (.not. associated(node%parent)) then
            close(1)
        end if
    end subroutine KDTree_out

!----------------------------------------------------------------------
!     recursive subroutine nearest_search(Tar, node, nearest, i)
!     !Find the nearest triangle of space point Tar in the KDTree begin with node
!!     .. Input Arguments ..
!        type(typpoint)                                    :: Tar
!        type(KDT_node), pointer              :: node
!!     .. Output Arguments ..
!        type(KDT_node), pointer              :: nearest
!!     .. Local Arguments ..
!        type(KDT_node), pointer              :: near, far
!        real(R8)                                       :: dist, dist1, dist2, aa
!        integer                                          :: i, split
!
!        i = i + 1
!        !write(102,*) node%level
!        aa = 1.d0
!        if ( .not. associated(node%left) .and. .not.associated(node%right) ) then
!            dist = distance(Tar, nearest%the_data%p(4))
!            dist2 = distance(Tar, node%the_data%p(4))
!            if ( dist2 - dist < 0.d0 ) then
!                nearest => node
!            end if
!            return
!        end if
!
!        split = node%splitaxis
!        dist1 = Tar%p(split) - node%the_data%p(4)%p(split)
!        aa = 1.d0
!        if ( dist1 <= 0.d0 ) then
!            near => node%left
!            far => node%right
!        else
!            near => node%right
!            far => node%left
!        end if
!
!        if ( ASSOCIATED(near) ) then
!            call nearest_search(Tar, near, nearest, i)
!        end if
!
!
!        dist = distance(Tar, nearest%the_data%p(4))
!        if ( dist <= abs(dist1) ) then
!            return
!        else if( dist > dist1 .and. dist < 1.d0) then
!        else if ( dist > abs(dist1) ) then
!            dist2 = distance(Tar, node%the_data%p(4))
!            if ( dist2 - dist < 0.d0 ) then
!                nearest => node
!                dist = dist2
!            end if
!            if ( ASSOCIATED(far) ) then
!                call nearest_search(Tar, far, nearest, i)
!            end if
!        end if
!        aa = 1.d0
!
!    end subroutine nearest_search

!----------------------------------------------------------------------
    recursive subroutine nearest_search(Tar, node, nearest, i)
    ! Find the nearest triangle of space point Tar in the KDTree begin with node
    ! .. Input Arguments ..
    real(R8),INTENT(IN)     :: Tar(3)
    type(KDT_node), pointer :: node
    ! .. Output Arguments ..
    type(KDT_node), pointer :: nearest
    ! .. Local Arguments ..
    type(KDT_node), pointer :: near, far
    real(R8)                :: dist, dist1, dist2
    integer                 :: i, split
    real(R8)                :: cy,cz,y1,z1,y2,z2,y0,z0,r
    logical                 :: SSC

    i = i + 1
    if ( .not. associated(node%left) .and. &
         .not. associated(node%right) ) then
        dist = distance(Tar, nearest%the_data%p(4)%P(1:3))
        dist2 = distance(Tar, node%the_data%p(4)%P(1:3))
        if ( dist2 < dist ) then
            nearest => node
        end if
        return
    end if

    split = node%splitaxis
    dist1 = Tar(split) - node%the_data%p(4)%p(split)
    if ( dist1 <= 0.d0 ) then
        near => node%left
        far => node%right
    else
        near => node%right
        far => node%left
    end if

    if ( ASSOCIATED(near) ) call nearest_search(Tar, near, nearest, i)

    dist = distance(Tar, nearest%the_data%p(4)%P(1:3))
    if ( dist <= abs(dist1) ) then
        return
    else if (split == 1) then
        cy = Tar(2)
        cz = Tar(3)
        r  = sqrt(dist*dist - dist1*dist1)
        y1 = node%box(2)
        z1 = node%box(3)
        y2 = node%box(5)
        z2 = node%box(6)
        y0 = (node%box(2) + node%box(5))/2
        z0 = (node%box(3) + node%box(6))/2
        SSC= SphereSplitsurfaceCross(cy,cz,y1,z1,y2,z2,y0,z0,r)
        if(dist > abs(dist1) .and. SSC)then
            dist2 = distance(Tar, node%the_data%p(4)%P(1:3))
            if ( dist2 - dist < 0.d0 ) then
                nearest => node
                dist = dist2
            end if
            if ( ASSOCIATED(far) ) then
                call nearest_search(Tar, far, nearest, i)
            end if
        endif
    elseif(split == 2)then
        cy = Tar(1)
        cz = Tar(3)
        r  = sqrt(dist*dist - dist1*dist1)
        y1 = node%box(1)
        z1 = node%box(3)
        y2 = node%box(4)
        z2 = node%box(6)
        y0 = (node%box(1) + node%box(4))/2
        z0 = (node%box(3) + node%box(6))/2
        SSC= SphereSplitsurfaceCross(cy,cz,y1,z1,y2,z2,y0,z0,r)
        if(dist > abs(dist1) .and. SSC)then
            dist2 = distance(Tar, node%the_data%p(4)%P(1:3))
            if ( dist2 - dist < 0.d0 ) then
                nearest => node
                dist = dist2
            end if
            if ( ASSOCIATED(far) ) then
                call nearest_search(Tar, far, nearest, i)
            end if
        endif
    elseif(split == 3)then
        cy = Tar(1)
        cz = Tar(2)
        r  = sqrt(dist*dist - dist1*dist1)
        y1 = node%box(1)
        z1 = node%box(2)
        y2 = node%box(4)
        z2 = node%box(5)
        y0 = (node%box(1) + node%box(4))/2
        z0 = (node%box(2) + node%box(5))/2
        SSC= SphereSplitsurfaceCross(cy,cz,y1,z1,y2,z2,y0,z0,r)
        if(dist > abs(dist1) .and. SSC)then
        dist2 = distance(Tar, node%the_data%p(4)%P(1:3))
        if ( dist2 - dist < 0.d0 ) then
            nearest => node
            dist = dist2
        end if
        if ( ASSOCIATED(far) ) then
            call nearest_search(Tar, far, nearest, i)
        end if
        endif
    end if

    end subroutine nearest_search
!----------------------------------------------------------------------
    logical function SphereSplitsurfaceCross(cy,cz,y1,z1,y2,z2,y0,z0,r)
        real(R8)    :: cy,cz,y1,z1,y2,z2,y0,z0,r
        real(R8)    :: minx,miny

        minx = min(abs(y1-cy),abs(y2-cy))
        miny = min(abs(z1-cz),abs(z2-cz))
        if(minx*minx+miny*miny<r*r)then
            SphereSplitsurfaceCross= .true.
            return
        endif

        if((abs(y0-cy)<(abs(y2-y1)/2+r)) .and. abs(z0-cz)<abs(z2-z1)/2)then
            SphereSplitsurfaceCross= .true.
            return
        endif

        if((abs(z0-cz)<(abs(z2-z1)/2+r)) .and. abs(y0-cy)<abs(y2-y1)/2)then
            SphereSplitsurfaceCross= .true.
            return
        endif
        SphereSplitsurfaceCross= .false.
        return
    end function SphereSplitsurfaceCross
!----------------------------------------------------------------------
    function distance(a, b)
      real(R8)                                    :: distance
        ! type(typPoint)                                   :: a, b
        real(R8)                                    :: a(3), b(3)
        integer                                       :: i

        distance = 0.d0
        do i = 1, 3
            distance = distance + (a(i)-b(i)) * (a(i)-b(i))
        end do
        distance = dsqrt(distance)
        return

    end function distance
!----------------------------------------------------------------------
    real(R8) function dot(a, b)
        type(typPoint)              :: a, b
        dot =  a%P(1)*b%P(1) + a%P(2)*b%P(2) + a%P(3)*b%P(3)
        return
    end function dot
!----------------------------------------------------------------------
    type(typPoint) function rot(a, b)
        type(typPoint)              :: a, b
        rot%P(1) =  a%P(2)*b%P(3) - a%P(3)*b%P(2)
        rot%P(2) =  a%P(3)*b%P(1) - a%P(1)*b%P(3)
        rot%P(3) =  a%P(1)*b%P(2) - a%P(2)*b%P(1)
        return
    end function rot
    end module geometry_mod2
!======================================================================
!   1. 3D Geometry read from POINTWISE file (*.facet) and change into TECPLOT file (*.plt)
!   2. store Geometry surface information into k-D trees
!   3. Find nearest triangle in the Geometry surface for an arbitrary point in 3D space
!   Copyright@BI LIN��2019.06.21
!======================================================================
    subroutine ReadFacet
    use ModPrecision
    use ModInpGlobal
    use ModGeometry
    use geometry_mod2
    use ModTypDef
    use ModKDTree
    implicit none

    integer                     :: i, j, k, l
    integer                     :: vertex(3), ng, nsp, nse, a
    character(80)               :: FileName
    character(20)               :: dd
    logical                     :: file_exist
    type(typPoint), allocatable :: sp(:)
    type(typKDTtree), pointer   :: tp => null()
    real(R8)                    :: box(6), aa, bb
    real(R8), allocatable       :: dis(:, :, :)
    type(typPoint), allocatable :: att(:, :, :), att2
    type(KDT_node), pointer     :: nearest
    real(R8)                    :: tstart,tend

    ! write(*,"(a, I3)") ' Number of seperated objects in this simulation is: ', nGeometry
    ! write(*,*)'-----------------------------------------------------------------------------'
    ! write(*,"(a, 1x, a, 1x, a, 1x)") '| Name of Objects |','Number of surface points |', 'Number of surface elements |'
    ! write(*,"(a, a)")'|-----------------|--------------------------|----------------------------', '|'

    allocate(body(nGeometry), kdtree(nGeometry))

    open(22, file='./Data/Geometry.dat', status='unknown', form='formatted')

    do ng = 1, nGeometry
        if (ng == 0) then
        elseif (ng .LT. 10) then
            write(dd,'(g0)') ng
            dd = '00'//adjustl(dd)
        elseif (ng .LT. 100) then
            write(dd,*) ng
            dd = '0'//adjustl(dd)
        else
            write(dd,*) ng
            dd = adjustl(dd)
        end if

        FileName = trim(GeometryName)//trim(dd)//'.facet'
        print*,"Read file: ", FileName
        inquire( file = FileName , exist = file_exist )
        if(.not. file_exist) stop "Error====> Error geometry file name"

        open(1, file=FileName, status='old', form='formatted')
        read(1,'(///)')
        read(1,*) nsp
        allocate(sp(nsp))
        body(ng)%nsp = nsp
        do i = 1, body(ng)%nsp
            read(1,*) sp(i)%P(1), sp(i)%P(2), sp(i)%P(3)
            sp(i)%label = i
            if ( i == 1) then
                box(1) = sp(i)%P(1)
                box(2) = sp(i)%P(2)
                box(3) = sp(i)%P(3)
                box(4) = sp(i)%P(1)
                box(5) = sp(i)%P(2)
                box(6) = sp(i)%P(3)
            else
                box(1) = min(box(1), sp(i)%P(1))
                box(2) = min(box(2), sp(i)%P(2))
                box(3) = min(box(3), sp(i)%P(3))
                box(4) = max(box(4), sp(i)%P(1))
                box(5) = max(box(5), sp(i)%P(2))
                box(6) = max(box(6), sp(i)%P(3))
            end if
        end do
        read(1,*)
        read(1,*)
        read(1,*) nse
        allocate(body(ng)%se3d(nse))
        body(ng)%nse = nse
        do i = 1, body(ng)%nse
            read(1,*) vertex(1), vertex(2), vertex(3)
            body(ng)%se3d(i)%p(1) = sp(vertex(1))
            body(ng)%se3d(i)%p(2) = sp(vertex(2))
            body(ng)%se3d(i)%p(3) = sp(vertex(3))
            body(ng)%se3d(i)%p(4)%P(:) = (body(ng)%se3d(i)%p(1)%P(:) + &
                                          body(ng)%se3d(i)%p(2)%P(:) + &
                                          body(ng)%se3d(i)%p(3)%P(:))/3.
            body(ng)%se3d(i)%p(4)%label = i
        end do
        body(ng)%box(:) = box(:)

        write(22, *)'VARIABLES="X",        "Y",        "Z"'
        write(22, *)'ZONE N=', nsp, 'E=', nse, 'F=FEPOINT,', 'ET=TRIANGLE'
        do i = 1, nsp
            write(22,*) sp(i)%P(1), sp(i)%P(2), sp(i)%P(3)
        end do
        do i = 1, nse
            write(22,*) body(ng)%se3d(i)%p(1)%label, body(ng)%se3d(i)%p(2)%label, body(ng)%se3d(i)%p(3)%label
        end do

        ! write(*,"(a, 5x, a, 5x, a, 3x, a, 11x, a, 6x, a, 10x, a)") '|', head//dd, '|', nsp, '|', nse, '|'
        ! if (ng .NE. nGeometry) then
        !     write(*,"(a, a)")' ----------------- -------------------------- ----------------------------', ' '
        ! else
        !     write(*,*)'-----------------------------------------------------------------------------'
        ! end if

        close(1)
        deallocate(sp)
    end do
    close(22)

    nsp = 0
    nse = 0
    do ng = 1, nGeometry
        nsp = nsp + body(ng)%nsp
        nse = nse + body(ng)%nse
    end do
    write(*,*) 'Total surface   points number of all objects is : ', nsp
    write(*,*) 'Total surface elements number of all objects is : ', nse

    tp => kdtree(1)
    call CPU_TIME(tStart)
    call create_KDT_tree_for_body_i(tp, 1)
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "creat tree time:",tEnd-tStart

    call CPU_TIME(tStart)
    call KDTree_out(tp%root)
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "output tree time:",tEnd-tStart
    DEALLOCATE(body)

    return
    end subroutine ReadFacet
!======================================================================
    subroutine ReadStl
    use ModPrecision
    use ModInpGlobal
    use ModGeometry
    use geometry_mod2
    use ModTypDef
    use ModKDTree
    implicit none
    integer                     :: i, j, k
    integer                     :: ng  ! Geometry number
    integer                     :: nsp ! Point number
    integer(I4)                 :: nse ! Element number
    character(80)               :: FileName
    character(20)               :: dd
    logical                     :: file_exist
    type(typPoint), allocatable :: sp(:)
    type(typKDTtree), pointer   :: tp => null()
    real(R8)                    :: box(6)
    character(80)               :: Header
    real(R4)                    :: Normal(3)
    real(R4)                    :: Vertex(3)
    integer(I2)                 :: unknown
    real(R8)                    :: tstart,tend
    real(R8)                    :: tmp(3)

    allocate(body(nGeometry), kdtree(nGeometry))

    open(22, file='./Data/Geometry.dat', status='unknown', form='formatted')

    do ng = 1, nGeometry
        if (ng == 0) then
        elseif (ng .LT. 10) then
            write(dd,'(g0)') ng
            dd = '00'//adjustl(dd)
        elseif (ng .LT. 100) then
            write(dd,*) ng
            dd = '0'//adjustl(dd)
        else
            write(dd,*) ng
            dd = adjustl(dd)
        end if

        FileName = trim(GeometryName)//trim(dd)//'.stl'
        print*,"Read file: ", FileName
        inquire( file = FileName , exist = file_exist )
            if(.not. file_exist) stop "Error====> Error geometry file name"
        open(1, file=FileName, status='old', form='unformatted', &
                access='stream')

        read(1) Header
        print*,'Header=',Header
        read(1) nse
        print*,'elements=',nse

        allocate(body(ng)%se3d(nse))
        body(ng)%nse = nse
        body(ng)%nsp = nse*3

        read(1) (Normal(k), k=1,3)
        do i = 1, 3
            do k = 1, 3
                read(1) Vertex(k)
                body(ng)%se3d(1)%P(i)%p(k) = Vertex(k)
            enddo
        enddo
        read(1) unknown
        box(1) = min(body(ng)%se3d(1)%p(1)%P(1), &
                     body(ng)%se3d(1)%p(2)%P(1), &
                     body(ng)%se3d(1)%p(3)%P(1))
        box(2) = min(body(ng)%se3d(1)%p(1)%P(2), &
                     body(ng)%se3d(1)%p(2)%P(2), &
                     body(ng)%se3d(1)%p(3)%P(2))
        box(3) = min(body(ng)%se3d(1)%p(1)%P(3), &
                     body(ng)%se3d(1)%p(2)%P(3), &
                     body(ng)%se3d(1)%p(3)%P(3))
        box(4) = box(1)
        box(5) = box(2)
        box(6) = box(3)
        body(ng)%se3d(1)%p(4)%P(:) = (body(ng)%se3d(1)%p(1)%P(:) + &
                                      body(ng)%se3d(1)%p(2)%P(:) + &
                                      body(ng)%se3d(1)%p(3)%P(:))/3.
        body(ng)%se3d(1)%p(4)%label = 1
        if (checkvec(body(ng)%se3d(1),real(Normal,R8))) then
            tmp(:) = body(ng)%se3d(1)%p(2)%P(:)
            body(ng)%se3d(1)%p(2)%P(:) = body(ng)%se3d(1)%p(3)%P(:)
            body(ng)%se3d(1)%p(3)%P(:) = tmp(:)
        endif

        do i = 2, nse
            read(1) (Normal(k), k=1,3)
            do j = 1, 3
                do k = 1, 3
                    read(1) Vertex(k)
                    body(ng)%se3d(i)%P(j)%p(k) = Vertex(k)
                enddo
                box(1) = min(box(1), body(ng)%se3d(i)%p(j)%P(1))
                box(2) = min(box(2), body(ng)%se3d(i)%p(j)%P(2))
                box(3) = min(box(3), body(ng)%se3d(i)%p(j)%P(3))
                box(4) = max(box(4), body(ng)%se3d(i)%p(j)%P(1))
                box(5) = max(box(5), body(ng)%se3d(i)%p(j)%P(2))
                box(6) = max(box(6), body(ng)%se3d(i)%p(j)%P(3))
            enddo
            read(1) unknown
            body(ng)%se3d(i)%p(4)%P(:) = (body(ng)%se3d(i)%p(1)%P(:) + &
                                          body(ng)%se3d(i)%p(2)%P(:) + &
                                          body(ng)%se3d(i)%p(3)%P(:))/3.
            body(ng)%se3d(i)%p(4)%label = i
            if (checkvec(body(ng)%se3d(i),real(Normal,R8))) then
                tmp(:) = body(ng)%se3d(i)%p(2)%P(:)
                body(ng)%se3d(i)%p(2)%P(:) = body(ng)%se3d(i)%p(3)%P(:)
                body(ng)%se3d(i)%p(3)%P(:) = tmp(:)
            endif
        end do
        body(ng)%box(:) = box(:)

        write(22, *)'VARIABLES="X",        "Y",        "Z"'
        write(22, *)'ZONE N=', body(ng)%nsp, 'E=', nse
        write(22, *)'F=FEPOINT,', 'ET=TRIANGLE'
        do i = 1, nse
            do j = 1, 3
                write(22,*) (body(ng)%se3d(i)%p(j)%P(k), k=1,3)
            enddo
        end do
        do i = 0, nse*3, 3
            write(22,*) i+1, i+2, i+3
        end do
        close(1)
    end do
    close(22)

    nsp = 0
    nse = 0
    do ng = 1, nGeometry
        nsp = nsp + body(ng)%nsp
        nse = nse + body(ng)%nse
    end do
    write(*,*) 'Total surface   points number of all objects is : ', nsp
    write(*,*) 'Total surface elements number of all objects is : ', nse


    tp => kdtree(1)
    call CPU_TIME(tStart)
    call create_KDT_tree_for_body_i(tp, 1)
    call CPU_TIME(tEnd)
    write(*,'(1X,A,F10.2)') "creat tree time:",tEnd-tStart

    ! call CPU_TIME(tStart)
    ! call KDTree_out(tp%root)
    ! call CPU_TIME(tEnd)
    ! write(*,'(1X,A,F10.2)') "output tree time:",tEnd-tStart

    DEALLOCATE(body)
    return
    contains
!----------------------------------------------------------------------
        logical function checkvec(tri,Nin)
        ! T: Inconsistent
        use ModPrecision
        use ModTypDef
        use ModTools, only: CROSS_PRODUCT_3
        implicit none
        type(triangle)          :: tri
        real(R8),INTENT(IN)     :: Nin(3)
        real(R8)                :: Nnow(3)
        real(R8)                :: n(3), v1(3), v2(3)
        real(R8)                :: a

        checkvec = .false.
        v1 = tri%P(2)%P - tri%P(1)%P
        v2 = tri%P(3)%P - tri%P(2)%P
        Nnow  = CROSS_PRODUCT_3(v1,v2)
        a  = DOT_PRODUCT(Nin,Nnow)

        ! a = 0
        ! if (Nnow(1)>=0.and.Nin(1)>=0)then
        ! continue
        ! elseif (Nnow(1)<0.and.Nin(1)<0)then
        ! continue
        ! else
        !     a = a
        ! endif
        ! if (Nnow(2)>=0.and.Nin(2)>=0)then
        ! continue
        ! elseif (Nnow(2)<0.and.Nin(2)<0)then
        ! continue
        ! else
        ! continue
        !     a = a
        ! endif
        ! if (Nnow(3)>=0.and.Nin(3)>=0)then
        ! continue
        ! elseif (Nnow(3)<0.and.Nin(3)<0)then
        ! continue
        ! else
        !     a = a
        ! endif

        if (a<=0) then
            print*,'Inconsistent geometry vector has modified.'
            checkvec = .true.
        endif
        endfunction checkvec
    endsubroutine ReadStl
!======================================================================
    subroutine initCheckVector
    ! use ModPrecision
    use ModInpGlobal,only: nGeometry
    use ModKDTree
    implicit none
    integer                     :: i

    do i= 1, nGeometry
        call CheckVector(kdtree(i)%root)
    enddo
    contains
!----------------------------------------------------------------------
        recursive subroutine CheckVector(kdnode)
        use ModPrecision
        use ModGeometry
        use ModTypDef
        use ModTools, only : CROSS_PRODUCT_3
        implicit none
        type(KDT_node),pointer      :: kdnode
        real(R8)                    :: n(3), v1(3), v2(3)
        real(R8)                    :: p(3), p1(3)
        real(R8)                    :: dp, dp1

        if (.not.ASSOCIATED(kdnode)) return

        p1 = kdnode%the_data%P(1)%P
        v1 = kdnode%the_data%P(2)%P - kdnode%the_data%P(1)%P
        v2 = kdnode%the_data%P(3)%P - kdnode%the_data%P(2)%P
        n  = CROSS_PRODUCT_3(v1,v2)
        p  = p1 + n ! p(x,y,z)=p1(x,y,z)+n(x,y,z)
        dp = sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3)) ! Distance from np to orgin.
        dp1= sqrt(p1(1)*p1(1)+p1(2)*p1(2)+p1(3)*p1(3))
        if (dp<dp1) print*,"error geometry vector"
        call CheckVector(kdnode%left)
        call CheckVector(kdnode%right)

        endsubroutine CheckVector
    endsubroutine initCheckVector
!======================================================================
