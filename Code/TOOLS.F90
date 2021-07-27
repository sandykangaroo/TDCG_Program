!======================================================================
    module ModTools
    contains

        PURE function CROSS_PRODUCT_3(a,b)
        ! Back cross-product result for a(3) and b(3).
        use ModPrecision
        implicit none
        real(R8)           :: CROSS_PRODUCT_3(3)
        real(R8),INTENT(IN):: a(3), b(3)
        CROSS_PRODUCT_3 = [a(2)*b(3)-a(3)*b(2), &
                        a(3)*b(1)-a(1)*b(3), &
                        a(1)*b(2)-a(2)*b(1)]
        endfunction CROSS_PRODUCT_3
!----------------------------------------------------------------------
        PURE function PtF_normal(p,tri) ! Point to Face normal.
        ! PtF_Dist: normal vector from p -> tri.
        !           Note not the vector tri -> p.
        use ModPrecision
        use ModTypDef, only: triangle
        implicit none
        real(R8)                :: PtF_normal(3)
        real(R8),INTENT(IN)     :: p(3)
        type(triangle),pointer  :: tri
        real(R8)                :: n(3) ! TriFace normal vector.
        real(R8)                :: v0(3), v1(3), v2(3), dist
        real(R8)                :: module_v0, module_n, cos_alph

        V1 = tri%P(2)%P - tri%P(1)%P
        V2 = tri%P(3)%P - tri%P(2)%P
        n  = CROSS_PRODUCT_3 (V1, V2)
        v0 = p - tri%P(4)%P
        module_n = sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
        dist = DOT_PRODUCT(v0, n)/module_n
        PtF_normal = -n / module_n * dist

        endfunction PtF_normal
!----------------------------------------------------------------------
        function BBOX(p,box)
        ! =-1 outside bbox; = 0 on the box; = 1 inside bbox.
        use ModPrecision
        implicit none
        integer(I2)         :: BBOX
        real(R8),INTENT(IN) :: p(3)
        real(R8),INTENT(IN) :: box(6)

        if (p(1) < box(1) .or. p(2) < box(2) .or. p(3) < box(3) .or. &
            p(1) > box(4) .or. p(2) > box(5) .or. p(3) > box(6)) then
            BBOX = -1
        elseif (p(1)==box(1) .or. p(2)==box(2) .or. p(3)==box(3) .or. &
                p(1)==box(4) .or. p(2)==box(5) .or. p(3)==box(6)) then
            BBOX = 0
        else
            BBOX = 1
        endif
        endfunction BBOX
!----------------------------------------------------------------------
        function CrossPoint(Point,D,tri)
        use ModPrecision
        use ModTypDef, only: triangle
        implicit none
        real(R8)                 :: CrossPoint(4)
        real(R8)      ,INTENT(IN):: Point(3)   ! point(3) = x, y, z.
        real(R8)      ,INTENT(IN):: D(3)
        type(triangle),INTENT(IN):: tri
        real(R8):: V0(3), V1(3), V2(3)
        real(R8):: P(3), Q(3)
        real(R8):: t, u, v
        real(R8):: A

        V0 = Point      - tri%P(1)%P
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
        if (A>0 .and. u>=0 .and. v>=0 .and. u+v<=A .and. &
            t>=0 .and. t<=A) then
            CrossPoint(1:3) = point(:) + D(:) * t / A
            CrossPoint(4)   = 1
            return
        elseif (A<0 .and. u<=0 .and. v<=0 .and. u+v>=A .and. &
            t<=0 .and. t>=A) then
            CrossPoint(1:3) = point(:) + D(:) * t / A
            CrossPoint(4)   = 1
            return
        else
            CrossPoint(:)   = (/0,0,0,0/)
            return
        endif
        endfunction CrossPoint
!----------------------------------------------------------------------
        PURE function MollerTrumbore(Point,D,tri,ios)
        use ModPrecision
        use ModTypDef
        implicit none
        logical                  :: MollerTrumbore
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
        else
            if (A>0 .and. u>=0 .and. v>=0 .and. u+v<=A .and. t>=0) then
                MollerTrumbore=.true.; return
            elseif (A<0 .and. u<=0 .and. v<=0 .and. u+v>=A .and. t<=0) then
                MollerTrumbore=.true.; return
            else
                MollerTrumbore=.false.; return
            endif
        endif
        endfunction MollerTrumbore
!----------------------------------------------------------------------
    endmodule ModTools
!======================================================================
!======================================================================
