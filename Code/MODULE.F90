!======================================================================
    module ModPrecision
        implicit none
        integer,parameter::R4=4
        integer,parameter::R8=8
        integer,parameter::I1=1
        integer,parameter::I2=2
        integer,parameter::I4=4
    end module ModPrecision
!======================================================================
!======================================================================
    module ModInpGlobal
        use ModPrecision
        implicit none
        logical :: Aniso
        integer :: nStep
        integer :: nSave
        integer :: nAdaptStep
        integer :: Debug
        Real(R8):: CFL
        logical :: MeshOnly
        logical :: NRR
        logical :: Chimera
        integer :: Limiter
        logical :: Restart
        integer :: DimensionALL ! 20, 2D; 30, 3D
        integer :: nGeometry ! Number of body
        character(50):: GeometryName
        character(50):: GeometryFormat
        character(50):: OutputName
        character(10):: OutputFormat
    end module ModInpGlobal
!----------------------------------------------------------------------
    module ModInpMesh
        use ModPrecision
        implicit none
        integer :: nCell(3)
        real(R8):: DomainMin(3), DomainMax(3)
        integer :: InitRefineLVL
        integer :: AdaptRefineLVL
        integer :: cIntersectMethod
        logical :: useKDT
    end module ModInpMesh
!----------------------------------------------------------------------
    module ModInpInflow
        use ModPrecision
        implicit none
        real(R8):: Alpha
        real(R8):: Beta
        real(R8):: Re
        real(R8):: T00
        real(R8):: Ma00
        real(R8):: Gama00
        real(R8):: Rgas    ![J/(kg.K)]
    end module ModInpInflow
!----------------------------------------------------------------------
    module ModInpNRRset
        use ModPrecision
        implicit none
        integer :: NRRSeed
        real(R8):: NRRLength
        real(R8):: NRRTheta
    end module ModInpNRRset
!======================================================================
    module ModGlobal
        use ModPrecision
        implicit none
        integer :: ProcessLocation
    end module ModGlobal
!======================================================================
! Define dynamic data structure
    module ModTypDef
        use ModPrecision
        implicit none
! One dimension dynamic real array
!         type Real1DArray
!             real(R8)::R1D(:)
!         end type Real1DArray

! ! Two dimension dynamic real array
!         type Real2DArray
!             real(R8)::R2D(:,:)
!         end type Real2DArray

! ! One dimension dynamic integer array
!         type Integer1DArray
!             integer::I1D(:)
!         end type Integer1DArray

! ! Two dimension dynamic integer array
!         type Integer2DArray
!             integer::I2D(:,:)
!         end type Integer2DArray

! Cartesian OctCell's point coordinate
        type typPoint
            real(R8):: P(3)
            integer :: label
        end type typPoint

        type segment
            type(typPoint):: P(2)  !p(1) = p1, p(2) = p2
        end type segment

        type triangle
            type(typPoint):: P(4)
            !p(1) = p1, p(2) = p2, p(3) = p3, p(4) = center
        end type triangle

! Cartesian grid data structure
        ! nBGCell    = number of the back-ground cells (root node)
        ! nCell      = number of the cells (leaf node)
        ! levelx     = x-direct level
        ! levely     = y-direct level
        ! levelz     = z-direct level
        ! cross      relationship between Cell and the object surface.
        !            = -5 Initial
        !            = -4 Not intersect Cell.
        !            = -3 Intersect Cell.
        !            =  0 Cell outside the object surface.
        !            =  1 Intersect while Cell center outside the object surface
        !            =  2 Intersect while Cell center inside the object surface
        !            =  3 Cell inside the object surface
        !            = -1 NRR: ray region
        !            = -2 NRR: region between ray-ray
        ! fSplitTyp = 0 Cell is isotropical
        !            = 1 Cell is obtained by refined in x-direction
        !            = 2 Cell is obtained by refined in y-direction
        !            = 3 Cell is obtained by refined in z-direction
        !            = 4 Cell is obtained by refined in xy-direction
        !            = 5 Cell is obtained by refined in xz-direction
        !            = 6 Cell is obtained by refined in yz-direction
        ! Location   = 0 Cell does not have a son
        !            = 1,2,3,...8 Cell's location among siblings
        ! Node       number of the eight vertexs
        ! U          1 rou*u
        !            2 rou*v
        !            3 rou*w
        !            4 rou
        !            5 T
        ! Mark       Mark in subroutine SmoothMesh
        ! Center     1 x
        !            2 y
        !            3 z
        ! Neighbor   1 Xminus
        !            2 Xplus
        !            3 Yminus
        !            4 Yplus
        !            5 Zminus
        !            6 Zplus
        ! nCrossTri  Size of CrossTri(:)
        ! CrossTri   KDT pointer for the triangle cross the cell
        type FttCell
            integer(I4)             :: nBGCell(3)
            integer(I4)             :: nCell
            integer(I2)             :: Location
            integer(I2)             :: LVL(3)
            integer(I2)             :: fSplitTyp
            integer(I2)             :: Cross
            logical                 :: Mark(3)
            real(R8)                :: Center(3)
            real(R8)                :: U(5)
            ! real(R8)                :: walldistance
            type(FttCell),pointer   :: Father
            type(FttOct), pointer   :: Octson
            type(tCrossTri),pointer :: CrossTri
        end type FttCell

        type FttOct
            integer(I2)             :: nSon
            type(FttCell),pointer   :: Neighbor1, Neighbor2, Neighbor3, &
                                       Neighbor4, Neighbor5, Neighbor6
            type(FttCell)           :: son(8)
        end type FttOct

        type tCrossTri
            type(KDT_node), pointer :: tri
            type(tCrossTri),pointer :: next
        end type tCrossTri

        ! type BlockCell

        type typStructCell
            integer :: cross
            real(R8):: Center(3)
            real(R8):: U(5)
        end type typStructCell

        type KDT_node
            type(triangle), pointer:: the_data
            integer                :: SplitAxis ! The dimension of split.
                                        !  =1 x_axis; =2 y_axis; =3 z_axis
            real(R8)               :: box(6) ! Bounding box of The_data.
                ! box(1:3) = xmin, ymin, zmin; box(4:6) = xmax, ymax, zmax
            integer                :: level ! depth of KDTree
            type(KDT_node), pointer:: left, right, parent
        end type KDT_node

    end module ModTypDef
!======================================================================
! Define globally share constants
    module ModGlobalConstants
        use ModPrecision
        implicit none
! ! constants used in sutherlan'law, =110.3 in NSMB5.0
! !        real(R8),parameter:: C00=100.4
! !        real(R8),parameter:: gama=1.4, gama1=gama-1.0
! !        real(R8),parameter:: pr=0.72, prt=0.9
!         real(R8),parameter:: PI=3.14159265358979
!         integer,parameter:: SchemeNND2=1,SchemeWENO3=3
!         integer,parameter:: FluxRoe=1,Fluxcentral=2
!         integer,parameter:: BCWall=2, BCSymmetry=3, BCFarfield=4
!         integer,parameter:: TimeRK3=1,TimeLUSGS=0
!         integer,parameter:: TurSA=1,TurSST=2,TurKW=3
        real(R4),parameter :: epsR4 = 0.00001
        real(R8),parameter :: epsR8 = 0.00000000000001
    end module ModGlobalConstants
!======================================================================
    module ModKDTree
        use ModPrecision
        use ModTypDef
        implicit none

        type typKDTtree
            type(KDT_node), pointer:: root=>null()
        end type typKDTtree

        type(typKDTtree), allocatable, target   :: KDTree(:)

    endmodule ModKDTree
!======================================================================
    module ModGeometry
    ! Define geometry discreted points and elements.
        use ModPrecision
        use ModTypDef
        implicit none

        type geom
            integer                          :: nsp, nse
            type(triangle) , allocatable     :: se3d(:)
            real(R8)                         :: box(6)
        end type geom

        type(geom), allocatable, target     :: body(:)

    end module ModGeometry
!======================================================================
    module ModMesh
        use ModPrecision
        use ModTypDef
        implicit none
        integer :: nBGCells     ! Number of the background Cells.
        integer :: nCells     ! Number of total Cells.
        real(R8):: BGCellSize(3)    ! Step size for background OctCell.
        type(FTTCell),pointer :: OctCell(:,:,:)
        ! type(typStructCell),pointer:: StrCell(:,:,:,:)
    endmodule ModMesh
!======================================================================
    module ModSolve
        use ModPrecision
        implicit none
        integer :: step
        real(R8):: TimeStep
    endmodule ModSolve
!======================================================================
    module ModOutput
        use ModPrecision
        use ModTypDef
        implicit none
        real(R8),ALLOCATABLE :: Nodes(:,:)      ! Nodes coordinate
        real(R8),ALLOCATABLE :: cVariables(:,:) ! OctCell variables value
        integer ,ALLOCATABLE :: cNodes(:,:)     ! OctCell nodes number
        integer :: nNodes
    endmodule ModOutput
!======================================================================
!======================================================================

        ! type(FTTCell),pointer :: cs
        ! integer:: is

        ! if(ASSOCIATED(c%Octson))then
        !     do is=1,c%Octson%nSon
        !         cs=>c%Octson%son(is)
        !         call (cs)
        !     enddo
        !     return
        ! endif
