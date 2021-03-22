!======================================================================
!
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
        logical :: NRR
        logical :: Chimera
        integer :: Limiter
        logical :: Restart
        character(50):: OutputNameStr
        character(50):: GeometryName
        character(10):: OutputFormat
    end module ModInpGlobal
!======================================================================
    module ModInpMesh
        use ModPrecision
        implicit none
        integer :: nCell(3)
        real(R8):: Domain1(3), Domain2(3)
        integer :: InitRefineLVL
        integer :: AdaptRefineLVL
        integer :: cIntersectMethod
    end module ModInpMesh
!======================================================================
    module ModInpInflow
        use ModPrecision
        implicit none
        real(R8):: Alpha
        real(R8):: Beta
        real(R8):: ReyNum
        real(R8):: T00
        real(R8):: Mach00
        real(R8):: Gama00
        real(R8):: Rgas    ![J/(kg.K)]
    end module ModInpInflow
!======================================================================
    module ModInpNRRset
        use ModPrecision
        implicit none
        integer :: NRRSeed
        real(R8):: NRRLength
        real(R8):: NRRTheta
    end module ModInpNRRset
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
    
! Cartesian Cell's point coordinate
        type typPoint
            real(R8):: P(3)
        end type typPoint
    
! ! Triangle's three vertexs��odered by right-hand rule
!         type triangle
!             type(typPoint)::tripoint1, tripoint2, tripoint3
!         end type triangle
    
! Cartesian grid data structure
! nBGCell    = number of the back-ground cells (root node) 
! nCell      = number of the cells (leaf node) 
! levelx     = x-direct level
! levely     = y-direct level
! levelz     = z-direct level
! cross      relationship between cell and the object surface.
!            = 0 Cell outside the object surface.
!            = 1 Intersect while cell center outside the object surface
!            = 2 Intersect while cell center inside the object surface
!            = 3 Cell inside the object surface
!            = -1 NRR: ray region
!            = -2 NRR: region between ray-ray
! fSplitType = 0 Cell is isotropical 
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
! Center     1 x
!            2 y
!            3 z
! Neighbor   1 plus
!            2 minus
        type typCell
            integer :: nBGCell, nCell
            integer :: lvl(3)
            integer :: cross
            integer :: fSplitType, Location
            integer :: Node(8)
            real(R8):: Center(3)
            real(R8):: U(5)
            type(typCell),pointer :: Father
            type(typCell),pointer :: son1, son2, son3, son4,    &
                                     son5, son6, son7, son8
            type(typCell),pointer :: NeighborX1, NeighborX2,    &
                                     NeighborY1, NeighborY2,    &
                                     NeighborZ1, NeighborZ2

        end type typCell
    end module ModTypDef
!======================================================================
! Define globally share constants
    module ModGlobalConstants
        use ModPrecision
        implicit none
! constants used in sutherlan'law, =110.3 in NSMB5.0
!        real(R8),parameter:: C00=100.4
!        real(R8),parameter:: gama=1.4, gama1=gama-1.0
!        real(R8),parameter:: pr=0.72, prt=0.9
        real(R8),parameter:: PI=3.14159265358979
        integer,parameter:: SchemeNND2=1,SchemeWENO3=3
        integer,parameter:: FluxRoe=1,Fluxcentral=2
        integer,parameter:: BCWall=2, BCSymmetry=3, BCFarfield=4
        integer,parameter:: TimeRK3=1,TimeLUSGS=0
        integer,parameter:: TurSA=1,TurSST=2,TurKW=3  
    end module ModGlobalConstants
!======================================================================
!======================================================================
!======================================================================