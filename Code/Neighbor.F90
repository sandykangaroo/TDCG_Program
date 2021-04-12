    module ModNeighbor
        use ModMesh
        use ModTypDef 
        use ModInpMesh
        implicit none
    contains

! X- neighbor 
    recursive function NeighborX1(c)
        implicit none
        type(octCell),pointer :: c,cf,cintial,cfn,NeighborX1
    end function NeighborX1
    
! X+ neighbor 
    recursive function NeighborX2(c)
        implicit none
        type(octCell),pointer :: c,cf,cinitial,cfn,NeighborX2 
    end function NeighborX2
! X- neighbor 
    recursive function NeighborY1(c)
        implicit none
        type(octCell),pointer :: c,cf,cintial,cfn,NeighborY1
    end function NeighborY1
    
! X+ neighbor 
    recursive function NeighborY2(c)
        implicit none
        type(octCell),pointer :: c,cf,cinitial,cfn,NeighborY2 
    end function NeighborY2
! X- neighbor 
    recursive function NeighborZ1(c)
        implicit none
        type(octCell),pointer :: c,cf,cintial,cfn,NeighborZ1
    end function NeighborZ1
    
! X+ neighbor 
    recursive function NeighborZ2(c)
        implicit none
        type(octCell),pointer :: c,cf,cinitial,cfn,NeighborZ2 
    end function NeighborZ2
    end module ModNeighbor
