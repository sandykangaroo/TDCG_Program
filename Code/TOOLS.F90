!======================================================================
    module ModTools
    contains

    function CROSS_PRODUCT_3(a,b)
    use ModPrecision
    implicit none
    real(R8)           :: CROSS_PRODUCT_3(3)
    real(R8),INTENT(IN):: a(3), b(3)
    CROSS_PRODUCT_3 = [a(2)*b(3)-a(3)*b(2), &
                       a(3)*b(1)-a(1)*b(3), &
                       a(1)*b(2)-a(2)*b(1)]
    endfunction CROSS_PRODUCT_3
    endmodule ModTools
!======================================================================
