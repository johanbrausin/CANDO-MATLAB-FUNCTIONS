function dy = odeexcCL(t,y,A)
    
%==========================================================================
%AUTHOR : E. Boulais, LCBB MIT
%DATE   : 04/04/2013
%INPUT  : t      : Symbolic time variable
%         y      : Symbolic solution variable
%         A      : Transition matrix
%OUTPUT : dy    ==> Symbolic solution variable
%
%==========================================================================

dy=A*y;
end

