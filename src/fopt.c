/*
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA 02111-1307 USA


  This is the C-Fortran Interface for a toy model representing a box bounded  
  one-dimensioanl optimizer written in Fortran which requires only a function 
  call to the objective function, the Fortran subroutine fun(n, x, value). 
  The Fortran subroutine fun() calls the C function cfun() which makes the
  function call to the R function through 'eval(objective, environment)'. Here
  'objective' and 'environment are global SEXP variables (although I know that
  global variables are not an optimal solution).

  The SEXP C function \code{mfopt()} is the main C function which is called 
  from the R function \code{mfopt()} through '.Call("mfopt", fn, rho)', where 
  'fn' and 'rho' denote the objective function and the embedding environment.

  Diethelm Wuertz
  ETH Zurich
  www.rmetrics.org
*/


#include <Rinternals.h>
#include <R_ext/RS.h>


SEXP objective ;
SEXP environment ;


/*
--------------------------------------------------------------------------------
*/


void cfunc(int *n, double x[], double value[])
{
    SEXP PAR ;
    int i ;
    PROTECT(PAR = findVarInFrame(environment, install(".par"))) ;
    for (i = 0; i < *n; i++) REAL(PAR)[i] = x[i] ;
    value[0] = asReal(eval(objective, environment)) ;
    UNPROTECT(1) ;
}


/*
--------------------------------------------------------------------------------
*/


void F77_SUB(cfunc)(int *n, double x[], double value[])
{ 
    cfunc(n, x, value) ;    
}  


/*
--------------------------------------------------------------------------------
*/


SEXP mfopt(SEXP fn, SEXP rho)
{
    
    SEXP PAR, XLO, XUP, VAL, RPAR, N;
   
    objective = fn ;
    environment = rho ;
    
    PROTECT(PAR  = findVarInFrame(rho, install(".par"))) ;
    PROTECT(XLO  = findVarInFrame(rho, install(".xlo"))) ;
    PROTECT(XUP  = findVarInFrame(rho, install(".xup"))) ;
    PROTECT(VAL  = findVarInFrame(rho, install(".val"))) ;
    PROTECT(RPAR = findVarInFrame(rho, install(".rpar"))) ;
    PROTECT(N    = findVarInFrame(rho, install(".n"))) ;
        
    /* 
    Call the Fortran Solver:
    */
    
    F77_CALL(fopt)(REAL(PAR), INTEGER(N), REAL(XLO), REAL(XUP), REAL(VAL), 
        INTEGER(N), REAL(RPAR), INTEGER(N)) ;
    
    UNPROTECT(6) ;
    return R_NilValue;
}


/*
--------------------------------------------------------------------------------
*/


