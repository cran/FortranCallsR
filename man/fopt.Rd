\name{toyOptim}


\alias{toyOptim}


\title{Toy Model Optimizer}


\description{

    This is an educational package which shows a simple solution 
    how to implement a function optimizer written in Fortran into 
    the R environment for box bounded one-dimensional optimization 
    calling the objective function from R.
    
    The example can be used as a prototype template for implementing
    Fortran subroutines which need to access a R function from R.
            
}


\usage{
toyOptim(start, objective, lower = -Inf, upper = Inf, trace = TRUE, 
    control = list(), \dots)
}


\arguments{
     
    \item{start}{
        a numerical value, a guess where the optimizer should start.\cr
        (In this toy example you can give any real number, it will 
        be overwritten by an optimal choice of the solver.)
        }
    \item{objective}{
        the objective function.
        }
    \item{lower, upper}{
        lower and upper values of the searchable interval.\cr
        (In this toy example the default values are -1e99 and 1e99 
        if lower and upper are specified as -Inf and +Inf.)
        }
    \item{trace}{
        a logical flag, should the optimization be traced?
        }
    \item{control}{
        list of control parameters.\cr
        (In this toy example there is only one control parameter,
        named \code{tol}. Its default value is 1.0e-6.)
        }
    \item{\dots}{
        optional arguments passed to the objective function.
        }
        
}


\details{
  
    \bold{The Fortran Solver: File src/foptsrc.f}
    
    To demonstrate how to call a R function from Fortran we consider a 
    toy model for a box bounded one-dimensioanl optimizer which requires 
    only a function call to the objective function fun(n, x, value). The
    optimzer itself is written in Fortran, the objective function will be
    passed from the R environment via a C wrapper.

    So this example can be considered as a starting point and prototype
    template for similar tasks. These are tasks where we have Fortran 
    Code which optimizes, finds zeroes, integrates or even is good for 
    other problems, and the objective function should be made available 
    through the R environment.

    The toy model considered here is Brent's one dimensional optimzer (1973). 
    We have slightly modified the original code which we copied from the
    GAMS 'golden oldies' section. The optimizer and objective function are 
    Fortran subroutines. Note, that the strategy presented here can be 
    used to implement other Fortran Routines in the same way.
    
    
    \bold{The Fortran-C Interface: File src/fopt.c}
    
    This is the Fortran-C Interface for the toy model representing Brent's
    optimizer written in Fortran. The Fortran subroutine \code{fun()} calls \
    the C function \code{cfun()} which makes the function call to the R 
    function through the call \code{eval(objective, environment)}. Here
    \code{objective} and \code{environment} are global SEXP variables 
    (although we are aware of the fact that global variables may not be
    an optimal solution).

    The SEXP C function \code{mfopt()} is the main C function which is 
    called from the R function \code{mfopt()} through  
    \code{.Call("mfopt", fn, rho)}, where \code{fn} and \code{rho} 
    denote the objective function and the embedding environment.
    
    
    \bold{The R Function Script: File R/fopt.R}
    
    The R function \code{toyOptim()} runs the optimization. The argument
    \code{start} holds a guess where to start the optimization, 
    \code{objective} defines the objective written in R, \code{lower}
    and \code{upper} define the interval bounds as constrains, \code{trace}
    allows to trace the iteration path, and the \code{control} list
    may define optional control parameters required by the solver. 

}


\value{

    The function returns a list with the results obtained by the solver.
    \cr 
    
    For our toy model this is the parameter \emph{x}, named \code{$par},
    where the function \code{f(x)}, named \code{$objective} takes its 
    minimum, the value of the objective function at the minimum, and 
    the tolerance, named \code{$tol}.
    
}


\note{

    Note, the prototype template presented here uses global SEXP 
    variables for the function to be accessed from R and the embedding
    environment. Another example for a Fortran/R interface is the
    \code{nlminb} optimizer using a different approach named "reverse 
    communications". We recommend to inspect the source code. 
    Unfortunately, most of the Fortran subroutines don't allow to apply 
    this technique in a straightforward way, and then the most simple 
    approach is the method to work with global variables. Furthermore, 
    we like to remark that a (minor) drawback of the use of global 
    variables ist that it will not work in (the very seldom case of) 
    multithreaded systems.

}


\author{

    Diethelm Wuertz for this R port, \cr
    Richard Brent for the box-bounded one-dimensional Fortran solver FMIN.
    
}


\examples{
## toyOptim - 
   
   # Simple Test Function 1:
   toyFun1 = function(x) { 
     fun = x^2 
     # .trace is a global variable 
     if (.trace) print(c(x = x, fun = fun))
     fun
   }
   toyOptim(start = 0, objective = toyFun1, lower = -1, upper = 1)
   
   
   # Simple Test Function 2 - with addtitional Parameter a:
   toyFun2 = function(x, a = 1) { 
     fun = -sin(a*x) 
     if (.trace) print(c(x = x, fun = fun))
     fun
   } 
   toyOptim(start = 0, objective = toyFun2, lower = 0, upper = pi)
}


\keyword{programming}

