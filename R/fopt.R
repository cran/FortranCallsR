
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


################################################################################
# FUNCTION:         DESCRIPTION:
#  toyOptim          A simple box bounded one-dimensional optimizer
################################################################################


toyOptim = function(start, objective, lower = -Inf, upper = Inf, 
trace = TRUE, control = list(), ...)
{	# A function implemented by Diethelm Wuertz

	# Description:
	#	A simple box bounded one-dimensional optimizer
	
	# Arguments:
	#	start - a numerical value, a guess where the optimizer should 
	#		start. (In this toy example you can give any real number, it 
	#		will be overwritten by an optimal choice of the solver.)
	# 	objective - the objective function.
	#	lower, upper - lower and upper values of the searchable interval.
	#		(In this toy example the default values are -1e99 and 1e99 
	#		if lower and upper are specified as -Inf and +Inf.)
	#	trace - a logical flag, should the optimization be traced?
	#	control - list of control parameters. (In this toy example there 
	#		is only one control parameter, named 'tol'. Its default value 
	#		is 1.0e-6.)
	#	... - optional arguments passed to the objective function.
		
	# FUNCTION:
	
	# Parameter Settings:
	val = 0
	BIG <- 1.e9
	con <- list(tol = 1.0e-8)
	con[(namc <- names(control))] <- control
	.trace <<- trace

	# Environment Settings:
	n <- length(X <- c(as.double(start), as.double(max(lower, -BIG)), 
		as.double(min(BIG, upper)), as.double(val), 
		as.double(con$tol), as.integer(length(start))))
	fn <- quote(objective(.par, ...))	
	rho = new.env(parent = environment())
    assign(".par",  X[1], envir = rho) 
    assign(".xlo",  X[2], envir = rho)
    assign(".xup",  X[3], envir = rho)
    assign(".val",  X[4], envir = rho)
    assign(".rpar", X[5], envir = rho)
    assign(".n",    X[6], envir = rho)
    
    # Call Solver:
    if (.trace) cat("\nIteration Path:\n\n")
    .Call("mfopt", fn, rho, PACKAGE = "FortranCallsR")
	
    # Extract Results:
    ans = list(
    	par = get(".par", envir = rho), 
    	objective = get(".val", envir = rho),
    	tol = con$tol)
    if (.trace) cat("\n\nFinal Result:\n\n")
	
	# Return Value:
	ans
}		


# ------------------------------------------------------------------------------

