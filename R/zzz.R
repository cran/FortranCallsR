
	
.First.lib =  
function(lib, pkg)
{   # A function implemented by Diethelm Wuertz
    
  
    # Load dll:
    library.dynam("FortranCallsR", pkg, lib)
}
