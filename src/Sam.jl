# Sam.jl

module Sam

import BGZF

include("SamMeta.jl")
include("Alignments.jl")
include("SamFile.jl")
include("BamFile.jl")
#include("Bai.jl")

end # module Sam
    
