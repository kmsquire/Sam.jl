# SamMeta.jl

using DataStructures
import Base: push!, readline, close, eof, isequal


# Cigar strings
cigar_ops = "MIDNSHP=X"
cigar_op_num(x::Char) = search(cigar_ops, x)
cigar_op_num(x::Uint8) = search(cigar_ops.data, x)
cigar_op(n::Uint32) = cigar_ops[n+1]

function bam_cigar2str(bc::Vector{Uint32})
    cigar = Uint8[]
    for c in bc
        op_len = c >> 4
        op = cigar_op(c & 0x00000F)
        append!(cigar, string(op_len).data)
        push!(cigar, uint8(op))
    end
    bytestring(cigar)
end

# Sequence

nucs = "=ACMGRSVTWYHKDBN"
nuc_num(x::Char) = search(nucs, x)


### Aliases for handling samtools metadata
typealias SamData Dict{String, String, Ordered}
typealias SamHeaderData Dict{String, Vector{SamData}, Ordered}

get_header_refs(header::SamHeaderData) = [RefSeq(s["SN"], int32(s["LN"])) for s in header["SQ"]]

# Reference Sequence

immutable RefSeq
    rname::String
    size::Int32
end

isequal(a::RefSeq, b::RefSeq) = (isequal(a.rname, b.rname) && a.size == b.size)
    
type SamMeta
    version::String
    sort_order::String
    header::SamHeaderData
    refs::Array{RefSeq}
end

function SamMeta(;sam_version=1.0, sort_order="unknown")
    HD = SamHeaderData(("HD", SamData(("VN", "SO"), (string(sam_version), sort_order))))
    SamMeta(string(sam_version), sort_order, HD, RefSeq[])
end

function push!(sm::SamHeaderData, kv::Tuple)
    (tag, value) = kv
    if has(sm, tag)
        push!(sm[tag], value)
    else
        sm[tag] = [value]
    end
    kv
end



    
