

immutable Chunk
    chunk_beg::Uint64
    chunk_end::Uint64
end

immutable Bin
    bin::Uint32
    n_chunk::Uint32
    chunks::Vector{Chunk}
end

immutable Interval_16
    ioffset::Uint64
end

immutable RefIndex
    n_bin::Int32
    bins::Vector{Bin}
    n_intv::Int32
    intv::Vector{Interval_16}
end

immutable BamIndex
    n_ref::Int32
    refs::Vector{RefIndex}
end

const BAI_MAGIC="BAI\1"

function BamIndex(fn::String)
    if !isfile(fn)
        error(fn, " not found!")
    end
        
    open(fn) do io
        BamIndex(io)
    end
end

function BamIndex(io::IO)
    if read(io, Array(Uint,4)) != BAI_MAGIC.data
        error("Not a bam index")
    end

    n_ref = read(io, Int32)
    refbins = Array(RefIndex, n_ref)
    for ref = 1:nref
        n_bin = read(io, Int32)
        bins = Array(Bin, n_bin)
        for bin = 1:n_bin
            binnum = read(io, Uint32)
            n_chunk = read(io, Int32)
            chunks = read(io, Array(Chunk, n_chunk))

            bins[bin] = Bin(binnum, n_chunk, chunks)
        end
        
        n_intv = read(io, Int32)
        intv = read(io, Array(Interval_16, n_intv))

        refbins[ref] = RefIndex(n_ref, n_bin, bins, n_intv, intv)
    end
    BamIndex(refbins)
end

function region2bin(idx::BamIndex, chr::Int, r_beg::Int, r_end::Int)
    r_beg -= 1
    r_end -= 1
    if (r_beg >> 14 == r_end >> 14) return div(((1 << 15)-1),7) + (r_beg >> 14) + 1 end
    if (r_beg >> 17 == r_end >> 17) return div(((1 << 12)-1),7) + (r_beg >> 17) + 1 end
    if (r_beg >> 20 == r_end >> 20) return div(((1 <<  9)-1),7) + (r_beg >> 20) + 1 end
    if (r_beg >> 23 == r_end >> 23) return div(((1 <<  6)-1),7) + (r_beg >> 23) + 1 end
    if (r_beg >> 26 == r_end >> 26) return div(((1 <<  3)-1),7) + (r_beg >> 26) + 1 end
    return 1
end

const MAX_BIN=div(((1<<18)-1),7)

function region2bins(r_beg::Int, r_end::Int)
    bins = Array(Bins,0)
    size_hint(bin_list, MAX_BIN)

    r_beg -= 1
    r_end -= 1

    push!(bins, 1)
    for k =    1 + (r_beg>>26) :    1 + (r_end>>26);  push!(bins, k+1);  end
    for k =    9 + (r_beg>>23) :    9 + (r_end>>23);  push!(bins, k+1);  end
    for k =   73 + (r_beg>>20) :   73 + (r_end>>20);  push!(bins, k+1);  end
    for k =  585 + (r_beg>>17) :  585 + (r_end>>17);  push!(bins, k+1);  end
    for k = 4681 + (r_beg>>14) : 4681 + (r_end>>14);  push!(bins, k+1);  end

    bins
end
