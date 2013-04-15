# Sam.jl

module Sam

using StrPack
using GZip

import OrderedCollections.OrderedDict
import Base.push!, Base.readline, Base.close, Base.eof

export SamFile, 
       BamFile, 
       Alignment, 
       BamAlignment, 
       SamAlignment, 
       readline, 
       close,
       eof

### Aliases for handling samtools metadata
typealias SamData OrderedDict{String, String}
typealias SamHeaderData OrderedDict{String, Vector{SamData}}

# function push!(sm::SamData, kv::(String, Any))
#     (tag, value) = kv
#     sm[tag] = value
#     kv
# end

function push!(sm::SamHeaderData, kv::(String, Any))
    (tag, value) = kv
    if has(sm, tag)
        push!(sm[tag], value)
    else
        sm[tag] = [value]
    end
    kv
end

# sam files
abstract AbstractSamFile

type SamFile <: AbstractSamFile
    io::IO
    header::SamHeaderData
    sort_order::String
end
SamFile(io::IO, header::SamHeaderData) = SamFile(io, header, header["HD"][1]["SO"])

# sam file header

const SAM_MAGIC = "@HD\t"

function parse_header_kvs(str::String)
    d = OrderedDict{String,String}()
    for kv in split(str, "\t")
        k = kv[1:2]
        v = kv[4:end]
        d[k] = v
    end
    d
end

function parse_header_line(line::String)
    linetag = line[2:3]
    tagdict = parse_header_kvs(line[5:end])
    (linetag, tagdict)
end

function parse_header_lines{S<:String}(header_lines::Vector{S})
    header = SamHeaderData()
    for line in header_lines
        push!(header, parse_header_line(line))
    end
    header
end

function read_sam_header(io::GZipStream)
    # Check magic string
    ## currently checked before calling 
    #read(io, Array(Uint8,4)) == SAM_MAGIC.data || error("Not a sam file.")
    
    header_lines = [SAM_MAGIC * rstrip(readline(io))]
    
    c = GZip.gzgetc(io)
    while c == '@'
        GZip.gzungetc(c,io)
        push!(header_lines, rstrip(readline(io)))
        try
            c = GZip.gzgetc(io)
        catch e
            # Occurs, e.g., when there are no reads in a file
            return parse_header_lines(header_lines)
        end
    end
    GZip.gzungetc(c,io)

    parse_header_lines(header_lines)
end

# sam file alignments

abstract Alignment

type SamAlignment <: Alignment
    qname::String
    flag::Uint16
    rname::String
    pos::Int32
    mapq::Uint8
    cigar::String
    rnext::String
    pnext::Int32
    tlen::Int32
    seq::String
    qual::String
    aux::String
end

SamAlignment(qname::String, flag::String, rname::String, pos::String, mapq::String, cigar::String, rnext::String, pnext::String, tlen::String, seq::String, qual::String, aux::String) = 
    SamAlignment(qname, uint16(flag), rname, int32(int32(pos)), uint8(mapq), cigar, rnext, int32(int32(pnext)), int32(int32(tlen)), seq, qual, aux)

read_alignment(s::SamFile) = SamAlignment(split(readline(s.io), '\t', 12)...)
readline(sf::AbstractSamFile) = read_alignment(sf)


# bam files

type RefSeq
    ref::String
    size::Int32
end

type BamFile <: AbstractSamFile
    io::IO
    header::SamHeaderData
    refs::Array{RefSeq}
    sort_order::String
end
BamFile(io::IO, header::SamHeaderData, refs::Array{RefSeq}) = BamFile(io, header, refs, header["HD"][1]["SO"])

# bam file header

BAM_MAGIC="BAM\1"

function read_bam_header(io::IO)
    # Check magic string
    ## currently checked before calling 
    #read(io, Array(Uint8,4)) == BAM_magic || error("Not a bam file.")
    
    # Read header
    l_text = read(io, Int32)
    header = bytestring(read(io, Array(Uint8, l_text))) | rstrip | x->split(x,"\n") | parse_header_lines        
end

# bam file alignments

@struct type AlignmentInfo
    refID::Int32
    pos::Int32
    l_readname::Uint8
    mapq::Uint8
    bin::Uint16
    n_cigar_op::Uint16
    flag::Uint16
    #l_seq::Int32
    rlen::Int32
    next_refID::Int32
    next_pos::Int32
    tlen::Int32
end

const strpack_asize=StrPack.STRUCT_REGISTRY[AlignmentInfo].asize
const ainfo_size = sizeof(AlignmentInfo)

type BamAlignment <: Alignment
    info::AlignmentInfo
    readname::String
    cigar::Vector{Uint32}
    seq::Vector{Uint8}
    qual::Vector{Uint8}
    aux::Vector{Uint8}
end

function read_bam_refs(io::IO)
    n_ref = read(io, Int32)

    refs = RefSeq[]

    for i = 1:n_ref
        l_name = read(io, Int32)
        name = rstrip(bytestring(read(io, Array(Uint8,l_name))),"\0")
        l_ref = read(io, Int32)
        push!(refs, RefSeq(name, l_ref))
    end

    refs
end

function read_alignment(b::BamFile)
    blocksize = read(b.io, Uint32)
    info = unpack(b.io, AlignmentInfo, strpack_asize, align_packed, :LittleEndian)

    readname = rstrip(bytestring(read(b.io, Array(Uint8, info.l_readname))), "\0")
    cigar = read(b.io, Array(Uint32, info.n_cigar_op))
    seq = read(b.io, Array(Uint8, (info.rlen+1)>>1))
    qual = read(b.io, Array(Uint8, info.rlen))
               
    bytesread = ainfo_size + 
                info.l_readname + 
                info.n_cigar_op<<2 + 
                (info.rlen+1)>>1 +
                info.rlen
    # TODO: parse
    aux = read(b.io, Array(Uint8, blocksize-bytesread))

    BamAlignment(info, readname, cigar, seq, qual, aux)
end


function samopen(io::GZipStream)
    magic = read(io, Array(Uint8, 4))
    if magic == BAM_MAGIC.data
        BamFile(io, read_bam_header(io), read_bam_refs(io))
    elseif magic == SAM_MAGIC.data
        SamFile(io, read_sam_header(io))
    else
        error("File does not seem to be a sam or bam file")
    end
end

samopen(fn::String) = samopen(GZip.open(fn))
samopen(io::IOStream) = samopen(GZip.gzdopen(fd(io)))
open = samopen

close(sf::AbstractSamFile) = close(sf.io)
eof(sf::AbstractSamFile) = eof(sf.io)

end # module Sam
    
