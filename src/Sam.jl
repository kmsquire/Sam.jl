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
typealias SamData OrderedDict{String, String}
typealias SamHeaderData OrderedDict{String, Vector{SamData}}


# Reference Sequence

type RefSeq
    rname::String
    size::Int32
end

type SamMeta
    header::SamHeaderData
    refs::Array{RefSeq}
end

function SamMeta(;sam_version=1.0, sort_order="unknown")
    HD = SamHeaderData(("HD", SamData(("VN", "SO"), (string(sam_version), sort_order))))
    SamMeta(HD, RefSeq[])
end

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
    meta::SamMeta
end

# sam file header

const SAM_MAGIC = "@HD\t"

function parse_headerline(line::String)
    linetag = line[2:3]
    tagdict = OrderedDict{String,String}()
    for kv in split(str, "\t")
        k = kv[1:2]
        v = kv[4:end]
        tagdict[k] = v
    end
    tagdict
    (linetag, tagdict)
end

function compose_headerline(linetag::String, tagdict::OrderedDict{String,String})
    headerline = ["@"*linetag]
    for (k,v) in d
        push!(headerline, k*":"*v)
    end
    join(headerline, "\t")
end

function parse_headerlines{S<:String}(headerlines::Vector{S})
    header = SamHeaderData()
    for line in headerlines
        push!(header, parse_header_line(line))
    end
    header
end

function compose_headerlines(header::SamHeaderData)
    headerlines = String[]
    for (linetag, tds) in header
        for td in tds
            push!(headerlines, compose_headerline(linetag, td))
        end
    end
    headerlines
end

function read_sam_header(io::GZipStream)
    # Check magic string
    ## currently checked before calling 
    #read(io, Array(Uint8,4)) == SAM_MAGIC.data || error("Not a sam file.")
    
    headerlines = [SAM_MAGIC * rstrip(readline(io))]
    
    c = GZip.gzgetc(io)
    while c == '@'
        GZip.gzungetc(c,io)
        push!(headerlines, rstrip(readline(io)))
        try
            c = GZip.gzgetc(io)
        catch e
            # Occurs, e.g., when there are no reads in a file
            return parse_headerlines(headerlines)
        end
    end
    GZip.gzungetc(c,io)

    parse_headerlines(headerlines)
end

function write_sam_header(io::IO, header::SamHeaderData)
    if first(header)[1] != "HD"
        error("write_sam_header: header must start with HD line")
    end
    for line in compose_headerlines(header)
        write(io, line, "\n")
    end
end

# sam file alignments

abstract Alignment

type SamAlignment <: Alignment
    qname::String
    flag ::Uint16
    rname::String
    pos  ::Int32
    mapq ::Uint8
    cigar::String
    rnext::String
    pnext::Int32
    tlen ::Int32
    seq  ::String
    qual ::String
    aux  ::String
end

SamAlignment(qname::String,
             flag::String,
             rname::String,
             pos::String,
             mapq::String,
             cigar::String,
             rnext::String,
             pnext::String,
             tlen::String,
             seq::String,
             qual::String,
             aux::String) =
    SamAlignment(       qname, 
                 uint16(flag), 
                        rname,
                 int32 (pos),
                 uint8 (mapq),
                        cigar,
                        rnext,
                 int32 (pnext),
                 int32 (tlen),
                        seq,
                        qual,
                        aux)

read_alignment(s::SamFile) = SamAlignment(split(readline(s.io), '\t', 12)...)
readline(sf::AbstractSamFile) = read_alignment(sf)

# bam files

type BamFile <: AbstractSamFile
    io::IO
    meta::SamMeta
end

function BamFile(io::IO; header=nothing, refs=nothing)
    local sort_order
    if header == nothing
        header = SamHeaderData()
        sort_order = "unknown"
    else
        sort_order = header["HD"][1]["SO"]
    end
    if refs == nothing
        refs = RefSeq[]
    end
    BamFile(io, header, refs, sort_order)
end

# bam file header

BAM_MAGIC="BAM\1"

write_bam_magic(io::IO) = write(io, BAM_MAGIC)

function read_raw_bam_header(io::IO)
    l_text = read(io, Int32)
    read(io, Array(Uint8, l_text))
end

function write_raw_bam_header(io::IO, header::Array{Uint8})
    write(io, int32(length(header)))
    write(io, header)
end

function read_bam_header(io::IO)
    # Check magic string
    ## currently checked before calling 
    #read(io, Array(Uint8,4)) == BAM_magic || error("Not a bam file.")
    
    # Read header
    header = read_raw_bam_header(io) | bytestring | rstrip | x->split(x,"\n") | parse_headerlines        
end

function write_bam_header(io::IO, header::SamHeaderData)
    header = [h*"\n" for h in compose_headerlines(header)] | join
    write_raw_bam_header(io, header.data)
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
    rlen::Int32           # was: l_seq::Int32
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

convert(::Type{SamAlignment}, ba::BamAlignment, meta::SamMeta) = 
    SamAlignment(ba.readname,
                 ba.info.flag,
                 ba.info.refID < 0 ? "*" : meta.refs[ba.info.refID+1].rname,
                 ba.info.pos,
                 ba.info.mapq,
                 bam_cigar2str(ba.cigar),
                 ba.info.next_refID < 0 ? "*" : meta.refs[ba.info.next_refID+1].rname,
                 ba.info.next_pos,
                 ba.info.tlen,
                 unpack_seq(ba.seq),
                 bytestring(ba.qual),
                 unpack_aux(ba.aux))

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

function write_bam_refs(io::IO, refs::Vector{RefSeq})
    write(io, int32(length(refs)))

    for ref in refs
        write(io, int32(length(ref.ref)))
        write(io, ref.ref, "\0")
        write(io, ref.size)
    end
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

function write_alignment(b::BamFile, r::BamAlignment)
    blocksize::Int32 = sizeof(AlignmentInfo) + 
                       length(readname) + 1 +
                       sizeof(cigar) + 
                       sizeof(seq) + 
                       sizeof(qual) + 
                       sizeof(aux)
    write(b.io, blocksize)
    pack(b.io, r.info, strpack_asize, align_packed, :LittleEndian)
    write(b.io, r.readname, "\0")
    write(b.io, r.cigar)
    write(b.io, r.seq)
    write(b.io, r.qual)
    write(b.io, r.aux)
end

get_header_refs(header::SamHeaderData) = [RefSeq(s["SN"], s["LN"]) for s in header["SQ"]]

function read_meta(io::IO)
    header = read_bam_header(io)
    if first(header)[1] != "HD"
        error("read_meta: header does not start with HD line")
    end

    refs = read_bam_refs(io)

    # Make sure refs match
    if has(header, "SQ")
        header_refs = get_header_refs(header)
        if length(header_refs) != length(refs)
            error("Different number of reference sequences in header and ref list... what happened?!")
        elseif any(header_refs .!= refs)
            pos = find(header_refs .!= refs)
            error("""Reference sequence mismatch between header and ref list at position(s) $pos:
                  header: 
                     $([pos header_refs[pos]])
                  ref list: 
                     $([pos refs[pos]])
                  """ )
        end
    end
    return SamMeta(header, refs)
end

function read_meta(io::IO)
    header = read_sam_header(io)
    if first(header)[1] != "HD"
        error("read_meta: header does not start with HD line")
    end

    refs = get_header_refs(header)
    return SamMeta(header, refs)
end


function samopen(io::GZipStream, mode="r"; meta=nothing)
    if mode == "r"
        magic = read(io, Array(Uint8, 4))
        if magic == BAM_MAGIC.data
            BamFile(io, read_meta(io))
        elseif magic == SAM_MAGIC.data
            SamFile(io, read_meta(io))
        else
            error("File does not seem to be a sam or bam file")
        end
    elseif mode == "w"
        if meta == nothing
            error("")
        end
        if splitext(io.name)[2] == ".bam"
            bam = BamFile(io, meta)

        else
            #Initialize writeable same
        end
    else
        error("mode must be \"w\" or \"r\"")
    end
end

samopen(fn::String, mode="r"; header=nothing, refs=nothing) = samopen(GZip.open(fn, mode), mode, header=header, refs=refs)
samopen(io::IOStream, mode="r"; header=nothing, refs=nothing) = samopen(GZip.gzdopen(fd(io)), mode, header=header, refs=refs)
const open = samopen

function samopen(f::Function, args...)
    io = samopen(args...)
    try
        f(io)
    finally
        close(io)
    end
end

#function getheader(fn::String)
#    io = samopen(

close(sf::AbstractSamFile) = close(sf.io)
eof(sf::AbstractSamFile) = eof(sf.io)

end # module Sam
    
