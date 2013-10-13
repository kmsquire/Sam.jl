# Bam.jl

export BamReader, read_alignment

# bam files

type BamReader
    io::IO
    meta::SamMeta
end

const BAM_MAGIC="BAM\1"

function BamReader(io::IO)
    io = BGZF.BGZFReader(io)
    
    magic = read(io, Array(Uint8, 4))
    if magic != BAM_MAGIC.data
        error("Not a bam file")
    end
    
    header = read_bam_header(io)
    if first(header)[1] != "HD"
        error("Bam header does not start with HD line")
    end

    sort_order = header["HD"][1]["SO"]
    version = header["HD"][1]["VN"]

    refs = read_bam_refs(io)

    # Make sure refs match
    if haskey(header, "SQ")
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
    BamReader(io, SamMeta(version, sort_order, header, refs))
end

function read_raw_bam_header(io::IO)
    l_text = read(io, Int32)
    read(io, Array(Uint8, l_text))
end

function read_bam_header(io::IO)
    # Read header
    header = read_raw_bam_header(io) |> bytestring |> rstrip |> x->split(x,"\n") |> parse_headerlines        
end

function read_bam_refs(io::IO)
    n_ref = read(io, Int32)

    refs = RefSeq[]

    for i = 1:n_ref
        l_name = read(io, Int32)
        name = rstrip(bytestring(read(io, Array(Uint8,l_name))),'\0')
        l_ref = read(io, Int32)
        push!(refs, RefSeq(name, l_ref))
    end

    refs
end

function read_alignment(b::BamReader)
    blocksize = read(b.io, Uint32)
    info = read(b.io, Array(BamAlignmentInfo,1))[1]

    readname = rstrip(bytestring(read(b.io, Array(Uint8, info.l_readname))), '\0')
    cigar = read(b.io, Array(Uint32, info.n_cigar_op))
    seq = read(b.io, Array(Uint8, (info.rlen+1)>>1))
    qual = read(b.io, Array(Uint8, info.rlen))
               
    bytesread = sizeof(info) + 
                info.l_readname + 
                info.n_cigar_op<<2 + 
                (info.rlen+1)>>1 +
                info.rlen
    # TODO: parse
    aux = read(b.io, Array(Uint8, blocksize-bytesread))

    BamAlignment(info, readname, cigar, seq, qual, aux, b.meta)
end

close(br::BamReader) = close(br.io)
eof(br::BamReader) = eof(br.io)


# bam writer

type BamWriter
    io::IO
    meta::SamMeta
end

write_bam_magic(io::IO) = write(io, BAM_MAGIC)

function write_raw_bam_header(io::IO, header::Array{Uint8})
    write(io, int32(length(header)))
    write(io, header)
end

function write_bam_header(io::IO, header::SamHeaderData)
    header = [h*"\n" for h in compose_headerlines(header)] |> join
    write_raw_bam_header(io, header.data)
end

function write_bam_refs(io::IO, refs::Vector{RefSeq})
    write(io, int32(length(refs)))

    for ref in refs
        write(io, int32(length(ref.ref)))
        write(io, ref.ref, "\0")
        write(io, ref.size)
    end
end

function write_alignment(b::BamWriter, r::BamAlignment)
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

close(bw::BamWriter) = close(bw.io)
eof(bw::BamWriter) = eof(bw.io)

