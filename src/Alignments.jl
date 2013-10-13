
import Base: convert, show

export convert, show

# sam file alignments

type SamAlignment
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


##################

# bam file alignments

immutable BamAlignmentInfo
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

type BamAlignment
    info::BamAlignmentInfo
    readname::String
    cigar::Vector{Uint32}
    seq::Vector{Uint8}
    qual::Vector{Uint8}
    aux::Vector{Uint8}
    meta::SamMeta
end

function unpack_seq(cseq::Vector{Uint8}, rlen::Int32)
    seq = Array(Uint8, 0)
    sizehint(seq, rlen)
    count = 0
    for b in cseq
        push!(seq, nucs.data[b>>4+1])
        count += 1
        if count == rlen;  break;  end
        push!(seq, nucs.data[b&0x0F+1])
        count += 1
    end
    bytestring(seq)
end

function unpack_aux(seq::Vector{Uint8})
    return ""
    #string(seq)
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
                 unpack_seq(ba.seq, ba.info.rlen),
                 bytestring(ba.qual+33),
                 unpack_aux(ba.aux))

show(io::IO, ba::BamAlignment) = show(io, convert(SamAlignment, ba, ba.meta))
