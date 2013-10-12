module BGZip

using StrPack
import Zlib

import Base: read, eof, close #, nb_available

export read, nb_available, eof, close


const FTEXT    = 0b00000001
const FHCRC    = 0b00000010
const FEXTRA   = 0b00000100
const FNAME    = 0b00001000
const FCOMMENT = 0b00010000


@struct type GZipHeaderFixed
    ID1   ::Uint8
    ID2   ::Uint8
    CM    ::Uint8
    FLG   ::Uint8
    MTIME ::Uint32
    XFL   ::Uint8
    OS    ::Uint8
end

const strpack_gzheader_asize = StrPack.STRUCT_REGISTRY[GZipHeaderFixed].asize
const gzheader_size = sizeof(GZipHeaderFixed)

read(io::IO, ::Type{GZipHeaderFixed}) = unpack(io, GZipHeaderFixed, strpack_gzheader_asize, align_packed, :LittleEndian)

###

@struct type GZipExtraHeader
    SI1   ::Uint8
    SI2   ::Uint8
    LEN   ::Uint16
end

const strpack_gzextra_asize = StrPack.STRUCT_REGISTRY[GZipExtraHeader].asize
const gzextra_size = sizeof(GZipExtraHeader)

read(io::IO, ::Type{GZipExtraHeader}) = unpack(io, GZipExtraHeader, strpack_gzextra_asize, align_packed, :LittleEndian)

###

type GZipExtra
    header::GZipExtraHeader
    data  ::Array{Uint8}
end

function read(io::IO, ::Type{GZipExtra})
    header = read(io, GZipExtraHeader)
    data = read(io, Array(Uint8, header.LEN))
    GZipExtra(header, data)
end

###

type GZipHeader
    header ::GZipHeaderFixed
    xlen   ::Uint16
    extra  ::Array{GZipExtra}
    fname  ::String
    comment::String
    fhcrc  ::Uint16
end

function read(io::IO, ::Type{GZipHeader})
    header = read(io, GZipHeaderFixed)
    xlen   = uint16(0)
    extra  = GZipExtra[]

    if header.FLG & FEXTRA > 0
        xlen = read(io, Uint16)
        extra_buffer = IOBuffer(read(io, Array(Uint8,xlen)))
        while !eof(extra_buffer)
            push!(extra, read(extra_buffer, GZipExtra))
        end
    end

    fname   = header.FLG & FNAME    > 0 ? readuntil(io, '\0') : ""
    comment = header.FLG & FCOMMENT > 0 ? readuntil(io, '\0') : ""
    fhcrc   = header.FLG & FHCRC    > 0 ? read(io, Uint16)    : uint16(0)

    GZipHeader(header, xlen, extra, fname, comment, fhcrc)
end

###

type GZStream <: IO
    header::GZipHeader
    io::IO
end

#nb_available(io::GZStream) = nb_available(io.io)


function open(fname::String, mode::String="r")
    io = Base.open(fname, mode)

    # Only handle read-only streams for now
    if mode != "r";  return io;  end

    ## Check for the GZip ID

    # mark(io)
    ID = read(io, Array(Uint8, 2))
    # reset(io)
    seek(io, 0)

    if ID != [31, 139];  return io;  end

    ## Read in the header and return a GZStream

    header = read(io, GZipHeader)
    GZStream(header, Zlib.Reader(io))
end

read(this::GZStream, ::Type{Uint8}) = Base.read(this.io, x)
read(this::GZStream, x::Array) = Base.read(this.io, x)
read(this::GZStream, x::BitArray) = Base.read(this.io, x)
read(this::GZStream, x::AbstractArray) = Base.read(this.io, x)

eof(s::GZStream) = eof(s.io)

close(s::GZStream) = close(s.io)

end
