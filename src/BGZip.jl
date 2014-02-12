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


###

type GZipExtra
    SI1   ::Uint8
    SI2   ::Uint8
    LEN   ::Uint16
    data  ::Array{Uint8}
end

function read(io::IO, ::Type{GZipExtra})
    si1 = read(io, Uint8)
    si2 = read(io, Uint8)
    len = read(io, Uint16)
    data = read(io, Array(Uint8, len))
    GZipExtra(si1, si2, len, data)
end

###

immutable GZipHeader
    ID1   ::Uint8
    ID2   ::Uint8
    CM    ::Uint8
    FLG   ::Uint8
    MTIME ::Uint32
    XFL   ::Uint8
    OS    ::Uint8
    xlen   ::Uint16
    extra  ::Array{GZipExtra}
    fname  ::String
    comment::String
    fhcrc  ::Uint16
end

function read(io::IO, ::Type{GZipHeader})
    id1   = read(io, Uint8)
    id2   = read(io, Uint8)
    cm    = read(io, Uint8)
    flg   = read(io, Uint8)
    mtime = read(io, Uint32)
    xfl   = read(io, Uint8)
    os    = read(io, Uint8)

    xlen  = uint16(0)
    extra = GZipExtra[]

    if flg & FEXTRA > 0
        xlen = read(io, Uint16)
        extra_buffer = IOBuffer(read(io, Array(Uint8,xlen)))
        while !eof(extra_buffer)
            push!(extra, read(extra_buffer, GZipExtra))
        end
    end

    fname   = flg & FNAME    > 0 ? readuntil(io, '\0') : ""
    comment = flg & FCOMMENT > 0 ? readuntil(io, '\0') : ""
    fhcrc   = flg & FHCRC    > 0 ? read(io, Uint16)    : uint16(0)

    GZipHeader(id1, id2, cm, flg, mtime, xfl, os, xlen, extra, fname, comment, fhcrc)
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
    GZStream(header, Zlib.Reader(io, true))
end

read(this::GZStream, args...) = read(this.io, args...)

eof(s::GZStream) = eof(s.io)

close(s::GZStream) = close(s.io)

end
