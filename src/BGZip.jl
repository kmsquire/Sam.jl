module BGZip

using StrPack
using Zlib

import Base: read, readuntil, readline, readall, nb_available, eof, open, close

export read, readuntil, readline, readall, nb_available, eof, close


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
        extra_buffer = IOBuffer(read(io, Uint8[XLEN]))
        while !eof(extra_buffer)
            push(extra, read(extra_buffer, GZipExtra))
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
    buffer::IOBuffer

    readsize::Int
    rawbuf::Array{Uint8}
end

nb_available(io::GZStream) = nb_available(io.buffer)


function open(fname::String, mode::String="r", readsize = 65536)
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
    GZStream(header, io, PipeBuffer(), readsize, Uint8[])
end


function read(this::GZStream, ::Type{Uint8})
    buf = this.buffer
    wait_readnb(this, 1)
    read(buf, Uint8)
end


function read{T}(this::GZStream, a::Array{T})
    assert(isbits(T),"Read from Buffer only supports bits types or arrays of bits types")
    nb = length(a)*sizeof(T)
    buf = this.buffer
    @assert buf.seekable == false
    @assert buf.maxsize >= nb
    wait_readnb(this,nb)
    read(buf, a)
    return a
end

function readline(this::GZStream)
    wait_readline(this)
    readline(this.buffer)
end

function readuntil(this::GZStream,c::Uint8)
    buf = this.buffer
    @assert buf.seekable == false
    wait_readbyte(this,c)
    readuntil(buf,c)
end

function readall(this::GZStream)
    buf = this.buffer
    @assert buf.seekable == false
    wait_readall(this)
    takebuf_string(this.buffer)
end

eof(s::GZStream) = eof(s.buffer) && eof(s.io)

close(s::GZStream) = (close(s.io); close(s.buffer))

function read_decompress(x::GZStream, nb::Int)
    resize!(x.rawbuf, nb)
    read(x.io, x.rawbuf)

    sz = length(x.buffer.data)
    decompress(x.rawbuf, true, x.buffer.data)
    x.buffer.size += length(x.buffer.data) - sz
end


function wait_readnb(x::GZStream, nb::Int)
    while !eof(x.io) && (nba = nb_available(x.buffer)) < nb
        bytes_to_read = max((nb-nba), min(x.readsize, nb_available(x.io)))
        read_decompress(x, bytes_to_read)
    end
end

function wait_readbyte(x::GZStream, c::Uint8)
    while !eof(x.io) && search(x.buffer,c) <= 0
        bytes_to_read = min(x.readsize, nb_available(x.io))
        read_decompress(x, bytes_to_read)
    end
end

wait_readline(x) = wait_readbyte(x, uint8('\n'))

function wait_readall(x::GZStream)
    while !eof(x.io)
        bytes_to_read = min(x.readsize, nb_available(x.io))
        read_decompress(x, bytes_to_read)
    end
end

end
