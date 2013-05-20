using StrPack

import Base.read

const FTEXT    = 0b00000001
const FHRCR    = 0b00000010
const FEXTRA   = 0b00000100
const FNAME    = 0b00001000
const FCOMMENT = 0b00010000

@struct type GZipHeader
    ID1   ::Uint8
    ID2   ::Uint8
    CM    ::Uint8
    FLG   ::Uint8
    MTIME ::Uint32
    XFL   ::Uint8
    OS    ::Uint8
end

const strpack_gzheader_asize = StrPack.STRUCT_REGISTRY[GZipHeader].asize
const gzheader_size = sizeof(GZipHeader)

@struct type GZipExtraHeader
    SI1   ::Uint8
    SI2   ::Uint8
    LEN   ::Uint16
end

type GZipExtra
    header::GZipExtraHeader
    data  ::Array{Uint8}
end

const strpack_gzextra_asize = StrPack.STRUCT_REGISTRY[GZipHeader].asize
const gzextra_size = sizeof(GZipExtra)

type GZipFile
    header ::GZipHeader
    xlen   ::Uint16
    extra  ::Array{GZipExtra}
    fname  ::String
    comment::String
    fhcrc  ::Uint16
    crc32  ::Uint32
    isize  ::Uint32
end

read(io::IO, GZipHeader) = unpack(io, GZipHeader, strpack_gzheader_asize, align_packed, :LittleEndian)
read(io::IO, GZipExtra) = unpack(io, GZipExtra, strpack_gzextra_asize, align_packed, :LittleEndian)

function read_header(io::IO)
    header = read(io, GZipHeader)
    if header.FLG & FEXTRA > 0
        local XLEN = read(io, Uint16)
        local xlen_buffer = IOBuffer(read(io, Uint8[XLEN]))
    end
    if header.FLG & FNAME > 0
        local filename = readuntil(io, '\0')
    end
    if header.FLG & FCOMMENT > 0
        local comment = readuntil(io, '\0')
    end
    if header.FLG & FHCRC > 0
        local crc16 = read(io, Uint16)
    end
end
