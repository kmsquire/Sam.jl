# Sam.jl

using OrderedCollections

# sam files

type SamFile
    io::IO
    meta::SamMeta
end

# sam file header

const SAM_MAGIC = "@HD\t"

function parse_headerline(line::String)
    linetag = line[2:3]
    tagdict = Dict{String,String,Ordered}()
    for kv in split(line[5:end], "\t")
        k = kv[1:2]
        v = kv[4:end]
        tagdict[k] = v
    end
    tagdict
    (linetag, tagdict)
end

function compose_headerline(linetag::String, tagdict::Dict{String,String})
    headerline = ["@"*linetag]
    for (k,v) in d
        push!(headerline, k*":"*v)
    end
    join(headerline, "\t")
end

function parse_headerlines{S<:String}(headerlines::Vector{S})
    header = SamHeaderData()
    for line in headerlines
        (k,v) = parse_headerline(line)
        if k in keys(header)
            push!(header[k], v)
        else
            header[k] = [v]
        end
        #push!(header, parse_headerline(line))
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

function read_sam_header(io::IO)
    # Check magic string
    ## currently checked before calling 
    #read(io, Array(Uint8,4)) == SAM_MAGIC.data || error("Not a sam file.")
    
    headerlines = [SAM_MAGIC * rstrip(readline(io))]
    
    c = Base.peek(io)
    while c == '@'
        push!(headerlines, rstrip(readline(io)))
        try
            c = Base.peek(io)
        catch e
            # Occurs, e.g., when there are no reads in a file
            break
        end
    end

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

read_alignment(s::SamFile) = SamAlignment(split(readline(s.io), '\t', 12)...)


function samopen(io::IO, mode="r"; meta=nothing)
    if mode == "r"
        magic = read(io, Array(Uint8, 4))
        if magic == BAM_MAGIC.data
            BamFile(io, header=read_raw_bam_header(io))
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

close(sf::SamFile) = close(sf.io)
eof(sf::SamFile) = eof(sf.io)

    
