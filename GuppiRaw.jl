module GuppiRaw
	export Guppi
	export parseValueString
	export readHeader
	export Int4
	export readComplex
	export readGuppiBlock
	export writeComplex!
	export writeHeader!
	export calcBlockSize
	export calcBlockShape

	GuppiHeaderDict = Dict{String, Union{Int, Float32, String}}

	struct Guppi
		filepath::String
		headers::Array{GuppiHeaderDict}
		blockByteOffsets::Array{Unsigned}
	end
	
	function Guppi(filepath)
		headers::Array{GuppiHeaderDict} = []
		blockByteOffsets::Array{Unsigned} = []

		open(filepath, "r") do fio
			while !eof(fio) && length(headers) < 128
				push!(headers, readHeader!(fio))
				blockSize = headers[1]["BLOCSIZE"]
				align512 = (haskey(headers[1], "DIRECTIO") ? headers[1]["DIRECTIO"] == 1 : false)
				dataSize = GuppiRaw.calcBlockSize(headers[end])
				@assert blockSize >= dataSize "$blockSize >= $dataSize"

				if align512; seekAligned!(fio, 512) end
				push!(blockByteOffsets, position(fio))
				skip(fio, blockSize)
				if align512; seekAligned!(fio, 512) end
			end
			println()
			return Guppi(filepath, headers, blockByteOffsets)
		end
	end

	function seekAligned!(fio, alignBoundary)
		pos = position(fio)
		offset = (pos+(alignBoundary-1)) & ~(alignBoundary-1) - pos
		skip(fio, offset)
	end
	
	function readHeader!(fio::IO)
		headerDict = GuppiHeaderDict()
		while true
			str = string(Char.(read(fio, 80))...)
			parts = split(str, "=")
			if contains(str, "END"*" "^77)
				break
			end
			if length(parts) != 2
				return false
			end
			
			headerDict[strip(parts[1])] = parseValueString(strip(parts[2]))
		end
		headerDict
	end

	function parseValueString(str)
		ret = tryparse(Int, str)
		ret = (isnothing(ret) ? tryparse(Float32, str) : ret)
		ret = (isnothing(ret) ? string(str) : ret)
		ret
	end
	
	function readGuppiBlock(gp, bnum::Integer)
		blockShape = calcBlockShape(gp.headers[bnum])
		nbits = gp.headers[bnum]["NBITS"]
		open(gp.filepath, "r") do fio
			seek(fio, gp.blockByteOffsets[bnum])
			return [readComplex!(fio, nbits) for i in CartesianIndices(blockShape)]
		end
	end

	function calcBlockSize(header::GuppiHeaderDict)
		prod(length.(calcBlockShape(header)))*header["NBITS"]*2/8 # assumes each sample is re+im
	end 

	function calcBlockShape(header::GuppiHeaderDict)
		dimensionKeys = [
			"NPOL",
			"PKTNTIME",
			"PKTNCHAN",
			"PIPERBLK",
			"NSTRM",
			"NANTS"
		]

		if all(haskey(header, key) for key in dimensionKeys) == false
			println.(i for i in header)
		end

		dimensions = [1:header[key] for key in dimensionKeys]
		dimensions[4] = 1:div(length(dimensions[4]), header["PKTNTIME"])
		Tuple(dimensions)
	end

	function Int4(val::Integer)::Int8
		return (val&0xf) > 0x7 ? Int8(-8+(val&0x7)) : (val&0x7)
	end

	function UInt4(val::Integer)::Int8
		return val&0xf
	end

	function readComplex!(fio::IO, nbits; unsigned = true)
		if nbits == 4
			byte = read(fio, UInt8)
			return unsigned ? 
					UInt4(byte>>4) + UInt4(byte)im :
					Int4(byte>>4) + Int4(byte)im
		elseif nbits == 8
			return unsigned ? 
					UInt8(read(fio, Int8)) + UInt8(read(fio, Int8))im :
					Int8(read(fio, Int8)) + Int8(read(fio, Int8))im
		elseif nbits == 16
			return unsigned ? 
					UInt16(read(fio, Int16)) + UInt16(read(fio, Int16))im :
					Int16(read(fio, Int16)) + Int16(read(fio, Int16))im
		else
			return nothing
		end
	end

	function writeComplex!(fio::IO, c::Complex, nbits; cast4bit=false)
		if nbits == 4
			byte = write(fio, UInt8((real(c)&0xf) << 4 + (imag(c)&0xf)))
		elseif nbits == 8
			mask = cast4bit ? 0x0f : 0xff
			return write(fio, UInt8(real(c)&mask), UInt8(imag(c)&mask))
		elseif nbits == 16
			return write(fio, UInt16(real(c)&0xffff), UInt16(imag(c)&0xffff))
		else
			return nothing
		end
	end

	function writeHeader!(fio::IO, header::GuppiHeaderDict)
		for kv in header
			write(fio, rpad(kv[1]*"="*string(kv[2]), 80))
		end
		write(fio, rpad("END", 80))
	end

end