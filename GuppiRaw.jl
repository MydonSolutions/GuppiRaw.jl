module GuppiRaw
	export Guppi
	export parseValueString
	export readHeader
	export Int4
	export readComplex
	export readGuppiBlock
	export writeComplex!
	export writeHeader!

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
			push!(headers, readHeader!(fio))
			blockSize = headers[1]["BLOCSIZE"] # should read per header

			while !eof(fio) && length(headers) < 128
				push!(blockByteOffsets, position(fio))
				skip(fio, blockSize)
				push!(headers, readHeader!(fio))
			end
			println()
			return Guppi(filepath, headers, blockByteOffsets)
		end
	end

	function readGuppiBlock(gp, bnum::Integer)
		dimensionKeys = [
			"NPOL",
			"PKTNTIME",
			"PKTNCHAN",
			"NSTRM",
			"NANTS"
		]

		npkts = div(gp.headers[bnum]["PIPERBLK"], gp.headers[bnum]["PKTNTIME"])
		dimensions = [1:gp.headers[bnum][key] for key in dimensionKeys]
		push!(dimensions, 1:npkts)
		nbits = gp.headers[bnum]["NBITS"]
		open(gp.filepath, "r") do fio
			seek(fio, gp.blockByteOffsets[bnum])
			return [readComplex!(fio, nbits) for i in CartesianIndices(Tuple(dimensions))]
		end
	end

	function parseValueString(str)
		ret = tryparse(Int, str)
		ret = (ret == nothing ? tryparse(Float32, str) : ret)
		ret = (ret == nothing ? string(str) : ret)
		ret
	end

	function readHeader!(fio::IO)
		headerDict = GuppiHeaderDict()
		while true
			str = string(Char.(read(fio, 80))...)
			parts = split(str, "=")
			if length(parts) != 2 || contains(str, "END"*" "^77)
				break
			end
			headerDict[strip(parts[1])] = parseValueString(strip(parts[2]))
		end
		headerDict
	end

	function Int4(val::Integer)::Int8
		return (val&0xf) > 0x7 ? Int8(-8+(val&0x7)) : (val&0x7)
	end

	function readComplex!(fio::IO, nbits)
		if nbits == 4
			byte = read(fio, UInt8)
			return Int4(byte>>4) + Int4(byte)im
		elseif nbits == 8
			return Int8(read(fio, Int8)) + Int8(read(fio, Int8))im
		elseif nbits == 16
			return Int16(read(fio, Int16)) + Int16(read(fio, Int16))im
		else
			return nothing
		end
	end

	function writeComplex!(fio::IO, c::Complex, nbits)
		if nbits == 4
			byte = write(fio, UInt8((real(c)&0xf) << 4 + imag(c)&0xf))
		elseif nbits == 8
			return write(fio, Int8(real(c)&0xff), Int8(imag(c)&0xff))
		elseif nbits == 16
			return write(fio, Int16(real(c)&0xffff), Int16(imag(c)&0xffff))
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