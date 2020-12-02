using Printf
include("./GuppiRaw.jl/GuppiRaw.jl")

gp4 = GuppiRaw.Guppi("/mnt/buf0/dmpauto_wael/GUPPI/guppi_59179_74555_000000_Unknown_0001.0000.raw")
gp4 = GuppiRaw.Guppi("/mnt/buf0/dmpauto_wael/GUPPI/guppi_59184_32766_153781_Unknown_0001.0000.raw")

idx = [gp4.headers[i]["PKTIDX"] for i in 1:10]


[header["NDROP"] for header in gp4.headers]
gp4.headers[3]["PIPERBLK"]/gp4.headers[3]["PKTNTIME"]


fp8 = "./upsamples/guppi_59179_74555_000000_Unknown_0001.0000.raw"

fio8 = open(fp8,"w")
for bIdx in 1:length(gp4.headers)
	header = copy(gp4.headers[bIdx])
	header["NBITS"] = 8
	header["BLOCSIZE"] *= 2
	header["DIRECTIO"] = 0
	GuppiRaw.writeHeader!(fio8, header)
	for cmplx in GuppiRaw.readGuppiBlock(gp4, bIdx, unsigned=false)
		GuppiRaw.writeComplex!(fio8, cmplx, 8)
	end
end
close(fio8)

blk4 = GuppiRaw.readGuppiBlock(gp4, 1, unsigned=false)

gp8 = GuppiRaw.Guppi(fp8)
blk8 = GuppiRaw.readGuppiBlock(gp8, 1, unsigned=false)

blk4 == blk8




blk4bit = GuppiRaw.readGuppiBlock(gp4, 1, unsigned=false)

blk4bit[1, 1,  1:32, 1, 1]

chan1 = blk4bit[1, :,  1, 1, 1]

unique(chan1)
findfirst(x->x== 0+2im, chan1)
length(chan1)/4

gp4 = GuppiRaw.Guppi("/mnt/buf0/timecontig/GUPPI/./guppi_40587_00013_000002_Unknown_0001.0000.raw")
gp4 = GuppiRaw.Guppi("/mnt/buf0/dmpauto_wael/GUPPI/./guppi_59172_55282_001055_Unknown_0001.0007.raw")

blk4bit = GuppiRaw.readGuppiBlock(gp4, 1, unsigned=false)


### Outdated automated verification

average([header["NDROP"]/pktperblk for header in gp.headers])

blk = GuppiRaw.readGuppiBlock(gp, 1)

blk[:, :, :, 1, 1, 1]

average(v) = sum(v)/length(v)

for blockI in 1:length(gp.blockByteOffsets)
	blk = GuppiRaw.readGuppiBlock(gp, blockI, unsigned=false)
	println(size(blk))
	breakBool = false
	for fengI in 1:size(blk, 6)
		for strmI in 1:size(blk, 5)

			for pktI in 1:size(blk, 4)
				print((pktI, strmI, fengI, blockI), "\t\t\t\t\r")
				for chanN in 1:size(blk, 3)
					c = GuppiRaw.Int4(fengI-1) + GuppiRaw.Int4(strmI-1 + (chanN%2 == 1 ? 0 : 8))*1im

					a = average(Complex{Int}.(blk[:,:,chanN, pktI, strmI, fengI]))
					if c != a
					# if any(blk[:,:,chanN, pktI, strmI, fengI] .!= c)
						println("\nChan "*string(chanN)*"\t"*string(a)*"!="*string(c))
						breakBool = true
						break
					end
				end
				if breakBool; break end
			end
			# if breakBool; break end
		end
		# if breakBool; break end
	end
	# if breakBool; break end
end
