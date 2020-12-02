include("./GuppiRaw.jl/GuppiRaw.jl")
# include("./GuppiRaw.jl")

fp4 = "/homelocal/sonata/mydonsol/rawspec/inputs_4bit/guppiSigFile4bit-random.0000.raw"
fp8 = "/homelocal/sonata/mydonsol/rawspec/inputs/guppiSigFile8bit-4ant-8chan.0000.raw"

header4 = GuppiRaw.GuppiHeaderDict(
	"BLOCSIZE"=> 2^27,
	"PIPERBLK"=> 2^18,
	"PKTNTIME"=> 16,
	"NPOL"		=> 2,
	"NBITS"		=> 4,
	"NSTRM"		=> 1,
	"NANTS"		=> 1,
	"PKTNCHAN"=> 256,
	"OBSNCHAN"=> 256,
	"OBSFREQ" => 1420,
	"OBSBW" 	=> 128,
	"TBIN" 		=> 2.0f-6,
	"PKTIDX" 	=> 0
)
header8 = GuppiRaw.GuppiHeaderDict(
	# "BLOCSIZE"=> 2^28,
	# "PIPERBLK"=> div(2^28, 2*2*1*1024),
	"BLOCSIZE"=> 2^27,
	"PIPERBLK"=> 2^20,
	"PKTNTIME"=> 16,
	"NPOL"		=> 2,
	"NBITS"		=> 8,
	"NSTRM"		=> 1,
	"NANTS"		=> 4,
	"PKTNCHAN"=> 8,
	"OBSNCHAN"=> 1024,
	"OBSFREQ" => 1420,
	"OBSBW" 	=> 128,
	"TBIN" 		=> 2.0f-6,
	"PKTIDX" 	=> 0
)

GuppiRaw.calcBlockShape(header8)
log(2.0, GuppiRaw.calcBlockSize(header8))

@assert header4["BLOCSIZE"] == GuppiRaw.calcBlockSize(header4)
@assert header8["BLOCSIZE"] == GuppiRaw.calcBlockSize(header8)

# dimIdx has dimensions ["NPOL", "PIPERBLK", "PKTNCHAN", "NSTRM", "NANTS"]
function genComplex(blkDimIdx, realMag=1)
	# r = (blkDimIdx[1]%2==0 ? -1 : 1)*realMag
	# i = (blkDimIdx[3]%2==0 ? -1 : 1) * blkDimIdx[2]
	# if blkDimIdx[1]%2 == 0 || blkDimIdx[2] != 1 || blkDimIdx[4] != 1
	# if blkDimIdx[2] != 1 || blkDimIdx[4] != 1
	# if blkDimIdx[5] != 2 || blkDimIdx[2] != 1
	if blkDimIdx[2] != 1
		return 0 + 0im
	end
	# return rand(Int8) + rand(Int8)im
	return 1+1im#r + i*1im
end

buffio4 = IOBuffer(maxsize=2^27)
buffio8 = IOBuffer(maxsize=2^28)

print("generating ")
for dimIdx in CartesianIndices(GuppiRaw.calcBlockShape(header8))#iterate through the BLOCSIZE, through all of the axes
	# dimIdx has dimensions ["NPOL", "PKTNTIME", "PKTNCHAN", "PKTPERBLK", "NSTRM", "NANTS"]
	cmplx = genComplex(dimIdx, 1)
	cmplx4 = GuppiRaw.Int4(real(cmplx)) + GuppiRaw.Int4(imag(cmplx))*1im
	# GuppiRaw.writeComplex!(buffio4, cmplx4, 4)
	GuppiRaw.writeComplex!(buffio8, cmplx4, 8)
		# GuppiRaw.writeComplex!(fio,c_sig[i],4)
end
print("block\n")

buffReadio4 = IOBuffer(take!(buffio4), write=false, read=true)
buffReadio8 = IOBuffer(take!(buffio8), write=false, read=true)

# peek some byte values
# sample3 = take!(buffReadio4)[1:10]
# sample3_8 = take!(buffReadio8)[1:20]

fio4 = open(fp4,"w")
fio8 = open(fp8,"w")

# @time begin
		# fio = Mmap.mmap(fd, UInt8, (128*(15*80 + 128*1024*1024),), 0)
		for bIdx in 1:1
				# GuppiRaw.writeHeader!(fio4,header4)
				GuppiRaw.writeHeader!(fio8,header8)

				# write(fio4,take!(buffReadio4))
				write(fio8,take!(buffReadio8))
				
				# header4["PKTIDX"] += header4["PIPERBLK"]
				header8["PKTIDX"] += header8["PIPERBLK"]
		end
	# end
# end

close(fio8)
close(fio4)
