--[[
	A file for loading spectra that we can use to test variable DDA chunking algorithms.
	There is also a method for evaluating a chunking job
	
	LoadSpectra() : gets a set of spectra over the elution time for a peak.  Takes about 1 second to load
	EvaluateAChunking( spectrum, Chunker ) : evaluates your chunking job for one spectrum
	EvaluateFullChunkingJob( spectra, Chunker ) : the final test, evaluates your chunking job for all spectra
--]]
local mat = require "mat"
require "ut_Math"

local chunkEval = {}
local snStr, rtStr = "%d+", "%d+%.%d+"
local snSearchStr = "S%s+(" .. snStr .. ")%s+.+"
local rtSearchStr = "S%s+" .. snStr .. "%s+("..rtStr..")%s+.+"

local function _GetSpectrum( line )
	local sn = string.match(line, snSearchStr)
	if not sn then return end
	local rt = string.match(line, rtSearchStr)
	assert(rt, "Could find scan number but not rt : " .. line)
	sn, rt = tonumber(sn), tonumber(rt)
	return {sn = sn, rt = rt}	
end

local deisotopeStr, chargeStr, intensityStr, basePeakStr = "%d+%.%d+", "%d", "%d+", "%d+%.%d+"
local deisotopeSearchStr = "P%s+(" .. deisotopeStr .. ")%s+.+"
local chargeSearchStr = "P%s+" .. deisotopeStr .. "%s+(" .. chargeStr .. ")%s+.+"
local intensitySearchStr = "P%s+" .. deisotopeStr .. "%s+" .. chargeStr .. "%s+(" .. intensityStr .. ")%s+.+"
local basePeakSearchStr = "P%s+" .. deisotopeStr .. "%s+" .. chargeStr .. "%s+" .. intensityStr .. "%s+(" .. basePeakStr .. ")%s+.+"

local function _GetPeak(line)
	local deisotopedMass = string.match(line, deisotopeSearchStr)
	-- on spectrum lines this will fail
	if not deisotopedMass then return end
	
	local charge = string.match(line, chargeSearchStr)
	assert(charge,"Could find mass but not charge : " .. line)
	local intensity = string.match(line,intensitySearchStr)
	assert(intensity,"Could find mass and charge but not intensity : " .. line)	
	local basePeak = string.match(line, basePeakSearchStr)
	assert(basePeak,"Could find mass, charge, and intensity but not base peak : " .. line)
	deisotopedMass,charge,intensity, basePeak = tonumber(deisotopedMass), tonumber(charge), tonumber(intensity), tonumber(basePeak)
	
	return { mono = (deisotopedMass + charge) / charge, charge = charge, intensity = intensity, mass = basePeak}
end

--[[
	Load a hardklor file and extract a certain portion of it.  The format is 
	spectra = {
		{sn = scan number, rt = run time (minutes), {mono = monoisotopic mass, charge = charge, intensity = intensity, mass = base peak mass},
		...
	}
	\param fileName a hardklor file name
	\param startSpectrumNumber : the scan number at which to start loading spectra.  must be a valid MS1 spectrum number.
	\param stopSpectrumNumber : the scan number at which to stop loading spectra, inclusive.  must be a valid MS1 spectrum number.
	
--]]
function chunkEval.GetHardKlorFileSpectraSubSet( fileName, startSpectrumNumber, stopSpectrumNumber )
	assert(startSpectrumNumber <= stopSpectrumNumber,"start scan number must be <= stop scan number")
	
	local file,err = io.open(fileName,"r")
	assert(file,err)
	local spectra = {}
	local foundStartLine
	for line in file:lines() do 
		-- just looking for spectrum lines that start with S, 
		local spectrum = _GetSpectrum( line )
		if spectrum and spectrum.sn == startSpectrumNumber then 
			foundStartLine = true			
			spectra[#spectra+1] = spectrum
			break
		end
	end
	assert(foundStartLine,"Couldn't find the starting spectrum line")
	
	local foundStopLine
	for line in file:lines() do
		local spectrum = _GetSpectrum( line )		
		-- a new spectrum has started
		if spectrum then 
			spectra[#spectra+1] = spectrum
		else
			-- this is a peptide line
			local peak = _GetPeak(line)
			assert(peak, "Couldn't parse a peak on this line : " .. line)			
			table.insert(spectra[#spectra], peak)
		end
		
		if spectrum and spectrum.sn == stopSpectrumNumber then 
			foundStopLine = true
			break
		end
	end
	assert(foundStopLine,"Couldn't find the stop spectrum line")
	
	-- sort all the spectra so that later functions can use binary search
	for i,spectrum in ipairs(spectra) do 
		table.sort(spectrum, function(a,b) return a.mass < b.mass end)
	end
	
	return spectra
end

local _firstMass, _lastMass = 375.0, 1500.0
local _massRange = _lastMass - _firstMass
function chunkEval.GetFirstMass() return _firstMass end
function chunkEval.GetLastMass() return _lastMass end

--[[
	Returns a function that makes random chunks.
	This can be used to evaluate the chunk evaluation functions
--]]
local _minRandWidth, _maxRandWidth = 0.4, 50.0
function chunkEval.RandomChunkGenerator()
	local chunks = {}
	local startMass, endMass = _firstMass
	while true do 
		-- generate a random width
		local width = math.random() * (_maxRandWidth - _minRandWidth) + _minRandWidth
		endMass = startMass + width
		if endMass > _lastMass then 
			break 
		end
		local center = (startMass + endMass) / 2.0
		table.insert(chunks, {mass = center, width = width})
		
		startMass = endMass
	end
	endMass = math.min(endMass, _lastMass)
	return chunks
end

local _BinarySearchForMass = math.MakeABinarySearchFunction( function( idx, spectrum ) return spectrum[idx].mass end )

function chunkEval.GetStartIndex( mass, spectrum )
	local idx = _BinarySearchForMass( mass, spectrum )		
	-- if we are greater, then go lower
	while idx >= 1 and spectrum[idx].mass > mass do 
		idx = idx - 1
		moved = true
	end	
	
	-- if we are less than the mass of interest, go higher until we are greater
	while idx >= 1 and idx <= #spectrum and spectrum[idx].mass < mass do 
		idx = idx + 1
	end	
	
	idx = math.max(1, math.min(idx, #spectrum))
	return idx
end

function chunkEval.GetStopIndex( mass, spectrum )
	local idx = _BinarySearchForMass( mass, spectrum )
		
	-- if we were less than the mass of interest, go higher until we are greater
	local moved 
	while idx <= #spectrum and spectrum[idx].mass < mass do 
		idx = idx + 1
		moved = true
	end	
	
	-- if we are greater than the mass of interest, try to go lower until we aren't greater
	while idx <= #spectrum and idx >= 1 and spectrum[idx].mass > mass do 
		idx = idx - 1
	end	
	
	idx = math.max(1, math.min(idx, #spectrum))
	return idx
end

function chunkEval.GetChunkBoundaries( chunk )
	local w2 = chunk.width / 2.0
	return chunk.mass - w2, chunk.mass + w2
end

function chunkEval.GetSpectrumIndicesFromChunk( chunk, spectrum )
	local m0,m1 = chunkEval.GetChunkBoundaries( chunk )
	return  chunkEval.GetStartIndex( m0, spectrum ), chunkEval.GetStopIndex( m1, spectrum )
end

--[[
	Sum the intensity of the peaks in a chunk
--]]
function chunkEval.GetArea( chunk, spectrum )
	local i0, i1 = chunkEval.GetSpectrumIndicesFromChunk( chunk, spectrum )
	local sum = 0.0
	for i=i0, i1 do 
		sum = sum + spectrum[i].intensity
	end
	return sum
end

--[[
	Count the number of peaks per Da in a chunk
--]]
function chunkEval.GetPeaksPerDa( chunk, spectrum )
	local i0, i1 = chunkEval.GetSpectrumIndicesFromChunk( chunk, spectrum )
	-- the number of peaks is the difference in indices + 1
	local peaks = i1 - i0 + 1
	return peaks / chunk.width
end

--[[
	Compute how much of the feasible range, from _firstMass to 
	_lastMass is covered by the chunks.  Assume the chunks are
	sorted in ascending order by mass, done in EvaluateAChunking
--]]
function chunkEval.GetChunksCoverage( chunks )
	local coverage = 0.0
	for i,chunk in ipairs(chunks) do
		coverage = coverage + chunk.width
		if i > 1 then
			local m0 = chunkEval.GetChunkBoundaries( chunk )
			local _,pm1 = chunkEval.GetChunkBoundaries( chunks[i-1] )
			-- subtract any overlap between the chunks
			coverage = coverage - math.max(0.0, pm1 - m0)			
		end
	end
	return coverage / _massRange
end

--[[
	
	\param spectrum : \sa GetHardKlorFileSpectraSubSet for the format of this table
	\param Chunker : a function(spectrum) -> chunks
		For an input spectrum, produces a chunks table with the following format 
		chunks = {
			{mass = center of isolation window, width = full width of the isolation window},
			...			
		}
--]]
function chunkEval.EvaluateAChunking( spectrum, Chunker )
	local chunks = Chunker( spectrum )
	-- make sure the chunks are sorted
	table.sort(chunks, function(a,b) return (a.mass - a.width) < (b.mass - b.width) end)
	
	local areas, peaksPerDa = {}, {}
	for i,chunk in ipairs(chunks) do 
		areas[i] = chunkEval.GetArea( chunk, spectrum )
		peaksPerDa[i] = chunkEval.GetPeaksPerDa( chunk, spectrum )	
	end
	
	local areaStdev,areaMean = mat.StdDev( areas )
	local pksDaStdDev,pksDaMean = mat.StdDev( peaksPerDa )
	local coverageFraction = chunkEval.GetChunksCoverage( chunks )
	
	print(areaMean, pksDaMean, coverageFraction )
end

--[[
	Calls EvaluateAChunking on a set of spectra, and produces some metrics
--]]
function chunkEval.EvaluateFullChunkingJob( spectra, Chunker )
	
	
end

--[[
	Load some spectra from a hardklor file over a certain range.
--]]
function chunkEval.LoadSpectra()
	return chunkEval.GetHardKlorFileSpectraSubSet("20160919_ClassicOTOTfast_HCD_1ug_70min_01_hardklor.hk", 2869, 3122)
end

--local spectra = chunkEval.LoadSpectra()
--chunkEval.EvaluateAChunking( spectra[1], chunkEval.RandomChunkGenerator )

return chunkEval
