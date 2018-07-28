--[[
	Load file and get rid of extra blank lines
--]]
local mat = assert( require("mat"), "Failed to load mat.lua")
local loader = {}

-- figure out how many columns there are
function _DetermineNumColumns(line, filePattern)
	local count=0
	for match in string.gmatch(line,filePattern) do
		count = count + 1
	end
	
	return count
end

local function _DetermineRowsAndColumns( fileName, filePattern, numHeaderLines, RowValidator )
	local file,errStr = io.open(fileName,"r")
	if not file then
		error(errStr)
	else
		-- skip the header lines
		for i=1,numHeaderLines do file:read() end
		
		local rowCount,colCount = 0, 0
		for line in file:lines() do		
			-- don't count bad rows
			if not RowValidator or (RowValidator and RowValidator(line,filePattern)) then					
				if rowCount == 0 then
					colCount = _DetermineNumColumns(line, filePattern)
				end
				rowCount = rowCount + 1
			end
		end
		file:close()
		return rowCount, colCount
	end
end

local function _AllocateVector( fileName, filePattern, numHeaderLines, RowValidator)
	local rows,cols = _DetermineRowsAndColumns( fileName, filePattern, numHeaderLines, RowValidator )
	local matrix = mat.New( rows, cols)
	matrix.headerLines = {}
	assert(matrix,"Error : No data were allocated for " .. fileName)
	return matrix
end

function loader.LoadFile(args)

	local fileName 		= args.fileName
	--local filePattern = args.pattern or "([-]?%d+%.?%d*[Ee]?[+-]?%d*)"
	-- we can use different conversion functions for different columns
	local converter		= args.converter or function(colIdx, x) return tonumber(x) end
	-- this pulls out any value between commas.  There is a comma before and a comma after, because
	-- in xcel .csv files, the last column does not have a comma after it
	local filePattern 	= args.pattern or ",-([%d%a%s%+%-_%.%(%)%.%/%<%>%%%.%:%[%]]*),?"
	local showDataLoad 	= args.showDataLoad
	local numHeaderLines= args.numHeaderLines
	local RowValidator  = args.RowValidator
	local skipColCountAssert = args.skipColCountAssert

	local data = _AllocateVector( fileName, filePattern, numHeaderLines, RowValidator )
	local numCols = data:Cols()
	local file,errStr = io.open(fileName,"r")
	if not file then
		return false, errStr
	else
		local row = 1
		for line in file:lines() do
			-- save the number of header lines that were requested
			if #data.headerLines < numHeaderLines then
				table.insert(data.headerLines,line)
			else
				if showDataLoad then io.write(row .. "\t") end
				if not RowValidator or (RowValidator and RowValidator(line,filePattern)) then					
					-- parse the line for the given filePattern
					local col=1
					for match in string.gmatch(line,filePattern) do					
						local num = converter(col, match)
						if showDataLoad then io.write(tostring(num) .. "\t") end
						if not skipColCountAssert then 
							assert(col <= numCols,string.format("Row %d has more columns, %d,  than the first line, which had %d",row,col,numCols))	
						end											
						data[row][col] = num
						col = col + 1
					end
					if showDataLoad then io.write("\n") end	
					row = row + 1
				end				
			end
		end
		file:close()
	end
	return data
end

--[[
local data = LoadFile{fileName = "TSQ Endura SRM Table.csv", numHeaderLines = 1, showDataLoad = true, 
	converter = function(col, x)
		if col == 1 or col == 4 then
			return x
		else
			return tonumber(x)
		end
	end}
data:Print()
--]]

return loader