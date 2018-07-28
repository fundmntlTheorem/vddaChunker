--[[

--]]

local utils = assert(require("utils"))
local fun = assert(require("fun"))

local mat = {}

local mtbl = {
	__call = function(tbl,...)
		-- allow to retrive a value using the syntax vec(n) instead of vec[n][1]
		local args = table.pack(...)
		if #args == 1 then 
			return tbl[ args[1] ][1]
		elseif #args == 2 then
			return tbl[ args[1] ][ args[2] ]
		else
			error("mat( i, [j] ) should specify row or row,col dimensions")
		end		
	end,
	__eq = function(lhs, rhs)
		return mat.Equals(lhs,rhs)
	end,
	__lt = function(lhs, rhs)
		
		-- the compiler can switch the order of arguments
		-- for < versus > if you put a 'not' in front of an expression
		if type(lhs) == "number" then
			return mat.GreaterThan(rhs, lhs)
		else
			return mat.LessThan(lhs,rhs)
		end
	end,
	__le = function(lhs, rhs)
		if type(lhs) == "number" then
			return mat.GreaterThanEquals(rhs, lhs)
		else
			return mat.LessThanEquals(lhs, rhs)
		end
	end,
	__add = function(a, b)
		return mat.Add(a, b)
	end,
	__sub = function(a, b)
		return mat.Sub(a,b)
	end,
	-- look up the methods of the mat table
	__index = function(tbl, key)
		return mat[key]
	end,	
}

--[[
	mat.Resize
	Resize the matrix to the indicated dimensions.  All current data are cleared
					+-------------------- matrix object
					|		+------------ first matrix dimension	
					|		|	+------- (optional) second dimension of matrix
					|		|	|	+--- (optional) the initial value to set everyon to
					|		|	|	|
					v		v	v	v	
--]]
function mat.Resize(matrix, m, n, initVal)
	
	utils.CheckArgs(matrix, {"table"},
					m, {"number", function(arg) return arg > 0,"rows \'m\' must be > 0" end},
					n, {function(arg) return arg == nil or (type(arg) == "number" and arg > 0),
							"cols '\n\' must be > 0 or nil" end},
					initVal, {function(arg) return arg == nil 
								or type(arg) == "number","initial value must be a number" end})
	n = n or 1
	initVal = initVal or 0
	
	local mi,ni = mat.Size( matrix )
	
	-- clear rows we don't need
	for rIdx=m+1,mi do matrix[rIdx] = nil end
	
	for rIdx=1,m do
		
		matrix[rIdx] = {}
		
		for cIdx=1,n do
			matrix[rIdx][cIdx] = initVal
		end			
		
		-- clear unneeded rows
		for cIdx=n+1,ni do matrix[rIdx][cIdx] = nil end
	end
	
	return matrix
end

--[[
	mat.New
	Create a new matrix
	
	one argument
		number : create a vector with specified size
		table : copy constructor
	two arguments
		number x number : create a matrix with specified size	
--]]
function mat.New(...)
	
	-- give this table the mat class metatable
	local o = {{}}
	
	local args = table.pack(...)
	if #args == 0 then
    
		-- no arguments passed, just make a valid matrix
		mat.Resize(o, 1)
    
	elseif #args == 1 then
		
		if type(args[1]) == "number" then
			-- with just one dimension specified, we are a vector	
			mat.Resize(o, args[1])
			
		elseif type(args[1]) == "table" then
      
			-- here we have passed some numbers, so we are either copying a table to 
			-- a matrix object, or copying a matrix object to another matrix object
			-- make it a valid size
			o = mat.Copy(o, args[1] )
			
		else
			error("mat.New : first argument has unallowed type : \'"..type(args[1]).."\'")
		end
	elseif #args == 2 then
		
		-- with two dimensions specified, we are creating a matrix	
		if type(args[1]) == "number" and type(args[2]) == "number" then
			
			mat.Resize(o, args[1], args[2])
		else
			
			error("mat.New : unallowed argument types : \'"..type(args[1]) .. "\' , \'" .. type(args[2]).."\'")
		end		
	elseif #args > 2 then
		
		error("mat.New : unallowed number of arguments : \'"..#args.."\'")
	end
		
	setmetatable(o, mtbl)	
	return o
end

--[[mat.IsMatrix
	A table is an object of the mat type if it has the same metatable.
	Also return an error string that can be used to throw with error()
--]]
function mat.IsMatrix(matrix)
	return getmetatable(matrix) == mtbl, "argument is not a matrix object"
end

--[[mat.IsColVector
	check that the column dimension is equal to 1
--]]
function mat.IsColVector(matrix)
	
	if not mat.IsMatrix(matrix) then return false end
	return mat.Cols(matrix) == 1
end

--[[mat.IsRowVector
	check that the row dimension is equal to 1
--]]
function mat.IsRowVector(matrix)
	
	if not mat.IsMatrix(matrix) then return false end
	return mat.Rows(matrix) == 1
end

--[[
	Check for if the argument is a column vector or a 1D table
--]]
function mat.IsVectorOr1DTable( x ) 
	return mat.IsColVector(x) or 
			(type(x)=="table" and type(x[1])=="number"),
			"Argument is not a 1D table or matrix" 
end

--[[
	Check for if the argument is a 2D matrix or table
--]]
function mat.Is2DMatrixOr2DTable( x )
	return mat.Cols( x ) > 1 and mat.Rows( x ) > 0, 
		"Argument is not a 2D table or 2D matrix"	
end

--[[
	This can be useful, since both matrix and 2D table are 2D tables
--]]
function mat.IsMatrixOr2DTable( x )
	return mat.IsMatrix( x ) or mat.Cols( x ) > 1,
		"Argument is not a 2D table or matrix"
end

--[[
	Check if argument is either a mat object or a table
--]]
function mat.IsMatrixOrTable( x )
	return mat.IsMatrix(x) or type(x) == "table",
		"Argument is not a table or a matrix"
end

function mat.Is1DTable( x )
	return type(x) == "table" and type(x[1]) == "number"
end

--[[
	mat.Size
	Get the dimensions of the matrix
	both dimensions should exist because of how
	mat.Resize works
--]]
function mat.Size(matrix)
	
	if mat.IsMatrix( matrix ) then
		return #matrix,#matrix[1]		
	else
		utils.CheckArgs( matrix, {"table"})
		
		if type(matrix[1]) == "table" then
			-- return maximum column size
			return #matrix, fun.Reduce( math.max, fun.Map( function(x) return #x end, matrix ) )			
		else 
			-- this is a 1D table
			return #matrix,1
		end		
	end	
end

--[[
  mat.Total
  total number of elements in the matrix
--]]
function mat.Total(matrix)
	local m,n = mat.Size( matrix )
	return m * n
end

--[[mat.AreSameSize	
	Compare two matrix sizes
--]]
function mat.AreSameSize(lhs, rhs)
	
	utils.CheckArgs(lhs, {"table"},
					rhs, {"table"})
				
	local ml,nl = mat.Size(lhs)
	local mr,nr = mat.Size(rhs)
	
	return ml == mr and nl == nr,"matrices aren't the same size : \'"..ml.."\' ~= \'"..mr
									.."\' or \'"..nl.."\' ~= \'"..nr.."\'"
end

--[[mat.Copy
	copy the values from a mat object or 1D or 2D table of numbers
	into the lhs object
					+--------------- mat object to copy numbers into
					|	+--- mat object or table to copy numbers from
					|	|
					v	v
--]]
function mat.Copy(lhs, rhs)
	
	utils.CheckArgs(lhs, {"table"},
					rhs, {"table"})
	
	if mat.IsMatrix(rhs) then
		
		-- copy a matrix to a matrix, make sure the sizes are the same though
		if not mat.AreSameSize(lhs,rhs) then 
			mat.Resize(lhs,mat.Size( rhs )) 
		end
		
		for i,row in ipairs(rhs) do
			for j,val in ipairs(row) do
				lhs[i][j] = val
			end
		end		
    
	else
    
		-- in this case, the rhs is just a table, it should be either 1D or 2D
		if type(rhs[1]) == "table" then
      
			-- make sure all the dimensions of the tables are the same
			local sizes = fun.Map( function(x) return #x end, rhs )
			assert( fun.All( function(x) return x == #rhs[1] end, sizes ), "All rows of 2D table must be the same size" )
      
			-- resize the lhs to be the same size as the rhs
			mat.Resize( lhs, #rhs, sizes[1] )
      
			-- copying 2D table
			for i,row in ipairs(rhs) do
				
				assert(type(row) == "table","\nrow \'"..i.."\' of rhs is not a table, is a \'"..type(row).."\'")
				
				for j,val in ipairs(row) do
					assert(type(val) == "number","\nitem at index \'"..i.." , "
						..j.."\' is not a number, is a \'"..type(val).."\'")
					
					lhs[i][j] = val
				end
			end			
			
		elseif type(rhs[1]) == "number" then 
			
			-- copying a vector
			local m,n = mat.Size(lhs)
			if m ~= #rhs or n ~= 1 then mat.Resize(lhs, #rhs, 1) end
						
			for i,val in ipairs(rhs) do
				assert(type(val) == "number","\nindex \'"..i.."\' of the rhs is not a number, is a \'"..type(val).."\'")
				lhs[i][1] = val
			end	
			
		end	
	end
	return lhs
end

--[[mat.All 
	true if all items in a matrix meet a condition,
	otherwise false
					+--------------- mat object
					|		+------- function to operate on each item in matrix
					|		|
					|		|
					v		v
--]]
function mat.All(matrix, func)
	
	utils.CheckArgs( func, {"function"} )
	
	if mat.IsMatrixOr2DTable( matrix ) then
		
		-- outer Reduce the inner results
		return fun.All( function(a) 
				
				-- apply to each value in the row
				return fun.All( function(b) 
						
							return func( b )
						end,a)		
			end, matrix)		
	else
		
		utils.CheckArgs( matrix, {mat.Is1DTable} )
		
		return fun.All( function(b) 
		
			return func( b )
		end, matrix)			
	end
end

--[[mat.Any
	If any of the items in the matrix meet a condition
	
					+----------- mat object
					|	+------- function to operate on each item in matrix
					|	|
					|	|	
					v	v	
--]]
function mat.Any(matrix,func)
	
	utils.CheckArgs( func, {"function"} )
	
	if mat.IsMatrixOr2DTable( matrix ) then
		
		-- outer checks all the inner results
		return fun.Any( function(a) 
				
				-- apply to each value in the row
				return fun.Any( function(b) 
						
							return func( b )
						end,a)		
			end,matrix)		
	else
		
		utils.CheckArgs( matrix, {mat.Is1DTable} )
		
		return fun.Any( function(b) 
		
			return func( b )
		end, matrix)			
	end
end

--[[
	Return 1D table version of 2D structure
	ex {{1, 2}, {3, 4}} -> {1, 2, 3, 4}
	This can be convenient if you wanted to zip
	two 2D tables together.
--]]
function mat.Flatten( matrix )

	utils.CheckArgs( matrix, {mat.IsMatrixOrTable} )

	local v = {}
	
	if mat.IsMatrixOr2DTable( matrix ) then
	
		fun.Map( function(a)
				fun.Map( function(b)
					table.insert(v, b)
				end,a)
			end,
			matrix)
	
	else
	
		utils.CheckArgs( matrix, {mat.Is1DTable} )
		
		fun.Map( function(b)
					table.insert(v, b)
				end,matrix)
	end

	return v
end

--[[mat.CompareAll
	General purpose comparing function that checks
	the arguments and allows to compare each value of
	a matrix against a number or against each value
	of another matrix
--]]
function mat.CompareAll(lhs, rhs, compareFunc)
	
	utils.CheckArgs( lhs, {mat.IsMatrixOrTable} )
	
	compareFunc = compareFunc or operator.eq
	
	local areAll
	if type(rhs) == "number" then
		
		utils.CheckArgs(lhs, {mat.IsMatrixOrTable})
		
		-- every value in lhs must be the number rhs
		areAll = mat.All( lhs, function(x)
			return compareFunc(x, rhs)
		end)
		
	else
		
		utils.CheckArgs({lhs, rhs}, {function(arg) return mat.AreSameSize(arg[1], arg[2]) end})		
		
		-- every value in lhs and rhs, must be the same
		areAll = fun.All( function( tbl )
			return compareFunc( tbl[1], tbl[2] )
		end,
		fun.Zip( mat.Flatten( lhs ), mat.Flatten( rhs ) ))
	
	end
	
	return areAll,"lhs and rhs are not the same"	
end

--[[mat.CompareAny
	returns true if any value meets condition, otherwise
	returns false.  This version assures that the arguments
	meet the conditions of lhs is mat object, rhs is a number
	or a mat object
--]]
function mat.CompareAny(lhs, rhs, compareFunc)
	
	utils.CheckArgs( lhs, {mat.IsMatrixOrTable} )
	
	compareFunc = compareFunc or operator.eq
	
	local isAny
	if type(rhs) == "number" then
	-- every value in lhs must be the number rhs
	
		isAny = mat.Any(lhs, function(x)
			return compareFunc(x, rhs)
		end)
		
	else
		
		utils.CheckArgs( {lhs, rhs}, {function(arg) return mat.AreSameSize(arg[1], arg[2]) end})		
	
		-- every value in lhs and rhs must be the same
		isAny = fun.Any( function( tbl )
			return compareFunc( tbl[1], tbl[2] )
		end,
		fun.Zip( mat.Flatten( lhs ), mat.Flatten( rhs ) ))
	
	end
	
	return isAny, "lhs and rhs share no values"
end

--[[mat.Equals
	return true if all elements of two matrices are the same,
	otherwise return false
--]]
function mat.Equals(lhs, rhs)
	return mat.CompareAll(lhs, rhs, function(l, r) return l == r end)
end

function mat.LessThan(lhs, rhs)
	return mat.CompareAll(lhs, rhs, function(l, r) return l < r end)
end

function mat.GreaterThan(lhs, rhs)
	return mat.CompareAll(lhs, rhs, function(l, r) return l > r end)
end

function mat.LessThanEquals(lhs, rhs)
	return mat.CompareAll(lhs, rhs, function(l, r) return l <= r end)
end

function mat.GreaterThanEquals(lhs, rhs)
	return mat.CompareAll(lhs, rhs, function(l, r) return l >= r end)
end

--[[mat.SetConst
	Set the matrix to a constant value
						+----------- a mat object
						|		+--- a number to set the matrix to
						|		|
						v		v
--]]
function mat.SetConst(matrix, val)
	
	utils.CheckArgs(matrix, {mat.IsMatrix},
					val, {"number"})
	
	for i,row in ipairs(matrix) do
		for j,_ in ipairs(row) do
			matrix[i][j] = val
		end
	end
	return matrix
end

--[[mat.Col
	Copy all the values in the given column into a new mat object
--]]
function mat.Col(matrix,colIdx)
	
	utils.CheckArgs(matrix, {mat.IsMatrix},
					colIdx, {"number", function(arg) return utils.Frac(arg) == 0 end,
									   function(arg) 
											local m,n = matrix:Size()
											return colIdx <= n
										end})
								
	-- make a new vector that fits all the rows
	local m = matrix:Rows()
	local newMat = mat.New(m)
	
	for i=1,m do
		newMat[i][1] = matrix[i][colIdx]
	end
	
	return newMat
end

--[[mat.Cols
	return the number of columns in the matrix
--]]
function mat.Cols(matrix)
	local m,n = mat.Size( matrix )
	return n	
end

--[[mat.Row
	Copy all the values in the given row into a new mat object
--]]
function mat.Row(matrix,rowIdx)
	
	utils.CheckArgs(matrix, {mat.IsMatrix},
					rowIdx, {"number", function(arg) return utils.Frac(arg) == 0 end,
									   function(arg) 
											local m,n = matrix:Size()
											return rowIdx <= m
										end})
								
	-- make a new vector that fits all the columns
	local n = matrix:Cols()
	local newMat = mat.New(1,n)
	
	for j=1,n do
		newMat[1][j] = matrix[rowIdx][j]
	end
	
	return newMat	
end

--[[mat.Rows
	return the number of rows in the matrix
--]]
function mat.Rows(matrix)
	local m,n = mat.Size( matrix )
	return m
end

--[[mat.SetRandomUniformInt
	set each value to a random integer value between
	optional arguments lower and upper
--]]
function mat.SetRandomUniformInt(matrix,lower,upper)
	
	utils.CheckArgs(matrix, {mat.IsMatrix})
	if lower then utils.CheckArgs(lower, {"number"}) end
	if upper then utils.CheckArgs(upper, {"number"}) end
	
	for i,row in ipairs(matrix) do
		for j,_ in ipairs(row) do
			matrix[i][j] = math.random(lower, upper)
		end
	end
  
	return matrix
end

function mat.SetRandomUniformReal(matrix,...)
	
end

function mat.SetRandomDistribution(matrix,...)
	
end

--[[
  Set to a series of numbers.
  Give start, stop, number of steps, or
       start, stop, step size
  resize into a column vector of appropriate size
--]]
function mat.SetSeries(matrix, args)
	
	utils.CheckArgs(matrix, {mat.IsMatrix},
                  args, {function(x) return args.steps or args.stepSize and not (args.steps and arg.stepSize), 
                           "Either steps or stepSize must be specified" end},
                  args.start, {"number"}, args.stop, {"number"})
  
	local start,stop = args.start, args.stop
	local steps, stepSize = args.steps, args.stepSize
	if steps then
		-- the stepsize should ensure that we start and stop on the proper values
		stepSize = ( stop - start ) / (steps - 1)
	else
		steps = ( stop - start ) / stepSize + 1
	end  

	matrix:Resize( steps )

	for i=1,steps do
		matrix[i][1] = start + stepSize * (i-1)
	end

	return matrix
end

--[[mat.MathOperation
	apply a function to each item in rhs, where that
	function takes the lhs item, and either the rhs item
	if it is a mat object, or else it is a number
--]]
function mat.MathOperation(lhs, rhs, op)
	
	utils.CheckArgs(lhs, {mat.IsMatrix},
					op,  {"function"})
		
	local result = mat.New(lhs:Size())
	
	if type(rhs) == "number" then
		
		for i,row in ipairs(lhs) do
			for j,val in ipairs(row) do
				result[i][j] = op(val, rhs)
			end
		end			
	else
		utils.CheckArgs(rhs, {mat.IsMatrix},
						{lhs, rhs}, {function(arg) return mat.AreSameSize(arg[1], arg[2]) end})
		
		
		for i,row in ipairs(lhs) do
			for j,val in ipairs(row) do
				result[i][j] = op(val, rhs[i][j])
			end
		end		
	end
	return result
end

function mat.Add(lhs, rhs)
	return mat.MathOperation(lhs, rhs, function(a, b) return a + b end)
end

function mat.Sub(lhs, rhs)
	return mat.MathOperation(lhs, rhs, function(a, b) return a - b end)
end

function mat.Mult(lhs, rhs)
	
end

function mat.Div(lhs, rhs)
	
end

--[[
	compute the l-norm of the matrix values
--]]
function mat.Norm(matrix, l)
	l = l or 2
	
	utils.CheckArgs( 	matrix, {"table"},
						l, {"number", function(x) return x > 0,"norm l must be > 0" end} )
	
	local norm = 0
	if mat.IsMatrix( matrix ) or type(matrix[1]) == "table" then
		
		for i,row in ipairs(matrix) do
			for j,val in ipairs(row) do
				norm = norm + val^l
			end
		end	
	else
		
		if type(matrix[1]) == "number" then
			for i,val in ipairs(matrix) do
				norm = norm + val^l
			end
		else
			error("mat.Norm expects a mat object or a 1D or 2D table of numbers")
		end		
	end
	
	return norm^(1/l)
end

--[[
	Dot product of two vectors.  Arguments must be 
	1D tables or column vector mat objects, or a mix of these
--]]
function mat.Dot( lhs, rhs )
	
	utils.CheckArgs( lhs, {mat.IsVectorOr1DTable}, rhs, {mat.IsVectorOr1DTable}, 
		{lhs, rhs}, {function(x) return #x[1] == #x[2],"lhs and rhs must have the same size" end} )
	
	local function Getter( tbl )
		return mat.IsMatrix(tbl) and function(i) return tbl[i][1] end 
			or function(i) return tbl[i] end
	end
	
	local getLHS, getRHS = Getter(lhs), Getter(rhs)
	
	local dot = 0
	for i=1,#lhs do
		dot = dot + getLHS(i) * getRHS(i)
	end
	
	return dot
end

--[[
	calculates a 0-1 number such that 0 is orthogonal vectors,
	and 1 is the same vector
--]]
function mat.CosineSimilarity( lhs, rhs )
	
	utils.CheckArgs( lhs, {IsVectorOr1DTable}, rhs, {IsVectorOr1DTable}, 
		{lhs, rhs}, {function(x) return #x[1] == #x[2],"lhs and rhs must have the same size" end} )
	
	return mat.Dot(lhs, rhs) / ( mat.Norm( lhs ) * mat.Norm( rhs ) )	
end

--[[
	sum up all the values in the matrix
--]]
function mat.Sum(matrix)
  
	local sum = 0
	if mat.IsMatrixOr2DTable( matrix ) then
		for _,row in ipairs(matrix) do
			for _,val in ipairs(row) do
				sum = sum + val
			end
		end
	elseif mat.IsVectorOr1DTable( matrix ) then
		for _,val in ipairs(matrix) do
			sum = sum + val
		end
	else
		error("Unrecognized argument type " .. type(matrix))
	end
	
	return sum  
end

--[[
	compute the average of all the values
--]]
function mat.Average(matrix)
  	
	local total = mat.Total( matrix )
	return total == 0 and 0 or mat.Sum( matrix ) / total
end

--[[
	Compute the standard deviation of all the values
	and return the mean as well, since it is
	calculated as a part of this
--]]
function mat.StdDev(matrix)
   
	local var, mu = 0, 0
	local n = mat.Total( matrix )
	if n <= 0 then return var, mu end
		
	if mat.IsMatrixOr2DTable( matrix ) then
 
		for _,row in ipairs(matrix) do
			for _,x in ipairs(row) do
				var = var + x^2
				mu = mu + x
			end
		end 
	elseif mat.IsVectorOr1DTable( matrix ) then
		
		for _,x in ipairs(matrix) do
			var = var + x^2
			mu = mu + x
		end
	else
		error("Unrecognized argument type " .. type(matrix))
	end

	mu = mu / n
	var = var / n

	return math.sqrt( (var - mu^2) ), mu
end
	
function mat.StdDevSample(matrix)
   
	local var, mu = 0, mat.Average( matrix )
	local n = mat.Total( matrix )
	if n <= 0 then return var, mu end
		
	if mat.IsMatrixOr2DTable( matrix ) then
 
		for _,row in ipairs(matrix) do
			for _,x in ipairs(row) do
				var = var + (x - mu)^2
			end
		end 
	elseif mat.IsVectorOr1DTable( matrix ) then
		
		for _,x in ipairs(matrix) do
			var = var + (x - mu)^2
		end
	else
		error("Unrecognized argument type " .. type(matrix))
	end

	var = var / (n-1)

	return math.sqrt( var ), mu
end

--[[
	Return the median of all the values in the matrix
--]]
function mat.Median( matrix )
	
	assert( mat.Total( matrix ) > 0 )
	
	local flat = mat.Flatten( matrix )
	table.sort(flat)
	
	return flat[ math.floor( #flat / 2 + 0.5) ]	
end

--[[
	Helper function for the Min and Max methods
--]]
local function _MinOrMax( func, matrix )
	
	if mat.IsMatrixOr2DTable( matrix ) then
	
		local minOrMax, mr, mc = matrix[1][1], 1, 1
		for i,row in ipairs(matrix) do
			for j,x in ipairs(row) do
				if func(x, minOrMax) then
					minOrMax, mr, mc = x, i, j
				end
			end
		end
	
		return minOrMax, mr, mc
	elseif mat.Is1DTable(matrix) then
	
		local minOrMax, mr = matrix[1], 1
		for i,x in ipairs( matrix ) do
			if func(x , minOrMax) then
				minOrMax, mr = x, i
			end
		end
		
		return minOrMax, mr
	end	
end

--[[
	return minimum value in matrix
	and the indices where it occurs
--]]
function mat.Min( matrix )
	
	return _MinOrMax( operator.lt, matrix )
end

--[[
	return maximum value in matrix
	and the indices where it occurs
--]]
function mat.Max( matrix )
	
	return _MinOrMax( operator.gt, matrix )
end

--[[
	Add a value to a vector, or a row to a matrix
--]]
function mat.Push( matrix, val )
	
	if mat.Is1DTable( matrix ) then
		
		utils.CheckArgs(val, {"number"})
		
		table.insert(matrix, val)
		
	elseif mat.IsColVector( matrix ) then
		
		utils.CheckArgs(val, {"number"})
		
		table.insert(matrix, {val})
		
	elseif mat.Is2DMatrixOr2DTable( matrix ) then
		
		utils.CheckArgs(val, {mat.Is1DTable})
		local m,n = mat.Size( matrix )
		local mi,ni = mat.Size( val )
		assert(n == mi,"Pushed item has cols \'" .. ni .. "\' ~= \'" .. n .. "\'")
		table.insert(matrix, val)
		
	else
		error("unrecognized argument type \'" .. type(matrix) .. "\'")
	end
end

--[[
	Return the subset of the matrix
	If the matrix is actually a 1D table, eg {1,2,3},
	then i and j are the start and stop indices at which to copy.
	
	If the matrix is 2D ( table or mat ), then i, j are the start indices,
	and rows, cols are the size of the subsection to copy
--]]
function mat.SubSet( matrix, i, j, rows, cols )
	
	utils.CheckArgs( i, {"number"}, j, {"number"})
	
	local t = {}
	if mat.Is1DTable( matrix ) then
		
		assert(i > 0 and i < #matrix and j > 0 and j >= i and j<= #matrix,
			"SubSet indices must be inside array bounds")
		
		-- i and j are start and stop
		local m=1
		for r=i,j do
			t[m] = matrix[r]
			m = m + 1
		end
		
	elseif mat.IsMatrixOrTable( matrix ) then
		
		utils.CheckArgs( rows, {"number"}, cols, {"number"} )
		
		local mf,nf = i + rows, j + cols
		local m,n = mat.Size( matrix )
		assert(mf <= m and nf <= n,"start + subsize must be < size of matrix")		
		
		t = mat.New(rows, cols)
		
		m,n = 1,1
		for r = i,rows do
			for c = j,cols do
				t[m][n] = matrix[r][c]
				n=n+1
			end
			m=m+1
		end
		
	else
		error("unrecognized argument type \'" .. type(matrix) .. "\'")
	end
	
	return t
end

function mat.Print(matrix,fmt)	
	
	fmt = fmt or "%10.3f"
	
	local writers = {	string = function(val) return val end,
						number = function(val) return string.format(fmt,val) end}
	
	if mat.IsMatrix( matrix ) or mat.IsMatrixOr2DTable( matrix ) then
		
		for i,row in ipairs(matrix) do
			for j,val in ipairs(row) do
				io.write(writers[type(val)](val) .. "\t")
			end
			io.write("\n")
		end	
		
	elseif mat.Is1DTable( matrix ) then
		
		for _,val in ipairs(matrix) do
			print(writers[type(val)](val))
		end
		
	end
end

do
math.randomseed(os.time())

local ret, err

---------------------------------------------
-- test the mat.New
do
	-- pass 0 args
	ret,err = utils.ProtectedRun(mat.New)
	assert(ret)
		
	-- pass 1 arg, make a vector
	ret,err = utils.ProtectedRun(mat.New,5)
	assert(ret(3) == 0)
	-- 1D table
	ret,err = utils.ProtectedRun(mat.New,{1,2,3})
	assert(ret)
	ret,err = utils.ProtectedRun(mat.New,{{1,2},{3,4}})
	assert(ret)
	-- shouldn't be garbage input
	ret,err = utils.ProtectedRun(mat.New,"bad")
	assert(not ret)
	ret,err = utils.ProtectedRun(mat.New,"bool")
	assert(not ret)
	
	-- 3 args is too many
	ret,err = utils.ProtectedRun(mat.New,1,2,3)
	assert(not ret,err)

	-- 2nd arg has to be a number
	ret,err = utils.ProtectedRun(mat.New,1,{})
	assert(not ret,err)
	ret,err = utils.ProtectedRun(mat.New,1,"bad")
	assert(not ret,err)

	-- 2 args with numbers
	ret,err = utils.ProtectedRun(mat.New,5,5)
	assert(ret,err)
	assert(ret[3][3] == 0)
end

---------------------------------------------
-- test mat.Resize
do
local m,n
local vec = mat.New(5) 
-- have to pass a table
ret, err = utils.ProtectedRun(mat.Resize)
assert(not ret)
ret, err = utils.ProtectedRun(mat.Resize,3)
assert(not ret,err)

-- passing a 1D table
ret, err = utils.ProtectedRun(mat.Resize,{},3)
assert(ret,err)
m,n = mat.Size(ret)
assert(m==3 and n==1)

ret, err = utils.ProtectedRun(mat.Resize,vec,3)
m,n = mat.Size(vec)
assert(m==3 and n==1)

ret, err = utils.ProtectedRun(mat.Resize,{},3,2)
assert(ret,err)
m,n = mat.Size(ret)
assert(m==3 and n==2)
ret, err = utils.ProtectedRun(mat.Resize,vec,3,2)
m,n = mat.Size(vec)
assert(m==3 and n==2)

ret, err = utils.ProtectedRun(mat.Resize,{{},{}},3)
assert(ret,err)
m,n = mat.Size(ret)
assert(m==3 and n==1)

ret, err = utils.ProtectedRun(mat.Resize,{{},{}},3,5)
assert(ret,err)
m,n = mat.Size(ret)
assert(m==3 and n==5)
ret, err = utils.ProtectedRun(mat.Resize,vec,3,5)
m,n = mat.Size(vec)
assert(m==3 and n==5)

ret, err = utils.ProtectedRun(mat.Resize,vec,-3)
assert(not ret,err)

ret, err = utils.ProtectedRun(mat.Resize,vec,3)
assert(ret,err)

ret, err = utils.ProtectedRun(mat.Resize,vec,3,3)
assert(ret,err)

ret, err = utils.ProtectedRun(mat.Resize,vec,3,-3)
assert(not ret,err)

ret, err = utils.ProtectedRun(mat.Resize,vec,3,"badarg")
assert(not ret,err)

ret, err = utils.ProtectedRun(mat.Resize,vec,3,3,1)
assert(ret,err)

ret, err = utils.ProtectedRun(mat.Resize,vec,3,3,"badarg")
assert(not ret,err)

end

---------------------------------------------
-- test the mat.Size
do
	local m1 = mat.New(3,2)
	local m2 = mat.New(m1:Size())
	local m,n = m2:Size()
	assert(m == 3 and n == 2)
	
	m,n = mat.Size{1}
	assert(m==1 and n==1)
	
	m,n = mat.Size{{1},{1}}
	assert(m==2 and n==1)
	
	m,n = mat.Size{{1,2},{1}}
	assert(m==2 and n==2)
	
	m,n = mat.Size{1,2}
	assert(m==2 and n==1)
end

---------------------------------------------
-- test the mat.Total
do
	assert(mat.Total( {1, 2} ) == 2)
	assert(mat.Total{ {1, 2}, {2, 3} } == 4)
end

---------------------------------------------
-- test the mat.AreSameSize
do
	assert(mat.AreSameSize( mat.New(5), {{1},{2},{3},{4},{5}}))
	assert(mat.AreSameSize( mat.New(2,2), {{1,1},{2,2}}))
end

---------------------------------------------
-- test the mat.Copy
do
	
local m1 = mat.New(5):SetConst(8)
local t1 = {1,2,3}

ret,err = utils.ProtectedRun(mat.Copy,m1,-1)
assert(not ret,err)

ret,err = utils.ProtectedRun(mat.Copy,-1,m1)
assert(not ret,err)

ret,err = utils.ProtectedRun(mat.Copy,m1,t1)
assert(ret,err)
assert(m1(1) == 1 and m1(2) == 2 and m1(3) == 3)

local m2 = mat.New(5):SetConst(6)
ret,err = utils.ProtectedRun(mat.Copy,m1,m2)
assert(mat.Equals(m1,m2))

local m3 = mat.New{6, 6, 6, 6, 6}
assert(mat.Equals(m2, m3))

local m7 = mat.New()
mat.Copy(m7, {{6}, {6}, {6}, {6}, {6}})
assert(mat.Equals(m7, m3))

m7:Copy(mat.New{{6, 6, 6, 6, 6}})
assert(mat.Equals(m7, {{6, 6, 6, 6, 6}}))

local m4 = mat.New(2,2):SetConst(4)
local m5 = mat.New{ {4, 4}, {4, 4} }
assert(mat.Equals(m4, m5))

ret,err = utils.ProtectedRun(mat.New, { {4, 4}, {4} } )
assert( not ret, err )

end

---------------------------------------------
-- test the mat.Avg, mat.StdDev, mat.Norm, mat.Dot, mat.CosineSimilarity
do
	local m = mat.New(5):SetConst(1)
	assert(m:Sum() == 5)
	
	m = mat.New():SetSeries{start=0,stop=10,steps=11}
	assert(m:Average() == 5)	
--	assert(utils.TruncateDigits(m:StdDev(),1e12)  == 2.872281323269)
	assert(utils.TruncateDigits(m:Norm(), 1e5) == 19.62141)
	assert(utils.TruncateDigits(m:Dot(m), 1e5) == 385)
	
	-- test with combinations of matrices and tables
	local t = {0,1,2,3,4,5,6,7,8,9,10}
	assert(utils.TruncateDigits(mat.Dot(t,m), 1e5) == 385)
	assert(utils.TruncateDigits(mat.Dot(m,t), 1e5) == 385)
	assert(utils.TruncateDigits(mat.Dot(t,t), 1e5) == 385)
	assert(utils.TruncateDigits(mat.Norm(t), 1e5) == 19.62141)
	
	assert(utils.TruncateDigits(mat.CosineSimilarity( t, t),1e5) == 1)
	assert(utils.TruncateDigits(mat.CosineSimilarity( m, t),1e5) == 1)
	assert(utils.TruncateDigits(mat.CosineSimilarity( t, m),1e5) == 1)
end

---------------------------------------------
-- test the mat.Median
do
	assert(mat.Median( mat.New{1,2,3,4,5} ) == 3 )
	assert(mat.Median( {1,2,3,4,5} ) == 3 )
	assert(mat.Median( mat.New{1,2,3,4} ) == 2 )
	
	assert(mat.Median( mat.New{3,1,5,2,4} ) == 3 )	
	assert(mat.Median( {3,1,5,2,4} ) == 3 )
	
	assert(mat.Median( mat.New{{3,1,5,2,4}} ) == 3 )	
	assert(mat.Median( {{3,1,5,2,4}} ) == 3 )
	
	assert(mat.Median( mat.New{{3,1,5,2,4}, {1,2,3,4,5}} ) == 3 )	
end

---------------------------------------------
-- test the mat.SetConst
do
	local m = mat.New(5)
	ret,err = utils.ProtectedRun(mat.SetConst,m,-1)
	assert(ret,err)
	for i=1,5 do assert(m(i) == -1,"SetConst of a vec doesn't work") end

	m = mat.New(5,5)
	ret,err = utils.ProtectedRun(mat.SetConst,m,-1)
	assert(ret,err)
	for i=1,5 do 
		for j=1,5 do
			assert(m[i][j] == -1,"SetConst of a mat doesn't work")
		end
	end
end

---------------------------------------------
-- test the mat.SetASeries
do
	local m = mat.New():SetSeries{start=1,stop=10,stepSize=1}
	local mh = mat.New{1,2,3,4,5,6,7,8,9,10}
	assert(m == mh)

	m = mat.New():SetSeries{start=1,stop=10,steps=10}
	assert(m == mh)
end

---------------------------------------------
-- test the mat.Equals
do
	
	local m1,m2 = mat.New(5):SetConst(10),mat.New(5):SetConst(10)
	ret,err = utils.ProtectedRun(mat.Equals,m1,m2)
	assert(ret,"Can't tell that two matrices are the same")
		
	m1,m2 = mat.New(5):SetConst(10),mat.New(5):SetConst(5)
	ret,err = utils.ProtectedRun(mat.Equals,m1,m2)
	assert(not ret,"Can't tell that two matrices are the different")

	assert(not (m1 == m2),"metamethod == doesn't work")
	assert(m1 ~= m2,"metamethod ~= doesn't work")

	ret,err = utils.ProtectedRun(mat.Equals,m1,10)
	assert(ret, err)

	ret,err = utils.ProtectedRun(mat.Equals,m1,1)
	assert(not ret,err)

end

---------------------------------------------
-- mat.Cols, and mat.Rows on tables
do
	local t = {{1,2,3,4,5,6,7,8,9,10}}
	assert(mat.Cols(t) == 10)
	assert(mat.Cols({1,2,3,4,5,6,7,8,9,10}) == 1)
	assert(mat.Cols{{1,1}, {2,2}} == 2)
	assert(mat.Cols{{1,1}, {2}} == 2)
end

---------------------------------------------
-- test the mat.GreaterThan and mat.LessThan
do
	local m1, m2 = mat.New(5):SetConst(10), mat.New(5):SetConst(8)

	ret,err = utils.ProtectedRun(mat.GreaterThan, m1, m2)
	assert(ret, err)
	assert(m1 > m2,"can't tell m1 > m2")

	ret,err = utils.ProtectedRun(mat.LessThan, m1, m2)
	assert(not ret, err)
	assert(not (m1 < m2),"can't tell m1 < m2")

	ret,err = utils.ProtectedRun(mat.GreaterThan, m1, 5)
	assert(ret, err)

	ret,err = utils.ProtectedRun(mat.GreaterThan, m1, 15)
	assert(not ret, err)

	assert(m1 < 15)
	assert(not (m1 < 0))
	assert(m1 > 0)
	assert(not (m1 > 100))

end

---------------------------------------------
-- test the mat.GreaterThanEquals and mat.LessThanEquals
do
	
	local m1, m2 = mat.New(5):SetConst(10), mat.New(5):SetConst(8)
	ret,err = utils.ProtectedRun(mat.LessThanEquals, m1, m2)
	assert(not ret, err)
	assert(not (m1 <= m2),"can't tell m1 <= m2")
	assert(m1 <= 100)
	assert(m1 <= 10)
	assert(not (m1 <= 0))

	ret,err = utils.ProtectedRun(mat.GreaterThanEquals, m1, m2)
	assert(ret, err)
	assert(m1 >= m2,"can't tell m1 >= m2")
	assert(m1 >= 0)
	assert(not (m1 >= 100))

	local m1, m2 = mat.New(5):SetConst(8), mat.New(5):SetConst(8)
	ret,err = utils.ProtectedRun(mat.LessThanEquals, m1, m2)
	assert(ret, err)
	assert(m1 <= m2,"can't tell m1 <= m2")

	ret,err = utils.ProtectedRun(mat.LessThanEquals, m2, m1)
	assert(ret, err)
	assert(m2 <= m1,"can't tell m2 <= m1")

	ret,err = utils.ProtectedRun(mat.GreaterThanEquals, m1, m2)
	assert(ret, err)
	assert(m2 >= m1,"can't tell m2 >= m1")

end

---------------------------------------------
-- test the mat.Any
do
	local m1 = mat.New(5):SetConst(5)
	ret,err = utils.ProtectedRun(mat.Any, m1, function(arg) return arg > 0 end)
	assert(ret, err)
	
	local m2 = mat.New{ 1, 2, 3, 4 }
	assert( not mat.CompareAny(m2, 0) )
	assert( mat.CompareAny(m2, 2) )
	
	local m3 = mat.New{ {1, 2}, {3, 4} }
	assert( not mat.CompareAny(m3, 0) )
	assert( mat.CompareAny(m3, 2) )	
	
	ret,err = utils.ProtectedRun(mat.CompareAny, m2, m3 )
	assert(not ret, err)
	
	local m4 = mat.New{ {1, 2}, {3, 5} }
	assert(not mat.CompareAll( m3, m4 ))
	assert(mat.CompareAny( m3, m4 ))
	
	local t1 = { {1, 2}, {3, 4} }
	local t2 = { {1, 2}, {3, 5} }
	assert(not mat.CompareAll( t1, t2 ))
	assert(mat.CompareAny( t1, t2 ))	
	
	assert(mat.CompareAll( {1, 2, 3}, {1, 2, 3} ))
	
end

---------------------------------------------
-- test the mat.Col and mat.Row
do
	local m1 = mat.New(5,5):SetRandomUniformInt(1,10)
	local v1 = m1:Col(2)
	for i,_ in ipairs(v1) do
		assert(v1[i][1] == m1[i][2])
	end	
	
	local v2 = m1:Row(2)
	for j,_ in ipairs(v2) do
		assert(v2[1][j] == m1[2][j])
	end
end

---------------------------------------------
-- test mat.Add
do
	local m1 = mat.New(5):SetConst(1)
	local m2 = mat.New(5):SetConst(2)
	
	local m3 = m1 + m2
	assert(mat.All(m3, function(arg) return arg == 3 end))
	assert(mat.All(m1, function(arg) return arg == 1 end))
	
	local m4 = m1 + 10
	assert(mat.All(m4, function(arg) return arg == 11 end))	  
end

---------------------------------------------
-- test mat.Sub
do
	local m1 = mat.New(5):SetConst(1)
	local m2 = mat.New(5):SetConst(5)
	local m3 = m1 - m2
	
	assert(mat.All(m3, function(arg) return arg == -4 end))
	local m4 = m2 - 10
	assert(mat.All(m4, function(arg) return arg == -5 end))
end

---------------------------------------------
-- test mat.Min and mat.Max
do
	local min,max,i,j
	local m = mat.New(5):SetSeries{start=1,stop=10,steps=10}
	min,i,j = m:Min()
	assert(min == 1 and i == 1 and j == 1)	
	max,i,j = m:Max()
	assert(max == 10 and i == 10 and j == 1)
end

---------------------------------------------
-- test mat.Push
do
	local m = mat.New(5):SetConst( 3 )
	m:Push( 5 )
	
	assert( m:Equals{3, 3, 3, 3, 3, 5} )
	assert( mat.Equals( m, mat.New{3,3,3,3,3,5} ))
	ret,err = utils.ProtectedRun(mat.Push,m,{1})
	assert(not ret, err)
	
	local t = {2, 2, 2, 2}
	mat.Push(t, 3)
	assert( mat.Equals(t, {2, 2, 2, 2, 3}) )
	assert( mat.Equals(t, mat.New{2, 2, 2, 2, 3}) )
	
	local n = mat.New(2, 2):SetConst(1)
	n:Push{3, 3}
	assert( n:Equals{{1, 1}, {1, 1}, {3, 3}} )
	
	ret,err = utils.ProtectedRun(mat.Push,n,{1})
	assert(not ret, err)
end



end

return mat