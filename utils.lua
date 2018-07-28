local fun = assert( require("Fun"), "Failed to load Fun.lua")

--[[
	Some functions that will check if an argument has the given type
--]]
local knownStringChecks = {
	["string"] = function(arg) return type(arg) == "string", "arg is not a \'string\', but a \'"..type(arg).."\'" end,
	["function"]= function(arg) return type(arg) == "function", "arg is not a \'function\', but a \'"..type(arg).."\'" end,
	["number"] = function(arg) return type(arg) == "number", "arg is not a \'number\', but a \'"..type(arg).."\'" end,
	["nil"] = function(arg) return type(arg) == "nil", "arg is not a \'nil\', but a \'"..type(arg).."\'" end,
	["boolean"] = function(arg) return type(arg) == "boolean", "arg is not a \'boolean\', but a \'"..type(arg).."\'" end,
	["table"] = function(arg) return type(arg) == "table", "arg is not a \'table\', but a \'"..type(arg).."\'" end,
}

local utils = {}

--[[
	Pass a list of arguments and tables to check the arguments against
	i.e.
	arg1, {check1, check2, ...}, arg2, {check1, check2, ...}, ...
	
	each argument will be verified that it meets the certain specifications.
	if the checks are a string, then we'll check that the argument has that type, ex "string".
	if the checks are a function, that function should have interface function(arg) return true/false/nil,"errorString" end
	where the optional error string will be used to pass to error.
	if the check is a table, then the argI must match at least one of the items in the table
--]]
function utils.CheckArgs(...)
	
	local args = table.pack(...) or {}
	
	-- number of arguments has to be even
	assert(#args > 0 and #args % 2 == 0,"\nThe number of arguments has to be an even number > 0, but \'"..#args.."\' arguments were passed")
	
	for i=1,#args,2 do
		
		local argNum = math.floor(i/2+1)
		local argument,conditions = args[i], args[i+1]
		assert(type(conditions) == "table","The check conditions for argument \'" .. argNum .. "\' must be a table"
			.. "\nNot a \'"..type(conditions).."\'")
		
		for cIdx,condition in ipairs(conditions) do
			
			assert(type(condition) == "function" or type(condition) == "string" or type(condition) == "table",
				"Condition \'"..cIdx.."\' for argument \'" .. argNum .. " is \'"
				..type(condition).."\', but should be a function, a type, or a table of allowed values")
			
			local checkFunc
			if type(condition) == "string" then
			
				if not knownStringChecks[condition] then
					error("No type \'"..condition.."\' exists")
				else
					checkFunc = knownStringChecks[condition]
				end
				
			elseif type(condition) == "function" then
				
				checkFunc = condition
				
			elseif type(condition) == "table" then
				
				checkFunc = function(arg)
					
					for _,val in pairs(condition) do
						if val == arg then return true end
					end
					return false,"Argument \'"..argNum.."\' didn't match any of the required values of condition \'"..cIdx.."\'"
				end
			end
			
			local result, err = checkFunc(argument)
			assert(result, err or "Argument \'"..argNum.."\' has failed its condition \'"..cIdx.."\'")
		end
	end
	
	return true
end
	
--[[
	Make a list that associates the values with the keys
	of a table, instead of the other way around
--]]
function utils.Enumify( tbl )
	
	local list = {}
	for k,v in ipairs( tbl ) do
		list[v] = k
	end
	return list
end

function utils.NewEnumIntegers( args )
	
	-- associate each thing in args with an integer
	local o = {}
	for k,v in pairs( args ) do
		o[v] = k
		o[k] = v
	end	
	--[[
		set the metatable so that we can compare
		other variables with this one.  ex
		local a = enum.New{ "thing1", "thing2" }
		local b = 1
		assert( not a.CheckType(b))
		assert(a.thing1 == 1 and a.thing2 == 2)

		local c = enum.New{ "thing3" }
		assert(not c.CheckType(a))
	--]]	
	setmetatable(o, {})
	
	function o.CheckType( x )
		return getmetatable(x) == getmetatable(o)
	end
	
	return o
end

--[[
	return the decimal fraction of a number
--]]
function utils.Frac(number)
	return number - math.floor(number)
end
	
--[[
	Remove trailing digits from a number, 
	ex TruncateDigits( 5.234, 100 ) = 5.23
--]]
function utils.TruncateDigits( num, trunc )
	return math.floor(num*trunc) / trunc	
end
	
--[[utils.AddReadOnly
	add the values in the childTable to parentTable in such
	a way that they are read-only and can't be altered.  This is
	done by not actually storing the childTable values in parentTable,
	but rather in a 'secret' table in parentTable's metatable, so that
	each time the protected variables are accessed, we hit the
	__index and __newindex methods, and can check to see if the variables
	exist in the secret table or not.
	
								+--------------- some existing table
								|			+--- table with values to be nominally added to parentTable
								|			|
								v			v
--]]
function utils.AddReadOnly(parentTable, childTable)
    local mt = getmetatable (parentTable)

    if mt then
        if mt._ROTable == nil then
			-- creating a readonly meta field
            mt._ROTable = {}              
        end
    else    -- no metatable has been set

        mt = {}
        mt._ROTable = {}
        mt.__index = function (t,k)
			
			local value = mt._ROTable[k]
            if value then return value end
			
            error("Unknown key: " .. k)
        end

        mt.__newindex = function (t,k,v)
					
            if mt._ROTable[k] then
                error(k .. " is read-only")
            else
                error("Unknown key: " .. k)
            end
        end

        setmetatable(parentTable, mt)
    end

    if type(childTable) ~= "table" then return end

	-- add all the values from the child table to the read only table
	-- these values won't exist in the parent table, so that all calls to the variables
	-- have to go through __index and __newindex
    for k, v in pairs (childTable) do
        mt._ROTable[k] = v
    end	
	
	return parentTable
end
	
function utils.Save2DTableToFile( fileName, dataTable )

	fun.Map( function(x) assert( type(x) == "table", "dataTable should be 2D") end, dataTable)

	local file,error = io.open(fileName,"w")
	if not file then
		print(error)
	else
			
		local maxSize = 0
		for i,tbl in ipairs(dataTable) do
			maxSize = math.max(maxSize, #tbl)
		end

		-- write the data
		for row=1,maxSize do
				
			-- iterate through the headers
			for i,tbl in ipairs(dataTable) do
											
				-- write the value or a tab for a place-holder
				local val = tbl[row] or "\t"
				file:write(val)

				if i<#dataTable then file:write("\t") else file:write("\n") end
			end
		end

		-- close the file
		file:flush() file:close()
	end
end

--[[
	Call a function with some arguments, in a protected mode, so that
	even if an error is thrown, we don't stop execution.  This can be
	used to run unit tests of functions that would otherwise have error thrown
--]]
function utils.ProtectedRun(func,...)
	
	local ret, err = pcall(func,...)
	if ret then 
		return err
	else 
		return ret, err
	end
end

-- This function returns a deep copy of a given table. 
-- The function below also copies the metatable to the new table if there is one, 
-- so the behaviour of the copied table is the same as the original. 
-- But the 2 tables share the same metatable, 
-- one can avoid this by changing this 'getmetatable(object)' to 'l_copy( getmetatable(object) )'. 
function utils.DeepCopy(object)
    local l_lookupTable = {}
    local function l_copy(object)
        if type(object) ~= "table" then
            return object
        elseif l_lookupTable[object] then
            return l_lookupTable[object]
        end
        local l_newTable = {}
        l_lookupTable[object] = l_newTable
        for index, value in pairs(object) do
            l_newTable[l_copy(index)] = l_copy(value)
        end
        return setmetatable(l_newTable, l_copy(getmetatable(object)))
    end
    return l_copy(object)
end


-- unit tests of the CheckArgs function
local ret,err

ret,err = utils.ProtectedRun(utils.CheckArgs,1)
assert(not ret,"\nWe can't tell that there's only 1 arg")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,2,3)
assert(not ret,"\nWe can't tell that there's odd args")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,2)
assert(not ret,"\nWe can't tell that even args isn't a table")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,{2})
assert(not ret,"\nWe can't tell that check condition isn't a string or function")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,{"number"})
assert(ret,"\nWe can't tell that things are good for one complete argument, check pair")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,{"string"})
assert(not ret,"\nWe can't tell that number is not a string")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,{"number",function(arg) return arg == 1 end})
assert(ret,"\nWe can't tell that 1 is a number and is 1")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,{"number",function(arg) return arg == 2 end})
assert(not ret,"\nWe can't tell that 1 is a number and is not 2")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,{"number",function(arg) return arg == 1 end},true,{"boolean",function(arg) return arg==true end})
assert(ret,"\nWe can't tell that 1 is 1 and true is true")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,{"number",function(arg) return arg == 1 end},true,{"boolean",function(arg) return arg==false end})
assert(not ret,"\nWe can't tell that 1 is 1 and true is not false : ")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,{"number",function(arg) return arg == 2,"Argument is not 1" end})
assert(not ret,"\nWe can't tell that 1 is 1 and true is not false : ")
assert(string.match(err,"Argument is not 1"),"\nWe aren't able to pass function error messages")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,{{2,5,1}})
assert(ret,"\nWe can't tell that 1 matches a required condition")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,{{2,5,num=1}})
assert(ret,"\nWe can't tell that 1 matches a required condition with key")

ret,err = utils.ProtectedRun(utils.CheckArgs,1,{{2,5,6}})
assert(not ret,err)

ret,err = utils.ProtectedRun(utils.CheckArgs,"test",{{2,"test",6}})
assert(ret,err)

local a,b = {},{"item1","item2"}
utils.AddReadOnly(a,b)

ret, err = utils.ProtectedRun(function() a[1] = "itemBad" return true end)
assert(not ret, err)

return utils