--[[
	I was going to make my own set of functional programming functions,
	but then Googled 'map function lua', and found this other stuff already made.
--]]

local fun = {}

 -- Functional Library
 --
 -- @file    functional.lua
 -- @author  Shimomura Ikkeitable.unpack
 -- @date    2005/05/18
 --
 -- @brief    porting several convenience functional utilities form Haskell,Python etc..
 -- map(function, table)
 -- e.g: map(double, {1,2,3})    -> {2,4,6}
 function fun.Map(func, tbl)
     local newtbl = {}
     for i,v in pairs(tbl) do
         newtbl[i] = func(v)
     end
     return newtbl
 end

 -- filter(function, table)
 -- e.g: filter(is_even, {1,2,3,4}) -> {2,4}
 function fun.Filter(func, tbl)
     local newtbl= {}
     for i,v in pairs(tbl) do
         if func(v) then
	     newtbl[i]=v
         end
     end
     return newtbl
 end
 
 -- ifilter(function, table)
 -- e.g: filter(is_even, {1,2,3,4}) -> {2,4}
 function fun.iFilter(func, tbl)
     local newtbl= {}
     for i,v in pairs(tbl) do
         if func(v) then
			table.insert(newtbl,v)
         end
     end
     return newtbl
 end
 
 -- head(table)
 -- e.g: head({1,2,3}) -> 1
 function fun.Head(tbl)
     return tbl[1]
 end

 -- tail(table)
 -- e.g: tail({1,2,3}) -> {2,3}
 --
 -- XXX This is a BAD and ugly implementation.
 -- should return the address to next porinter, like in C (arr+1)
 function fun.Tail(tbl)
     if #tbl < 1 then
         return nil
     else
         local newtbl = {}
         local tblsize = #tbl
         local i = 2
         while (i <= tblsize) do
             table.insert(newtbl, i-1, tbl[i])
             i = i + 1
         end
        return newtbl
     end
 end
 
 -- foldr(function, default_value, table)
 -- e.g: foldr(operator.mul, 1, {1,2,3,4,5}) -> 120
 function fun.Foldr(func, val, tbl)
     for i,v in pairs(tbl) do
         val = func(val, v)
     end
     return val
 end

 -- reduce(function, table)
 -- e.g: reduce(operator.add, {1,2,3,4}) -> 10
 function fun.Reduce(func, tbl)
     return fun.Foldr(func, fun.Head(tbl), fun.Tail(tbl))
 end
 
--[[
	zip together values of multiple tables
	ex. fun.Zip( {1,2,3}, {4,5,6} ) = {{1,4}, {2,5}, {3,6}}
--]] 
function fun.Zip( ... )
	
	local args = table.pack(...)
	args.n = nil
	
	local size = #args[1]		
	fun.Map( function(x) 
				assert( #x == size, 
				"fun.Zip : a table has size " .. #x .. ", all tables should have same size " .. #args[1]) 
			end, args)
	
	
	local z = {}
	for i=1,size do
		z[i] = fun.Map( function(x) return x[i] end, args)
	end

	return z
end

function fun.ZipWithKeys( ... )
	
	local allArgs = table.pack(...)
	assert(#allArgs %2 == 0,"Number of arguments must be even, each table associated with a key")
	
	local tables, keys = {}, {}
	for i,arg in ipairs(allArgs) do
		if i % 2 == 1 then
			assert(type(arg)=="table","argument " .. i .. " must be a table")
			table.insert(tables, arg)
		else
			assert(type(arg)=="string","argument " .. i .. " must be a string")
			table.insert(keys, arg)
		end
	end

	local size = #tables[1]	
	for i,x in ipairs(tables) do
		assert( #x == size, 
				"fun.Zip : the \'"..keys[i] .. "\'has size " .. #x .. ", all tables should have same size " .. #tables[1]) 
	end	
	local zippedArgs = fun.Zip( keys, tables )
	
	local z = {}
	for i=1,size do
		z[i] = {}
		fun.Map( function(x) 
			z[i][ x[1] ] = x[2][i]
		end, zippedArgs)
	end

	return z
end

--[[
  Check if everything in tbl meets a condition
--]]
function fun.All( func, tbl )
    --return fun.Reduce( func, tbl )
	return fun.Foldr( function(a, b) return a and func(b) end, true, tbl )
end
  
--[[
  Check if any item meets the condition
  false is the initial condition, so we have 
  return false or func(tbl[1]) or ... func(tbl[n])
--]]
function fun.Any( func, tbl )
	return fun.Foldr( function(a, b) return a or func(b) end, false, tbl )
end  
  
 -- curry(f,g)
 -- e.g: printf = curry(io.write, string.format)
 --          -> function(...) return io.write(string.format(table.unpack(...))) end
 function fun.Curry(f,g)
     return function (...)
         return f(g(unpack(...)))
     end
 end
 
 -- bind1(func, binding_value_for_1st)
 -- bind2(func, binding_value_for_2nd)
 -- @brief
 --      Binding argument(s) and generate new function.
 -- @see also STL's functional, Boost's Lambda, Combine, Bind.
 -- @examples
 --      local mul5 = bind1(operator.mul, 5) -- mul5(10) is 5 * 10
 --      local sub2 = bind2(operator.sub, 2) -- sub2(5) is 5 -2
 function fun.Bind1(func, val1)
     return function (val2)
         return func(val1, val2)
     end
 end
 
 function fun.Bind2(func, val2) -- bind second argument.
     return function (val1)
         return func(val1, val2)
     end
 end
  
 -- is(checker_function, expected_value)
 -- @brief
 --      check function generator. return the function to return boolean,
 --      if the condition was expected then true, else false.
 -- @example
 --      local is_table = is(type, "table")
 --      local is_even = is(bind2(math.mod, 2), 1)
 --      local is_odd = is(bind2(math.mod, 2), 0)
 fun.Is = function(check, expected)
     return function (...)
         if (check(table.unpack(...)) == expected) then
             return true
         else
             return false
         end
     end
 end
 
 -- operator table.
 -- @see also python's operator module.
 operator = {
     mod = math.mod;
     pow = math.pow;
     add = function(n,m) return n + m end;
     sub = function(n,m) return n - m end;
     mul = function(n,m) return n * m end;
     div = function(n,m) return n / m end;
     gt  = function(n,m) return n > m end;
     lt  = function(n,m) return n < m end;
     eq  = function(n,m) return n == m end;
     le  = function(n,m) return n <= m end;
     ge  = function(n,m) return n >= m end;
     ne  = function(n,m) return n ~= m end;

 }
 -- enumFromTo(from, to)
 -- e.g: enumFromTo(1, 10) -> {1,2,3,4,5,6,7,8,9}
 -- TODO How to lazy evaluate in Lua? (thinking with coroutine)
 fun.EnumFromTo = function (from,to)
     local newtbl = {}
     local step = bind2(operator[(from < to) and "add" or "sub"], 1)
     local val = from
     while val <= to do
         table.insert(newtbl, #newtbl+1, val)
         val = step(val)
     end
     return newtbl
 end
 
 -- make function to take variant arguments, replace of a table.
 -- this does not mean expand the arguments of function took,
 -- it expand the function's spec: function(tbl) -> function(...)
 function fun.Expand_args(func)
     return function(...) return func(...) end
 end
	
return fun