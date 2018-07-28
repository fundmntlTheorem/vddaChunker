--ut_Math
require "ut_LUSolve"

--[[
	The 1D line search is an important tool for optimization.  We use it to search
	along a direction for the minimum of an objective function.  In our signal processing
	application, the objective function is the sum of squared differences between data and model.
	But this function allows to pass in any object that has the method Calc_F that returns a scalar.
	So what is the Golden Ratio?  This is the special number that the Greeks found exists in many
	places in nature, and is pleasing to the eye when things are made according to it.  The ratio holds
	for two numbers a,b, a < b, when a / b = b / (a+b).  It is useful in line searches for the following reason:

	Let's say we have two starting points a0 and b0 that bracket a minimum.  We can subdivide the line that
	connects a0 and b0 by choosing two new points a1 and b1.  if f(a1) < f(b1), then we choose the new endpoints
	to be a0 and b1, knowing that the minimum lies between these two.  If rho is the golden ratio, and we originally 
	chose a1 and b1 to be a1 = a0 + rho*(b0-a0) and b1 = a0 + (1-rho)*(b0-a0), then a1 will be perfectly placed in the
	new segment a0 a1 b1, such that we only have to pick one more spot between a1 and b1 to evaluate.  Thus we only
	have to evaluate F once on each iteration after the first one.

	The other nice thing we do here is convert a multidimension problem into a 1D problem, by way of the objective function
	value.  Thus we are finding the minimizer of the function Calc_F(Z(alpha)), where Z = x0 + dk*alpha, Z,x0,dk in R^n, alpha in R^1.
	x0 = a0, and x0 + dk = x2 = b0. So alpha is	the scalar minimizer we are finding, that gives us the final vector minimizer xMin = x0 + alpha*dk.
	
	FofX is the function that returns a scalar value to be minimized, given parameters x
	x0 and x2 are Vectors with the parameters at the start and stop of the bounds
	note that they can also be scalars because the Vector arithmetic overloads the operators
	graph needs to know how to plot x and y, so you can pass him in for testing things
--]]
function math.GoldenSectionSearch(args)
	-- in the optimization they use PHI for f(x) sometimes, so we'll be cool too
	local PHI,x0,x2 = args.FofX,args.x0,args.x2
	local GOLDEN_TOLERANCE,MAX_GOLDEN_ITERATIONS = args.tol or 0.01,args.iter or 100
	local graph = args.graph
	
	local rho = (3-math.sqrt(5.0))/2 -- golden ratio
	local oneMinusRho = 1-rho
	-- the four scalar points that define our search
	local a0,a1,b0,b1
	-- and the function value points
	local f_a1,f_b1
	-- dk is the direction vector for the search, starting at x0, it goes to x2 
	local dk = x2 - x0

	-- the initial bounds of alpha
	a0,b0 = 0,1
	
	local function Z(alpha)
		return x0 + dk * alpha
	end
	
	-- take the first step
	a1 = rho*b0;
	b1 = oneMinusRho*b0
		f_a1 = PHI(Z(a1))
		f_b1 = PHI(Z(b1))
	
	-- step until the tolerance is reached or we exceed the max iterations
	local k=0;
	while math.abs(b0-a0) > GOLDEN_TOLERANCE and k < MAX_GOLDEN_ITERATIONS do
		if f_a1 < f_b1 then -- shift the search region to the LEFT
			b0,b1 = b1,a1
			f_b1 = f_a1
			-- evaluate the new point
			a1 = a0 + rho*(b0-a0)
				f_a1 = PHI(Z(a1))
				
			if graph then graph(Z(a1),f_a1) end
		else -- shift the search region to the RIGHT
			a0,a1 = a1,b1
				f_a1 = f_b1;
			-- evaluate the new point
			b1 = a0 + oneMinusRho*(b0-a0)
				f_b1 = PHI(Z(b1))
				
			if graph then graph(Z(b1),f_b1) end
		end
		k = k + 1
	end
	local xMin
	-- take the minimum to be the lower point
	if f_a1 < f_b1 then
		xMin = Z(a1)
	else
		xMin = Z(b1)
	end
	return xMin
end

function math.PolyFit(xValues,yValues,polyOrder,weights,firstX,lastX)

	firstX = firstX or 1;
	lastX = lastX or #xValues;
--	use DirectPolyFit which takes a list of exponents
--	first construct the exponent list (starting from 0)

	if not (#xValues == #yValues) then
		print ("number of xValues must equal number of yValues");
		return;
	end
	
	if not (#xValues > polyOrder) then
		print ("number of data points must exceed polynomial order");
		return;
	end
	
	local expList = {};
	for i = 0, polyOrder do
		expList[i+1] = i;
	end
	local coeffList = math.DirectPolyFit(xValues,yValues,expList,weights,firstX,lastX);
--	now reindex the array, map the first element in the result to index 0
	local result = {};
	for i = 0, polyOrder do
		result[i] = coeffList[i+1];
	end
	return result;
end


--[[
	The 1D line search is an important tool for optimization.  We use it to search
	along a direction for the minimum of an objective function.  In our signal processing
	application, the objective function is the sum of squared differences between data and model.
	But this function allows to pass in any object that has the method Calc_F that returns a scalar.
	So what is the Golden Ratio?  This is the special number that the Greeks found exists in many
	places in nature, and is pleasing to the eye when things are made according to it.  The ratio holds
	for two numbers a,b, a < b, when a / b = b / (a+b).  It is useful in line searches for the following reason:

	Let's say we have two starting points a0 and b0 that bracket a minimum.  We can subdivide the line that
	connects a0 and b0 by choosing two new points a1 and b1.  if f(a1) < f(b1), then we choose the new endpoints
	to be a0 and b1, knowing that the minimum lies between these two.  If rho is the golden ratio, and we originally 
	chose a1 and b1 to be a1 = a0 + rho*(b0-a0) and b1 = a0 + (1-rho)*(b0-a0), then a1 will be perfectly placed in the
	new segment a0 a1 b1, such that we only have to pick one more spot between a1 and b1 to evaluate.  Thus we only
	have to evaluate F once on each iteration after the first one.

	The other nice thing we do here is convert a multidimension problem into a 1D problem, by way of the objective function
	value.  Thus we are finding the minimizer of the function Calc_F(Z(alpha)), where Z = x0 + dk*alpha, Z,x0,dk in R^n, alpha in R^1.
	x0 = a0, and x0 + dk = x2 = b0. So alpha is	the scalar minimizer we are finding, that gives us the final vector minimizer xMin = x0 + alpha*dk.
	
	FofX is the function that returns a scalar value to be minimized, given parameters x
	x0 and x2 are Vectors with the parameters at the start and stop of the bounds
	note that they can also be scalars because the Vector arithmetic overloads the operators
	graph needs to know how to plot x and y, so you can pass him in for testing things
--]]
function math.GoldenSectionSearch(l_args)
	-- in the optimization they use PHI for f(x) sometimes, so we'll be cool too
	local PHI,x0,x2 = l_args.FofX,l_args.x0,l_args.x2
	local GOLDEN_TOLERANCE,MAX_GOLDEN_ITERATIONS = l_args.tol or 0.01,l_args.iter or 100
	local l_graph = l_args.graph
	
	local rho = (3-math.sqrt(5.0))/2 -- golden ratio
	local oneMinusRho = 1-rho
	-- the four scalar points that define our search
	local a0,a1,b0,b1
	-- and the function value points
	local f_a1,f_b1
	-- dk is the direction vector for the search, starting at x0, it goes to x2 
	local dk = x2 - x0

	-- the initial bounds of alpha
	a0,b0 = 0,1
	
	local function Z(alpha)
		return x0 + dk * alpha
	end
	
	-- take the first step
	a1 = rho*b0;
	b1 = oneMinusRho*b0
		f_a1 = PHI(Z(a1))
		f_b1 = PHI(Z(b1))
	
	-- step until the tolerance is reached or we exceed the max iterations
	local k=0;
	while math.abs(b0-a0) > GOLDEN_TOLERANCE and k < MAX_GOLDEN_ITERATIONS do
		if f_a1 < f_b1 then -- shift the search region to the LEFT
			b0,b1 = b1,a1
			f_b1 = f_a1
			-- evaluate the new point
			a1 = a0 + rho*(b0-a0)
				f_a1 = PHI(Z(a1))
				
			if l_graph then l_graph(Z(a1),f_a1) end
		else -- shift the search region to the RIGHT
			a0,a1 = a1,b1
				f_a1 = f_b1;
			-- evaluate the new point
			b1 = a0 + oneMinusRho*(b0-a0)
				f_b1 = PHI(Z(b1))
				
			if l_graph then l_graph(Z(b1),f_b1) end
		end
		k = k + 1
	end
	local xMin
	-- take the minimum to be the lower point
	if f_a1 < f_b1 then
		xMin = Z(a1)
	else
		xMin = Z(b1)
	end
	return xMin
end

function math.round(num)
	local frac = num % 1
	if frac > 0.5 then 
		num = math.ceil(num) 
	else
		num = math.floor(num)
	end
	return num
end

--[[
	
	Grub's statistic is G = max[i=1...N](abs(yi - yavg) / stddev)
	so find the biggest one, then test if it is past a threshold value.
	The threshold value requires student-t values, and I don't want to 
	do that right now, maybe later.
	Use a recursive call to get rid of outliers until there are no more,
	or list has size 1.  So only one outlier is removed per call maximum,
	and we are assuming a normal distribution
--]]
function math.GrubsOutlierTest(l_args)
	local l_array = l_args.array
	local l_thresh= l_args.thresh or 2
	if #l_array <=1 then
		return l_array
	end
	
	local l_std,l_avg = GetArrayStdDev(l_array)
	local l_maxG,l_maxIndex = 0,0
	for i,val in ipairs(l_array) do
		local l_g = math.abs(val - l_avg) / l_std
		print(val .. " : " .. l_g)
		if l_g > l_maxG then 
			l_maxG = l_g
			l_maxIndex = i
		end
	end
	-- now test against the threshold
	if l_maxG > l_thresh then
		-- remove the outlier
		table.remove(l_array,l_maxIndex)
		return math.GrubsOutlierTest{array=l_array,thresh=l_thresh}
	else
		-- no outliers
		return l_array
	end
end

--[[
	see how much the mean moves as you remove each point, probably this is
	kind of slow.
--]]
function math.PiercesOutlierTest(l_args)
	local l_array = l_args.array
	-- use this if there are other arrays that have data that
	-- depend on the one being filtered, for example if one was looking
	-- at errors in one array, and associated with that was mass values
	local l_associatedArrays = l_args.associatedArrays
	local l_thresh= l_args.thresh or 3
	local l_originalMean = GetArrayAverage(l_array)
	local l_removers = {}
	-- two loops, inner loop finds mean with one point removed
	local function l_newAverage(index)
		local l_tempArray = {}
		for j,val in ipairs(l_array) do
			if j~=index then
				table.insert(l_tempArray,val)
			end
		end	
		return GetArrayAverage(l_tempArray)
	end
	
	-- function to remove the values at the list of indecies
	local function l_removeValues(l_tempArray)
		for i=#l_removers,1,-1 do	
			table.remove(l_tempArray,l_removers[i])
		end	
	end	

	for i=1,#l_array do
		-- compute the mean with one point removed
		local l_tempMean = l_newAverage(i)
		if math.abs(l_tempMean - l_originalMean)/l_tempMean > l_thresh then
			table.insert(l_removers,i)
		end
	end

	-- now remove the bad ones
	l_removeValues(l_array)
	-- and also remove from the associated arrays
	if l_associatedArrays then
		for _,arr in ipairs(l_associatedArrays) do
			l_removeValues(arr)
		end
	end
	
	return l_array
end

--[[
	Check how much the line coefficients move as you remove each point.  if
	they move too much, remove the point.  remove only up to one point
--]]
function math.FilterLinePoints(l_args)
	local l_xData,l_yData = l_args.xData,l_args.yData
	local l_xyPair = Vector.NewXY(l_xData,l_yData)
	local l_thresh = l_args.thresh or 6
	-- get the original coeffient fit
	local l_originalFit = l_xyPair:PolyFit().coeffs
	if not l_originalFit then return false end
	local function l_getNewFit(idx)
		local l_newXYPair = Vector.NewXY()
		-- build arrays without this point
		for i,_ in ipairs(l_xyPair.x) do
			if i~=idx then
				l_newXYPair.x:Push(l_xyPair.x[i])
				l_newXYPair.y:Push(l_xyPair.y[i])
			end
		end
		local l_newFit = l_newXYPair:PolyFit().coeffs
		if not l_newFit then return false end
		local l_slopeDiff = math.abs((l_newFit[1] - l_originalFit[1]) / l_originalFit[1])
		local l_intDiff   = math.abs((l_newFit[0] - l_originalFit[0]) / l_originalFit[0])
		return l_slopeDiff > l_thresh or l_intDiff > l_thresh,math.sqrt(l_slopeDiff^2 + l_intDiff^2)
	end
	
	-- go through and get the fit for with the point removed
	local l_maxIndex
	local l_maxError = 0
	for i,_ in ipairs(l_xData) do
		local l_bad,l_magnitude = l_getNewFit(i)
		if l_bad and l_magnitude > l_maxError then
			l_maxIndex = i
			l_maxError = l_magnitude
		end		
	end
	
	-- check if there was a bad one
	if l_maxIndex then
		print("   Outlier Detected : removing point x " .. l_xData[l_maxIndex] .. " y " .. l_yData[l_maxIndex])
		table.remove(l_xyPair.x,l_maxIndex)
		table.remove(l_xyPair.y,l_maxIndex)	
		l_originalFit = l_xyPair:PolyFit().coeffs
	end
	
	return l_originalFit
end

function math.Median(values)
	table.sort(values);
	local midindex = math.floor(((#values)+1)/2);
	if (#values%2==0) then return (values[midindex] + values[midindex+1])/2;
	else return values[midindex];
	end
end

--auxiliary functions: LineThroughTwoPoints, MassToQuadDAC, MassToTrapDAC, MassToSRIG
--MassToQuadDAC converts a mass value into a MASS_DAC value.
--Four constants (y2Mass, y6Mass, y2QuadDAC, y6QuadDAC) are provided above
--and used to define a linear relationship between mass and the MASS_DAC

function math.LineThroughTwoPoints(x1,y1,x2,y2)
	local dx = x2 - x1;
	local slope = (y2-y1)/dx;
	local intercept = (x2*y1 - x1*y2)/dx;
	return slope,intercept
end

--[[
	Gauss-Jordan Elimination Method for solving matrix system Ax=b
	return matrix beta now containing x, and also matrix alpha, which now
	contains the inverse of the input matrix alpha
--]]
function math.SolveGaussJordan(alphaMatrix,betaMatrix)
	local indexSwappedColumns 	= {}	-- book keeping array for what was the pivot column
	local indexSwappedRows		= {}	-- book keeping array for what was the pivot row
	local indexPivots			= {}	-- flags for whether pivot point was found
	local numberAlphaRows		= #alphaMatrix
	local numberBetaCols		= 1 or #betaMatrix[1]
	local irow,icol
	-- initialize flags
	for i=1,numberAlphaRows do
		indexPivots[i] = 0
	end
	-- main loop over columns to be reduced
	for i=1,numberAlphaRows do

		local big = 0
		for j=1,numberAlphaRows do 		-- searching for pivot point
			if indexPivots[j] ~= 1 then
				for k=1,numberAlphaRows do
					if indexPivots[k] == 0 then
						local datum = math.abs(alphaMatrix[j][k])
						-- searching for the biggest value
						if datum >=big then
							big = datum
							irow = j
							icol = k
						end
					end
				end
			end
		end
		if not icol then return "Error: Singlular Matrix" end		-- this is a failure
		indexPivots[icol] = indexPivots[icol] + 1	-- set this columns flag
		-- now we have pivot element, interchange rows if needed to put the pivot
		-- element on the diagonal
		--print("pivot at : " .. irow .. ", " .. icol .. " value " .. alphaMatrix[irow][icol] )
		if irow ~= icol then
			for index=1,numberAlphaRows do
				local temp = alphaMatrix[irow][index]
				alphaMatrix[irow][index] = alphaMatrix[icol][index]
				alphaMatrix[icol][index] = temp
			end
			if numberBetaCols > 1 then
				for index=1,numberBetaCols do
					local temp = betaMatrix[irow][index]
					betaMatrix[irow][index] = betaMatrix[icol][index]
					betaMatrix[icol][index] = temp
				end
			else
				local temp = betaMatrix[irow]
				betaMatrix[irow] = betaMatrix[icol]
				betaMatrix[icol] = temp
			end
		end
		-- save pivot row,col so we can unscramble alphaMatrix later
		indexSwappedRows[i] = irow
		indexSwappedColumns[i] = icol
		-- now divide the pivot row by the pivot element at irow, icol
		if alphaMatrix[icol][icol] == 0 then return "Error: Singlular Matrix" end
		-- scaling factor
		local pivinv = 1 / alphaMatrix[icol][icol]
		alphaMatrix[icol][icol] = 1.0
		for index=1,numberAlphaRows do
			alphaMatrix[icol][index] = alphaMatrix[icol][index]*pivinv
		end
		if numberBetaCols > 1 then
			for index=1,numberBetaCols do
				betaMatrix[icol][index] = betaMatrix[icol][index]*pivinv
			end
		else
			betaMatrix[icol] = betaMatrix[icol]*pivinv
		end
		-- now reduce the rows
		for redRow=1,numberAlphaRows do
			if redRow ~= icol then
				local temp = alphaMatrix[redRow][icol]
				alphaMatrix[redRow][icol] = 0.0
				for redCol=1,numberAlphaRows do
					alphaMatrix[redRow][redCol] = alphaMatrix[redRow][redCol] - alphaMatrix[icol][redCol]*temp
				end
				if numberBetaCols > 1 then
					for redCol=1,numberBetaCols do
						betaMatrix[redRow][redCol] = betaMatrix[redRow][redCol] - betaMatrix[icol][redCol]*temp
					end
				else
					betaMatrix[redRow] = betaMatrix[redRow] - betaMatrix[icol]*temp
				end
			end
		end
	end -- end of main loop over the columns

	-- now unscramble the solution
	for i=numberAlphaRows,1,-1 do
		if indexSwappedRows[i] ~= indexSwappedColumns[i] then
			for k=1,numberAlphaRows do
				local temp = alphaMatrix[k][indexSwappedRows[i]]
				alphaMatrix[k][indexSwappedRows[i]] = alphaMatrix[k][indexSwappedColumns[i]]
				alphaMatrix[k][indexSwappedColumns[i]] = temp
			end
		end
	end
	-- return the solution and inverse of alphaMatrix too
	return betaMatrix,alphaMatrix
end

--[[	
	example function to show how to
	make a test set of matrices A*x=b and solve for x
--]]
function math.TestingGaussJordan()

	local alpha = {}
	local beta  = {}

	alpha[1] = {} alpha[2] = {} alpha[3] = {}
	alpha[1][1] = 1
	alpha[1][2] = -1
	alpha[1][3] = 1
	alpha[2][1] = 2
	alpha[2][2] = -3
	alpha[2][3] = 4
	alpha[3][1] = -2
	alpha[3][2] = -1
	alpha[3][3] = 1

	beta[1] = 0
	beta[2] = -2
	beta[3] = 7

	local solution = math.SolveGaussJordan(alpha,beta)

	for index,soln in ipairs(beta) do
		print(index .. ", " .. soln)
	end
end

function math.PlotPolyFit2D(data,coeffs)
	GRAPH.Init(0,0,0,0);
	local lastXValue = nil;
	local minYValue,maxYValue;
	local traceIndex = 0;
--	loop through data points
	for index,dataPoint in ipairs(data) do
		local xValue = dataPoint[1];
		local yValue = dataPoint[2];
		local zValue = dataPoint[3];
		if (xValue~=lastXValue) then
			if (traceIndex>0) then
				GRAPH.Legend(traceIndex,lastXValue);
				traceIndex = traceIndex + 1;
				PlotPolySlice(lastXValue,coeffs,minYValue,maxYValue,1000,traceIndex);
			end
			traceIndex = traceIndex + 1;
			minYValue = yValue;
			maxYValue = yValue;
			lastXValue = xValue;
		else
			minYValue = math.min(minYValue,yValue);
			maxYValue = math.max(maxYValue,yValue);
		end
		GRAPH.Plot(yValue,zValue,traceIndex,traceIndex);
	end
	GRAPH.Legend(traceIndex,lastXValue);
	traceIndex = traceIndex + 1;
	PlotPolySlice(lastXValue,coeffs,minYValue,maxYValue,1000,traceIndex);
end
		
function PlotPolySlice(x,coeffs,minY,maxY,numberOfSteps,traceNum)
	GRAPH.Legend(traceNum,x);
	local yStep = (maxY-minY)/numberOfSteps;
	for y = minY, maxY, yStep do
		local z = ComputePoly2D(x,y,coeffs);
		GRAPH.Plot(y,z,traceNum,traceNum);
	end
end

function ComputePoly2D(x,y,coeffs)
	local powersOfX = {};
	local maxPowerOfX = #coeffs;
	powersOfX[0] = 1;
	powersOfX[1] = x;
	local maxPowerOfY = math.max(#coeffs[0],#coeffs[1]);
	for xPower = 2, maxPowerOfX do
		powersOfX[xPower] = powersOfX[xPower-1] * x;
		maxPowerOfY = math.max(#coeffs[xPower],maxPowerOfY);
	end

	local powersOfY = {};
	powersOfY[0] = 1;
	powersOfY[1] = y;
	for yPower = 2, maxPowerOfY do
		powersOfY[yPower] = powersOfY[yPower-1] * y;
	end

	local sum = 0;
	for xPower = 0, maxPowerOfX do
		local xFactor = powersOfX[xPower];
		for yPower = 0, #coeffs[xPower] do
			local yFactor = powersOfY[yPower];
			local term = coeffs[xPower][yPower] * xFactor * yFactor;
			sum = sum + term;
		end
	end
	return sum;
end

function math.QuadraticEquationSolver(coeffTable)
	local a = coeffTable[2];
	local b = coeffTable[1];
	local c = coeffTable[0];
	if (a==0) then return -c/b; end
	local disc = b*b - 4*a*c;
	if (disc<0) then return nil; end
	local term1 = -b/(2*a);
	if (disc==0) then return term1; end
	local term2 = math.sqrt(disc)/(2*a);
	local root1 = term1 + term2;
	local root2 = term1 - term2;
	if (disc>0) then return root1, root2; end
end

function math.RoundDown(value, roundTo)
	return math.floor(value/roundTo)*roundTo
end

function math.RoundUp(value, roundTo)
	return math.ceil(value/roundTo)*roundTo
end

--[[
	Let's say I want a certain precision number.  Say a ppt.
	So I'd want 100.1234 to be truncated to 100.1, and 10.12345
	to be truncated to 10.12.
--]]
function math.TruncateAtPrecision( x, precision )
	if x == 0.0 then
        return x
	end
	local m = math
    local ex = m.floor(m.log10(m.abs(x))) - precision + 1
    local div = m.pow(10, ex)
    return m.floor(x / div + 0.5) * div
end

function math.InnerProduct(a,b)
	local sum = 0;
	for i,av in pairs(a) do
		local bv = b[i];
		if bv==nil then return nil; end
		sum = sum + av*bv;
	end
	return sum;
end

function math.QuadraticForm(a,M,b)
	local N = #a;
	if #b~=N or #M~=N then return nil; end
	local sum = 0;
	for i = 1, N do
		local ai = a[i];
		local Mi = M[i];
		for j = 1,N do
			sum = sum + ai*Mi[j]*b[j];
		end
	end
	return sum;
end

function math.ZeroVector(N)
	local result = {}
	for i=1, N do
		result[i] = 0
	end
	return result
end

function math.ZeroMatrix(N,M)
	local result = {}
	for i=1, N do
		result[i] = math.ZeroVector(M)
	end
	return result
end

function math.CoordinateAxis(n,N)
	local result = math.ZeroVector(N)
	result[n] = 1
	return result
end

function math.GramSchmidt(v1)
	local basis = {};
	basis[1] = math.NormalizeVector(v1);
	local N = #v1;
	for n = 2, N do
		local v = math.CoordinateAxis(n,N);
		for m = 1, n-1 do
			local projection = math.ScaleVector(math.InnerProduct(basis[m],v),basis[m]);
			v = math.SubtractVector(v,projection);
		end
		basis[n] = math.NormalizeVector(v);
	end
	return basis;
end

-- project a onto b
function ProjectVector(a,b)
	return math.ScaleVector(math.InnerProduct(a,b),math.InnerProduct(b,b));
end

function math.AddVector(a,b)
	local N = #a;
	if #b~=N then return nil; end
	local sum = {};
	for i = 1, N do
		sum[i] = a[i] + b[i];
	end
	return sum;
end

function math.PrintVector(v)
	for i = 1, #v do
		print(v[i]);
	end
end

function math.ScaleVector(c,v)
	local cv = {};
	for i,vv in pairs(v) do
		cv[i] = c*vv;
	end
	return cv;
end

function math.NormalizeVector(v)
	local length = math.sqrt(math.InnerProduct(v,v));
	return math.ScaleVector(1/length,v);
end

function math.SubtractVector(a,b)
	return math.AddVector(a,math.ScaleVector(-1,b));
end

function math.DirectPolyFit(xValues,yValues,expValues,weights,firstX,lastX)

	firstX = firstX or 1;
	lastX = lastX or #xValues;
	local X,W = {},{}
	
	local N = lastX - firstX + 1;
	for row = 1, N do
		X[row] = {}
		W[row] = {}
		for col = 1,N do 
			W[row][col] = 0
		end
	end
	
	local Y = {};
	for sampleIndex = firstX, lastX do
		local x = xValues[sampleIndex];
		local shiftedIndex = sampleIndex - firstX + 1;
		Y[shiftedIndex] = yValues[sampleIndex];
		W[shiftedIndex][shiftedIndex] = weights[sampleIndex]
		for expIndex, exponent in ipairs(expValues) do
			X[shiftedIndex][expIndex] = x ^ exponent;
		end		
	end
	
	local XtW = math.MatrixMult(math.Transpose(X), W)
	return LUSolve(math.MatrixMult(XtW,X),math.MatrixVectorMult(XtW,Y));
end

function math.Transpose(a)

	local aTranspose = {};
	for row = 1, #a[1] do
		aTranspose[row] = {};
	end
	
	for rowIndex, rowVector in ipairs(a) do
		for colIndex, entry in ipairs(rowVector) do
			aTranspose[colIndex][rowIndex] = entry;
		end
	end
	
	return aTranspose;
end

function math.PrintMatrix(a)
	for rowIndex,rowVector in ipairs(a) do
		print ("row: ",rowIndex);
		math.PrintVector(rowVector);
	end
end

function math.MatrixMult(a,b)
	local result = {};
	local bT = math.Transpose(b);
	for rowIndex,aRowVector in ipairs(a) do
		result[rowIndex] = {};
		for colIndex,bColVector in ipairs(bT) do
			result[rowIndex][colIndex] = math.InnerProduct(aRowVector,bColVector);
		end
	end
	return result;
end

function math.MatrixVectorMult(a,b)
	local result = {};
	for rowIndex,aRowVector in ipairs(a) do
		result[rowIndex] = math.InnerProduct(aRowVector,b);
	end
	return result;
end

function math.DiscreteHistogram(data)
    local numberOfSamples = #data;
	if numberOfSamples==0 then return end;
	table.sort(data);
	local histogram = {};
	histogram.N = numberOfSamples;
	histogram.minValue = data[1];
	histogram.x = {};
	histogram.y = {};
	local lastValue = histogram.minValue;
	local binIndex = 1;
	while binIndex <= #data do
		local count = 0;
		while data[binIndex]==lastValue do
			binIndex = binIndex + 1;
			count = count+1;
		end
		histogram.x[#(histogram.x)+1] = lastValue;
		histogram.y[#(histogram.y)+1] = count;
		lastValue = data[binIndex];
	end

--	convert counts to count density
	histogram.dx = {};
	if #(histogram.x) == 1 then
		histogram.y[1] = 1;
		histogram.dx[1] = 0;
		return histogram;
	end
	
	local leftBoundary = (3*histogram.x[1] - histogram.x[2])/2; 
	for histIndex = 1, #histogram.x do
		local rightBoundary;
		if (histIndex<#histogram.x) then
			rightBoundary = (histogram.x[histIndex] + histogram.x[histIndex+1])/2;
		else
			rightBoundary = (3*histogram.x[histIndex] - histogram.x[histIndex-1])/2;
		end
		local binWidth = rightBoundary-leftBoundary;
		histogram.dx[histIndex] = binWidth;
		histogram.y[histIndex] = histogram.y[histIndex]/binWidth;		--convert counts to count density
		leftBoundary = rightBoundary;
	end
	
--	normalize histogram to unit area (i.e. a probability density function)
	for histIndex,yValue in ipairs(histogram.y) do
		histogram.y[histIndex] = yValue/numberOfSamples;
	end
	return histogram;
end
		
function math.Histogram(data,numberOfBins)
	numberOfBins = numberOfBins or 100;
	local numberOfSamples = #data;
	if numberOfSamples==0 then return;
	elseif numberOfSamples<numberOfBins then
		numberOfBins = numberOfSamples;
	end
	table.sort(data);
	local histogram = {};
	histogram.N = numberOfSamples;
	histogram.minValue = data[1];
	histogram.x = {};
	histogram.y = {};
	histogram.dx = {};
	local previousXValue = histogram.minValue;
	local weightPerBin = 1/numberOfBins;											--area under curve will integrate to 1
	for binIndex=1, numberOfBins do
		local sampleIndex = math.floor(binIndex*numberOfSamples/numberOfBins);		--each bin contains same number of samples
		local currentXValue = data[sampleIndex];									--right boundary of current bin
		local xRange = currentXValue - previousXValue;								--bin size
		local density;
		if (xRange>0) then
			density = weightPerBin/xRange;											--# of counts in bin is fixed, density increases as bin size shrinks
		else
			local index = sampleIndex;
			while (data[index]==currentXValue) do index = index + 1; end
			if data[index] then
				currentXValue = data[index];
				xRange = currentXValue - previousXValue;
				density = weightPerBin/xRange;
			else
				density = 0;
			end
		end
		local xBinValue = 0.5*(previousXValue+currentXValue);						--center of bin
		histogram.x[binIndex] = xBinValue;
		histogram.y[binIndex] = density;
		histogram.dx[binIndex] = xRange;
		previousXValue = currentXValue;												--left boundary of next bin
	end

	return histogram;
end

function math.PlotHistogram(histogram,plotIndex)
	if not histogram then return end;
	plotIndex = plotIndex or 1;
	local numberOfBins = #histogram.x;
	if plotIndex==1 then GRAPH.ActiveGraph(0); GRAPH.Init(0,0,0,0); end
	GRAPH.PlotTable(histogram.x,histogram.y,plotIndex,plotIndex)
end

function math.PlotCumulative(histogram,plotIndex,scale)
	plotIndex = plotIndex or 1;
	scale = scale or 1;
	local minValue = histogram.minValue;
	local binSize = histogram.binSize;
	if (plotIndex==1) then GRAPH.Init(minValue,0,0,0) end
	local sum = 0;
	local histIndex = 1;
	repeat
		local binValue = minValue + binSize*(histIndex-0.5);
		local value = histogram[histIndex];
		sum = sum + value;
		GRAPH.Plot(binValue,sum*scale,plotIndex,plotIndex);
		histIndex = histIndex + 1;
	until sum>0.995 or histIndex>10000
end

function math.Mean(data)
	local sum = 0;
	for dataIndex = 1, #data do
		local value = data[dataIndex];
		sum = sum + value;
	end
	return sum/#data;
end

function math.StDev(data,mean)
	mean = mean or math.Mean(data);
	local sum = 0;
	for dataIndex = 1, #data do
		local value = data[dataIndex];
		sum = sum + value*value;
	end
	local variance = sum/#data - mean*mean;
	if (variance<0) then return 0; end	-- just to catch errors in precision
	return math.sqrt(variance);
end

function math.Variance(data,mean) --This is a somewhat naive algorithm and should be updated for any serious use for something which handles floating point error a bit better.
	mean = mean or math.Mean(data);
	local sum = 0;
	for dataIndex = 1, #data do
		local value = data[dataIndex];
		sum = sum + value*value;
	end
	local variance = sum/#data - mean*mean;
	if (variance<0) then return 0; end	-- just to catch errors in precision
	return variance;
end

-- search array for a certain key.  assume array is sorted.
-- if exact value is found, second return argument is true
-- if not, then return false.  this will probably be the case
-- if searching for a real number (instead of integer).  in this
-- case you may then wish to interpolate
function math.BinarySearch(array,key)
	local l_first = 1
	local l_last = #array
	local l_mid
	while l_first <= l_last do
		l_mid = math.floor((l_first+l_last)/2)
		if key > array[l_mid] then
			l_first = l_mid+1
		elseif key < array[l_mid] then
			l_last = l_mid-1
		else
			return l_mid,true
		end
	end
	return l_mid,false
end
	
--[[
	The normal binary search function requires that the array that
	is searched be an array of numbers.  It is also useful to be
	able to search a table of something else, so we can create a 
	closure that searches a table by some other means.  In this case
	that other means is by using the get function, which should
	reference the table of interest.  Ex., where map was a table of
	{ {rt = :}, {rt = :}, ...}
	
	local func = math.MakeABinarySearchFunction( function( key )
		return map[key].rt
	end, #map)

			+--<< returns a function which searches for the item 'key'
			|								+--- a function with signature ( int, array ) -> value
			|								|	
			|								|		
			v								v		
--]]
function math.MakeABinarySearchFunction( get )

	return function(key, array)
	
		local first = 1
		local last = #array
		local mid
		while first <= last do
			mid = math.floor((first+last)/2)
			
			-- here you see we use get( mid ) to access the 
			-- values we are searching for
			if key > get( mid, array ) then
				first = mid+1
			elseif key < get( mid, array ) then
				last = mid-1
			else
				return mid,true
			end
		end
		return mid,false
	end
end	
	
--[[
	multiply by rotation matrix
--]]	
function math.RotateXYPoint(x,y,angleRad)
	local l_cos,l_sin = math.cos(angleRad),math.sin(angleRad)
	return l_cos*x-l_sin*y,l_sin*x+l_cos*y
end	
	
--[[
	Eh, I thought I'd make my own histogram the way I like them.
--]]	
function HistogramP(l_args)
	l_args = l_args or {}
	local l_histogram = {}
	l_histogram.minBin 	= l_args.minBin
	l_histogram.maxBin 	= l_args.maxBin
	l_histogram.numBins	= l_args.numBins or 256
	l_histogram.binSize = (l_histogram.maxBin-l_histogram.minBin)/l_histogram.numBins
	l_histogram.graph  = l_args.graph
	-- whether to clamp values to range and thus add them to end bins
	-- if not defined, then they won't be recorded at all
	l_histogram.clamp 	= l_args.clamp	
	l_histogram.fixedY	= l_args.fixedY
	-- initialize the bins
	l_histogram.x = Vector.New(l_histogram.numBins)
	l_histogram.y = Vector.New(l_histogram.numBins)
	l_histogram.sum = 0
	for bin=1,l_histogram.numBins do
		local l_xValue = l_histogram.minBin + l_histogram.binSize*(bin-1)
		-- plot in middle of bin
		l_histogram.x[bin] = l_xValue + l_histogram.binSize/2
		l_histogram.y[bin] = 0
	end
	
	-- function to update the histgram counts
	-- lots of times I like to plot the histogram every say 100 scans or something
	-- to give me an idea what things are looking like
	function l_histogram:Update(data,start)
		local l_start = start or 1
		for i=l_start,#data do
			local l_binIndex = math.floor((data[i]-self.minBin)/self.binSize)
			if self.clamp==true then
				l_binIndex = math.min(self.numBins,l_binIndex)
				l_binIndex = math.max(1,l_binIndex)
			end
			-- add if within bounds
			if l_binIndex >= 1 and l_binIndex <=self.numBins then
				self.y[l_binIndex] = self.y[l_binIndex]+1
				-- update count so plot can normalize y axis
				self.sum = self.sum + 1
			end
		end
		return self -- allow to call this function on initialization HistogramP{}:Update(data)
	end	
	-- plot the data
	function l_histogram:Plot(l_args)
		l_args = l_args or {}
		local l_reInit = l_args.reInit
		local l_noNormalize = l_args.noNormalize
		local l_norm = self.sum
		if l_noNormalize then
			l_norm = 1
		end
	
		-- clear the graph if requested
		if l_reInit==true then
			-- if a graph object exists, use it
			if self.graph then
				self.graph:Init()
			else -- or just dirty init the graph
				GRAPH.Init(0,0,0,0)
			end
		end
		-- do the plotting
		local l_trace = l_args.trace or 1
		for index,_ in ipairs(self.x) do
			if self.graph then
				self.graph:Plot(self.x[index],self.y[index]/l_norm,l_trace)
			else
				GRAPH.Plot(self.x[index],self.y[index]/l_norm,l_trace)
			end
		end	
		return self
	end
	
	--[[
		A neat histogram that is rotated 90 degrees
	--]]
	function l_histogram:PlotRotated(l_args)
		if #self.y ==0 then
			return 
		end
		l_args = l_args or {}
		local l_trace = l_args.trace
		local l_maximumY = self.fixedY or l_args.fixedY or 1
		local l_xOffset  = l_args.xOffset or 0
		-- we want a consistent height of the graph
		local _,l_normY = FindArrayMax(self.y)
		if l_args.activeGraph then
			GRAPH.ActiveGraph(l_args.activeGraph)
		end
		
		-- plot rotated points
		for i,x in ipairs(self.x) do
			local l_y = self.y[i]/l_normY*l_maximumY
			if l_y > 0 then
				local l_rx,l_ry = math.RotateXYPoint(x,l_y,math.pi/2)
				GRAPH.Plot(l_rx+l_xOffset,l_ry,l_trace)
			end
		end
		return self
	end
	
	-- return an average value
	-- SUM[i](xi * pi) = 1/N*SUM[i](xi *yi)
	function l_histogram:Average()
		local l_sumY = self.sum
		local l_avg = 0
		for index,_ in ipairs(self.x) do
			l_avg = l_avg + (self.x[index]*self.y[index])
		end
		-- finish the average
		return l_avg/l_sumY,l_sumY
	end
	
	-- var(x) = SUM[i](xi^2 * pi) - avg^2 = 1/N*SUM[i](xi^2 *yi) - avg^2
	function l_histogram:StdDev()
		local l_avg,l_sum = self:Average()
		local l_std = 0
		for index,_ in ipairs(self.x) do
			l_std = l_std + self.x[index]*self.x[index]*self.y[index]
		end
		-- complete the variance value
		l_std = l_std / l_sum
		return math.sqrt(l_std - l_avg^2)
	end
	
	function l_histogram:GetCDF()
		local l_y = {}
		local l_val = 0
		for index,datum in ipairs(self.y) do
			l_val = l_val + datum/self.sum
			table.insert(l_y,l_val)
		end
		return l_y
	end	
	
	-- return the xvalue at the pct percentile of the distribution
	-- ex. I want the xvalue at 0.95
	function l_histogram:GetPercentile(pct)
		local l_cdf = self:GetCDF()
		local l_val=self.x[1]
		local l_index = 1
		for index,datum in ipairs(l_cdf) do
			if datum > pct then
				l_val = self.x[index]
				l_index = index
				break
			end
		end	
		return l_val,l_index
	end
	
	-- See if there's data in the histogram.
	function l_histogram:HasData()
		local l_sumY = 0
		for index,value in ipairs(self.y) do 
			l_sumY = l_sumY + value
		end
		local l_result
		if l_sumY == 0 then 
			l_result = false
		else
			l_result = true
		end
		return l_result
	end
	
	-- return this little class
	return l_histogram
end

--	TestDataAnalysis is the main program
--	It creates a fake data set (if requested), fits it to a model function, and plots the data and the fit.
function math.RobFit(fitDataArray,I0,tau,T,dt,plotData)
	--	FakeData
	-- 	returns an array of data points {(tn,kn): n = 1, N}, tn = (n-1)*dt, N = 1 + T/dt
	--	kn is generated by calling ReactionRateConstant
	local function l_RobFit_FakeData(I0,tau,T,dt)
		--	ReactionRateConstant
		--	This is the analytic function for generating fake data.
		--	It is qualitatively like the observed time-dependence of reaction rate constant.
		local function l_Rob_ReactionRateConstant(t,I0,tau)
			return I0*(1-math.exp(-t/tau));
		end

		local dataArray = {};
		for t = 0, T, dt do
			local dataPoint = {};
			dataPoint.t = t;
			dataPoint.k = l_Rob_ReactionRateConstant(t,I0,tau);
			dataArray[#dataArray+1] = dataPoint;
		end
		return dataArray;
	end
	
	--	PlotData
	local function l_Rob_PlotData(dataArray)
		GRAPH.Layout(2);										--use trace 1 for data (will use
		GRAPH.ActiveGraph(0);
		GRAPH.Title("Data");
		GRAPH.Axes("Target","Rate Constant (msec-1)");
		local minTime = dataArray[1].t;
		local maxTime = dataArray[#dataArray].t;
		GRAPH.Init(minTime,0,maxTime,0);
		local plotIndex = 1;
		for _,dataPoint in ipairs(dataArray) do
			local t = dataPoint.t;
			local k = dataPoint.k;
			GRAPH.Plot(t,k,plotIndex,plotIndex);
		end
	end
	
	--	SquaredError
	--	returns sum of squared differences between model and data
	--	The model is hard-coded and has two linear segments, resulting in two sums.
	local function l_Rob_SquaredError(dataArray,t0Index,I0)
		local t0 = dataArray[t0Index].t;
		local errorSum1 = 0;
		for sampleIndex = 1, t0Index do
			local dataPoint = dataArray[sampleIndex];
			local t = dataPoint.t;
			local k_obs = dataPoint.k;
			local k_calc = I0 * t/t0;							-- linear segment of model
			local delta_k = k_obs - k_calc;
			local squaredErrorTerm = delta_k * delta_k;
			errorSum1 = errorSum1 + squaredErrorTerm;
		end
		
		local errorSum2 = 0;
		for sampleIndex = t0Index+1, #dataArray do
			local dataPoint = dataArray[sampleIndex];
			local t = dataPoint.t;
			local k_obs = dataPoint.k;
			local k_calc = I0;									-- flat segment of model
			local delta_k = k_obs - k_calc;
			local squaredErrorTerm = delta_k * delta_k;
			errorSum2 = errorSum2 + squaredErrorTerm;
		end

		local squaredError = errorSum1 + errorSum2;
		return squaredError;
	end
	
	--	FitI0AndT0
	--	Model is a piecewise linear function
	--	line 1: (0,0) to (t0,I0)
	--	line 2: (t0,I0) to (T,I0)
	--	There are two parameters t0,I0.  T is arbitrary (defined by data).
	--	The fit minimizes the sum of squared differences.
	--	t0 is stepped over the data points.  A closed form for I0 is calculated at each t0.
	local function l_Rob_FitI0AndT0(dataArray,plotData)

		--	For a fixed t0, the least-squares value I0 has a closed-form solution.
		--	It is the ratio of two pairs of two sums.
		--	i.e. bestI0 = (N1 + N2) / (D1 + D2)
		--	The two sums correspond to the two linear segments.	
		local function l_Rob_FitI0GivenT0(dataArray,t0Index)
			local t0 = dataArray[t0Index].t;
			local numeratorSum1 = 0;
			local denominatorSum1 = 0;
			for sampleIndex = 1, t0Index do
				local dataPoint = dataArray[sampleIndex];
				local t = dataPoint.t;
				local k_obs = dataPoint.k;
				local k_calc_norm = t/t0;								-- linear segment of model(divided by I0)
				local numeratorTerm = k_obs * k_calc_norm;
				local denominatorTerm = k_calc_norm * k_calc_norm;		
				numeratorSum1 = numeratorSum1 + numeratorTerm;
				denominatorSum1 = denominatorSum1 + denominatorTerm;
			end
			
			local numeratorSum2 = 0;
			local denominatorSum2 = 0;
			for sampleIndex = t0Index+1, #dataArray do
				local dataPoint = dataArray[sampleIndex];
				local t = dataPoint.t;
				local k_obs = dataPoint.k;
				local k_calc_norm = 1;									-- flat segment of model(divided by I0)
				local numeratorTerm = k_obs;
				local denominatorTerm = k_calc_norm;
				numeratorSum2 = numeratorSum2 + numeratorTerm;
				denominatorSum2 = denominatorSum2 + denominatorTerm;
			end

			local numerator = numeratorSum1 + numeratorSum2;			-- add sums for two segments
			local denominator = denominatorSum1 + denominatorSum2;
			local bestI0 = numerator/denominator;
			return bestI0;
		end	
	
	--	Calculate best I0 for each of multiple t0 values, where t0 is set to each observed data point.
	--	Also calculate the squared error for each choice of t0.
		local squaredErrorArray = {};
		local i0Array = {};
		for t0Index = 1, #dataArray do
			local I0 = l_Rob_FitI0GivenT0(dataArray,t0Index);
			local squaredError = l_Rob_SquaredError(dataArray,t0Index,I0);
			--print(t0Index,I0,squaredError);
			squaredErrorArray[t0Index] = squaredError;
			i0Array[t0Index] = I0;
		end

		local errorIndex, errorMax = FindArrayMax(squaredErrorArray,2);
		local errorIndex, errorMin = FindArrayMin(squaredErrorArray,2);
		if plotData then 
			GRAPH.ActiveGraph(1);
			GRAPH.Init(0,0,0,0)
			local plotIndex = 3;
			GRAPH.Axes("Target","squared error");
			for t0Index = 2, #dataArray do
				local t0 = dataArray[t0Index].t;
				GRAPH.Plot(t0,squaredErrorArray[t0Index],plotIndex,plotIndex);
			end
		end
	--	Now select (and return) the model from among these t0 values with the minimum squared error. 
	--	The error value at t0=1 is bogus, so start with second element.
		local startIndex = 2;
		local bestT0Index = FindArrayMin(squaredErrorArray,startIndex);
		local bestI0 = i0Array[bestT0Index];
		local bestT0 = dataArray[bestT0Index].t;
		return bestI0,bestT0,squaredErrorArray;
	end

	--	Plot Model
	--	The model (two linear segments) is superposed on the data plot (trace = 2)
	--	A legend displays the optimal parameter values.
	local function l_Rob_PlotModel(I0,t0,T)
		GRAPH.ActiveGraph(0);
		local plotIndex = 2;
		local legendString = string.format("Model: (k = %6.2f, Target = %6.2f)",I0,t0);
		GRAPH.Legend(plotIndex,legendString);
		GRAPH.Plot(0,0,plotIndex,plotIndex);
		GRAPH.Plot(t0,I0,plotIndex,plotIndex);
		GRAPH.Plot(T,I0,plotIndex,plotIndex);
	end

	local I0 = I0 or 100									--intensity at the plateau (asymptote)
	local tau = tau or 10									--exponential time constant (for fake data)
	local dt = dt or 1										--uniform time step between data points
	local T1 = 100											--Set T to 100 for fake data set. Overwrite later for real data set from routine.
	local fitDataArray = fitDataArray or l_RobFit_FakeData(I0,tau,T1,dt)		--pass routine data or create a fake data set
	local tempArray = {}
	for i,dataPoint in ipairs(fitDataArray) do
		tempArray[i] = dataPoint.t
	end
	local maxIndex, max = FindArrayMax(tempArray)
	T1 = max																		--duration of data stream
	if plotData then l_Rob_PlotData(fitDataArray) end								--plot the data
	local I0_est,t0_est,squaredErrorArray = l_Rob_FitI0AndT0(fitDataArray,plotData) --find best model (I0,t0)
	if plotData then l_Rob_PlotModel(I0_est,t0_est,T1) end							--plot the model (superimpose on data)
	return t0_est, I0_est, T1, squaredErrorArray
end

--[[
	General purpose bubble sort
	unsortedTable is the one we want to sort.  
	otherTables are tables that need to be indexed the same as unsortedTable
	sortFunc is how to sort the table
--]]
function math.Sort(unsortedTable,otherTables,sortFunc)
	otherTables = otherTables or {}
	-- by default we sort from small to large
	sortFunc = sortFunc or  function (x,y)
								return x > y
							end
	for i=1,#unsortedTable do
		for j=1,#unsortedTable-1 do
			if sortFunc(unsortedTable[j],unsortedTable[j+1]) then
				-- swap values with neat Lua method
				unsortedTable[j],unsortedTable[j+1] = unsortedTable[j+1],unsortedTable[j]
				-- swap any other arrays that belong to unsortedTable
				for _,tbl in ipairs(otherTables) do
					tbl[j],tbl[j+1] = tbl[j+1],tbl[j]
				end
			end
		end
	end
end

--[[
	evaluate a polynomial, given x and the coefficients.
	the coefficients are in increasing order, 1 to n
--]]
function math.PolyEval(x,coeff)
	-- check to see if the array is zero based.  the polyfit routine
	-- is zerobased, but most lua tables aren't made that way
	local l_start,l_end = #coeff
	if coeff[0] then
		l_end = 0
	else
		l_end = 1
	end
	-- do the evaluation
	local l_f = coeff[l_start]
	for i=l_start-1,l_end,-1 do
		l_f = l_f*x + coeff[i]
	end
	return l_f		
end

--[[
	Solve a complex polynomial using Newton's method
--]]
function math.NewtonMethod(l_args) 

	-- Import the data to analyze and initialize the output
	local l_coefficients = l_args.coefficients
	local l_initialValue = l_args.initialValue or 0
	local l_iterations = l_args.iterations or 500
	local l_minIterations = l_args.minIterations or 10
	local l_maxIterations = l_args.maxIterations or 20
	local l_minDeltaValue = l_args.minDeltaValue or 0.0001
	
	-- Initialize a few more local variables and functions
	local l_oldValue = l_initialValue
	local l_currentValue = l_initialValue
	local l_deltaValue
	
	local function l_function(xValue,coefficients)
		-- Treat every funtion as a simply x-y function
		local l_xValue = xValue
		local l_coefficients = coefficients
		local l_yValue = 0
		-- Calculate the y-value of the main function
		for i = 1,#coefficients do
			l_yValue = l_yValue + coefficients[i]*l_xValue^(i-1)
		end
		return l_yValue
	end	
	
	local function l_derivativeFunction(xValue,coefficients)
		-- Treat every funtion as a simply x-y function
		local l_xValue = xValue
		local l_coefficients = coefficients
		local l_yValue = 0
		-- Calculate the y-value of the main function
		for i = 1,#coefficients do
			l_yValue = l_yValue + coefficients[i]*(i-1)*l_xValue^(i-2)
		end
		return l_yValue
	end	
	
	-- Enter the main loop of the gradient descent procedure
	for i = 1,l_iterations do
		-- Calculate the yValue of the function and its derivative
		local l_yValue = l_function(l_currentValue,l_coefficients)
		local l_yValueDerivative = l_derivativeFunction(l_currentValue,l_coefficients)
		-- Update the currrent value
		l_currentValue = l_currentValue - l_yValue/l_yValueDerivative
		-- Calculate the change in value and update the old value
		l_deltaValue = math.abs((l_currentValue-l_oldValue)/l_oldValue)
		l_oldValue = l_currentValue	
		if l_minIterations < i then	
			-- If the change in cost was less than the minimum threshold then bail out and don't waste any more time
			if l_deltaValue < l_minDeltaValue then break end
		end
		if l_maxIterations < i then	break end
	end
	
	return l_currentValue
end

--[[
	modified bessel function of the first kind,
	from NR
--]]
function math.BesselMod1(xvalue)
	-- coefficients for evaluating some polynomials
	local l_p = {5.000000000000000e-1,6.090824836578078e-2,
				2.407288574545340e-3,4.622311145544158e-5,5.161743818147913e-7,
				3.712362374847555e-9,1.833983433811517e-11,6.493125133990706e-14,
				1.693074927497696e-16,3.299609473102338e-19,4.813071975603122e-22,
				5.164275442089090e-25,3.846870021788629e-28,1.712948291408736e-31}
	local l_q = {4.665973211630446e-1,1.677754477613006e-3,
					2.583049634689725e-6,2.045930934253556e-9,7.166133240195285e-13}
	local l_pp = {1.286515211317124e-1,1.930915272916783e-1,
					6.965689298161343e-2,7.345978783504595e-3,1.963602129240502e-4}
	local l_qq = {3.309385098860755e-1,4.878218424097628e-1,
					1.663088501568696e-1,1.473541892809522e-2,1.964131438571051e-4,
					-1.034524660214173e-6}
	-- there are two regimes, when x < 15 and x >= 15
	local l_absX = math.abs(xvalue)
	if l_absX < 15 then
		local l_y = xvalue*xvalue
		return xvalue*math.PolyEval(l_y,l_p)/math.PolyEval(225-l_y,l_q)
	else
		local l_z = 1.0-15/l_absX
		local l_ans = math.exp(l_absX)*math.PolyEval(l_z,l_pp)/math.PolyEval(l_z,l_qq)/math.sqrt(l_absX)
		if l_ans <= 0 then l_ans = -l_ans end
		return l_ans
	end
end

--[[
	Round value to the specified digits
	procedure is to multiply by 10^digits place, so
	that the specified digit is at the ones place, then use floor,
	then multiply back
		ex 10.743
		round to hundreths place
		10.743 * 10^2 = 1074.3
		floor(1074.3) = 1074
		1074 * 10^-2 = 10.74
	The convention here is that fractional numbers have 10^-digits,
	columns greater than 1s column are positive sign
--]]
function math.RoundDigits(number,digits)
	return math.floor(number*10^-digits)*10^digits
end

function math.PolynomialModel1D(...)
	local coefficients={...}
	return function (x)
		local retval=0
		for i=1,#coefficients do --This could be done with map and reduce but this implementation is lighter.
			retval=retval+coefficients[i]*x^(i-1)
		end
		return retval
	end
end

--Check for nan (not a number)
function math.isNan(number)
	return number ~= number 
end
