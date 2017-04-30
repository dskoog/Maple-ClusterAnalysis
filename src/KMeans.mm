KMeans := proc(
	X::{rtable,DataFrame},
	{clusters::posint := 2}, #Number of clusters
	{initializationmethod::identical(Forgy, Random) := Forgy}, #Initialization Method
	{tolerance::positive:= 1e-4}, #Convergence Criterion
	{maxiterations::posint := 100},
	{output::{list,name} := [':-record'] }
)
	local clusterindex, ctotal, i, j, l_error, l_output, ldata, ncols, notoutput, nrows, outresults, r_col, r_row, resample, results, U;

	if hastype(output, 'list') then
		l_output := output;
	else
		l_output := [output];
	end if;

	kernelopts(opaquemodules=false):

	ldata := Statistics:-PreProcessData(X, 2, 'copy');
	nrows, ncols := upperbound(ldata);

	if clusters > nrows then
		error "number of clusters cannot be greater than number of points";
	end if;

	results := Record('center', 'data', 'initcenter', 'distance', 'cluster', 'numiterations', 'tss', 'tally');

	l_error := tolerance;
	r_row := rand(1..nrows);
	r_col := rand(1..clusters);

	U := Matrix(nrows, clusters, 'datatype' = float);
	results:-data := evaln( X );
	results:-distance := Matrix(nrows, clusters, 'datatype' = float);
	results:-initcenter := Matrix(clusters, ncols, 'datatype' = float);
	clusterindex := Matrix(nrows, clusters, 'fill' = 0);

	if initializationmethod = 'Random' then

		resample := true;

		while resample = true do

			resample := false;
			for i to nrows do  
				clusterindex[i, r_col()] := 1;
			end do;

			for j to clusters do
				ctotal := add(clusterindex[..,j]);
				if ctotal <> 0 then #Guard against no points in a cluster
					results:-initcenter[j] := add(ldata[[ListTools:-SearchAll(1, convert(clusterindex[..,j], 'list'))]][i], i = 1 .. ctotal) / ctotal;
				else
					resample := true;
					break;
				end if;
			end do;

		end do;

	elif initializationmethod = 'Forgy' then

		for j to clusters do
			results:-initcenter[j] := ldata[r_row()]; #Should guard against multiple same rows
		end do;

		for i to nrows do  
			for j to clusters do
				results:-distance[i,j] := sqrt(add(x^2.0, x = ldata[i] - results:-initcenter[j]));
			end do;
			clusterindex[i, min[index](results:-distance[i])] := 1;
		end do;

	end if;

	results:-center := copy(results:-initcenter);	

	for results:-numiterations to maxiterations while l_error >= tolerance do
		ArrayTools:-Copy(results:-distance, U);
		for j to clusters do
			ctotal := add(clusterindex[..,j]);
			if ctotal <> 0 then #Guard against no points in a cluster
				results:-center[j] := add(ldata[[ListTools:-SearchAll(1, convert(clusterindex[..,j], 'list'))]][i], i = 1 .. ctotal) / ctotal;
			end if;
		end do;
		clusterindex := Matrix(nrows, clusters, 'fill' = 0);
		for i to nrows do  
			for j to clusters do
				results:-distance[i,j] := sqrt(add(x^2.0, x = ldata[i] - results:-center[j]));
			end do;
			clusterindex[i, min[index](results:-distance[i])] := 1;
		end do;
        LinearAlgebra:-MatrixAdd(U, results:-distance, 1, -1, 'inplace');
        l_error := max(map['evalhf', 'inplace'](abs, U));
	end do;

	results:-numiterations := results:-numiterations - 1;

	if results:-numiterations >= maxiterations then  
	      error "did not converge: error = %1; try increasing maxiterations or tolerance", l_error;
	end if;
	
	results:-cluster := Array(1..nrows, 'fill' = 0);
	for i to nrows do
		results:-cluster(i) := max[index](clusterindex[i, ..])
	end do;

	results:-tss := add(min(results:-distance[i]), i =  1 .. nrows); #Add the minimum value (distance to closest point) in each row of results:-distance
	results:-tally := Statistics:-Tally(results:-cluster);
	userinfo(1, KMeans, "Iterations: ", results:-numiterations, "Total Sum of Squares: ", results:-tss);

	outresults := Array([]);
	notoutput := Array([]);

	if has(l_output, ':-record') then
		return results;
	else
		for i in l_output do
			if member(i, ['center', 'initcenter', 'distance', 'numiterations', 'cluster', 'tss', 'tally']) then
				ArrayTools:-Append(outresults, results[i]);
			else
				ArrayTools:-Append(notoutput, i);
			end if;
		end do;

		if numelems(notoutput) > 0 then
			error "unrecognized output option: %1", convert(notoutput, 'list');
		end if;

		return seq(outresults);
	end if;

end proc: