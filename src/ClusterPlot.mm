ClusterPlot := proc( X::{rtable,DataFrame,record}, # n x 2 data set
					{center::listlist:=NULL},						
					{cluster::{rtable,list,DataSeries}:=NULL},  
					{style::identical(convexhull, line, ellipse, voronoi, none) := 'convexhull'}
					)
	local datawcluster, diffvalues, i, lcenter, ldata, ldatasplit, metadata, numdf, p1;

	kernelopts(opaquemodules=false):

	if type(X, ':-record') then
		try
			ldata := Statistics:-PreProcessData(X:-data, 2, 'copy');
			datawcluster := Matrix(<ldata | convert(X:-cluster, Vector[column])>);
			lcenter := X:-center;
			metadata := Record( 'numiterations' = X:-numiterations, 'tss' = X:-tss );
		catch:
			error("invalid datatype"):
		end try:
	elif (cluster <> NULL and center <> NULL) 
	  or (hastype(X, {'rtable', 'DataFrame'}) and cluster <> NULL) then
		ldata := Statistics:-PreProcessData(X, 2, 'copy');
		if LinearAlgebra:-ColumnDimension(ldata) <> 2 or LinearAlgebra:-RowDimension(ldata) <> numelems(cluster) then
			error("incompatible data set size"):
		elif center <> NULL and numelems(center) <> numelems(convert(cluster,'set')) then
			error("number of specified center points does not match number of cluster values");
		end if;
		datawcluster := Matrix(<ldata | convert(cluster, Vector[column])>);
		lcenter := center; #TODO - needs to find center for undefined center
		metadata := NULL;
	elif hastype(X, {'rtable', 'DataFrame'}) and cluster = NULL then
		error("option cluster is required for DataFrame and rtable datatypes")
	else
		error("invalid datatype"):
	end if;

	if LinearAlgebra:-ColumnDimension(datawcluster) <> 3 then
		error("incompatible data set size"):
	end if;

	ldata := Statistics:-SplitByColumn(datawcluster, 3);
	numdf := nops(ldata);
	p1 := Array(1..numdf);
	if style = 'convexhull' then
		ldatasplit := Array(1..numdf);
		for i to numdf do
			ldatasplit(i) := ldata[i];
			p1(i) := plottools:-polygon(simplex:-convexhull(convert(ldatasplit[i][..,1..2], 'list', 'nested')), 'color' = plots:-setcolors()[i], 'transparency' = 0.75);
		end do;
	elif style = line then
		ldatasplit := Array(1..numdf);
		for i to numdf do
			ldatasplit(i) := ldata[i];
			p1(i) := CentralPointPlot(ldatasplit[i][..,1..2], lcenter[i,..], 'color' = plots:-setcolors()[i]);
		end do;
	elif style = ellipse then
		ldatasplit := Array(1..numdf);
		for i to numdf do
			ldatasplit(i) := ldata[i];
			p1(i) := MVEE(ldatasplit[i][..,1..2], 'output' = 'plot', 'color' = plots:-setcolors()[i], _rest);
		end do;
	elif style = voronoi then
		if lcenter = NULL then
			error("values for center required with style = voronoi");
		end if;
		p1 := [ComputationalGeometry:-VoronoiDiagram(lcenter, 'colorregions' = 'false', _rest)];
	else
		p1 := [NULL];
	end if;

	if type(X, ':-record') and cluster <> NULL then
		diffvalues := convert(ConfusionMatrix(X:-cluster, cluster, difference), 'list');

		return plots:-display(
				plots:-pointplot(datawcluster[[ListTools:-SearchAll(true,diffvalues)],1..2], 'colorscheme' = ["valuesplit", cluster[[ListTools:-SearchAll(true,diffvalues)]]], 'symbol' = 'solidcircle', 'symbolsize' = 15),
				plots:-pointplot(datawcluster[[ListTools:-SearchAll(false,diffvalues)],1..2], 'colorscheme' = ["valuesplit", cluster[[ListTools:-SearchAll(false,diffvalues)]]], 'symbol' = 'diagonalcross', 'symbolsize' = 15),
                `if`(lcenter <> NULL, plots:-pointplot(lcenter[..,1..2], 'colorscheme' = ["valuesplit", [seq(1..numdf)]], 'symbol' = 'diamond', 'symbolsize' = 30), NULL),
                seq(i, i in p1),
 			'size' = [0.75,0.35], 'axes' = 'boxed', `if`(metadata <> NULL, 'caption' = sprintf("Iterations: %d, TSS: %0.4f", metadata:-numiterations, metadata:-tss), NULL), _rest);
	else
		return plots:-display(
				plots:-pointplot(datawcluster[..,1..2], 'colorscheme' = ["valuesplit", datawcluster[..,3]], 'symbol' = 'solidcircle', 'symbolsize' = 15), 
                `if`(lcenter <> NULL, plots:-pointplot(lcenter[..,1..2], 'colorscheme' = ["valuesplit", [seq(1..numdf)]], 'symbol' = 'diamond', 'symbolsize' = 30), NULL),
                seq(i, i in p1),
 			'size' = [0.75,0.35], 'axes' = 'boxed', `if`(metadata <> NULL, 'caption' = sprintf("Iterations: %d, TSS: %0.4f", metadata:-numiterations, metadata:-tss), NULL), _rest);
	end if;
end proc:

#Local for ClusterPlot - Generates a connected PointPlot where all points are connected to the center
CentralPointPlot := proc(XY, x)
    local i, ldata, ncols, nrows, p1;

	kernelopts(opaquemodules=false):

	ldata := Statistics:-PreProcessData(XY, 2, 'copy');
	nrows, ncols := upperbound(XY);

    p1 := Array(1..nrows);
    for i to nrows do
    	p1(i) := plots:-pointplot([convert(XY[i,..],'list'), convert(x, 'list')], 'style' = 'line', 'symbol' = 'solidcircle', 'symbolsize' = 15, 'transparency' = 0.75, _rest);
    end do;
    return plots:-display(seq(i, i in p1));
end proc: