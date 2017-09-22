#Minimum Volume Enclosing Ellipsoid
MVEE := proc(XY, 
			  {tolerance::positive:= 1e-4}, #Convergence Criterion
			  {maxiterations::posint := 100},
			  {output::{identical(data,plot),list(identical(data,plot))} := data},
              {filled::truefalse := false} 
			)

    local alpha, evalues, evectors, i, l_error, ldata, ldataext, M, maxvalindex, n, ncols, nrows, p1, semiaxes, stepsize, U, U1, x, X, y;
    local A, center, l_output; #Output

    if hastype(output, 'list') then
		l_output := output;
	else
		l_output := [output];
	end if;

	kernelopts(opaquemodules=false):

	ldata := Statistics:-PreProcessData(XY, 2, 'copy');

	nrows, ncols := upperbound(ldata);
	ldataext := Matrix([ldata, Vector[column](nrows, ':-fill' = 1)], 'datatype = float');

	if ncols <> 2 then
		error "expected 2 columns of data, got %1", ncols;
	end if;

	l_error := 1;

	U := Vector[column](1..nrows, 'fill' = 1/nrows);

	##Khachiyan Algorithm##
	for n to maxiterations while l_error >= tolerance do

		X := LinearAlgebra:-Transpose(ldataext) . LinearAlgebra:-DiagonalMatrix(U) . ldataext;
		M := LinearAlgebra:-Diagonal(ldataext . LinearAlgebra:-MatrixInverse(X) . LinearAlgebra:-Transpose(ldataext));
		maxvalindex := max[index](map['evalhf', 'inplace'](abs, M));
		stepsize := (M[maxvalindex] - ncols - 1)/((ncols + 1) * (M[maxvalindex] - 1));
		U1 := (1 - stepsize) * U;
		U1[maxvalindex] := U1[maxvalindex] + stepsize;
		l_error := LinearAlgebra:-Norm(LinearAlgebra:-DiagonalMatrix(U1 - U));
		U := U1;

	end do;

	A := (1/ncols) * LinearAlgebra:-MatrixInverse(LinearAlgebra:-Transpose(ldata) . LinearAlgebra:-DiagonalMatrix(U) . ldata - (LinearAlgebra:-Transpose(ldata) . U) . LinearAlgebra:-Transpose((LinearAlgebra:-Transpose(ldata) . U)));
	center := LinearAlgebra:-Transpose(ldata) . U;
	evalues, evectors := LinearAlgebra:-Eigenvectors(A);
	evectors := evectors(.., sort[index](1 /~ (sqrt~(Re~(evalues))), `>`, ':-output' = ':-permutation'));
	semiaxes := sort(1 /~ (sqrt~(Re~(evalues))), `>`);
	alpha := arctan(Re(evectors[2,1]) / Re(evectors[1,1]));

	if l_output = [':-data'] then
		return A, center, evectors, evalues;
	elif has( l_output, ':-plot' ) then
			x := t -> center[1] + semiaxes[1] * cos(t) * cos(alpha) - semiaxes[2] * sin(t) * sin(alpha);
			y := t -> center[2] + semiaxes[1] * cos(t) * sin(alpha) + semiaxes[2] * sin(t) * cos(alpha);
			if filled then
				p1 := plots:-display(subs(CURVES=POLYGONS, plot([x(t), y(t), t = 0..2*Pi], ':-transparency' = 0.95, _rest)));
			else
				p1 := plot([x(t), y(t), t = 0..2*Pi], _rest);
			end if;
		return p1, `if`( has(l_output, ':-data'), op([A, center, evectors, evalues]), NULL );
	end if;

end proc: