using Gaston

# helper
function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    (repmat(vx, m, 1), repmat(vy, 1, n))
end

# Create system of linear equations
# for given problem
# column-major encoding of matrix
function make_system(A, fun, h, border)

	(m, n) = size(A);
	m = m - 2; # vrstice (visina)
	n = n - 2; # stolci (sirina)
	nm = n * m;

	# main diagonal
	i = [1:nm]
	j = [1:nm]
	s = ones(nm) .* -4
	# 1st up
	i = [i, 1:nm-1]
	j = [j, 2:nm]
	tmp = ones(nm-1);
	idx = find(mod([1:length(tmp)], m) .== 0);
	tmp[idx] = 0;
	s = [s, tmp]
	# 1st down
	i = [i, 2:nm]
	j = [j, 1:nm-1]
	s = [s, tmp]
	# 2nd up
	i = [i, 1:nm-m]
	j = [j, m+1:nm]
	s = [s, ones(nm-m)]
	# 2nd down
	i = [i, m+1:nm]
	j = [j, 1:nm-m]
	s = [s, ones(nm-m)]

	Z = sparse(i, j, s)

	xx = linspace(border[1, 1], border[1, 2], n);
	yy = linspace(border[2, 1], border[2, 2], m);
	(xxx, yyy) = meshgrid(xx, yy);
	b = 2 .* h .* fun(xxx, yyy);
	b = vec(reshape(b, nm, 1));

	b[1:m] = b[1:m] - A[2:m+1, 1]; #levi rob
	b[nm-m+1:nm] = b[nm-m+1:nm] - A[2:m+1, n+2]; #desni rob
	idx_up = find(mod([1:nm], m) .== 1);
	idx_down = find(mod([1:nm], m) .== 0);
	b[idx_up] = b[idx_up] - A[1, 2:n+1]';
	b[idx_down] = b[idx_down] - A[m+2, 2:n+1]';

	(Z, b)

end

#	CG method for solving Ax=b
# 	ARGs:
#		A: matrix
#		b: vector
#	RETURN: (x, err)
#		x: solution vector
#		err: error history
function cg(A, b)
	cg(A, b, zeros(size(b)))
end


#	CG method for solving Ax=b
# 	ARGs:
#		A: matrix
#		b: vector
#		x: initial approximation of solution (can be zero)
#	RETURN: (x, err)
#		x: solution vector
#		err: error history
function cg(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, x::Array{Float64,1})
	r=b-A*x #residual
	p=r #search direction

    rsold=(r'*r)[1] 
    err = zeros(1000)
	seps = eps(Float32)*10e4

    for i=1:1000
        Ap = vec(A*p) # just mul
        alpha = (rsold/(p'*Ap))[1]
        x = x + alpha * p # new solution
        r = r - alpha * Ap # new residual
        rsnew = (r'*r)[1];
        err[i] = rsnew;
        if rsnew < seps
              break;
        end
        p = r + rsnew / rsold * p
        rsold = rsnew;
    end
    (x, nonzeros(err))
end

function sor3(A, b, w)


	(m, n) = size(A)
	PING = copy(A)
	PONG = copy(A)

	for it = 1:1000
	  	#for all black (i,j) grid points
	  	for i = 2:2:m-1
	  		for j = 2:2:n-1
	     		PONG[i,j] = PING[i,j] + w * (PING[i-1,j] + PING[i+1,j] + PING[i,j-1] + PING[i,j+1] + b[i,j] - 4*PING[i,j])/4
	     	end
	  	end

	  	#for all red (i,j) grid points
		for i = 3:2:m-1
	  		for j = 3:2:n-1
	    		PONG[i,j] = PING[i,j] + w * (PONG[i-1,j] + PONG[i+1,j] + PONG[i,j-1] + PONG[i,j+1] + b[i,j] - 4*PING[i,j])/4
	  		end
	  	end
#	  	if norm(PING-PONG) < eps(Float32)*10e4
#	  		break
#	  	end
	  	PING = copy(PONG)

  	end

  	PONG
end

#	SOR method for solving Ax=b
# 	ARGs:
#		A: matrix
#		b: vector
#		w: relaxation factor (0<w<2)
#		x: initial approximation of solution (can be zero)
#	RETURN: (x, err)
#		x: solution vector
#		err: error history
function sor2(A, b, w, x0)
	d = diag(A)
	D = sparse([1:length(d)], [1:length(d)], d)
	AD = A-D
	L = -tril(AD)
	U = -triu(AD)

	DL = D - w*L
	UE = (1-w)*D + w*U
	bE = w*b

	err = zeros(1000)
	seps = eps(Float32)*10e4

	x = copy(x0)

	for it = 1:1000
		rhv = vec(UE*x0 + bE)
		x = DL\rhv
		e = norm(x-x0, Inf)
		err[it] = copy(e)
		if e < seps
			break
		end
		x0 = copy(x)
	end

	(x, nonzeros(err))
end

#	SOR method for solving Ax=b
# 	ARGs:
#		A: matrix
#		b: vector
#		w: relaxation factor (0<w<2)
#		x: initial approximation of solution (can be zero)
#	RETURN: (x, err)
#		x: solution vector
#		err: error history
function sor(A, b, w, x)
	
	(m, n) = size(A)
	xold = copy(x);
	err = zeros(1000);
	seps = eps(Float32) * 10e4

	for k = 1:1000
		for i = 1:n
			sigma = 0
			for j = 1:i-1
				sigma += A[i, j] * x[j]
			end
			for j = i+1:m
				sigma += A[i, j] * x[j]
			end
			x[i] = (1-w) * x[i] + (w/A[i, i]) * (b[i] - sigma);
		end

		e = norm(x - xold, Inf);
		err[k] = e;
		if e < seps
			break;
		end
		xold = copy(x);
	end

	(x, err)
end

# find optimal relaxation factor for matrix A
#	INPUT:
#		A: matrix
#	RETURN: (wy, wx)
#		wy: relaxation factor in y-dimension
#		wx: relaxation factor in x-dimension

function opt_relax_fact(A)
	(sizey, sizex) = size(A);
	wy = 2/(1+pi/sizey)
	wx = 2/(1+pi/sizex)
	(wy, wx)
end

function solve(A, fun)
	border = [0 1; 0 1];
	solve(A, fun, border);
end

function solve(A, fun, border)
	solve(A, fun, border, "cg");
end

function solve(A, fun, border, method)
	(sizey, sizex) = size(A);
	hx = (border[1, 2] - border[1,1]) / sizex;
	hy = (border[2, 2] - border[2,1]) / sizey;
	assert(abs(hx-hy) < eps(hx), "h should be equal in both dimensions!");
	solve(A, fun, border, hx, method);
end

function solve(A, fun, border, h, method)
	(sizey, sizex) = size(A);
	hx = (border[1, 2] - border[1,1]) / sizex;
	hy = (border[2, 2] - border[2,1]) / sizey;
	assert(abs(hx-hy) < eps(hx), "h should be equal in both dimensions!");
	assert(abs(hx-h) < eps(h), "Specified h differs from computed. Wrong borders?!");
	(Z, b) = make_system(A, fun, h, border);

	init = vec(reshape(A[2:sizey-1, 2:sizex-1], length(b), 1))

	if method == "sor"
		(wy, wx) = opt_relax_fact(A)
		#println("uporabljam SOR-fatkor $wy")
		(x, err) = sor(Z, b, max(wx,wy), init);
	elseif method == "sor2"
		(wy, wx) = opt_relax_fact(A)
		#println("uporabljam SOR-fatkor $wy")
		(x, err) = sor2(Z, b, max(wx,wy), init);
	else
		(x, err) = cg(Z, b, init);
	end
	
	X = reshape(x, sizey-2, sizex-2);
	A[2:sizey-1, 2:sizex-1] = X;

	(A, err);
end

function plot_solution(A, border)
	plot_solution(A, border, "Solution")
end

function plot_solution(A, border, title)
	(sizey, sizex) = size(A);
	xx = linspace(border[1, 1], border[1,2], sizex);
	yy = linspace(border[2, 1], border[2,2], sizey);
	figure()
	surf(xx, yy, A, "title", title);
end

function plot_err(err)
	plot_err(err, "Error")
end

function plot_err(err, title)
	figure()
	plot(linspace(0, length(err), length(err)), err, "axis", "semilogy", "title", title);
end

function border1(sizex)
	sizey = 2* sizex;
	border = [0 1*pi; 0 2*pi];
	xx = linspace(border[1, 1], border[1,2], sizex);
	yy = linspace(border[2, 1], border[2,2], sizey);

	A = zeros(sizey, sizex);
 	A[1, :] = sin(xx);
  	A[sizey, :] = sin(xx);
  	A[:, 1] = cos(yy).-1;
  	A[:, sizex] =  cos(yy).-1;

  	(A, border);
end

function border0(sizex, sizey)
	border = [0 1; 0 sizey/sizex];
	A = zeros(sizey, sizex);
	(A, border);
end

function function1(x, y)
	6.*x.*y.*(1.-y)-2.*x.^3;
end

function function2(x, y)
	x.*0;
end

function function3(x, y)
	x.*0 + 1;
end

function demo()

	(A, border) = border0(50, 50)
	(Z2, err2) = solve(A, function1, border, "sor2")
	(A, border) = border0(50, 50)
	(Z1, err1) = solve(A, function1, border, "cg")

	plot_err(err1, "Napaka pri CG #1")
	plot_err(err2, "Napaka pri SOR #1")
	plot_solution(Z1, border, "Restive CG #1")
	plot_solution(Z2, border, "Resitev SOR #1")

	d1 = norm(Z1-Z2)
	println("Razlika je $d1")

	(A, border) = border1(150)
	(Z2, err2) = solve(A, function3, border, "sor2")
	(A, border) = border1(150)
	(Z1, err1) = solve(A, function3, border, "cg")

	d1 = norm(Z1-Z2)
	println("Razlika je $d1")

	plot_err(err1, "Napaka pri CG #2")
	plot_err(err2, "Napaka pri SOR #2")
	plot_solution(Z1, border, "Restive CG #2")
	plot_solution(Z2, border, "Resitev SOR #2")

end

function bench()

	# Test case 0
	println("Warm up. JIT fun!")
	(A, border) = border0(10, 10);
	print("CG: ")
	@time solve(A, function1, border, "cg");
	(A, border) = border0(10, 10);
	print("SOR: ")
	@time solve(A, function1, border, "sor2");

	# Test case 1
	println("Test case #1: 50x50 grid")
	(A, border) = border0(50, 50);
	print("CG: ")
	@time solve(A, function1, border, "cg");
	(A, border) = border0(50, 50);
	print("SOR: ")
	@time solve(A, function1, border, "sor2");

	# Test case 2
	println("Test case #2: 250x50 grid")
	(A, border) = border0(250, 50);
	print("CG: ")
	@time solve(A, function1, border, "cg");
	(A, border) = border0(250, 50);
	print("SOR: ")
	@time solve(A, function1, border, "sor2");

	# Test case 3
	println("Test case #3: 50x250 grid")
	(A, border) = border0(50, 250);
	print("CG: ")
	@time solve(A, function1, border, "cg")
	(A, border) = border0(50, 250);
	print("SOR: ")
	@time solve(A, function1, border, "sor2");

	# Test case 4
	println("Test case #4: 250x250 grid")
	(A, border) = border0(250, 250);
	print("CG: ")
	@time solve(A, function1, border, "cg");
	(A, border) = border0(250, 250);
	print("SOR: ")
	@time solve(A, function1, border, "sor2");
	
	# Test case 5
	println("Test case #5: 500x1000 grid")
	(A, border) = border0(500, 1000);
	print("CG: ")
	@time solve(A, function1, border, "cg");
	(A, border) = border0(500, 1000);
	print("SOR: ")
	@time solve(A, function1, border, "sor2");

	# Test case 5
	println("Test case #5: 1000x1000 grid")
	(A, border) = border0(1000, 1000);
	print("CG: ")
	@time solve(A, function1, border, "cg");
	(A, border) = border0(1000, 1000);
	print("SOR: ")
	@time solve(A, function1, border, "sor2");

	()
end