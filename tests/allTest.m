%% Main function to generate tests
function tests = allTest
    tests = functiontests( ...
        {@chapter1test, @chapter2test, @chapter3test, @chapter4test, ...
        @chapter5test, @chapter6test, @chapter8test, @chapter9test, ...
        @chapter10test, @chapter11test, @chapter13test});
end

function chapter1test(testCase)
    assertTrue(testCase, isapprox(horner([1,-3,3,-1], 1.6), 0.6^3));
end

function chapter2test(testCase)
    A = [ 1 2 3 0; -1 1 2 -1; 3 1 2 4; 1 1 1 1 ];
    [L, U] = lufact(A);
    assertTrue(testCase, isapprox(L*U, A));
    assertTrue(testCase, isapprox(tril(L), L));
    assertTrue(testCase, isapprox(triu(U), U));
    assertTrue(testCase, isapprox(diag(L), ones(4, 1)));
    b = [1;10;0;-1] / 5;
    assertTrue(testCase, isapprox(L\b, forwardsub(L, b)));
    assertTrue(testCase, isapprox(U\b, backsub(U, b)));
    [L, U, p] = plufact(A);
    assertTrue(testCase, isapprox(L*U, A(p, :)));
    assertTrue(testCase, isapprox(tril(L), L));
    assertTrue(testCase, isapprox(triu(U), U));
    assertTrue(testCase, isapprox(diag(L), ones(4, 1)));
end

function chapter3test(testCase)
    A = [3 4 5;-1 0 1;4 2 0; 1 1 2; 3 -4 1];
    b = (5:-1:1)';
    assertTrue(testCase, isapprox(lsnormal(A, b), A\b))
    assertTrue(testCase, isapprox(lsqrfact(A, b), A\b))
    [Q, R] = qrfact(A);
    assertTrue(testCase, isapprox(Q*R, A))
    assertTrue(testCase, isapprox(Q*Q', eye(5)))
    assertTrue(testCase, isapprox(triu(R), R))
    assertEqual(testCase, size(R), [5, 3])
end

function chapter4test(testCase)
    for c = [2 4 7.5 11]
        f = @(x) exp(x) - x - c;
        dfdx = @(x) exp(x) - 1;
        x = newton(f, dfdx, 1.0);  r = x(end);
        assertTrue(testCase, isapprox(f(r), 0, 0, 1e-13))
    end
    
    for c = [2 4 7.5 11]
        f = @(x) exp(x) - x - c;
        dfdx = @(x) exp(x) - 1;
        x = secant(f, 1.0, 2.0);  r = x(end);
        assertTrue(testCase, isapprox(f(r), 0, 0, 1e-13))
    end
    
        function [f, J] = nlfun(x)
            f = zeros(3,1);
            f(1) = exp(x(2)-x(1)) - 2;
            f(2) = x(1)*x(2) + x(3);
            f(3) = x(2)*x(3) + x(1)^2 - x(2);
            J = zeros(3,3);
            J(1, :) = [-exp(x(2)-x(1)) exp(x(2)-x(1))  0];
            J(2, :) = [x(2) x(1) 1];
            J(3, :) = [2*x(1) x(3)-1 x(2)];
        end
    
    x = newtonsys(@nlfun, [0,0,0]);
    r = x(:, end);
    assertTrue(testCase, isapprox(nlfun(r), [0;0;0], 0, 1e-12))
    x = newtonsys(@nlfun, [1,2,3]);
    r = x(:, end);
    assertTrue(testCase, isapprox(nlfun(r), [0;0;0], 0, 1e-12))
    x = levenberg(@nlfun, [4,-4,-3]);
    r = x(:, end);
    assertTrue(testCase, isapprox(nlfun(r), [0;0;0], 0, 1e-12))
end


function chapter5test(testCase)
    f = @(t) cos(5*t);
    [Q, t] = intadapt(f,-1,3,1e-8);
    assertTrue(testCase, isapprox(Q, (sin(15)+sin(5))/5, 1e-5))
    [T,~] = trapezoid(f,-1,3,820);
    assertTrue(testCase, isapprox(T, (sin(15)+sin(5))/5, 1e-4))
    
    t = [-2,-0.5,0,1,1.5,3.5,4]/10;
    w = fdweights(t-0.12, 2);
    f = @(x) cos(3*x);
    assertTrue(testCase, isapprox(dot(w,f(t)), -9*cos(0.36), 1e-3))
    
    t = [-2,-0.5,0,1,1.5,3.5,4]/10;
    
    H = hatfun(t, 5);
    assertTrue(testCase, isapprox(H(0.22), (0.22-t(5))/(t(6) - t(5))))
    assertTrue(testCase, isapprox(H(.38), (t(7)-.38)/(t(7) - t(6))))
    assertEqual(testCase, H(0.06), 0)
    assertEqual(testCase, H(t(6)), 1)
    assertEqual(testCase, H(t(7)), 0)
    
    p = plinterp(t, f(t));
    x = [-0.14 -0.01 0 0.03 0.11 0.15 0.25 0.3 0.4];
    y = interp1(t, f(t), x, "linear");
    assertTrue(testCase, isapprox(p(x), y))
    
    S = spinterp(t, exp(t));
    x = [-.17, -0.01, 0.33, .38];
    assertTrue(testCase, isapprox(S(x), exp(x), 1e-5))
    assertTrue(testCase, isapprox(S(t), exp(t), 1e-10))
end

function chapter6test(testCase)
    f = @(t, u, p) u + p*t^2;
    u_ex = exp(1.5) - 2*(-2 + 2*exp(1.5) - 2*1.5 - 1.5^2);
    ivp = ode(ODEFcn=f);
    ivp.InitialValue = 1;
    ivp.Parameters = -2;
    sol = solve(ivp, 0, 1.5);
    [t, u] = eulerivp(ivp, 0, 1.5, 5000);
    assertTrue(testCase, isapprox(u(end), sol.Solution(end), 5e-3))
    
    [t, u] = am2(ivp, 0, 1.5, 4000);
    assertTrue(testCase, isapprox(u(end), sol.Solution(end), 1e-4))
    
    ivp = ode(ODEFcn=@(t, u, p) [t+p-sin(u(2)); u(1)]);
    ivp.InitialValue = [-1; 4];
    ivp.Parameters = -6;
    ivp.InitialTime = 1;
    ivp.RelativeTolerance = 1e-7;
    ivp.AbsoluteTolerance = 1e-7;
    sol = solve(ivp, 1, 2);
    
    [t, u] = eulerivp(ivp, 1, 2, 5000);
    assertTrue(testCase, isapprox(u(:,end), sol.Solution(:,end), 5e-3))
    
    [t, u] = ie2(ivp, 1, 2, 2000);
    assertTrue(testCase, isapprox(u(:,end), sol.Solution(:,end), 1e-4))
    
    [t, u] = rk4(ivp, 1, 2, 1000);
    assertTrue(testCase, isapprox(u(:,end), sol.Solution(:,end), 1e-6))
    
    [t, u] = ab4(ivp, 1, 2, 1000);
    assertTrue(testCase, isapprox(u(:,end), sol.Solution(:,end), 1e-6))
    
    [t, u] = rk23(ivp, 1, 2, 1e-5);
    assertTrue(testCase, isapprox(u(:,end), sol.Solution(:,end), 2e-5))
    
    [t, u] = am2(ivp, 1, 2, 2000);
    assertTrue(testCase, isapprox(u(:,end), sol.Solution(:,end), 1e-4))
end

function chapter8test(testCase)
	V = randn(4, 4);
	D = diag([-2,0.4,-0.1,0.01]);
	A = V*D/V;

	[beta,x] = poweriter(A, 30);
	assertTrue(testCase, isapprox(beta(end), -2, 1e-10))
    d = dot(x,V(:,1)) / (norm(V(:,1)) * norm(x));
	assertTrue(testCase, isapprox(abs(d), 1, 1e-10))

    [beta,x] = inviter(A, 0.37, 15);
	assertTrue(testCase, isapprox(beta(end), 0.4, 1e-10))
    d = dot(x,V(:,2)) / (norm(V(:,2)) * norm(x));
	assertTrue(testCase, isapprox(abs(d), 1, 1e-10))

	[Q, H] = arnoldi(A, ones(4,1), 4);
    assertTrue(testCase, isapprox(A*Q(:,1:4), Q*H))
	
    b = ones(4,1);
	[x, res] = arngmres(A, b, 3);
    assertTrue(testCase, isapprox(norm(b - A*x), res(end)))
	[x, res] = arngmres(A, b, 4);
	assertTrue(testCase, isapprox(b, A*x, 1e-10))
end

function chapter9test(testCase)
	f = @(x) exp(sin(x) + x.^2);
	t = -cos((0:40)'*pi/40);
	p = polyinterp(t, f(t));
    assertTrue(testCase, isapprox(p(-0.12345), f(-0.12345)))

	f = @(x) exp(sin(pi*x));
	n = 30;
	t = 2*(-n:n)'/(2*n+1) ;
	p = triginterp(t,f(t));
	assertTrue(testCase, isapprox(p(-0.12345) , f(-0.12345)))
	t = (-n:n-1)'/n;
	p = triginterp(t,f(t));
	assertTrue(testCase, isapprox(p(-0.12345) , f(-0.12345)))

	F = @(x) tan(x/2 - 0.2);
	f = @(x) 0.5*sec(x/2-0.2).^2;
	assertTrue(testCase, isapprox(ccint(f,40), F(1)-F(-1)))
    assertTrue(testCase, isapprox(glint(f,40), F(1)-F(-1)))

	f = @(x) 1./(32 + 2*x.^4);
	assertTrue(testCase, isapprox(intinf(f, 1e-9), sqrt(2)*pi/32, 1e-5))

	f = @(x) (1-x) ./ ( sin(x).^0.5 );
	assertTrue(testCase, isapprox(intsing(f,1e-8), 1.34312, 1e-5))
end

function chapter10test(testCase)
	lambda = 0.6;
	phi = @(r,w,dwdr) lambda./w^2 - dwdr./r;
	a = eps;  b = 1;
	ga = @(u,du) du;    
	gb = @(u,du) u-1;  

	[r,w,dwdx] = shoot(phi,a,b,ga,gb,[0.8;0], 1e-6);
	assertTrue(testCase, isapprox(w(1), 0.78776, 1e-4))

	f = @(x) exp(x.^2-3*x);
	df = @(x) (2*x-3).*f(x);
	ddf = @(x) ((2*x-3).^2+2).*f(x);

	[t,D,DD] = diffmat2(400,[-0.5,2]);
	assertTrue(testCase, isapprox(df(t), D*f(t), 1e-3))
    assertTrue(testCase, isapprox(ddf(t), DD*f(t), 1e-3))
	[t,D,DD] = diffcheb(80,[-0.5,2]);
	assertTrue(testCase, isapprox(df(t), D*f(t), 1e-7))
    assertTrue(testCase, isapprox(ddf(t), DD*f(t), 1e-7))

	exact = @(x)  exp(sin(x));
	p = @(x)  -cos(x);
	q = @sin;
	r = @(x) 0*x;
	[x,u] = bvplin(p,q,r,0,pi/2,1,exp(1),300);
	assertTrue(testCase, isapprox(u, exact(x), 1e-3))

	phi = @(t,theta,om) -0.05*om - sin(theta);
	g1 = @(u,du) u-2.5;
	g2 = @(u,du) u+2;
	init = linspace(2.5,-2,101)';

	[t,theta] = bvp(phi,0,5,g1,g2,init);
	assertTrue(testCase, isapprox(theta(7), 2.421850016880724,1e-10))

	c = @(x)  x.^2;
	 q = @(x)  4*ones(size(x));
	f = @(x)  sin(pi*x);
	[x,u] = fem(c,q,f,0,1,100);
	assertTrue(testCase, isapprox(u(33), 0.1641366907307196, 1e-10))
end

function chapter11test(testCase)
	s = @(x)  sin(pi*(x-0.2));
	c = @(x)  cos(pi*(x-0.2));
	f = @(x)  1 + s(x).^2;
	df = @(x)  2*pi*s(x).*c(x);
	ddf = @(x)  2*pi^2*(c(x).^2 - s(x).^2);

	[t, D, DD] = diffper(400,[0,2]);
	assertTrue(testCase, isapprox(df(t), D*f(t), 1e-3))
	assertTrue(testCase, isapprox(ddf(t), DD*f(t), 1e-3))

	phi = @(t,x,u,ux,uxx) uxx + t*u;
	g1 = @(u,ux) ux;
	g2 = @(u,ux) u-1;
	init = @(x) x.^2;
	[x, u] = parabolic(phi, [0,1], 40, g1, g2, [0, 2], init);
    uu = u(0.5);
	assertTrue(testCase, isapprox(uu(21), 0.845404, 1e-3))
    uu = u(1);
	assertTrue(testCase, isapprox(uu(end), 1, 1e-4))
    uu = u(2);
	assertTrue(testCase, isapprox(uu(1), 2.45692, 1e-3))
end

function chapter13test(testCase)
	f = @(x,y) -sin(3*x.*y-4*y) .* (9*y.^2 + (3*x-4).^2);
	g = @(x,y) sin(3*x.*y - 4*y);
	xspan = [0,1];  yspan = [0,2];
	[X, Y, U] = poissonfd(f, g, 60, xspan, 60, yspan);
    assertTrue(testCase, isapprox(g(X,Y), U, 1e-3))

	lambda = 1.5;
	function F = pde(X, Y, U, Ux, Uxx, Uy, Uyy)
		F = Uxx + Uyy - lambda./(U+1).^2;
	end
	g = @(x,y) x+y ;   
	u = elliptic(@pde, g, 30, [0,2.5], 24, [0,1]);
	assertTrue(testCase, isapprox(u(1.25, 0.5), 1.7236921361, 1e-6))
	assertTrue(testCase, isapprox(u(1,0), 1, 1e-6))
end

function result = isapprox(a, b, rtol, atol)
    if nargin < 4
        atol = 0;
        if nargin < 3
            rtol = 1e-10;
        end
    end
    result = norm(a - b) <= atol + rtol * max(norm(a), norm(b));
end