///////////// Computational Methods ///////////////
dim = 3
npts = 4
eps = 1e-20

function x = bisection(f_, a, b, l0, cc)
    c = (a+b)*0.5
    
    i = 0;
    for i=1:6
        fc = f_(c, l0, cc)
        i = i+1
        y = (a+c)*0.5
        fy = f_(y, l0, cc)
        if(fy <= fc) 
            b = c
            c = y
        else
            z = (c+b)*0.5
            fz = f_(z, l0, cc)
            if(fc <= fz)
                a = y
                b = z
            else
                a = c
                c = z
            end
        end
    end
    x = c
endfunction


function [x, c, s1, s2] = grad_desc(x0, l0, f_, lf_, df_, alpha, c, maxit)
    err = 1
    i=0
    s1 = zeros(maxit+1,1)
    s2 = zeros(maxit+1,1)
    last_grad = zeros(x)
    while(err > eps & i < maxit)
        s1(i+1)=f_(normalize(x0))
        s2(i+1)=lf_(x0, l0, c)
        grad = df_(x0, l0, c)
        x = bisection(lf_, x0, x0 - alpha * grad, l0, c)
        last_grad = grad
        
        if(sum(isnan(x))|sum(isinf(x)))
            maxit=i
            x = x0
        else
            err = norm(grad)
            x0 = x
            i = i + 1
        end
    end
    s1 = s1(1:i)
    s2 = s2(1:i)
    x = normalize(x)
endfunction

function [x, l, s1, s2, s3] = optimize(x0, l0, f_, lf_, hf_, df_, alpha1, alpha2, c, maxit, smaxit)
    s1 = []
    s2 = []
    i = 0
    err = 1
    x = x0
    l = l0
    s3 = list()
    while(err > eps & i < maxit)

        if(i < maxit/2)
            [x, c, ts1, ts2] = grad_desc(x, l, f_, lf_, df_, alpha1, c, smaxit)
        else
            t = (i - maxit/2)/(maxit/2)
            alpha  = (1-t)*alpha1 + t*alpha2
            [x, c, ts1, ts2] = grad_desc(x, l, f_, lf_, df_, alpha, c, smaxit)
        end
        l = l + c * hf_(x)
        i = i + 1
        
        err = norm(df_(x))
        s3($+1)=x;
        s1 = [s1; ts1]
        s2 = [s2; ts2]
    end
    
endfunction

/////////////F(x) ///////////////

function [y, ii, jj] = f(x)
    y = -1000
    ii = 1
    jj = 2
    for i=1:npts
        xi = getv(x,i)
        for j=i+1:npts
            xj = getv(x,j)
            if(xi'*xj > y)
                y = xi'*xj
                ii = i
                jj = j
            end
        end
    end
endfunction

function y = lf(x, l0, c)
    h = hf(x);
    y = f(x) + l0'*h + c.*h'*h * 0.5;
endfunction

function dx = df(x, l, c)
    dx = zeros(x)
    dl = zeros(x)
    [y, j, i] = f(x)
    
    xj = getv(x,j)
    xi = getv(x,i)
    dx = setv(dx,xj,i)
    dx = setv(dx,xi,j)
    
    for i=1:npts
        xi = getv(x,i)
        dxi = 2 * l(i) * xi
        nxi =(xi'*xi)-1
        dl = setv(dl,dxi + c.*xi.*nxi, i)
    end
    dx = dx + dl
endfunction

function hx = hf(x)
    hx = zeros(npts,1)
    for i=1:npts
        xi = getv(x,i)
        hx(i) = (xi'*xi)-1
    end
endfunction

///////////// G(x) ///////////////


function y = g(x)
    y=0
    for i=1:npts
        xi = getv(x,i)
        for j=i+1:npts
            xj = getv(x,j)
            xij = xi-xj
            y = y + 1/sqrt(xij'*xij)
        end
    end
endfunction

function dx = dg(x, l, c)
    dx = zeros(x)
    dl = zeros(x)
    
    for i=1:npts
        xi = getv(x,i)
        tmp = zeros(xi)
        for j=1:npts
            xj = getv(x,j)
            if(i~=j)
                xij = xi-xj
                nij = xij'*xij
                tmp = tmp -2*xij/(nij*nij)
            end
        end
        dx = setv(dx, tmp, i)
    end
    
    for i=1:npts
        xi = getv(x,i)
        dxi = 2 * l(i) * xi
        nxi = (xi'*xi)-1
        dl = setv(dl,dxi+c*xi*nxi, i)
    end
    dx = dx + dl
endfunction

function hx = hg(x)
    hx = zeros(npts,1)
    for i=1:npts
        xi = getv(x,i)
        hx(i) = (xi'*xi)-1
    end
endfunction

function y = lg(x, l0, c)
    h = hg(x)
    y = g(x) + l0'*h + c*h'*h * 0.5
endfunction

///////////// Utils ///////////////

function x = normalize(x)
    for i=1:npts
        xi = getv(x,i)
        xi = xi/norm(xi)
        x = setv(x,xi,i)
    end
endfunction

function xi = getv(x,i)
    si = (i-1)*dim + 1
    xi = x(si:si+dim-1)
endfunction

function x = setv(x,xi,i)
    si = (i-1)*dim + 1
    x(si:si+dim-1) = xi
endfunction



function a = min_angle(x)
    a = acos(f(x)).*180/%pi
endfunction



///////////// Initialization ///////////////

function x0 = initEq()
    x0 = zeros((dim)*npts, 1)
    pi2 = %pi/4
    for i=1:npts
        theta = 2*%pi*i/npts
        xi = [cos(pi2) 0 -sin(pi2);0 1 0; sin(pi2) 0 cos(pi2)]*[sin(theta); cos(theta); 0]
        x0 = setv(x0, xi, i)
    end
    x0 = normalize(x0)
endfunction

function x0 = initRand()
    x0 = rand((dim)*npts, 1)-0.5
    x0 = normalize(x0)
endfunction

function x0 = initRM(i, f)
    x0 = rand((dim)*npts, 1)-0.5
    x0 = normalize(x0)
    fx0 = f(x0)
    for k=1:i
        x = rand((dim)*npts, 1)-0.5
        x = normalize(x)
        fx = f(x)
        if(fx < fx0)
            x0 = x
            fx0 = fx
        end
    end
endfunction

function x0 = initRM2(i, f)
    x0 = rand((dim)*int(npts/2), 1)-0.5
    x0 = [x0; -x0]
    if(modulo(npts, 2) == 1) 
        x0 = [x0 ; rand(dim, 1)-0.5]
    end
    x0 = normalize(x0)
    fx0 = f(x0)
    for k=1:i
        x = rand((dim)*npts, 1)-0.5
        x = normalize(x)
        fx = f(x)
        if(fx < fx0)
            x0 = x
            fx0 = fx
        end
    end
endfunction

/////////////// Run ///////////////

result =[]
l=list()

form = ['b', 'b--', 'r', 'r--', 'r-.', 'g', 'g--', 'g-.']
leg = []
cf = 2
cg = 0.5
af = 0.025
ag = 0.1
ai = 0.0

for i=4:40
    tic()
    npts = i
    
   //Use same seed for comparison
    rand('seed',1)
    l0 = rand(npts,1)
    rand('seed',1)
    x0f = initRM2(5000, f)
    rand('seed',1)
    x0g = initRM2(5000, g)
    [x1, l1, s11, s12] = optimize(x0f, l0, f, lf, hf, df, af, ai, cf, 500, 5);
    [x2, l2, s21, s22] = optimize(x0g, l0, g, lg, hg, dg, ag, ai, cg, 500, 5);
    r = [min_angle(x1), min_angle(x2), g(x1), g(x2), min_angle(x0f), min_angle(x0g),g(x0f), g(x0g), norm(df(x1,l1,c)), norm(dg(x2,l2,c))]
    disp(r, "result for #" + string(i) + " points")
    result = [result; r]
    l($+1)=list(x0f, x0g, x1, x2, l1, l2, s11, s21, s12, s22)

    disp(min_angle(x1), 'f(x): ')
    disp(g(x2), 'g(x): ')
    disp(toc(), "time taken: ")
end



function plotp3d(x)
    l = length(x)
    x_ = x(1:3:l)
    y_ = x(2:3:l)
    z_ = x(3:3:l)
    param3d(x_, y_, z_)
    e=gce()
    e.line_mode="off"
    e.line_style=3
    e.mark_mode="on"
    e.mark_style=3
    a=gca()
    a.isoview="on"
endfunction

function plotTraj(l)
    lt = list()
    for i=1:npts
        m = []
        ll = length(l)
        for j=1:ll
            m = [m, getv(l(j),i)]
        end
        lt($+1)=m
        param3d(m(1,:), m(2,:), m(3,:))
        e=gce()
        e.line_mode="on"
        //e.line_style=3
        e.thickness=2
        param3d(m(1,ll), m(2,ll), m(3,ll))
        e=gce()
        e.line_mode="off"
        e.line_style=3
        e.mark_mode="on"
        e.mark_style=3
    end
    a=gca()
    a.isoview="on"
    
endfunction



/////////////// Stochastic Heuristic for comparison ///////////////

function p = normalizep(p)
    s = size(p)
    for i=1:s(2)
        p(:,i) = normalize(p(:,i))
    end
endfunction

function [x, fv] = de(fx, pop, init, cr, fr, maxit)
    p = rand(dim*npts, pop)
    p = normalizep(p)
    fv = zeros(pop,1)
    for i=1:pop
        fv(i) = fx(p(:,i))
    end
    disp(min(fv))
    for i=1:init
        tmp = normalize(rand(dim*npts,1))
        ft = fx(tmp)
        j = find(fv >= max(fv))
        if(ft < fv(j))
            p(:,j) = tmp
            fv(j)=ft
        end
    end
    disp(min(fv))
    i=0
    
    while(i < maxit & (max(fv) - min(fv)) > eps)
        i = i+1
        pp = p
        for k=1:pop
            idx=[]
            while(length(idx) < 3)
                idx = unique(grand(3,1,'uin',1,pop))
            end
            
            tmp = pp(:, idx(1)) + fr * (pp(:, idx(2)) - pp(:, idx(3)))
            rr = grand(1,1,'uin',1,dim*npts)
            r = rand(tmp)
            r(rr)=0
            for l=1:dim*npts
                if(r(l)>cr)
                    tmp(l) = p(l,k)
                end
            end
            
            tmp = normalize(tmp)
            ft = fx(tmp)
            if(ft < fv(k))
                fv(k)=ft
                pp(:,k)=tmp
            end
        end
        p = pp
     end
     disp(min(fv))
     disp("iterations: " + string(i))
     
     j = find(fv <= min(fv))
     x = p(:,j(1))
endfunction


/////////////// Test function ///////////////

function dx = dft(x, l0)
    dx = zeros(length(x))
    for i=1:length(x)/2
        dx(2*i-1) = 100 * 2 * x(2*i-1) * 2 * (x(2*i-1)^2 - x(2*i)) + 2 * (x(2*i-1) - 1)
        dx(2*i) = -100 * 2 * (x(2*i-1)^2 - x(2*i))
    end
endfunction

function hf = hft(x)
     hf = 0
endfunction

function y=ft(x)
    y = 0
    for i=1:length(x)/2
        y = y + 100 *(x(2*i-1)^2 - x(2*i))^2 + (x(2*i-1) - 1)^2
    end
endfunction




