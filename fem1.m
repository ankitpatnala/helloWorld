nel = 5;
nnel = 2;
ndof = 1;
nnode = 6;
sdof = nnode*ndof;
gcoord = [ 0.0 0.2 0.4 0.6 0.8 1.0 ];
nodes = [ 1 2 ; 2 3 ; 3 4 ; 4 5 ; 5 6 ];
acoef = 1;
bcoef = 3; 
ccoef = 2;
bcdof(1) = 1 ;
bcval(1)= 0.0;
bcdof(2) = 6;
bcval(2) = 0.0;

ff = zeros(sdof,1);
kk = zeros(sdof,sdof);
index = zeros(nnel*ndof,1);

for iel = 1 : nel 
    nl = nodes(iel,1); nr = nodes(iel,2);
    xl = gcoord(nl); xr = gcoord(nr);
    eleng = xr -xl;
    index = feeldof1(iel,nnel,ndof);
    
    k = feode21(acoef,bcoef,ccoef,eleng);
    f = fe1l(xl,xr);
    [kk,ff] = feasmbl2(kk,ff,k,f,index);
end

[kk , ff] = feaplyc2(kk,ff,bcdof,bcval);
    
fsol = kk\ff;

c1 = 0.5/exp(1);
c2 = 0.5*(1+1/exp(1));
for i = 1:nnode
    x = gcoord(i);
    esol(i) = c1*exp(2*x)+c2*exp(x)+1/2;
end

num = 1:1:sdof;
results = [num' fsol esol'];
plot(num , fsol);







