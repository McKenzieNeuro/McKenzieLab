
function c = redblue_local(n)
    r=[linspace(0.20,1.00,n/2) linspace(1.00,0.85,n/2)]';
    g=[linspace(0.40,0.90,n/2) linspace(0.90,0.20,n/2)]';
    b=[linspace(0.85,0.90,n/2) linspace(0.90,0.10,n/2)]';
    c=[r g b];
end