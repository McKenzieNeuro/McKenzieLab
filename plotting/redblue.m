function c = redblue(n)
    r = [linspace(0.8,1,n/2) linspace(1,0.8,n/2)]';
    g = [linspace(0.8,0.2,n/2) linspace(0.2,0.8,n/2)]';
    b = [linspace(1,0.8,n/2) linspace(0.8,1,n/2)]';
    c = [r g b];
    c = flipud(c);
end