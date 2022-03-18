
function cost = EMD(spiketrain1,spiketrain2,ti,tf)


if ~isempty(spiketrain1)
    X = [ti,repelem(spiketrain1',2),tf];
    n = 1/length(spiketrain1);
    Y = zeros(1,length(spiketrain1)-1);
    for i = 1:(length(spiketrain1)-1)
        Y(i) = i*n;
    end
    Y = [0,0,repelem(Y,2),1,1];
    Area1 = trapz(X,Y);
end


if ~isempty(spiketrain2)
    X = [ti,repelem(spiketrain2',2),tf];
    n = 1/length(spiketrain2);
    Y = zeros(1,length(spiketrain2)-1);
    for i = 1:(length(spiketrain2)-1)
        Y(i) = i*n;
    end
    Y = [0,0,repelem(Y,2),1,1];
    Area2 = trapz(X,Y);
end


if ~isempty(spiketrain1) & ~isempty(spiketrain2)
    cost = abs(Area1-Area2);
elseif ~isempty(spiketrain1) & isempty(spiketrain2)
    cost = Area1;
elseif isempty(spiketrain1) & ~isempty(spiketrain2)
    cost = Area2;
else
    cost = 0;
end

end