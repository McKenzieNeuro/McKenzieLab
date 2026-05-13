function out = mysoftmax(in)
    % number of columns of in is the number of classes
    K = size(in,2);
    in = (1/(K-1)).*in;
    inmax = max(in,[],2);
    in = in-inmax;
    % in = bsxfun(@minus,in,inmax);
    numerator = exp(in);
    denominator = sum(numerator,2);
    denominator(denominator == 0) = 1;
    out= numerator./denominator;
end