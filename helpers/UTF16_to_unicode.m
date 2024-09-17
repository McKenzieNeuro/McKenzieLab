function unicode=UTF16_to_unicode(UTF16)
%Convert UTF-16 to the code points stored as uint32
%
%See https://en.wikipedia.org/wiki/UTF-16
%
% 1 word (U+0000 to U+D7FF and U+E000 to U+FFFF):
%  xxxxxxxx_xxxxxxxx
% 2 words (U+10000 to U+10FFFF):
%  110110xx_xxxxxxxx 110111xx_xxxxxxxx

persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
UTF16=uint32(UTF16);

multiword= UTF16>55295 & UTF16<57344; %0xD7FF and 0xE000
if ~any(multiword)
    unicode=UTF16;return
end

word1= find( UTF16>=55296 & UTF16<=56319 );
word2= find( UTF16>=56320 & UTF16<=57343 );
try
    d=word2-word1;
    if any(d~=1)
        error('trigger error')
    end
catch
    error('input is not valid UTF-16 encoded')
end

%Binary header:
% 110110xx_xxxxxxxx 110111xx_xxxxxxxx
% 00000000 01111111 11122222 22222333
% 12345678 90123456 78901234 56789012
header_bits='110110110111';header_locs=[1:6 17:22];
multiword=UTF16([word1.' word2.']);
multiword=unique(multiword,'rows');
S2=mat2cell(multiword,ones(size(multiword,1),1),2);
unicode=UTF16;
for n=1:numel(S2)
    bin=dec2bin(double(S2{n}))';
    
    if ~strcmp(header_bits,bin(header_locs))
        error('input is not valid UTF-16 encoded')
    end
    bin(header_locs)='';
    if ~isOctave
        S3=uint32(bin2dec(bin  ));
    else
        S3=uint32(bin2dec(bin.'));%Octave needs an extra transpose.
    end
    S3=S3+65536;% 0x10000
    %Perform actual replacement.
    unicode=PatternReplace(unicode,S2{n},S3);
end
end
function out=bsxfun_plus(in1,in2)
%Implicit expansion for plus(), but without any input validation.
try
    out=in1+in2;
catch
    try
        out=bsxfun(@plus,in1,in2);
    catch
        sz1=size(in1);                    sz2=size(in2);
        in1=repmat(in1,max(1,sz2./sz1));  in2=repmat(in2,max(1,sz1./sz2));
        out=in1+in2;
    end
end
end
function out=PatternReplace(in,pattern,rep)
%Functionally equivalent to strrep, but extended to more data types.
out=in(:)';
if numel(pattern)==0
    L=false(size(in));
elseif numel(rep)>numel(pattern)
    error('not implemented (padding required)')
else
    L=true(size(in));
    for n=1:numel(pattern)
        k=find(in==pattern(n));
        k=k-n+1;k(k<1)=[];
        %Now k contains the indices of the beginning of each match.
        L2=false(size(L));L2(k)=true;
        L= L & L2;
        if ~any(L),break,end
    end
end
k=find(L);
if ~isempty(k)
    for n=1:numel(rep)
        out(k+n-1)=rep(n);
    end
    if numel(rep)==0,n=0;end
    if numel(pattern)>n
        k=k(:);%Enforce direction.
        remove=(n+1):numel(pattern);
        idx=bsxfun_plus(k,remove-1);
        out(idx(:))=[];
    end
end
end