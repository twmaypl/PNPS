function w = chebfft1(v,n)
w = zeros(size(v));
if n==1
    w = chebfft(-v);
else
    for i = 1:n
        w(:,i)=chebfft(-v(:,i));
    end
end
end