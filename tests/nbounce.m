function nbounce(p,q)

if mod(p,2)==1 && mod(q,2)==2
    k=1;
else
    k=1/2;
end

n = k*(2*p+min(2*p,6*q))

end