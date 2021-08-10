function c = bisectionMethod(f,a,b,error)%f=@(x)x^2-3; a=1; b=2; (ensure change of sign between a and b) error=1e-4
c=(a+b)/2;
while abs(f(c))>error
    if f(c)<0&&f(a)<0
        a=c;
    else
        b=c;
    end
    c=(a+b)/2;
end