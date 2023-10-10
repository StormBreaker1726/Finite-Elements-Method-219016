function sh = shl(nen,nint)
if(nint == 2)
    pt(1) = -sqrt(3.)/3.;
    pt(2) = sqrt(3.)/3.;
    w(1) = 1.;
    w(2) = 1.;
end

for l=1:nint
    t=pt(l);
    if(nen==2)
        sh(1,l) = (1.0-t)/2.0;
        sh(2,l) = (1.0+t)/2.0;
    end
end
end