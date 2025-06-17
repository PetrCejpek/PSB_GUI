function I=TrapezoidIntegral(x,y)
I=sum((x(2:end)-x(1:end-1)).*(y(2:end)+y(1:end-1))/2);