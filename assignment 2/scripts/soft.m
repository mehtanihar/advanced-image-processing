function [x]=soft(y,lambda)
    x = sign(y).*( max( 0, abs(y)-lambda ) );

end
