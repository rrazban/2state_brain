function Xnew= Isingify2(t,X)

Xnew = X(2:t,:)-X(1:t-1,:);
Xnew( Xnew < 0 ) = -1;
Xnew( Xnew > 0 ) = 1;

end
