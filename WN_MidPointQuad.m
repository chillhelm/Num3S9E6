% quick and dirty implementation of the Edge Mid Point Quadrature rule
% Parameters are: 
% Triangle ... The triangle in a matrix of 3 column vectors
% f ... the function to be integrated
function integral=WN_MidPointQuad(Triangle,f)
    integral=0;
    integral=integral+f((Triangle(:,1)+Triangle(:,2))/2);
    integral=integral+f((Triangle(:,1)+Triangle(:,3))/2);
    integral=integral+f((Triangle(:,2)+Triangle(:,3))/2);
    Volume=abs(det([Triangle(:,2)-Triangle(:,1) Triangle(:,3)-Triangle(:,1)])/2.0);
    integral=integral*Volume/3;
end

