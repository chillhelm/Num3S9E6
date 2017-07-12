%   Copyright 2017, Wilhelm Neubert
%   This file is part of Num3S9E6.
%   
%   Num3S9E6 is free software: you can redistribute it and/or modify
%   it under the terms of the GNU Lesser Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%   
%   Num3S9E6 is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser Public License for more details.
%   
%   You should have received a copy of the GNU Lesser Public License
%   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.



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

