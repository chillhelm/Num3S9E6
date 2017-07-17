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
%   along with Num3S9E6.  If not, see <http://www.gnu.org/licenses/>.



function err=WN_ComputeErrors(Grid, GridNodes,uh,solution,solgrad)

    GridSize = size(Grid,2);
    NodeCount=size(GridNodes,2);

    H1 = 0;
    L2 = 0;

    for Element=Grid
        %compute transformation to current grid element from reference triangle
        LocalTransformJacobian = [GridNodes(:,Element(2))-GridNodes(:,Element(1)) GridNodes(:,Element(3))-GridNodes(:,Element(1))];
        LocalTransformJacobianInv = inv(LocalTransformJacobian);
        
        % compute the gradient of Uh on this triangle
        GradUh=[uh(Element(2))-uh(Element(1)) uh(Element(3))-uh(Element(1))]*LocalTransformJacobianInv;
        % compute H1 semi-norm on this triangle
        h1 = WN_MidPointQuad(GridNodes(:,Element'),@(x)dot(GradUh-solgrad(x),GradUh-solgrad(x)));
        H1 = H1+h1;

        % compute L2 norm on this triangle
        uh_interp = @(x)uh(Element(1))*(1-x(1)-x(2))+uh(Element(2))*x(1)+uh(Element(3))*x(2);
        l2 = WN_MidPointQuad(GridNodes(:,Element'),@(x)(solution(x)-uh_interp(x))^2);
        L2 = L2+l2;
    end

    % take square roots
    err=[sqrt(L2);sqrt(H1)];

end

