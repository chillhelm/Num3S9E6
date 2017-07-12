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



function b=WN_RHSAssembler(Grid, GridNodes, rhs_function, boundary_terms)

    NodeCount=size(GridNodes,2);

    % the local basis on the reference triangle
    LocalBasisFunction1=@(x)1-x(1)-x(2);
    LocalBasisFunction2=@(x)x(1);
    LocalBasisFunction3=@(x)x(2);

    b=zeros(NodeCount,1);

    
    for Element=Grid
        
        % transformation to the current grid element
        LocalTransformJacobian = [GridNodes(:,Element(2))-GridNodes(:,Element(1)) GridNodes(:,Element(3))-GridNodes(:,Element(1))];
        LocalTransformJacobianInv = inv(LocalTransformJacobian);

        % integrate phi_1*rhs over the triangle
        FunctionToIntegrate=@(x)LocalBasisFunction1(LocalTransformJacobianInv*(x-GridNodes(:,Element(1))))*rhs_function(x);
        local_update = WN_MidPointQuad(GridNodes(:,Element'),FunctionToIntegrate);
        
        %update the rhs vector
        b(Element(1)) = b(Element(1)) + local_update;

        %do the same thing for phi_2*rhs and phi_3*rhs
        
        FunctionToIntegrate=@(x)LocalBasisFunction2(LocalTransformJacobianInv*(x-GridNodes(:,Element(2))))*rhs_function(x);

        local_update = WN_MidPointQuad(GridNodes(:,Element'),FunctionToIntegrate);
        b(Element(2)) = b(Element(2)) + local_update;

        FunctionToIntegrate=@(x)LocalBasisFunction3(LocalTransformJacobianInv*(x-GridNodes(:,Element(3))))*rhs_function(x);

        local_update = WN_MidPointQuad(GridNodes(:,Element'),FunctionToIntegrate);
        b(Element(3)) = b(Element(3)) + local_update;
    end

    % special stuff for boundary nodes
    for j=(sqrt(NodeCount)-2)^2+1:NodeCount
        b(j)=boundary_terms(GridNodes(:,j));
    end
end

