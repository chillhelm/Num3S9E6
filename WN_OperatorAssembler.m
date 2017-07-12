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



function A=WN_OperatorAssembler(Grid, GridNodes)

    GridSize = size(Grid,2);
    NodeCount=size(GridNodes,2);

    LocalBasisGradients=[[-1 -1];[1 0];[0 1]];

    % 5 entries per row might be enough? I wasnt sure so I went with a bit extra ...
    A=spalloc(NodeCount,NodeCount,NodeCount*7);

    %iterate over all grid elements
    for Element=Grid
        % compute the transformation to this triangle ...
        LocalTransformJacobian = [GridNodes(:,Element(2))-GridNodes(:,Element(1)) GridNodes(:,Element(3))-GridNodes(:,Element(1))];
        % and it's inverse
        LocalTransformJacobianInv = inv(LocalTransformJacobian);
        
        % Dimension of local basis (1 for each node)
        LocalSize=size(Element,1);
        for j=1:LocalSize
            % compute the gradient of phi_j
            GradPhiJ=LocalBasisGradients(j,:)*LocalTransformJacobianInv;
            for i=1:LocalSize
                %compute the gradient of phi_i
                GradPhiI=LocalBasisGradients(i,:)*LocalTransformJacobianInv;
                %integrate their product over the triangle in question
                local_update = WN_MidPointQuad(GridNodes(:,Element'),@(x)GradPhiI*GradPhiJ');
                %update the matrix
                A(Element(j),Element(i)) = A(Element(j),Element(i)) + local_update;
            end
        end
    end

    % the boundary elements have to be treated special
    for j=(sqrt(NodeCount)-2)^2+1:NodeCount
        A(j,:)=zeros(1,size(A,2));
        A(j,j)=1;
    end

    %todo: sparsenes of A
end

