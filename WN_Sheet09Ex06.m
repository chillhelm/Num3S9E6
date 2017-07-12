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




boundary_terms=@(x)x(1)^4*x(2)^5-17*sin(x(1)*x(2)); %this is the correct solution as well
rhs_function=@(x) 20*x(1)^4*x(2)^3 + 12*x(1)^2*x(2)^5 + 17*x(1)^2*sin(x(1)*x(2)) + 17*x(2)^2*sin(x(1)*x(2)); % the laplace of the solution
solgrad=@(x)[4*x(1)^3*x(2)^5-17*x(2)*cos(x(1)*x(2)) 5*x(1)^4*x(2)^4-17*x(1)*cos(x(1)*x(2)) ]; % gradient of the solution (needed for error computation
Errors=[];
for n=2:8
    h=2^-n;
    %Generate a grid with the appropriate width
    [Nodes,Elements]=WN_GenerateGrid(h);
    %Assemble the linear FEM operator
    A=WN_OperatorAssembler(Elements,Nodes);
    %Assemble the right hand side
    b=WN_RHSAssembler(Elements,Nodes, rhs_function, boundary_terms);
    %Solve
    uh=A\b;
    %Compute errors
    Errors = [Errors [n;h; WN_ComputeErrors(Elements,Nodes,uh,boundary_terms, solgrad)]];
    %Compute convergence rate
    if size(Errors,2)>1
      rate_L2 = (log(Errors(3,end))-log(Errors(3,end-1)))/(log(Errors(2,end))-log(Errors(2,end-1)));
      rate_H1 = (log(Errors(4,end))-log(Errors(4,end-1)))/(log(Errors(2,end))-log(Errors(2,end-1)));
    end
end

%Show a table with the errors
Errors



s = sprintf('The order of convergence in the L2 based on finest meshes is %g',rate_L2);
disp(s);
t = sprintf('The order of convergence in the H1 semi norm based on finest meshes is %g',rate_H1);
disp(t);

