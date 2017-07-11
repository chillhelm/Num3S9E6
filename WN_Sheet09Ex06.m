
boundary_terms=@(x)x(1)^4*x(2)^5-17*sin(x(1)*x(2)); %this is the correct solution as well
rhs_function=@(x) 20*x(1)^4*x(2)^3 + 12*x(1)^2*x(2)^5 + 17*x(1)^2*sin(x(1)*x(2)) + 17*x(2)^2*sin(x(1)*x(2)); % the laplace of the solution
solgrad=@(x)[4*x(1)^3*x(2)^5-17*x(2)*cos(x(1)*x(2)) 5*x(1)*4*x(2)^4-17*x(1)*cos(x(1)*x(2)) ] % gradient of the solution (needed for error computation
Errors=[]
for n=2:8
    h=2^-n
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
    if n>2
      rate_L2 = (log(Errors(3,end))-log(Errors(3,end-1)))/(log(Errors(2,end))-log(Errors(2,end-1)));
      rate_H1 = (log(Errors(4,end))-log(Errors(4,end-1)))/(log(Errors(2,end))-log(Errors(2,end-1)));
    end
end

%Show a table with the errors
Errors



s = sprintf('The order of convergence in the sup-norm based on finest meshes is %g',rate_L2);
disp(s);
t = sprintf('The order of convergence in the Euclidean norm based on finest meshes is %g',rate_H1);
disp(t);