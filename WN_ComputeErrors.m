function err=WN_ComputeErrors(Grid, GridNodes,uh,solution,solgrad)

    GridSize = size(Grid)(2);
    NodeCount=size(GridNodes)(2);

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

