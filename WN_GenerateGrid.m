% Generates a grid on (0,1)^2 with the given grid width.
% Returns a list of Nodes (in the 1st return value)
% and then a list of triangles, presented as column vectors
% where each column is a triangle, by giving the indices of the 3 nodes
function [Nodes,Elements] = WN_GenerateGrid(width)

    NoOfNodesPerRow = 1/width+1;

    Nodes=zeros(2,NoOfNodesPerRow^2);
    counter=1;
    % create the nodes of the central square (the inner nodes)
    for i = 1:NoOfNodesPerRow-2
        for j=1:NoOfNodesPerRow-2
            Nodes(:,counter) = [j*width;i*width];
            counter=counter+1;
        end
    end

    % create the nodes of the bottom row
    for i=1:NoOfNodesPerRow
        Nodes(:,counter) = [(i-1)*width;0];
        counter = counter + 1;
    end
    % create the nodes of the right column
    for i=2:NoOfNodesPerRow
        Nodes(:,counter) = [1;(i-1)*width];
        counter = counter + 1;
    end
    % create the nodes of the top row
    for i=NoOfNodesPerRow-1:-1:1
        Nodes(:,counter) = [(i-1)*width;1];
        counter = counter + 1;
    end
    % create the nodes of the left column
    for i=NoOfNodesPerRow-1:-1:2
        Nodes(:,counter) = [0; (i-1)*width];
        counter = counter + 1;
    end

    % we find the triangles by iterating over the squares and then adding to triangles for each square
    Elements=zeros(3,(NoOfNodesPerRow-1)^2*2);
    counter = 1;
    for i=0:NoOfNodesPerRow-2
        for j=0:NoOfNodesPerRow-2
            p1=[width*i;width*j];
            p2=p1+[width;0];
            p3=p1+[width;width];
            p4=p1+[0;width];
            % now we look for the indices of these points. Couldn't find a built-in matlab function that does this
            % so I loop over the columns of Nodes by hand and compare
            p1_index=0;
            p2_index=0;
            p3_index=0;
            p4_index=0;
            for k=1:size(Nodes)(2)
                if(Nodes(:,k) == p1)
                    p1_index=k;
                end
                if(Nodes(:,k) == p2)
                    p2_index=k;
                end
                if(Nodes(:,k) == p3)
                    p3_index=k;
                end
                if(Nodes(:,k) == p4)
                    p4_index=k;
                end
                if(p1_index!=0 && p2_index!=0 && p3_index!=0 && p4_index!=0)
                    break
                end
            end
            Elements(:,counter) = [p1_index;p2_index;p3_index];
            counter=counter+1;
            Elements(:,counter) = [p1_index;p4_index;p3_index];
            counter=counter+1;
        end
    end
end