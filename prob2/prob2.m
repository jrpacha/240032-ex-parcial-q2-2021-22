clearvars
close all

p=[1.15,0.66];
T=@(x,y) 16*cosh(x.*y);

eval('fluidMesh2');

numNodes=size(nodes,1);
numElem=size(elem,1);

temp=T(nodes(:,1),nodes(:,2));

for e=1:numElem
    nods=elem(e,:);
    vertexs=nodes(nods,:);
    [alphas,isInside]=baryCoord(vertexs,p);
    if isInside == 1
        fprintf('PART A\n')
        fprintf('The mean of the x-component of the three nodes of\n')
        fprintf('the element containing p is: %.4e\n',...
            sum(vertexs(:,1))/3)
        fprintf('Hint. The y-component of the third node of the\n')
        fprintf('element containing p is: %.4e\n',vertexs(3,2))

        fprintf('PART B\n')
        fprintf('The value of Psi3(p) is: %.4e\n', alphas(3))

        fprintf('PART C\n')
        fprintf('absErr:=|tempInterp(p)-temp(p)| = %.4e\n',...
            abs(alphas*temp(nods)-T(p(1,1),p(1,2))))
        fprintf('Hint. The maximum value of the temperature at the\n')
        fprintf('element''s nodes to which the point belongs is: %.4e\n',...
            max(temp(nods)))
        break;
    end
end



