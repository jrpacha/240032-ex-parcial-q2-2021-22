function varargout = funcioProb3(dudx,a,alpha)
% If alpha is present in the input argument list
%   varargout{1} = minU
%   varargout{2} = interpU
%   varargout{3} = u
% or
%   varargout{1} = alpha s.t. u(27)-u(70) = 2*u(52)
%   varargout{2} = u
%   varargout{3} = ''
% 

numDiv=100;
x=0.337;
x0=0;
xL=1;
a1=1.0; beta=1.0; a0=3.0;

h=(xL-x0)/numDiv;
nodes=(x0:h:xL)';
elem=[1:numDiv; 2:numDiv+1]';
numNodes=size(nodes,1);
numElem=size(elem,1);

K=zeros(numNodes);
F=zeros(numNodes,1);
Q=zeros(numNodes,1);
u=zeros(numNodes,1);

Ke1 = beta*[1,-1;-1,1]/h + a0*h*[2,1;1,2]/6;

for e=1:numElem
    rows=[elem(e,1), elem(e,2)];
    cols=rows;
    x1=nodes(rows(1,1),1);
    x2=nodes(rows(1,2),1);
    Ke = 0.5*a1*(x1+x2)*[1,-1;-1,1]/h+Ke1;
    K(rows,cols)=K(rows,cols)+Ke;
    if a ~= 0
        Fe = a*h*[2*x1+x2;x1+2*x2]/6;
        F(rows) = F(rows) + Fe;
    end
end


%B.C 
fixedNods=1;
freeNods=setdiff(1:numNodes,fixedNods);

%Nautural
Q(numNodes)=dudx*(a1*nodes(numNodes)+beta);

switch nargin
    case 3
        %Essential B.C.
        u(1)=alpha;

        %Reduced system
        Qm=Q(freeNods)+F(freeNods)-K(freeNods,fixedNods)*u(fixedNods);
        Km=K(freeNods,freeNods);

        %Solve the reduced system
        um=Km\Qm;
        u(freeNods)=um;
        varargout{1}=min(u);
        varargout{2}=interp1(nodes,u,x);
        varargout{3}=u;
    case 2        
        Km = zeros(numNodes); %If you want to keep the original K
        Qm = zeros(numNodes,1); %If you want to keep the original Q
    
        %"Reduced" system
        Qm(1:numNodes-1)=Q(freeNods)+F(freeNods);
        Km(1:numNodes-1,:)=K(freeNods,:);
        Km(numNodes,27)=1;
        Km(numNodes,52)=-2;
        Km(numNodes,70)=-1;    
    
        %Solve the reduced system
        u=Km\Qm;
        varargout{1}=u(1);
        varargout{2}=u;
        varargout{3}='';
    otherwise
        error('You must pass 2 or 3 arguments')
end
    
end 