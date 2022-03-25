function [xf, fobj, info, VF, VFmean, optTime] = Hip_Implant(varargin)
%-------------------------------------------------------------------------%
% inputs = filtr,kfiltr,VF,MaxVF,MinDensEle,MaxDensEle
% filtr = 0 No filter
%         1 Density filter
% kfiltr = filter size parameter (filter radius = kfiltr*(element size))
% VF = Global volume fraction
% MaxVF = Maximum global volume fraction
% MinDensEle = Minimum volume fraction per element (non-linear constraint)
% MaxDensEle = Maximum volume fraction per element (non-linear constraint)
%-------------------------------------------------------------------------%
%% FEM
% Geometry
model = mphload('Hip_Implant.mph');
[~,meshdata] = mphmeshstats(model,'mesh2'); % Mesh with 14894 elements
% [~,meshdata] = mphmeshstats(model,'mesh3'); % Mesh with 848247 elements
nodes = meshdata.vertex';
XYZ = nodes*1e3;
elements = double(meshdata.elem{2}+1)';
CONEC(:,1) = elements(:,3);
CONEC(:,2) = elements(:,1);
CONEC(:,3) = elements(:,2);
CONEC(:,4) = elements(:,4);
nele = size(CONEC, 1);
nnod = size(XYZ, 1);
P1 = XYZ(CONEC(:,1),:);
P2 = XYZ(CONEC(:,2),:);
P3 = XYZ(CONEC(:,3),:);
P4 = XYZ(CONEC(:,4),:);
V = 1/6*abs(dot(P2-P1,cross(P3-P1,P4-P1,2),2));


% Boundary nodes
bnodes = nodes(1:max(max(meshdata.elem{3}))+1,:);


% Boundary conditions
LZ = 0;
base = find(bnodes(:,3)<LZ);
MFIXDOF = [3*base(:)-2; 3*base(:)-1; 3*base(:)];
MFIXVAL = MFIXDOF*0;


% Loads
LZ = 0.102;
Fext = zeros(3*nnod,1);

pnod = find(bnodes(:,3)>LZ);
P = 1;
pdof = 3*pnod;
Fext(pdof) = -P/length(pnod);


% Material
mater = struct('ElasticModuli',@(x) DEMOD(x),'ConstitutiveLaw', 'LinearElastic');


% Nodes and elements
EleType = 'TETRA4';
NOD = struct('Coordinates',XYZ,'ExternalForce',Fext,'UFixedDOFs',MFIXDOF,'UFixedValues',MFIXVAL);
ELE = struct('Type',EleType,'Connectivity',CONEC,'Material',mater,'Volumes',V);

ELE.dom1 = find(meshdata.elementity{1,2}==1);
ELE.dom2 = find(meshdata.elementity{1,2}==2);
ELE.dom3 = find(meshdata.elementity{1,2}==3);


%% Constraints
% Bounds
constr.xmin = zeros(7*length(ELE.dom2),1);
constr.xmin(1:7:end) = 0.06;
constr.xmin(2:7:end) = 0.06;
constr.xmin(3:7:end) = 0.06;
constr.xmin(4:7:end) = 0.6;
constr.xmin(5:7:end) = 0;
constr.xmin(6:7:end) = 0;
constr.xmin(7:7:end) = 0;
constr.xmax = zeros(7*length(ELE.dom2),1);
constr.xmax(1:7:end) = 0.95;
constr.xmax(2:7:end) = 0.95;
constr.xmax(3:7:end) = 0.95;
constr.xmax(4:7:end) = 1.1;
constr.xmax(5:7:end) = pi;
constr.xmax(6:7:end) = pi;
constr.xmax(7:7:end) = pi;

% Volume Fraction
if nargin > 2
    if varargin{3} > 0
        constr.AvrgDens = varargin{3};
    end
end
if ~isfield(constr,'AvrgDens')
    constr.AvrgDens = -inf;
end
if nargin > 3
    if varargin{4} > 0
        constr.MaxAvrgDens = varargin{4};
    end
end
if ~isfield(constr,'MaxAvrgDens')
    constr.MaxAvrgDens = -inf;
end
if nargin > 4
    if varargin{5} > 0
        constr.MinDensEle = varargin{5};
    end
end
if ~isfield(constr,'MinDensEle')
    constr.MinDensEle = -inf;
end
if nargin > 5
    if varargin{6} > 0
        constr.MaxDensEle = varargin{6};
    end
end
if ~isfield(constr,'MaxDensEle')
    constr.MaxDensEle = -inf;
end


%% Filter
if nargin > 0
    filtr = varargin{1};
else
    filtr = 0;
end
if nargin > 1
    kfiltr = varargin{2};
else
    filtr = 0;
    kfiltr = 0;
end
H = CalcSideLength(NOD.Coordinates,ELE.Connectivity,ELE.Type);
hmax = max(H,[],2);
rfiltr = kfiltr*hmax;
Vars2Filtr = 1:7;


%% Optimization
% Seed
xi = zeros(7*length(ELE.dom2),1);
xi(1:7:end) = 0.95*0.9999;
xi(2:7:end) = 0.95*0.9999;
xi(3:7:end) = 0.95*0.9999;
xi(4:7:end) = 1.1*0.9999;
xi(5:7:end) = pi*0.9999;
xi(6:7:end) = pi*0.9999;
xi(7:7:end) = pi*0.9999;

% Termination options
termOpt.maxIter = 10000;
termOpt.desTol = 1e-5;
termOpt.accTol = 1e-4;
termOpt.accIter = 10;

% Output file
printlevel = 5;
file = ['Hip_Implant_NEle',num2str(nele)];
if constr.AvrgDens > 0
    file = [file,'_VF',num2str(constr.AvrgDens*100)];
end
if constr.MaxAvrgDens > 0
    file = [file,'_MaxVF',num2str(constr.MaxAvrgDens*100)];
end
if constr.MinDensEle > 0
    file = [file,'_MinVFEle',num2str(constr.MinDensEle*100)];
end
if constr.MaxDensEle > 0
    file = [file,'_MaxVFEle',num2str(constr.MaxDensEle*100)];
end
if filtr
    file = [file,'_kf',num2str(kfiltr,'%0.1f')];
end
file = [file,'_IPOpt'];
save([file '.mat'],'ELE','NOD','termOpt','constr','xi');

% Run
VarElems = [];
opt = 0;
[xf, fobj, info, VF, VFmean, optTime] = OptMain(opt, termOpt, file, printlevel, ELE, NOD, constr, filtr, rfiltr, Vars2Filtr, VarElems, xi);
fid = fopen([file '.txt'],'w');
fprintf(fid,'Optimized compliance = %0.4f [1/MPa]\n', fobj);
fprintf(fid,'Volume fraction = %0.2f\n', VFmean);
fprintf(fid,'Total time = %0.0f [s]\n', optTime);
fprintf(fid,'Number of iterations = %d\n', info.iter);
fclose(fid);
