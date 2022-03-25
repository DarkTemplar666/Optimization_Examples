function [xf, fobj, info, VF, VFmean, optTime] = LShape(varargin)
%-------------------------------------------------------------------------%
% inputs = filtr,kfiltr,VF,MaxVF
% filtr = 0 No filter
%         1 Density filter
% kfiltr = filter size parameter (filter radius = kfiltr*(element size))
% VF = Global volume fraction
% MaxVF = Maximum global volume fraction
%-------------------------------------------------------------------------%
%% FEM
% Geometry
LX = 6e3;
LY = 2e3;
LZ = 6e3;
NX = 60;
NY = 20;
NZ = 60;
model = mphload('Lshape.mph');
[~,meshdata] = mphmeshstats(model);
nodes = meshdata.vertex';
XYZ = nodes*1e3;
elements = double(meshdata.elem{2}+1)';
CONEC(:,1) = elements(:,7);
CONEC(:,2) = elements(:,3);
CONEC(:,3) = elements(:,4);
CONEC(:,4) = elements(:,8);
CONEC(:,5) = elements(:,5);
CONEC(:,6) = elements(:,1);
CONEC(:,7) = elements(:,2);
CONEC(:,8) = elements(:,6);
nele = size(CONEC, 1);
nnod = size(XYZ, 1);


tol = 1e-3;


% Boundary conditions
clampnod = find(XYZ(:,3)>LZ-tol);
MFIXDOF = [3*clampnod(:)-2; 3*clampnod(:)-1; 3*clampnod(:)];
MFIXVAL = MFIXDOF*0;


% Loads
P = -100;
Fext = zeros(3*nnod,1);
pnod = find(XYZ(:,1)>LX-tol & XYZ(:,2)>LY/2-LY/NY-tol & XYZ(:,2)<LY/2+LY/NY+tol & XYZ(:,3)>LY/2-LY/NY-tol & XYZ(:,3)<LY/2+LY/NY+tol);
pdof = 3*pnod;
Fext(pdof) = 1;
Fext = Fext*(P/sum(Fext));


% Material
mater = struct('ElasticModuli',@(x) DEMOD(x),'ConstitutiveLaw', 'LinearElastic');


% Nodes and elements
NOD = struct('Coordinates',XYZ,'ExternalForce',Fext,'UFixedDOFs',MFIXDOF,'UFixedValues',MFIXVAL);
V = ones(size(CONEC,1),1);
EleType = 'HEXA8';
ELE = struct('Type',EleType,'Connectivity',CONEC,'Material',mater,'Volumes',V);


%% Constraints
% Bounds
constr.xmin = zeros(7*nele,1);
constr.xmin(1:7:end) = 0.1375;
constr.xmin(2:7:end) = 0.1375;
constr.xmin(3:7:end) = 0.1375;
constr.xmin(4:7:end) = 0.6;
constr.xmin(5:7:end) = 0;
constr.xmin(6:7:end) = 0;
constr.xmin(7:7:end) = 0;
constr.xmax = zeros(7*nele,1);
constr.xmax(1:7:end) = 0.1375;
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
eleSize = LX/NX;
rfiltr = kfiltr*eleSize;
Vars2Filtr = 1:7;


%% Optimization
% Seed
xi = zeros(7*nele,1);
xi(1:7:end) = 0.1375;
xi(2:7:end) = 0.95*0.9999;
xi(3:7:end) = 0.95*0.9999;
xi(4:7:end) = 1.4*0.9999;
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
file = ['LShape_NEle',num2str(nele)];
if constr.AvrgDens > 0
    file = [file,'_VF',num2str(constr.AvrgDens*100)];
end
if constr.MaxAvrgDens > 0
    file = [file,'_MaxVF',num2str(constr.MaxAvrgDens*100)];
end
file = [file,'_MinVFEle6_MaxVFEle41'];
if filtr
    file = [file,'_kf',num2str(kfiltr,'%0.1f')];
end
file = [file,'_IPOpt'];
save([file '.mat'],'ELE','NOD','termOpt','constr','xi');

% Run
VarElems = 1:nele;
opt = 1;
[xf, fobj, info, VF, VFmean, optTime] = OptMainFixTc(opt, termOpt, file, printlevel, ELE, NOD, constr, filtr, rfiltr, Vars2Filtr, VarElems, xi);
fid = fopen([file '.txt'],'w');
fprintf(fid,'Optimized compliance = %0.4f [1/MPa]\n', fobj);
fprintf(fid,'Volume fraction = %0.2f\n', VFmean);
fprintf(fid,'Total time = %0.0f [s]\n', optTime);
fprintf(fid,'Number of iterations = %d\n', info.iter);
fclose(fid);
