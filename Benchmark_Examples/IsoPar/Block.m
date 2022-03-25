function [xf, fobj, info, VF, VFmean, optTime] = Block(varargin)
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
LX = 2e3;
LY = 2e3;
LZ = 2e3;
NX = 20;
NY = 20;
NZ = 20;
EleType = 'HEXA8';
[XYZ, CONEC] = GenPrismaticMesh(LX, LY, LZ, NX, NY, NZ, EleType);
nele = size(CONEC, 1);
nnod = size(XYZ, 1);


tol = 1e-3;


% Boundary conditions
esquina = find(XYZ(:,1)>(LX-tol-LX/NX) & XYZ(:,2)<(tol+LY/NY) & XYZ(:,3)<tol);
carader = find(XYZ(:,2)>LY-tol);
carapos = find(XYZ(:,1)<tol);
MFIXDOF = [3*esquina(:); 3*carader(:) - 1; 3*carapos(:) - 2];
MFIXVAL = MFIXDOF*0;


% Loads
P = -100/4;
Fext = zeros(3*nnod,1);
centro = find(XYZ(:,3)<tol & XYZ(:,1)<(LX/NX+tol) & XYZ(:,2)>(LY-LY/NY-tol));
Fext(3*centro) = 1;
centro = find(XYZ(:,3)<tol & XYZ(:,1)<(LX/NX+tol) & XYZ(:,2)>(LY-tol));
Fext(3*centro) = 0.5;
centro = find(XYZ(:,3)<tol & XYZ(:,1)<(tol) & XYZ(:,2)>(LY-LY/NY-tol));
Fext(3*centro) = 0.5;
centro = find(XYZ(:,3)<tol & XYZ(:,1)<(tol) & XYZ(:,2)>(LY-tol));
Fext(3*centro) = 0.25;
Fext = Fext*(P/sum(Fext));


% Material
mater = struct('ElasticModuli',@(x) DEMOD(x),'ConstitutiveLaw', 'LinearElastic');


% Nodes and elements
NOD = struct('Coordinates',XYZ,'ExternalForce',Fext,'UFixedDOFs',MFIXDOF,'UFixedValues',MFIXVAL);
V = ones(size(CONEC,1),1);
ELE = struct('Type',EleType,'Connectivity',CONEC,'Material',mater,'Volumes',V);


%% Coefficients of parametrization
load('NCOEFFS2.mat','NCOEFFS')


%% Constraints
% Bounds
constr.xmin = zeros(7*nele,1);
constr.xmin(1:7:end) = -1;
constr.xmin(2:7:end) = -1;
constr.xmin(3:7:end) = -1;
constr.xmin(4:7:end) = 0.6;
constr.xmin(5:7:end) = 0;
constr.xmin(6:7:end) = 0;
constr.xmin(7:7:end) = 0;
constr.xmax = zeros(7*nele,1);
constr.xmax(1:7:end) = 1;
constr.xmax(2:7:end) = 1;
constr.xmax(3:7:end) = 1;
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
xi(1:7:end) = 0.9999;
xi(2:7:end) = 0.9999;
xi(3:7:end) = 0.9999;
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
file = ['Block_NEle',num2str(nele)];
if constr.AvrgDens > 0
    file = [file,'_VF',num2str(constr.AvrgDens*100)];
end
if constr.MaxAvrgDens > 0
    file = [file,'_MaxVF',num2str(constr.MaxAvrgDens*100)];
end
file = [file,'_MinVFEle5_MaxVFEle45'];
if filtr
    file = [file,'_kf',num2str(kfiltr,'%0.1f')];
end
file = [file,'_IPOpt'];
save([file '.mat'],'ELE','NOD','termOpt','constr','xi');

% Run
VarElems = 1:nele;
opt = 1;
[xf, fobj, info, VF, VFmean, optTime] = OptMainPar2(opt, termOpt, file, printlevel, ELE, NOD, NCOEFFS, constr, filtr, rfiltr, Vars2Filtr, VarElems, xi);
fid = fopen([file '.txt'],'w');
fprintf(fid,'Optimized compliance = %0.4f [1/MPa]\n', fobj);
fprintf(fid,'Volume fraction = %0.2f\n', VFmean);
fprintf(fid,'Total time = %0.0f [s]\n', optTime);
fprintf(fid,'Number of iterations = %d\n', info.iter);
fclose(fid);
