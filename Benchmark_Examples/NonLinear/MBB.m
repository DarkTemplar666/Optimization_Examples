function [xf, fobj, info, VF, VFmean, optTime] = MBB(varargin)
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
LX = 6e3;
LY = 1e3;
LZ = 2e3;
NX = 60;
NY = 10;
NZ = 20;
EleType = 'HEXA8';
[XYZ, CONEC] = GenPrismaticMesh(LX, LY, LZ, NX, NY, NZ, EleType);
nele = size(CONEC, 1);
nnod = size(XYZ, 1);


tol = 1e-3;


% Boundary conditions
faceXf = find(XYZ(:,1)>(LX-tol));
faceYf = find(XYZ(:,2)>(LY-tol));
corner1 = find(XYZ(:,1)<(tol+LX/NX) & XYZ(:,2)<(tol+LY/NY) & XYZ(:,3)<tol);
MFIXDOF = [3*faceXf(:)-2; 3*faceYf(:)-1; 3*corner1(:)];
MFIXVAL = MFIXDOF*0;


% Loads
P = -100/4;
Fext = zeros(3*nnod,1);
corner2 = find(XYZ(:,1)>(LX-tol-LX/NX) & XYZ(:,2)>(LY-tol-LY/NY) & XYZ(:,3)>(LZ-tol));
Fext(3*corner2) = 1;
corner2 = find(XYZ(:,1)>(LX-tol) & XYZ(:,2)>(LY-tol-LY/NY) & XYZ(:,3)>(LZ-tol));
Fext(3*corner2) = 0.5;
corner2 = find(XYZ(:,1)>(LX-tol-LX/NX) & XYZ(:,2)>(LY-tol) & XYZ(:,3)>(LZ-tol));
Fext(3*corner2) = 0.5;
corner2 = find(XYZ(:,1)>(LX-tol) & XYZ(:,2)>(LY-tol) & XYZ(:,3)>(LZ-tol));
Fext(3*corner2) = 0.25;
Fext = Fext*(P/sum(Fext));


% Material
mater = struct('ElasticModuli',@(x) DEMOD(x),'ConstitutiveLaw', 'LinearElastic');


% Nodes and elements
NOD = struct('Coordinates',XYZ,'ExternalForce',Fext,'UFixedDOFs',MFIXDOF,'UFixedValues',MFIXVAL);
V = ones(size(CONEC,1),1);
ELE = struct('Type',EleType,'Connectivity',CONEC,'Material',mater,'Volumes',V);


%% Constraints
% Bounds
constr.xmin = zeros(7*nele,1);
constr.xmin(1:7:end) = 0.07;
constr.xmin(2:7:end) = 0.07;
constr.xmin(3:7:end) = 0.07;
constr.xmin(4:7:end) = 0.6;
constr.xmin(5:7:end) = 0;
constr.xmin(6:7:end) = 0;
constr.xmin(7:7:end) = 0;
constr.xmax = zeros(7*nele,1);
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
eleSize = LX/NX;
rfiltr = kfiltr*eleSize;
Vars2Filtr = 1:7;


%% Optimization
% Seed
xi = zeros(7*nele,1);
xi(1:7:end) = 0.95*0.9999;
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
file = ['MBB2_NEle',num2str(nele)];
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
VarElems = 1:nele;
opt = 1;
[xf, fobj, info, VF, VFmean, optTime] = OptMain(opt, termOpt, file, printlevel, ELE, NOD, constr, filtr, rfiltr, Vars2Filtr, VarElems, xi);
fid = fopen([file '.txt'],'w');
fprintf(fid,'Optimized compliance = %0.4f [1/MPa]\n', fobj);
fprintf(fid,'Volume fraction = %0.2f\n', VFmean);
fprintf(fid,'Total time = %0.0f [s]\n', optTime);
fprintf(fid,'Number of iterations = %d\n', info.iter);
fclose(fid);
