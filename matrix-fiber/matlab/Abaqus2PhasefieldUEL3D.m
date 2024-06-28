function Abaqus2PhasefieldUEL3D(inputPath,MatProp,style,enable_reduced)
% style = single or double; 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!------------------ jing Lee----------!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!------------------2024/06/20-----------!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!! --------- BETA Version --------- !!!!!!!!!!!!!!!!!!!!
%  When this script is used, it is the user's responsibility to check and
%                         interpret all results.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
%     The matlab script is modifeid according to molnar,
%
%     More on the theory and implementation of the UEL to model diffuse 
%     fracture with phase-field can be found in our recent paper:
%       G. Moln, A. Gravouil, 2D and 3D Abaqus implementation of a
%        robust staggered phase-field solution for modeling brittle
%        fracture, Finite Elements in Analysis and Design, 130 pp. 27-38, 2017.
%
%     https://www.sciencedirect.com/science/article/pii/S0168874X16304954
%
%    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     The script is not able to convert all types of input files, only the
%     ones created the right order discussed in these tutorials (e.g. the
%     mesh needs to be defined on the parts):
%        http://www.molnar-research.com/tutorials.html
%    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%               
%
%  -----  INPUT variables are the following -------------------------------
% inputPath - is a string which point to the input file created by ABAQUS,
%				if left blank [], the script opens an inquiry
%
% MatProp - is an array of (x,22) values containing materials parameters
%           both for the stress and phase field elements, where x is the
%           number of different materials used:
%
%  PROPS(1) = Density (rho)
%  PROPS(2) = Longitudinal Young¡¯s modulus (E1) 
%  PROPS(3) = Transverse Young¡¯s modulus (E2,E3)
%  PROPS(4) = Shear modulus (G12,G13)
%  PROPS(5) = Shear modulus (G23)
%  PROPS(6) = Poisson¡¯s ratio (V12,V13)
%  PROPS(7) = Poisson¡¯s ratio (V23)
%  PROPS(8) = Transverse tensile strength (ST1) 
%  PROPS(9) = Transverse shear strength (ST2)
%  PROPS(10) = In-plane shear strength (STL) 
%  PROPS(11) = Longitudinal tensile strength(YL)
%  PROPS(12) = Transverse tensile critical energy release rate(GcT1)
%  PROPS(13) = Transverse shear critical energy release rate (GcT2)
%  PROPS(14) = In-plane shear critical energy release rate (GcTL)
%  PROPS(15) = Longitudinal tensile critical energy release rate(GcL)
%  PROPS(16) = matrix phase field internal length scales (l0m)
%  PROPS(17) = fiber phase field internal length scales(l0L)
%  PROPS(18) = artificial viscosity parameters for matrix(etam)
%  PROPS(19) = artificial viscosity parameters for fiber(etaL)
%  PROPS(20) =Material direction vector component (nx)
%  PROPS(21) =Material direction vector component (ny)
%  PROPS(22) =Material direction vector component (nz)
% -------------------------------------------------------------------------
%
% BASIC structure of the script:
% 0) Loading input variables
% 1) Start copying original input file into the new one
% 2) When the beginning of the element definition section is found
%    (*Element) it stops copying and starts uploading elements into ElMat
%    array.
% 3) Parallel it writes the first layer of finite elements
% 4) When ready it replicates the two additional layers
% 5) Creates element sets and assigns material properties
% 6) Continues copying the original input file until the end of the assembly
%	 part. Here the script adds an additional element set summing up all
%    the UMAT elements.
% 7) Then it finishes the with the simulation steps, adding always a field
%    output request for the SDVs
% -------------------------------------------------------------------------
%
%node offset
offset = 0.01;
narginLoc=nargin;

if narginLoc<2
    if narginLoc==0
% ---------------- Asking for the input file -------------------
        [FileName,PathName,FilterIndex] = uigetfile('*.inp');
        inputPath=[PathName FileName];
    else
    end
    if isa(inputPath,'char')~=1
        error('Provide a valid input file.')
    end
% ---------------- Asking for material properties ---------------
    nMatini=input('Number of materials? - ');
        for i=1:nMatini
            MatProp(i,1)=input(['Density (rho)  (' num2str(i) ')? - ']);
            MatProp(i,2)=input(['Longitudinal Young¡¯s modulus (E1) (' num2str(i) ')? - ']);
            MatProp(i,3)=input(['Transverse Young¡¯s modulus (E2,E3)  (' num2str(i) ')? - ']);
            MatProp(i,4)=input(['Shear modulus (G12,G13) (' num2str(i) ')? - ']);
            MatProp(i,5)=input(['Shear modulus (G23) (' num2str(i) ')? - ']);
            MatProp(i,6)=input(['Poisson¡¯s ratio (V12,V13) (' num2str(i) ')? - ']);
            MatProp(i,7)=input(['Poisson¡¯s ratio (V23) (' num2str(i) ')? - ']);
            MatProp(i,8)=input(['Transverse tensile strength (ST1)  (' num2str(i) ')? - ']);
            MatProp(i,9)=input(['Transverse shear strength (ST2) (' num2str(i) ')? - ']);
            MatProp(i,10)=input(['In-plane shear strength (STL)  (' num2str(i) ')? - ']);
            MatProp(i,11)=input(['Longitudinal tensile strength(YL) (' num2str(i) ')? - ']);
            MatProp(i,12)=input(['Transverse tensile critical energy release rate(GcT1) (' num2str(i) ')? - ']);
                       
            MatProp(i,13)=input(['Transverse shear critical energy release rate (GcT2) (' num2str(i) ')? - ']);
            MatProp(i,14)=input(['In-plane shear critical energy release rate (GcTL) (' num2str(i) ')? - ']);
            MatProp(i,15)=input(['Longitudinal tensile critical energy release rate(GcL) (' num2str(i) ')? - ']);
            MatProp(i,16)=input(['matrix phase field internal length scales (l0m) (' num2str(i) ')? - ']);
            MatProp(i,17)=input(['fiber phase field internal length scales(l0L) (' num2str(i) ')? - ']);
            MatProp(i,18)=input(['artificial viscosity parameters for matrix(etam) (' num2str(i) ')? - ']);
            MatProp(i,19)=input(['artificial viscosity parameters for fiber(etaL) (' num2str(i) ')? - ']);
            MatProp(i,20)=input(['Material direction vector component (nx) (' num2str(i) ')? - ']);
            MatProp(i,21)=input(['Material direction vector component (ny) (' num2str(i) ')? - ']);
            MatProp(i,22)=input(['Material direction vector component (nz) (' num2str(i) ')? - ']);
        end
end

if length(MatProp(1,:))<22
    error('Not enought material properties are given.')
end
if length(MatProp(1,:))>22
    error('Too much material properties are given.')
end
if isempty(inputPath)==1
    [FileName,PathName,~] = uigetfile('*.inp');
    inputPath=[PathName FileName];
end

% ---------------- Uploaing material properties -------------------
rho=MatProp(:,1);
E1=MatProp(:,2);
E2=MatProp(:,3);
G12=MatProp(:,4);
G23=MatProp(:,5);
V12=MatProp(:,6);
V23=MatProp(:,7);
ST1=MatProp(:,8);
ST2=MatProp(:,9);
STL=MatProp(:,10);
SYL=MatProp(:,11);
GcT1=MatProp(:,12);
GcT2=MatProp(:,13);
GcTL=MatProp(:,14);
GcYL=MatProp(:,15);
l0m=MatProp(:,16);
l0L=MatProp(:,17);
etam=MatProp(:,18);
etaL=MatProp(:,19);
nx=MatProp(:,20);
ny=MatProp(:,21);
nz=MatProp(:,22);



% ---------------- Initializing output file -------------------
inputName=inputPath;
output=[inputPath(1:end-4) '_UEL.inp'];

fid=fopen(inputName,'rt');
fout1 = fopen(output, 'wt');
a=fgets(fid);

cnt=0;

% ---------------- Sart reading input file -------------------
cntE=1;
while(ischar(a))

% ------ search nodes ------
while(ischar(a))
    if strfind(a,'*Node')~=0
        break
    end
    fprintf(fout1,a);    
    a=fgets(fid);
end
a=fgets(fid);
% ------ save nodes ------
cnt=0;
while(ischar(a))
    if strfind(a,'*Element')~=0
        break
    end
    cnt=cnt+1;
    Nodes(cnt,:)=str2num(a);    
    a=fgets(fid);
end
% ------ Repeating nodes ------
fprintf(fout1,['*Node\n']);
nnodes = size(Nodes,1);
for j=1:nnodes
    fprintf(fout1,[num2str(j) ', ' num2str(Nodes(j,2:end-1),'%d, ') ' ' num2str(Nodes(j,end)) '\n']);
end
if strcmp(style,"double")
% an extra layer for the second phase
% the offset is for the z coordiated, should be choosed according the model
% geometry
    for j=1:nnodes
    fprintf(fout1,[num2str(nnodes+j) ', ' num2str(Nodes(j,2:end-1)-offset,'%d, ') ' ' num2str(Nodes(j,end)-offset) '\n']);
    end
end
a=fgets(fid);
% ------ Defintion of UEL types ------
if length(str2num(a))==5
   % triangular elements 
    fprintf(fout1,['*********************** TETRAHEDRAL ****************************\n']);
    fprintf(fout1,['*User element, nodes=4, type=VU1, properties=22, coordinates=3, VARIABLES=26\n']);
    fprintf(fout1,['1,2,3,11\n']);
    NnodeE(cntE)=4;
else
   % rectangular elements 
    fprintf(fout1,['*********************** BRICK ***************************\n']);
    fprintf(fout1,['*User element, nodes=8, type=VU2, properties=22, coordinates=3, VARIABLES=208\n']);
    fprintf(fout1,['1,2,3,11\n']);
    NnodeE(cntE)=8;
end

% ------ Uploading original mesh ------
cnt=0;
while isempty(strfind(a,'*Nset,')) && isempty(strfind(a,'*Element,'))
    cnt=cnt+1;
    ElMat{cntE}(cnt,:)=str2num(a);
    a=fgets(fid);
end
nElem(cntE)=length(ElMat{cntE}(:,1));

fprintf(fout1,['***************************************************************\n']);

if isempty(strfind(a,'*Element,'))==0
    cntE=cntE+1;
else
    break
end
end

% ------ All elements ------
nElemAll=sum(nElem);

% ------ Replicating UELs ------
% ------ for double-phase field modelling, an extra layer of cells need to be replicated 
for j=1:cntE
    if NnodeE(j)==4
        ej=1;
    elseif NnodeE(j)==8
        ej=2;
    end
    fprintf(fout1,['***************************************************************\n']);
    fprintf(fout1,['*Element, type=VU' num2str(ej) '\n']);
    for i=1:nElem(j)
        fprintf(fout1,[num2str(ElMat{j}(i,1)) ', ' num2str(ElMat{j}(i,2:end-1),'%d, ') ' ' num2str(ElMat{j}(i,end)) '\n']);
    end
    
    if strcmp(style,"double")
        fprintf(fout1,['***************************************************************\n']);
        fprintf(fout1,['*Element, type=VU' num2str(ej) '\n']);
        for i=1:nElem(j)
            fprintf(fout1,[num2str(ElMat{j}(i,1)+2*nElemAll) ', ' num2str(ElMat{j}(i,2:end-1)+nnodes,'%d, ') ' ' num2str(ElMat{j}(i,end)+nnodes) '\n']);
        end
    end
end

% ------ Sets and material properties ------
fprintf(fout1,['***************************************************************\n']);
fprintf(fout1,['********************** ASSIGNING MATERIAL PROP ****************\n']);

cnt2=1;
while isempty(strfind(a,'** Section'))
   while isempty(strfind(a,'*Elset,'))
       a=fgets(fid);
   end
   bin=strfind(a,'elset=');
   Ename=a(bin+6:end-1);
   if isempty(strfind(Ename,','))==0
      bin=strfind(Ename,',');
      Ename=Ename(1:bin-1);
   end
   ElstName{cnt2}=Ename;
   fprintf(fout1,a);
   ElSeT{1,cnt2}=a;
   a=fgets(fid);
   cnt3=2;
   while isempty(strfind(a,'** Section:')) && isempty(strfind(a,'*Nset,'))
       ElSeT{cnt3,cnt2}=a;
       fprintf(fout1,a);
       a=fgets(fid);
       cnt3=cnt3+1;
   end
   Nline(cnt2)=cnt3-1;
   cnt2=cnt2+1;
end
fprintf(fout1,['***************************************************************\n']);

nMat=length(ElstName);
if length(E1)<nMat
    error('Model contains more materials. Provide more material properties!')
end

for k=1:nMat
    fprintf(fout1,['*Uel property, elset=' ElstName{k} '\n']);
% eight data per line!
        fprintf(fout1,[num2str(rho(k)) ', ' num2str(E1(k)) ', ' num2str(E2(k)) ', ' num2str(G12(k)) ', ' num2str(G23(k)) ', ' num2str(V12(k)) ', ' num2str(V23(k)) ', ' num2str(ST1(k)) ',\n' num2str(ST2(k)) ', ' num2str(STL(k)) ', ' num2str(SYL(k)) ', ' ...
            num2str(GcT1(k)) ', ' num2str(GcT2(k)) ', ' num2str(GcTL(k)) ', ' num2str(GcYL(k)) ', ' num2str(l0m(k)) ',\n' num2str(l0L(k)) ', ' num2str(etam(k)) ', ' num2str(etaL(k)) ', ' num2str(nx(k)) ', ' num2str(ny(k)) ', ' num2str(nz(k)) '\n']);
end
fprintf(fout1,['**rho, E1, E2, G12, G23, V12, V23, ST1, ST2, STL, SYL, GcT1, GcT2, GcTL, GcYL, l0m, l0L, etam, etaL, nx, ny, nz\n']);

fprintf(fout1,['***************************************************************\n']);
if strcmp(style,"double")
for k=1:nMat
    b=ElSeT{1,k};
    c=length(ElstName{k});
    fprintf(fout1,[b(1:14) ElstName{k} '_P2' b(14+c+1:end)]);
    if isempty(strfind(b,'generate'))
        for j=2:Nline(k)-1
            b=str2num(ElSeT{j,k});
            fprintf(fout1,'%d, ',b+2*nElemAll);
            fprintf(fout1,'\n');
        end
        b=str2num(ElSeT{Nline(k),k});
        fprintf(fout1,'%d, ',b(1:end-1)+2*nElemAll);
        fprintf(fout1,'%d',b(end)+2*nElemAll);
        fprintf(fout1,'\n');
        
    else
        b=str2num(ElSeT{2,k});
        fprintf(fout1,[num2str(b(1)+2*nElemAll) ', ' num2str(b(2)+2*nElemAll) ', ' num2str(b(3)) '\n']);
    end
end


fprintf(fout1,['***************************************************************\n']);
for k=1:nMat
    fprintf(fout1,['*Uel property, elset=' ElstName{k} '_P2\n']);
        fprintf(fout1,[num2str(rho(k)) ', ' num2str(E1(k)) ', ' num2str(E2(k)) ', ' num2str(G12(k)) ', ' num2str(G23(k)) ', ' num2str(V12(k)) ', ' num2str(V23(k)) ', ' num2str(ST1(k)) ',\n' num2str(ST2(k)) ', ' num2str(STL(k)) ', ' num2str(SYL(k)) ', ' ...
            num2str(GcT1(k)) ', ' num2str(GcT2(k)) ', ' num2str(GcTL(k)) ', ' num2str(GcYL(k)) ', ' num2str(l0m(k)) ',\n' num2str(l0L(k)) ', ' num2str(etam(k)) ', ' num2str(etaL(k)) ', ' num2str(nx(k)) ', ' num2str(ny(k)) ', ' num2str(nz(k)) '\n']);
end
fprintf(fout1,['**rho, E1, E2, G12, G23, V12, V23, ST1, ST2, STL, SYL, GcT1, GcT2, GcTL, GcYL, l0m, l0L, etam, etaL, nx, ny, nz\n']);
end


% ------ Replicating dummy UMAT elements ------
fprintf(fout1,['***************************************************************\n']);

for j=1:cntE
    if (enable_reduced&&NnodeE(j)>4)
    fprintf(fout1,['*Element, type=C3D' num2str(NnodeE(j)) 'R\n']);
    else
    fprintf(fout1,['*Element, type=C3D' num2str(NnodeE(j)) '\n']);
    end        
    for i=1:nElem(j)
        fprintf(fout1,[num2str(ElMat{j}(i,1)+nElemAll) ', ' num2str(ElMat{j}(i,2:end-1),'%d, ') ' ' num2str(ElMat{j}(i,end)) '\n']);
    end
    fprintf(fout1,['***************************************************************\n']);
end

while isempty(strfind(a,'material='))
    a=fgets(fid);
end

UMATname=a(strfind(a,'material=')+9:end);

fprintf(fout1,['*Elset, elset=umatelem, generate\n']);
fprintf(fout1,[num2str(nElemAll+1) ', ' num2str(nElemAll*2)   ', 1\n']);

fprintf(fout1,['***************************************************************\n']);
fprintf(fout1,['*Solid Section, elset=umatelem, material=' UMATname '1.0\n']);


while isempty(strfind(a,'*End Part'))
    a=fgets(fid);
end

while(ischar(a))
    if strfind(a,'*Instance, name=')~=0
        break
    end
    fprintf(fout1,a);    
    a=fgets(fid);
end

InstanceName=a(strfind(a,'*Instance, name=')+16:strfind(a,', part=')-1);

while(ischar(a))
    if strfind(a,'*End Assembly')~=0
        break
    end
    fprintf(fout1,a);    
    a=fgets(fid);
end

fprintf(fout1,['**\n']);

fprintf(fout1,['*Elset, elset=umatelem, instance=' InstanceName ', generate\n']);
fprintf(fout1,[num2str(nElemAll+1) ', ' num2str(nElemAll*2) ', \n']);

fprintf(fout1,['**\n*End Assembly\n']);
fprintf(fout1,['***************************************************************\n']);

a=fgets(fid);

while(ischar(a))

    while(ischar(a))
        if strfind(a,'*Output, history')~=0
            break
        end
        fprintf(fout1,a);    
        a=fgets(fid);
    end
        fprintf(fout1,['*element output, elset=umatelem\n']);
        fprintf(fout1,['SDV\n']);
    
    fprintf(fout1,['**\n']);
    a=fgets(fid);
    fprintf(fout1,a);
    a=fgets(fid);

end
fclose('all');

% --------------- Start modifying the UEL fortran file -------------------

% determining if the OS if Lunix, Windows of Mac
comp = computer;
if strfind(comp,'WIN')~=0
    comptype=1;     % windows
elseif strfind(comp,'LNX')~=0
    comptype=2;     % linux
elseif strfind(comp,'MAC')~=0
    comptype=3;     % mac
else
    comptype=input('Please provide OS type: 1 - Microsoft Windows; 2 - Linux; 3 - Apple Mac: ');     % asking user
    if comptype==1 || comptype==2 || comptype==3
    else
        comptype=input('Wrong type.\n Please provide OS type: 1 - Microsoft Windows; 2 - Linux; 3 - Apple Mac: ');     % asking user
        if comptype==1 || comptype==2 || comptype==3
        else
        error('Wrong type. Try again.')
    end
    end
end

fullp=mfilename('fullpath');
if isempty(fullp)~=1
    input2=[fullp(1:end-30) '\subroutine\doublephase.for'];
else
    fullp=matlab.desktop.editor.getActiveFilename;
    input2=[fullp(1:end-32) '\subroutine\doublephase.for'];
end

if comptype==1
    output=[inputName(1:end-4) '_UEL.for'];
else
    output=[inputName(1:end-4) '_UEL.f'];
end

fid=fopen(input2,'rt');
fout1 = fopen(output, 'wt');
a=fgets(fid);

while(ischar(a))
    if strfind(a,'N_ELEM=')~=0
        strs=strfind(a,'N_ELEM=');
        ends=strs+6;
        a=[a(1:ends) num2str(nElemAll) a(ends+2:end)];
    end
    if strfind(a,'logical :: enable_c3d8r')~=0
        strs=strfind(a,'enable_c3d8r = ');
        ends=strs+14;        
        if (enable_reduced)
            a=[a(1:ends) '.true.\n'];
        else
            a=[a(1:ends) '.false.\n'];
        end
    end
    fprintf(fout1,a);    
    a=fgets(fid);
end

fclose('all');

return