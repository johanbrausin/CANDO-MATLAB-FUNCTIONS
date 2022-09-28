function [A,rhoj,T,rho,eff,avg_TE,K,J, nR] = FRETopt(D,n,R,t_end,R0,tau,mol_ext,kind,ss,random,pos_don,pos_acc,par)

%===================================================================================================================================
% AUTHOR : E. BOULAIS, LCBB, JULY 2013
% DESCRIPTION : Calculates energy transfer in the FRET limit among an
% assembly of dyes. Function has been designed to take advantage of Matlab
% functionality and avoid double for-loops which are usually computationally inefficient.
% Computation time for 10,000 dyes is approximately one hour using 13 GB RAM. See Supporting Information for additional solution times.

%nt==> Number of types of dyes
%%       --Symbol--    --Name--      --Dimension--      --Units--          --Description
%INPUT :    D        Dipole matrix       n*3              N.A.          Norm. dipole orientation
%           n        Number of dyes       1               N.A.          Total number of dyes
%           R        Dye position        n*3              [m]           Dipole center position
%          t_end     Time of sim.         1               [s]           Total time of simulation
%           R0       Foerster distance   nt*nt            [m]           Matrix of Foerster radius
%                                                                       R0(i,j) is R0 for i-->j
%           tau      Lifetime            1*nt             [s]           Fluor. Lifetime
%         mol_ext    Molar extinction    1*nt          [M^1*cm^-1]      Molar extinction
%          kind      Type of dye         1*n              N.A.          Type of each dye
%           ss       Steady-state flag    1               N.A.          If ss=1, steady-state calculation only  
%          random    Random flag         1*x              N.A.          Position of dyes with random dip. orient.
%          pos_don   Donor position      1*x              N.A.          Index of dyes that can be considered as donors
%          pos_acc   Acceptor position   1*x              N.A.          Index of dyes that can be considered as acceptors
%           par      Parallel flag        1               N.A.          If par==1, then explicit parallelisation of matrix assembly
%
%OUTPUT    avg_TE    Transm. efficiency   1               N.A.          Transmission efficiency
%          eff       Donor(s) efficiency  1               N.A.          Quenching efficiency of donors 
%          rho       Exciton population  n*1(ss)/n*ts(td) N.A.         Exciton population of each dye in the steady-state (ss) or
%                                                                       time-dependent(td)
%            T       Time steps           ts               s            Time-steps of simulation
%          rhoj      Initial condition   n*1              N.A.          Initial exciton population
%            A       Transition matrix   n*n              s^-1          Total transition matrix
%DEPENDANCE : odeexcCL.m
%===================================================================================================================================
%%
% Turn on explicit parallelisation of assembly on flag par
if par==1
    matlabpool open  %before it wae ####parpool local###
end


%==========================================================================
%% Generation of the transition matrix J
%==========================================================================
%Interaction matrix based on the dipole-dipole approximation

%==========================================================================
% Generation of matrix Rd, calculating the dye difference of positions
%==========================================================================
%construct matrix R

Rm1=zeros(3*n,n);

Rt=R';
rowIdx=(1:size(Rt,1))'; % this part of code thos the same as Rm2=repmat(R',n,1)
colIdx=(1:size(Rt,2))';   % but is faster
Rm2=Rt(rowIdx(:, ones(n,1)), colIdx(:, ones(1,1)));
nR=zeros(n,n);
for i=1:n
    RtR=R(i,:)';
    rowIdxR=(1:size(RtR,1))'; %this part of code does the same as 
    colIdxR=(1:size(RtR,2))'; %Rm1(1+3*(i-1):3*i,:)=repmat(R(i,:)',1,n);
    Rm1buf=RtR(rowIdxR(:, ones(1,1)), colIdxR(:, ones(n,1)));
    Rm1(1+3*(i-1):3*i,:)=Rm1buf;
    
    
end
% Rd is a matrix containing Rj-Ri in the form   dx11 dx12 dx13...dx1n
%                                                dy11 dy12 dy13 ... dy1n
%                                                dz11 dz12 dz13 ... dz1n
%                                                dx21 dx22 dx23 ... dx2n
%                                                dy21 dy22 dy23 ... dy1n
%                                                dz21 dz22 dz23 ... dz1n
%                                               ...  ...  ...  ... ...
%                                                dzn1 dzn2 dzn3 ... dznn
Rd=zeros(3*n,n); 
Rd=Rm2-Rm1;
Rm1=[];
Rm2=[];
mult=ones(n,n);
parfor i=1:n
    %calculate the norm of Rd
    nR(i,:)=max(5e-10,sqrt(Rd(1+3*(i-1),:).^2+Rd(2+3*(i-1),:).^2+Rd(3+3*(i-1),:).^2));
    %Calculate (Di.Rij)(Dj.Rij)
    buff=ones(n,n);
    buff(i,:)=D(i,:)*Rd(3*(i-1)+1:3*i,:);
    buff(:,i)=(D(i,:)*reshape(Rd(:,i),3,n))';
    mult=mult.*buff;

end
%Clear useless matrix
buff=[];
J=zeros(n,n); 
K=zeros(n,n);
%Interaction matrix
J=D*D'./nR.^3-3.*mult./nR.^5;

%Fix interaction with randomly oriented dyes

if isempty(random)==0 %if some dyes are randomly oriented 
%Fix interaction with randomly oriented dyes
%random-random ==> kappa^2=2/3
%random-fixed ==> kappa^2=1/3+(D.R)^2
    for j=1:length(random)
        for k=1:n
            if (isempty(find(random==k,1))==1) %dye k is not random
                 Dpr=D(k,:)*Rd(3*(k-1)+1:3*k,random(j));
                 Dpr=Dpr/nR(random(j),k);
                 J(random(j),k)=sqrt(1/3+Dpr^2)./nR(random(j),k).^3;
                 J(k,random(j))=J(random(j),k);
            else % dye k is also random
                 J(random(j),k)=sqrt(2/3)./nR(random(j),k).^3;
                 J(k,random(j))=J(random(j),k);
             end
         end
              
     end
 end

tauK=tau(kind);
rowIdxtau=(1:size(tauK,1))'; % Same as repmat(tauK,n,1) but faster
colIdxtau=(1:size(tauK,2))';
tauKr=tauK(rowIdxtau(:, ones(n,1)), colIdxtau(:, ones(1,1)));

%Transition rate matrix
K=R0(kind,kind)'.^6.*J.*conj(J)./(tauKr)*1.5;

%Fix infinite self-interaction
for i=1:n
    K(i,i)=0;
end
%K=min(K,1e10);
%max(max(K))
%Clear useless matrix
%J=[];
%nR=[];
Rd=[];
mult=[];

%==========================================================================
%% Generation of the Transition Matrix (K in the main manuscript)
%==========================================================================

    DI=-sum(K)-1./tauK;
    A=K+diag(DI,0);

    %======================================================================    
%% Solving the system
%==========================================================================
    %#ode45 is inefficient for stiff problems. ode23t is efficient and 
    %prevents artificial numerical damping.
if ss~=1
   
    %time-dependant resolution
    tsave=[0 t_end];
    rhoinit=zeros(1,n);
    rhoinit(pos_don)=mol_ext(kind(pos_don));
    rhoinit=rhoinit./sum(rhoinit);
    rhoj=rhoinit;
    options = odeset('InitialStep',1e-15);

    [T,rho]=ode23t(@(t,y)odeexcCL(t,y,A),tsave,rhoinit,options);

    avg_TEbuf=zeros(size(pos_acc));

    for i=1:length(pos_acc)
        avg_TEbuf(i)=trapz(T,rho(:,pos_acc(i))/tauK(pos_acc(i)));
    end
    
    avg_TE=sum(avg_TEbuf);
    eff_buf=zeros(size(pos_don));

    for i=1:length(pos_don)
        eff_buf(i)=trapz(T,rho(:,pos_don(i))/tauK(pos_don(i)));
    end
    eff=1-sum(eff_buf);

    
else
    %Steady-state resolution, assuming constant input at donors
    T=0;
    rhoj=zeros(n,1);
    rhoj(pos_don)=mol_ext(kind(pos_don));
    rhoj=rhoj./sum(rhoj);
    rho=A\(-rhoj);
    avg_TEbuf=zeros(size(pos_acc));
    eff_buf=zeros(size(pos_don));
    for i=1:length(pos_acc)
        avg_TEbuf(i)=rho(pos_acc(i))/tauK(pos_acc(i));
    end
    avg_TE=sum(avg_TEbuf);
    for i=1:length(pos_don)
        eff_buf(i)=rho(pos_don(i))/tauK(pos_don(i));
    end
    eff=1-sum(eff_buf);
end

  
if par==1
    matlabpool close   %before it was ####matlabpool###
end

end

