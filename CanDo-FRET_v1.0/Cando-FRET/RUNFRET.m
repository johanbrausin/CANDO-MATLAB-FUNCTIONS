clear
clc

n = input('please enter the number of molecules in the system->');
nt = input('enter the types of different dyes in the system->');
kind = [1:n];
ss = 0 ;
R0 = zeros(n,n);
for i=1:n
    for j=1:n
        R0(i,j)= input('please enter the Forster Radius of the molecule i with molecule j and check out the matrix (each value is multiplied by 10^-10 to have meters scale)->')*10^-10
    end
end
R = zeros(n,3);
for k=1:n
   R(k,:)= input('enter each dye coordinates as a vector [x,y,z] and check the matrix each value is multiplied by 10^-10 to have meters scale) ->')*10^-10 
end
pos_don = input('enter which kind of mol are DONORS depending the row of the previous matrix R  as an array [x1 x2 ...]->');
pos_acc = input('enter which kind of mol are ACCEPTORS depending the row of the previous matrix R as an array [x1 x2 ...]->');
a = input('are you considering the transition dipole of each dye? (1/0)->');
switch a
    case 1
        D = zeros(n,3);
        for z=1:n
            D(z,:)= input('if you are considering the transition dipole matrix input each unitary vector as [x y z] and check the matrix->')
            nr = input('eneter the number of molecules with random dipole orientation ->');
            random = zeros(1,nt);
            random(z)= input('enter which dyes are random (1,2,...)->')
        end 
    case 0
        D = [0 0 1 ;0 1 0; 1 0 0; 0 0 1];
        random = [1:n];
end 

t_end = input('enter the time of simulation (this will be multiplied by 10^-9, to have seconds scale)->')*10^-9;

%tau = [3.1 3.4 3.2 0.7 ].*10^-9;
tau = zeros(1,nt);
for s=1:nt
    
    tau(s) = input ('fluorescence life time of the dyes multiplied by 10^-9 (depends what number you gave each type of dye )->')*10^-9
end
mol_ext = zeros(1,nt);
for t=1:nt
    mol_ext(t) = input ('molar extinction in [M^1*cm^-1] unities ->')
end
par = 0; %%its better having 0 meanwhile
[A,rhoj,T,rho,eff,avg_TE, K,J, nR] = FRETopt(D,n,R,t_end,R0,tau,mol_ext,kind,ss,random,pos_don,pos_acc,par);


%%TEST VALUES 
%%_________________________________________________________________________________________________________________________
% par = 0; %%its better having 0 meanwhile
% ss = 0 ;
% n=4;
% t_end = 10*10^-9;
% kind=[1 2 3 4];
% random =[1 2 3 4];
% pos_don =[2];
% pos_acc = [3 4];
% R = [210.628 206.962 -1.370 ;148.366 144.530 -2.574; 103.563 98.914 6.607; 47.802 43.959 -2.427].*10^-10;
% R0 = [0 62 51 45; 26 58 67 63; 19 10 0 79; 45 63 79 0].*10^-10;
% D = [0 0 10 ;0 1 0; 1 0 0; 0 0 0];
% tau = [3.1 3.4 3.2 0.7].*10^-9;
% random = [1 2 3 4];
% mol_ext = [90000 120000 150000 240000];
% [A,rhoj,T,rho,eff,avg_TE, K,J, nR] = FRETopt(D,n,R,t_end,R0,tau,mol_ext,kind,ss,random,pos_don,pos_acc,par);

%%________________________________________________________________________________________________

disp('transition matrix')
disp(A)
disp('Donors efficency')
disp(eff)
disp('Transmission efficency')
disp(avg_TE)

%x1 = tau*abs(A);
%X = x1.^(-1/6);
%f1 = K*tau';
%f2 = f1.^(-1/6);
%phitden = 1+(f2.^6);
%phit = phitden.^-1;
%figure
%scatter(f2,phit,'r')
%xlabel('$\frac{r}{R_0}$','Interpreter','latex','FontSize',15)
%ylabel('efficency transfer','Interpreter','latex','FontSize',15)
%title('Variations in the transfer efficency','Interpreter','latex','FontSize',15)
%grid
figure
imagesc(A);
title('Transition Matrix','Interpreter','latex','FontSize',15)
colorbar
figure 
mesh(rho);
title('Exciton Population for each Dye Surface graph','Interpreter','latex','FontSize',15)
colorbar
figure 
imagesc(rho);
title('Exciton Population for each Dye Matrix graph','Interpreter','latex','FontSize',15)
colorbar
grid
%figure
% plot(T,rho(:,3))
% xlabel('$Time(sec)$','Interpreter','latex','FontSize',15)
% ylabel('Exciton population Dye n','Interpreter','latex','FontSize',15)
% title('Variations in the transfer efficency','Interpreter','latex','FontSize',15)
% grid
for i = 1:n
   figure
   plot(T,rho(:,i))
   xlabel('$Time(sec)$','Interpreter','latex','FontSize',15)
   ylabel('Exciton population Dye n','Interpreter','latex','FontSize',15)
   title('Variations in the transfer efficency','Interpreter','latex','FontSize',15)
   grid
end

% s = size(K);
% X = zeros(s(1),s(2)); 
% for i = 1:length(tau)
%     for j = 1:s(2)
%         X(i,j) = (K(i,j)*tau(i))^(-1/6);  
%     end
% end
%BOYD NON LINEAR OPTIC BOOK

