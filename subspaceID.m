function [A,B,C,D,x0,sv] = subspaceID(u,y,s,n,method)
%% Instructions:
% Implement your subspace ID methods here.
% Avoid duplicate code! 
% Write the method specific code into the switch case!
% Use the following function inputs and outputs.

% Function INPUT 
% u         system input (matrix of size N x m)
% y         system output (matrix of size N x l)
% s         block size (scalar)
% n         model order (scalar)
% method    method (string e.g. 'moesp')
%
% Function OUTPUT
% A         System matrix A (matrix of size n x n)
% B         System matrix B (matrix of size n x m)
% C         System matrix C (matrix of size l x n)
% D         System matrix D (matrix of size l x m)
% x0        Initial state (vector of size n x one)
% sv        Singular values (vector of size n x one)
    m=1;%inputs
    l=1;%outputs
    switch method
        case 'moesp'
            %Step 1 - LQ factorization
            Han_u = hankel(u(1:l*s),u(s*l:end));%block Hankel matrix of the input 
            Y = hankel(y(1:l*s),y(s*l:end));%block Hankel matrix of the output sequence y
            L = (triu(qr([Han_u ; Y]')))';
            
            %Step 2 - SVD of L22-determinition of n,A,C
            L22 = L(s+1:2*s,s+1:2*s);
            [Un,Sn,Vn] = svd(L22);
            sv = diag(Sn);%Singular values
            A = Un(1:s*l-l,1:n)\Un(l+1:end,1:n);
            C = Un(1:l,1:n);
            
            %Step3 - Determine x0,B,D
            I = eye(l); 
            Phi = zeros(length(y),n+n*m +l*m);
            for k =0:(length(y)-1)
                SUM =zeros(1,n);
                for j =0:k-1
                    SUM = SUM + kron(u(j+1),C*A^(k-1-j));
                end
                Phi(k+1,:) = [C*A^k  SUM  kron(u(k+1,1),I)];
            end
            BD = Phi\y;
            x0 = BD(1:n);
            B = BD(n+1:2*n);
            D = BD(2*n+1);
        case 'pi-moesp'
            Up = hankel(u(1:m*s),u(s*m:length(y)-s));%past input 
            Uf = hankel(u(s+1:2*s), u(2*s*m:end)); %future input
            Yf = hankel(y(s+1:2*s),y(2*s:end)); %future output
            R =(triu(qr([Uf;Up;Yf]')))'; %LQ factorization
            R32 = R(2*s+1:2*s*m+s*l,s*m+1:s*m+s*l); %take R32 element
            
            [Un,Sn,Vn] = svd(R32);%SVD of R32
            sv = diag(Sn);%Singular values
            A = Un(1:s*l-l,1:n)\Un(l+1:end,1:n);
            C = Un(1:l,1:n);
            
            %Step3 - Determine x0,B,D
            I = eye(l); 
            Phi = zeros(length(y),n+n*m +l*m);
            for k =0:(length(y)-1)
                SUM =zeros(1,n);
                for j =0:k-1
                    SUM = SUM + kron(u(j+1),C*A^(k-1-j));
                end
                Phi(k+1,:) = [C*A^k  SUM  kron(u(k+1,1),I)];
            end
            BD = Phi\y;
            x0 = BD(1:n);
            B = BD(n+1:2*n);
            D = BD(2*n+1);
            
        case 'po-moesp'
            Up = hankel(u(1:m*s),u(s*m:length(y)-s));
            Uf = hankel(u(s+1:2*s), u(2*s*m:end));
            Yp = hankel(y(1:l*s),y(s*l:length(y)-s));
            Yf = hankel(y(s+1:2*s),y(2*s:end));  
            
            %-Step1 - LQ factorization
            Z =[Up;Yp];
            R =(triu(qr([Uf;Z;Yf]')))';
            R32 = R(3*s+1:3*s*m+s*l,s*m+1:s*m+s*l);% take R32
            
            %Step 2 - SVD
            [Un,Sn,Vn] = svd(R32);
            sv = diag(Sn);%Singular values
            A = Un(1:s*l-l,1:n)\Un(l+1:end,1:n);
            C = Un(1:l,1:n);
            
            %Step3 - Determine x0,B,D
            I = eye(l); 
            Phi = zeros(length(y),n+n*m +l*m);
            for k =0:(length(y)-1)
                SUM =zeros(1,n);
                for j =0:k-1
                    SUM = SUM + kron(u(j+1),C*A^(k-1-j));
                end
                Phi(k+1,:) = [C*A^k  SUM  kron(u(k+1,1),I)];
            end
            BD = Phi\y;
            x0 = BD(1:n);
            B = BD(n+1:2*n);
            D = BD(2*n+1);
    end
       
end