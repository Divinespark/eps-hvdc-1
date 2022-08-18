basemva=100; accuracy=0.000001; accel=1.6; maxiter=100;

%% Step 1: Busdata of test case
global busdata linedata
%        Bus Bus  Voltage Angle   ---Load---- -------Generator----- Static Mvar
%        No  code Mag.    Degree  MW    Mvar  MW  Mvar Qmin Qmax    +Qc/-Ql
busdata=[
                1   1    1.03    0         0        0     0   0   0   0  0
                2   0    1.05    0        80      25    0   0   0   0  0
                3   2    1.00    0         0        0     0   0   0   0  0];


%% Step 2: Linedata of test case
%         Bus bus   R      X     1/2 B   = 1 for lines
%         nl  nr  p.u.   p.u.   p.u.     > 1 or < 1 tr. tap at bus nl
linedata=[
1   2   0.0    0.15    0     1
2   3   0.0    0.06    0     1.079/1]; 

%% Step 3: Obtain the Bus Admittance Matrix for power flow solution
j=sqrt(-1); i = sqrt(-1);
nl = linedata(:,1); nr = linedata(:,2); R = linedata(:,3);
X = linedata(:,4); Bc = j*linedata(:,5); a = linedata(:, 6);
nbr=length(linedata(:,1)); nbus = max(max(nl), max(nr));
Z = R + j*X; y= ones(nbr,1)./Z;        %branch admittance
for n = 1:nbr
if a(n) <= 0  
    a(n) = 1; else end
Ybus=zeros(nbus,nbus);     % initialize Ybus to zero
               % formation of the off diagonal elements
for k=1:nbr;
       Ybus(nl(k),nr(k))=Ybus(nl(k),nr(k))-y(k)/a(k);
       Ybus(nr(k),nl(k))=Ybus(nl(k),nr(k));
    end
end
              % formation of the diagonal elements
for  n=1:nbus
     for k=1:nbr
         if nl(k)==n
         Ybus(n,n) = Ybus(n,n)+y(k)/(a(k)^2) + Bc(k);
         elseif nr(k)==n
         Ybus(n,n) = Ybus(n,n)+y(k) +Bc(k);
         else, end
     end
end
clear Pgg


%% Step 4: Power flow solution by Newton-Raphson method
ns=0; ng=0; Vm=0; delta=0; yload=0; deltad=0;
nbus = length(busdata(:,1));
kb=[];Vm=[]; delta=[]; Pd=[]; Qd=[]; Pg=[]; Qg=[]; Qmin=[]; Qmax=[];  
Pk=[]; P=[]; Qk=[]; Q=[]; S=[]; V=[]; 
for k=1:nbus
    n=busdata(k,1);
    kb(n)=busdata(k,2); Vm(n)=busdata(k,3); delta(n)=busdata(k, 4);
    Pd(n)=busdata(k,5); Qd(n)=busdata(k,6); Pg(n)=busdata(k,7); Qg(n) = busdata(k,8);
    Qmin(n)=busdata(k, 9); Qmax(n)=busdata(k, 10);
    Qsh(n)=busdata(k, 11);
    if Vm(n) <= 0  
        Vm(n) = 1.0; V(n) = 1 + j*0;
    else delta(n) = pi/180*delta(n);
         V(n) = Vm(n)*(cos(delta(n)) + j*sin(delta(n)));
         P(n)=(Pg(n)-Pd(n))/basemva;
         Q(n)=(Qg(n)-Qd(n)+ Qsh(n))/basemva;
         S(n) = P(n) + j*Q(n);
    end
end
for k=1:nbus
    if kb(k) == 1
        ns = ns+1; 
    else
    end
    if kb(k) == 2
        ng = ng+1; 
    else
    end
    ngs(k) = ng;
    nss(k) = ns;
end
Ym=abs(Ybus); t = angle(Ybus);
m=2*nbus-ng-2*ns;
maxerror = 1; converge=1;
iter = 0;
%%%% added for parallel lines
mline=ones(nbr,1);
for k=1:nbr
      for m=k+1:nbr
         if((nl(k)==nl(m)) && (nr(k)==nr(m)));
            mline(m)=2;
         elseif ((nl(k)==nr(m)) && (nr(k)==nl(m)));
         mline(m)=2;
         else
         end
      end
   end
%%%   end of statements for parallel lines   

%% Step 5: Start of iterations
clear A  DC   J  DX
Voltage2(1,1)=1;
while maxerror >= accuracy && iter <= maxiter % Test for max. power mismatch
    for ii=1:m
        for k=1:m
            A(ii,k)=0;      %Initializing Jacobian matrix
        end
    end
iter = iter+1;
for n=1:nbus
nn=n-nss(n);
lm=nbus+n-ngs(n)-nss(n)-ns;
J11=0; J22=0; J33=0; J44=0;
   for ii=1:nbr
   	if mline(ii)==1   % Added to include parallel lines
      	if nl(ii) == n || nr(ii) == n
         	if nl(ii) == n 
                l = nr(ii); 
            end
         	if nr(ii) == n 
                l = nl(ii); 
            end
         J11=J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
         J33=J33+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        		if kb(n)~=1
                    J22=J22+ Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
                    J44=J44+ Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
                else
                end
        		if kb(n) ~= 1  && kb(l) ~=1
        			lk = nbus+l-ngs(l)-nss(l)-ns;
        			ll = l -nss(l);
      			% off diagonalelements of J1
                    A(nn, ll) =-Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
                    if kb(l) == 0  % off diagonal elements of J2
                        A(nn, lk) =Vm(n)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
                    end
                    if kb(n) == 0  % off diagonal elements of J3
                        A(lm, ll) =-Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n)+delta(l));
                    end
                    if kb(n) == 0 && kb(l) == 0  % off diagonal elements of  J4
                        A(lm, lk) =-Vm(n)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
                    end
                else
                    end
        else
        end
    else
    end
   end
   Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33;
   Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;
   if kb(n) == 1 
       P(n)=Pk; Q(n) = Qk; 
   end   % Swing bus P
     if kb(n) == 2  
         Q(n)=Qk;
         if Qmax(n) ~= 0
           Qgc = Q(n)*basemva + Qd(n) - Qsh(n);
           if iter <= 7                  % Between the 2th & 6th iterations
              if iter > 2                % the Mvar of generator buses are
                if Qgc  < Qmin(n),       % tested. If not within limits Vm(n)
                Vm(n) = Vm(n) + 0.01;    % is changed in steps of 0.01 pu to
                elseif Qgc  > Qmax(n),   % bring the generator Mvar within
                Vm(n) = Vm(n) - 0.01;
                end % the specified limits.
              else
              end
           else
           end
         else
         end
     end
   if kb(n) ~= 1
     A(nn,nn) = J11;  %diagonal elements of J1
     DC(nn) = P(n)-Pk;
   end
   if kb(n) == 0
     A(nn,lm) = 2*Vm(n)*Ym(n,n)*cos(t(n,n))+J22;  %diagonal elements of J2
     A(lm,nn)= J33;        %diagonal elements of J3
     A(lm,lm) =-2*Vm(n)*Ym(n,n)*sin(t(n,n))-J44;  %diagonal of elements of J4
     DC(lm) = Q(n)-Qk;
   end
end
DX=A\DC';
for n=1:nbus
  nn=n-nss(n);
  lm=nbus+n-ngs(n)-nss(n)-ns;
    if kb(n) ~= 1
        delta(n) = delta(n)+DX(nn); 
    end
    if kb(n) == 0
        Vm(n)=Vm(n)+DX(lm); 
    end
   if length(busdata(:,1))~= linedata(2,2)
  clear; clc;
  end
 end
%maxerror=max(abs(DC)); %taking power as convergence condition
Voltage2(iter+1,1)=Vm(1,2);
maxerror=abs(Voltage2(iter+1,1)-Voltage2(iter,1)); %taking voltage as convergence condition
     if iter == maxiter && maxerror > accuracy 
         fprintf('\nWARNING: Iterative solution did not converged after ')
         fprintf('%g', iter), fprintf(' iterations.\n\n')
         fprintf('Press Enter to terminate the iterations and print the results \n')
         converge = 0; 
         pause
     else
     end
end

if converge ~= 1
   tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); else, 
   tech=('                   Power Flow Solution by Newton-Raphson Method');
end   
V = Vm.*cos(delta)+j*Vm.*sin(delta);
deltad=180/pi*delta;
i=sqrt(-1);
k=0;
for n = 1:nbus
     if kb(n) == 1
         k=k+1;
         S(n)= P(n)+j*Q(n);
         Pg(n) = P(n)*basemva + Pd(n);
         Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
         Pgg(k)=Pg(n);
         Qgg(k)=Qg(n);
     elseif  kb(n) ==2
         k=k+1;
         S(n)=P(n)+j*Q(n);
         Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
         Pgg(k)=Pg(n);
         Qgg(k)=Qg(n);  
  end
yload(n) = (Pd(n)- j*Qd(n)+j*Qsh(n))/(basemva*Vm(n)^2);
end
busdata(:,3)=Vm';
busdata(:,4)=deltad';
Pgt = sum(Pg); 
Qgt = sum(Qg); 
Pdt = sum(Pd);
Qdt = sum(Qd); 
Qsht = sum(Qsh);


%% Step 6:  prints the power flow solution in a tabulated form on the screen.
disp(tech)
fprintf('                      Maximum Power Mismatch = %g \n', maxerror)
fprintf('                             No. of Iterations = %g \n\n', iter)
head =['    Bus  Voltage  Angle    ------Load------    ---Generation---   Injected'
       '    No.  Mag.     Degree     MW       Mvar       MW       Mvar       Mvar '
       '                                                                          '];
disp(head)
for n=1:nbus
     fprintf(' %5g', n), fprintf(' %7.3f', Vm(n)),
     fprintf(' %8.3f', deltad(n)), fprintf(' %9.3f', Pd(n)),
     fprintf(' %9.3f', Qd(n)),  fprintf(' %9.3f', Pg(n)),
     fprintf(' %9.3f ', Qg(n)), fprintf(' %8.3f\n', Qsh(n))
end
    fprintf('      \n'), fprintf('    Total              ')
    fprintf(' %9.3f', Pdt), fprintf(' %9.3f', Qdt),
    fprintf(' %9.3f', Pgt), fprintf(' %9.3f', Qgt), fprintf(' %9.3f\n\n', Qsht)%  This program is used in conjunction with lfgauss or lf Newton
%  for the computation of line flow and line losses.

SLT = 0;
fprintf('\n')
fprintf('                           Line Flow and Losses \n\n')
fprintf('     --Line--  Power at bus & line flow    --Line loss--  Transformer\n')
fprintf('     from  to    MW      Mvar     MVA       MW      Mvar      tap\n')
for n = 1:nbus
busprt = 0;
   for L = 1:nbr;
       if busprt == 0
       fprintf('   \n'), fprintf('%6g', n), fprintf('      %9.3f', P(n)*basemva)
       fprintf('%9.3f', Q(n)*basemva), fprintf('%9.3f\n', abs(S(n)*basemva))

       busprt = 1;
       else, end
       if nl(L)==n      k = nr(L);
       In = (V(n) - a(L)*V(k))*y(L)/a(L)^2 + Bc(L)/a(L)^2*V(n);
       Ik = (V(k) - V(n)/a(L))*y(L) + Bc(L)*V(k);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       elseif nr(L)==n  k = nl(L);
       In = (V(n) - V(k)/a(L))*y(L) + Bc(L)*V(n);
       Ik = (V(k) - a(L)*V(n))*y(L)/a(L)^2 + Bc(L)/a(L)^2*V(k);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       else, end
         if nl(L)==n | nr(L)==n
         fprintf('%12g', k),
         fprintf('%9.3f', real(Snk)), fprintf('%9.3f', imag(Snk))
         fprintf('%9.3f', abs(Snk)),
         fprintf('%9.3f', real(SL)),
             if nl(L) ==n & a(L) ~= 1
             fprintf('%9.3f', imag(SL)), fprintf('%9.3f\n', a(L))
             else, fprintf('%9.3f\n', imag(SL))
             end
         else, end
  end
end
SLT = SLT/2;
A=SLT;
fprintf('   \n'), fprintf('    Total loss                         ')
fprintf('%9.3f', real(SLT)), fprintf('%9.3f\n', imag(SLT))
clear Ik In SL SLT Skn Snk