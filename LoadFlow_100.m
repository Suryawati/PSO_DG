function [RugiRugi,ProblemSistem,VMin,VMax,Qg,bahayaGenerator]=LoadFlow_100(busdata,linedata,P,Q,Position,jum_DG,k,GenRestric)
% ===========load flow ==============================%
basemva = 100;  accuracy = 0.001; accel = 1.8; maxiter = 100;
%        Bus Bus  Voltage Angle   ---Load---  -----Generator-----   Shunt
%        No  code Mag.    Degree  MW    Mvar   MW   Mvar Qmin Qmax   Mvar
%        IEEE 30-BUS TEST SYSTEM
%        Bus Bus  Voltage Angle     ---Load---- -------Generator------------Static Mvar
%        No  code Mag.    Degree    MW    Mvar     MW       Mvar    Qmin  Qmax  +Qc/-Ql
% busdata=[1	  1	   1.05	     0	     0	     0	  260.2    -16.1	  0	   0	0
%          2	  2	   1.043    -5.48   21.7    12.7   40	    50	    -40	  50	0
%          3	  0	   1.021	-7.96	 2.4	 1.2	0	     0	      0	   0	0
%          4	  0	   1.012	-9.62	 7.6	 1.6	0	     0	      0	   0	0
%          5	  2	   1.01	   -14.37   94.2	19	    0	    37	    -40	  40	0
%          6	  0	   1.01	   -11.34	 0	     0	    0	     0	      0	   0	0
%          7	  0	   1.002   -13.12	22.8	10.9	0	     0	      0	   0	0
%          8	  2	   1.01	   -12.1	30	    30	    0	    37.3	-10	  40	0
%          9	  0	   1.051   -14.38	 0	     0	    0        0	      0	   0	0
%         10	  0	   1.045   -15.97	 5.8	 2	    0	     0	      0	   0	0.19
%         11	  2	   1.062	-14.39	 0	     0	    0	    16.2	-10	  24    -8.5
%         12	  0	   1.057	-15.24	11.2	7.5	    0	     0	      0	   0	0
%         13	  2	   1.05 	-15.24	 0	     0	    0	    10.6	-10	  24    -10
%         14	  0	   1.042	-16.13	6.2	    1.6	    0	     0	      0	   0	0
%         15	  0	   1.038	-16.22	8.2	    2.5	    0	     0	      0	   0	0
%         16	  0	   1.045	-15.83	3.5	    1.8	    0	     0	      0	   0	0
%         17	  0	   1.04	    -16.14	9	    5.8	    0	     0	      0	   0	0
%         18	  0	   1.028	-16.82	3.2	    0.9	    0	     0	      0	   0	0
%         19	  0	   1.026	-17	    9.5	    3.4	    0	     0	      0	   0	0
%         20	  0	   1.03	    -16.8	2.2	    0.7	    0	     0	      0	   0	0
%         21	  0	   1.033	-16.42	17.5	11.2	0	     0	      0	   0	0
%         22	  0	   1.033	-16.41	0	    0	    0	     0	      0	   0	0
%         23	  0	   1.027	-16.61	3.2	    1.6	    0	     0	      0	   0	0
%         24	  0	   1.021	-16.78	8.7	    6.7	    0	     0	      0	   0	0.043
%         25	  0	   1.017	-16.35	0	    0	    0	     0	      0	   0	0
%         26	  0	   1	    -16.77	3.5	    2.3	    0	     0	      0	   0	0
%         27	  0	   1.023	-15.82	0	    0	    0	     0	      0	   0	0
%         28	  0	   1.007	-11.97	0	    0	    0	     0	      0	   0	0
%         29	  0	   1.003	-17.06	2.4	    0.9	    0	     0	      0	   0	0
%         30	  0	   0.992	-17.94	10.6	1.9	    0	     0	      0	   0	0];
%     
    busdata(:,5:6)=busdata(:,5:6)*1;
    
for a=1:jum_DG
    busdata((Position(k,a)),7)=busdata((Position(k,a)),7)+P(k,a);
    busdata((Position(k,a)),8)=busdata((Position(k,a)),8)+Q(k,a);
end
      %                          Line code

%                                        Line code
%         Bus bus   R      X     1/2 B   = 1 for lines
%         nl  nr  p.u.   p.u.   p.u.     > 1 or < 1 tr. tap at bus nl
% linedata=[1   2   0.0192   0.0575   0.02640    1
%           1   3   0.0452   0.1852   0.02040    1
%           2   4   0.0570   0.1737   0.01840    1
%           3   4   0.0132   0.0379   0.00420    1
%           2   5   0.0472   0.1983   0.02090    1
%           2   6   0.0581   0.1763   0.01870    1
%           4   6   0.0119   0.0414   0.00450    1
%           5   7   0.0460   0.1160   0.01020    1
%           6   7   0.0267   0.0820   0.00850    1
%           6   8   0.0120   0.0420   0.00450    1
%           6   9   0.0      0.2080   0.0        1
%           6  10   0         .5560   0          1
%           9  11   0         .2080   0          1
%           9  10   0         .1100   0          1
%           4  12   0         .2560   0          1
%          12  13   0         .1400   0          1
%          12  14    .1231    .2559   0          1
%          12  15    .0662    .1304   0          1
%          12  16    .0945    .1987   0          1
%          14  15    .2210    .1997   0          1
%          16  17    .0524    .1923   0          1
%          15  18    .1073    .2185   0          1
%          18  19    .0639    .1292   0          1
%          19  20    .0340    .0680   0          1
%          10  20    .0936    .2090   0          1
%          10  17    .0324    .0845   0          1
%          10  21    .0348    .0749   0          1
%          10  22    .0727    .1499   0          1
%          21  22    .0116    .0236   0          1
%          15  23    .1000    .2020   0          1
%          22  24    .1150    .1790   0          1
%          23  24    .1320    .2700   0          1
%          24  25    .1885    .3292   0          1
%          25  26    .2544    .3800   0          1
%          25  27    .1093    .2087   0          1
%          28  27     0       .3960   0          1
%          27  29    .2198    .4153   0          1
%          27  30    .3202    .6027   0          1
%          29  30    .2399    .4533   0          1
%           8  28    .0636    .2000   0.0214     1
%           6  28    .0169    .0599   0.065      1];
%       
        


%---------------------------------------------------------------------------Pembentukan Matrix Ybus
%j=sqrt(-1); i = sqrt(-1);
nl = linedata(:,1);nr = linedata(:,2); R = linedata(:,3);
X = linedata(:,4); Bc = j*linedata(:,5); a = linedata(:, 6);
nbr=length(linedata(:,1));nbus = max(max(nl), max(nr));   
Z = R + j*X; y= ones(nbr,1)./Z; % admitansi cabang     
for n = 1:nbr
if a(n) <= 0,  a(n) = 1; else end
Ybus=zeros(nbus,nbus);     % inisialisasi Ybus  

% pembentukan elemen off diagonal 
for k=1:nbr;
       Ybus(nl(k),nr(k))=Ybus(nl(k),nr(k))-y(k)/a(k);
       Ybus(nr(k),nl(k))=Ybus(nl(k),nr(k));
    end
end
% pembentukan elemen diagonal
for  n=1:nbus
     for k=1:nbr
         if nl(k)==n
         Ybus(n,n) = Ybus(n,n)+y(k)/(a(k)^2) + Bc(k);
         elseif nr(k)==n
         Ybus(n,n) = Ybus(n,n)+y(k) +Bc(k);
         else, end
     end
end
Ybus;

%-------------------------------------------------------------------Load Flow dengan Newton-Raphson

ns=0; ng=0; Vm=0; delta=0; yload=0; deltad=0;
nbus = length(busdata(:,1));
for k=1:nbus
n=busdata(k,1);
kb(n)=busdata(k,2);Vm(n)=busdata(k,3);delta(n)=busdata(k, 4);
Pd(n)=busdata(k,5); Qd(n)=busdata(k,6); Pg(n)=busdata(k,7); Qg(n) = busdata(k,8);
Qmin(n)=busdata(k, 9); Qmax(n)=busdata(k, 10);
Qsh(n)=busdata(k, 11);
    if Vm(n) <= 0  Vm(n) = 1.0; V(n) = 1 + j*0;
    else delta(n) = pi/180*delta(n);
         V(n) = Vm(n)*(cos(delta(n)) + j*sin(delta(n)));
         P(n)=(Pg(n)-Pd(n))/basemva;
         Q(n)=(Qg(n)-Qd(n)+ Qsh(n))/basemva;
         S(n) = P(n) + j*Q(n);
    end
end
for k=1:nbus
if kb(k) == 1, ns = ns+1; else, end
if kb(k) == 2 ng = ng+1; else, end
ngs(k) = ng;
nss(k) = ns;
end
Ym=abs(Ybus); t = angle(Ybus);
m=2*nbus-ng-2*ns;
maxerror = 1; converge=1;
iter = 0;
% Mulai Iterasi
clear A  DC   J  DX
while maxerror >= accuracy & iter <= maxiter % Tes untuk Maks. Daya yang Tidak Sesuai
for i=1:m
for k=1:m
   A(i,k)=0;      %Inisialisasi Matrix Jacobian
end, end
iter = iter+1;
for n=1:nbus
nn=n-nss(n);
lm=nbus+n-ngs(n)-nss(n)-ns;
J11=0; J22=0; J33=0; J44=0;
   for i=1:nbr
     if nl(i) == n | nr(i) == n
        if nl(i) == n,  l = nr(i); end
        if nr(i) == n,  l = nl(i); end
        J11=J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        J33=J33+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        if kb(n)~=1
        J22=J22+ Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        J44=J44+ Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        else, end
        if kb(n) ~= 1  & kb(l) ~=1
        lk = nbus+l-ngs(l)-nss(l)-ns;
        ll = l -nss(l);
      % element off diagonal J1
        A(nn, ll) =-Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
              if kb(l) == 0  % element off diagonal J2
              A(nn, lk) =Vm(n)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));end
              if kb(n) == 0  % element off diagonal J3
              A(lm, ll) =-Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n)+delta(l)); end
              if kb(n) == 0 & kb(l) == 0  % element off diagonal J4
              A(lm, lk) =-Vm(n)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));end
        else end
     else , end
   end
   Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33;
   Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;
   if kb(n) == 1 P(n)=Pk; Q(n) = Qk; end   % Swing bus P
     if kb(n) == 2  Q(n)=Qk;
         if Qmax(n) ~= 0
           Qgc = Q(n)*basemva + Qd(n) - Qsh(n);
           if iter <= 7                  % Antara Iterasi Ke-2 dan Ke-6
              if iter > 2                % MVAR Bus-Bus Generator dites.
                if Qgc  < Qmin(n),       % Jika tidak dalam batas Vm(n)
                Vm(n) = Vm(n) + 0.01;    % Maka Dirubah ke Langkah 0.01pu
                elseif Qgc  > Qmax(n),   % untuk Memberi MVAR dengan Batas 
                Vm(n) = Vm(n) - 0.01;end % yang telah Ditentukan
              else, end
           else,end
         else,end
     end
   if kb(n) ~= 1
     A(nn,nn) = J11;  %element diagonal J1
     DC(nn) = P(n)-Pk;
   end
   if kb(n) == 0
     A(nn,lm) = 2*Vm(n)*Ym(n,n)*cos(t(n,n))+J22;  %element diagonal J2
     A(lm,nn)= J33;        %element diagonal J3
     A(lm,lm) =-2*Vm(n)*Ym(n,n)*sin(t(n,n))-J44;  %element diagonal J4
     DC(lm) = Q(n)-Qk;
   end
end
DX=A\DC';
for n=1:nbus
  nn=n-nss(n);
  lm=nbus+n-ngs(n)-nss(n)-ns;
    if kb(n) ~= 1
    delta(n) = delta(n)+DX(nn); end
    if kb(n) == 0
    Vm(n)=Vm(n)+DX(lm); end
 end
  maxerror=max(abs(DC));
     if iter == maxiter & maxerror > accuracy 
   fprintf('\nPERINGATAN : Solusi Iteratif tidak Konvergen Setelah')
   fprintf('%g', iter), fprintf(' iterasi.\n\n')
   fprintf('Tekan Enter untuk Mengakhiri Iterasi dan Cetak Hasil Komputasi\n')
   converge = 0; pause, else, end
   
end

if converge ~= 1
   tech= ('                      SOLUSI ITERATIF TIDAK KONVERGEN'); else, 
   tech=('                   Solusi Aliran Daya dengan Metode Newton-Raphson');
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
     Qgg(k)=Qg(n);     %june 97
     elseif  kb(n) ==2
     k=k+1;
     S(n)=P(n)+j*Q(n);
     Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
     Pgg(k)=Pg(n);
     Qgg(k)=Qg(n);  % June 1997
  end
yload(n) = (Pd(n)- j*Qd(n)+j*Qsh(n))/(basemva*Vm(n)^2);
end
busdata(:,3)=Vm'; busdata(:,4)=deltad';
Pgt = sum(Pg);  Qgt = sum(Qg); Pdt = sum(Pd); Qdt = sum(Qd); Qsht = sum(Qsh);

%clear A DC DX  J11 J22 J33 J44 Qk delta lk ll lm
%clear A DC DX  J11 J22 J33  Qk delta lk ll lm


% %---------------------------------------------------------------------------------------Data BusOut
% disp(tech)
% fprintf('                      Maksimum Daya Tidak Sesuai = %g \n', maxerror)
% fprintf('                             Nomor Iterasi = %g \n\n', iter)
% head =['    Bus  Voltage  Angle    ------Load------    ---Generation---   Injected'
%        '    No.  Mag.     Degree     MW       Mvar       MW       Mvar       Mvar '
%        '                                                                          '];
% disp(head)
% for n=1:nbus
%      fprintf(' %5g', n), fprintf(' %7.3f', Vm(n)),
%      fprintf(' %8.3f', deltad(n)), fprintf(' %9.3f', Pd(n)),
%      fprintf(' %9.3f', Qd(n)),  fprintf(' %9.3f', Pg(n)),
%      fprintf(' %9.3f ', Qg(n)), fprintf(' %8.3f\n', Qsh(n))
% end
%     fprintf('      \n'), fprintf('    Total              ')
%     fprintf(' %9.3f', Pdt), fprintf(' %9.3f', Qdt),
%     fprintf(' %9.3f', Pgt), fprintf(' %9.3f', Qgt), fprintf(' %9.3f\n\n', Qsht)
%     
%----------------------------------------------------------------Line Flow

%============================CONSTRAIN TEGANGAN=================
%voltageProblem=0;
VMin=min(Vm);
VMax=max(Vm);
if VMin<0.9 || VMax>1.1
    voltageProblem=1;
else
    voltageProblem=0;    
end

an=length(GenRestric);
bahayaGenerator(1)=0;
for n=2:an
    %bahayaGenerator(n)=0+bahayaGenerator(n-1);
    if Qg(GenRestric(n))<busdata(GenRestric(n),9) || Qg(GenRestric(n))>busdata(GenRestric(n),10)
        bahayaGenerator(n)=1;
    else
        bahayaGenerator(n)=0;
    end
end

%===============================================================
SLT = 0;
% fprintf('\n')
% fprintf('                           Line Flow and Losses \n\n')
% fprintf('     --Line--  Power at bus & line flow    --Line loss--  Transformer\n')
% fprintf('     from  to    MW      Mvar     MVA       MW      Mvar      tap\n')
for n = 1:nbus
busprt = 0;
   for L = 1:nbr;
       if busprt == 0
       P(n)*basemva;%fprintf('   \n'), fprintf('%6g', n), fprintf('      %9.3f', P(n)*basemva)
       Q(n)*basemva;%abs(S(n)*basemva);fprintf('%9.3f', Q(n)*basemva), fprintf('%9.3f\n', abs(S(n)*basemva))

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
         %fprintf('%12g', k),
         real(Snk);imag(Snk);%fprintf('%9.3f', real(Snk)), fprintf('%9.3f', imag(Snk))
         abs(Snk);%fprintf('%9.3f', abs(Snk)),
         real(SL);%fprintf('%9.3f', real(SL)),
             if nl(L) ==n & a(L) ~= 1
             imag(SL);a(L);%fprintf('%9.3f', imag(SL)), fprintf('%9.3f\n', a(L))
             else,imag(SL); %fprintf('%9.3f\n', imag(SL))
             end
         else, end
  end
end
 SLT = SLT/2;
% fprintf('   \n'), fprintf('    Total loss                         ')
% real(SLT);imag(SLT);fprintf('%9.4f', real(SLT)), fprintf('%9.4f\n', imag(SLT))

%===============================resul constrain========================
if voltageProblem==1 || sum(bahayaGenerator)>0
    ProblemSistem=1;
else
    ProblemSistem=0;
end
%======================================================================
RugiRugi=real(SLT);
% %clear Ik In SL SLT Skn Snk