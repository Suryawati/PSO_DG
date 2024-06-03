  %program konfigurasi peletakan DG menggunakan PSO%
  clc;
  clear all;
  close all;
  tic
  
  %inisialisasi parameter PSO
   iter_max=50;    % iterasi maksimum
   c1=0.7;         % acceleration 1
   c2=0.9;         % acceleration 2
   w=0.40;         % weigth
   wmax=3;         % weigth maximum
   wmin=1;         % weigth minimum
   swarm=10;       % jumlah swarm (jika di GA kita familiar dengan populasi)
   jum_DG=3;       % jumlah DG yang akan kita letakan
   Pmin=0;              % isikan P minimum DG
   Pmax=10;             % isikan P maksimum DG
   Qmin=-2;             % isikan Q minimum DG
   Qmax=-0.001;         % isikan Q maksimum DG
   jum_bus=30; %data (30 0r 14 Bus)
   stmbus=jum_bus;
   [busdata,linedata,GenRestric]=DataIEEE(stmbus);
   
   for i= 1 : iter_max
        wt(i)=wmax-((wmax-wmin)/iter_max)*i ;
    end
        
    P=(Pmax-(rand(swarm,jum_DG)*(Pmax-Pmin)));
    Q=(Qmax-(rand(swarm,jum_DG)*(Qmax-Qmin)));
    V=rand(swarm,jum_DG*3)*0;
    Position=round(jum_bus-(rand(swarm,jum_DG)*(jum_bus-3)));
    for a=1:swarm
         for b=1:jum_DG
             if ismember(Position(a,b),GenRestric)==1,
                 Position(a,b)=10;                                 
             elseif Position(a,b)>jum_bus 
                 Position(a,b)=jum_bus;
             elseif Position(a,b)<1 
                 Position(a,b)=10;
             end
         end
     end             
    X=[P Q Position];     
      
    
    Rugi =zeros(swarm,1);           %pembentukan tabel rugi dalam satu kolom    
    
    %masuk ke load flow untuk menghitung rugi-rugi daya    
    
    for k=1:swarm;
    [RugiRugi,ProblemSistem,VMin,Vmax,Qg]=LoadFlow_100(busdata,linedata,P,Q,Position,jum_DG,k,GenRestric);  %perhitungan rugi2 
    ProbSis(k)=ProblemSistem;
    Rugi_1(k)=RugiRugi;               %rugi daya di setiap particle
    end

     for k=1:swarm
       %Rugi(k)=Rugi100(k); %rugi rata2
       if ProbSis(k)==0
           Rugi(k)=Rugi_1(k);
       elseif ProbSis(k)==1
           Rugi(k)=100;
       end
    end
    %----------------------------------------------------------------------
   
    [Fbest(1),C]=min(Rugi);  %C menunjukan indeks dan Fbest adalah rugi terbaik
    Pbest(1,:)=X(C,:);       %partikel terbaik
    
    Fgraph(1)=Rugi(C,:);     
     
    [Fgbest,Iterbest]=min(Fbest);  %menentukan rugi terbaik secara global dari seluruh iterasi
    GlobalBest=Pbest(Iterbest,:);  %Global best
    FglobalBest(1)=Fgbest;
    
    %update velocity
    for a=1:swarm
    V(a,:)=w*V(a,:)+c1*rand*(Pbest(1,:)-X(a,:))+c2*rand*(GlobalBest-X(a,:));
    end
    %update posisi
    X=(X+V);
    
%---Buat Grafik
hfig = figure;
hold on
title('Grafik Optimisasi Peletakan DG dengan PSO');
set(hfig, 'position', [50,40,600,300]);
set(hfig, 'DoubleBuffer', 'on');
hbestplot = plot(1:iter_max,zeros(1,iter_max));
htext1 = text(0.6*iter_max,30,sprintf('Rugi daya : %f', 0.0));
xlabel('Iterasi');
ylabel('Fungsi Objektif');
hold off
drawnow;
% 
   it=1;
   while it<=iter_max,
        it=it+1;
       
    %----------------------------------------------------------------------
    %masuk ke load flow untuk menghitung rugi-rugi daya      
   
     P=X(:,(1:jum_DG));
     Q=X(:,(jum_DG+1:jum_DG*2));
     Position=round(X(:,((jum_DG*2)+1):jum_DG*3));      
         
     for a=1:swarm
         for b=1:jum_DG
             if ismember(Position(a,b),GenRestric)==1,
                 Position(a,b)=10;                                 
             elseif Position(a,b)>jum_bus 
                 Position(a,b)=jum_bus;
             elseif Position(a,b)<1 
                 Position(a,b)=10;
             end
         end
     end             
         
     for a=1:swarm
         for b=1:jum_DG
             if P(a,b)<Pmin
                 P(a,b)=Pmin;
             elseif P(a,b)>Pmax
                 P(a,b)=Pmin+(0.25*Pmin);
             end
         end
     end
     for a=1:swarm
         for b=1:jum_DG
             if Q(a,b)<Qmin
                 Q(a,b)=Qmin;
             elseif Q(a,b)>Qmax
                 Q(a,b)=Qmax+(0.25*Qmin);
             end
         end
     end
  
     X=[P Q Position];
    
      
    for k=1:swarm;
    [RugiRugi,ProblemSistem,VMin,Vmax,Qg]=LoadFlow_100(busdata,linedata,P,Q,Position,jum_DG,k,GenRestric);  %perhitungan rugi2
    ProbSis(k)=ProblemSistem;
    Rugi100(k)=RugiRugi;               %rugi daya di setiap particle
    end

     for k=1:swarm
       %Rugi(k)=Rugi100(k);
       if ProbSis(k)==0
           Rugi(k)=Rugi100(k);
       elseif ProbSis(k)==1
           Rugi(k)=100;
       end
     end
     
     %----------------------------------------------------------------------
   
      [Fbest(it),C]=min(Rugi);
      Pbest(it,:)=X(C,:);      
      Fgraph(it)=Rugi(C,:);
     
    [Fgbest,Iterbest]=min(Fbest);
    GlobalBest=Pbest(Iterbest,:);
    FglobalBest(it)=Fgbest;
    
        
    %update velocity
    for a=1:swarm
        V(a,:)=w*V(a,:)+c1*rand*(Pbest(it,:)-X(a,:))+c2*rand*(GlobalBest-X(a,:));    
    end
    %update Position
    X=(X+V);
    
    %mutasi
            for a=1:swarm
                p=rand;
                if p<0.3
                d=round((jum_DG*3)*rand()+0.5);
                  if d<=8
                    X(a,d)=Pmax-(X(a,d)-Pmin);
                  elseif d>8 & d<=16,
                   X(a,d)=Qmax-(X(a,d)-Qmin);
                  elseif d>16
                   X(a,d)=jum_bus-(X(a,d)-2);
                  end
                end            
            end
        
    
    plotvector = get(hbestplot,'YData');
    plotvector(it-1) = FglobalBest(it-1);
    set(hbestplot,'YData',plotvector);
    set(htext1,'String',sprintf('Rugi daya optimal: %f', FglobalBest(it-1)));
    hold on;
    drawnow 
    
end
   
 
 clc;
 toc
 [RugiRugi]=LF_nonDG(busdata,linedata);
 Rugi_TanpaDG=RugiRugi;
 fprintf('\n')
 busPQ=GlobalBest;
 [RugiRugi]=LF(busdata,linedata,busPQ,jum_DG);
 Rugi_withDG=RugiRugi;
 fprintf('\n')
 
 
 busP=[busPQ(1:jum_DG)',busPQ((jum_DG*2)+1:jum_DG*3)'];
 busP=sortrows(busP,2);
 busQ=[busPQ(jum_DG+1:jum_DG*2)',busPQ((jum_DG*2)+1:jum_DG*3)'];
 busQ=sortrows(busQ,2);
 
 fprintf('\n')
 fprintf('Rugi Tanpa DG : %f \n',Rugi_TanpaDG);
 fprintf('\n')
 fprintf('Rugi setelah ada DG : %f \n',Rugi_withDG);
 fprintf('\n')
 fprintf(' Daya Aktif(P)     no Bus       Kapasitas\n')
 for a=1:jum_DG
     fprintf('%10.3f',busP(a,1)); fprintf('          %1.f',busP(a,2));      fprintf('          %1.f MW\n',Pmax);                
 end
 fprintf('\n')
 fprintf(' Daya Reaktif(Q)   no Bus       Kapasitas\n')
 for a=1:jum_DG
     fprintf('%10.3f',busQ(a,1));  fprintf('         %1.f',busQ(a,2));      fprintf('          %1.f MVar\n',Qmax);          
 end
