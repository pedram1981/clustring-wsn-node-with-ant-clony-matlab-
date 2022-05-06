function ant_wsn()

clc;
clear;
close all;
global  E_elec E_amp node_count screen paket_size position_sensor station ant count_cluster  cluster_center
clc;
E_elec=50;
E_amp=10;
node_count=200;
screen=200;
paket_size=4000;
station=[floor(screen/2) floor(screen/2)];
count_cluster=10;
position_sensor=[];
nVar=node_count;
set_position();

MaxIt=100;      

nAnt=50;        

Q=1;

tau0=10;        

alpha=0.3;      

rho=0.1;        

tau=tau0*ones(nVar,count_cluster);

BestCost=zeros(MaxIt,1); 

empty_ant.Tour=[];
empty_ant.Cost=[];

ant=repmat(empty_ant,nAnt,1);

BestSol.Cost=inf;

for it=1:MaxIt
   clc;
   it
    for k=1:nAnt
        
        ant(k).Tour=zeros(1,nVar);
        
        for l=1:nVar
            
            P=tau(l,:).^alpha;
            
            P=P/sum(P);
            
            j=RouletteWheelSelection(P);
            
            ant(k).Tour(l)=j;
            
        end
        
        ant(k).Cost=cost(ant(k).Tour);
        
        if ant(k).Cost<BestSol.Cost
            BestSol=ant(k);
            BestCost(it)=ant(k).Cost;
        end
        
    end
    if(BestCost(it)==0)
        BestCost(it)=BestCost(it-1);
    end
   
    for k=1:nAnt
        
        tour=ant(k).Tour;
        
        for l=1:nVar
            
            tau(tour(l),:)=tau(tour(l),:)+Q/ant(k).Cost;
            
        end
        
    end
    
    tau=(1-rho)*tau;
    ant_=1+ceil((nAnt-1)*rand(1,1));
    Auction(BestSol,ant(ant_),ant_);

        
end
BestCost(it)
figure;
plot(BestCost,'b','lineWidth', 1.5);
hold on
 
xlabel('Number of iteration','FontSize',14,'FontWeight','bold');
ylabel('Enerjy','FontSize',14,'FontWeight','bold');

legend('ant clony with Auction');
set(gca,'FontName','Tahoma','fontsize',14)
grid on;

end

function set_position()
global  position_sensor node_count
for i=1:node_count
    position_sensor(i,1)=floor(100*rand(1,1));
    position_sensor(i,2)=floor(100*rand(1,1));
end
end

function y=cost(ant1)
global E_elec E_amp node_count screen paket_size position_sensor station  count_cluster  cluster_center 
set_center(ant1);
p1=[];
l=paket_size/100;
 sum_energy=0;
for i=1:count_cluster
   for j=1:node_count
        if(ant1(j)==i && j~=cluster_center(i) && cluster_center(i)~=0)
           d=sqrt(((position_sensor(cluster_center(i),1)-position_sensor(j,1))^2)+((position_sensor(cluster_center(i),2)-position_sensor(j,2))^2));
           Er=(E_elec*l)+(E_amp*l*d*d);     
           ER=(E_elec*l);
           sum_energy=Er+ER+sum_energy;
        end
    end
end

for i=1:count_cluster
        if(cluster_center(i)~=0)
           d=sqrt(((position_sensor(cluster_center(i),1)-(screen/2))^2)+((position_sensor(cluster_center(i),2)-(screen/2))^2));
           Er=(E_elec*l)+(E_amp*l*d*d);     
           ER=(E_elec*l);
           sum_energy=Er+ER+sum_energy;
        end
end
 
y=sum_energy;

end

function set_center(ant1)
global  node_count  position_sensor  count_cluster  cluster_center
p1=[];

cluster_center=[];
for i=1:count_cluster
    count=0;
    for j=1:node_count
       if(ant1(j)==i)
           count=count+1;
        p1(count,1)=j;   
        
       end
    end
   index=0;
    temp=10000000;
    for i1=1:count
       sum_=0;
         for j1=1:count
        sum_=sqrt(((position_sensor(p1(i1),1)-position_sensor(p1(j1),1))^2)+((position_sensor(p1(i1),2)-position_sensor(p1(j1),2))^2))+sum_;
      end
     
      if(sum_<temp)
          temp=sum_;
          index=p1(i1);
      end
    end
   
    cluster_center(i)=index;
end

end

function Auction(best_ant,ant_,id)
global ant node_count
mask=round(rand(1,node_count));
ant1=zeros(1,node_count+1);
ant2=zeros(1,node_count+1);
for i=1:node_count
    if(mask(i)==1)
        ant1(i)=best_ant.Tour(i);
        ant2(i)=ant_.Tour(i);
    else
        ant2(i)=best_ant.Tour(i);
        ant1(i)=ant_.Tour(i);
    end
end
            cost_new1 =cost(ant1);
            cost_new2 = cost(ant2);
            if (cost_new1 < ant_.Cost)
              ant(id).Tour=ant1;
              ant(id).Cost=cost_new1;
            end
            if (cost_new2 < cost_new1)
                 ant(id).Tour=ant2;
                 ant(id).Cost=cost_new2;
            end
            
end