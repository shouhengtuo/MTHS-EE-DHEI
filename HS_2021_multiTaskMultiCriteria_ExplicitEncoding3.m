 function [Task,NC,flag,DFes] = HS_2021_multiTaskMultiCriteria_ExplicitEncoding3(data,dim_epi,HMS,max_iter,CX,objFunction,K0,K,F)

%input--------------------------------------------------------------------
% data-----------------input dataset
% epi_dim--------------the epistasis dimension
% HMS--------------the size of harmony memory(HM)

%%-------------------------------------------------------------------------
% initial arguments
HMCR=0.98;
D=length(CX);
DFes = HMS;
PAR = 0.7;
TP = 0.2;
TP2 = 0.6;
n=size(data,2);
State=data(:,n);
% K = dim_epi - 1; % 任务数 2-->epi_dim 

CriteriaNum = 4;
%% ---------------------------------------------------------------

SNPs=n-1;  %% 总SNP个数
flag = 0;
ki = 0; 
taskNum = K - K0 + 1;
Orders = K0:K;
bestNum = 10;
for k = 1:taskNum
     ki = ki + 1;
    %% 初始化
   dim = Orders(ki); % 第个任务的维度(阶数）为k+1
    X = zeros(HMS,dim);
    HM = X;
    snp=[];
    for i=1:HMS
        snp(1)=ceil(rand*SNPs);        
        for j=2:dim % 检查每个个体中是否有重复的SNP 位点
          snp(j)=ceil(rand*SNPs); 
          while ismember(snp(j),snp(1:j-1)) 
             snp(j)=ceil(rand*SNPs);        
          end
        end
        temp=snp;
        snp2=sort(snp(1:dim));        
        while ismember(snp,X,'rows')
            j=ceil(rand*dim);
            snp(j)=ceil(rand*SNPs); 
            temp=snp;
            snp2=sort(snp(1:dim));
        end

        X(i,:)=snp2;   %% X中存放有序的解        
        HM(i,:)=temp;  %% HM中相应存放无序解      
        Scores(i,:) = objFunction(data(:,HM(i,:)),State); 
        snp=[];
    end
    
    for j = 1:CriteriaNum
        Criteria(j).X = X;
        Criteria(j).HM = HM;
        Criteria(j).Score = Scores(:,j);
         
        [bestScores,bestIds] = sort(Criteria(j).Score);
                
         Criteria(j).bestScores = bestScores(1:bestNum);
          Criteria(j).bestIds = bestIds(1:bestNum); 
           Criteria(j).bestXs = Criteria(j).X(bestIds(1:bestNum),:);
    end
    Task(ki).Criteria = Criteria;

end
% for i = 1:taskNum
%    Task(i).Criteria(1).bestX
% end


NC=HMS*taskNum;
 LT=0;
%%-------------------------------------------------------------------------
tic;
while NC <= max_iter 
    
    
   for k = 1 : taskNum
       if mod(NC,10) == 0
            for j = 1:CriteriaNum 
                 [bestScores,bestIds] = sort(Task(k).Criteria(j).Score);
                
                 Task(k).Criteria(j).bestScores = bestScores(1:bestNum);
                 Task(k).Criteria(j).bestIds = bestIds(1:bestNum); 
                Task(k).Criteria(j).bestXs = Task(k).Criteria(j).X(bestIds(1:bestNum),:);
                
            end
       end
      
       
       Xnew = [];
       Xtemp = [];
       dim = Orders(k);  %% 从第k个任务探索
       
       %transfer = ceil(rand*5 - 1);
       if rand < TP %transfer == 1 
           %%
           Flag = 1;          
               kr = ceil(rand*taskNum);
               while kr==k
                  kr = ceil(rand*taskNum);
               end
            
                CJ = ceil(rand*CriteriaNum);
               if k > kr %% 从阶数低的任务 向 阶数高的任务 整体迁移
                     a = ceil(rand*HMS);
                    
                     Xnew = Task(kr).Criteria(CJ).X(a,:);
                     i = kr+2;
                     while i <= dim
                         if rand < HMCR
                             a = ceil(rand*HMS);
                             kr0 = ceil(rand*taskNum);
                             j = ceil(rand*bestNum);
                             b = ceil(rand*(kr0+1));
                             x = Task(kr0).Criteria(CJ).HM(a,b);
                             if rand < PAR
                                 x = Task(kr0).Criteria(CJ).bestXs(j,b);
                             end
                         else
                             x = ceil(rand*SNPs);
                         end
                         if ~ismember(x,Xnew)
                             Xnew(i) = x;
                             i = i + 1;
                         end
                     end
                     Xtemp = Xnew;
                     Xnew = sort(Xnew);
               else  %% 从阶数高的任务选择迁移
                  selIndex = randperm(Orders(kr));
                   if rand < TP2  %  从相同标准的记忆库中选取
                        a = ceil(rand*HMS);
                        
                        Xnew = Task(kr).Criteria(CJ).X(a,selIndex(1:dim));
                   else
                        cj = ceil(rand*CriteriaNum);
                        a = ceil(rand*HMS);
                        Xnew = Task(kr).Criteria(cj).X(a,selIndex(1:dim));
                   end
                   Xtemp = Xnew;
               end
                    
                      while ( ismember(Xnew ,Task(k).Criteria(1).X,'rows')||ismember(Xnew ,Task(k).Criteria(2).X,'rows')||...
                              ismember(Xnew ,Task(k).Criteria(3).X,'rows')||ismember(Xnew ,Task(k).Criteria(4).X,'rows'))
                          j=ceil(rand*dim);
                          r=ceil(rand*SNPs);
                          while ismember(r,Xnew)
                              r=ceil(rand*SNPs);
                          end
                          Xnew(j)=r;
                          Xtemp=Xnew;
                          Xnew =sort(Xnew);
                      end
       else
           if rand < TP2  %% 从相同评价标准的 记忆库中组合学习
               CJ = ceil(rand*CriteriaNum);
             %% ----------------------------------------------
              i=1;
              while i<= dim
                      if rand<HMCR  
                             a = ceil(rand*HMS); 
                             b = ceil(rand*dim); 
        %                      cr = ceil(rand*CriteriaNum);
                             Xnew(i) = Task(k).Criteria(CJ).X(a,b);                     
                             if rand<PAR
                                 sPar = ceil(rand * 4);
                                 k0 = ceil(rand*taskNum);
                                         while k0 == k
                                              k0 = ceil(rand*taskNum);
                                         end
                                 switch sPar


                                     case 1
                                         
        %                                  c = ceil(rand*CriteriaNum);
                                         d = ceil(rand*HMS);
                                         j = ceil(rand*(k0+1));
                                         Xnew(i) = Task(k0).Criteria(CJ).X(d,j);
                                     case 2
                                         c = ceil(rand*CriteriaNum);
                                         d = ceil(rand*HMS);
                                         j = ceil(rand*(k+1));

                                              L = Task(k).Criteria(CJ).X(d,j) - Xnew(i);
                                         Xnew(i) = Xnew(i) + round(F*rand *(Task(k).Criteria(CJ).X(d,j) - Xnew(i)));
                                        if  Xnew(i) > SNPs
                                             Xnew(i) = SNPs - max(0,round(normrnd(0,min(5,max(L/10,1)))));
                                         elseif Xnew(i) < 1
                                             Xnew(i) = 1 + max(0,round(normrnd(0,min(5,max(L/10,1)))));
                                         end
                                     case 3
                                           
                                         c = ceil(rand*bestNum);
                                         d = ceil(rand*HMS);
                                         j = ceil(rand*(k0+1));
                                        

                                        L = Task(k0).Criteria(CJ).bestXs(c, j) - Xnew(i);
                                         Xnew(i) = Xnew(i) + round(F*rand *(Task(k0).Criteria(CJ).bestXs(c, j) - Xnew(i)));
                                         
                                         if  Xnew(i) > SNPs
                                             Xnew(i) = SNPs - max(0,round(normrnd(0,min(5,max(L/10,1)))));
                                         elseif Xnew(i) < 1
                                             Xnew(i) = 1 + max(0,round(normrnd(0,min(5,max(L/10,1)))));
                                         end
                                 end   
                             end
                      else
                             Xnew(i)=ceil(rand*SNPs);
                      end


                      if i==1 || ~ismember(Xnew(i),Xnew(1:i-1))
                          i = i + 1;         
                      end

                  if ( i-1 == dim ) 
                      Xtemp=Xnew;
                      Xnew=sort(Xnew);
                      while ( ismember(Xnew ,Task(k).Criteria(1).X,'rows')||ismember(Xnew ,Task(k).Criteria(2).X,'rows')||...
                              ismember(Xnew ,Task(k).Criteria(3).X,'rows')||ismember(Xnew ,Task(k).Criteria(4).X,'rows'))
                          j=ceil(rand*dim);
                          r=ceil(rand*SNPs);
                          while ismember(r,Xnew)
                              r=ceil(rand*SNPs);
                          end
                          Xnew(j)=r;
                          Xtemp=Xnew;
                          Xnew =sort(Xnew);
                      end
                  end
              end
           else %% 从不同标准的记忆库中组合优化
                %% ----------------------------------------------
              i=1;
              while i<= dim
                      if rand<HMCR  
                             a = ceil(rand*HMS); 
                             b = ceil(rand*dim); 
                             cr = ceil(rand*CriteriaNum);
                             Xnew(i) = Task(k).Criteria(cr).HM(a,b);                     
                             if rand<PAR
                                 sPar = ceil(rand * 3);
                                 switch sPar
                                     case 1                                                
                                         d = ceil(rand*HMS);
                                         j = ceil(rand*(k+1));
                                         Xnew(i) = Task(k).Criteria(cr).X(d,j);                                    
                                     case 2
                                         d = ceil(rand*HMS);
                                         c = ceil(rand*bestNum);
                                         j = ceil(rand*(k+1));
                                         Xnew(i) = Task(k).Criteria(cr).bestXs(c,j);      
                                     case 3
                                         c = ceil(rand*CriteriaNum);
                                         d = ceil(rand*HMS);
                                         j = ceil(rand*(k+1));
                                         L = Task(k).Criteria(cr).X(d,j) - Xnew(i);
                                         Xnew(i) = round(Xnew(i) + F*rand*L);
%                                          if  Xnew(i) > SNPs
%                                              Xnew(i) = SNPs - max(0,round(normrnd(0,2)));
%                                          elseif Xnew(i) < 1
%                                              Xnew(i) = 1 + max(0,round(normrnd(0,2)));
%                                          end
                                         if  Xnew(i) > SNPs
                                             Xnew(i) = SNPs - max(0,round(normrnd(0,min(5,max(L/10,1)))));
                                         elseif Xnew(i) < 1
                                             Xnew(i) = 1 + max(0,round(normrnd(0,min(5,max(L/10,1)))));
                                         end
                                 end   
                             end
                      else
                             Xnew(i)=ceil(rand*SNPs);
                      end


                      if i==1 || ~ismember(Xnew(i),Xnew(1:i-1))
                          i = i + 1;         
                      end

                  if ( i-1 == dim ) 
                      Xtemp=Xnew;
                      Xnew=sort(Xnew);
                      while ( ismember(Xnew ,Task(k).Criteria(1).X,'rows')||ismember(Xnew ,Task(k).Criteria(2).X,'rows')||...
                              ismember(Xnew ,Task(k).Criteria(3).X,'rows')||ismember(Xnew ,Task(k).Criteria(4).X,'rows'))
                          j=ceil(rand*dim);
                          r=ceil(rand*SNPs);
                          while ismember(r,Xnew)
                              r=ceil(rand*SNPs);
                          end
                          Xnew(j)=r;
                          Xtemp=Xnew;
                          Xnew =sort(Xnew);
                      end
                  end
              end
           end
       end
    
     % XnewScore = Gtest_score(data(:,Xnew),State); 
   
     XnewScores = objFunction(data(:,Xnew),State); 
      NC=NC+1;
      if k+1 == D
          DFes = DFes + 1;
      end

   %%
   for cj = 1:CriteriaNum
       if ~ismember(Xnew, Task(k).Criteria(cj).X,'rows')
              [Task(k).Criteria(cj).fworst, Idworst(cj)] = max(Task(k).Criteria(cj).Score);
              if XnewScores(cj) < Task(k).Criteria(cj).fworst
                  if cj < 4 || (cj==4 && rand < (1 - NC / max_iter))
                      Task(k).Criteria(cj).Score(Idworst(cj)) = XnewScores(cj);
                      Task(k).Criteria(cj).X(Idworst(cj),:) = Xnew;
                      Task(k).Criteria(cj).HM(Idworst(cj),:) = Xtemp;
                  end
              end
       end
   end
       
     %%  The program is terminted if the Xnew is the solution. 
    
       if (length(Xnew)==length(CX)) 
           if (Xnew == CX)
               flag = 1;
               
%                fprintf(' successfully! dim=%d, NC = %d,score=(%f,%f,%f,%f),Xbest=(%d,',dim, NC, XnewScores,Xnew(1));
%                for j = 2:dim
%                    fprintf('%d ',Xnew(j));
%                end
%                fprintf(')\n');
                  
              
           end
       end
       if flag == 1
           break;
       end
%        if mod(NC,1501)==0
%             fprintf('dim=%d, Xbest=(%d',k+1,Task(k).EliteSet(4).X(1,1));
%             for j = 2:length(Task(k).EliteSet(4).X(1,:))-1
%                fprintf(',%d',Task(k).EliteSet(4).X(1,j));              
%             end
%              fprintf(',%d),score=%12.8f\n', Task(k).EliteSet(1).X(1,end),Task(k).EliteSet(4).Scores(1));
%        end
     
   end
       if flag == 1
           break;
       end
end









