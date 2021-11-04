function[meet]= intersectingPoints(samples,point3,fillers,prototypes,threshold1,threshold2)
occupied=zeros(length(point3),1); 
meet=zeros(length(point3),3);
count=1; %beta=1;
%ca2d = @(cosa,beta) (exp(-beta*cosa+beta)-1)/(exp(2*beta)-1);
ca2P = @(cosa,beta) (exp( beta*cosa+beta)-1)/(exp(2*beta)-1);
%finding points of intersection of more than 2 domains
for idx=1:samples
           [val1, ind1] = sort(ca2P(fillers(idx,:)*prototypes',1)./sum(ca2P(fillers(idx,:)*prototypes',1),2),'descend');
           if (abs(val1(1)-val1(2))<=threshold1 & abs(val1(1)-val1(3))<=threshold2)
             if occupied(find(sum(unique(ind1(1:3))==(point3),2)==3))==0
                meet(find(sum(unique(ind1(1:3))==point3,2)==3),:)=fillers(idx,:);
                %meet(count,:)=fillers(idx,:);
                occupied(find(sum(unique(ind1(1:3))==point3,2)==3))=1;
               count=count+1;
             else
                 continue;
             end
           end
end


