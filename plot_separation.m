function[curve]= plot_separation(point3,meet)
count=1;
for idx=1:length(point3)
    meet1=meet(idx,:); 
    meet2=meet(cell2mat(arrayfun(@(c) length(intersect(point3(idx,:),point3(c,:)))==2,1:length(point3),'uni',0)),:);
    for idk=1:size(meet2,1)
        p2=meet2(idk,:)'; p1=meet1'; r=norm(p1);
        v3=cross(cross(p1,p2),p1);  v3=r*v3/norm(v3);
        t = linspace(0,atan2(norm(cross(p1,p2)),dot(p1,p2)));
        v = p1*cos(t)+v3*sin(t); % v traces great circle path, relative to center
        curve{count}=v; count=count+1;
    end
    meet2=[];
end
end