function [nodes,t,e] = DelMesh(pv,h)

% Delaunay Mesh generator on polygon with nodes given by pv and mesh-sizing
% function h; need h > 0

% Typical choices of pv
% pv = [0,0; 1,0; 1,.5 ; .2 , .5 ; .2 .6; 1,.6; 1 , 1 ; 0,1; 0,0];
% pv = [0,0; 1,0; .5,.5 ; 1 , 1 ; 0,1; 0,0];

% Typical choices of h
% h = @(x) .01+norm(x)/10;
% h = @(x) (1 + x(1))*(1 + x(2))/30;


tic

bnd = size(pv,1)-1;
Proj = zeros(2,2,bnd);
for k = 1:bnd
    diff = pv(k+1,:) - pv(k,:);
    Proj(:,:,k) = diff'*diff/norm(diff)^2 - eye(2);
end


L = size(pv,1);
nodes = pv(1:end-1,:);

for i = 1:L-1
    x0 = pv(i,:);
    x1 = pv(i+1,:);
    
    while 1
    p = x1 - x0; 
    dis = norm(p);
    p = p/dis;
    
    hpr = (h(x0+0.01*p) - h(x0))/0.01;
    delx = h(x0)/(1 - hpr);
    
    if 2*delx < dis
        x0 = x0 + delx*p;
        nodes(end+1,:) = x0;
    else
        break
    end
    end
    
end


triangles = delaunayn(nodes);
t = [];
for i = 1:size(triangles,1)
    p = sum(nodes(triangles(i,:),:))/3;
    if inpolygon(p(1),p(2),pv(:,1),pv(:,2))
        t(end+1,:) = triangles(i,:);
    end
end



count = 0;

while 1

tplot(nodes,t); 
% count = count + 1;

bnd_dists = zeros(1,bnd);

S = 0; j = 0;
for i = 1:size(t,1)
    A = sum([norm(nodes(t(i,1),:) - nodes(t(i,2),:)) , norm(nodes(t(i,1),:) - nodes(t(i,3),:)) , ...
        norm(nodes(t(i,3),:) - nodes(t(i,2),:))])/sum([h(nodes(t(i,1),:)) , h(nodes(t(i,2),:)) , ...
        h(nodes(t(i,3),:))]);
    
    if A > S 
        S = A;
        j = i;
    end
end

if S > 1
    g = size(nodes,1);
    [cent,r] = CCentre(nodes(t(j,1),:),nodes(t(j,2),:),nodes(t(j,3),:));
    
    for bs = 1:bnd
        bnd_dists(bs) = norm(Proj(:,:,bs)*(cent - pv(bs,:))');
    end
    
    bnd_dists = bnd_dists/r;
    
    if inpolygon(cent(1),cent(2),pv(:,1),pv(:,2)) && all(bnd_dists > .5)
        
    else        
        
        dists = [norm(nodes(t(j,1),:) - nodes(t(j,2),:)) , norm(nodes(t(j,2),:) - nodes(t(j,3),:)) , ...
        norm(nodes(t(j,3),:) - nodes(t(j,1),:))];
        
        hvals = [.5*(h(nodes(t(j,1),:)) + h(nodes(t(j,2),:))) , ...
        .5*(h(nodes(t(j,2),:)) + h(nodes(t(j,3),:))) , ...
        .5*(h(nodes(t(j,3),:)) + h(nodes(t(j,1),:)))];
    
        eval = dists./hvals;

        [~,ind1] = max(eval);
        
        ind2 = mod(ind1,3) + 1;
        
        cent = nodes(t(j,ind1),:) + nodes(t(j,ind2),:);
        cent = .5*cent;
    
    end
    
    nodes = [nodes; cent];
    
    
else
    break
end


t_old = [];
nodes_new = [];

for k = 1:size(t,1)
    [centk, rk] = CCentre(nodes(t(k,1),:),nodes(t(k,2),:),nodes(t(k,3),:));
    
    if norm(cent - centk) < rk
        nodes_new(end+1:end+3,:) = t(k,:)';
    else t_old(end+1,:) = t(k,:);
    end
end



nodes_new = unique(nodes_new); 
nodes_new(end+1) = g+1;


nodes_del = nodes(nodes_new,:);

tr_n = delaunayn(nodes_del);


for l = 1:size(tr_n,1)
    tr_n(l,:) = nodes_new(tr_n(l,:))';
end

[gon_x, gon_y] = hul(nodes_del(:,1),nodes_del(:,2));

triangles = [t_old ; tr_n];


t = t_old;
for i = size(t_old,1) + 1:size(triangles,1)
    
    p = [];
    for a = 1:3
        for b = 1:3;
            if a ~= b
                p(end+1,:) = .9*nodes(triangles(i,a),:) + .1*nodes(triangles(i,b),:);
            end
        end
    end

    if [inpolygon(p(:,1),p(:,2),gon_x,gon_y) ; inpolygon(p(:,1),p(:,2),pv(:,1),pv(:,2))]
        t(end+1,:) = triangles(i,:);
    end
end


end


e = boundary_nodes(t);

tplot(nodes,t)

toc

end















