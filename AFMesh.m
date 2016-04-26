function [nodes,t,e] = AFMesh(pv,h)

% Advancing Front Mesh generator on polygon with nodes given by pv and mesh-sizing
% function h; need h > 0

% Typical choices of pv
% pv = [0,0; 1,0; 1,.5 ; .2 , .5 ; .2 .6; 1,.6; 1 , 1 ; 0,1; 0,0];
% pv = [0,0; 1,0; .5,.5 ; 1 , 1 ; 0,1; 0,0];
% pv = [0,0; 1,0 ; 1 , 1 ; 0,1; 0,0];

% Typical choices of h
% h = @(x) .01 + norm(x)/10;
% h = @(x) (1 + x(1))*(1 + x(2))/30;


tic

L = size(pv,1);
nodes = pv(1:end-1,:);

fr = 1;

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
        fr(end+1) = size(nodes,1);
        
    else
        break
    end
    end
    
    fr(end+1) = i+1;
end
fr(end) = 1;
t = [];


Fr = {fr};

while ~isempty(Fr) 

    len_fr = length(Fr);    
for ct = 1:len_fr
fr = Fr{ct};    
    
while length(fr) > 3
    
    % length(fr) 
    
    dists = [];
    for i = 1:length(fr)-1
        dists(i) = norm(nodes(fr(i),:) - nodes(fr(i+1),:));
    end
    [len, ind] = max(dists);
    
    p = [nodes(fr(ind),2)-nodes(fr(ind+1),2) , nodes(fr(ind+1),1)-nodes(fr(ind),1)];
    p = p/norm(p);
    
    x0 = .5*(nodes(fr(ind),:) + nodes(fr(ind+1),:)) + len*sqrt(3)/2*p;
    x0_eps = .5*(nodes(fr(ind),:) + nodes(fr(ind+1),:)) + 0.001*len*p;
    
    hpr = (h(x0+0.01*p) - h(x0))/0.01;
        
    dels = (len - h(x0))/(hpr + (nodes(fr(ind+1),:) - x0)*p'/len);
    
    x0 = x0 + dels*p;
    
    r0 = norm(nodes(fr(ind+1),:) - x0);
    r = r0/2;
    
    if ind == 1
        v1 = nodes(fr(end-1),:) - nodes(fr(end),:);
    else
        v1 = nodes(fr(ind-1),:) - nodes(fr(ind),:);
        
    end
    
    
    if ind == length(fr)-1
        v2 = nodes(fr(2),:) - nodes(fr(1),:);
    else v2 = nodes(fr(ind+2),:) - nodes(fr(ind+1),:);
    end
    
    v1 = v1/norm(v1); v2 = v2/norm(v2);
    
    w1 = x0 - nodes(fr(ind),:);
    w2 = x0 - nodes(fr(ind+1),:);
    w1 = w1/norm(w1); w2 = w2/norm(w2);
    
    Angs = [v1*w1' , v2*w2'];
    candists = [];
    
    
    
    while 1
        cand = [];
        k_cand = [];
        
        for i = 1:length(fr)-1
            if i == ind || i == ind+1
                continue;
            elseif ind == length(fr)-1 && i == 1
                continue;
            elseif norm(nodes(fr(i),:) - x0) < r
                cand(end+1) = fr(i);
                k_cand(end+1) = i;

                
            end
        end
        
        if isempty(cand)
            t_tent = [nodes(fr(ind),:) ; nodes(fr(ind+1),:) ; x0];
            ps = [];
            for a = 1:3
            for b = 1:3;
                if a ~= b
                    ps(end+1,:) = .9*t_tent(a,:) + .1*t_tent(b,:);
                    ps(end+1,:) = .6*t_tent(a,:) + .4*t_tent(b,:);
                end
            end
            end
            
            if [all(Angs < 0.8) ; inpolygon(x0(1),x0(2),nodes(fr,1),nodes(fr,2)) ; ...
                    inpolygon(ps(:,1),ps(:,2),nodes(fr,1),nodes(fr,2))]
               
                nodes(end+1,:) = x0;
                t(end+1,:) = [fr(ind) , fr(ind + 1) , size(nodes,1)];
                fr = [fr(1:ind) , size(nodes,1) , fr(ind+1:end)];
                break;
            else r = r*1.2;
            end
        else
            for k = 1:length(cand)
                candists(k) = norm(nodes(cand(k),:) - x0);
            end
            
            cou = 0;
            br1 = 0; br2 = 0;
            
            
        
            while br2 == 0
                
            [~, kind] = min(candists);
            t_tent = [nodes(fr(ind),:) ; nodes(fr(ind + 1),:) ; nodes(cand(kind),:) ; nodes(fr(ind),:)];
            
            ps = [];
            for a = 1:3
            for b = 1:3;
                if a ~= b
                    ps(end+1,:) = .9*t_tent(a,:) + .1*t_tent(b,:);
                    ps(end+1,:) = .6*t_tent(a,:) + .4*t_tent(b,:);
                end
            end
            end
            
            if [inpolygon(x0_eps(1),x0_eps(2),t_tent(:,1),t_tent(:,2)); ...
                    inpolygon(ps(:,1),ps(:,2),nodes(fr,1),nodes(fr,2))]
                t(end+1,:) = [fr(ind) , fr(ind + 1) ,  cand(kind)];
                br2 = 1;
            else
                candists(kind) = 100; 
                cou = cou + 1;
                 
            end
            
            if cou == length(cand)
                r = r*1.2;
                br1 = 1; br2 = 1;
            end
            
            end
            
            if br1
                continue;
            end
            
            lind = k_cand(kind);
            
            if ind == 1
                s1 = length(fr) - lind;
                s2 = lind - 2;
                if s1 < s2
                    gr = [fr(1) fr(lind:end)];
                    fr = [fr(2:lind) fr(2)];
                    
                    
                else
                    gr = [fr(2:lind) fr(2)];
                    fr = [fr(1) fr(lind:end)];
                end
            elseif ind == length(fr) - 2 || ind == length(fr) - 1 
                s1 = lind - 1;
                s2 = length(fr) - lind + 1;
                if s1 < s2
                    
                    gr = [fr(1:lind) fr(ind+1:end)];
                    fr = [fr(lind:ind) fr(lind)];
                else
                    gr = [fr(lind:ind) fr(lind)];
                    fr = [fr(1:lind) fr(ind+1:end)];
                end
            elseif lind > ind
                gr = [fr(ind+1:lind) fr(ind+1)];
                fr = [fr(1:ind) fr(lind:end)];
                
            else gr = [fr(lind:ind) fr(lind)];
                fr = [fr(1:lind) fr(ind+1:end)];
            end
            
           
            Fr{end+1} = gr;
                
            
            break;
        end
        
        
    end
    axis([-.5 1.5 -.5 1.5])   
    grid on
    tplot(nodes,t); 
    
end

Fr{ct} = fr;

end

Gr = Fr;
Fr = {};
for jj = 1:length(Gr)
    if length(Gr{jj}) == 4
        u = Gr{jj};
        
        t(end+1,:) = u(1:3);
    elseif length(Gr{jj}) > 4
        Fr{end+1} = Gr{jj};
    end
end

tplot(nodes,t); 

end

    
t = sort(t,2);
t = unique(t,'rows');


T = size(t,1);
Int_angs = zeros(T,3);
Vecs = zeros(3,2);

List = [];
while 1
for i = 1:T
    for a = 1:3
        b = mod(a,3) + 1;
        u = nodes(t(i,b),:) - nodes(t(i,a),:);
        Vecs(a,:) = u/norm(u);
    end
    
    for a = 1:3
        b = mod(a,3) + 1;
        Int_angs(i,a) = -Vecs(a,:)*Vecs(b,:)';
    end
end


i0 = -1; j0 = -1;
for i = 1:T
    if any(List == i)
        continue;
    end
    
     if any(Int_angs(i,:) < -0.5)
        [~ , a] = min(Int_angs(i,:));
        a = mod(a+1,3) + 1;
        n1 = t(i,a); n2 = t(i,mod(a,3)+1);
        m1 = t(i,mod(a+1,3)+1);
        i0 = i;
     end
     
    
end

if i0 ~= -1
for j = 1:T
    if any(t(j,:) == n1) && any(t(j,:) == n2) && j ~= i0 
        j0 = j;
        break;
    end
end

else
    i0;
    break;
end


if j0 == -1 
List(end+1) = i0;
continue;    
%     j0
%     break
end 
    
m2 = sum([t(i0,:) t(j0,:)]) - 2*n1 - 2*n2 - m1;

t(i0,:) = [m1 n1 m2];
t(j0,:) = [m1 n2 m2];

tplot(nodes,t)


end



    
e = boundary_nodes(t);
tplot(nodes,t)

toc


end













