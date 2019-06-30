function [dr, is_inside] = points2polygon(rp, pp)
%
% rp - точки (nx2)
% pp - n точек (nx2) вершины многоугольника с обходом против часовой стрелки

% Количество вершин многоугольника
nv = size(pp,1);
np = size(rp,1);
% Добавляем в конец списка первую вершину
pp = repmat([pp;pp(1,:)],1,1,np);
rp = permute(repmat(rp,1,1,size(pp,1)-1),[3,2,1]);
% Вектора ребер
edges = diff(pp,1,1);
% Единичные вектора ребер
ee    = edges./repmat(sqrt(sum(edges.^2,2)),1,2,1);
% Единичные вектора нормалей ребер, направленные внутрь, если вершины
% обходятся против часовой стрелке
en    = [-ee(:,2,:) ee(:,1,:)];
% От вершины до точки
d     = rp - pp(1:end-1,:,:);
% Расстояние от каждого ребра до точки
dist  = sum(d.*en,2);
% Если все больше нуля, то точка rp внутри многоугольника
is_inside = prod(dist>0,1);
% Определяем минимальное
[mind,imin] = min(dist,[],1);

imin = squeeze(imin);
mind = squeeze(mind);
is_inside = squeeze(is_inside);
en   = en(:,:,1);
% Направляем вектор наружу
dr = -repmat(mind,1,2).*en(imin,:).*repmat(is_inside,1,2);

end

