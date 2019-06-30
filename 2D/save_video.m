%
%  
% Сохранение видео
%
v = VideoWriter('motion_fairing_test.avi');
v.Quality = 98; 
open(v);

figure;
set(gca,'FontSize',20);
axis([min(r0(1,:))-10.0 max(r0(1,:)) -15 15]);
daspect([1 1 1]);
set(gcf, 'Units','pixels' , 'OuterPosition', [10, 10, 1280, 720+79]);
hold on;
for frame = 1:2:size(t,1)
    [~, F0, Fp, ncon] = dqdt_polygon(t, q(frame,:)', params);
    A = [cos(q(frame,3)) -sin(q(frame,3));sin(q(frame,3)) cos(q(frame,3))];
    cla;

    yi=q(:,5:2:params.n*2+3);
    xi=q(:,4:2:params.n*2+3);
    [maxy,imaxy] = min(yi,[],1);
    patch(diag(xi(imaxy,:)),maxy,[0.95 0.95 0.95],'EdgeColor','w');

    text(-5,14,sprintf('t=%5.2f c | F_p^{max}=%6.2f кН | F_x^b=%5.1f кН | F_y^b=%5.1f кН',t(frame), max(sqrt(sum(Fp.*Fp,2)))*0.001,A'*F0*0.001),...
    'FontSize',18,'FontName','FixedWidth','Color','b');    
    line([-25 -15 -15],[1 1 -2],'LineWidth',4);
	plot(q(frame,4:2:params.n*2+2),q(frame,5:2:params.n*2+3),'r.-','LineWidth',3.0);
    plot(q(frame,4:2:params.n*2+2),q(frame,5:2:params.n*2+3),'w.');
    
    r0i   = q(frame,1:2);
    phi0 = q(frame,3);
    pp   = (repmat(r0i',1,size(params.pp,1)) + [cos(phi0) -sin(phi0);sin(phi0) cos(phi0)]*params.pp')';
    patch([pp(:,1);pp(1,1)],[pp(:,2);pp(1,2)],[98 141 65]/255.0,'linewidth',1.5);    

    rc = q(frame,1:2)';
    vc = q(frame,params.n*2+4:params.n*2+5)';
    quiver(rc(1),rc(2),vc(1)*0.7,vc(2)*0.7,0,'Color',[0.8,0.8,0.8],'LineWidth',2);
   

	f = getframe(gcf);
    if saveVideo 
       writeVideo(v,f); 
    end
end

close(v);
close(gcf);
