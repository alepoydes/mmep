dir='fields/HRBVW2/';
d=dlmread([dir,'batch.txt']);
b=0.02:0.004:0.116;
e=d-d(:,end);

figure(1);
plot(linspace(0,1,size(d,2)),d','-');
ylabel('Energy (J)');
xlabel('Reaction coordinate');
title('Energy along MEP from skyrmion to fermionic state');
grid on
mark=@(n,q) text(size(q,2),q(n,end),['H=',num2str(b(n)),'(J)'],'horizontalalignment','right');
mark(1,d);
mark(13,d);
mark(size(d,1),d);
print([dir,'figure1.pdf'],'-dpdf');

figure(2);
plot(linspace(0,1,size(d,2)),a','-');
ylabel('E-E_FM (J)');
xlabel('Reaction coordinate');
title('Energy along MEP from skyrmion to fermionic state');
grid on
mark=@(n,q) text(1,q(n,1),['H=',num2str(b(n)),'(J)'],'horizontalalignment','left');
mark(1,a);
mark(13,a);
mark(size(a,1),a);
print([dir,'figure2.pdf'],'-dpdf');

figure(3);
plot(b,a(:,1),'-',b,max(a,[],2),'-');
ylabel('E-E_FM (J)');
xlabel('Magnetic field H (J)');
title('Energy of skyrmion and barrier as functions of the magnetic field');
legend('Skyrmion','Barrier','location','southeast');
grid on
print([dir,'figure3.pdf'],'-dpdf');