clear, clc
n = 100; % # vertices
m = 20; % # of samples
K = 200; % # trials
Q = rand(n); % probability matrix
P=Q;
P(P<0.9)=0;

%%

for k=1:K
    A=nan(n,n,m);
    for i=1:m
        for u=1:n
            for v=1:n
                A(u,v,i)=rand<P(u,v);
            end
        end
    end
    
    %%
    Abar=mean(A,3);
    Abars=Abar;
    Abars(Abars<0.1)=0;
    
    loss(k)=sum(sum((Abar-P).^2));
    loss_chance(k)=sum(sum((rand(n)-P).^2));
    losss(k)=sum(sum((Abars-P).^2));
    loss0(k)=sum(sum((zeros(n)-P).^2));
    
end
%%
risk=mean(loss)
risks=mean(losss)
risk_chance=mean(loss_chance)
risk0=mean(loss0)

%%
figure(1), clf
imagesc(P), colorbar, colormap('gray')

figure(2), clf
imagesc(A(:,:,1)), colormap('gray')
title('sample graph')

figure(3), clf
imagesc((A(:,:,1)-P).^2), colorbar, colormap('gray')
title('A-P')

figure(4),
subplot(121), imagesc(Abar)
subplot(122), imagesc((Abar-P).^2), colormap('gray'), colorbar
title('Abar-P')


[foo, bar]=hist(loss);
[foo2, bar2]=hist(loss_chance);
[foos, bars]=hist(losss);
[foo0, bar0]=hist(loss0);

figure(5),clf, hold on
plot(bar,foo,'r')
plot(bar2,foo2,'k')
plot(bars,foos,'b')
plot(bar0,foo0,'b')
legend('Abar','chance','sparse','zero')
set(gca,'xscale','log')

figure(6), clf
imagesc(Abar-Abars), colorbar

figure(7), clf
imagesc(Abars<0.05)
