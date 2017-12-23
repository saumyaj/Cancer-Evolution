alphar= 0.20;
betar= 1e-10;
alphas= 0.20;
betas= 1e-10;

deltapr=0.01925;
deltacr=0.017325;
deltaps=0.01925;
deltacs=0.017325;

uprcs=1e-11;
upscr=1e-11;
uprps=1e-10;
upspr=1e-10;

i_betar = 0.5; 
i_kc = 0;
i_kp = 0.98;
i_betas = 0.5;
i_upspr = 5e-9;

global lambda;
lambda= 0.00396;
k = 5e14;

n=2000;
iteration=n+1;

Pr=zeros(1,iteration);
Ps=zeros(1,iteration);
Ps(1)=1;
Cr=zeros(1,iteration);
Cs=zeros(1,iteration);
t1 = 600;
t2 = 842;
p = 21; 
q = 2;
CycleLen = 6; %cycle length
for i=2:iteration
    Pr(i) = Pr(i-1) + ((1-alphar-uprcs-uprps)*Pr(i-1)+upspr*Ps(i-1))*f(Pr(i-1), Ps(i-1), Cr(i-1), Cs(i-1), k) + betar*Cr(i-1) - deltapr*Pr(i-1); 
    Ps(i) = Ps(i-1) + ((1-alphas - upscr - upspr)*Ps(i-1) + uprps*Pr(i-1))* f(Pr(i-1), Ps(i-1), Cr(i-1), Cs(i-1), k) + betas*Cs(i-1) - deltaps*Ps(i-1);
    Cr(i) = Cr(i-1) + (alphar*Pr(i-1) + upscr*Ps(i-1))*f(Pr(i-1), Ps(i-1), Cr(i-1), Cs(i-1), k) - betar*Cr(i-1) - deltacr*Cr(i-1) ;
    Cs(i) = Cs(i-1) + (alphas*Ps(i-1) + uprcs*Pr(i-1))*f(Pr(i-1), Ps(i-1), Cr(i-1), Cs(i-1), k) - betas*Cs(i-1) - deltacs*Cs(i-1) ;
    if ((Pr(i)+Ps(i)+Cr(i)+Cs(i))>=10e12)
        i
        break;
    end;
    if((mod(i-t1,p)==0 && ((i-t1)/p <= (CycleLen-1)) && i>=t1) || ((mod(i-t2,p)==0 && ((i-t2)/p <= (CycleLen-1)) && i>=t2)))
        Pr(i) = Pr(i) + i_upspr*Ps(i-1) + i_betar*Cr(i-1);
        Ps(i) = Ps(i) - i_upspr*Ps(i-1) + i_betas*Cs(i-1)-i_kp*Ps(i-1);
        Cr(i)=Cr(i)-i_betar*Cr(i-1);
        Cs(i)=Cs(i) - (i_betas+i_kc)*Cs(i-1);
    end;
end

plot(log10(Pr+Ps+Cr+Cs));

