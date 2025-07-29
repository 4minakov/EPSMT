%% MT2tauk
function [tau,k] = MT2tauk(fm)

M = diag(sort(eig(fm),'descend'));

M_Iso=(1/3)*trace(M)*eye(3);

M_Dev=(M-M_Iso);

if M_Dev(2,2)>0
    k = (trace(M)/3)/(abs(trace(M)/3) - M_Dev(3,3));
    T = -2*M_Dev(2,2)/M_Dev(3,3);
elseif M_Dev(2,2)==0
    k = (trace(M)/3)/(abs(trace(M)/3) + M_Dev(1,1));
    T = 0;
elseif M_Dev(2,2)<0
    k = (trace(M)/3)/(abs(trace(M)/3) + M_Dev(1,1));
    T = 2*M_Dev(2,2)/M_Dev(1,1);
end

tau = T*(1 - abs(k));