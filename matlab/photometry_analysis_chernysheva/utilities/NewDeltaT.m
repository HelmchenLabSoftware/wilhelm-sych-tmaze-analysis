function [dt,sm_dt,sd_dt]=NewDeltaT(Timevector1,Timevector2)
%time difference between time points

dt=[];
for k=1:size(Timevector1)
    dt(k)=Timevector2(k)-Timevector1(k);
end
sm_dt=median(dt);
sd_dt=var(dt,1);
end
