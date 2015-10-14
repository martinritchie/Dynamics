function [SIT] = gillespie_SIS(A, initial, tau, gamma, dt, Tend)
N = length(A);
NN = 1:N;
state = zeros(1,N);
S(1) = N-initial;
I(1) = initial;
T(1) = 0;
for i = 1:initial;
    state(i)=1;
end
state = state(randperm(N));
t=2;
rate_vector = zeros(1,N);
for i = 1:N
    if state(i)==0
        neigh = find(A(i,:));
        for j = 1:length(neigh)
            if state(neigh(j))==1
                rate_vector(i) = rate_vector(i) + tau;
            end
        end
    end
    if state(i) == 1
        rate_vector(i)= gamma;
    end
end

time = 0;
while time<=Tend+0.5
    I(t) = I(t-1);
    S(t) = S(t-1);
    rate = sum(rate_vector);
    Cum=cumsum(rate_vector);
    indexr = 1:length(rate_vector); 
    if rate>0.000001
        tstep = log(1-rand)/(-rate);
        T(t) = T(t-1)+tstep;
        
        Event=indexr(Cum>rand*rate);
        Event = Event(1);
        temp = NN(A(Event,:)>0);
        s_contacts = temp(state(temp)==0);
        switch state(Event)
            case 0
                S(t) = S(t)-1;
                I(t) = I(t)+1;
                state(Event) = 1;
                rate_vector(Event) = gamma;
                rate_vector(s_contacts) =  rate_vector(s_contacts) + tau;
            case 1
                I(t) = I(t)-1;
                S(t) = S(t)+1;
                state(Event) = 0;
                rate_vector(Event) = tau*length(find(state(temp)==1));
                rate_vector(s_contacts) =  rate_vector(s_contacts) - tau;
        end
        time = T(t);
        t=t+1;
    else
        time = T(t-1);
        while time<= Tend + 0.5
            S(t)=S(t-1);
            I(t)=I(t-1);
            T(t)=T(t-1)+0.5;
            time = T(t);
            t = t+1;
        end
        break
    end
end


%%%%%%%%%%%%%%%%
% interpolation
M = Tend/dt + 1;
%T, I, S, R, motif1tevol
Ti=linspace(0,Tend,M);
Si=zeros(1,M);
Ii=zeros(1,M);
indexv = 1:length(T);
for jj=1:M
    k = T<=Ti(jj);
    k = indexv(k);
    k = max(k); 
%     k=find(T<=(jj-1)*dt,1,'last');
    Si(jj)=S(k);
    Ii(jj)=I(k);
    Ti(jj)=T(k);
end
SIT = [Si; Ii; Ti];
