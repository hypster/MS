%% shopSimulation
%  Simulation of a shop with fresh goods. The shop has only a single
%  employee that works 8 hours a day. Since they
%  are selling fresh goods, they will lose profitCustomerTooLate money if
%  a customer has to wait longer than max_delay. If the customer is served
%  within max_delay time, they will make profitCustomerInTime.
%  To regulate the earnings, they can choose between operating a FIFO or a
%  LIFO queue. Also they can set a maximum queue capacity and customers
%  arriving to the system, while at maximum queue capacity, will be turned
%  away, which is cost neutral.
%  Parameters
%  - queueCapacity Maximum capacity in the queue. If a customer arrives at
%    the system, while there are this many customers in the queue, the
%    arriving customer will be send away
%  - FIFO Boolean parameter describing the queue strategy. If true, the
%  queue operates FIFO, if false, the queue operates LIFO
%  Output
%  - averageProfit how much money was earned per hour
%  - averageDelay Average delay  of the customers accepted in the system
%  - fractionRejected Proportion of customers that was rejected
function [averageProfit,averageDelay,fractionRejected] = shopSimulation(queueCapacity,FIFO)

if(nargin<1)
    queueCapacity=5; % optional parameter with queue capacity
end

if(nargin<2)
    FIFO=true;
end

service_mu=1.5; % average service time in minutes
service_sigma=0.8; % standard deviation of the service time
arrival_lambda=0.60; % arrival rate of customers per minute
% System state:
serverStatus=0; % 0 is idle, 1 is busy
nQueue=0;
tArrivals=[]; % times of arrival of the customers waiting in the queue
clock=0;
tArrivalEvent=generateArrivalTime(arrival_lambda);
tDepartureEvent=inf;
% statistical counters:
nArrivals=0; % number of customers arrived at the system
nTurnedAway=0; % number of customers turned away due to a full queue
nDelay=0; %number of people who enter the service, including 0 delay, updated when customer enters into service
sDelay=0; %total duration of delay, updated when customer leaves the queue
aQueue=0;
aBusy=0;
% profit system
maxDelay=7.0; % max time in the queue that customers can wait in minutes
profitCustomerInTime=10; % how much will you earn from a customer that was served in time
profitCustomerTooLate=-4; % how much will you loose from a customer that was served too late
profit=0; % profit made

maxTime=8*60; % max time in minutes to run the simulation


%     % Display initial status
%     disp('===============================================');
%     disp(['Time: ' num2str(clock)]);
%     disp(['Server status: ' num2str(serverStatus)]);
%     disp(['Next arrival event: ' num2str(tArrivalEvent)]);
%     disp(['Next departure event: ' num2str(tDepartureEvent)]);
%     disp(['Number in queue: ' num2str(nQueue)]);
%     disp(['Arrival times in queue: ' num2str(tArrivals)]);
%     disp(['Number of customers with completed delays: ' num2str(nDelay)]);
%     disp(['Total duration of delays: ' num2str(sDelay)]);
%     disp(['Delay time so far: ' num2str(aQueue)]);
%     disp(['Busy time so far: ' num2str(aBusy)]);
%     disp('===============================================');

    while clock<maxTime
        % check whether the next event is an arrival or a departure
        prevClock=clock;
        if tArrivalEvent<tDepartureEvent
            clock=tArrivalEvent;
            nArrivals=nArrivals+1; % register an arrival

            if(nQueue<queueCapacity) % check whether we want to receive this customer

                % update waiting time
                aQueue=aQueue+(clock-prevClock)*nQueue;
                % update busy time
                aBusy=aBusy+(clock-prevClock)*serverStatus;

                % Check whether  customer can proceed to server
                if serverStatus==0
                    nDelay=nDelay+1;
                    serverStatus=1;
                    profit=profit+profitCustomerInTime; % we're making money!
                    tDepartureEvent = clock + generateServiceTime(service_mu,service_sigma); % schedule departure
                    tArrivalEvent=clock+generateArrivalTime(arrival_lambda); % schedule new arrival
                else
                % or enqueue him
                    if(FIFO)
                        tArrivals=[tArrivals clock];% store time of arrival FIFO
                    else
                        tArrivals=[clock tArrivals];% store time of arrival LIFO
                    end
                    nQueue=nQueue+1;
                    tArrivalEvent=clock+generateArrivalTime(arrival_lambda); % schedule new arrival
                end
            else
                nTurnedAway=nTurnedAway+1;
                tArrivalEvent=clock+generateArrivalTime(arrival_lambda);
            end
        else
        % we have a departure
            clock = tDepartureEvent;

            % update waiting time
            aQueue=aQueue+(clock-prevClock)*nQueue;
            % update busy time
            aBusy=aBusy+(clock-prevClock)*serverStatus;

            % take customer out of service
            serverStatus=0;
            tDepartureEvent = inf; % don't forget this one

            % check whether there is someone waiting in queue
            if nQueue>0
                % register his delay
                sDelay = sDelay + clock-tArrivals(1);
                nDelay = nDelay + 1;
                % check whether we made money
                if(clock-tArrivals(1)<=maxDelay)
                    profit=profit+profitCustomerInTime;
                else
                    profit=profit+profitCustomerTooLate;
                end
                % calculate departure
                serverStatus=1;
                tDepartureEvent = clock + generateServiceTime(service_mu,service_sigma);
                %take him out of queue
                tArrivals(1)=[];
                nQueue=nQueue-1;
            end
        end
%         % Display status
%         disp(['Time: ' num2str(clock)]);
%         disp(['Server status: ' num2str(serverStatus)]);
%         disp(['Next arrival event: ' num2str(tArrivalEvent)]);
%         disp(['Next departure event: ' num2str(tDepartureEvent)]);
%         disp(['Number in queue: ' num2str(nQueue)]);
%         disp(['Arrival times in queue: ' num2str(tArrivals)]);
%         disp(['Number of customers with completed delays: ' num2str(nDelay)]);
%         disp(['Total duration of delays: ' num2str(sDelay)]);
%         disp(['Delay time so far: ' num2str(aQueue)]);
%         disp(['Busy time so far: ' num2str(aBusy)]);
%         disp('===============================================');
    end
    
    % generate outputs
    averageProfit=60*profit/maxTime;
    averageDelay=sDelay/nDelay;
    fractionRejected=nTurnedAway/nArrivals;
end

function S = generateServiceTime(mu,sigma)
    % generate a strictly positive service time
    S=-1;
    while(S<=0)
        S=mu+randn(1)*sigma;
    end
end

function A = generateArrivalTime(lambda)
    A=(-1/lambda)*log(rand());
end
