function[errChange]=UpdateExperiment(verbose)
%Generated a Random Erdos-Reyni graph and ran the hals 
%upates using NNLS
%Graph initialization
n = 2000;
A = rand(n,n);
A = triu(A,1);
A = A + A';
thres_conn = log(n)/n;
thres_giant = 1/(n-1);
t1 = linspace(0, thres_giant, 50);
t2 = linspace(thres_giant, thres_conn, 100);
t = [t1 ones(1,40)*thres_giant t2];
G = @(p) A < p;
A_init=G(t(2));
%NMF initializations
k=50;
updateTolerance=1e-3;
%run initial hals
[W_init,H_init]=hals(A_init,rand(n,k),rand(k,n),k,1e-6,100,0);
initGradH=norm(H_init,'fro');
updatedH=H_init;
prevError=initGradH+1000;
if(verbose==1)
    errChange=zeros(numel(t)-2,3);
end
currentIteration=1;
for (tt = t(3:end))
    %get updated graph
    currentGraph=G(tt);
    %call to one update using NNLS.    
    updatedH=nnlsm_blockpivot(W_init,currentGraph,0,updatedH); 
    %compute the currentError
    currentError=norm(updatedH,'fro')/initGradH;
    if(verbose==1)
        errChange(currentIteration,1)=(currentError-prevError);
        errChange(currentIteration,2)=norm(A-W_init*updatedH,'fro');
    end
    %if we are deviating from the previous iteration more than tolerance
    %recompute W_init and H_init
    if((currentError-prevError)>updateTolerance)
        %if this is the case we have to update W.
        if(verbose==1)
            errChange(currentIteration,3)=1;
        end
        [W_init,H_init]=hals(currentGraph,W_init,updatedH,k,1e-6,100,0);
        initGradH=norm(H_init,'fro');
        updatedH=H_init;
        prevError=initGradH+1e-2;
    end
    prevError=currentError;
    currentIteration=currentIteration+1;
end
end
