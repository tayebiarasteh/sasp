function states=hmm2_viterb(hparm,lbtj,handles)
axes(handles.dispViterbi);

cla
hold on;
[nobs,n]=size(lbtj);
xlim([1 20])
ylim([0.7 n+0.3])
plot(handles.states,'r*')
xlabel('Time index')
ylabel('State')
set(gca,'YTick',[1 2 3])

DEL=1.0e-10;

state_est = zeros(nobs,n);
delta = zeros(nobs,n);
states=zeros(nobs,1);

%
%        Initialization
%
delta(1,:) = log(max(hparm.Pi',DEL)) +  lbtj(1,:);
state_est(1,:) = 1;


%
%    Iteration
%
for t=2:nobs,
	for  j=1:n,
		tmp = delta(t-1,:)' + log(max(hparm.A(:,j),DEL));
		[mx,imax]=max(tmp);
		delta(t,j)=mx+lbtj(t,j);
		state_est(t,j)=imax;
		plot([t-1,t],[state_est(t,j),j],'Color',[.8 .8 .8]);
		
		if t>16 && t<nobs-3
			xlim([t-16, t+3])
		elseif t>=nobs-3
			xlim([nobs-19, nobs])
		else
			xlim([1 20])
		end
		if ~get(handles.fast,'Value')
			pause(0.1)
		end
	end;
end;


%
%     termination
%
[mx,imax]=max(delta(end,:));
states(end) = imax;

if ~get(handles.fast,'Value')
	pause(1)
end

%
%     backtracking
%
for t= nobs-1 : -1 : 1,
	states(t) = state_est(t,states(t+1));
	plot([t-1,t],[states(t),states(t+1)],'b','LineWidth',2);
	if t>16 && t<nobs-3
		xlim([t-16, t+3])
	elseif t>=nobs-3
		xlim([nobs-19, nobs])
	else
		xlim([1 20])
	end
	if ~get(handles.fast,'Value')
		pause(0.3)
	end
end;
hold off
xlim([1 nobs])
return
