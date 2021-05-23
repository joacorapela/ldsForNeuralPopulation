
do_plot = 0;

RF = cell(1,size(params.model.C,1)); % number of neurons
for n = 1:size(RF,2) % n is number of neurons
    
    for in = 1:4 % input
        for i=1:11 % 10 lags
            RF{n}(i,in) = params.model.C(n,:)*((params.model.A)^(i-1))*params.model.B(:,in);
        end
    end
    
    if do_plot
        figure;
        hold on; plot(RF{n}(:,1),'Color',[0 0.7 0 ]);% vis go
        hold on; plot(RF{n}(:,2),'Color',[0.7 0 0 ]);% vis nogo
        hold on; plot(RF{n}(:,3),'g');% laser go
        hold on; plot(RF{n}(:,4),'r');% laser nogo
        legend({'vg','vn','lg','ln'})
    end
end