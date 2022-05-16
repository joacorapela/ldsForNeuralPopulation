

for i = 1:8
    figure;plot(Xval.params{2}.model.B(8+i,2:9))
    hold on;plot(Xval.params{1}.model.B(8+i,2:9))
    hold on;plot(Xval.params{3}.model.B(8+i,2:9))
    hold on;plot(Xval.params{4}.model.B(8+i,2:9))
    hold on;plot(Xval.params{5}.model.B(8+i,2:9))
end