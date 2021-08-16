rootdir = '/mnt/data/Mitra/cache/repos/figures/16082021';

a = dir(fullfile(rootdir,'/**/LM-*.fig'));
for i=1:length(a)
    i
    cd(a(i).folder)
    
    uiopen(fullfile(a(i).folder,a(i).name),1)
    f=gcf;
    set(f,'Color','w')
  
    %saveas(f,[a(i).name(1:end-4),'.svg'])   
    %export_fig('combined.pdf','-pdf', '-m2','-append')
    export_fig('LM-combined.pdf','-pdf','-append','-r200','-nocrop','-p10')
    %print(f,'-dpsc2','-noui','-append','-r300','combined.pdf');
    close(f);
    pause(1)
end




