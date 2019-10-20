clear
clc

setenv('MKL_NUM_THREADS','8')
setenv('MKL_SERIAL','YES')
setenv('MKL_DYNAMIC','NO')  

%/mnt/lustre/working/lab_christing/jingleiL/EEG-fMRI


addpath(genpath('/mnt/lustre/home/jingleiL/bin/spm12'))
addpath(genpath('/mnt/lustre/home/jingleiL/bin/Tensorlab'))
addpath '/mnt/lustre/home/jingleiL/bin/eeglab13_6_5b'
addpath '/mnt/lustre/home/jingleiL/bin/eeglab13_6_5b/functions/sigprocfunc'
addpath '/mnt/lustre/home/jingleiL/bin/eeglab13_6_5b/functions/guifunc'
addpath '/mnt/lustre/home/jingleiL/bin/eeglab13_6_5b/functions/adminfunc'
addpath '/mnt/lustre/home/jingleiL/bin/eeglab13_6_5b/functions/miscfunc'
addpath '/mnt/lustre/home/jingleiL/bin/eeglab13_6_5b/functions/popfunc'
mycmd='module load fsl/5.0.9_eddy';
system(mycmd);
 addpath '/mnt/lustre/home/jingleiL/bin/spams-matlab/build'
 
load /mnt/lustre/working/lab_christing/jingleiL/EEG-fMRI/all_chan_struct_array.mat
inpath='/mnt/lustre/working/lab_christing/jingleiL/EEG-fMRI';
opath='/mnt/lustre/working/lab_christing/jingleiL/EEG-fMRI/tensorlab_Avalon/Test_4D_CPD_A2_MICCAI/A2_nnegabcd_R60_2_45_1.5';
cmpnum=60;
lambda=1.5;       % Sparsity parameter
bhrf=1;
startf=2;
fbandnum=45;
fftnum=500;
fftoverlapnum=450;


mycmd=['mkdir ',opath];
system(mycmd);
chids=1:61;
chan=allchanl(chids);
trlength=[215,530];


neegfall=[];
subnum=0;
for scnid={'A'};
    for runid=2

        for subid=[2 3 4 5 7 8 9 10 11 12 13 14 16 17 18 19 20]
            
            if runid==1
               fmricolnum=215;
               eegcolnum=236501;
            else
               fmricolnum=530;
               eegcolnum=583001;
            end

            if subid==19&&strcmp(scnid,'B')&&runid==2
                fprintf('Skip 19-B-run2.....\n');
            else
                
                fname=[inpath,'/cEEG_interpolate/',num2str(subid),cell2mat(scnid),'.run',num2str(runid),'_data.mat']
                eeg=load(fname);
                eeg=eeg.tttt;
                eeg=eeg';
                neeg=eeg;
                neegf=zeros(fbandnum,61,fmricolnum);

                for i=1:fmricolnum
                    eegseg=neeg(int32((i-1)*2.2*500+1):int32(i*500*2.2),:);
                    eegseg=eegseg-repmat(mean(eegseg),[size(eegseg,1) 1]);
                    eegseg=eegseg./repmat(std(eegseg),[size(eegseg,1) 1]);
                    [pxx,f] = pwelch(eegseg,fftnum,fftoverlapnum,fftnum,500,'onesided','psd');
                    neegf(:,:,i)=pxx(startf:startf+fbandnum-1,:);               
                end
                clear eeg neeg;
                subnum=subnum+1;
                neegfall(:,:,:,subnum)=neegf;
                
            end
        end
    end
end
        
        R=cmpnum;
        model = struct;
        model.variables.a = randn(size(neegfall,1),R);
        model.variables.b = randn(size(neegfall,2),R);
        model.variables.c = randn(size(neegfall,3),R);
        %model.variables.c = randn(I*R-0.5*R*(R-1),1);
        model.variables.d = randn(size(neegfall,4),R);
        %model.variables.e = randn(1,R);
        
        model.factors.A = {'a',@struct_nonneg};
        model.factors.B = {'b',@struct_nonneg};
        %orthq = @(z,task)struct_orth(z,task,[I R]);
        model.factors.C = {'c',@struct_nonneg};
        model.factors.D = {'d',@struct_nonneg};
        %model.factors.E = {'e',@struct_nonneg};
        model.factorizations.tensor.data =neegfall; %fmt(neegfall,'TolSparsity',2);
        model.factorizations.tensor.cpd  = {'A', 'B', 'C', 'D'};%'E'
        
        sdf_check(model, 'print');
        [sol, output] = sdf_nls(model, 'Display', 10, 'MaxIter', 1000,'CGMaxIter',200,'TolFun',eps^2,'TolX',eps);    
         fname=[opath,'/neegf.mat'];
         save(fname,'neegf');
         fname=[opath,'/sol.mat']
         save(fname,'sol');

         fname=[opath,'/R.mat']
         save(fname,'R');
         fname=[opath,'/f.mat']
         save(fname,'f');
         
%         qq=sol.factors.E;
%         fname=[opath,'/E.mat']
%              save(fname,'qq'); 
%          [Y,I]=sort(qq,'descend');
%          
        mm=sol.factors.B;
        nn=sol.factors.A;
        pp=sol.factors.D;
        dd=sol.factors.C;
%         mm=mm(:,I);
%         nn=nn(:,I);
%         pp=pp(:,I);
%         dd=dd(:,I);
       
        dd=dd-repmat(mean(dd),[size(dd,1) 1]);
        dd=dd./repmat(std(dd),[size(dd,1) 1]);
        if bhrf==1
            ndd=zeros(size(dd));

            for jj=1:size(dd,2)
                x=dd(:,jj);
                BF.dt=2.2;
                BF.name='hrf';
                BF=spm_get_bf(BF);
                U.u=x;
                U.name={'test'};
                convx=spm_Volterra(U,BF.bf);
                ndd(:,jj)=convx;
            end
        else
            ndd=dd;
        end
 
        
          fname=[opath,'/ndd.mat']
             save(fname,'ndd');
             mycmd=['mkdir ',opath,'/EVs'];
            system(mycmd);
            for i=1:R
                fname=[opath,'/EVs/ev',num2str(i),'.txt'];
                dlmwrite(fname,dd(:,i),'\t');
            end

        
       for i=1:R
           figure;
           set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,4]);
           topoplot(mm(:,i),chan);
           colorbar
           title('EEG spatial map');
           fname=[opath,'/',num2str(i),'_eeg.jpeg'];
           print(gcf,'-djpeg',fname,'-r100');
           close
           figure;
           set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,4]);
           bar(f(startf:startf+fbandnum-1),nn(:,i));
           title('Frequency distribution');
           xlabel('HZ')
           fname=[opath,'/',num2str(i),'_frequency.jpeg'];
           print(gcf,'-djpeg',fname,'-r100');
           close
           figure;
           set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,4]);
           trs=1:fmricolnum;
           trs=trs.*2.2;
           plot(trs,ndd(:,i));
           title('Temporal Weights (TRs)');
           xlabel('s')
           fname=[opath,'/',num2str(i),'_temperal.jpeg'];
           print(gcf,'-djpeg',fname,'-r100');
           close
           figure;
           set(gcf,'PaperUnits','inches','PaperPosition',[0,0,4,4]);
           bar([1:subnum],pp(:,i));
           title('Subject distribution');
           xlabel('subject')
           fname=[opath,'/',num2str(i),'_subject.jpeg'];
           print(gcf,'-djpeg',fname,'-r100');
           close
       end
       
            fname=[opath,'/trs.mat']
            save(fname,'trs');
   

for scnid={'A'};
    for runid=2

      
        for subid=[2 3 4 5 7 8 9 10 11 12 13 14 16 17 18 19 20]
                       
            if runid==1
               fmricolnum=215;
               eegcolnum=236501;
            else
               fmricolnum=530;
               eegcolnum=583001;
            end

            if subid==19&&strcmp(scnid,'B')&&runid==2
                fprintf('Skip 19-B-run2.....\n');
            else
                fname=[inpath,'/FMRI_signals/',num2str(subid),cell2mat(scnid),'_run',num2str(runid),'_3mm_retrend_whole_b_signals.txt']
                fid=fopen(fname);
                fmri=textscan(fid,[repmat('%f',[1,fmricolnum])]);
                fclose(fid);
                fmri=cell2mat(fmri);
                fmri=fmri';
                fmri=fmri-repmat(mean(fmri),[size(fmri,1) 1]);
                fmri=fmri./repmat(std(fmri),[size(fmri,1) 1]);


                nopath=[opath,'/',num2str(subid),cell2mat(scnid),'_run',num2str(runid)];
                mycmd=['mkdir ',nopath];
                system(mycmd);

                fprintf('Starting mexLasso.....\n');
                param.mode=2;
                param.lambda=lambda;
                %param.pos=true;
                param.lambda2=0;
                param.numThreads=-1;
                tic
                A=mexLasso(fmri,ndd,param);
                t=toc;
                fprintf('time of computation for sparsecoding: %f\n',t)
                A=full(A);
                A(isnan(A))=0;

                fname=[nopath,'/fmricomp.txt'];
                dlmwrite(fname,A,'\t');

                obase=[nopath,'/test']
                mycmd=['/mnt/lustre/home/jingleiL/bin/Map_reference2brain /mnt/lustre/working/lab_christing/jingleiL/EEG-fMRI/MNI152_3mm_mask.nii.gz ',fname,' ',num2str(R),' ',obase]
                system(mycmd);

                mycmd=['cp /mnt/lustre/working/lab_christing/jingleiL/EEG-fMRI/tensorlab_Avalon/Test_frequency_hrf/map_network.sh ',nopath,'/']
                system(mycmd);

                mycmd=['chmod 755 ',nopath,'/map_network.sh']
                system(mycmd); 

                mycmd=[nopath,'/map_network.sh ',num2str(R),' ',nopath,...
                    ' ',num2str((median(median(abs(A)))+1.5*abs(median(median(abs(A)))-mean(mean(abs(A)))))),' ',num2str((median(median(abs(A)))+0.5*abs(max(max(abs(A)))-median(median(abs(A))))))];
                system(mycmd); 

                mycmd=['cp /mnt/lustre/working/lab_christing/jingleiL/EEG-fMRI/tensorlab_Avalon/Test_frequency_hrf/cmps.html ',nopath,'/'];
                system(mycmd);
                
            end
        end       
    end
end


