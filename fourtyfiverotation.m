%%%lossless 45� rotation MZ 2016

close all
clear oldx M
tic,
standardintp=0;
    interpmethode =  'linear'; %'spline'; %'nearest' % 
    
%data
if 0 
    uiopen('C:\Users\gi2570\Desktop\40_by_40_thumbnail_of_''Green_Sea_Shell''.png',1)
else
    uiopen('./160_by_160_thumbnail_of_''Green_Sea_Shell''.png',1);
end

try
  Green_Sea_Shell_160=x160_by_160_thumbnail_of_0x27Green_Sea_Shell0x27;
catch
   Green_Sea_Shell_160= x40_by_40_thumbnail_of_0x27Green_Sea_Shell0x27;
end

gss_org=(sum(double(Green_Sea_Shell_160(:,:,:)),3));
gss_org=interpft(interpft(gss_org,1*max(size(gss_org)),1),1*max(size(gss_org)),2);
gss_org=padarray(gss_org,[2 2], 0 ,'both'); %mean(gss_org(:)) ,'both'); %double
%gss_org=zeros(40,40);
%gss_org(2,2)=700; gss_org(38,38)=700; gss_org(2,38)=700; gss_org(38,2)=700;

if 0 %%% spatial domain
    %%gss_org=padarray(gss_org,[max(size(gss_org))/2 max(size(gss_org))/2], mean(gss_org(:)) ,'both'); %double
    gss=gss_org;

    iteration=0;
    errorl=[];
    while iteration< 100
        iteration=iteration+1

        if  mod(iteration,2)
            gss_p=padarray(gss,[max(size(gss))/2+1 max(size(gss))/2+1], 0 ,'both'); %double

            %gss_p=ifft2(fft2(gss_p),'symmetric');
            %figure; imagesc((gss_p))
            %gss_p=ifft2(fft2(gss_p).*fftshift(padarray(ones(size(gss_org)),[max(size(gss))/2 max(size(gss))/2])));

            siz=2*max(size(gss));
            if iteration==1
                siz=siz*2;
            end
            siz=siz+1;
        else
            gss_p=gss;
            siz=2*max(size(gss))+1;
        end
        gss_pi=interpft(interpft(gss_p,siz,1),siz,2);

        %figure(3); imagesc(gss_pi)

        %mask=(cat(2,triu(ones(siz/2,siz/2),1),fliplr(triu(ones(siz/2,siz/2),1))));
        %mask=cat(1,fliplr(flipud(mask)),mask);

        %checkk=repmat([0,1],[1 ceil((siz+1)^2)]);
        %checkk=(reshape(checkk(1:(siz+1)^2)',[siz+1, siz+1 ])); checkk=checkk(1:siz,1:siz);
        %figure(2); imagesc((logical(mask.*logical(checkk))))

        oldx=[]; oldy=[];  newxx=[]; newyy=[];
        for newy=1:floor(siz/2)+1, for newx=1:floor(siz/2)+1 oldx=[oldx floor(siz/2)+newx-newy+1];oldy=[oldy newy+newx-1]; newxx=[newxx newx]; newyy=[newyy newy];
            end,end
        % figure; plot3(oldx,oldy,repmat(0,[length(oldx) 1]),'.')

        iiind=sub2ind([siz,siz],oldx,oldy);
        figure(4); imagesc(reshape(gss_pi(iiind),[floor(siz/2)+1 floor(siz/2)+1]))

        gss=reshape(gss_pi(iiind),[floor(siz/2)+1 floor(siz/2)+1]);
        gss=gss(1:siz/2,1:siz/2);
        figure(5);imagesc(gss(1:2:end,1:2:end));
        figure(6);imagesc(flipud(fliplr(gss_org)));
        if mod(iteration,8)==0
            errorl=[errorl sum(sum(abs(gss(1:2:end,1:2:end)-gss_org)))];
        end

    end
    plot(errorl)

else %%try 2
    %%gss_org=padarray(gss_org,[max(size(gss_org))/2 max(size(gss_org))/2], mean(gss_org(:)) ,'both'); %double
    gss=gss_org;
 
    iteration=0;
    errorl=[];
    while iteration < 160
        iteration=iteration+1

        if iteration==1 %1 % mod(iteration,2)
            gss_p=padarray(gss,[floor(max(size(gss))/2) floor(max(size(gss))/2)], mean(gss(:)) ,'post'); %double
            gss_p=padarray(gss_p,[ceil(max(size(gss))/2) ceil(max(size(gss))/2)], mean(gss(:)) ,'pre'); %double
            %gss_p=ifft2(fft2(gss_p),'symmetric');
            %figure; imagesc((gss_p))
            %gss_p=ifft2(fft2(gss_p).*fftshift(padarray(ones(size(gss_org)),[max(size(gss))/2 max(size(gss))/2])));
            %gss_p=gss;
            siz=max(size(gss));
            %if iteration==1
            %    siz=siz*2;       
            %end
           siz
        else
            gss_p=gss;
            siz=2*max(size(gss));
            siz=siz;
        end
        %gss_pi=gss_p;
        
        %reduce again to high freQ!
        if 1 %mod(iteration,2)
            
             if standardintp
                    [x,y] = meshgrid(1:size(gss_p,1),1:size(gss_p,1));
                    [xn,yn] = meshgrid(linspace(1,size(gss_p,1),2*size(gss_p,1)),linspace(1,size(gss_p,1),2*size(gss_p,1)));
                    gss_pi=interp2(x,y,gss_p,xn,yn,interpmethode);
                    siz=max(size(gss_p));
             else
                 gss_pi=interpft(interpft(gss_p,2*size(gss_p,1),1),2*size(gss_p,1),2);
                 siz=max(size(gss_p));
             end
            %if iteration==1
            %    siz=siz*2;       
            %end
           siz
        else
           gss_pi=gss_p;
           siz=max(size(gss));
        end
        %figure(10); imagesc(gss_pi(1:4:end,1:4:end))
        %figure(3); imagesc(gss_pi)

        %mask=(cat(2,triu(ones(siz/2,siz/2),1),fliplr(triu(ones(siz/2,siz/2),1))));
        %mask=cat(1,fliplr(flipud(mask)),mask);

        %checkk=repmat([0,1],[1 ceil((siz+1)^2)]);
        %checkk=(reshape(checkk(1:(siz+1)^2)',[siz+1, siz+1 ])); checkk=checkk(1:siz,1:siz);
        %figure(2); imagesc((logical(mask.*logical(checkk))))
   
        %siz=siz/0.5
        if ~exist('oldx','var') || isempty(oldx)
            oldx=zeros(siz,siz); oldy=zeros(siz,siz); newxx=[]; newyy=[];
            for newy=1:1:floor(siz), 
                for newx=1:1:floor(siz)
                    oldx(newx,newy)=ceil(siz)+newx-newy+1;
                    oldy(newx,newy)=newy+newx-1;
                    %oldx=[oldx ceil(siz)+newx-newy+1]; oldy=[oldy newy+newx-1]; newxx=[newxx newx]; newyy=[newyy newy];
                end,
            end
            oldx=oldx(:); oldy=oldy(:);
            % figure; plot3(oldx,oldy,repmat(0,[length(oldx) 1]),'.')

            iiind=sub2ind([2*(siz),2*(siz)],oldx,oldy);
        end

            %figure(9); imagesc(abs(reshape(gss_pi(iiind(1:1:end)+0),[siz siz])))
            
            
            if mod(iteration,2)
                if standardintp
                    [x,y] = meshgrid(1:2*siz,1:2*siz);
                    gss=reshape(interp2(x,y,gss_pi,oldy,oldx,interpmethode),[siz siz]);
                    %[x,y] = meshgrid(1:(max(size(gss)))\(siz):siz,1:(max(size(gss)))\(siz):siz);
                    %gss=reshape(interp2(x,y,gss,oldy,oldx-siz/2,interpmethode),[siz siz]);
                else
                    gss=reshape(gss_pi(iiind(1:1:end)+0),[siz siz]);
                end
            else
                if standardintp
                    [x,y] = meshgrid(1:siz*2,1:siz*2);
                    gss=zeros(2*siz,2*siz);
                    gss(ceil(siz/2):ceil(siz/2)+siz-1,ceil(siz/2):ceil(siz/2)+siz-1)=((reshape(interp2(x,y,gss_pi,reshape(oldy,[siz siz]),reshape(oldx,[siz siz]),interpmethode),[siz siz])));
                 
                    %[x,y] = meshgrid(1:siz,1:siz);
                    %gss_t=zeros(2*siz,2*siz);
                    %gss_t(ceil(siz/2):ceil(siz/2)+siz-1,ceil(siz/2):ceil(siz/2)+siz-1)=((reshape(interp2(x,y,gss,reshape(oldy,[siz siz]),reshape(oldx,[siz siz]),interpmethode),[siz siz])));
                    %gss=gss_t;
                else
                    gss=zeros(2*siz,2*siz);
                    gss(ceil(siz/2):ceil(siz/2)+siz-1,ceil(siz/2):ceil(siz/2)+siz-1)=reshape(gss_pi(iiind),[siz siz]);
                end
            end
           
           if ~mod(iteration,2)
              % gss=gss(21:60,21:60); 
              gss=gss(2:2:end,2:2:end); 
           else
             %gss=gss(21:60,21:60);   
             %gss=gss(siz/2-1:siz/2+siz,siz/2-1:siz/2+siz);
           end
            %gssl=ifftshift(fft2(gss_pi));
            %gssl=circshift(fftshift(reshape(gss(iiind),[siz siz])),[0 0]);
       %     gss=zeros(1*(siz),1*(siz));
      
      % if mod(iteration,2)
      %      gss(iiind)=gssl;        
      % else
      %    gss(ceil(siz/2):ceil(siz/2)+siz-1,ceil(siz/2):ceil(siz/2)+siz-1)=reshape(gssl(iiind),[siz siz]);        
      %  end
        
        %gss=circshift(ifftshift(gss(1:2:end,1:2:end)),[0 -1]);
        %figure(2); imagesc((abs(gss(1:1:end,1:1:end))))
        
        %gssl=circshift(fftshift(reshape(gss(iiind),[siz/2 siz/2])),[0 1]);
        %gss=fftshift(ifft2(gssl)); %enforce hermetian symmetry
        %gss=fftshift(abs(ifft2((gssl),'symmetric'))./(1)); %enforce hermetian symmetry
        %gss=fftshift(ifft2(gss(1:1:end,1:1:end),'symmetric'));
        %  if iteration ==1 %mod(iteration,2)
        %      gss=circshift(gss,[0 2]);
        %  end
        
        %figure(7); imagesc(fftshift(((gss(:,:)))))
        
        % figure(4); imagesc(log(abs((gss)))./2);
        subplot(2,2,2),
          figure(4);
        imagesc(((gss))); title(['roation number ' num2str(iteration)]); colorbar;
        %figure(4); imagesc(fftshift(log(abs(fft2(gss)))));
        %gss=gss(1:siz/2,1:siz/2);
        %figure(5);imagesc(gss(1:2:end,1:2:end)); colorbar; drawnow;
        subplot(2,2,1),
        imagesc(gss_org);colorbar; title('original image')
        
        if mod(iteration,8)==0
            %%%normalize -> workaround FIXME
            %    gss=gss.*(max(gss_org(:))/max(gss(:)));
            gss_t=gss(siz/4+1:siz*3/4,siz/4+1:siz*3/4);
            % gss_t=gss_t(2:2:end,2:2:end);
            subplot(2,2,4),
            errorl=[errorl sum(sum(abs(gss_t(1:1:end,1:1:end)-gss_org)))];
            plot(errorl); title(['mean error in percent to mean value' num2str(100.*mean(mean(abs(gss_t(1:1:end,1:1:end)-gss_org)))/mean(gss_org(:)))]);

            subplot(2,2,3),
            imagesc(gss_t(1:1:end,1:1:end)-gss_org); colorbar; title('difference image')
           drawnow;
        end
        %saveas(gcf,'output','png')
        M(iteration) = getframe(gcf);

    end
toc,
movie(M)
    %plot(errorl)
end
%edit fourtyfiverotation.mcls