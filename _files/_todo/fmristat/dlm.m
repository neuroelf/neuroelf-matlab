function pdlm=dlm(ts,fwhm_file,df,mask_file,mask_thresh,nconj,nvar);

%DLM finds P-values using the discrete local maxima method.
%
% PDLM = DLM( TS , FWHM_FILE [, DF [, MASK_FILE [, MASK_THRESH 
%                [, NCONJ [, NVAR ]]]]])
%
% TS is a vector of statistic values, the first one must be > 1.
%
% For other parameters see STAT_SUMMARY. 
%
% If FWHM_FILE is a scalar or 3-vector of fwhms, they must be in
% voxels, not mm.
%
% If MASK_FILE is a scalar, it is the number of voxels in the mask. 
%
% PDLM is the row vector of DLM P-values.

if nargin<3;  df=Inf; end
if nargin<4;  mask_file=[];  end
if nargin<5;  mask_thresh=[]; end
if nargin<6;  nconj=1;  end
if nargin<7;  nvar=1;  end

pt=stat_threshold(0,1,0,df,ts,[],[],nconj,nvar);
zs=stat_threshold(0,1,0,Inf,pt);
ecz=stat_threshold([0 0 0 1],-Inf,1,Inf,zs);
ect=stat_threshold([0 0 0 1],-Inf,1,df,ts,[],[],nconj,nvar);
rs=ect./ecz;
nzs=length(zs);
rhobar=zeros(nzs,3);

if ~isempty(mask_file) & isempty(mask_thresh)
   mask_thresh=fmri_mask_thresh(mask_file);
end
if isstr(mask_file)
   mask_thresh1=mask_thresh(1);
   if length(mask_thresh)>=2
      mask_thresh2=mask_thresh(2);
   else
      mask_thresh2=Inf;
   end
   d=fmris_read_image(mask_file,0,0);
   d=fmris_read_image(mask_file,1:d.dim(3),1);
   mask=d.data;
   mask=double(mask>mask_thresh1 & mask<=mask_thresh2);
end
if isempty(mask_file) & isstr(fwhm_file)
   d=fmris_read_image(fwhm_file,0,0);
   mask=ones(d.dim(1),d.dim(2),d.dim(3));
end

if ~isempty(mask)
   vs=find(mask);
   N=length(vs);
   nbrs=[];
   for k=1:d.dim(3)
      n0=mask(:,:,k);
      n1=conv2(n0,ones(3,1),'same');
      n2=conv2(n0,ones(1,3),'same');
      n3=n0;
      if k>1
         n3=n3+mask(:,:,k-1);
      end
      if k<d.dim(3)
         n3=n3+mask(:,:,k+1);
      end
      nbr=1+n1+n2*4+n3*16;
      nbrs=[nbrs; nbr(find(n0))];
   end
   nn=hist(nbrs,1:64);
else
   nn=[zeros(1,64-length(mask_file)) mask_file];
end

if isstr(fwhm_file)
   irhos=unique(round(((1:1000)-0.5)/1000*(N-1)))+1;
   for l=1:3
      rho_vol=fmris_read_image(fwhm_file,1:d.dim(3),l+2);
      rhos=sort(rho_vol.data(vs));
      rho=rhos(irhos);
      for iz=1:nzs
         rhobar(iz,l)=1-(mean(sqrt((1-abs(rho).^rs(iz).*sign(rho))),1)).^2;
      end
   end
else
   fwhms=repmat(fwhm_file,1,4-length(fwhm_file));
   rho=1-2*log(2)./fwhms.^2;
   rho=rho.*(rho>0);
   for iz=1:nzs
      rhobar(iz,:)=rho.^rs(iz);
   end
end
   
m=200;
w=((1:m)-0.5)/m;
nz=200;
pz=((1:nz)-0.5)/nz*0.001;
z=stat_threshold(0,1,0,Inf,pz);

Q1=zeros(nzs,3,nz);
Q2=zeros(nzs,3,nz);
for irho=1:nzs
   for l=1:3
      rho=rhobar(irho,l);
      hz=sqrt((1-rho)/(1+rho))*z;
      y0=exp(-(z/(1+rho)).^2);
      y=w'*y0;
      s=-2*log(y+(y==0)/2);
      f=(y>0)./(s.*sqrt(s-ones(m,1)*hz.^2));
      Q1(irho,l,:)=1-erfc(hz/sqrt(2))/2;
      Q2(irho,l,:)=1-erfc(hz/sqrt(2))+mean(f).*(y0.*hz)/pi;
   end
end

dpz=1/nz*0.001;
pzz=(0:nz)/nz*0.001;
pdlm=zeros(1,nzs);
for iz=1:length(zs)
   PQs=0;
   for ni=1:64
      if nn(ni)>0
         [n1,n2,n3]=ind2sub([4 4 4],ni);
         ns=[n1 n2 n3]-1;
         PQ=1;
         for l=1:3
            if ns(l)==2
               PQ=PQ.*Q1(iz,l,:);
            end
            if ns(l)==3
               PQ=PQ.*Q2(iz,l,:);
            end
         end
         PQs=PQs+PQ*nn(ni);
      end
   end
   pdlmv=[0 cumsum(PQs(:))'*dpz];
   pdlm(iz)=interp1(pzz,pdlmv,pt(iz));
end

return
