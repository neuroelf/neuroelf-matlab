function df = multistat(input_files_ef, input_files_sd, ...
   df_data, fwhm_data, X, contrast, output_file_base, ...
   which_stats, fwhm_varatio, niter, df_limit)

%MULTISTAT fits a mixed effects linear model.
%
% Combines effects (E) and their standard errors (S) using a linear mixed 
% effects model:     E = X b + e_fixed + e_random,     where
%    b is a vector of unknown coefficients,
%    e_fixed  is normal with mean zero, standard deviation S,
%    e_random is normal with mean zero, standard deviation sigma (unknown).
% The model is fitted by REML using the EM algorithm with NITER iterations. 
%
% Gives the conjunction (minimum) of the random effects T stats for the data.
%
% Estimates the effective FWHM in mm of the residuals from the linear model,
% as if the residuals were white noise smoothed with a Gaussian filter 
% whose fwhm was FWHM. Applies a first-order correction for small degrees 
% of freedom and small FWHM relative to the voxel size. Great care has been
% taken to make sure that FWHM is unbiased, particularly for low df, so that
% if it is smoothed spatially then it remains unbiased. The bias in FWHM
% is less than 1/10 of the voxel size for FWHM > 2*voxel size. However the 
% continuity correction for small FWHM breaks down if FWHM < voxel size,
% in which case FWHM is generally too large. If FWHM > 50, FWHM = 50.
%
% WARNING: Multistat is very slow if the number of columns of X is more than 
% 1 and less than the number of input files, and INPUT_FILES_SD is not empty,
% since it loops over voxels, rather than doing calculations in parallel.
% 
% DF = MULTISTAT( INPUT_FILES_EF , INPUT_FILES_SD , DF_DATA [, FWHM_DATA 
%         [, X [, CONTRAST [, OUTPUT_FILE_BASE [, WHICH_STATS 
%           [, FWHM_VARATIO [, NITER [, DF_LIMIT ]]]]]]]] )
% 
% INPUT_FILES_EF is the input fmri effect files, the dependent variables,
% usually the _ef.img or _ef.mnc files, padded with extra blanks if necessary; 
% they will be removed before use.
% 
% INPUT_FILES_SD is the input fmri sdeffect files, the standard
% deviations of the dependent variables, usually the _sd.img or _sd.mnc files,
% padded with extra blanks if necessary; they will be removed before use.
% If INPUT_FILES_SD=[], then INPUT_FILES_SD is assumed to be
% zero for all voxels, DF_DATA is set to Inf, and FWHM_VARATIO now
% smoothes the voxel sd. This allows multistat to duplicate DOT for 
% analysing PET data (with smoothed voxel sd, rather than pooled sd).
% 
% DF_DATA is the row vector of degrees of freedom of the input files.
% If empty (default), these are read from the headers of INPUT_FILES_SD.
%
% FWHM_DATA is the fwhm in mm of INPUT_FILES_EF. It is only  
% used to calculate the degrees of freedom, printed out at the beginning. 
% If empty (default), it is estimated from the least-squares residuals,
% or it is read from the headers of INPUT_FILES_EF if available.
%
% X is the design matrix, whose rows are the files, and columns
% are the explanatory (independent) variables of interest. 
% Default is X=[1; 1; 1; ..1] which just averages the files. If the
% rank of X equals the number of files, e.g. if X is square, then 
% the random effects cannot be estinmated, but the fixed effects
% sd's can be used for the standard error. This is done very quickly.
% 
% CONTRAST is a matrix whose rows are contrasts for the statistic images.
% Default is [1 0 ... 0], i.e. it picks out the first column of X.
% 
% OUTPUT_FILE_BASE: matrix whose rows are the base for output statistics,
% one base for each row of CONTRAST, padded with extra blanks if 
% necessary; they will be removed before use.
%
% WHICH_STATS: character matrix inidicating which statistics for output,
% one row for each row of CONTRAST. If only one row is supplied, it is used 
% for all contrasts. The statistics are indicated by strings, which can
% be anywhere in WHICH_STATS, and outputed in OUTPUT_FILE_BASEstring.ext, 
% depending on the extension of INPUT_FILES_EF. The strings are: 
% _t       T statistic image =ef/sd.
% _ef      effect (b) image for magnitudes.
% _sd      standard deviation of the effect for magnitudes. 
% _rfx     ratio of random to fixed effects standard deviation. Note that
%          rfx^2 is the F statistic for testing for random effects.
% _conj    not recommended - use conjunction.m instead. 
%          conjunction (minimum) of the T statistics for the data, i.e. 
%          min(INPUT_FILES_EF/sd) using a mixed effects sd. 
% _resid   the residuals from the model, only for non-excluded frames.
% _wresid  the whitened residuals from the model normalized by dividing
%          by their root sum of squares, only for non-excluded frames.
% _fwhm    FWHM information.
%          Frame 1: effective FWHM in mm of the whitened residuals,
%          as if they were white noise smoothed with a Gaussian filter 
%          whose fwhm was FWHM. FWHM is unbiased so that if it is smoothed  
%          spatially then it remains unbiased. If FWHM > 50, FWHM = 50.
%          Frame 2: resels per voxel, again unbiased.
%          Frames 3,4,5: correlation of adjacent resids in x,y,z directions.
% e.g. WHICH_STATS='try this: _t _ef _sd and_fwhm blablabla' 
% will output t, ef, sd and fwhm.
% You can still use 1 and 0's for backwards compatiability with previous 
% versions - see help from previous versions. 
% If empty (default), only DF is returned. The strings are: 
% 
% FWHM_VARATIO is the fwhm in mm of the Gaussian filter used to smooth the
% ratio of the random effects variance divided by the fixed effects variance.
%  - 0 will do no smoothing, and give a purely random effects analysis;
%  - Inf will do complete smoothing to a global ratio of one, giving a 
%    purely fixed effects analysis. 
% The higher the FWHM_VARATIO, the higher the ultimate degrees of
% freedom DF of the tstat image (printed out at the beginning of the run), 
% and the more sensitive the test. However too much smoothing will
% bias the results. The program prints and returns DF as its value.
% Alternatively, if FWHM_VARATIO is negative, it is taken as the desired
% df, and the fwhm is chosen to get as close to this as possible (if fwhm>50,  
% fwhm=Inf). Default is -100, i.e. the fwhm is chosen to achieve 100 df. 
%
% NITER is the number of iterations of the EM algorithm. Default is 10.
%
% DF_LIMIT controls which method is used for estimating FWHM. If DF_RESID > 
% DF_LIMIT, then the FWHM is calculated assuming the Gaussian filter is 
% arbitrary. However if DF is small, this gives inaccurate results, so
% if DF_RESID <= DF_LIMIT, the FWHM is calculated assuming that the axes of
% the Gaussian filter are aligned with the x, y and z axes of the data. 
% Default is 4. 
%
% DF.t is the effective mixed-effects df of the sd and T statistics.
% DF.resid is the df of a random effects analysis.
% DF.fixed is the df of a fixed effects analysis.
% DF.rfx is the numerator df of _rfx (denominator is DF.fixed).
% Note that DF.resid <= DF.t <= DF.fixed. 

%############################################################################
% COPYRIGHT:   Copyright 2002 K.J. Worsley,
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center, 
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada. 
%              worsley@math.mcgill.ca, liao@math.mcgill.ca
%
%              Permission to use, copy, modify, and distribute this
%              software and its documentation for any purpose and without
%              fee is hereby granted, provided that the above copyright
%              notice appear in all copies.  The author and McGill University
%              make no representations about the suitability of this
%              software for any purpose.  It is provided "as is" without
%              express or implied warranty.
%############################################################################

% Defaults:

if nargin<3; df_data=[]; end
if nargin<4; fwhm_data=[]; end
if nargin<5; X=[]; end
if nargin<6; contrast=[]; end
if nargin<7; output_file_base=[]; end
if nargin<8; which_stats=[]; end
if nargin<9; fwhm_varatio=[]; end
if nargin<10; niter=[]; end
if nargin<11; df_limit=[]; end

parent_file=deblank(input_files_ef(1,:));
if isempty(X); X=ones(size(input_files_ef,1),1); end
if isempty(contrast); contrast=[1 zeros(1,size(X,2)-1)]; end
if isempty(output_file_base); output_file_base=parent_file(1:(max(findstr('.',parent_file))-1)); end
if isempty(fwhm_varatio); fwhm_varatio=-100; end
if isempty(niter); niter=10; end
if isempty(df_limit); df_limit=4; end

% Open images:

n=size(input_files_ef,1)
d=fmris_read_image(parent_file);
d.dim
numslices=d.dim(3);
numys=d.dim(2);
numxs=d.dim(1);
numpix=numys*numxs;
Steps=d.vox
numcolX=size(contrast,2)
D=2+(numslices>1)

% Open files for writing:

[base,ext]=fileparts2(input_files_ef(1,:));
out.parent_file=parent_file;
numcontrasts=size(contrast,1)

if ~isempty(which_stats)
   if size(which_stats,1)==1
      which_stats=repmat(which_stats,numcontrasts,1);
   end
   if isnumeric(which_stats)
      which_stats=[which_stats zeros(numcontrasts,9-size(which_stats,2))];
   else
      ws=which_stats;
      which_stats=zeros(numcontrasts,9);
      fst=['_t      ';'_ef     ';'_sd     ';'_rfx    ';'_conj   ';'_resid  ';'_wresid ';'not used';'_fwhm   '];
      for i=1:numcontrasts
         for j=1:3
            which_stats(i,j)= ~isempty(strfind(ws(i,:),deblank(fst(j,:)))); 
         end
      end
      for j=4:9
         which_stats(1,j)= ~isempty(strfind(ws(1,:),deblank(fst(j,:))));
      end
   end
   which_stats
end

X2=X'.^2;
p=rank(X)
df.resid=n-p;

if isempty(input_files_sd)
   df_data=Inf
else
   if isempty(df_data)
      for ifile=1:n
         d=fmris_read_image(input_files_sd(ifile,:),0,0);
         if isfield(d,'df')
            df_data(ifile)=d.df;
         end
      end
   end
end
if length(df_data)==1
   df_data=ones(1,n)*df_data;
end
df.fixed=sum(df_data);

if isempty(which_stats)
   df_data
   df
   return
end

if df.resid>0
   
   % Degrees of freedom is greater than zero, so do mixed effects analysis:
   
   Y=zeros(n,numpix);
   S=ones(n,numpix);
   varfix=ones(1,numpix);
   varatio_vol=zeros(numpix, numslices);
   
   if fwhm_varatio<0
      % find fwhm to achieve target df:
      df_target=-fwhm_varatio;
      if df_target<=df.resid
         fwhm_varatio=0
      elseif df_target>=df.fixed
         fwhm_varatio=Inf
      end
   end
   
   if fwhm_varatio<Inf
      
      % Now loop over voxels to get variance ratio, varatio:
      
      Sreduction=0.99
      for slice=1:numslices
         First_pass_slice=slice            
         for ifile=1:n
            d=fmris_read_image(input_files_ef(ifile,:),slice,1);
            Y(ifile,:)=reshape(d.data,1,numpix);
         end
         sigma2=sum((Y-X*(pinv(X)*Y)).^2,1)/df.resid;
         
         if ~isempty(input_files_sd)
            for ifile=1:n
               d=fmris_read_image(input_files_sd(ifile,:),slice,1);
               S(ifile,:)=reshape(d.data,1,numpix).^2;
            end
            varfix=df_data*S/df.fixed;  
         end
         
         if isempty(fwhm_data) & (fwhm_varatio~=0)
            % Quick fwhm of input files:
            if ~isempty(input_files_sd)
               sdd=(varfix>0)./sqrt(varfix*df.resid+(varfix<=0));
            else
               sdd=(sigma2>0)./sqrt(sigma2*df.resid+(sigma2<=0));
            end
            wresid_slice=((Y-X*(pinv(X)*Y)).*repmat(sdd,n,1))';
            if slice==1
               D=2+(numslices>1);
               sumr=0;
               i1=1:(numxs-1);
               j1=1:(numys-1);
               nxy=conv2(ones(numxs-1,numys-1),ones(2));
               u=reshape(wresid_slice,numxs,numys,n);
               if D==2
                  ux=diff(u(:,j1,:),1,1);
                  uy=diff(u(i1,:,:),1,2);
                  axx=sum(ux.^2,3);
                  ayy=sum(uy.^2,3);
                  if df.resid>df_limit
                     axy=sum(ux.*uy,3);
                     detlam=(axx.*ayy-axy.^2);
                  else
                     detlam=axx.*ayy;
                  end
                  r=conv2((detlam>0).*sqrt(detlam+(detlam<=0)),ones(2))./nxy;
               else
                  r=zeros(numxs,numys);
               end
               
               if ~isempty(input_files_sd)
                  mask=reshape(sqrt(varfix),numxs,numys);
                  Df=df.fixed;
               else
                  mask=reshape(sqrt(sigma2),numxs,numys);
                  Df=df.resid;
               end
               tot=sum(sum(mask));
               
               alphar=1/2;
               dr=df.resid/Df;
               dv=df.resid-dr-(0:D-1);
               if df.resid>df_limit
                  biasr=exp(sum(gammaln(dv/2+alphar)-gammaln(dv/2)) ...
                     +gammaln(Df/2-D*alphar)-gammaln(Df/2))*dr^(-D*alphar);
               else
                  biasr=exp((gammaln(dv(1)/2+alphar)-gammaln(dv(1)/2))*D+ ...
                     +gammaln(Df/2-D*alphar)-gammaln(Df/2))*dr^(-D*alphar);
               end

            else 
               uz=reshape(wresid_slice,numxs,numys,n)-u;
               ux=diff(u(:,j1,:),1,1);
               uy=diff(u(i1,:,:),1,2);
               uz=uz(i1,j1,:);
               axx=sum(ux.^2,3);
               ayy=sum(uy.^2,3);
               azz=sum(uz.^2,3);
               if df.resid>df_limit
                  axy=sum(ux.*uy,3);
                  axz=sum(ux.*uz,3);
                  ayz=sum(uy.*uz,3);
                  detlam=(axx.*ayy-axy.^2).*azz-(axz.*ayy-2*axy.*ayz).*axz-axx.*ayz.^2;
               else
                  detlam=axx.*ayy.*azz;
               end
               mask1=mask;
               if ~isempty(input_files_sd)
                  mask=reshape(sqrt(varfix),numxs,numys);
               else
                  mask=reshape(sqrt(sigma2),numxs,numys);
               end
               tot=tot+sum(sum(mask));
               r1=r;
               r=conv2((detlam>0).*sqrt(detlam+(detlam<=0)),ones(2))./nxy;
               sumr=sumr+sum(sum((r1+r)/(1+(slice>2)).*mask1));
               u=reshape(wresid_slice,numxs,numys,n);
            end
            if slice==numslices
               sumr=sumr+sum(sum(r.*mask));
               fwhm_data=sqrt(4*log(2))*(prod(abs(Steps(1:D)))*tot/sumr*biasr)^(1/3)
               clear u ux uy uz wresid_slice 
            end
         end
         
         if ~isempty(input_files_sd)
            minS=min(S)*Sreduction;
            Sm=S-repmat(minS,n,1);
            if size(X,2)==1
               % When X is a vector, calculations can be done in parallel:
               for iter=1:niter
                  Sms=Sm+repmat(sigma2,n,1);
                  W=(Sms>0)./(Sms+(Sms<=0));
                  X2W=X2*W;
                  XWXinv=(X2W>0)./(X2W+(X2W<=0));
                  betahat=XWXinv.*(X'*(W.*Y));
                  R=W.*(Y-X*betahat);
                  ptrS=p+sum(Sm.*W,1)-(X2*(Sm.*W.^2)).*XWXinv;
                  sigma2=(sigma2.*ptrS+(sigma2.^2).*sum(R.^2,1))/n; 
               end
               sigma2=sigma2-minS;
            else
               % Otherwise when X is a matrix, we have to loop over voxels:
               for pix=1:numpix
                  sigma2_pix=sigma2(pix);
                  Sm_pix=Sm(:,pix);
                  Y_pix=Y(:,pix);
                  for iter=1:niter
                     Sms=Sm_pix+sigma2_pix;
                     W=(Sms>0)./(Sms+(Sms<=0));
                     Whalf=diag(sqrt(W));
                     WhalfX=Whalf*X;
                     pinvX=pinv(WhalfX);
                     WhalfY=Whalf*Y_pix;
                     betahat=pinvX*WhalfY;
                     R=WhalfY-WhalfX*betahat;
                     SW=diag(Sm_pix.*W);
                     ptrS=p+sum(Sm_pix.*W,1)-sum(sum((SW*WhalfX).*pinvX'));   
                     sigma2_pix=(sigma2_pix.*ptrS+ ...
                        (sigma2_pix.^2).*sum(W.*(R.^2),1))/n; 
                  end
                  sigma2(pix)=sigma2_pix-minS(pix);
               end         
            end
         end
         varatio_vol(:,slice)=(sigma2./(varfix+(varfix<=0)).*(varfix>0))';  
         
      end
      
      if fwhm_varatio<0
         % find fwhm to achieve target df:
         df_target=-fwhm_varatio;
         fwhm_varatio_limit=50;
         ff=0:1:100;
         dfs=[];
         for f=ff;
            df=regularized_df(f,D,Steps,numslices,df,fwhm_data);
            dfs=[dfs df.t];
         end
         fwhm_varatio=interp1(dfs,ff,df_target);
         if isnan(fwhm_varatio) | fwhm_varatio > fwhm_varatio_limit
            fwhm_varatio=Inf;
            varatio_vol=zeros(numpix,numslices);
         end
         fwhm_varatio
      end
      
      [df, ker_x, ker_y, K]=regularized_df(fwhm_varatio,D,Steps,numslices,df,fwhm_data);
      df.rfx=round(df.rfx);
      df.t=round(df.t);
      
      if fwhm_varatio>0 & fwhm_varatio<Inf
         % Smoothing varatio in slice is done using conv2 with a kernel ker_xy.
         for slice=1:numslices
            varatio_slice=reshape(varatio_vol(:,slice),numxs,numys);
            varatio_vol(:,slice)=reshape(conv2(ker_x,ker_y,varatio_slice,'same'),numpix,1);
         end
         % Smoothing betwen slices is done by straight matrix multiplication
         % by a toeplitz matrix K normalized so that the column sums are 1.
         varatio_vol=varatio_vol*K;
      end
      
   else
      df.rfx=Inf;
      df.t=df.fixed;
   end
   df_data
   df
   
   pinvX=pinv(X);
   ncpinvX=sqrt(sum((contrast*pinvX).^2,2));
   
   if isempty(fwhm_data)
      for ifile=1:n
         d=fmris_read_image(input_files_ef(ifile,:),0,0);
         if isfield(d,'fwhm')
            fwhm_data(ifile)=d.fwhm;
         end
      end
      if ~isempty(fwhm_data)
         fwhm_data=mean(fwhm_data(fwhm_data>0))
      end
   end
   if ~isempty(fwhm_data)
      out.fwhm=fwhm_data;
   end
   
   % Second loop over slices to get statistics:
   
   for slice=1:numslices
      Second_pass_slice=slice            
      for ifile=1:n
         d=fmris_read_image(input_files_ef(ifile,:),slice,1);
         Y(ifile,:)=reshape(d.data,1,numpix);
      end
      if ~isempty(input_files_sd)
         for ifile=1:n
            d=fmris_read_image(input_files_sd(ifile,:),slice,1);
            S(ifile,:)=reshape(d.data,1,numpix).^2;
         end
         varfix=df_data*S/df.fixed;   
         sigma2=varatio_vol(:,slice)'.*varfix;
         Sigma=S+repmat(sigma2,n,1);
         % To ensure that Sigma is always positive:
         if fwhm_varatio>0
            Sigma=max(Sigma,S*0.25);
         end
         W=(Sigma>0)./(Sigma+(Sigma<=0));
         if size(X,2)==1
            % Do calculations in parallel:
            X2W=X2*W;
            XWXinv=(X2W>0)./(X2W+(X2W<=0));
            betahat=XWXinv.*(X'*(W.*Y));
            sdeffect_slice=contrast*sqrt(XWXinv);
         else
            % Loop over pixels:
            betahat=zeros(size(X,2),numpix);
            sdeffect_slice=zeros(numcontrasts,numpix);
            for pix=1:numpix
               Whalf=diag(sqrt(W(:,pix)));
               pinvX=pinv(Whalf*X);
               betahat(:,pix)=pinvX*Whalf*Y(:,pix);
               sdeffect_slice(:,pix)=sqrt(sum((contrast*pinvX).^2,2));
            end  
         end
      else
         betahat=pinvX*Y;
         sigma2=varatio_vol(:,slice)';            
         sdeffect_slice=ncpinvX*sqrt(sigma2);
         Sigma=repmat(sigma2,n,1);
         W=(Sigma>0)./(Sigma+(Sigma<=0));
      end
      effect_slice=contrast*betahat;
      tstat_slice=effect_slice./(sdeffect_slice+(sdeffect_slice<=0)) ...
         .*(sdeffect_slice>0);
      sdratio_slice=sqrt(varatio_vol(:,slice)+1);
      
      % output:
      
      for k=1:numcontrasts
         if which_stats(k,1)
            out.data=[]; if isfield(out,'nconj'); out=rmfield(out,'nconj'); end
            out.file_name=[deblank(output_file_base(k,:)) '_t' ext];
            out.dim=[numxs numys numslices 1];
            out.data=reshape(tstat_slice(k,:),numxs,numys);
            out.df=df.t;
            fmris_write_image(out,slice,1);
         end
         if which_stats(k,2)
            out.data=[]; if isfield(out,'nconj'); out=rmfield(out,'nconj'); end
            out.file_name=[deblank(output_file_base(k,:)) '_ef' ext];
            out.dim=[numxs numys numslices 1];
            out.data=reshape(effect_slice(k,:),numxs,numys);
            out.df=n;
            fmris_write_image(out,slice,1);
         end
         if which_stats(k,3)
            out.data=[]; if isfield(out,'nconj'); out=rmfield(out,'nconj'); end
            out.file_name=[deblank(output_file_base(k,:)) '_sd' ext];
            out.dim=[numxs numys numslices 1];
            out.data=reshape(sdeffect_slice(k,:),numxs,numys);
            out.df=df.t;
            fmris_write_image(out,slice,1);
         end
      end
      if which_stats(1,4)
         out.data=[]; if isfield(out,'nconj'); out=rmfield(out,'nconj'); end
         out.file_name=[deblank(output_file_base(1,:)) '_rfx' ext];
         out.dim=[numxs numys numslices 1];
         out.data=reshape(sdratio_slice,numxs,numys);
         out.df=[df.rfx df.fixed];
         if isfield(out,'fwhm'); out.fwhm=[out.fwhm fwhm_varatio]; end
         fmris_write_image(out,slice,1);
         if isfield(out,'fwhm'); out.fwhm=out.fwhm(1); end
      end
      if which_stats(1,5)
         out.file_name=[deblank(output_file_base(1,:)) '_conj' ext];
         out.dim=[numxs numys numslices 1];
         out.data=reshape(min(sqrt(W).*Y,[],1),numxs,numys);
         out.df=df.t;
         out.nconj=n;
         fmris_write_image(out,slice,1);
      end
      if which_stats(1,6)
         out.data=[]; if isfield(out,'nconj'); out=rmfield(out,'nconj'); end
         out.file_name=[deblank(output_file_base(1,:)) '_resid' ext];
         out.dim=[numxs numys numslices n];
         out.data=reshape((Y-X*betahat)',numxs,numys,n);
         out.df=df.resid;
         fmris_write_image(out,slice,1:n);
      end  
      if which_stats(1,7) | which_stats(1,9)
         wresid_slice=(sqrt(W/df.resid).*(Y-X*betahat))';
      end
      if which_stats(1,7) 
         out.data=[]; if isfield(out,'nconj'); out=rmfield(out,'nconj'); end
         out.file_name=[deblank(output_file_base(1,:)) '_wresid' ext];
         out.dim=[numxs numys numslices n];
         out.data=reshape(wresid_slice,numxs,numys,n);
         out.df=[df.resid df.t];
         fmris_write_image(out,slice,1:n);
      end  
      
      if which_stats(1,9)
         
         % Finds an estimate of the fwhm for each of the 8 cube corners surrounding
         % a voxel, then averages. 
         
         if slice==1
            
            % setup for estimating the FWHM:
            I=numxs;
            J=numys;
            IJ=I*J;
            Im=I-1;
            Jm=J-1;
            nx=conv2(ones(Im,J),ones(2,1));
            ny=conv2(ones(I,Jm),ones(1,2));
            nxy=conv2(ones(Im,Jm),ones(2));
            f=zeros(I,J);
            r=zeros(I,J);
            Azz=zeros(I,J);
            ip=[0 1 0 1];
            jp=[0 0 1 1];
            is=[1 -1  1 -1];
            js=[1  1 -1 -1];
            D=2+(numslices>1);
            alphaf=-1/(2*D);
            alphar=1/2;
            Step=abs(prod(Steps(1:D)))^(1/D);
            Df=df.t;
            dr=df.resid/Df;
            dv=df.resid-dr-(0:D-1);
            if df.resid>df_limit
               % constants for arbitrary filter method:
               biasf=exp(sum(gammaln(dv/2+alphaf)-gammaln(dv/2)) ...
                  +gammaln(Df/2-D*alphaf)-gammaln(Df/2))*dr^(-D*alphaf);
               biasr=exp(sum(gammaln(dv/2+alphar)-gammaln(dv/2)) ...
                  +gammaln(Df/2-D*alphar)-gammaln(Df/2))*dr^(-D*alphar);
            else
               % constants for filter aligned with axes method:
               biasf=exp((gammaln(dv(1)/2+alphaf)-gammaln(dv(1)/2))*D+ ...
                  +gammaln(Df/2-D*alphaf)-gammaln(Df/2))*dr^(-D*alphaf);
               biasr=exp((gammaln(dv(1)/2+alphar)-gammaln(dv(1)/2))*D+ ...
                  +gammaln(Df/2-D*alphar)-gammaln(Df/2))*dr^(-D*alphar);
            end
            consf=(4*log(2))^(-D*alphaf)/biasf*Step;
            consr=(4*log(2))^(-D*alphar)/biasr;
            
            u=reshape(wresid_slice,I,J,n);
            ux=diff(u,1,1);
            uy=diff(u,1,2);
            Axx=sum(ux.^2,3);
            Ayy=sum(uy.^2,3);
            dxx=([Axx; zeros(1,J)]+[zeros(1,J); Axx])./nx;
            dyy=([Ayy  zeros(I,1)]+[zeros(I,1)  Ayy])./ny;
            if D==2
               for index=1:4
                  i=(1:Im)+ip(index);
                  j=(1:Jm)+jp(index);
                  axx=Axx(:,j);
                  ayy=Ayy(i,:);
                  if df.resid>df_limit
                     axy=sum(ux(:,j,:).*uy(i,:,:),3)*is(index)*js(index);
                     detlam=(axx.*ayy-axy.^2);
                  else
                     detlam=axx.*ayy;
                  end
                  f(i,j)=f(i,j)+(detlam>0).*(detlam+(detlam<=0)).^alphaf;
                  r(i,j)=r(i,j)+(detlam>0).*(detlam+(detlam<=0)).^alphar;
               end
            end
         else 
            uz=reshape(wresid_slice,I,J,n)-u;
            dzz=Azz;
            Azz=sum(uz.^2,3);
            dzz=(dzz+Azz)/(1+(slice>1));
            % The 4 upper cube corners:
            for index=1:4
               i=(1:Im)+ip(index);
               j=(1:Jm)+jp(index);
               axx=Axx(:,j);
               ayy=Ayy(i,:);
               azz=Azz(i,j);
               if df.resid>df_limit
                  axy=sum(ux(:,j,:).*uy(i,:,:),3)*is(index)*js(index);
                  axz=sum(ux(:,j,:).*uz(i,j,:),3)*is(index);
                  ayz=sum(uy(i,:,:).*uz(i,j,:),3)*js(index);
                  detlam=(axx.*ayy-axy.^2).*azz-(axz.*ayy-2*axy.*ayz).*axz-axx.*ayz.^2;
               else
                  detlam=axx.*ayy.*azz;
               end
               f(i,j)=f(i,j)+(detlam>0).*(detlam+(detlam<=0)).^alphaf;
               r(i,j)=r(i,j)+(detlam>0).*(detlam+(detlam<=0)).^alphar;
            end
            f=consf/((slice>2)+1)*f./nxy;
            r=consr/((slice>2)+1)*r./nxy;
            out.data=[]; if isfield(out,'nconj'); out=rmfield(out,'nconj'); end
            out.file_name=[deblank(output_file_base(1,:)) '_fwhm' ext];
            out.dim=[numxs numys numslices 5];
            out.data=f.*(f<50)+50*(f>=50);
            out.df=[df.resid df.t];
            fmris_write_image(out,slice-1,1);
            out.data=r;
            fmris_write_image(out,slice-1,2);
            out.data=1-dxx/2;
            fmris_write_image(out,slice-1,3);
            out.data=1-dyy/2;
            fmris_write_image(out,slice-1,4);
            out.data=1-dzz/2;
            fmris_write_image(out,slice-1,5);
            
            f=zeros(I,J);
            r=zeros(I,J);
            u=reshape(wresid_slice,I,J,n);
            ux=diff(u,1,1);
            uy=diff(u,1,2);
            Axx=sum(ux.^2,3);
            Ayy=sum(uy.^2,3);
            dxx=([Axx; zeros(1,J)]+[zeros(1,J); Axx])./nx;
            dyy=([Ayy  zeros(I,1)]+[zeros(I,1)  Ayy])./ny;
            % The 4 lower cube corners:
            for index=1:4
               i=(1:Im)+ip(index);
               j=(1:Jm)+jp(index);
               axx=Axx(:,j);
               ayy=Ayy(i,:);
               azz=Azz(i,j);
               if df.resid>df_limit
                  axy=sum(ux(:,j,:).*uy(i,:,:),3)*is(index)*js(index);
                  axz=-sum(ux(:,j,:).*uz(i,j,:),3)*is(index);
                  ayz=-sum(uy(i,:,:).*uz(i,j,:),3)*js(index);
                  detlam=(axx.*ayy-axy.^2).*azz-(axz.*ayy-2*axy.*ayz).*axz-axx.*ayz.^2;
               else
                  detlam=axx.*ayy.*azz;
               end
               f(i,j)=f(i,j)+(detlam>0).*(detlam+(detlam<=0)).^alphaf;
               r(i,j)=r(i,j)+(detlam>0).*(detlam+(detlam<=0)).^alphar;
            end
         end
         if slice==numslices
            f=consf*f./nxy;
            r=consr*r./nxy;
            out.data=f.*(f<50)+50*(f>=50);
            fmris_write_image(out,slice,1);
            out.data=r;
            fmris_write_image(out,slice,2);
            out.data=1-dxx/2;
            fmris_write_image(out,slice,3);
            out.data=1-dyy/2;
            fmris_write_image(out,slice,4);
            out.data=1-Azz/2;
            fmris_write_image(out,slice,5);
         end
      end % of if which_stats(1,9)
      
   end
else
   % If degrees of freedom is zero, estimate effects by least squares,
   % and use the standard errors to estimate the sdeffect.
   
   for ifile=1:n
      d=fmris_read_image(input_files_ef(ifile,:),0,0);
      if isfield(d,'fwhm')
         fwhm_data(ifile)=d.fwhm;
      end
   end
   if ~isempty(fwhm_data)
      fwhm_data=mean(fwhm_data(fwhm_data>0))
      out.fwhm=fwhm_data;
   end

   Y=zeros(n,numpix);
   S=ones(n,numpix);
   cXinv=contrast*pinv(X);
   cXinv2=cXinv.^2;    
   fprintf('Note: df.resid=0, so doing a fixed effects analysis.');
   if ~isempty(df_data)
      df.t=round(sum(cXinv2,2).^2./((cXinv2.^2)*(1./df_data')))
   else
      df.t=Inf
   end
   
   if isempty(which_stats)
      return
   end
   
   for slice=1:numslices
      slice            
      for ifile=1:n
         d=fmris_read_image(input_files_ef(ifile,:),slice,1);
         Y(ifile,:)=reshape(d.data,1,numpix);
      end
      if ~isempty(input_files_sd)
         for ifile=1:n
            d=fmris_read_image(input_files_sd(ifile,:),slice,1);
            S(ifile,:)=reshape(d.data,1,numpix).^2;
         end
      end
      effect_slice=cXinv*Y;
      sdeffect_slice=sqrt(cXinv2*S);
      tstat_slice=effect_slice./(sdeffect_slice+(sdeffect_slice<=0)) ...
         .*(sdeffect_slice>0);
      for k=1:numcontrasts
         if which_stats(k,1)
            out.file_name=[deblank(output_file_base(k,:)) '_t' ext];
            out.dim=[numxs numys numslices 1];
            out.data=reshape(tstat_slice(k,:),numxs,numys);
            out.df=df.t;
            fmris_write_image(out,slice,1);
         end
         if which_stats(k,2)
            out.file_name=[deblank(output_file_base(k,:)) '_ef' ext];
            out.dim=[numxs numys numslices 1];
            out.data=reshape(effect_slice(k,:),numxs,numys);
            out.df=n;
            fmris_write_image(out,slice,1);
         end
         if which_stats(k,3)
            out.file_name=[deblank(output_file_base(k,:)) '_sd' ext];
            out.dim=[numxs numys numslices 1];
            out.data=reshape(sdeffect_slice(k,:),numxs,numys);
            out.df=df.t;
            fmris_write_image(out,slice,1);
         end
      end
   end
   
end

return



function [df, ker_x, ker_y, K]=regularized_df(fwhm_varatio,D,Steps,numslices,df,fwhm_data);

if fwhm_varatio>0
   fwhm_x=fwhm_varatio/abs(Steps(1));
   ker_x=exp(-(-ceil(fwhm_x):ceil(fwhm_x)).^2*4*log(2)/fwhm_x^2);
   ker_x=ker_x/sum(ker_x);
   fwhm_y=fwhm_varatio/abs(Steps(2));
   ker_y=exp(-(-ceil(fwhm_y):ceil(fwhm_y)).^2*4*log(2)/fwhm_y^2);
   ker_y=ker_y/sum(ker_y);
   fwhm_z=fwhm_varatio/abs(Steps(3));
   ker_z=exp(-(0:(numslices-1)).^2*4*log(2)/fwhm_z^2);
   K=toeplitz(ker_z);
   K=K./(ones(numslices)*K);
   df_indep=df.resid/(sum(ker_x.^2)*sum(ker_y.^2)*sum(sum(K.^2))/numslices);
   df_correl=df.resid*(2*(fwhm_varatio/fwhm_data)^2+1)^(D/2);
   df.rfx=min(df_indep,df_correl);
   df.t=1/(1/df.rfx+1/df.fixed);
else
   ker_x=1;
   ker_y=1;
   K=eye(numslices);
   df.rfx=df.resid;
   df.t=df.resid;
end

return

