function rec_cg = TrendSeasonalFit_COLDLine_ori(dir_l,n_rst,ncols,nrows,T_cg,Tmax_cg,conse,num_c,nbands,B_detect)
% CCDC-COLD 13.4 version - Zhe Zhu, University of Connecticut, Storrs
% Continuous Change Detection (CCD) designed for Continuous mOnitoring of Land Distubance (COLD)
%% Revisions: $ Date: 9/03/2020 $ Copyright: Zhe Zhu
%  Version 13.04  Update model for every three observations (09/03/2020)
%  Version 13.03  Modified model intitializatin test (10/29/2018)
%  Version 13.02  Add included angle to exlcude false positive change (10/05/2018)
%  Version 13.01  Do not need clear observations more than 25 percent (03/28/2018)
%% Version 13.0:  Optimized parameters for monitoring disturbance - COLD (03/20/2018) 
%  Version 12.36  Reduce 75 percent of memory needed for line_t (11/22/2017)
%  Version 12.35  Adjust T_cg based on conse (10/03/2017)
%  Version 12.34  Add data density requirement (09/20/2017)
%  Version 12.33  Adjust conse based on delta time (09/19/2017)
%  Version 12.32  Fix bud for not enough data (09/17/2017)
%  Version 12.31  Read input varibles in txt file (07/05/2016)
%  Version 12.30  Fixed a bug for pixels without minimum observations (11/20/2015)
%  Version 12.29  Modified fit for perennial snow and Fmask failed pixels (09/20/2015)
%  Version 12.29  Do not fit disturbed time period (09/18/2015)
%  Version 12.28  Fixed a bug for missing values in land cover maps (09/16/2015)
%  Version 12.27  Fixed bugs for persistent snow and falied Fmask pixels (06/17/2015)
%  Version 12.26  Connected time for all models (05/28/2015)
%  Version 12.25  Bug fixed in snow percent (05/19/2015)
%  Version 12.24  Change T_const in Tmask (03/31/2015)
%  Version 12.23  Update iteratively before 24 observations (03/22/2015)
%  Version 12.22  Adjust mini RMSE based on temporal variability (03/22/2015)
%  Version 12.21  Add more categories and update i_start in the end (03/14/2015)
%  Version 12.20  Convert BT to from K to C before analysis (03/12/2015)
%  Version 12.19  Fit for permanent snow if is more than 75% (03/12/2015)
%  Version 12.18  No change detection if clear observation less than 25% (03/12/2015)
%  Version 12.17  Use median value for very simple model & change magnitude (02/24/2015)
%  Version 12.16  Finding changes in all water pixels (02/24/2015)
%  Version 12.15  Use the original multitemporal cloud mask (02/15/2015)
%  Version 12.14  Do not need example_img in images folder (02/09/2015)
%  Version 12.13: More infromation in "category" (11/10/2014)
%  This version (12.13) is used for the third round of the LCMS project.
%  Command: TrendSeasonalFit_v12Plot(N_row,N_col,min=0.5,T_cg=0.99,n_times=3,conse=6,B_detect=2:6)
%  Version 12.12: Fit for pixels where Fmask fails (11/09/2014)
%  Version 12.11: Bug fixed in num_fc (11/09/2014)
%  Version 12.10: Better multietmporal cloud detection at the beginning (11/06/2014)
%  Version 12.9:  Detect change only for land pixel (water/snow speical case) (10/31/2014)
%  Version 12.8:  Speed up by reducing time for RMSE and model computing (10/17/2014)
%  Version 12.7:  mini rmse should be larger than 10% of the mean (10/13/2014)
%  Version 12.6:  Fit model again when there are a 33.3% more data (10/08/2014)
%  Version 12.5:  Use subset of bands (2-6) for detecting surface change (10/01/2014)
%  Version 12.4:  Only apply multitemporal cloud masking during model initialization (09/29/2014)
%  Version 12.3:  Use subset of bands (3-5) to balance change in diferent dimensions (09/01/2014)
%  This version (12.3) is used for the second round of the LCMS project.
%  Command: TrendSeasonalFit_v12Plot(N_row,N_col,min=1,T_cg=0.99,n_times=3,conse=5,B_detect=3:6)
%  Version 12.2:  Bug fixed in model intialization (08/14/2014)
%  Version 12.1:  Use subset of bands (3-6) to avoid atmosphere influences (08/04/2014)
%% Version 12.0   Detecting change based on probability (07/19/2014)
%  Version 11.6:  No need to change folder name & faster in speed (by Christ Holden 06/06/2014)
%  Version 11.5:  Improved calculation of temporally adjusted RMSE (04/23/2014)
%  Version 11.4:  Revise "rec_cg.category" to better seperate different fit processes (04/01/2014)
%  This version (11.4) is used for generating synthetic data for ACRE project and
%  detecting change for LCMS project.
%  Command: TrendSeasonalFit_v11Plot(N_row,N_col,min=1,T_cg=2,n_times=3,conse=6,B_detect=1:6)
%  Version 11.3:  Add "rec_cg.magnitude" as indicator of change magnitude (04/01/2014)
%  Version 11.2:  Change very simple fit with mean value for start and end of timeseries (04/01/2014)
%  Version 11.1:  Do not need metadata in the image folder to run CCDC (03/25/2014)
%% Version 11.0:  Use change vector magnitude as threshold for detecting change (03/25/2014)
%  Version 10.13: Use subset of bands (1-6) to avoid atmosphere influences (01/31/2014)
%  Version 10.12: More accurate number of days per year "num_yrs" (01/30/2014)
%  Version 10.11: RMSE up% agriculture activity (9)dates with time series fit (01/26/2014)
%  Version 10.10: Update temperature extreme in recent studies (01/16/2014)
%  Version 10.9:  Find break in max value in any of the band (01/08/2014)
%  Version 10.8:  Add very simple fit with median value for start and end of timeseries (10/21/2013)
%  This version (10.8) is used for generating synthetic data for the LCMS project.
%  Command: TrendSeasonalFit_v10Plot('stack',N_row,N_col,mini=0.5,T_cg=3,n_times=3,conse=6,B_detect=2:6)
%  Version 10.7:  Better multitemporal cloud detection (10/19/2013)
%  Version 10.6:  Add "Tmax_cg" for last step noise removal (10/18/2013)
%  Version 10.5:  Use subset of bands (2-6) to avoid atmosphere influences (10/18/2013)
%  Version 10.4:  Let dynamic fitting for pixels at the beginning (09/23/2013)
%  Version 10.3:  Able to detect change at the verying beginning (09/06/2013)
%  Version 10.2:  Add mini years "mini_yrs" in model intialization (09/03/2013)
%  Version 10.1:  Reduce time for calcuating "v_dif" (09/02/2013)
%% Version 10.0:  Fit for beginning and end of the time series (08/31/2013)
%  Version 9.9:   Only fit more than 50% of Landat images overlap area (08/28/2013)
%  Version 9.8:   Force model fit for persistent snow pixels (08/27/2013)
%  Version 9.7:   Add "rec_cg.category" as indicator of fitting procudure (08/20/2013)
%                 Add rec_cg.change_prob as indicator of change probability (08/20/2013)
%                 Add rec_cg.num_obs ad indicator of number of observations (08/20/2013)
%  Version 9.6:   Remove mininum rmse "mini" and minimum years "mini_yrs" (08/16/2013)
%  Version 9.5:   Model gets more coefficients with more observations (08/16/2013)
%  Version 9.4:   Bug fixed in calculating temporally adjusted rmse (08/01/2013)
%  Version 9.3:   Fit curve again after one year (03/28/2013)
%  This version (9.3) is used for mapping land cover for the IDS project.
%  Command: TrendSeasonalFit_v9Plot('stack',N_row,N_col,T_cg=2,n_times=3,conse=4)
%  Version 9.2:   Use "mini = T_const/T_cg" for small rmse cases (03/26/2013)
%  Version 9.1:   Remove out of range pixels before time series analysis (02/09/2013)
%% Version 9.0:   Using 8 coefficients and lasso fit (02/01/2013)
%  Version 8.4:   Use "max v_slope" instead of "average v_slope" (01/16/2013)
%  Version 8.3:   Start initialization when "time_span" > 1 year (01/16/2013)
%  Version 8.2:   Bug fix% agriculture activity (9)ed in not fitting models at the begining (01/16/2013)
%  Version 8.1:   Bug fixed in counting "i" and "i_span"(01/13/2013)
%% Version 8.0:   Temporally changing RMSE (01/09/2013)
%% Version 7.3:   Continuous Change Detection and Classification (CCDC) (07/11/2012)
%  This version (7.3) is explained by Zhu, Z. & Woodcock, C.E., Continuous Change
%  Detection and Classification (CCDC) of land cover using all available
%  Landsat data, Remote Sensing of Environment (2014).
%  Command: TrendSeasonalFit_v7Plot('stack',N_row,N_col,T_cg=3,n_times=3,conse=3)
%% Version 1.0:   Continous Monitoring of Forest Disturbance Algorithm (CMFDA) (07/13/2010)
%  This version (1.0) is explained by Zhu, Z., Woodcock, C.E., Olofsson, P.,
%  Continuous monitoring of forest disturbance using all available Landsat
%  data, Remote Sensing of Environment (2012).
%
%% Inputs:
% stk_n='stack'; stack image name
% ncols = 8021; % number of pixels processed per line
% nrows=1; % the nrowsth lines
% for example    1 2 3 4 5
%                6 7 8 9 10
%
%% Outputs:
%
% rec_cg RECord information about all curves between ChanGes
% rec_cg(i).t_start record the start of the ith curve fitting (julian_date)
% rec_cg(i).t_end record the end of the ith curve fitting (julian_date)
% rec_cg(i).t_break record the first observed break time (julian_date)
% rec_cg(i).coefs record the coefficients of the ith curve
% rec_cg(i).pos record the position of the ith pixel (pixel id)
% rec_cg(i).magnitude record the change vector of all spectral bands
% rec_cg(i).category record what fitting procudure and model is used
% cateogry category 5x: persistent snow    4x: Fmask fails
% cateogry category 3x: modified fit       2x: end fit
% category category 1x: start fit           x: normal procedure
% cateogry category x1: mean value         x4: simple model
% category category x6: advanced model     x8: full model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  defining variables
%% Constants
% maximum number of coefficient required
% 2 for tri-modal; 2 for bi-modal; 2 for seasonality; 2 for linear;
min_num_c = 4;
mid_num_c = 6;
max_num_c = 8;
% number of clear observation / number of coefficients
n_times = 3;
% initialize NUM of Functional Curves
num_fc = 0;
% number of days per year
num_yrs = 365.25;
% number of bytes: int16
num_byte = 2;
% Band for multitemporal cloud/snow detection (Green)
num_B1 = 1;
% Band for multitemporal shadow/snow shadow detection (SWIR)
num_B2 = 6;
% Threshold for cloud, shadow, and snow detection.
T_const = 4.42;
% minimum year for model intialization
mini_yrs = 1;
% no change detection for permanent snow pixels
t_sn = 0.75;
% threshold (degree) of mean included angle 
nsign = 45;
% Fmask fails threshold
% t_clr = 0.25;
% agriculture activity (9)
% get num of total folders start with "L"
imf = dir('L*'); % folder names
% filter for Landsat folders
imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
imf = [imf{:}];
imf = vertcat(imf{:});
% sort according to yeardoy
yeardoy = str2num(imf(:, 10:16)); 
[~, sort_order] = sort(yeardoy);
imf = imf(sort_order, :);
% number of folders start with "L"
num_t = size(imf,1);
% consecutive number
def_conse = conse;

% Tmasking of noise
Tmax_cg = chi2inv(Tmax_cg,length(B_detect));
% adjust threshold based on chi-squared distribution
def_pT_cg = T_cg;
def_T_cg = chi2inv(def_pT_cg,length(B_detect));

% initialize the struct data of RECording of ChanGe (rec_cg)
rec_cg = struct('t_start',[],'t_end',[],'t_break',[],'coefs',[],'rmse',[],...
    'pos',[],'change_prob',[],'num_obs',[],'category',[],'magnitude',[]);

% % mask for study area (1 fit, 0 no fit)
% fit_mask = enviread('GZ_Mask');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Get ready for Xs & Ys
%% Read in Xs & Ysq
% transforming to serial date number (0000 year)
sdate = zeros(num_t,1); % Xs
line_t = zeros(num_t,nbands*ncols); %Ys

for i = 1:num_t
    im_dir = dir(imf(i, :));
    im = '';
    for f = 1:size(im_dir, 1)
        % use regular expression to match:
        %   'L(\w*)'    Any word begining with L that has any following chars
        %   stk_n       includes stack name somewhere after L
        %   '$'           ends with the stack name (e.g., no .hdr, .aux.xml)
        if regexp(im_dir(f).name, ['L(\w*)', 'stack']) == 1
            im = ['C:/Master_Thesis/COLD/Data','/',imf(i, :), '/', im_dir(f).name];
            break
        end
    end
    % Check to make sure we found something
    if strcmp(im, '')
        error('Could not find stack image for directory %s\n', imf(i));
    end
    % Find date for folder imf(i)
    yr = str2num(imf(i, 10:13)); %#ok<*ST2NM>
    doy = str2num(imf(i, 14:16));
    sdate(i) = datenum(yr, 1, 0) + doy;
    dummy_name = im;
    fid_t = fopen(dummy_name,'r'); % get file ids
    fseek(fid_t,num_byte*(nrows-1)*ncols*nbands,'bof');
    line_t(i,:) = fread(fid_t,nbands*ncols,'uint16','ieee-le'); % get Ys
end
fclose('all'); % close all files

for i_ids = 1:ncols
    % get default conse & T_cg
    conse = def_conse;
    T_cg = def_T_cg;

    %     % Only run CCDC for places where more than 50% of images has data
    %     idexist = line_m < 255;
    %     overlap_pct = sum(idexist)/num_t;
    %     if overlap_pct < 0.5
    %         continue;
    %     end
    

    
    % clear pixel should have reflectance between 0 and 1
    % brightness temperature should between -93.2 to 70.7 celsius degree
    idrange = line_t(:,nbands*(i_ids-1)+1)>0&line_t(:,nbands*(i_ids-1)+1)<10000&...
        line_t(:,nbands*(i_ids-1)+2)>0&line_t(:,nbands*(i_ids-1)+2)<10000&...
        line_t(:,nbands*(i_ids-1)+3)>0&line_t(:,nbands*(i_ids-1)+3)<10000&...
        line_t(:,nbands*(i_ids-1)+4)>0&line_t(:,nbands*(i_ids-1)+4)<10000&...
        line_t(:,nbands*(i_ids-1)+5)>0&line_t(:,nbands*(i_ids-1)+5)<10000&...
        line_t(:,nbands*(i_ids-1)+6)>0&line_t(:,nbands*(i_ids-1)+6)<10000;
    
  
            idgood = idrange;
            
            % Xs & Ys for computation
            clrx = sdate(idgood);
            % bands 1-5,7,6
            clry = line_t(idgood,(nbands*(i_ids-1)+1):(nbands*(i_ids-1)+nbands));
            clry = double(clry);
            % find repeated ids
            [clrx,uniq_id,~] = unique(clrx);
            % mean of repeated values
            tmp_y = zeros(length(clrx),nbands);
            % get the mean values
            for i = 1:nbands-1
                tmp_y(:,i) = clry(uniq_id,i);
            end
            clry = tmp_y;
            
            idclr = clry(:,num_B1) < median(clry(:,num_B1)) + 400;
            n_clr = sum(idclr);
            
            if n_clr < n_times*min_num_c % not enough clear pixels
                continue
            else
                % Xs & Ys for computation
                clrx = clrx(idclr);
                clry = clry(idclr,:);
                
                % the first observation for TSFit
                i_start = 1;
                % identified and move on for the next curve
                num_fc = num_fc + 1; % NUM of Fitted Curves (num_fc)
                
                % defining computed variables
                fit_cft = zeros(max_num_c,nbands);
                % rmse for each band
                rmse = zeros(nbands,1);
                % Fmask fail qa = 40
                qa = 40;
                
                for i_B = 1:nbands
                    % fit basic model for all within range snow pixels
                    [fit_cft(:,i_B),rmse(i_B)] = autoTSFit(clrx,clry(:,i_B),min_num_c);
                end
                
                % record time of curve start
                rec_cg(num_fc).t_start = clrx(i_start);
                % record time of curve end
                rec_cg(num_fc).t_end = clrx(end);
                % record break time
                rec_cg(num_fc).t_break = 0;
                % record postion of the pixel
                rec_cg(num_fc).pos = (nrows-1)*ncols+i_ids;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = length(clrx);
                % record fit category
                rec_cg(num_fc).category = qa + min_num_c;
                % record change magnitude
                rec_cg(num_fc).magnitude = zeros(1,nbands);


            end % end of if sum(idgood) statement    
end % end of for i_ids loop
% save([dir_l,'/',n_rst,'/','record_change',num2str(nrows)],'rec_cg');
% when saving matlab file, the computer is shut off. The matlab file may be broken.
file_name_mat = ['record_change',num2str(nrows),'.mat'];
file_name_mat_part = [file_name_mat, '.part'];
% save it first as name with part
save([dir_l,'/',n_rst,'/', file_name_mat_part],'rec_cg');
% rename it
movefile([dir_l,'/',n_rst,'/',file_name_mat_part], [dir_l,'/',n_rst,'/',file_name_mat]);

end % end of function

% function to caculate included angle between ajacent pair of change vectors
function y = angl(v_dif)

[row,~] = size(v_dif);
y = zeros(row-1,1);

if row > 1
    for i = 1:row-1
        a = v_dif(i,:);
        b = v_dif(i+1,:);
        % y measures the opposite of cos(angle)
        y(i) = acos(a*b'/(norm(a)*norm(b)));
    end
else
    y = 0;
end

% convert angle from radiance to degree
y = y*180/pi;

end