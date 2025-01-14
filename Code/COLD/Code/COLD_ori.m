clc
clear
close all


cd 'C:/Master_Thesis/'
root_dir = 'C:/Master_Thesis/';

data_dir = sprintf('%COLD/Data/rstudio-export/', root_dir);
addpath(data_dir)
code_dir = sprintf('%sGitHub_code/COLD',root_dir);
addpath(code_dir)

Read the data
opts = detectImportOptions("COLD/Data/rstudio-export/Lansat_8_SRB2.csv")
opts = setvartype(opts, 'double')
LandsatSR1 = readtable("COLD/Data/rstudio-export/Lansat_8_SRB2.csv", opts)
LandsatSR1(:,1:2) = []
opts = detectImportOptions("COLD/Data/rstudio-export/Lansat_8_SRB3.csv")
opts = setvartype(opts, 'double')
LandsatSR2 = readtable("COLD/Data/rstudio-export/Lansat_8_SRB3.csv", opts)
LandsatSR2(:,1:2) = []
opts = detectImportOptions("COLD/Data/rstudio-export/Lansat_8_SRB4.csv")
opts = setvartype(opts, 'double')
LandsatSR3 = readtable("COLD/Data/rstudio-export/Lansat_8_SRB4.csv", opts)
LandsatSR3(:,1:2) = []
opts = detectImportOptions("COLD/Data/rstudio-export/Lansat_8_SRB5.csv")
opts = setvartype(opts, 'double')
LandsatSR4 = readtable("COLD/Data/rstudio-export/Lansat_8_SRB5.csv", opts)
LandsatSR4(:,1:2) = []
opts = detectImportOptions("COLD/Data/rstudio-export/Lansat_8_SRB6.csv")
opts = setvartype(opts, 'double')
LandsatSR5 = readtable("COLD/Data/rstudio-export/Lansat_8_SRB6.csv", opts)
LandsatSR5(:,1:2) = []
opts = detectImportOptions("COLD/Data/rstudio-export/Lansat_8_SRB7.csv")
opts = setvartype(opts, 'double')
LandsatSR6 = readtable("COLD/Data/rstudio-export/Lansat_8_SRB7.csv", opts)
LandsatSR6(:,1:2) = []


root_dir = 'C:/Master_Thesis/';

data_dir = sprintf('%COLD/Data/', root_dir);
addpath(data_dir)

%save("Landsatsixband_stack.mat", "LandsatSR1", "LandsatSR2", "LandsatSR3", "LandsatSR4", "LandsatSR5", "LandsatSR6")

LandsatSR1 = table2array(LandsatSR1)
LandsatSR2 = table2array(LandsatSR2)
LandsatSR3 = table2array(LandsatSR3)
LandsatSR4 = table2array(LandsatSR4)
LandsatSR5 = table2array(LandsatSR5)
LandsatSR6 = table2array(LandsatSR6)

stack bip
% Preallocate a 3D matrix to store the stacked data
stack = zeros(11773, 186, 6, 'double');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\Landsat8sixband\Landsatsixband_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = LandsatSR1;
stack(:,:,2) = LandsatSR2;
stack(:,:,3) = LandsatSR3;
stack(:,:,4) = LandsatSR4;
stack(:,:,5) = LandsatSR5;
stack(:,:,6) = LandsatSR6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'double');

% Close the file after writing
fclose(fileID);


% Define dimensions and data type
rows = 11818;
cols = 186;
bands = 6;
datatype = 'double'; % Same data type used for writing

% Open the BIP file for reading
fileID = fopen('C:\Master_Thesis\COLD\Data\Landsat8sixband\Landsat_stack\stacked_matrices.bip', 'r');

% Read the data as a vector
data = fread(fileID, rows * cols * bands, datatype);

% Close the file
fclose(fileID);

% Reshape the data into a 3D matrix (rows x cols x bands)
stack = reshape(data, [rows, cols, bands]);

LandsatSR1 = "LandsatSR1.mat"
save(LandsatSR1, "LandsatSR1")
LandsatSR2 = "LandsatSR2.mat"
save(LandsatSR2, "LandsatSR2")
LandsatSR3 = "LandsatSR3.mat"
save(LandsatSR3, "LandsatSR3")
LandsatSR4 = "LandsatSR4.mat"
save(LandsatSR4, "LandsatSR4")
LandsatSR5 = "LandsatSR5.mat"
save(LandsatSR5, "LandsatSR5")
LandsatSR6 = "LandsatSR6.mat"
save(LandsatSR6, "LandsatSR6")

ymd = readtable('C:\Master_Thesis\DRMAT\Data\ymd_ndvi_cal.csv');
ymd(:,1) = []
ymd = table2array(ymd)

sdate = zeros(186,1);
for i = 1:size(ymd,1)
    yr = ymd(i,1);
    mm = ymd(i,2);
    dd = ymd(i,3);
    doy = datenummx(yr,mm,dd)-datenummx(yr,1,0);
    sdate(i) = datenum(yr, 1, 0) + doy;
end

sdate = reshape(sdate, 1, 186)

%sdate_rep = repmat(sdate, 1, 6)
fileName = 'C:\Master_Thesis\COLD\Data\sdate.txt';
fid = fopen(fileName, 'wt');
[rows, columns] = size(sdate)
for row = 1:rows
    for col = 1: columns
        fprintf(fid, '%d ', sdate(row,col))
    end
    fprintf(fid, '\n')
end
fclose(fid)

fileName = 'C:\Master_Thesis\COLD\Data\sdate.txt'; % Path to your file
fid = fopen(fileName, 'rt'); % Open the file in read mode (text format)

if fid == -1
    error('Failed to open file. Check if the path is correct.');
end

% Read the data into a matrix
sdate_rep_read = fscanf(fid, '%d '); 

fclose(fid); % Close the file

% Reshape the data to match its original structure
rows = 1; % Original row size
columns = 186 ; % Number of columns based on the original replication
sdate_rep_read = reshape(sdate_rep_read, rows, columns);

disp(sdate_rep_read); % Display the read matrix

Run main
main_COLD(1,1)

test resut
root_dir = 'C:/Master_Thesis/COLD/Data/'
data_dir = sprintf('%TSFitMap/', root_dir);
addpath(data_dir)

x = load('C:/Master_Thesis/COLD/Data/TSFitMap/record_change1001.mat')
y = load('C:/Master_Thesis/COLD/Data/TSFitMap/record_change1002.mat')
x_table = struct2table(x.rec_cg, AsArray= true)


combine all the result of SR cal
root_path = 'C:/Master_Thesis/COLD/Data/TSFitMap/';
file_list = dir(fullfile(root_path,'*.mat'));
tic
results = cell(numel(file_list), 1); 
parfor k = 1:numel(file_list)
    A  =  fullfile(root_path, file_list(k).name);
    data = load(A);
        % Convert the structure to a table
    A_table = struct2table(data.rec_cg, 'AsArray', true);
    results{k} = A_table;
end
toc
% Combine all tables into a single table
result = vertcat(results{:});
SR_cal_COLD_output = sortrows(result,"ID");


t_break = SR_cal_COLD_output(:,"t_break")
t_break_date = cell(numel(t_break), 1)
for k = 10:numel(t_break)

t_break_date{k} = datetime(t_break.t_break(k,:),'convertfrom','juliandate','Format','yyy-MM-dd')
end









% Specify the folder containing .mat files
folder = 'C:/Master_Thesis/COLD/Data/TSFitMap/';  % Update with your folder path
matFiles = dir(fullfile(folder, '*.mat'));  % Get all .mat files in the folder
% for file_nr = 1:length(matFiles)
% Initialize a cell array to store file names with meaningful variables
filesWithNonEmptyContent = {};

% Loop through all .mat files
for i = 1:length(matFiles)
    % Full path to the current .mat file
    matFilePath = fullfile(folder, matFiles(i).name);
    
    % Load all variables from the .mat file
    loadedData = load(matFilePath);
    if ~isempty(extractfield(loadedData.rec_cg, 't_break'))| sum(extractfield(loadedData.rec_cg, 't_break')) ~= 0
        filesWithNonEmptyContent = {i};
    end
    
%     % Check each variable for meaningful content
%     variableNames = fieldnames(loadedData);  % Get variable names
%     hasMeaningfulVariable = false;  % Flag to track meaningful variables
%     for j = 1:length(variableNames)
%         variable = loadedData.(variableNames{j});  % Access the variable
% 
% 
% 
%         % Initialize an empty cell array to store non-empty variable names
% nonEmptyVariableNames = {};
% 
% % Get all field names of the struct
% fieldNames = fieldnames(variable);
% % Loop through each field to check for non-emptiness
% for j = 1:length(fieldNames)
%     fieldName = fieldNames{j};
%     if ~isempty(variable.(fieldName))
%         % Add non-empty variable name to the list
%         nonEmptyVariableNames{end+1} = fieldName; %#ok<AGROW>
%         filesWithNonEmptyContent{end+1} = matFiles(i).name; %#ok<AGROW>
%     end
% end
% 
% 
end
    

% end

% Display the results
if isempty(filesWithNonEmptyContent)
    disp('No files contain variables with meaningful values.');
else
    disp('Files with variables that contain meaningful values:');
    disp(filesWithNonEmptyContent);
end

fileID = fopen('nine.bin','w');
fwrite(fileID,[1:1000],'int16');
fclose(fileID);
fileID = fopen('nine.bin');
A = fread(fileID,70,'uint16','ieee-le')

% Define dimensions
nRows = 11811;  % Number of rows
nCols = 186;    % Number of columns
nLayers = 6;    % Number of layers

% Create data
data = reshape(1:(nRows * nCols * nLayers), [nRows, nCols, nLayers]);

% Open a binary file for writing
fileID = fopen('six_layers.bin', 'w');

% Write the data to the binary file layer by layer
for layer = 1:nLayers
    fwrite(fileID, data(:, :, layer), 'int16'); % Write each layer as 'int16'
end

% Open the binary file for reading
fileID = fopen('six_layers.bin', 'r');

% Read the data
readData = fread(fileID, [nRows * nCols * nLayers], 'int16');

% Reshape the data to the original dimensions
readData = reshape(readData, [nRows, nCols, nLayers]);

% Close the file
fclose(fileID);

% Display the values from one layer (e.g., the first layer)
disp('Values in the first layer:');
disp(readData(:, :, 1));
disp(readData(:, :, 2));
disp(readData(:, :, 3));

% Example: Display a specific row from the first layer
rowIndex = 1; % Choose the row you want to view
disp(['Values in row ', num2str(rowIndex), ' of the first layer:']);
disp(readData(rowIndex, :, 1));

Read example of ARD tif file
img1 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241007_20241019_02\LC08_CU_017007_20241007_20241019_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241007_20241019_02\LC08_CU_017007_20241007_20241019_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241007_20241019_02\LC08_CU_017007_20241007_20241019_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241007_20241019_02\LC08_CU_017007_20241007_20241019_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241007_20241019_02\LC08_CU_017007_20241007_20241019_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241007_20241019_02\LC08_CU_017007_20241007_20241019_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)

stack
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC08_CU__20241007_20241019_02\Landsat_stack\LC08_CU__20241007_20241019_02.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);
% Define dimensions and data type
rows = 5000;
cols = 5000;
bands = 6;
datatype = 'uint16'; % Same data type used for writing

% Open the BIP file for reading
fileID = fopen('C:\Master_Thesis\COLD\Data\LC08_CU__20241007_20241019_02\Landsat_stack\LC08_CU__20241007_20241019_02.bip', 'r');

% Read the data as a vector
data = fread(fileID, rows * cols * bands, datatype);

% Close the file
fclose(fileID);

% Reshape the data into a 3D matrix (rows x cols x bands)
stack = reshape(data, [rows, cols, bands]);

img1 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241014_20241025_02\LC08_CU_017007_20241014_20241025_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241014_20241025_02\LC08_CU_017007_20241014_20241025_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241014_20241025_02\LC08_CU_017007_20241014_20241025_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241014_20241025_02\LC08_CU_017007_20241014_20241025_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241014_20241025_02\LC08_CU_017007_20241014_20241025_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241014_20241025_02\LC08_CU_017007_20241014_20241025_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC08_CU__20241014_20241025_02\Landsat_stack\LC08_CU__20241014_20241025_02.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);


img1 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241101_20241108_02\LC08_CU_017007_20241101_20241108_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241101_20241108_02\LC08_CU_017007_20241101_20241108_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241101_20241108_02\LC08_CU_017007_20241101_20241108_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241101_20241108_02\LC08_CU_017007_20241101_20241108_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241101_20241108_02\LC08_CU_017007_20241101_20241108_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241101_20241108_02\LC08_CU_017007_20241101_20241108_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC08_CU__20241101_20241108_02\Landsat_stack\LC08_CU__20241101_20241108_02.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);




img1 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241108_20241117_02\LC08_CU_017007_20241108_20241117_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241108_20241117_02\LC08_CU_017007_20241108_20241117_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241108_20241117_02\LC08_CU_017007_20241108_20241117_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241108_20241117_02\LC08_CU_017007_20241108_20241117_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241108_20241117_02\LC08_CU_017007_20241108_20241117_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241108_20241117_02\LC08_CU_017007_20241108_20241117_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC08_CU__20241108_20241117_02\Landsat_stack\LC08_CU__20241108_20241117_02.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);



img1 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241115_20241122_02\LC08_CU_017007_20241115_20241122_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241115_20241122_02\LC08_CU_017007_20241115_20241122_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241115_20241122_02\LC08_CU_017007_20241115_20241122_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241115_20241122_02\LC08_CU_017007_20241115_20241122_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241115_20241122_02\LC08_CU_017007_20241115_20241122_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241115_20241122_02\LC08_CU_017007_20241115_20241122_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC08_CU__20241115_20241122_02\Landsat_stack\LC08_CU__20241115_20241122_02.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);


img1 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241117_20241130_02\LC08_CU_017007_20241117_20241130_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241117_20241130_02\LC08_CU_017007_20241117_20241130_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241117_20241130_02\LC08_CU_017007_20241117_20241130_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241117_20241130_02\LC08_CU_017007_20241117_20241130_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241117_20241130_02\LC08_CU_017007_20241117_20241130_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241117_20241130_02\LC08_CU_017007_20241117_20241130_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC08_CU__20241117_20241130_02\Landsat_stack\LC08_CU__20241117_20241130_02.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);


img1 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241124_20241130_02\LC08_CU_017007_20241124_20241130_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241124_20241130_02\LC08_CU_017007_20241124_20241130_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241124_20241130_02\LC08_CU_017007_20241124_20241130_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241124_20241130_02\LC08_CU_017007_20241124_20241130_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241124_20241130_02\LC08_CU_017007_20241124_20241130_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC08_CU__20241124_20241130_02\LC08_CU_017007_20241124_20241130_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC08_CU__20241124_20241130_02\Landsat_stack\LC08_CU__20241124_20241130_02.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);

new ard
img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC8_CU1__20230106_20230114_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);

img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230207_20230221_02\LC08_CU_017007_20230207_20230221_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230207_20230221_02\LC08_CU_017007_20230207_20230221_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230207_20230221_02\LC08_CU_017007_20230207_20230221_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230207_20230221_02\LC08_CU_017007_20230207_20230221_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230207_20230221_02\LC08_CU_017007_20230207_20230221_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230207_20230221_02\LC08_CU_017007_20230207_20230221_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20230207_20230221_02\LC8_CU1__20230207_20230221_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);

img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230327_20230410_02\LC08_CU_017007_20230327_20230410_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230327_20230410_02\LC08_CU_017007_20230327_20230410_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230327_20230410_02\LC08_CU_017007_20230327_20230410_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230327_20230410_02\LC08_CU_017007_20230327_20230410_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230327_20230410_02\LC08_CU_017007_20230327_20230410_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230327_20230410_02\LC08_CU_017007_20230327_20230410_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20230327_20230410_02\LC8_CU1__20230327_20230410_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);




img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231005_20231014_02\LC08_CU_017007_20231005_20231014_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231005_20231014_02\LC08_CU_017007_20231005_20231014_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231005_20231014_02\LC08_CU_017007_20231005_20231014_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231005_20231014_02\LC08_CU_017007_20231005_20231014_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231005_20231014_02\LC08_CU_017007_20231005_20231014_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231005_20231014_02\LC08_CU_017007_20231005_20231014_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20231005_20231014_02\LC8_CU1__20231005_20231014_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);



img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230428_20230512_02\LC08_CU_017007_20230428_20230512_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230428_20230512_02\LC08_CU_017007_20230428_20230512_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230428_20230512_02\LC08_CU_017007_20230428_20230512_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230428_20230512_02\LC08_CU_017007_20230428_20230512_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230428_20230512_02\LC08_CU_017007_20230428_20230512_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230428_20230512_02\LC08_CU_017007_20230428_20230512_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20230428_20230512_02\LC8_CU1__20230428_20230512_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);


img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230523_20230606_02\LC08_CU_017007_20230523_20230606_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230523_20230606_02\LC08_CU_017007_20230523_20230606_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230523_20230606_02\LC08_CU_017007_20230523_20230606_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230523_20230606_02\LC08_CU_017007_20230523_20230606_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230523_20230606_02\LC08_CU_017007_20230523_20230606_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230523_20230606_02\LC08_CU_017007_20230523_20230606_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20230523_20230606_02\LC8_CU1__20230523_20230606_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);


img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230615_20230626_02\LC08_CU_017007_20230615_20230626_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230615_20230626_02\LC08_CU_017007_20230615_20230626_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230615_20230626_02\LC08_CU_017007_20230615_20230626_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230615_20230626_02\LC08_CU_017007_20230615_20230626_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230615_20230626_02\LC08_CU_017007_20230615_20230626_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230615_20230626_02\LC08_CU_017007_20230615_20230626_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20230615_20230626_02\LC8_CU1__20230615_20230626_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);


img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230717_20230728_02\LC08_CU_017007_20230717_20230728_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230717_20230728_02\LC08_CU_017007_20230717_20230728_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230717_20230728_02\LC08_CU_017007_20230717_20230728_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230717_20230728_02\LC08_CU_017007_20230717_20230728_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230717_20230728_02\LC08_CU_017007_20230717_20230728_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230717_20230728_02\LC08_CU_017007_20230717_20230728_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20230717_20230728_02\LC8_CU1__20230717_20230728_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);



img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230903_20230915_02\LC08_CU_017007_20230903_20230915_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230903_20230915_02\LC08_CU_017007_20230903_20230915_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230903_20230915_02\LC08_CU_017007_20230903_20230915_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230903_20230915_02\LC08_CU_017007_20230903_20230915_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230903_20230915_02\LC08_CU_017007_20230903_20230915_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230903_20230915_02\LC08_CU_017007_20230903_20230915_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20230903_20230915_02\LC8_CU1__20230903_20230915_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);



img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231106_20231115_02\LC08_CU_017007_20231106_20231115_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231106_20231115_02\LC08_CU_017007_20231106_20231115_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231106_20231115_02\LC08_CU_017007_20231106_20231115_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231106_20231115_02\LC08_CU_017007_20231106_20231115_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231106_20231115_02\LC08_CU_017007_20231106_20231115_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231106_20231115_02\LC08_CU_017007_20231106_20231115_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20231106_20231115_02\LC8_CU1__20231106_20231115_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);



img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231208_20231218_02\LC08_CU_017007_20231208_20231218_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231208_20231218_02\LC08_CU_017007_20231208_20231218_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231208_20231218_02\LC08_CU_017007_20231208_20231218_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231208_20231218_02\LC08_CU_017007_20231208_20231218_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231208_20231218_02\LC08_CU_017007_20231208_20231218_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20231208_20231218_02\LC08_CU_017007_20231208_20231218_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20231208_20231218_02\LC8_CU1__20231208_20231218_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);



img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240109_20240126_02\LC08_CU_017007_20240109_20240126_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240109_20240126_02\LC08_CU_017007_20240109_20240126_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240109_20240126_02\LC08_CU_017007_20240109_20240126_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240109_20240126_02\LC08_CU_017007_20240109_20240126_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240109_20240126_02\LC08_CU_017007_20240109_20240126_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240109_20240126_02\LC08_CU_017007_20240109_20240126_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20240109_20240126_02\LC8_CU1__20240109_20240126_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);




img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240329_20240413_02\LC08_CU_017007_20240329_20240413_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240329_20240413_02\LC08_CU_017007_20240329_20240413_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240329_20240413_02\LC08_CU_017007_20240329_20240413_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240329_20240413_02\LC08_CU_017007_20240329_20240413_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240329_20240413_02\LC08_CU_017007_20240329_20240413_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240329_20240413_02\LC08_CU_017007_20240329_20240413_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20240329_20240413_02\LC8_CU1__20240329_20240413_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);



img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240516_20240524_02\LC08_CU_017007_20240516_20240524_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240516_20240524_02\LC08_CU_017007_20240516_20240524_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240516_20240524_02\LC08_CU_017007_20240516_20240524_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240516_20240524_02\LC08_CU_017007_20240516_20240524_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240516_20240524_02\LC08_CU_017007_20240516_20240524_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240516_20240524_02\LC08_CU_017007_20240516_20240524_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20240516_20240524_02\LC8_CU1__20240516_20240524_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);



img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240617_20240709_02\LC08_CU_017007_20240617_20240709_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240617_20240709_02\LC08_CU_017007_20240617_20240709_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240617_20240709_02\LC08_CU_017007_20240617_20240709_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240617_20240709_02\LC08_CU_017007_20240617_20240709_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240617_20240709_02\LC08_CU_017007_20240617_20240709_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240617_20240709_02\LC08_CU_017007_20240617_20240709_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20240617_20240709_02\LC8_CU1__20240617_20240709_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);



img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240820_20240903_02\LC08_CU_017007_20240820_20240903_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240820_20240903_02\LC08_CU_017007_20240820_20240903_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240820_20240903_02\LC08_CU_017007_20240820_20240903_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240820_20240903_02\LC08_CU_017007_20240820_20240903_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240820_20240903_02\LC08_CU_017007_20240820_20240903_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240820_20240903_02\LC08_CU_017007_20240820_20240903_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20240820_20240903_02\LC8_CU1__20240820_20240903_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);




img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240921_20241001_02\LC08_CU_017007_20240921_20241001_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240921_20241001_02\LC08_CU_017007_20240921_20241001_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240921_20241001_02\LC08_CU_017007_20240921_20241001_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240921_20241001_02\LC08_CU_017007_20240921_20241001_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240921_20241001_02\LC08_CU_017007_20240921_20241001_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240921_20241001_02\LC08_CU_017007_20240921_20241001_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20240921_20241001_02\LC8_CU1__20240921_20241001_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);




img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240930_20241009_02\LC08_CU_017007_20240930_20241009_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240930_20241009_02\LC08_CU_017007_20240930_20241009_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240930_20241009_02\LC08_CU_017007_20240930_20241009_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240930_20241009_02\LC08_CU_017007_20240930_20241009_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240930_20241009_02\LC08_CU_017007_20240930_20241009_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20240930_20241009_02\LC08_CU_017007_20240930_20241009_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20240930_20241009_02\LC8_CU1__20240930_20241009_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');

% Close the file after writing
fclose(fileID);

img1 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B2.TIF', 'r')
SR_B1 = read(img1)
img2 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B3.TIF', 'r')
SR_B2 = read(img2)
img3 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B4.TIF', 'r')
SR_B3 = read(img3)
img4 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B5.TIF', 'r')
SR_B4 = read(img4)
img5 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B6.TIF', 'r')
SR_B5 = read(img5)
img6 = Tiff('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC08_CU_017007_20230106_20230114_02_SR_B7.TIF', 'r')
SR_B6 = read(img6)
% Preallocate a 3D matrix to store the stacked data
stack = zeros(5000, 5000, 6, 'uint16');

% Open a file to write in BIP format
fileID = fopen('C:\Master_Thesis\COLD\Data\LC8_CU1__20230106_20230114_02\LC8_CU1__20230106_20230114_02_stack.bip', 'w');

% Assign Landsat bands to the respective layers in the stack
stack(:,:,1) = SR_B1;
stack(:,:,2) = SR_B2;
stack(:,:,3) = SR_B3;
stack(:,:,4) = SR_B4;
stack(:,:,5) = SR_B5;
stack(:,:,6) = SR_B6;

% Write data in BIP format (row-major order)
fwrite(fileID, stack, 'uint16');
