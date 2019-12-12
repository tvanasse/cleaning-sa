% written by BAR 1/13/15
% removed sample rate 4/13/16
% modified to use take one input 6/20/17
function samples = readalignmentraw(rawpathname,rawfilenames,start_sample,end_sample)
if nargin == 1
    [pathname,rawfile,ext] = fileparts(rawpathname);
    rawpathname = pathname;
    rawfilenames = {[rawfile,ext]};
end
num_raws = length(rawfilenames);

for i = 1:length(rawfilenames)
    fileID      = fopen([rawpathname,filesep,rawfilenames{i}]);
    
    if nargin < 3
        channel_samples = fread(fileID,'*single');
        num_samples = length(channel_samples); 
        
    elseif nargin == 4
        ftell(fileID);
        num_samples = end_sample - start_sample + 1;
        
        fseek(fileID, start_sample*4, 'bof');
        
        channel_samples = fread(fileID,[1,num_samples],'*single')'; %multiply sample to move bytes
    else
        error('wrong number of inputs');
    end
    
    if i == 1
        samples       = zeros(num_raws,num_samples,'single');
        samples(1,:)  = channel_samples;
    else
        samples(i,:)  = channel_samples;
    end
    
    fclose(fileID);
end