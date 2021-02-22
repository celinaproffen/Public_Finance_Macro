% ******************************************************************************
% readfile.m
% 
% Reads files in ascii format. First line is deleted. 
% n: gives dimension of input matrix from input file
% ******************************************************************************

function z = readfile(x, y, n)

fname = [x,y];
fid = fopen(fname);

fline = fgets(fid);     % first line is not important
z = fscanf(fid, '%g',[n, inf]);     % reads ascii, tab del., rows as columns

z = z(2 : n, :);
z = z';

fclose(fid);


