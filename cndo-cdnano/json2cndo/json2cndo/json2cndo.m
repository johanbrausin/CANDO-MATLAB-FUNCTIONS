function [] = json2cndo(jsonPATH, csvPATH, cndoPATH, latticeType)

% Convert the .json and .csv files to dnaInfo
dnaInfo = json2dnaInfo(jsonPATH, csvPATH, latticeType);
% Convert dnaInfo to the .cndo file
dnaInfo2cndo(dnaInfo, cndoPATH);

end