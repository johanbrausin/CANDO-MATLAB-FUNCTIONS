function seqInfo = csv2seq(filename)

seqInfo = [];

fid = fopen(filename);
fgetl(fid);                 % Read the header

tline = fgetl(fid);
while ischar(tline)
    C = strsplit(tline,',');
    assert(numel(C)==5);
    
    currSeqInfo.start = str2num(strrep(strrep(C{1},'[',' '),']',' '));
    currSeqInfo.end = str2num(strrep(strrep(C{2},'[',' '),']',' '));
    currSeqInfo.sequence = C{3};
    currSeqInfo.len = str2double(C{4});
    assert(numel(currSeqInfo.start)==2 && isempty(find(~isfinite(currSeqInfo.start),1)));
    assert(numel(currSeqInfo.end)==2 && isempty(find(~isfinite(currSeqInfo.end),1)));
    assert(numel(currSeqInfo.sequence)==currSeqInfo.len);
    
    seqInfo = cat(1,seqInfo,currSeqInfo);
    
    tline = fgetl(fid);
end

fclose(fid);

end