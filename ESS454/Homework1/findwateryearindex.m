function id=findwateryearindex(Data,yr)
% function searches the third column (local time) to find all in a
% specified water year (Oct 1 to sept. 30) and returns the indexes for that
% year.  JMB 2021

varname=Data.Properties.VariableNames{3};
id1=find(year(Data.(varname)(:))==yr-1 & month(Data.(varname)(:))==10 & day(Data.(varname)(:))==1  );
id2=find(year(Data.(varname)(:))==yr   & month(Data.(varname)(:))==9 &  day(Data.(varname)(:))==30  );
id=id1:id2;