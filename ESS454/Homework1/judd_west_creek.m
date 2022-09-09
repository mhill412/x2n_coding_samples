%% Matty Hill
%% ESS454
%% Vashon-Maury Island Water Budget
%% 01/12/22

%% Finding the water year indexes for each water year

 id2011=findwateryearindex(JuddWestPrecip,2011); 
 id2012=findwateryearindex(JuddWestPrecip,2012);
 id2013=findwateryearindex(JuddWestPrecip,2013);
 id2014=findwateryearindex(JuddWestPrecip,2014);
 id2015=findwateryearindex(JuddWestPrecip,2015);
 id2016=findwateryearindex(JuddWestPrecip,2016);
 id2017=findwateryearindex(JuddWestPrecip,2017);
 id2018=findwateryearindex(JuddWestPrecip,2018);
 id2019=findwateryearindex(JuddWestPrecip,2019);
 id2020=findwateryearindex(JuddWestPrecip,2020);
 id2021=findwateryearindex(JuddWestPrecip,2021);
 
 %% Analyzing the discharge and precipitation for each water year
 
 years = [2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021];
 for c = 1:length(years)
    analysis = AnalyzeJudd(JuddCreekDischarge,JuddWestPrecip,years(c))
    saveas(analysis, ['dis_precip_analysis' num2str(years(c))])
 end