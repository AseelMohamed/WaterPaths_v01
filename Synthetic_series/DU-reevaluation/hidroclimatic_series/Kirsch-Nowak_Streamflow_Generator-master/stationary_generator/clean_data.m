% This loads the mat files and combines the Q (cfs) time series into one
% txt file. Removes leap years by averaging with the prior date.
% MUDAR O CAMINHO DO DIRETORIO E O NOME DOS ARQUIVOS CONFORME O CASO

clc; clear all;
datadir = 'E:/UNB/Projeto/PRGroup/series_sinteticas/DU_reeval/streamflow_generator/Kirsch-Nowak_Streamflow_Generator-master/data/';
files = { 'qDescoberto_1971-2020.csv', ...
    'qSantaMaria_1971-2020.csv', 'qParanoa_1971-2020.csv', ...
    'qTorto_Bananal_1971-2020.csv',...
    'qCorumba_1971-2020.csv', 'evapDescoberto_1971-2020.csv', 'evapSantaMaria_1971-2020.csv','evapParanoa_1971-2020.csv',...
    'evapCorumba_1971-2020.csv'};
hist_data = [];

% load the historical data
for i=1:length(files)
    M{i} = load([datadir  files{i}]);
end

% find indices of leap years
% ALTERAR ABAIXO CONFORME O CASO
leaps = 60:365*3+366:365*(2020-1971+1)+ceil(2020-1971)/4;
all = 1:1:365*(2020-1971+1)+ceil(2020-1971)/4+1;
non_leaps = setdiff(all,leaps);

Qfinal = zeros(length(non_leaps),length(files));

for i=1:length(files)
    Q = M{i};
    Qfinal(:,i) = Q(non_leaps);
end

dlmwrite('E:/UNB/Projeto/PRGroup/series_sinteticas/DU_reeval/streamflow_generator/Kirsch-Nowak_Streamflow_Generator-master/data/Qdaily.txt', Qfinal, ' ');

% reshape into nyears x 365 and nyears x 12 for daily and monthly 
% statistical validation figures
Qfinal_monthly = convert_data_to_monthly(Qfinal);
% divide evaporation by 86400 (s/day) to get total monthly evap in in/month
%CODIGO ABAIXO COLOCADO COMO COMENTÁRIO PORQUE A APLICAÇÃO NÃO TEM
%EVAPORAÇÃO A PRINCIPIO
Qfinal_monthly{6} = Qfinal_monthly{6}/86400;
Qfinal_monthly{7} = Qfinal_monthly{7}/86400;
Qfinal_monthly{8} = Qfinal_monthly{8}/86400;
Qfinal_monthly{9} = Qfinal_monthly{9}/86400;


% create directories to write files to
mkdir('E:/UNB/Projeto/PRGroup/series_sinteticas/DU_reeval/streamflow_generator/Kirsch-Nowak_Streamflow_Generator-master/validation/historical');
for i=1:length(files)
   q_nx365 = reshape(Qfinal(:,i),365, [])';
   dlmwrite(['E:/UNB/Projeto/PRGroup/series_sinteticas/DU_reeval/streamflow_generator/Kirsch-Nowak_Streamflow_Generator-master/validation/historical/' files{i}(1:(length(files{i})-14)) '-daily.csv'], q_nx365);
   dlmwrite(['E:/UNB/Projeto/PRGroup/series_sinteticas/DU_reeval/streamflow_generator/Kirsch-Nowak_Streamflow_Generator-master/validation/historical/' files{i}(1:(length(files{i})-14)) '-monthly.csv'], Qfinal_monthly{i}); 
end
