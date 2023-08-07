
% This script shows how to use the streamflow generator available in the
% Kirsch-Nowak_Streamflow_Generator repository on a test dataset 
% (the Lower Susquehanna River basin).
%
% Copyright 2017 Matteo Giuliani, Jon Herman and Julianne Quinn
% 
% Post-doc Research Fellow at Politecnico di Milano
% matteo.giuliani@polimi.it
% http://giuliani.faculty.polimi.it
%
% Faculty Member at UC Davis
% jdherman8@gmail.com
%
% Postdoctoral Researcher at Cornell University
% jdq8@cornell.edu
%
% Please refer to README.txt for further information.
%
%
%     This code is free software: you can redistribute 
%     it and/or modify it under the terms of the GNU General Public License 
%     as published by the Free Software Foundation, either version 3 of the 
%     License, or (at your option) any later version.     
% 
%     This code is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.

%% prepare workspace
clear all
clc

% load multi-site observations of daily streamflow
Qdaily = load('E:/UNB/Projeto/PRGroup/series_sinteticas/streamflow_generator/Kirsch-Nowak_Streamflow_Generator-master/data/Qdaily.txt');
%CÓDIGO ABAIXO DEVE SER VERIFICADO EM CASO DE INCLUIR EVAPORAÇÃO
%Qdaily = Qdaily(:,1:4); % columns 4 and 5 the same; remove column 5

% make normally distributed evaporation log-normal like flows
% (monthly_gen.m takes the log of Qdaily to make all columns normally
% distributed)
%CODIGO ABAIXO DEVE SER VERIFICADO EM CASO DE INCLUIR EVAPORAÇÃO
Qdaily(:,6) = exp(Qdaily(:,6));
Qdaily(:,7) = exp(Qdaily(:,7));
Qdaily(:,8) = exp(Qdaily(:,8));
Qdaily(:,9) = exp(Qdaily(:,9));

sites = {'qDescoberto', 'qSantaMaria', 'qParanoa', 'qTorto_Bananal','qCorumba', 'evapDescoberto','evapSantaMaria','evapParanoa','evapCorumba'};
Nyears = size(Qdaily,1)/365;
Nsites = size(Qdaily,2);


%% Kirsch + Nowak generation
clc
% here you insert the number of realizations and the extension of each
% simulation
num_realizations = [1000];
num_years = [40];
dimensions = {'-1000x40'};
% directory to write output to
mkdir('E:/UNB/Projeto/PRGroup/series_sinteticas/streamflow_generator/Kirsch-Nowak_Streamflow_Generator-master/validation/synthetic');
for k=1:length(num_realizations)
    Qd_cg = combined_generator(Qdaily, num_realizations(k), num_years(k) );
    % back-transform evaporation
    %CODIGO ABAIXO DEVE SER AVALIADO EM CASO DE INCLUIR EVAPORAÇÃO
    Qd_cg(:,:,6) = log(Qd_cg(:,:,6));
    Qd_cg(:,:,7) = log(Qd_cg(:,:,7));
    Qd_cg(:,:,8) = log(Qd_cg(:,:,8));
    Qd_cg(:,:,9) = log(Qd_cg(:,:,9));
    
    % write simulated data to file
    for i=1:Nsites
        q_ = [];
        for j=1:num_realizations(k)
            qi = nan(365*num_years(k),1);
            qi(1:size(Qd_cg,2)) = Qd_cg(j,:,i)';
            q_ = [q_ reshape(qi,365,num_years(k))];
        end
        Qd2(:,i) = reshape(q_(:),[],1);
        saveQ = reshape(Qd2(:,i), num_years(k)*365, num_realizations(k))';
        dlmwrite(['E:/UNB/Projeto/PRGroup/series_sinteticas/streamflow_generator/Kirsch-Nowak_Streamflow_Generator-master/validation/synthetic/' sites{i} dimensions{k} '-daily.csv'], saveQ);
    end
    synMonthlyQ = convert_data_to_monthly(Qd2);
    %REAVALIAR O CÓDIGO ABAIXO EM CASO DE INCLUIR EVAPORAÇÃO - OBSERVAR O
    %NUMERO DA COLUNA QUE O DADO ESTA!!
    % divide evaporation by 86400 (s/day) to get total monthly evap in mm/month
    synMonthlyQ{6} = synMonthlyQ{6}/86400;
    synMonthlyQ{7} = synMonthlyQ{7}/86400;
    synMonthlyQ{8} = synMonthlyQ{8}/86400;
    synMonthlyQ{9} = synMonthlyQ{9}/86400;
        
    for i=1:Nsites
        saveMonthlyQ = reshape(synMonthlyQ{i}',12*num_years(k),num_realizations(k))';
        dlmwrite(['E:/UNB/Projeto/PRGroup/series_sinteticas/streamflow_generator/Kirsch-Nowak_Streamflow_Generator-master/validation/synthetic/' sites{i} dimensions{k} '-monthly.csv'], saveMonthlyQ);
    end
    dlmwrite(['E:/UNB/Projeto/PRGroup/series_sinteticas/streamflow_generator/Kirsch-Nowak_Streamflow_Generator-master/validation/synthetic/Qdaily' dimensions{k} '.csv'], Qd2);
    clear Qd2;
end
