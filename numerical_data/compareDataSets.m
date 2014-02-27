% compare vtu data sets

clear all;

path = './';
fileBase = 'testCase1';
tend = '2_000000';
tstep = '0_100000';
solv_list = ['t ';'ta';'b ';'i '];
nsol = length(solv_list);
solv_list = cellstr(solv_list);
regula_falsi = 'r';
fileExt = '.vtu';

numberOfPoints = 20;
error = zeros(numberOfPoints,1);
nProd = zeros(numberOfPoints,1);
rfProd = zeros(numberOfPoints,1);

fullName = @(solver,number) strcat(path,fileBase,'-s-',solver,'-T-',tend,'-t-',tstep,'-',number,fileExt);

for i = 1:nsol
    solv_char = solv_list{i};
    var = genvarname(['diff_s_' solv_char]);
    eval([var '= [];']);
end

for i = 0:numberOfPoints-1
    fileNumber = sprintf('%03d',i);
    rfFile = fullName('r',fileNumber);
    rf = readXmlDataArray(rfFile);
    for i = 1:nsol
        solv_char = solv_list{i};
        solverFile = fullName(solv_char,fileNumber);
        solverResults = readXmlDataArray(solverFile);
        
        var = genvarname(['diff_s_' solv_char]);
        eval([var '= [' var ' norm(solverResults-rf)];']);
        
        %error(i+1) = norm(solverResults(:,1)-rf(:,1));
        %nProd(i+1) = solverResults(end,1);
        %rfProd(i+1) = rf(end,1);
    end
end

for i = 1:nsol
    solv_char = solv_list{i};
    var = genvarname(['diff_s_' solv_char]);
    eval(['save ' var '.data ' var ' -ascii']);
end

%relErrProd = abs(nProd-rfProd); %./rfProd;