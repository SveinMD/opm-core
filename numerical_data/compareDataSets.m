% compare vtu data sets

path = './numerical_data/';
fileBase = 'testCase1';
newton = 'n';
regula_falsi = 'r';
trust_region = 't';
fileExt = '.vtu';

numberOfPoints = 20;
error = zeros(numberOfPoints,1);
nProd = zeros(numberOfPoints,1);
rfProd = zeros(numberOfPoints,1);

for i = 0:numberOfPoints-1

fileNumber = sprintf('%03d',i);
nFile = strcat(path,fileBase,'-',trust_region,'-',fileNumber,fileExt);
rfFile = strcat(path,fileBase,'-',regula_falsi,'-',fileNumber,fileExt);

n = readXmlDataArray(nFile);
rf = readXmlDataArray(rfFile);

%var = genvarname(['diff_t' fileNumber]);
%eval([var '= norm(n-rf)']);
error(i+1) = norm(n(:,1)-rf(:,1));
nProd(i+1) = n(end,1);
rfProd(i+1) = rf(end,1);

end

relErrProd = abs(nProd-rfProd); %./rfProd;