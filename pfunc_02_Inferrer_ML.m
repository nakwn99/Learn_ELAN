この関数の機能を説明してみる
参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説

この関数では、二値化したデータからイジングモデルのパラメーター　h, J を推定している
増田先生の解説の式（５）にあるように、h はそれぞれの活動パターンがもつエネルギーに対する重み
J は活動パターン同士の相互作用エネルギーに対する重み
これらのパラメーターを用いて計算したモデルの期待値が、実際のデータの平均値と等しくなる、
また計算したモデルの分散共分散行列（相関）の期待値が、実際のデータのそれらと等しくなる
ようにパラメーターを計算している
増田先生の解説では54ページに書かれている

与えられるデータは、行が、データの発生源であるノード１～nodenumber　
列が時間（時系列データ）



% 2016/11/11 by T. Ezaki 
% This function infers h and J by maximum likelihood estimation. 
% binarizedData: nodeNumber x tmax

function [h,J] = pfunc_02_Inferrer_ML(binarizedData)

[nodeNumber,dataLength] = size(binarizedData);
% tuning parameters
iterationMax = 5000000;
permissibleErr = 0.00000001;
dt = 0.2;


dataMean = mean(binarizedData,2);
dataCorrelation = (binarizedData*binarizedData')/dataLength;

h = zeros(nodeNumber,1);
J = zeros(nodeNumber);

for t=1:iterationMax
    [modelMean, modelCorrelation] = mfunc_ModelMeanCorrelation(h,J);
    dh = dt * (dataMean-modelMean);
    dJ = dt *  (dataCorrelation-modelCorrelation);
    dJ = dJ - diag(diag(dJ));% Jii = 0
    h = h + dh;
    J = J + dJ;
    
    if (sqrt(norm(dJ,'fro')^2 + norm(dh)^2)/nodeNumber/(nodeNumber+1) < permissibleErr)
        break;
    end
end
end
