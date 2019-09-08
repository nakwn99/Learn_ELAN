この関数の機能を説明してみる
資料：数理科学2019年6月号51ページ「エネルギー地形解析」増田直紀先生による解説　江崎先生による User's guide

この関数は、引数である二値化データ（実測値からの）を用いて、それに対応する確率、エントロピーを計算する
また pfunc_02 で計算したパラメータ h, J を用いて、モデルに対応する確率、エントロピー（-p log P）を計算する
mfunc_VectorList 関数を呼び出して、vectorList をまず作っている
この vectorList は



実測データからの確率、エントロピーの計算を vectorList 一行ごとに繰り返している
sampleVec を取り出し、行列に変換している　そこから
% 2016/11/11 by T. Ezaki 
% This function calculates accuracy indices:  
% r = (S1-S2)/(S1-SN), rD=(D1-D2)/D1
%
% probN: empirical prob dist (exact), 
% prob1: independent MEM 
% prob2; pairwise MEM
function [probN, prob1, prob2, rD, r] = pfunc_03_Accuracy(h, J, binarizedData)
[nodeNumber,dataLength] = size(binarizedData);
%% Create VectorList: nodeNumber x 2^nodeNumber
vectorList = mfunc_VectorList(nodeNumber);
%% Calculate empirical probability and its entropy (corresponding to N-th order model)
numVec = size(vectorList,2);
probN = zeros(numVec,1);
for iteVec = 1:numVec   
    sampleVec = vectorList(:,iteVec); % choose (iteVec)-th state
    sampleVecMatrix = sampleVec * ones(1, dataLength); 
    detector = sum(abs(sampleVecMatrix - binarizedData));
    numSample = size(find(detector == 0), 2);
    
    probN(iteVec,1) = numSample/size(binarizedData,2); % obtained empirical prob distribution
    
    % calculate entropy
    if numSample == 0
        SN_temp(iteVec,1) = 0;
    else
        SN_temp(iteVec,1) = -probN(iteVec,1) * log2(probN(iteVec,1));
    end
    
end

SN = sum(SN_temp);



%% Evaluate 1st-order model
%
% 1st-order model:
% E(V) = -sum h_i sigma_i
%  

probActive = mean((1+binarizedData)/2, 2); % probability of sigma being +1
probInactive = 1-probActive;

activeMatrix = probActive * ones(1,numVec);
inactiveMatrix = probInactive * ones(1,numVec);

probMatrix_temp = (vectorList+1)/2 .* activeMatrix + (1-vectorList)/2 .* inactiveMatrix;

probMatrix = prod(probMatrix_temp)';
prob1 = probMatrix;

S1 = - prob1' * log2(prob1);



%% Evaluate 2nd-order model

prob2 = mfunc_StateProb(h,J);
S2 = - prob2' * log2(prob2);



%%accuracy index 
r = (S1-S2)/(S1-SN);

%%accuracy index defined with the KL divergence

D1_Temp = zeros(numVec,1);
D2_Temp = zeros(numVec,1);

for iteVec = 1:numVec
    if probN(iteVec,1) == 0
        D1_Temp(iteVec,1) = 0;
        D2_Temp(iteVec,1) = 0;
    else
        D1_Temp(iteVec,1) = probN(iteVec,1) * log2(probN(iteVec,1)/prob1(iteVec,1));
        D2_Temp(iteVec,1) = probN(iteVec,1) * log2(probN(iteVec,1)/prob2(iteVec,1));
    end
end

D1 = sum(D1_Temp);
D2 = sum(D2_Temp);

rD = (D1-D2)/D1;

