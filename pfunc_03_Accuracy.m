Octave で確かめながらこの関数の機能を説明してみる
資料：数理科学2019年6月号51ページ「エネルギー地形解析」増田直紀先生による解説　江崎先生による User's guide

この関数は、引数である二値化データ（実測値からの）をすべて一度に用いて、それに対応する確率、エントロピーを計算する
また pfunc_02 で計算したパラメータ h, J を用いて、モデルに対応する確率、エントロピー（-p log P）を計算する　それには mfunc_StateProb(h,J)　を呼び出している
mfunc_VectorList 関数を呼び出して、vectorList をまず作っている
この vectorList は　例えば nodeNum = 3 なら
  -1   1  -1   1  -1   1  -1   1
  -1  -1   1   1  -1  -1   1   1
  -1  -1  -1  -1   1   1   1   1
のような行列を作る　各列が　0, 1, 2, 3, 4, 5, 6, 7 に相当する　増田先生の解説では図2の (C) の活動パターンに相当する
実測データからの確率、エントロピーの計算を vectorList 　一列　ごとに繰り返している（行ではない）
増田先生の解説では図2の (C) に相当する
vectorList から一列取り出し sampleVec とする　これを実測値の列の数だけ横に並べた行列に変換している　そこから vectorList を 
sum(abs(引き算)) している　これは「sampleVecMatrix と 二値化した実測値が完全に一致している」と 0 になり、一致しない部分が多いほど大きな値になる
find では、detector の値が == 0 になるインデックスを見つけている　size（リスト, 2）は列の数　 
だから numSample は、 「取り出した sampleVec と二値化した実測値が完全に一致している」回数になる
それを実測値の列の数で割ると「実測値から求めた、そのパターンの出現確率」P_emp(\sigma) になる
エントロピーは、情報エントロピーの基本的な式　-p log(p) で計算している
（つづく）
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

「N-th order model」に相当と書かれている　＝　ノードが N 個あり、それらが全てお互いに相互作用している状態 
実測値はそういうとても複雑なシステムから生じている可能性がある　全てをカバーするために、ここでは全ての実測値を一度に用いて計算する
（つづく）
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

この下でも実測値 binarizedData を用いて計算しているが、上の計算（実測値すべてを用いた計算）とは異なり
それぞれの行（ノード）のデータを独立に用いて確率を計算する　上には independent MEM と書かれている
ノードとノードの相互作用は考慮しない　＝　各ノード単独のエネルギー　これが全エネルギーよりずっと小さければ、「相互作用が強い」
binarizedData に 1 を足して 2 で割ると、-1 は 0 に、1 は 1 になる　mean(~, 2) は行方向に平均値を計算　この場合は、それぞれの行に +1 が
存在する確率（個数 / 列数）を計算することになる　probInactive は -1 の確率　　　　
各ノードごとに、　これらの値が計算され、縦に並んだベクトルになる
numVec は列の数　ones(1, numVec) を掛け算すると、縦ベクトルである probActive (or Inactive) を横向きに numVec だけ繰り返し並べることになる
probActive はノードの数と同じ数の要素を持つ縦ベクトルなので、binarizedData と同じサイズの行列ができる
(vectorList+1)/2 .* activeMatrix　では、vectorList で 1 だった部分に probActive の値が入る　
(1-vectorList)/2 .* inactiveMatrix　では、残りの -1 だった部分に probInactive の値が入る
prod(probMatrix_temp) は、各列の要素の積を求める　
その値は、「それぞれの行（ノード）のみのデータから求めた確率を元にして計算した、ある活動パターンが出現する確率」
最初のセクションで求めた値は、「データ全体」から「ある活動パターンが出現する確率」を求めているので、
二つの値の差は「行同士の相互作用」に相当する値になる
転置がついているので、数値が実測値の列の数だけ、縦に並んだベクトルが結果として得られる　それが prob1　これも -p log(p) でエントロピーに換算
（つづく）
%% Evaluate 1st-order model%
% 1st-order model:
% E(V) = -sum h_i sigma_i  
probActive = mean((1+binarizedData)/2, 2); % probability of sigma being +1
probInactive = 1-probActive;
activeMatrix = probActive * ones(1,numVec);
inactiveMatrix = probInactive * ones(1,numVec);
probMatrix_temp = (vectorList+1)/2 .* activeMatrix + (1-vectorList)/2 .* inactiveMatrix;
probMatrix = prod(probMatrix_temp)';
prob1 = probMatrix;
S1 = - prob1' * log2(prob1);

この下では、モデルのほうの確率 P_model(\sigma) を計算している　式 (5) は二つの項（単独、相互作用）がある二次モデルを示している　
mfunc_StateProb(h, J) を呼び出して、それを計算する　（つづく）
%% Evaluate 2nd-order model
prob2 = mfunc_StateProb(h,J);
S2 = - prob2' * log2(prob2);

SN はすべての実測値から計算した確率を用いた「すべての相互作用を含めて考慮したエントロピー」
S1 は、それぞれのノード（行）のみの実測値から計算した確率を用いた「相互作用を考慮しないエントロピー」　
S2 はモデル（式(5)）から計算したエントロピー
S1 との差を取ることで、相互作用の大きさを残している　それがどれくらい一致しているかを指標として示している
（つづく）
%%accuracy index 
r = (S1-S2)/(S1-SN);

ここでも、一致度を計算する　KL divergence という情報量と関係する指標を用いて計算している




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

%% Octave で確かめながらこの関数の機能を説明してみる
%% 資料：数理科学2019年6月号51ページ「エネルギー地形解析」増田直紀先生による解説　江崎先生による User's guide
