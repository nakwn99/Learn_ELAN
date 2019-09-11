%%% この関数の機能を説明してみる
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 江崎先生による User's guide

%%% pfunc_02 で尤度を最大にする h, J を求める際に必要な、モデルの期待値（平均）、分散共分散行列（パラメーターは h, J）を計算している
%%% prob は、解説の式 (4) で求められる値 P  　mfunc_StateProb で計算している
%%% 期待値は　（値 x その値に対応する確率）を全部足すことで求められる
%%% 一カラムずつ（各時間ごとに）計算して、さらに 2^nodeNumber をかけ算している　
%%% これは江崎先生による user's guide の Fig.3 に書かれているように、複数のノードの状態で構成される活動パターン (各要素は -1 か 1) を
%%% 一つの数値に換算するために行っている
%%% MATLAB の . ピリオドは、要素単位の演算を表す　

% This function calculates mean and correlation with respect to the pairwise MEM
% with given h and J.
function [modelMean, modelCorrelation] = mfunc_ModelMeanCorrelation(h,J)
nodeNumber = size(h,1);
vectorlist =  mfunc_VectorList(nodeNumber);
prob = mfunc_StateProb(h,J);

modelMean = mean((vectorlist.*(ones(nodeNumber,1)*prob')),2)*2^nodeNumber;
modelCorrelation = (vectorlist.*(ones(nodeNumber,1)*(sqrt(prob'))))*(vectorlist.*(ones(nodeNumber,1)*(sqrt(prob'))))';
end

%%% この関数の機能を説明してみる
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 江崎先生による User's guide
