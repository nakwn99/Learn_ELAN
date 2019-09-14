%%% この関数の機能を説明してみる
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説

%%% この関数では、行列として与えられたオリジナルデータ（時系列）を二値化 (+1, -1) している
%%% 行はデータの発生源になったそれぞれのノード、列は各時間での値
%%% mean(A,2) は各行の平均値をもつ列ベクトル　
%%% 各行について平均値を計算し、行ごとにそれぞれの要素の値から引き算した値の、正、負で二値化
%%% threshold 閾値を引数として与えている　この引数の値は０でもよい
%%% ones はすべての要素が 1 のベクトル・行列　average, threshold にこれをかけ算することで
%%% 元のデータと行数 nodeNumber、列数 datalength が等しい行列にしている
%%% eps は、MATLAB の数値の最小の刻み　これを足すことで sign の値が
%%% ０にならないようにする

%1 This function binarizes original data with respect to a threshold value.
function [binarizedData] = pfunc_01_Binarizer(originalData, threshold)
[nodeNumber,dataLength] = size(originalData);
average = mean(originalData,2);
binarizedData = sign(originalData - average*ones(1,dataLength) - (threshold+eps)*ones(nodeNumber,dataLength));%should not be equal to zero
end

%%% この関数の機能を説明してみる　　　
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 江崎先生による User's Guide
