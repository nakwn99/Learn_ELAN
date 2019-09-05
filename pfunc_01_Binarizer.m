この関数では、行列として与えられたオリジナルデータ（時系列）を二値化 (+1, -1) している
行はデータの元になったそれぞれのノード、列は各時間での値
mean(A,2) は各行の平均値をもつ列ベクトル　
各行について平均値を計算し、行ごとに引き算した値の正、負で二値化
threshold 閾値を引数として与えている　この引数の値は０でもよい
ones はすべての要素が 1 のベクトル・行列　threshold にこれをかけ算することで行列にしている
eps は、MATLAB の数値の最小の刻み　これを足すことで sign の値が０にならないようにする

%1 This function binarizes original data with respect to a threshold value.
function [binarizedData] = pfunc_01_Binarizer(originalData, threshold)
[nodeNumber,dataLength] = size(originalData);
average = mean(originalData,2);
binarizedData = sign(originalData - average*ones(1,dataLength) - (threshold+eps)*ones(nodeNumber,dataLength));%should not be equal to zero
end
