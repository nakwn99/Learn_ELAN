%%% モデルから推定した h, J を用いて、式（5）に基づいてエネルギーを計算する
%%% 資料：数理科学2019年6月号51ページ「エネルギー地形解析」増田直紀先生による解説　江崎先生による User's guide

function E = mfunc_Energy(h,J)
% Calculate energy

nodeNumber = size(h,1);
vectorList = mfunc_VectorList(nodeNumber);
numVec = size(vectorList,2);

E = -diag((0.5*J*vectorList)' * vectorList) ...
    - sum(((h * ones(1,numVec)) .* vectorList)', 2);

%%% モデルから推定した h, J を用いて、式（5）に基づいてエネルギーを計算する
%%% 資料：数理科学2019年6月号51ページ「エネルギー地形解析」増田直紀先生による解説　江崎先生による User's guide
