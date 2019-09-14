%%% この関数は、pfunc_02 で推定してあるパラメータ h, J を用いて、モデルのほうの確率 P_model(\sigma) を計算している
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 江崎先生による User's guide

%%% この関数は、pfunc_02 で推定してあるパラメータ h, J を用いて、モデルのほうの確率 P_model(\sigma) を計算している　
%%% h, J は実測データによって決まるので、ここで計算される値 probMEM （解説では式（４））も実測データを反映して決まる
%%% mfunc_VectorList で 活動パターン (= vectorList) を作る　出来上がった活動パターン (= vectorList) は　nodeNum = 3 なら
%%%   -1   1  -1   1  -1   1  -1   1
%%%   -1  -1   1   1  -1  -1   1   1
%%%   -1  -1  -1  -1   1   1   1   1
%%% になる
%%% Z は、式（4）の分母　分配関数と呼ばれる関数の計算になる　これで割り算することで確率になる
%%% probMEM は解説の式（4）　この値を返す

function probMEM = mfunc_stateProb(h,J)

nodeNumber = size(h,1);
vectorList = mfunc_VectorList(nodeNumber);
numVec = size(vectorList,2);

Z = sum(exp(-(-diag((0.5*J*vectorList)' * vectorList) - sum(((h * ones(1,numVec)) .* vectorList)', 2))),1);
probMEM = exp(-(-diag((0.5*J*vectorList)' * vectorList) - sum(((h * ones(1,numVec)) .* vectorList)', 2))) ./ Z;

%% この関数は推定した h, J を用いモデルの確率 P_model(\sigma) を計算
%%% この関数は、pfunc_02 で推定してあるパラメータ h, J を用いて、モデルのほうの確率 P_model(\sigma) を計算している
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 江崎先生による User's guide


