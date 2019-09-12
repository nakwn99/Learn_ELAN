%%% この関数の機能を説明してみる
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 江崎先生による User's guidefunction probMEM = mfunc_stateProb(h,J)

%%% この関数は、モデルのほうの確率 P_model(\sigma) を計算している　式 (5) は二つの項（単独、相互作用）がある二次モデルを示している　
%%% パラメータ h, J は pfunc_02 で推定してあるものを使う
%%% mfunc_VectorList で vectorList を作る　出来上がった vectorList は　nodeNum = 3 なら
%%%   -1   1  -1   1  -1   1  -1   1
%%%   -1  -1   1   1  -1  -1   1   1
%%%   -1  -1  -1  -1   1   1   1   1
%%% になる
%%% Z は、式（4）の分母　分配関数と呼ばれる関数の計算になる　これで割り算することで確率になる
%%% 


nodeNumber = size(h,1);
vectorList = mfunc_VectorList(nodeNumber);
numVec = size(vectorList,2);

Z = sum(exp(-(-diag((0.5*J*vectorList)' * vectorList) - sum(((h * ones(1,numVec)) .* vectorList)', 2))),1);
probMEM = exp(-(-diag((0.5*J*vectorList)' * vectorList) - sum(((h * ones(1,numVec)) .* vectorList)', 2))) ./ Z;

%%% この関数の機能を説明してみる
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 江崎先生による User's guide


