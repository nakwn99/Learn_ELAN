%%% この関数の機能を説明してみる
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 
%%% この関数では、二値化したデータからイジングモデルのパラメーター　h, J を推定している
%%% 増田先生の解説の式（５）にあるように、h はそれぞれの活動パターンがもつエネルギーに対する重み
%%% J は活動パターン同士の相互作用エネルギーに対する重み
%%% これらのパラメーターを用いて計算したモデルの期待値が、実際のデータの平均値と等しくなる、
%%% また計算したモデルの分散共分散行列（相関）の期待値が、実際のデータのそれらと等しくなる
%%% ようにパラメーターを計算している　（平均は、期待値で計算に使う確率が一様になった、特別な場合）
%%% 増田先生の解説では54ページに書かれている

%%% 与えられるデータは、行が、データの発生源であるノード１～nodenumber　
%%% 列が時間（時系列データ）
%%% dataMean=mean(binarizedData,2) は、それぞれの行の平均が並んでいるベクトルになる
%%% dataCorrelation は、分散共分散行列を計算する式に基づいて計算
%%% 普通の行列では　VAR(A) = \frac{1}{N-1} A^t %*% A　　しかし今回のデータ行列は、普通は列の方に割り当てる「データのラベル」が
%%% 行に割り当てられているので、普通と反対になり A %*% A^t　で計算している

%%% h と J を求めるために、最急降下法という方法を用いている。増田先生の解説では54ページに書かれている。
%%% まず h, J を 0 で初期化している。（６）、（７）の式に従って h, J の値を少しずつ更新することを繰り返す。
%%% h では、（データから計算した平均とモデルの期待値の差）* dt を足していく
%%% J では、（データから計算した分散共分散行列とモデルのそれとの差）* dt を計算する　対角成分
%%% （行列の要素の [1,1] [2,2] [3,3] ...）は 0 にする（自分自身との相互作用は考えないので 0）ので、diag(diag(dJ)) （=対角成分だけを残した行列）
%%% を引き算して、対角成分 = 0 にしている

%%% if (sqrt(norm(dJ,'fro')^2 + norm(dh)^2)/nodeNumber/(nodeNumber+1) < permissibleErr)　の部分で、データから計算した値と
%%% モデルの期待値との誤差を評価している。これが十分に小さくなれば OK
%%% norm はベクトル・行列の長さ、大きさのようなもの・長さ、大きさを一般化したもの　L1 ノルム、L2 ノルムとかたくさんの種類がある
%%% 'fro' を付けた場合、フロベニウスノルムと言うものを計算する　これは行列を対象として、すべての要素の二乗を足したものの平方根
%%% sqrt(~ + ~) の部分は、dJ の長さと dh の大きさを二つの辺にした三角形のもう一つの辺の長さを計算しているように解釈できる。
%%% だから dJ, dh がバランスよく十分に小さくなればよいことになる

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

%%% この関数の機能を説明してみる
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
