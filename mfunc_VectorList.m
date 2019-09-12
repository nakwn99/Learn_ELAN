%%% この関数は活動パターンを作る
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 江崎先生による User's guide
%%% この関数は、pfunc_03_Accuracy から呼び出され、活動パターン（解説の図2C）に相当する値 vectorList を作る
%%% vectorList は行列になり、行の数はノードの数　列の数は　2 ^ ノードの数
%%% dec2bin は小数 d の 2 進数表現を文字ベクトルとして返します。nodeNum が 3 なら、tmp は (000, 001, ..., 110, 111) を
%%% 転置した縦ベクトルになる（後ろに ' がついているので転置）　これを 8 行 3 列の行列と対応させる
%%% 次に vectorList の容器となる行列を作る　要素は全部 0　行数は tmp と同じ　列数は tmp の要素の文字数と同じ　000, .. ,111 なら ３
%%% vectorList(tmp == '0') = -1;　とすると、tmp の要素の文字で 0 の部分の要素を -1 に変換する
%%% vectorList(tmp == '1') =  1;　とすると、tmp の要素の文字で 1 の部分の要素を  1 に変換する
%%% flipud は、要素の -1 と +1 を入れ替える　' で転置するので、出来上がった vectorList は　nodeNum = 3 なら
%%%   -1   1  -1   1  -1   1  -1   1
%%%   -1  -1   1   1  -1  -1   1   1
%%%   -1  -1  -1  -1   1   1   1   1
%%% になる　解説の図2C の活動パターンと一致している

%This function creates a state vectorList: NumNode x 2^NumNode1
function vectorList = mfunc_VectorList(nodeNum)
node = 2^nodeNum;
tmp = dec2bin((0:node-1)');
vectorList = zeros(node, nodeNum, 'double');
vectorList(tmp == '0') = -1;
vectorList(tmp == '1') = 1;
vectorList = flipud(vectorList');

% vectorList = -ones(nodeNum, 1);
% 
% for ite = 1:nodeNum % ite = num of the active nodes
%     
%     activeNodeList = nchoosek(1:nodeNum,ite); % obtain a list of possible configurations
%     
%     vectorListTemp = - ones(nodeNum, size(activeNodeList,1));
%     
%     for ite2 = 1:size(activeNodeList,1)
%         vectorListTemp(activeNodeList(ite2,:),ite2) = 1;
%     end
%     
%     vectorList = [vectorList,vectorListTemp];
%     
% end
% clear ite ite2 activeNodeList vectorListTemp

%%% この関数は活動パターンを作る
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 江崎先生による User's guide
