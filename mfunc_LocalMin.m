%%% この関数の機能を解読してみる
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 江崎先生による User's guide
%% 解説では 54 ページに書かれている
%% mfunc_VectorList(nodeNumber)　で、活動パターンを作る　出来上がった活動パターン (= vectorList) は　nodeNum = 3 なら
%%%   -1   1  -1   1  -1   1  -1   1
%%%   -1  -1   1   1  -1  -1   1   1
%%%   -1  -1  -1  -1   1   1   1   1
%%% になる  列の数は 2^3 = 8
%% 解説にあるように、ある活動パターンに対して、それと一か所が変化（+1 -> -1 または -1 -> +1） したパターンを隣接パターンとする
%% ノードの数が N なら、一か所が変化する場合の数は N 通りあるので一つの活動パターンに N 個の隣接パターンがある
%% mfunc_VectorIndex では、隣接パターンを決定しやすくするための index として vectorIndex を作成している
%% EnergyIndex は、普通に左から 1, 2, 3, ... と順番につけた index
%% （つづく）
function [LocalMinIndex, BasinGraph, AdjacentList] = mfunc_LocalMin(nodeNumber, E)
% Find local minimum points
% Get vector list
vectorList = mfunc_VectorList(nodeNumber);
% Calculate vector index (See mfunc_VectorIndex())
[vectorIndex EnergyIndex] = mfunc_VectorIndex(vectorList);

%% ここでは、あるパターンと隣接するパターンのリストを作成している
%% 解説にあるように、ある活動パターンに対して、それと一か所が変化（+1 -> -1 または -1 -> +1） したパターンを隣接パターンとする
%% 隣接パターンは、ノードの数と同じ個数存在する
%% vectorIndex-1 と 1, 2, 4, ...  それぞれ XOR を取ると、それらの値が隣接パターンの index になる
%% それらの index は、EnergyIndex 左から順に数えた値になる
%% （つづく）
% Calculate local minimum points
% Make nearest neighbor index matrix
% Each data is one index different from original position
NeighborMatrix = vectorIndex;
Loc = 1;
for i=1:nodeNumber
    tmp = bitxor(vectorIndex-1, Loc)+1;
    NeighborMatrix = [NeighborMatrix, tmp];
    Loc = Loc * 2;  % 1 bit shift
end
%% N=3 なら
%% NeighborMatrix =
%%   1   2   3   5　　　パターン１と隣接するのはパターン２，３，５
%%   5   6   7   1　　　パターン５と隣接するのはパターン６，1，1
%%   3   4   1   7　　　パターン３と隣接するのはパターン4，1，7
%%   7   8   5   3　　　パターン７と隣接するのはパターン8，5，3
%%   2   1   4   6　　　パターン２と隣接するのはパターン1，4，6
%%   6   5   8   2　　　パターン６と隣接するのはパターン5，8，2
%%   4   3   2   8　　　パターン４と隣接するのはパターン3，2，8
%%   8   7   6   4　　　パターン８と隣接するのはパターン7，6，4

%% ここの処理で、Matrix をリストに変換？　値は変化しない(Octave では）
% Calculate adjacent list for adjacency matrix
AdjacentList = EnergyIndex(NeighborMatrix);

%%　隣接パターンが保持するエネルギーを取り出す　エネルギー E は N = 3 なら 2^3 = 8 個の値を持つ 
% Make energy list at each position including neighbors
NeighborEnergy = E(AdjacentList);
%%  E=[8,7,6,5,4,3,2,1] なら
%% >> E(AdjacentList)
%%   8   7   6   4　　　パターン１のエネルギーは８、隣接パターンのエネルギーは　７、６、４
%%   4   3   2   8
%%   6   5   8   2
%%   2   1   4   6
%%   7   8   5   3
%%   3   4   1   7
%%   5   6   7   1
%%   1   2   3   5　　　パターン８のエネルギーは１、隣接パターンのエネルギーは　２、３、５　　だから８は極小　
%% この場合パターン８だけが極小

% Calculate minimum energy and its index
[EnergyMin, tmp] = min(NeighborEnergy, [], 2);

% If the direction of BasinGraph is the same as its own, it is local
% minimum.
idx = (1:length(E))';
BasinGraph = AdjacentList(sub2ind(size(AdjacentList), idx, tmp));
LocalMinIndex = idx(idx == BasinGraph);
BasinGraph = [idx, BasinGraph];

% Calculate digraph
G = digraph(BasinGraph(:,1), BasinGraph(:,2));

% List up member of basin
% How many connectivity?
conn = conncomp(G, 'type', 'weak');
[dummy, memval] = sort(conn(LocalMinIndex));
member = zeros(length(idx), 1, 'double');
for i=1:length(LocalMinIndex)
    member(conn == i) = LocalMinIndex(memval(i));
end
BasinGraph = [BasinGraph, member];

% Show figures
for i=1:length(idx)
    Name{i} = num2str(idx(i));
end
% 2D Graph
figure(100);
h = plot(G, 'Layout', 'force');
h.NodeLabel = Name;
h.NodeCData = E;
h.LineWidth = 1;
c = colorbar;
c.Label.String = 'Energy';
c.Label.FontSize = 16;
ax = gca;
ax.XTickLabel = [];
ax.YTickLabel = [];
title('Local Minima and Basin (2D)', 'FontSize', 16);

% 3D Graph
figure(101);
h = plot(G, 'Layout', 'force');
h.NodeLabel = Name;
h.NodeCData = E;
h.LineWidth = 1;
layout(h, 'force3');
view(3);
h.ZData = E;
c = colorbar;
c.Label.String = 'Energy';
c.Label.FontSize = 16;
ax = gca;
ax.XTickLabel = [];
ax.YTickLabel = [];
title('Local Minima and Basin (3D)', 'FontSize', 16);

%%% この関数の機能を解読してみる
%%% 参考資料：　数理科学2019年6月号51ページ　「エネルギー地形解析」増田直紀先生による解説
%%% 江崎先生による User's guide
%% 解説では 54 ページに書かれている
