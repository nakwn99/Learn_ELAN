%%% Octave �œ��������������c�����t�@�C��
%%% ���̊֐��̋@�\����ǂ��Ă݂�
%%% �Q�l�����F�@�����Ȋw2019�N6����51�y�[�W�@�u�G�l���M�[�n�`��́v���c���I�搶�ɂ����
%%% �]��搶�ɂ�� User's guide
%% ����ł� 54 �y�[�W�ɏ�����Ă���
%% mfunc_VectorList(nodeNumber)�@�ŁA�����p�^�[�������@�o���オ���������p�^�[�� (= vectorList) �́@nodeNum = 3 �Ȃ�
%%%   -1   1  -1   1  -1   1  -1   1
%%%   -1  -1   1   1  -1  -1   1   1
%%%   -1  -1  -1  -1   1   1   1   1
%%% �ɂȂ�  ��̐��� 2^3 = 8
%% ����ɂ���悤�ɁA���銈���p�^�[���ɑ΂��āA����ƈꂩ�����ω��i+1 -> -1 �܂��� -1 -> +1�j �����p�^�[����אڃp�^�[���Ƃ���
%% �m�[�h�̐��� N �Ȃ�A�ꂩ�����ω�����ꍇ�̐��� N �ʂ肠��̂ň�̊����p�^�[���� N �̗אڃp�^�[��������
%% mfunc_VectorIndex �ł́A�אڃp�^�[�������肵�₷�����邽�߂� index �Ƃ��� vectorIndex ���쐬���Ă���
%% EnergyIndex �́A���ʂɍ����� 1, 2, 3, ... �Ə��Ԃɂ��� index
%% �i�Â��j

function [LocalMinIndex, BasinGraph, AdjacentList] = mfunc_LocalMin_for_oct2py(nodeNumber, E)
% Find local minimum points
% Get vector list
vectorList = mfunc_VectorList(nodeNumber);
% Calculate vector index (See mfunc_VectorIndex())
[vectorIndex EnergyIndex] = mfunc_VectorIndex(vectorList);

%% �����ł́A����p�^�[���Ɨאڂ���p�^�[���̃��X�g���쐬���Ă���
%% ����ɂ���悤�ɁA���銈���p�^�[���ɑ΂��āA����ƈꂩ�����ω��i+1 -> -1 �܂��� -1 -> +1�j �����p�^�[����אڃp�^�[���Ƃ���
%% �אڃp�^�[���́A�m�[�h�̐��Ɠ��������݂���
%% (vectorIndex-1) �� 1, 2, 4, ...  ���ꂼ�� XOR �����ƁA�����̒l���אڃp�^�[���� index �ɂȂ�悤�ɍI���ɍH�v����Ă���
%% �i�Â��j
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
%% N=3 �Ȃ�
%% NeighborMatrix =
%%   1   2   3   5�@�@�@�p�^�[���P�ivectorIndex ���P�j�Ɨאڂ���̂̓p�^�[���Q�C�R�C�T�i������̐����� vectorIndex �̂Q�C�R�C�T�j
%%   5   6   7   1�@�@�@�p�^�[���T�ivectorIndex ���j�Ɨאڂ���̂̓p�^�[���U�C�V�C1
%%   3   4   1   7�@�@�@�p�^�[���R�ivectorIndex ���j�Ɨאڂ���̂̓p�^�[��4�C1�C7
%%   7   8   5   3�@�@�@�p�^�[���V�ivectorIndex ���j�Ɨאڂ���̂̓p�^�[��8�C5�C3
%%   2   1   4   6�@�@�@�p�^�[���Q�ivectorIndex ���j�Ɨאڂ���̂̓p�^�[��1�C4�C6
%%   6   5   8   2�@�@�@�p�^�[���U�ivectorIndex ���j�Ɨאڂ���̂̓p�^�[��5�C8�C2
%%   4   3   2   8�@�@�@�p�^�[���S�ivectorIndex ���j�Ɨאڂ���̂̓p�^�[��3�C2�C8
%%   8   7   6   4�@�@�@�p�^�[���W�ivectorIndex ���j�Ɨאڂ���̂̓p�^�[��7�C6�C4

%% �����̏����ŁAMatrix �����X�g�ɕϊ��H�@�l�͕ω����Ȃ� EnergyIndex �͂ӂ��� (1,2,3,4,5,6,7,8) �ƕ��ׂĂ���
%% �i�Â��j
% Calculate adjacent list for adjacency matrix
AdjacentList = EnergyIndex(NeighborMatrix);

%%�@�אڃp�^�[�����ێ�����G�l���M�[�����o���@�G�l���M�[ E �� N = 3 �Ȃ� 2^3 = 8 �̃p�^�[�������ꂼ��G�l���M�[�l������ 
% Make energy list at each position including neighbors
NeighborEnergy = E(AdjacentList);
%%  E=[8,7,6,5,4,3,2,1] �Ȃ�
%% >> E(AdjacentList)
%%   8   7   6   4�@�@�@�p�^�[���P�̃G�l���M�[�͂W�A�אڃp�^�[���̃G�l���M�[�́@�V�A�U�A�S
%%   4   3   2   8
%%   6   5   8   2
%%   2   1   4   6
%%   7   8   5   3
%%   3   4   1   7
%%   5   6   7   1
%%   1   2   3   5�@�@�@�p�^�[���W�̃G�l���M�[�͂P�A�אڃp�^�[���̃G�l���M�[�́@�Q�A�R�A�T�@�@������W�͋ɏ��@
%% ���̏ꍇ�p�^�[���W�������ɏ���ԁ@�i�Â��j
 
%% EnergyMin �ɍŏ��̒l�Atmp �ɂ��̒l�ɑ΂��� index�@��̗�ł� tmp(8) = 1
%% �i�Â��j
% Calculate minimum energy and its index
[EnergyMin, tmp] = min(NeighborEnergy, [], 2);

%% Basingraph = �n�`�}�@
%% sub2ind([�s���A��], �s�A��j�́A���̍s��̗v�f�������i�s�A��j��P��̃C���f�b�N�X���l�ɕϊ�����
%% �l�̐U����́A���オ�P�ŁA�s�����i�������j�ɂQ�C�R�A�c�@�Ɛ����Ă����@��ԉ��̎��́A���̉E�̍s�̈�ԏ�
%% sub2ind(size(AdjacentList), idx, tmp) �ŁAAdjacentList �s��́A�e�s�ň�ԃG�l���M�[���Ⴂ�v�f���w������
%% BasnGraph �́A�e�s�ň�ԃG�l���M�[���Ⴂ�p�^�[�������� vectorIndex �̒l = ���̃p�^�[���́A���̍s�ɂ�����ɏ��@
%% LocalMinIndex = idx(idx == BasinGraph);�@����́ABasinGraph �̗v�f�̒l�� idx (��F1:8 �̏c�x�N�g���j�ƈ�v����
%% �s��I��ł���@��v���Ă��遁���̍s�̍Œ�l���אڃp�^�[���ł͂Ȃ��������g���ɏ����
%% �ɏ��ł���s�i�����p�^�[���j�͕������݂��Ă��悢
%% �i�Â��j

% If the direction of BasinGraph is the same as its own, it is local
% minimum.
idx = (1:length(E))';
BasinGraph = AdjacentList(sub2ind(size(AdjacentList), idx, tmp));
LocalMinIndex = idx(idx == BasinGraph);
BasinGraph = [idx, BasinGraph];

%%% Octave �œ��������������c�����t�@�C��
%%% ���̊֐��̋@�\����ǂ��Ă݂�
%%% �Q�l�����F�@�����Ȋw2019�N6����51�y�[�W�@�u�G�l���M�[�n�`��́v���c���I�搶�ɂ����
%%% �]��搶�ɂ�� User's guide
%% ����ł� 54 �y�[�W�ɏ�����Ă���