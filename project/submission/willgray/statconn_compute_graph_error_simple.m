function gErr = statconn_compute_graph_error_simple(lgTest, lgTruth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) [2014] The Johns Hopkins University / Applied Physics Laboratory All Rights Reserved. Contact the JHU/APL Office of Technology Transfer for any additional rights.  www.jhuapl.edu/ott
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lgTruthF = lgTruth+lgTruth';
lgTestF = lgTest+lgTest';

%using binary undirected graphs, no self loops for now TODO
lgTruthFB = single(lgTruthF > 0);
lgTestFB = single(lgTestF > 0);

gErr.scoreFrob = norm(lgTruthFB-lgTestFB,'fro');%  sum(sum(abs(logical(lgTruth)-logical(lgTest))));

%mtxDiff = lgTruthFB - lgTestFB;

% for i = 1:size(mtxDiff,1)
%     for j = 1:size(mtxDiff,2)
%         if i >= j
%             mtxDiff(i,j) = NaN;
%         end
%     end
% end
%nel = size(mtxDiff,1);

gErr.TP = sum(sum(lgTestFB == 1 & lgTruthFB == 1));
gErr.FP = sum(sum(lgTestFB == 1 & lgTruthFB == 0));
gErr.FN = sum(sum(lgTestFB == 0 & lgTruthFB == 1));
gErr.TN = sum(sum(lgTestFB == 0 & lgTruthFB == 0));
gErr.sum = gErr.FP + gErr.FN + gErr.TN + gErr.TP;
gErr.connFound = sum(lgTestFB(:));
gErr.connTrue = sum(lgTruthFB(:));
gErr.nSynapse = length(lgTest); %fixed
gErr.avgSynDegreeTest = mean(sum(lgTestF,2));
gErr.avgSynDegreeTruth = mean(sum(lgTruthF,2));
gErr.prec = gErr.TP/(gErr.TP + gErr.FP);
gErr.rec = gErr.TP/(gErr.TP + gErr.FN);
gErr.scoreF1 = 2 * (gErr.prec * gErr.rec) / (gErr.prec + gErr.rec);
