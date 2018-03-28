ClusterAnalysis := module()
	option package;
	description "a package for cluster analysis";

	export ClusterPlot, ConfusionMatrix, ConfusionTable, KMeans;
	local CentralPointPlot, MVEE;

$include "ClusterPlot.mm"
$include "ConfusionMatrix.mm"
$include "KMeans.mm"
$include "MVEE.mm"
end module: