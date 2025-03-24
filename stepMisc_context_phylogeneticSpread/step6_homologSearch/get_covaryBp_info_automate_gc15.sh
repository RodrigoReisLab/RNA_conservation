#! /bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --error=getPercentCovariation_gc15-noCacofold.err
#SBATCH --output=getPercentCovariation_gc15-noCacofold.out
#SBATCH --job-name=getPercentCovariation_gc15-noCacofold
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

echo -e "cluster\ttotalBP\tcovaryBP\tperCovary\n" > rscape_gc15_twoTest_covariation_results.txt;

for i in `find ./ -maxdepth 2 -name "*_gc15_rscape.power"`
do
	output=$(grep "^# BPAIRS" ${i} | sed -n '1p;$p')
	totalBP=$(echo "$output" | head -n1 | awk -F' ' '{print $NF}')
	covaryBP=$(echo "$output" | tail -n1 | awk -F' ' '{print $NF}')
	if [ $totalBP -gt 0 ]
	then
		perCovary=$(bc -l <<<"scale=2;$covaryBP/$totalBP*100")
		echo -e "$i\t$totalBP\t$covaryBP\t$perCovary";
	else
		totalBP=0;
		covaryBP=0;
		perCovary=0;
		echo -e "$i\t$totalBP\t$covaryBP\t$perCovary";
	fi
done >> rscape_gc15_twoTest_covariation_results.txt
