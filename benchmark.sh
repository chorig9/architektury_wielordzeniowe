
sum=0
max=0
min=10000000

bin=$1

for run in {1..10}; do
	 iter_result=$($1);

	echo $iter_result

         sum=$((sum + iter_result));

         max=$((max > iter_result ? max : iter_result));
         min=$((min < iter_result ? min : iter_result));
done

#echo "sum: $sum"
printf "mean: "
echo "$sum/10" | bc -l
echo "var: $((($max - $min)/2))"
echo "max: $max"
echo "min: $min"
